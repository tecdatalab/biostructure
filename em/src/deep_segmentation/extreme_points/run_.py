from datetime import datetime
from pathlib import Path

import fire
import torch

from torch.cuda.amp import GradScaler, autocast

import ignite
import ignite.distributed as idist
from ignite.contrib.engines import common
from ignite.engine import Engine, Events
from ignite.handlers import Checkpoint
from ignite.metrics import DiceCoefficient, IoU, Accuracy, Loss, ConfusionMatrix
from ignite.utils import manual_seed, setup_logger
from ignite.contrib.handlers.clearml_logger import (
    ClearMLLogger,
    ClearMLSaver,
    GradsHistHandler,
    GradsScalarHandler,
    WeightsHistHandler,
    WeightsScalarHandler,
    global_step_from_engine,
)

from clearml import Task
from matplotlib import pyplot as plt
plt.style.use(['science','no-latex'])
from torchsummary import summary

from SegmentationAgent import SegmentationAgent

import os
import glob



def training(local_rank, config):

    rank = idist.get_rank()
    manual_seed(config["seed"] + rank)
    device = idist.device()

    logger = setup_logger(name="3D-Unet-Training", distributed_rank=local_rank)

    log_basic_info(logger, config)

    output_path = config["output_path"]
    if rank == 0:
        if config["stop_iteration"] is None:
            now = datetime.now().strftime("%Y%m%d-%H%M%S")
        else:
            now = "stop-on-{}".format({config['stop_iteration']})

        folder_name = "3DUnet-{}-_backend-{}-{}_{}".format(config['depth'],idist.backend(),idist.get_world_size(),now)
        output_path = Path(output_path) / folder_name
        if not output_path.exists():
            output_path.mkdir(parents=True)
        config["output_path"] = output_path.as_posix()
        logger.info("Output path: {}".format(config['output_path']))

        if "cuda" in device.type:
            config["cuda device name"] = torch.cuda.get_device_name(local_rank)


        clearml_logger = ClearMLLogger(project_name="ProteinDeepSegmentation", task_name=folder_name)
        task = Task.current_task()
        # Log hyper parameters
        hyper_params = [
            "model",
            "depth",
            "pool_size",
            "extra_width",
            "batch_size",
            "momentum",
            "input_size",
            "depth",
            "weight_decay",
            "optim",
            "criterion",
            "with_amp",
            "num_epochs",
            "learning_rate",
        ]
        task.connect({k: config[k] for k in hyper_params})


    # Setup Segmentation Agent
    agent = SegmentationAgent(config['val_percentage'], config['test_ids'], config['num_classes'],config['optim'], config['criterion'],
                            config['pool_size'], int(config['extra_width']), config['batch_size'], (config['input_size'],config['input_size'],config['input_size']), 
                            config['data_path'], config['seed'], config['learning_rate'], config['momentum'], config['weight_decay'], 
                            config['depth'], device, config['with_amp'], idist)

    config["num_iters_per_epoch"] = len(agent.train_loader)
    print("Num iters per epoch: {}".format(config["num_iters_per_epoch"]))
    # Print model summary
    summary(agent.model, (2, config['input_size'], config['input_size'], config['input_size'])) 
    print(agent.model)
    # Create trainer for current task
    trainer = create_trainer(agent.model, agent.optimizer, agent.criterion, agent.train_loader.sampler, config, logger)
    #torch.autograd.set_detect_anomaly(True)
    # Let's now setup evaluator engine to perform model's validation and compute metrics
    cm = ConfusionMatrix(num_classes=config['num_classes'])
    metrics={
        'dice': DiceCoefficient(cm), 
        'iou': IoU(cm), 
        'accuracy': Accuracy(), 
        'loss': Loss(agent.criterion)
	}

    # We define two evaluators as they wont have exactly similar roles:
    # - `evaluator` will save the best model based on validation score
    evaluator = create_evaluator(agent.model, metrics=metrics, config=config)
    train_evaluator = create_evaluator(agent.model, metrics=metrics, config=config)
    test_evaluator = create_evaluator(agent.model, metrics=metrics, config=config)

    def run_validation(engine):
        epoch = trainer.state.epoch
        state = train_evaluator.run(agent.train_loader)
        log_metrics(logger, epoch, state.times["COMPLETED"], "Train", state.metrics)
        state = evaluator.run(agent.validation_loader)
        log_metrics(logger, epoch, state.times["COMPLETED"], "Val", state.metrics)

    def run_testing(engine):
        # Evaluate model on test data
        best_models = glob.glob(config["output_path"]+'/best*.pt')
        for i,m in enumerate(best_models):
            agent.model.load_state_dict(torch.load(m))
            state = test_evaluator.run(agent.test_loader)
            log_metrics(logger, config['num_epochs'], state.times["COMPLETED"], "Test model {}".format(os.path.basename(m)), state.metrics)
        

    trainer.add_event_handler(Events.EPOCH_COMPLETED(every=config["validate_every"]) | Events.COMPLETED, run_validation)

    if rank == 0:
        # Setup TensorBoard logging on trainer and evaluators. Logged values are:
        #  - Training metrics, e.g. running average loss values
        #  - Learning rate
        #  - Evaluation train/test metrics
        evaluators = {"training": train_evaluator, "test": evaluator}
        tb_logger = common.setup_tb_logging(output_path, trainer, agent.optimizer, evaluators=evaluators)
        clearml_logger.attach_opt_params_handler(
            trainer, event_name=Events.ITERATION_COMPLETED(every=100), optimizer=agent.optimizer
            )

        clearml_logger.attach(
            trainer, log_handler=WeightsScalarHandler(agent.model), event_name=Events.ITERATION_COMPLETED(every=100)
            )

        clearml_logger.attach(trainer, log_handler=WeightsHistHandler(agent.model), event_name=Events.EPOCH_COMPLETED(every=100))

        clearml_logger.attach(
            trainer, log_handler=GradsScalarHandler(agent.model), event_name=Events.ITERATION_COMPLETED(every=100)
            )

        clearml_logger.attach(trainer, log_handler=GradsHistHandler(agent.model), event_name=Events.EPOCH_COMPLETED(every=100))


    def score_function(engine):
        return engine.state.metrics['dice'][2].item()

    # Store 2 best models by validation accuracy starting from num_epochs / 2:
    best_model_handler = Checkpoint(
        {"model": agent.model},
        get_save_handler(config),
        #filename_prefix="best",
        n_saved=None,
        global_step_transform=global_step_from_engine(trainer),
        #score_name="test_accuracy",
        #score_function=score_function,
    )
    evaluator.add_event_handler(
        Events.COMPLETED(lambda *_: trainer.state.epoch % 2 == 0), best_model_handler
    )
    evaluator.add_event_handler(
        Events.COMPLETED(lambda *_: trainer.state.epoch == config["num_epochs"]), run_testing
    )
    # In order to check training resuming we can stop training on a given iteration
    if config["stop_iteration"] is not None:

        @trainer.on(Events.ITERATION_STARTED(once=config["stop_iteration"]))
        def _():
            logger.info("Stop training on {} iteration".format(trainer.state.iteration))
            trainer.terminate()

    try:
        trainer.run(agent.train_loader, max_epochs=config["num_epochs"])
    except Exception as e:
        logger.exception("")
        raise e

    if rank == 0:
        tb_logger.close()



def run(
    seed=42,
    model='3D-Unet',
    data_path="dataset/dataset_extreme_points.csv",
    output_path="results",
    num_classes=3,
    optim='SGD',
    criterion='CrossEntropy',
    depth=4,
    pool_size=1,
    extra_width=0,
    batch_size=1,
    test_ids = ["emd-8528", "emd-8682", "6uph", "emd-5275", "emd-11686", "4ez4", "emd-8336"],
    val_percentage=0.2,
    input_size=128,
    num_workers=1,
    num_epochs=100,
    learning_rate=0.001,
    momentum=0.9,
    weight_decay=0.0005,
    validate_every=1,
    checkpoint_every=1000,
    backend=None,
    resume_from=None,
    log_every_iters=15,
    stop_iteration=None,
    with_amp=False,
    **spawn_kwargs,
):
    """Main entry to train 3D Unet on Cryo-EM dataset.

    Args:
        seed (int): random state seed to set. Default, 42.
        data_path (str): input dataset csv. Default, "dataset/dataset.csv".
        output_path (str): output path. Default, "results".
        num_classes (int): number of classes. Default 3.
        depth (int): Unet model depth. Default 4.
        batch_size (int): total batch size. Default, 1.
        test_num (int): Number of maps to take appart for testing. Default, 10.
        input_size (int): Size of one side of the cube input 3D array. Default, 128.
        num_workers (int): number of workers in the data loader. Default, 12.
        num_epochs (int): number of epochs to train the model. Default, 24.
        learning_rate (float): learning rate. Default, 0.001.
        validate_every (int): run model's validation every ``validate_every`` epochs. Default, 1.
        checkpoint_every (int): store training checkpoint every ``checkpoint_every`` iterations. Default, 50.
        backend (str, optional): backend to use for distributed configuration. Possible values: None, "nccl". Default, None.
        resume_from (str, optional): path to checkpoint to use to resume the training from. Default, None.
        log_every_iters (int): argument to log batch loss every ``log_every_iters`` iterations.
            It can be 0 to disable it. Default, 15.
        stop_iteration (int, optional): iteration to stop the training. Can be used to check resume from checkpoint.
        with_amp (bool): if True, enables native automatic mixed precision. Default, False.
        **spawn_kwargs: Other kwargs to spawn run in child processes: master_addr, master_port, node_rank, nnodes
    """

    # Get running config
    config = locals()
    config.update(config["spawn_kwargs"])
    del config["spawn_kwargs"]

    with idist.Parallel(backend=backend, **spawn_kwargs) as parallel:
    	parallel.run(training, config)


def create_trainer(model, optimizer, criterion, train_sampler, config, logger):

    device = idist.device()

    # Setup Ignite trainer:
    # - let's define training step
    # - add other common handlers:
    #    - TerminateOnNan,
    #    - handler to setup learning rate scheduling,
    #    - ModelCheckpoint
    #    - RunningAverage` on `train_step` output
    #    - Two progress bars on epochs and optionally on iterations

    with_amp = config["with_amp"]
    scaler = GradScaler(enabled=with_amp)

    def train_step(engine, batch):

        x, y = batch[0], batch[1]

        if x.device != device:
            x = x.to(device, non_blocking=True)
            y = y.to(device, non_blocking=True)

        model.train()

        with autocast(enabled=with_amp):
            y_pred = model(x)
            loss = criterion(y_pred, y)

        optimizer.zero_grad()
        scaler.scale(loss).backward()
        scaler.step(optimizer)
        scaler.update()

        return {
            "batch loss": loss.item(),
        }

    trainer = Engine(train_step)
    trainer.logger = logger

    to_save = {"trainer": trainer, "model": model, "optimizer": optimizer}
    metric_names = [
        "batch loss",
    ]

    common.setup_common_training_handlers(
        trainer=trainer,
        train_sampler=train_sampler,
        to_save=to_save,
        save_every_iters=config["checkpoint_every"],
        save_handler=get_save_handler(config),
        output_names=metric_names if config["log_every_iters"] > 0 else None,
        with_pbars=False,
        clear_cuda_cache=False,
    )

    resume_from = config["resume_from"]
    if resume_from is not None:
        checkpoint_fp = Path(resume_from)
        assert checkpoint_fp.exists(), "Checkpoint '{}' is not found".format(checkpoint_fp.as_posix())
        logger.info("Resume from a checkpoint: {}".format(checkpoint_fp.as_posix()))
        checkpoint = torch.load(checkpoint_fp.as_posix(), map_location="cpu")
        Checkpoint.load_objects(to_load=to_save, checkpoint=checkpoint)

    return trainer


def create_evaluator(model, metrics, config, tag="val"):
    with_amp = config["with_amp"]
    device = idist.device()

    @torch.no_grad()
    def evaluate_step(engine: Engine, batch):
        model.eval()
        x, y = batch[0], batch[1]
        if x.device != device:
            x = x.to(device, non_blocking=True)
            y = y.to(device, non_blocking=True)

        with autocast(enabled=with_amp):
            output = model(x)
        return output, y

    evaluator = Engine(evaluate_step)

    for name, metric in metrics.items():
        metric.attach(evaluator, name)

    return evaluator


def log_metrics(logger, epoch, elapsed, tag, metrics):
    metrics_output = "\n".join([f"\t{k}: {v}" for k, v in metrics.items()])
    logger.info(f"\nEpoch {epoch} - Evaluation time (seconds): {elapsed:.2f} - {tag} metrics:\n {metrics_output}")


def log_basic_info(logger, config):
    logger.info(f"Train 3D Unet on Cryo-EM dataset of proteins")
    logger.info(f"- PyTorch version: {torch.__version__}")
    logger.info(f"- Ignite version: {ignite.__version__}")
    if torch.cuda.is_available():
        # explicitly import cudnn as
        # torch.backends.cudnn can not be pickled with hvd spawning procs
        from torch.backends import cudnn

        logger.info(f"- GPU Device: {torch.cuda.get_device_name(idist.get_local_rank())}")
        logger.info(f"- CUDA version: {torch.version.cuda}")
        logger.info(f"- CUDNN version: {cudnn.version()}")

    logger.info("\n")
    logger.info("Configuration:")
    for key, value in config.items():
        logger.info(f"\t{key}: {value}")
    logger.info("\n")

    if idist.get_world_size() > 1:
        logger.info("\nDistributed setting:")
        logger.info(f"\tbackend: {idist.backend()}")
        logger.info(f"\tworld size: {idist.get_world_size()}")
        logger.info("\n")

def log_weights_hist(parameters, logger, status='before'):
    fig = plt.figure(figsize=(15,5))
    w_index = 0
    for i in range(1,4):
        ax = fig.add_subplot(1,3,i)
        weights_list = parameters[w_index:w_index+3]
        for k,p in enumerate(weights_list):
            current_flat = torch.flatten(p[1].clone().detach().cpu())
            ax.hist(current_flat.numpy(),histtype='step',bins='auto',label=str(p[0]))
            ax.legend()
        ax.set_title('rom Layer : '+str(w_index+1)+' to '+str(w_index+3))
        w_index +=3
        plt.title('weights_'+str(status)+"_training.png")
        logger.report_matplotlib_figure(title='From Layer : '+str(w_index+1)+' to '+str(w_index+3), series=status, iteration=0, figure=plt ) 
        plt.clf()
    plt.close('all')

def get_save_handler(config):
    return ClearMLSaver(dirname=config["output_path"])


if __name__ == "__main__":
    fire.Fire({"run": run})
