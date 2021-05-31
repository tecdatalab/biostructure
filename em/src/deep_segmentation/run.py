from pathlib import Path

import numpy as np
import pandas as pd
import torch
from ignite.engine import Events, create_supervised_evaluator, \
    create_supervised_trainer
from ignite.metrics import DiceCoefficient, IoU, Accuracy, Loss, ConfusionMatrix
from ignite.handlers import Checkpoint, DiskSaver, global_step_from_engine
from matplotlib import pyplot as plt
plt.style.use(['science','no-latex'])
from torchsummary import summary

from SegmentationAgent import SegmentationAgent

VAL_PERCENTAGE = 0.2  # Amount of data to use for validation
TEST_NUM = 10  # Number of images to set aside for testing and visualization
NUM_CLASSES = 2  # Total number of classes in the dataset
BATCH_SIZE = 1  # Batch size for training
DIM_SIZE = 224  # The input size for model
CSV_PATH = Path('dataset/dataset.csv')  # Location of the dataset
SHUFFLE = True  # Shuffle the dataset before making the split
LR = 0.001  # Learning rate for the model
EPOCHS = 100  # Number of epochs to train the model
DEPTH = 4  # Depth of Unet model

DEVICE = 'cuda' if torch.cuda.is_available() else 'cpu'  # Device used to train

torch.multiprocessing.set_start_method('spawn')
cm = ConfusionMatrix(num_classes=NUM_CLASSES)
agent = SegmentationAgent(VAL_PERCENTAGE, TEST_NUM, NUM_CLASSES,
                          BATCH_SIZE, (DIM_SIZE,DIM_SIZE,DIM_SIZE), CSV_PATH, SHUFFLE, LR, DEPTH, DEVICE)

trainer = create_supervised_trainer(agent.model, agent.optimizer,
                                    agent.criterion)
evaluator = create_supervised_evaluator(agent.model, metrics={
        'dice': DiceCoefficient(cm, ignore_index=0), 'iou': IoU(cm, ignore_index=0), 'accuracy': Accuracy(), 'loss': Loss(agent.criterion)
})

def score_function(engine):
    return engine.state.metrics['dice'].item()

to_save = {'model': agent.model}
handler = Checkpoint(to_save, DiskSaver('./model', create_dir=True), n_saved=2,
                     filename_prefix='best', score_function=score_function, score_name="val_dice",
                     global_step_transform=global_step_from_engine(evaluator))

evaluator.add_event_handler(Events.COMPLETED, handler)

train_metrics = {'accuracy':[], 'dice':[], 'iou':[], 'loss':[]}
validation_metrics = {'accuracy':[], 'dice':[], 'iou':[], 'loss':[]}


@trainer.on(Events.EPOCH_COMPLETED)
def log_training_results(engine):
    """
    Print training accuracy and loss after each epoch
    """
    evaluator.run(agent.train_loader)
    metrics = evaluator.state.metrics
    avg_dice = metrics['dice'].item()
    avg_accuracy = metrics['accuracy']
    avg_iou = metrics['iou'].item()
    avg_loss = metrics['loss']
    global train_metrics 
    train_metrics['accuracy'].append(avg_accuracy)
    train_metrics['dice'].append(avg_dice)
    train_metrics['iou'].append(avg_iou)
    train_metrics['loss'].append(avg_loss)
    print(
            "Training - Epoch: {} Dice: {:.2f} Iou: {:.2f} Accuracy: {:.2f} Loss: {:.2f}".format(
                    engine.state.epoch, avg_dice, avg_iou, avg_accuracy, avg_loss))


@trainer.on(Events.EPOCH_COMPLETED)
def log_validation_results(engine):
    """
    Print validation accuracy and loss after each epoch
    """
    evaluator.run(agent.validation_loader)
    metrics = evaluator.state.metrics
    avg_dice = metrics['dice'].item()
    avg_accuracy = metrics['accuracy']
    avg_iou = metrics['iou'].item()
    avg_loss = metrics['loss']
    global validation_metrics
    validation_metrics['accuracy'].append(avg_accuracy)
    validation_metrics['dice'].append(avg_dice)
    validation_metrics['iou'].append(avg_iou)
    validation_metrics['loss'].append(avg_loss)
    print(
            "Validation - Epoch: {} Dice: {:.2f}  Iou: {:.2f}  Accuracy: {:.2f}  Loss: {:.2f}".format(
                    engine.state.epoch, avg_dice, avg_iou, avg_accuracy, avg_loss))

# Print summary of the model
summary(agent.model, (1, DIM_SIZE, DIM_SIZE, DIM_SIZE))
# Train model
trainer.run(agent.train_loader, max_epochs=EPOCHS)

pd.DataFrame(train_metrics).T.reset_index().to_csv('training.csv', header=False, index=False)

pd.DataFrame(validation_metrics).T.reset_index().to_csv('validation.csv', header=False, index=False)

plt.figure()
plt.plot(np.arange(len(train_metrics['dice'])) + 1, train_metrics['dice'], label='Train')
plt.plot(np.arange(len(validation_metrics['dice'])) + 1, validation_metrics['dice'],
         label='Validation')
plt.xlabel('Epochs')
plt.ylabel('Dice Coefficient')
plt.savefig('dice.svg')
plt.clf()

plt.figure()
plt.plot(np.arange(len(train_metrics['iou'])) + 1, train_metrics['iou'], label='Train')
plt.plot(np.arange(len(validation_metrics['iou'])) + 1, validation_metrics['iou'],
         label='Validation')
plt.xlabel('Epochs')
plt.ylabel('IoU')
plt.savefig('iou.svg')
plt.clf()

plt.figure()
plt.plot(np.arange(len(train_metrics['accuracy'])) + 1, train_metrics['accuracy'], label='Train')
plt.plot(np.arange(len(validation_metrics['accuracy'])) + 1, validation_metrics['accuracy'],
         label='Validation')
plt.xlabel('Epochs')
plt.ylabel('Accuracy')
plt.savefig('accuracy.svg')
plt.clf()


plt.figure()
plt.plot(np.arange(len(train_metrics['loss'])) + 1, train_metrics['loss'], label='Train')
plt.plot(np.arange(len(validation_metrics['loss'])) + 1, validation_metrics['loss'],
         label='Validation')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.savefig('loss.svg')


# Test the model on the test images
images, masks = next(iter(agent.test_loader))
with torch.no_grad():
    preds = agent.model(images)
cm = ConfusionMatrix(num_classes=NUM_CLASSES)
acc = DiceCoefficient(cm, ignore_index=0)
cm.update([preds, masks])
print('Test Dice: {}'.format(acc.compute().item()))

loss = agent.criterion(preds, masks)
loss = loss.cpu().detach().numpy()
print('Test Loss: {}'.format(loss))

images = images.cpu().numpy()
masks = masks.cpu().numpy()
preds = torch.argmax(preds, dim=1).cpu().numpy()

plt.figure(figsize=(30, 10))
for i in range(len(images)):
    plt.subplot(1, len(images), i + 1)
    #image = np.moveaxis(images[i], -1, 0)
    image = images[i]
    np.save('sample_{}.npy'.format(i), image)
    image = image * 255
    image = image.astype(int)
    plt.imshow(image[0][65])
    plt.axis('off')
plt.savefig('samples.png') 
plt.figure(figsize=(30, 10))
for i in range(len(masks)):
    plt.subplot(1, len(masks), i + 1)
    #image = np.moveaxis(masks[i], -1,0)
    image =masks[i]
    plt.imshow(image[65], cmap='Paired')
    plt.axis('off')
    np.save('mask_{}.npy'.format(i), image)
plt.savefig('masks.png')
plt.figure(figsize=(30, 10))
for i in range(len(preds)):
    plt.subplot(1, len(preds), i + 1)
    #image = np.moveaxis(preds[i], -1,0)
    image = preds[i]
    plt.imshow(image[65], cmap='Paired')
    plt.axis('off')
    np.save('pred_{}.npy'.format(i), image)
plt.savefig('preds.png')
