import matplotlib as mpl
import matplotlib.pyplot as plt
plt.style.reload_library()
plt.style.use(['science','ieee','std-colors'])
from sklearn.metrics import roc_curve,auc,precision_recall_curve,average_precision_score
from sklearn.preprocessing import label_binarize
from scipy import interp
import torch
import pandas as pd
import numpy as np

from SegmentationUNet import SegmentationUNet
from SegmentationDataset import SegmentationDataset
from torch.utils.data import DataLoader	



num_classes=3
folds = 3
depth=4
img_size=(64,64,64)
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
extra_width = 8
seed = 42
randg = torch.Generator()
randg.manual_seed(seed)

df = pd.read_csv('dataset/dataset_patches_final_new.csv', dtype=str)

test_list = ['4671', '5140', '6618', '8723', '22423', '8722', '4156', '4400', '20495', '13387', '20613', '11838', '10274', '24578', '25971', '3435', '8438', '25980', '3353', '6207', '4057', '11688', '22964', '3355', '30239', '22754', '21693', '22115', '25972', '23064', '8721']
# Get id list of EM maps in data
id_list = df.groupby('id')['id'].unique().index.tolist()
test_split = df[df['id'].isin(test_list)]
test_dataset = SegmentationDataset(test_split, num_classes, img_size, randg, device,augmentate=False, extra_width=extra_width, is_validation=True,is_nonoverlap_stride=True)

#models_patches =[['results/3DUnet-0x8fb3-_fold-0-1_20230606-155914/best_model_198_0.8113.pt','results/3DUnet-0xac2e-_fold-1-1_20230617-195714/best_model_198_0.5154.pt','results/3DUnet-0x454c-_fold-2-1_20230628-172120/best_model_196_0.8272.pt']]
models_patches =['results/3DUnet-0x8fb3-_fold-0-1_20230606-155914/best_model_198_0.8113.pt']

test_loader = DataLoader(test_dataset, batch_size=len(test_dataset)//2, shuffle=False)

for models in models_patches:
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    precision = dict()
    recall = dict()
    average_precision = dict()
    for i,fold in enumerate([models]):
        model = SegmentationUNet(num_classes, device, depth=depth, in_channels=2,batch_norm=True)
        model.load_state_dict(torch.load(fold, map_location=torch.device(device)))
        model.to(device)
        model.eval()
        preds = []
        targets = []
        with torch.no_grad():
            for inputs,y in test_loader:
                pred = model(inputs)
                prob = torch.nn.functional.softmax(pred, dim=1)
                y = label_binarize(y.cpu().detach().numpy().ravel(), classes=[0, 1, 2])
                preds.append(prob.cpu().detach().numpy())
                targets.append(y)
        fpr_i = dict()
        tpr_i = dict()
        roc_auc_i = dict()
        preds = np.concatenate(preds)
        targets = np.concatenate(targets)
        # Set up plot ROC
        #fig, ax = plt.subplots()
        #ax.plot([0, 1], [0, 1])
        for j in range(num_classes):
            fpr_i[j], tpr_i[j], _ = roc_curve(targets[:, j], preds[:, j].reshape(-1))
            roc_auc_i[j] = auc(fpr_i[j], tpr_i[j])
            print("AUC for fold {} class {}: {}".format(i+1, j, roc_auc_i[j]))
            #ax.plot(fpr[j], tpr[j], label="Class {} (AUC = {:.2f})".format(j,roc_auc[j]))
        fpr_i["micro"], tpr_i["micro"], _ = roc_curve(targets.ravel(), preds.ravel())
        roc_auc_i["micro"] = auc(fpr_i["micro"], tpr_i["micro"])
        print("AUC micro for fold {}: {}".format(i+1,roc_auc_i["micro"]))
        #ax.plot(fpr_i["micro"], tpr_i["micro"], label="Micro-average (AUC = {:.2f})".format(roc_auc_i["micro"]))
        #ax.set_xlim([0.0, 1.0])
        #ax.set_ylim([0.0, 1.05])
        #ax.set_xlabel("False Positive Rate")
        #ax.set_ylabel("True Positive Rate")
        #ax.legend(loc="lower right")
        #ax.autoscale(tight=True)
        #plt.savefig("roc_fold_{}_model_{}.eps".format(i+1,model_name))
        # PR plot
        precision_i = dict()
        recall_i = dict()
        average_precision_i = dict() 
        #fig,ax = plt.subplots()
        #no_skill = len(targets[targets==1]) / len(targets)
        #print("noskill {}".format(no_skill))
        #ax.plot([0, 1], [no_skill, no_skill])
        for j in range(num_classes):
            precision_i[j], recall_i[j], _ = precision_recall_curve(targets[:, j], preds[:, j].reshape(-1))
            average_precision_i[j] = average_precision_score(targets[:, j], preds[:, j].reshape(-1))
            print("AP for fold {} class {}: {}".format(i+1, j, average_precision_i[j]))
            #ax.plot(recall_i[j],precision_i[j], label='Class {} (AP = {:.2f})'.format(j,average_precision_i[j]))
        precision_i["micro"], recall_i["micro"], th = precision_recall_curve(targets.ravel(), preds.ravel())
        average_precision_i["micro"] = average_precision_score(targets, np.swapaxes(np.swapaxes(preds,0,1).reshape(3,-1),0,1), average="micro") 
        print("Micro Precision for fold {}:{}".format(i+1,precision_i['micro']))
        print("Micro recall for fold {}:{}".format(i+1,recall_i['micro']))
        print("AP micro for fold {}: {}".format(i+1, average_precision_i["micro"]))
        #ax.plot(recall_i['micro'],precision_i['micro'], label='Micro-average (AP = {:.2f})'.format(average_precision_i['micro']))
        #ax.set_xlim([0.0, 1.0])
        #ax.set_ylim([0.0, 1.05])
        #ax.set_xlabel("Recall")
        #ax.set_ylabel("Precision")
        #ax.legend(loc='lower left')
        #ax.autoscale(tight=True)
        #plt.savefig("pr_fold_{}_model_{}.eps".format(i+1,model_name))
        fpr[i] = fpr_i
        tpr[i] = tpr_i
        roc_auc[i] = roc_auc_i
        precision[i] = precision_i
        recall[i] = recall_i
        average_precision[i] = average_precision_i
    '''
    for j in range(num_classes):
        fig,ax = plt.subplots()
        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)
        for i,fold in enumerate(models):
            fpr_i = fpr[i]
            tpr_i = tpr[i]
            roc_auc_i = roc_auc[i]
            tprs.append(interp(mean_fpr, fpr_i[j], tpr_i[j]))
            tprs[-1][0] = 0.0
            aucs.append(roc_auc_i[j])
            ax.plot(fpr_i[j],tpr_i[j], label='ROC Fold {} (AUC = {:.2f})'.format(i,roc_auc_i[j]))
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        mean_auc = auc(mean_fpr, mean_tpr)
        std_auc = np.std(aucs)
        ax.plot(mean_fpr, mean_tpr, label=r'Mean ROC (AUC = {:.2f} $\pm$ {:.2f})'.format(mean_auc, std_auc))
        ax.set_xlim([0.0, 1.0])
        ax.set_ylim([0.0, 1.05])
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.legend(loc="lower right")
        ax.autoscale(tight=True)
        plt.savefig("roc_class_{}_model_{}.eps".format(j,model_name))
    '''
    class_names = ['Background', 'Outside subunit', 'Inside subunit']
    fig,ax = plt.subplots()
    for j in range(num_classes):
        prs = []
        aps = []
        mean_recall = np.linspace(0, 1, 100)
        for i,fold in enumerate([models]):
            precision_i = precision[i]
            recall_i = recall[i]
            average_precision_i = average_precision[i]
            prs.append(interp(mean_recall, precision_i[j], recall_i[j]))
            aps.append(average_precision_i[j])
            ax.plot(recall_i[j],precision_i[j], label=class_names[j])
        #mean_prs = np.mean(prs, axis=0)
        #mean_aps = np.mean(aps, axis=0)
        #std_aps = np.std(aps)
        #ax.plot(mean_prs, mean_recall, label=r'Mean PR (AP = {:.2f} $\pm$ {:.2f})'.format(mean_aps, std_aps))
    ax.set_xlim([0.0, 1.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel("Recall")
    ax.set_ylabel("Precision")
    ax.legend(loc='lower left')
    ax.autoscale(tight=True)
    plt.savefig("pr_best_patches_new.eps")

     

