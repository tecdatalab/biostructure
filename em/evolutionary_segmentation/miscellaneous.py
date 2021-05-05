import matplotlib
from matplotlib import pyplot as plt
plt.style.reload_library()
plt.style.use(['science','ieee'])
plt.rcParams.update({
    "font.size":12})     



def plot_evolutionary_graph(list_scores, list_segmentation_scores, plot_name, score_label): 
    fig, axs = plt.subplots(2)
    axs[0].plot(range(len(list_scores)),list_scores, label="Best score")
    axs[0].set_ylabel('Score')
    axs[1].plot(range(len(list_segmentation_scores)),list_segmentation_scores, label=score_name)
    axs[1].set_xlabel('Generations')
    axs[1].set_ylabel('Matching IoU')
    plt.savefig(plot_name+'.eps') 
    
