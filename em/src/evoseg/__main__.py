import os
import argparse
from em import molecule
from dataset.metrics import matching_iou
from evoseg import evolutionary_segmentation as evo_seg
from evoseg import miscellaneous as misc
import numpy as np
import pandas as pd
import pickle 
import copy
import glob

MAX_NUMBER_SEGMENTS=60


if __name__ == '__main__':

    ''' ARGUMENTS
    #   --input_path imput mrc file path
    #   --label_path input index mrc file path
    #   --output_dir output directory path to save output
    #   --initial_size initial population size
    #   --n_mates number of combinations
    #   --p_mates probability of crossovers
    #   --p_split probability of split mutation
    #   --p_merge probability of merge mutation
    #   --patience number of iterations without succesive improvements to stop
    #   --top_n save top n trough generations 
    #   --n_max max number of generations
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', required=True, help= 'Input mrc directory path')
    parser.add_argument('--output_path', required=True, help= 'Output directory path')
    parser.add_argument('--level', required=True, help= 'Isosurface level')
    #parser.add_argument('--label_path', required=True, help='Input mrc file path with labels')
    parser.add_argument('--init_size', required=False, default=30, help= 'Initial population size')
    parser.add_argument('--n_mates', required=False, default=25, help= 'Number of possible combinations')
    parser.add_argument('--p_mates', required=False, default=0.5, help= 'Probability of crossover ocurrence')
    parser.add_argument('--p_split', required=False, default=0.01, help= 'Probability of split mutation')
    parser.add_argument('--p_merge', required=False, default=0.01, help= 'Probability of merge mutation')
    parser.add_argument('--n_patience', required=False, default=30, help= 'Number of iterations without succesive improvements to stop')
    parser.add_argument('--n_max', required=False, default=200, help= 'Max number of generations')

    opt = parser.parse_args()

    if not os.path.exists(opt.output_path):
        os.mkdir(opt.output_path)
    print('init')
    classifier = pickle.load(open('classifier/classifier.pkl', 'rb'))
    print("Loading maps from {}".format(opt.input_path))
    filenames = glob.glob(opt.input_path)
    for f in filenames:
        mol = molecule.Molecule(f,0.01)
        #mol_gt = molecule.Molecule(opt.label_path, 0.01)
        data = mol.getDataAtContour(1) 
        #data_gt = mol_gt.getDataAtContour(1)
        #data[data_gt==0]=0
        mol.setData(data)
       
        pop = evo_seg.init_population(30, MAX_NUMBER_SEGMENTS, [1,4], [1,4], mol)
        save=False
        counter = 0
        run= True
        pop_fitness = None
        overall_best_score = 0
        patience = 30 #int(opt.n_patience)
        test_id = os.path.basename(opt.input_path).split('.')[0]
        
        score_list = []
        score_segmentation_list = []
        top_5 = [] 
        while(run):
            ma, pa = evo_seg.select_parents(pop, 25) #int(opt.n_mates))
            new_gen = evo_seg.mating(ma, pa, 0.5, 0.5, data)
            mutated_pop = evo_seg.mutate(new_gen, 0.01, 0.01, [1,4], [1,4], mol)
            pop_fitness = [ classifier.predict_proba(n['features'].reshape(1, -1))[0][1] for n in new_gen ]
            sorted_idx = np.argsort(pop_fitness)
            if len(top_5) > 0:
                for t in top_5:
                    mutated_pop.append(t)
            top_5 = [ copy.deepcopy(new_gen[sorted_idx[i]]) for i in range(-5,0)  ]
            pop_fitness = [ classifier.predict_proba(n['features'].reshape(1, -1))[0][1] for n in mutated_pop ]
            sorted_idx = np.argsort(pop_fitness)
            print("Population of gen {} fitness {}".format(counter,pop_fitness))
            best_individual_score = np.max(pop_fitness)
            print("***Optimizing {}, best score of generation {} is {}".format(test_id, counter, best_individual_score))
            save = True if best_individual_score > overall_best_score else False
            overall_best_score = best_individual_score if best_individual_score > overall_best_score else overall_best_score
            score_list.append(best_individual_score)
            current_segmentation = mutated_pop[sorted_idx[-1]]['labels']
            if save:
                print("     saving segmentation...")
                save_path = os.path.join(opt.output_path,'best_{0}_{1:.2f}.npy'.format(test_id,best_individual_score))
                np.save(save_path, current_segmentation)
                patience = 30
            else:
                patience -= 1
            pop = mutated_pop            
            run = False if ((counter >= 200) | (patience==0)) else True
            counter += 1

        # plot results
        #plot_evolutionary_graph(score_list, score_segmentation_list, os.path.join(opt.output_path,'evolutionary_result'), 'Matching IoU') 
