import os
import argparse
from em import molecule
from em.dataset.metrics import matching_iou
from em.evolutionary_segmentation import evolutionary_segmentation as evo_seg
from em.evolutionary_segmentation import miscellaneous as misc
import numpy as np
import pandas as pd
import pickle 
import copy


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
    parser.add_argument('--input_path', required=True, help= 'Input mrc file path')
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
    classifier = pickle.load(open('classifier.pkl', 'rb'))
    print(opt.input_path)
    mol = molecule.Molecule(opt.input_path,float(opt.level))
    #mol_gt = molecule.Molecule(opt.label_path, 0.01)
    data = mol.getDataAtContour(1)
    #data_gt = mol_gt.getDataAtContour(1)
    #data[data_gt==0]=0
    mol.setData(data)
   
    pop = evo_seg.init_population(int(opt.init_size), MAX_NUMBER_SEGMENTS, [1,4], [1,4], mol)
    print("Init population score:")
    for p in pop:
        print(classifier.predict_proba(p['features'].reshape(1, -1)))

    save=False
    counter = 0
    run= True
    pop_fitness = None
    overall_best_score = 0
    patience = int(opt.n_patience)
    test_id = os.path.basename(opt.input_path).split('.')[0]
    
    score_list = []
    score_segmentation_list = []
    
    while(run):
        ma, pa = evo_seg.select_parents(pop, int(opt.n_mates))
        new_gen = evo_seg.mating(ma, pa, float(opt.p_mates), 0.5, data)
        mutated_pop = evo_seg.mutate(new_gen, float(opt.p_split), float(opt.p_merge), [1,4], [1,4], mol)
        pop_fitness = [ classifier.predict_proba(n['features'].reshape(1, -1))[0][1] for n in mutated_pop ]
        print("Population of gen {} fitness {}".format(counter,pop_fitness))
        best_individual_score = np.max(pop_fitness)
        std_score_pop = np.std(pop_fitness)
   
        print("***Optimizing {}, best score of generation {} is {}".format(test_id, counter, best_individual_score))
        print("         population std of {}".format(std_score_pop))
        save = True if best_individual_score > overall_best_score else False
        sorted_idx = np.argsort(pop_fitness)
        overall_best_score = best_individual_score if best_individual_score > overall_best_score else overall_best_score
        score_list.append(best_individual_score)
        print('computing matching IoU..')
        current_segmentation = copy.deepcopy(mol)
        mol.setData(mutated_pop[sorted_idx[-1]]['labels'])
        score_segmentation_list.append(matching_iou(mol_gt,current_segmentation))
        if save:
            print("     saving segmentation...")
            save_path = os.path.join(output_path,'best_{0}_{1:.2f}.mrc'.format(test_id,best_individual_score))
            mol.save(save_path)
            patience = int(opt.n_patience)
        else:
            patience -= 1
        
        run = False if ((counter > int(opt.n_max)) | (patience==0)) else True
        counter += 1

    # plot results
    plot_evolutionary_graph(score_list, score_segmentation_list, os.path.join(opt.output_path,'evolutionary_result'), 'Matching IoU') 
