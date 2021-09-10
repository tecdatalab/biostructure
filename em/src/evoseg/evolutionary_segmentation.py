from em import molecule
from dataset import metrics
from skimage.measure import regionprops
import numpy as np
import json
from json import encoder
import pickle 
import copy
from sklearn.ensemble import RandomForestClassifier

def convert(o):
    if isinstance(o, np.generic): return o.item()
    raise TypeError

clf =  pickle.load(open('classifier/classifier.pkl', 'rb'))

SCALING_MIN = np.array([2.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0. \
,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-152.,-56.,-35.,-67.,-29.,-26.,-28.,-40.,-25.,-32.,-32.,-25.,-16.,-27.,-13.,-9.,-28.,-11.,-14.,-25.,-19.,-13.,-11. \
,-8.,-25.,-30.,-30.,-30.,-30.,-6.,-32.,-461.,-4.,-3.,-36.,-3.,-3.,-5.,-11.,-3.,-3.,-5.,-3.,-3.,-4.,-5.,-5.,-4.,-3.,-2.,-9.,-9.,-6.,-2.,-5.,-5.,-5.,0.,0.,0.])

SCALING_MAX = np.array([6.00000000e+01,9.99681351e-01,7.42918455e-01,9.29050051e-01 \
,5.21594140e-01,6.34306557e-01,9.29050051e-01,4.94328160e-01 \
,4.87711935e-01,4.98207008e-01,2.95542122e-01,4.35777611e-01 \
,3.01317516e-01,2.82279151e-01,4.58980084e-01,4.31133339e-01 \
,2.35659549e-01,2.22807083e-01,1.93021840e-01,1.35372520e-01 \
,3.57839678e-01,1.80313258e-01,1.72111026e-01,1.11779975e-01 \
,1.26686443e-01,1.25139822e-01,1.31950015e-01,1.80613564e-01 \
,1.85390428e-01,9.30613564e-02,9.18065153e-02,8.82728485e-02 \
,7.77530758e-01,1.10395215e-01,9.87761163e-02,8.84589572e-02 \
,1.12166133e-01,7.29273427e-02,7.05511270e-02,8.52065877e-02 \
,3.76810993e-01,8.33907366e-02,7.10520450e-02,3.76810993e-01 \
,3.26286377e-01,3.76810993e-01,5.67455715e-01,2.36081822e-01 \
,1.25500844e-01,2.25827710e-01,1.17118305e-01,8.52065877e-02 \
,1.21098693e-01,1.22311261e-01,1.36809363e-01,1.36809363e-01 \
,3.97220725e-01,7.42830030e-02,3.96936738e-02,3.94212183e-02 \
,5.42837809e-02,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00,1.00000000e+00,1.00000000e+00,1.00000000e+00 \
,1.00000000e+00])

MAX_NUMBER_SEGMENTS = 60

def create_vec_features(dictionary, max_subunits):
    matched_subunits_num = dictionary['matched_subunits']
    segments_volume = dictionary['voxels_assigned']
    euler_segments = dictionary['euler_segments']
    segments_volume = eval(segments_volume)
    segment_numbers = [ int(key) for key in segments_volume.keys()]
    segments_volume = [ segments_volume[key] for key in segments_volume.keys()]
    volume = np.sum(segments_volume)
    weigths = [ segment_volume/volume for segment_volume in segments_volume]
    hist = np.histogram(segment_numbers, weights=weigths, bins=max_subunits, range=(1,max_subunits+1))[0]
    euler_segments = eval(euler_segments)
    euler_seg_list = [ euler_segments[key] for key in euler_segments.keys()]
    ft_vector = np.zeros((2*max_subunits+1)) 
    ft_vector[0] = matched_subunits_num
    ft_vector[1:max_subunits+1] = hist
    euler_vect = np.zeros(max_subunits)
    euler_vect[:len(euler_seg_list)] = euler_seg_list
    ft_vector[max_subunits+1:] = euler_vect
    ft_vector = (ft_vector - SCALING_MIN)/ (SCALING_MAX-SCALING_MIN)
    return ft_vector
    


def init_population(pop_size, max_number_segments, steps_range, sigma_range, input_map):
    population = []
    for n in range(pop_size):
        number_segments = max_number_segments+1
        while((number_segments>max_number_segments) | (number_segments==1)):
            step = np.random.choice(range(steps_range[0],steps_range[1]+1),1)
            sigma = np.random.uniform(sigma_range[0], sigma_range[1],1)
            input_map.generateSegments(step[0], sigma[0])
            labels = input_map.labels.astype(np.int32)
            label_props = regionprops(labels)
            number_segments = len(label_props)
        segment_voxels_dict = {}
        segment_euler_dict = {}
        for l in label_props:
            segment_voxels_dict[l.label] = np.sum(labels == l.label)
            segment_euler_dict[l.label] = l.euler_number
        dict_attributes = {'matched_subunits':number_segments, 'voxels_assigned':json.dumps(segment_voxels_dict,default=convert), 'euler_segments':json.dumps(segment_euler_dict, default=convert)}
        ft_vector = create_vec_features(dict_attributes, max_number_segments)
        dict_to_append = {'features':ft_vector, 'labels':labels}
        population.append(dict_to_append) 
    return population


def fitness_function(individual):
    return clf.predict_proba(individual['features'].reshape(1, -1))[0][1]    


def select_parents(population, n_matings):
    ma = []
    pa = []
    sum_fitness = np.sum( [ fitness_function(individual) for individual in population] )
    for n in range(0, n_matings):
        r = np.random.uniform(0, sum_fitness, 1)
        s = 0
        for i in population:
            s += fitness_function(i)
            if s >= r:
                ma.append(i)
                break
    for n in range(0, n_matings):
        r = np.random.uniform(0, sum_fitness, 1)
        s = 0
        for i in population:     
            s += fitness_function(i)
            if s >= r:
                pa.append(i)
                break
    return ma, pa


def mating(ma, pa, crossover_proba, parent_proba, input_data):
    new_generation = []
    for m,p in zip(ma,pa):
        offspring = copy.deepcopy(input_data)
        marker = -10000
        offspring[offspring!=0] = marker
        #print("total voxels {}".format(np.sum(offspring!=0)))
        parents = [m['labels'],p['labels']]
        prob = np.random.uniform(0,1,1)
        if prob <= crossover_proba:
            list_indices = np.where(offspring == marker)
            list_indices = list(zip(list_indices[0], list_indices[1], list_indices[2]))
            number_voxels = np.sum(offspring==marker)
            
            current_parent = 0
            current_label = 1
            while(number_voxels!=0):
                random_voxel_ix = np.random.choice(range(len(list_indices)), size=1)[0]
                label = parents[current_parent][list_indices[random_voxel_ix]]
                for l in regionprops(parents[current_parent]):
                    if l.label==label:
                        #print("label {} found in parent {} with volume {}".format(l.label,current_parent,np.sum(parents[current_parent]==l.label)))
                        mask = parents[current_parent]==l.label
                        offspring[mask] = current_label
                        number_voxels = np.sum(offspring==marker)
                        current_parent = 1 if np.random.uniform(0,1,1)>0.5 else 0
                        #print("Tagged {} voxelx with label {}".format(np.sum(offspring==current_label), current_label))
                        #print("Voxels left: {}".format(number_voxels))
                        #print("Labels in offspring:{}".format(np.unique(offspring)))
                        current_label+=1
                        break
            if len(regionprops(offspring.astype(np.int32))) > 60:
                offspring = mating(m, p, 1, 0.5, 0.5)[0]
        elif np.random.uniform(0,1,1) <= parent_proba:
            offspring = parents[0]
        else:
            offspring = parents[1]
        offspring_voxels_dict = {}
        offspring_euler_dict = {}
        offspring_labelprops = regionprops(offspring.astype(np.int32))
        for i,l in enumerate(offspring_labelprops):
            offspring[offspring==l.label] = i+1
        offspring_labelprops = regionprops(offspring.astype(np.int32))
        for l in offspring_labelprops:
            #print("label {} volume {}".format(l.label, np.sum(offspring==l.label)))
            offspring_voxels_dict[l.label] = np.sum(offspring == l.label)
            offspring_euler_dict[l.label] = l.euler_number
        dict_attributes = {'matched_subunits':len(offspring_labelprops), 'voxels_assigned':json.dumps(offspring_voxels_dict,default=convert), 'euler_segments':json.dumps(offspring_euler_dict, default=convert)}
        offspring_ft_vector = create_vec_features(dict_attributes, MAX_NUMBER_SEGMENTS)
        new_generation.append({'features':offspring_ft_vector, 'labels':offspring.astype(np.int32)})
    return new_generation
        
        
def mutate(population, merge_prob, split_prob, steps_range, sigma_range, input_map):
    mutated_population = []
    for i,individual in enumerate(population):
        if np.random.uniform(0,1,1) <= split_prob:
            indiv_labels = individual['labels'].astype(np.int32)
            #print("preparing splitting")
            #for l in regionprops(indiv_labels):
            #    print("Label {} with {} voxels".format(l.label, np.sum(indiv_labels[indiv_labels==l.label])))
            segmented_count = 1
            while(segmented_count == 1):
                step = np.random.choice(range(steps_range[0],steps_range[1]+1),1, p=[0.4,0.4,0.1,0.1])
                sigma = np.random.uniform(sigma_range[0], sigma_range[1],1)
                input_map.generateSegments(step[0], sigma[0])
                labels = input_map.labels.astype(np.int32)
                label_props = regionprops(labels)
                segmented_count = len(label_props)
            indiv_label_props = regionprops(indiv_labels)
            indiv_label_list = [ l.label for l in indiv_label_props ]
            label_can_be_splitted = False
            count = 0
            while(label_can_be_splitted==False):
                label_to_be_splitted = np.random.choice(indiv_label_list)
                label_mask = (indiv_labels == label_to_be_splitted)
                labels_found = np.unique(labels[label_mask])
                number_segments = len(labels_found)
                if ((number_segments > 1) & (number_segments+len(indiv_label_list)-1<=60)):
                    #print("label {} can be splitted in {} segments for individual {}".format(label_to_be_splitted,number_segments,i))
                    label_can_be_splitted = True
                if count > len(indiv_label_list):
                    step = np.random.choice(range(steps_range[0],steps_range[1]+1),1, p=[0.4,0.4,0.1,0.1])
                    sigma = np.random.uniform(sigma_range[0], sigma_range[1],1)
                    #print("Recomputing segment mutation on individual {} with {} steps and {} sigma".format(i, step[0], sigma[0]))
                    input_map.generateSegments(step[0], sigma[0])
                    labels = input_map.labels.astype(np.int32)
                    count = 0
                count += 1
            #print("spliting label {} for individual {} in {} segments with labels {}".format(label_to_be_splitted,i,number_segments, labels_found))
            np.random.shuffle(labels_found)
            rename_label_dict = {}
            count = len(indiv_label_list)
            for l in labels_found:
                if count==len(indiv_label_list):
                    rename_label_dict[l]=label_to_be_splitted
                    count+=1
                else:
                    rename_label_dict[l] = count
                    count+=1
            #print("Rename label dict {}".format(rename_label_dict))
            
            new_labels = copy.deepcopy(indiv_labels)
            # Split and assign
            for key in np.sort(list(rename_label_dict.keys())):
                mask = np.logical_and(labels==key, new_labels==label_to_be_splitted)
                #print("Assigning label {} to {} voxels from gt".format(rename_label_dict[key],np.sum(mask)))
                new_labels[mask] = rename_label_dict[key]
            mutated_voxels_dict = {}
            mutated_euler_dict = {}
            mutated_labelprops = regionprops(new_labels.astype(np.int32))
            for i,l in enumerate(mutated_labelprops):
                new_labels[new_labels==l.label] = i+1
            mutated_labelprops = regionprops(new_labels.astype(np.int32))
            for l in mutated_labelprops:
                #print("label {} volume {}".format(l.label, np.sum(new_labels==l.label)))
                mutated_voxels_dict[l.label] = np.sum(new_labels == l.label)
                mutated_euler_dict[l.label] = l.euler_number
            dict_attributes = {'matched_subunits':len(mutated_labelprops), 'voxels_assigned':json.dumps(mutated_voxels_dict,default=convert), 'euler_segments':json.dumps(mutated_euler_dict, default=convert)}
            #print(dict_attributes)
            mutated_ft_vector = create_vec_features(dict_attributes, MAX_NUMBER_SEGMENTS)
            mutated_population.append({'features':mutated_ft_vector, 'labels':new_labels.astype(np.int32)})
        if np.random.uniform(0,1,1) <=  merge_prob:
            indiv_labels = individual['labels'].astype(np.int32)
            indiv_label_props = regionprops(indiv_labels)
            #print("preparing merging")
            #for l in regionprops(indiv_labels):
            #    print("Label {} with {} voxels".format(l.label, np.sum(indiv_labels[indiv_labels==l.label])))
            indiv_label_list = [ l.label for l in indiv_label_props ]
            label_to_be_merged = np.random.choice(indiv_label_list,1)[0]
            dict_label_centroids = {l.label:list(l.centroid) for l in indiv_label_props}
            distances_to_centroid_dict = {key:np.linalg.norm(np.array(dict_label_centroids[label_to_be_merged])-np.array(l_centroid), ord=2) for l_centroid,key in zip(dict_label_centroids.values(),dict_label_centroids.keys())}
            #print("Merging label {} with closest label in {}".format(label_to_be_merged, distances_to_centroid_dict))
            ordered_keys = [k for k,v in sorted(distances_to_centroid_dict.items(), key=lambda item: item[1])]
            label_to_merge_with = ordered_keys[1]
            #print("Selected label {}".format(label_to_merge_with))
            new_labels = copy.deepcopy(indiv_labels)
            new_label_id = len(indiv_label_list)+1
            new_labels[new_labels==label_to_be_merged] = new_label_id
            new_labels[new_labels==label_to_merge_with]  = new_label_id
            # Renumbering of labels
            for i,l in enumerate(regionprops(new_labels)):
                new_labels[new_labels==l.label] = i+1
            mutated_labelprops = regionprops(new_labels.astype(np.int32))
            mutated_voxels_dict = {}
            mutated_euler_dict = {}
            for l in mutated_labelprops:
                #print("label {} volume {}".format(l.label, np.sum(new_labels==l.label)))
                mutated_voxels_dict[l.label] = np.sum(new_labels == l.label)
                mutated_euler_dict[l.label] = l.euler_number
            dict_attributes = {'matched_subunits':len(mutated_labelprops), 'voxels_assigned':json.dumps(mutated_voxels_dict,default=convert), 'euler_segments':json.dumps(mutated_euler_dict, default=convert)}
            #print(dict_attributes)
            mutated_ft_vector = create_vec_features(dict_attributes, MAX_NUMBER_SEGMENTS)
            mutated_population.append({'features':mutated_ft_vector, 'labels':new_labels.astype(np.int32)})
        else: 
            mutated_population.append(individual)
    return mutated_population



#def reassign_labels(segmented_labels, gt_labels):
    





