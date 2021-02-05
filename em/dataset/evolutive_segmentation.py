from em import molecule
from em.dataset import metrics
from skimage.measure import regionprops
import numpy as np
import json
from json import encoder
import pickle 

from sklearn.ensemble import RandomForestClassifier

def convert(o):
    if isinstance(o, np.generic): return o.item()
    raise TypeError

clf =  pickle.load(open('classifier.pkl', 'rb'))

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

#test

mol = molecule.Molecule('models/aligned_emd_0044.mrc',5.45)
mol_gt = molecule.Molecule('models/aligned_emd_0044_gt.mrc',0.001)           
data = mol.getDataAtContour(1)
data_gt = mol_gt.getDataAtContour(1)
data[data_gt==0]=0
mol.setData(data)
pop = init_population(10, 60, [1,4], [1,4], mol)
classifier = pickle.load(open('classifier.pkl', 'rb'))
for p in pop:
    print(classifier.predict_proba(p['features'].reshape(1, -1)))

ma, pa = select_parents(pop, 5)

for m,p in zip(ma,pa):
    print(classifier.predict_proba(m['features'].reshape(1, -1)))
    print(classifier.predict_proba(p['features'].reshape(1, -1)))
