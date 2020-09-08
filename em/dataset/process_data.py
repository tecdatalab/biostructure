import pandas as pd
import os
import subprocess
from glob import glob
from xml.etree import ElementTree
import argparse
import sys
sys.path.append('..')
import molecule
import matplotlib 
import matplotlib.pyplot as plt
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from Bio.PDB import PDBParser
import numpy as np

# Excecute command
def execute_command(cmd):
    try:
        if os.system(cmd) != 0:
            raise Exception('Command "%s" does not exist' % command_emdb)
    except Exception as exc:
        raise RuntimeError('Command "%s" does not work' % cmd) from exc

# Sync headers to folder
def get_headers(header_path):
    if not os.path.exists(header_path):
        os.makedirs(header_path)
    print("Getting v1.9 headers...")
    command = 'rsync -rlpt --ignore-existing -v -z --delete --port=33444 --include "emd-*-v19.xml" "rsync.rcsb.org::emdb/structures/EMD-*/header/" '+ header_path
    execute_command(command)

# Scan headers and create a dataframe with all samples
def generate_dataframe(header_path):
    # Get list of xml headers
    all_xml_19 = glob(os.path.join(header_path,'*-v19.xml'))
    # Creates a dict with map path as key and map ID as value
    id_path_dict = {x:os.path.splitext(os.path.basename(x))[0][:-4] for x in all_xml_19}
    # Creates a dataframe to store data
    df = pd.DataFrame(all_xml_19, columns=["path"])
    # Creates map id column
    df["id"] = df["path"].map(id_path_dict.get)

    df['fitted_entries']=None
    df['method']= None
    df['resolution']=None
    df['contourLevel']=None

    #iterate over all xml files
    for i,meta in enumerate(all_xml_19):
        #parse xml content
        with open(meta, 'rt') as f:
            tree = ElementTree.parse(f)

            pdbe_fitted_entries = []
            for pdbe_entry in tree.findall('.//deposition/fittedPDBEntryIdList/fittedPDBEntryId'):
                pdbe_fitted_entries.append(pdbe_entry.text)
            
            processing_method_node = tree.find('.//processing/method')
            resolution_node = tree.find('.//processing/reconstruction/resolutionByAuthor')
            contour_node = tree.find('.//map/contourLevel')
        # If has several entries, use the first one.
        if len(pdbe_fitted_entries)>0:
            if len(pdbe_fitted_entries)>1:
                print("Model {} has more than one fitted structures, using the first one.".format(meta[-4]))
            df.loc[ df.id == id_path_dict[meta], ['fitted_entries'] ] =  pdbe_fitted_entries[0]
        # Get reported processing method
        if processing_method_node is not None:
            df.loc[ df.id == id_path_dict[meta], ['method'] ]  =  processing_method_node.text
        # Get reported resolution 
        if resolution_node is not None:
            df.loc[ df.id == id_path_dict[meta], ['resolution'] ] =  resolution_node.text
        # Get reported contour
        if contour_node is not None:
            df.loc[ df.id == id_path_dict[meta], ['contourLevel'] ] = contour_node.text 
        if i % 100:
            print("Processed {}".format(i))
    df.to_csv('dataset_metadata.csv', index=False)


# Process initial metadata, drop records that don't satisfy criteria
def processMetadata():
    df = pd.read_csv('./dataset_metadata.csv')
    # Get number of null elements to be removed
    non_fitted= len(df.index) - len(df.dropna(subset=['fitted_entries']).index)
    non_res = len(df.index) - len(df.dropna(subset=['resolution']).index)
    # Number of samples out of resolution gap
    outside_gap = len(df.index) - len(df[(df.resolution<4.5) & (df.resolution>10)].index)

    # Remove null values
    df = df.dropna(subset=['fitted_entries', 'resolution'])
    # Keep only samples within resolution gap
    df = df[(df.resolution>=4.5) & (df.resolution<=10)]

    print("{} entries without fitted structure".format(non_fitted))
    print("{} entries without reported resolution".format(non_res))
    print("{} entries out of resolution gap".format(outside_gap))
    print("{} candidates for dataset".format(len(df.index)))

    df["pdb_rsync"] = df["fitted_entries"].map(lambda pdb_id: "rsync.rcsb.org::ftp_data/structures/all/pdb/pdb"+pdb_id+".ent.gz")
    df["emdb_rsync"] = df["id"].map(lambda map_id: "rsync.rcsb.org::emdb/structures/"+str.upper(map_id)+"/map/"+map_id.replace("-","_")+".map.gz")

    df.to_csv('dataset_metadata.csv', index=False)

# Download candidate experimental models from database
def downloadModels(models_path):
    if not os.path.exists(models_path):
        os.makedirs(models_path)
    df = pd.read_csv('./dataset_metadata.csv')
    
    # Check existing models
    print("Check existing files to download new models only..")
    pdb_files = [f for f in glob(os.path.join(models_path,'*.ent'))]
    map_files = [f for f in glob(os.path.join(models_path,'*.map'))]

    pdb_filenames = [os.path.basename(f) for f in pdb_files]
    map_filenames = [os.path.basename(f) for f in map_files]

    pdb_id_path_dict = dict(zip(pdb_filenames,pdb_files))
    map_id_path_dict = dict(zip(map_filenames,map_files))

    print("{} pdbs and {} maps found".format(len(pdb_files), len(map_files)))

    #Get entries with missing pdb or map
    df['map_found'] = df["id"].map(lambda map_id: True if map_id.replace('-','_')+'.map' in map_filenames else False)
    df['pdb_found'] = df["fitted_entries"].map(lambda pdb_id: True if 'pdb'+pdb_id+'.ent' in pdb_filenames else False)

    # Choose only non existing files to download
    df_maps = df[df['map_found'] == False] 
    df_pdbs = df[df['pdb_found'] == False] 
    emdb_ftp_list = df_maps["emdb_rsync"].tolist()
    pdb_ftp_list = df_pdbs["pdb_rsync"].tolist()
    emdb_id_list = df_maps["id"].tolist()
    pdb_id_list = df_pdbs["fitted_entries"].tolist()
        
    for uri,name in zip(emdb_ftp_list, emdb_id_list):
        command_emdb = 'rsync -rlpt --ignore-existing -v -z --delete --port=33444 '+ uri  +' '+ models_path+'/'
        execute_command(command_emdb) 

    
    for uri,name in zip(pdb_ftp_list, pdb_id_list):
        command_pdb = 'rsync -rlpt --ignore-existing -v -z -L --delete --port=33444 '+  uri  +' '+ models_path+'/'
        execute_command(command_pdb)    
    
    command = 'gunzip -f ' +models_path+'/*.gz'
    execute_command(command)

# Download candidate experimental models from database
def downloadModelsPDB(pdb_path):
    if not os.path.exists(pdb_path):
        os.makedirs(pdb_path)

    # Copy pdb files from database
    print("Copying pdbs to generate simulated dataset...")
    command = 'rsync -rlpt --ignore-existing -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/pdb/'  +' '+ pdb_path
    #execute_command(command)
 
    # Unziping files 
    command = 'gunzip -f -r '+ pdb_path 
    #execute_command(command)
    
    # Move files from tree dir to pdb folder
    command = 'find '+pdb_path+' -name \'*.ent\' -exe mv {} '+pdb_path+' \;' 
    #execute_command(command)
    pdb_files = [f for f in glob(os.path.join(pdb_path,'*.ent'))]

    pdb_entries = [os.path.basename(f)[3:-4] for f in pdb_files]

    df_simulated = pd.DataFrame( list(zip(pdb_entries, pdb_files)), columns= ['entries','path'] )

    df_simulated.to_csv('metadata_synthetic.csv')

    print("{} pdbs downloaded".format(len(pdb_files)))
   

# Remove models without atom structure
def removeNonExistingModels(models_path):
    print("Scanning existing models..")
    pdb_files = [f for f in glob(os.path.join(models_path,'*.ent'))]
    map_files = [f for f in glob(os.path.join(models_path,'*.map'))]

    pdb_filenames = [os.path.basename(f) for f in pdb_files]
    map_filenames = [os.path.basename(f) for f in map_files]

    pdb_id_path_dict = dict(zip(pdb_filenames,pdb_files))
    map_id_path_dict = dict(zip(map_filenames,map_files))

    print("{} pdbs and {} maps found".format(len(pdb_files), len(map_files)))

    df = pd.read_csv('./dataset_metadata.csv')

    #Get entries with missing pdb or map
    df['map_found'] = df["id"].map(lambda map_id: True if map_id.replace('-','_')+'.map' in map_filenames else False)
    df['pdb_found'] = df["fitted_entries"].map(lambda pdb_id: True if 'pdb'+pdb_id+'.ent' in pdb_filenames else False)

    df_map_not_found = df[df["map_found"]==False]
    df_pdb_not_found = df[df["pdb_found"]==False]
    # Remove entries that doesn't have neither map or pdb
    indexes_to_remove = df[ (df['map_found'] == False) | (df['pdb_found'] == False) ].index
    #indexes_to_remove = df[ df['pdb_found'] == False ].index
    df = df.drop(indexes_to_remove)
    
    df['map_path'] = df["id"].map(lambda map_id: map_id_path_dict[map_id.replace('-','_')+'.map'])
    df['pdb_path'] = df["fitted_entries"].map(lambda pdb_id: pdb_id_path_dict['pdb'+pdb_id+'.ent'])
    
    print("{} maps not found".format(len(df_map_not_found.index)))
    print("{} pdbs not found".format(len(df_pdb_not_found.index)))

    df_map_not_found.to_csv('emdb_not_found.csv', index=False)
    df_pdb_not_found.to_csv('pdb_not_found.csv', index=False)

    df.to_csv('dataset_metadata.csv', index=False)

# Get chains in structure and compute number of subunits
def processfiles(models_path):
    
    # Define function to get number of chains in structure
    def pdb_num_chain_mapper(chain, filepath):
        parser = PDBParser(PERMISSIVE = True, QUIET = True)
        pbd_obj = parser.get_structure(chain, filepath)
        print("Pdb {} has {} chains".format(chain,len(list(pbd_obj.get_chains()))))
        return len(list(pbd_obj.get_chains()))

    print("Processing experimental samples")

    df = pd.read_csv('./dataset_metadata.csv')
    # Create dictionary from dataframe
    pdb_id_path_dict = pd.Series(df.pdb_path.values,index=df.fitted_entries).to_dict()
    # Compute number of subunits
    df['subunit_count'] = df["fitted_entries"].map(lambda pdb_id: pdb_num_chain_mapper(pdb_id, pdb_id_path_dict[pdb_id]))
    # Get number of models with less than 2 subunits
    non_enough_subunit = len(df[df.subunit_count<2].index)
    # Remove entries with less than 2 subunits
    df = df[df.subunit_count>=2]

    df.to_csv('dataset_metadata.csv', index=False)
    
    number_of_samples = len( df.index )
   
    print("{} experimental entries with less than 2 subunits".format(non_enough_subunit))
    print("{} experimental candidates for dataset".format(number_of_samples))

    print("Processing simulated samples")
    df_simulated = pd.read_csv('./metadata_synthetic.csv')
    # Create dictionary from dataframe
    pdb_id_path_dict = pd.Series(df_simulated.path.values,index=df_simulated.entries).to_dict()
    # Compute number of subunits
    df_simulated['subunit_count'] = df_simulated["entries"].map(lambda pdb_id: pdb_num_chain_mapper(pdb_id, pdb_id_path_dict[pdb_id]))
    # Get number of models with less than 2 subunits
    non_enough_subunit = len(df_simulated[df_simulated.subunit_count<2].index)
    # Remove entries with less than 2 subunits
    df_simulated = df_simulated[df_simulated.subunit_count>=2]
    number_of_samples = len( df_simulated.index )
   
    print("{} pdb entries with less than 2 subunits".format(non_enough_subunit))
    print("{} pdb candidates for simulated dataset".format(number_of_samples))

    df_simulated.to_csv('./metadata_synthetic.csv')




# Simulate map from pdb structure and calculate volume 
def simulateMapAndCompareVolume(index, df, sim_model_path):
    if not os.path.exists(sim_model_path):
        os.makedirs(sim_model_path)
    map_filename = df.at[df.index[index], 'map_path']
    pdb_filename = df.at[df.index[index], 'pdb_path']
    res = df.at[df.index[index], 'resolution']
    contourLevel = df.at[df.index[index], 'contourLevel']

    # Get map volume
    map_object = molecule.Molecule(map_filename, recommendedContour=contourLevel, cutoffRatios=[1])
    # Get dictionary of volumes, choose element with key 1 (corresponding to 100% recommended contour level)
    map_box = map_object.getGridSize()
    map_box_max_dim = np.max(map_box)
    if map_box_max_dim != np.min(map_box):
        # Change map extension to mrc due to compatibility with EMAN
        new_map_filename = map_filename.replace('.map','.mrc')
        df.at[df.index[index], 'map_path'] = new_map_filename
        command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2proc3d.py '+ map_filename + ' ' + new_map_filename +' --clip='+str(map_box_max_dim)+','+str(map_box_max_dim)+','+str(map_box_max_dim)
        print(command)
        execute_command(command)
    
    map_volume = map_object.getVolume()[1]
    # Get voxel size in Angstroms to generate simulated map 
    voxel_size = map_object.getVoxelSize()[0] #using only one value, supose its the same for all axes
    if voxel_size == 0:
        print("Map {} has voxel volume of 0, header is: \n {}".format(os.path.basename(map_filename), map_object.emMap.rawHeader))
    # Get map bounding box

    simulated_path = os.path.join(sim_model_path, 'sim_'+os.path.basename(map_filename).replace('.map','.mrc'))
    # Generate map
    command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2pdb2mrc.py -A=' + str(voxel_size)+ ' -R=' + str(res) + ' -B='+ str(map_box[2]) +','+ str(map_box[1])+ ','+ str(map_box[0]) +' '+  pdb_filename + ' ' + simulated_path
    print(command)
    execute_command(command)

    volume_result = dict()
    volume_result['index']=index
    volume_result['simulated_path'] = simulated_path
    volume_result['map_volume']=map_volume
    
    # Get simulated map volume
    if os.path.exists(simulated_path):
        map_object = molecule.Molecule(simulated_path, recommendedContour=contourLevel, cutoffRatios=[1])
        # Get dictionary of volumes, choose element with key 1 (corresponding to 100% recommended contour level)
        pdb_volume = map_object.getVolume()[1]
        volume_result['pdb_volume']=pdb_volume
    else:
        volume_result['pdb_volume']=-1

    return volume_result


# Function to generateSimulated samples from pdb
def generateSimulatedDataset(index, df, sim_model_path):
    if not os.path.exists(sim_model_path):
        os.makedirs(sim_model_path)
    pdb_filepath = df.at[df.index[index], 'path']
    pdb_entry = df.at[df.index[index], 'entries']

    # Choose random resolution over 5 and below 10 A
    upper_res = 10
    bottom_res = 5
    res = (upper_res-bottom_res) * np.random.random_sample() + bottom_res

    simulated_path = os.path.join(sim_model_path, 'sim_'+pdb_entry+'.mrc')

    result = dict()
    result['index'] = index
    result['resolution']=res
    result['path'] = simulated_path
    # Generate map
    try:
        command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2pdb2mrc.py -A=1.0 -R=' + str(res) + ' --center '+  pdb_filepath + ' ' + simulated_path
        print(command)
        execute_command(command)
    except:
        result['resolution']=-1
    return result

# Function to fit experimental map to match simulated map from pdb structure
def fitMaptoSim_Map(index, df, models_path):
    map_path = df.at[df.index[index], 'map_path']
    simulated_path = df.at[df.index[index], 'simulated_path']
    contourLvl = df.at[df.index[index], 'contourLevel']
    # Get filenames    
    map_filename = os.path.basename(map_path)
    mask_filename = 'mask_'+map_filename[:-4]+'.mrc'
    simulated_filename = os.path.basename(simulated_path)
    mask_sim_filename = 'mask_'+simulated_filename
    aligned_map_filename = 'aligned_'+map_filename[:-4]+'.mrc'
    mask_aligned_filename = 'mask_'+aligned_map_filename

    # Compute absolute paths
    aligned_map_path = os.path.join(models_path, aligned_map_filename)
    mask_sim_path = os.path.join(os.path.dirname(simulated_path), mask_sim_filename)
    mask_path = os.path.join(models_path, mask_filename)    
    mask_aligned_path = os.path.join(models_path, mask_aligned_filename)

    # Compute binary mask for simulated and original map
    command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2proc3d.py ' + map_path + ' '+ mask_path +' --process threshold.binary:value=' + str(contourLvl)
    print(command)
    execute_command(command)
    #command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2proc3d.py ' + simulated_path + ' '+ mask_sim_path +' --process threshold.binary:value=' + str(contourLvl)
    #print(command)
    #execute_command(command)
    
    # Compute similarity before alignment
    command =  ['/work/mzumbado/EMAN2/bin/python', '/work/mzumbado/EMAN2/bin/sximgstat.py', map_path, simulated_path, mask_path, '--ccc']
    print(' '.join(command))
    result = subprocess.run(command, stderr=subprocess.STDOUT, stdout=subprocess.PIPE, universal_newlines=True)
    print("Command output: {}".format(result.stdout))
    correlation_before = result.stdout.split()[-1] if result.stdout != None else 0
    
    # Compute alignment
    command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2proc3d.py  --align=rotate_translate_3d_tree:verbose=True --alignref='+simulated_path+' --multfile='+ mask_path +' '+map_path+' '+aligned_map_path
    print(command)
    try:
        execute_command(command)
    except RuntimeError as e:
        print(e)
    # Compute mask for aligned map
    command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2proc3d.py ' + aligned_map_path + ' '+ mask_aligned_path +' --process threshold.binary:value=' + str(contourLvl)
    print(command)
    try:
        execute_command(command)
    except RuntimeError as e:
        print(e)

    # Compute similarity for aligned map
    command =  ['/work/mzumbado/EMAN2/bin/python', '/work/mzumbado/EMAN2/bin/sximgstat.py', aligned_map_path, simulated_path, mask_path, '--ccc']
    print(' '.join(command))
    result = subprocess.run(command,  stderr=subprocess.STDOUT, stdout=subprocess.PIPE, universal_newlines=True)
    print("Command output: {}".format(result.stdout))
    correlation_after = result.stdout.split()[-1] if result.stdout != None else 0

    return {'index':index, 'similarity':correlation_before, 'aligned_path':aligned_map_path, 'similarity_aligned':correlation_after}
     

# Function to get symmetry information from pdb file
def getSymmetryFromPDB(models_path):
    import __main__
    __main__.pymol_argv = ['pymol','-qc']
    import pymol
    from pymol import cmd, stored
    pymol.finish_launching()

    def hasSymmetry(filename):
        cmd.reinitialize()
        cmd.load(filename)
        res = cmd.get_symmetry()
        if res is None:
            res = []
        return res
        
    list_pdb = df['pdb_path'] = df["fitted_entries"].map(lambda pdb_id:'/home/manzumbado/Development/Asistencia/BECA-TEC/biostructure/em/dataset/models/pdb'+pdb_id+'.ent').tolist()
    df['symm'] = df['pdb_path'].map(lambda x: True if (hasSymmetry(x)[0:-1] != [1.0, 1.0, 1.0, 90.0, 90.0, 90.0]) else False)

    df.to_csv('mv.csv')
    
# Function to select experimental data acording percentil distribution of volume capture ratio
def selectExperimentalDataset(result_dir, percentil_boundaries_tuple=None):
    # Open result csv from volume calculation
    df = pd.read_csv('dataset_volume.csv')
    # Replace missing volumes (produced by e2pdb2mrc failure to generate map from pdb)
    df['pdb_volume'].replace('', np.nan, inplace=True)
    df['map_volume'].replace('', np.nan, inplace=True)
    df.dropna(inplace=True)

    # Computing vol_capture which is the volume percent of the generated map respecto to the original one
    df['vol_capture'] = ((df['pdb_volume']*100)/df['map_volume'])/100
    # Cleaning infinite values (Some entries report 0 in volume, NEED TO CHECK)
    df['vol_capture'].replace(np.inf, np.nan, inplace=True)
    df.dropna(inplace=True)

    # Plot hist
    df.hist(column='vol_capture', grid=True, bins=30)
    plt.title('Ratio of simulated volume to em map volume at recommended contour')
    plt.xlabel('Bins')  
    plt.ylabel('Number of Samples')
    plt.savefig(os.path.join(result_dir,'volume_capture_hist.png'))
    plt.clf()
    # Computing basic stats
    num_samples = len(df)
    avg = df['vol_capture'].mean()
    med = df['vol_capture'].median()
    std = df['vol_capture'].std()
    print("Number of samples: {}".format(num_samples))
    print("Volume ratio avg: {}".format(avg))
    print("Volume ratio median: {}".format(med))
    print("Volume ratio std: {}".format(std))
    # Get column values to list
    values = df['vol_capture'].tolist()

    # Compute quartiles and plot box chart
    q1 = np.percentile(values, 25, interpolation = 'midpoint')
    q3 = np.percentile(values, 75, interpolation = 'midpoint')
    print("Quartile Q1: {}".format(q1))
    print("Quartile Q3: {}".format(q3))
    print("IQR: {}".format(q3-q1))

    low_boundary = q1-(q3-q1)*1.5 
    upper_boundary = q3+(q3-q1)*1.5 

    print("Lower boundary: {}".format(low_boundary))
    print("Upper boundary: {}".format(upper_boundary))

    print("Number of samples below lower boundary: {} and over upper boundary: {}".format(len( df[df.vol_capture < low_boundary] ), len( df[df.vol_capture >upper_boundary] )))
    
    # Printing percentile values
    for i in range(0,100):
        print("Percentile {} value: {}".format(i, df.vol_capture.quantile(i/100, interpolation = 'midpoint')))

    percentilA = df.vol_capture.quantile(int(percentil_boundaries_tuple[0])/100, interpolation = 'midpoint')
    percentilB = df.vol_capture.quantile(int(percentil_boundaries_tuple[1])/100,  interpolation = 'midpoint')

    if percentil_boundaries_tuple != None:
        selected_samples = df[ (df.vol_capture > percentilA) & (df.vol_capture < percentilB)]
        numMuestras = len( selected_samples)
        print("Number of samples between percentil {} and {}: {}".format(percentil_boundaries_tuple[0],percentil_boundaries_tuple[1], numMuestras))
    else: 
        selected_samples = df[(df.vol_capture < low_boundary) | (df.vol_capture >upper_boundary) ]
        print("Getting outliers following IQR...")
    # Save csv with selected samples
    selected_samples.to_csv('dataset_selected.csv')

    # Gerete boxplot
    res = df.boxplot(column='vol_capture', grid=True)
    plt.title('Ratio of simulated volume to em map volume at recommended contour')
    plt.xlabel('Volume ratio')  
    plt.ylabel('A/A')
    plt.savefig('volume_capture_barplot.png')

    


    

def main():
    global current_dir


    parser = argparse.ArgumentParser()
    parser.add_argument('--d', required=False, default=False, help='True if want to resync emdb data')
    parser.add_argument('--v', required=False, default=False, help='True if want to compare map volume with its simulated counterpart')
    parser.add_argument('--s', required=False, default=None, nargs=2, help='Please provide lower and upper experimental data percentil distribution to be included in final dataset')  
    parser.add_argument('--g', required=False, default=None, help='Please provide percent of simulated data to be included in final dataset')
    parser.add_argument('--f', required=False, default=False, help='True if want to align experimental data to simulated counterpart')
    parser.add_argument('--header_dir', default='header', required=False,help='output directory to store emdb headers') 
    parser.add_argument('--models_dir', default='models', required=False, help='output directory to store emdb models') 
    parser.add_argument('--simulated_dir', default='simulated', required=False, help='output directory to store simulated maps') 
    parser.add_argument('--result_dir', default='', required=False, help='output directory to store stats') 

    opt = parser.parse_args()
    current_dir = os.getcwd()
    header_path = os.path.join(current_dir, opt.header_dir)
    models_path = os.path.join(current_dir, opt.models_dir)
    simulated_path = os.path.join(current_dir, opt.simulated_dir)
    results_path = os.path.join(current_dir, opt.result_dir)
    pdb_path = os.path.join(models_path,'pdb')
    # Download and process data
    if int(opt.d):
        get_headers(header_path)
        generate_dataframe(header_path)
        processMetadata()
        downloadModels(models_path)
        downloadModelsPDB(pdb_path)
        removeNonExistingModels(models_path)
        processfiles(models_path)
    #Calculate volume
    if int(opt.v):
        df = pd.read_csv('./dataset_metadata.csv')
        print("Number of samples to process volume: ", len(df.index))
        df_volume = df[['id','map_path', 'fitted_entries', 'subunit_count', 'resolution','contourLevel']].copy()
        # Get index list to schedule processess 
        index_list = df_volume.index.tolist()
        print("Spawn procecess...")
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        with MPICommExecutor(comm, root=0, worker_size=size) as executor:
            if executor is not None:
                
                futures = []
                for i in index_list:
                    futures.append(executor.submit(simulateMapAndCompareVolume, i, df, simulated_path))
                for f in futures:
                    try:
                        d = f.result()
                        print("received vol dict: ",d)
                        res_index = d['index']
                        pdb_volume = d['pdb_volume']
                        map_volume = d['map_volume']
                        sim_path = d['simulated_path']
                        df_volume.loc[res_index, 'map_volume'] = map_volume
                        df_volume.loc[res_index, 'pdb_volume'] = pdb_volume
                        df_volume.loc[res_index, 'simulated_path'] = sim_path
                    except ValueError as error:
                        print("Error calculating volume for simulated {}: {}".format(df_volume.loc[i,'fitted_entries'],error))
                df_volume.to_csv('dataset_volume.csv', index=False)
    if opt.s != None:
        selectExperimentalDataset(results_path, percentil_boundaries_tuple=opt.s)
    # Generate simulated data to be included in final dataset
    if opt.g != None:
        # Open candidates to compute simulated map 
        df_synthetic = pd.read_csv('metadata_synthetic.csv')
        df_experimental = pd.read_csv('dataset_selected.csv')
        num_experimental = len(df_experimental)
        print("Number of candidate synthetic models: ", len(df_synthetic.index))
        proportion = int(opt.g)/100
        num_samples_to_pick = int(round(proportion*num_experimental/1-proportion))
        print("A proportion of {} synthetic data in final dataset requires {} simulated and {} experimental samples.".format(proportion,num_samples_to_pick,num_experimental ))
        df_synthetic_selected = df_synthetic.sample(num_samples_to_pick)
        # Get index list to schedule processess
        df_synthetic_selected.reset_index(inplace=True, drop=True)  
        index_list = df_synthetic_selected.index.tolist()
        print("Spawn procecess...")
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        with MPICommExecutor(comm, root=0, worker_size=size) as executor:
            if executor is not None:
                futures = []
                for i in index_list:
                    futures.append(executor.submit(generateSimulatedDataset, i, df_synthetic_selected, simulated_path))
                for f in futures:
                    try:
                        d = f.result()
                        print("received res dict: ",d)
                        res_index = d['index']
                        res = d['resolution']
                        path = d['path']
                        df_synthetic_selected.loc[res_index, 'resolution'] = res
                        df_synthetic_selected.loc[res_index, 'path'] = path
                    except ValueError as error:
                        print("Error calculating volume for simulated {}: {}".format(df_volume.loc[i,'fitted_entries'],error))
                df_synthetic_selected.to_csv('synthetic_selected.csv', index=False)
    
    if int(opt.f):
        # Fitting de los mapas
        print("Fitting experimental maps to match simulated counterpart")
        df_selected = pd.read_csv('dataset_selected.csv')
        index_list = df_selected.index.tolist()
        print("Spawn procecess...")
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        with MPICommExecutor(comm, root=0, worker_size=size) as executor:
            if executor is not None:
                futures = []
                for i in index_list:
                    futures.append(executor.submit(fitMaptoSim_Map, i, df_selected, models_path))
                for f in futures:
                    try:
                        d = f.result()
                        print("received res dict: ",d)
                        res_index = d['index']
                        corr = d['similarity']
                        corr_aligned = d['similarity_aligned']
                        path = d['aligned_path']
                        df_selected.loc[res_index, 'similarity'] = corr
                        df_selected.loc[res_index, 'similarity_aligned'] = corr_aligned
                        df_selected.loc[res_index, 'aligned_path'] = path
                    except ValueError as error:
                        print("Error fitting map {} to simulated counterpart: {}".format(df_volume.loc[i,'ID'],error))
                df_selected.to_csv('dataset_selected_sim.csv', index=False)
    '''
    if opt.p !=None:
        # Partir en pedazos, generar ground truth
        print("Tagging data")
    '''


       





       


       
        


if __name__ == '__main__':
    main()
