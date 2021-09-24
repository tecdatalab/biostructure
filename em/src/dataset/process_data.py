import pandas as pd
import os
from glob import glob
from xml.etree import ElementTree
import argparse
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from concurrent.futures import wait
from Bio.PDB import PDBParser, PDBIO, Select
import numpy as np

from miscellaneous import *
import em.molecule as molecule

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
def processMetadata(output_txt_path):
    df = pd.read_csv('./dataset_metadata.csv')
    total_num = len(df.index) 
    # Get number of null elements to be removed
    non_fitted= len(df.index) - len(df.dropna(subset=['fitted_entries']).index)
    non_res = len(df.index) - len(df.dropna(subset=['resolution']).index)
    # Number of samples out of resolution gap
    outside_gap = len(df.index) - len(df[(df.resolution<4.5) & (df.resolution>10)].index)

    # Remove null values
    df = df.dropna(subset=['fitted_entries', 'resolution'])
    # Keep only samples within resolution gap
    df = df[(df.resolution>=4.5) & (df.resolution<=10)]

    # Print stats to output file
    with open(output_txt_path, 'a') as out:
        out.write("{} total entries in EMDB\n".format(total_num))
        out.write("{} entries without fitted structure\n".format(non_fitted))
        out.write("{} entries without reported resolution\n".format(non_res))
        out.write("{} entries out of resolution gap\n".format(outside_gap))
        out.write("This leads to {} candidates for dataset\n".format(len(df.index)))

    df["pdb_rsync"] = df["fitted_entries"].map(lambda pdb_id: "rsync.rcsb.org::ftp_data/structures/all/pdb/pdb"+pdb_id+".ent.gz")
    df["emdb_rsync"] = df["id"].map(lambda map_id: "rsync.rcsb.org::emdb/structures/"+str.upper(map_id)+"/map/"+map_id.replace("-","_")+".map.gz")

    df.to_csv('dataset_metadata.csv', index=False)

# Download candidate experimental models from database
def downloadModels(models_path, output_txt_path):
    if not os.path.exists(models_path):
        os.makedirs(models_path)
    df = pd.read_csv('./dataset_metadata.csv')
    
    # Check existing models
    with open(output_txt_path, 'a') as out:
        out.write("Check existing files to download new models only..\n")
    pdb_files = [f for f in glob(os.path.join(models_path,'*.ent'))]
    map_files = [f for f in glob(os.path.join(models_path,'*.map'))]

    pdb_filenames = [os.path.basename(f) for f in pdb_files]
    map_filenames = [os.path.basename(f) for f in map_files]

    with open(output_txt_path, 'a') as out:
        out.write("{} pdbs and {} maps found\n".format(len(pdb_files), len(map_files)))

    #Get entries with missing pdb or map
    df['map_found'] = df["id"].map(lambda map_id: True if map_id.replace('-','_')+'.map' in map_filenames else False)
    df['pdb_found'] = df["fitted_entries"].map(lambda pdb_id: True if 'pdb'+pdb_id+'.ent' in pdb_filenames else False)

    # Choose only non existing files to download
    df_maps = df[df['map_found'] == False] 
    df_pdbs = df[df['pdb_found'] == False] 
    with open(output_txt_path, 'a') as out:
        out.write("Downloading {} pdbs and {} maps \n".format(len(df_pdbs), len(df_maps)))
    emdb_ftp_list = df_maps["emdb_rsync"].tolist()
    pdb_ftp_list = df_pdbs["pdb_rsync"].tolist()
    emdb_id_list = df_maps["id"].tolist()
    pdb_id_list = df_pdbs["fitted_entries"].tolist()
     
    for uri,name in zip(emdb_ftp_list, emdb_id_list):
        command = 'rsync -rlpt --ignore-existing -v -z --delete --port=33444 '+ uri  +' '+ models_path+'/'
        try:
            execute_command(command) 
        except Exception as e:
            with open('error.txt', 'a') as error_file:
                error_file.write("Error excecuting command {}\n".format(command))
    
    for uri,name in zip(pdb_ftp_list, pdb_id_list):
        command = 'rsync -rlpt --ignore-existing -v -z -L --delete --port=33444 '+  uri  +' '+ models_path+'/'
        try:
            execute_command(command)    
        except Exception as e:
            with open('error.txt', 'a') as error_file:
                error_file.write("Error excecuting command {}\n".format(command))
    
    print("Extracting donwloaded models from emdb")
    command = 'find '+models_path+' -maxdepth 1 -name \'*.gz\' -exec gunzip -f -n {}  \;'
    try:
        execute_command(command)
    except Exception as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Error excecuting command {}\n".format(command))
    


# Download candidate experimental models from database
def downloadModelsPDB(pdb_path, models_path, out_txt_path):
    if not os.path.exists(pdb_path):
        os.makedirs(pdb_path)

    # Copy pdb files from database
    print("Copying pdbs to generate simulated dataset...")
    command = 'rsync -rlpt --ignore-existing -v -z --delete --port=33444 rsync.rcsb.org::ftp_data/structures/divided/pdb/'  +' '+ pdb_path
    execute_command(command)
    print("Copying downloaded pdbs to models dir")
    # Copying downloaded files to models directory
    command = 'find '+pdb_path+' -name \'*.gz\' -exec cp -n {} '+models_path+' \;'
    execute_command(command)
    print("Unziping copyed models")
    # Unziping files 
    command = 'find '+models_path+' -maxdepth 1 -name \'*.ent.gz\' -exec gunzip -f -n {}  \;' 
    execute_command(command)
    
    pdb_files = [f for f in glob(os.path.join(models_path,'*.ent'))]

    pdb_entries = [os.path.basename(f)[3:-4] for f in pdb_files]

    df_simulated = pd.DataFrame( list(zip(pdb_entries, pdb_files)), columns= ['entries','path'] )

    df_simulated.to_csv('metadata_synthetic.csv')
    with open(out_txt_path, 'a') as out:
        out.write("{} pdbs downloaded\n".format(len(pdb_files)))
   

# Remove models without atom structure
def removeNonExistingModels(models_path, out_txt_path):

    with open(out_txt_path, 'a') as out:
        out.write("Detecting non downloaded models..\n")
    pdb_files = [f for f in glob(os.path.join(models_path,'*.ent'))]
    map_files = [f for f in glob(os.path.join(models_path,'*.map'))]

    pdb_filenames = [os.path.basename(f) for f in pdb_files]
    map_filenames = [os.path.basename(f) for f in map_files]

    pdb_id_path_dict = dict(zip(pdb_filenames,pdb_files))
    map_id_path_dict = dict(zip(map_filenames,map_files))
    
    with open(out_txt_path, 'a') as out:
        out.write("{} pdbs and {} maps found\n".format(len(pdb_files), len(map_files)))

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
    
    with open(out_txt_path, 'a') as out:
        out.write("{} maps are missing\n".format(len(df_map_not_found.index)))
        out.write("{} pdbs are missing\n".format(len(df_pdb_not_found.index)))

    df_map_not_found.to_csv('emdb_not_found.csv', index=False)
    df_pdb_not_found.to_csv('pdb_not_found.csv', index=False)

    df.to_csv('dataset_metadata.csv', index=False)

# Define function to split PDB file and save
def pdb_split_and_save(chain, filepath):
    parser = PDBParser(PERMISSIVE = True, QUIET = True)
    pdb_obj = parser.get_structure(chain, filepath)
    io = PDBIO()
    for chain in pdb_obj.get_chains():
        model = pdb_obj[0] 
        pdb_filename_chain = filepath.replace('.ent','_'+chain._id+'.ent')    
        io.set_structure(model[chain._id])
        io.save(pdb_filename_chain)
    print("{} splitted and saved".format(chain)) 

# Define function to get number of chains in structure
def pdb_num_chain_mapper(chain, filepath):
    parser = PDBParser(PERMISSIVE = True, QUIET = True)
    pbd_obj = parser.get_structure(chain, filepath)
    print("Pdb {} has {} chains".format(chain,len(list(pbd_obj.get_chains()))))
    return len(list(pbd_obj.get_chains()))

# Define function to get number of chains in structure
def pdb_get_chain_id_label_list(chain, filepath):
    parser = PDBParser(PERMISSIVE = True, QUIET = True)
    pbd_obj = parser.get_structure(chain, filepath)
    chain_id_list = [chain._id for chain in pbd_obj.get_chains()]
    chain_label_list = [i for i in range(1,len(chain_id_list)+1)]
    return (chain_id_list,chain_label_list)

# Get chains in structure and compute number of subunits
def processfiles(models_path, out_txt_path):
    
    with open(out_txt_path, 'a') as out:
        out.write("Processing experimental samples\n")

    df = pd.read_csv('./dataset_metadata.csv')
    # Create dictionary from dataframe
    pdb_id_path_dict = pd.Series(df.pdb_path.values,index=df.fitted_entries).to_dict()
    # Remove PDB extra info
    df["fitted_entries"].map(lambda pdb_id: pdb_split_and_save(pdb_id, pdb_id_path_dict[pdb_id]))
    # Compute number of subunits
    df['subunit_count'] = df["fitted_entries"].map(lambda pdb_id: pdb_num_chain_mapper(pdb_id, pdb_id_path_dict[pdb_id]))
    # Get number of models with less than 2 subunits
    non_enough_subunit = len(df[df.subunit_count<2].index)
    # Remove entries with less than 2 subunits
    df = df[df.subunit_count>=2]

    df.to_csv('dataset_metadata.csv', index=False)
     
    number_of_samples = len( df.index )
   
    with open(out_txt_path, 'a') as out:
        out.write("Removing {} experimental entries with less than 2 subunits\n".format(non_enough_subunit))
        out.write("{} experimental candidates for dataset\n".format(number_of_samples))
    
    with open(out_txt_path, 'a') as out:
        out.write("Processing simulated samples\n")
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
   
    with open(out_txt_path, 'a') as out:
        out.write("Removing {} pdb entries with less than 2 subunits\n".format(non_enough_subunit))
        out.write("{} pdb candidates for simulated dataset\n".format(number_of_samples))

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
    if map_box_max_dim != np.min(map_box) or map_box_max_dim % 2 != 0:
        # Change map extension to mrc due to compatibility with EMAN
        old_map_filename = map_filename
        map_filename = map_filename.replace('.map','.mrc')
        map_box_max_dim = map_box_max_dim if map_box_max_dim % 2 == 0 else map_box_max_dim+1
        command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2proc3d.py '+ old_map_filename + ' ' + map_filename +' --clip='+str(map_box_max_dim)+','+str(map_box_max_dim)+','+str(map_box_max_dim)
        print(command)
        try:
            execute_command(command)
        except Exception as e:
            with open('error.txt', 'a') as error_file:
                error_file.write("Error resizing map {}: {}\n".format(os.path.basename(old_map_filename), e))
        old_box = map_box
        map_object = molecule.Molecule(map_filename, recommendedContour=contourLevel, cutoffRatios=[1])
        map_box = map_object.getGridSize() 
        with open("rezised.txt", 'a') as txt:
            txt.write("Map {} with size {} has been resized to {} and saved as {}\n".format(os.path.basename(old_map_filename), old_box, map_box, os.path.basename(map_filename)))

    # Get map volume size in voxels
    map_volume = map_object.getVolume()[1]

    simulated_path = os.path.join(sim_model_path, 'sim_'+os.path.basename(map_filename).replace('.map','.mrc'))
    # Generate map
    #command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2pdb2mrc.py -A=1' + ' -R=' + str(res) + ' -B='+ str(map_box[2]) +','+ str(map_box[1])+ ','+ str(map_box[0]) +' '+  pdb_filename + ' ' + simulated_path
    command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2pdb2mrc.py --center --quiet -B='+ str(map_box[2]) +','+ str(map_box[1])+ ','+ str(map_box[0]) +' '+  pdb_filename + ' ' + simulated_path
    print(command)
    try:
        execute_command(command)
    except Exception as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Error simulating map {}: {}\n".format(os.path.basename(simulated_path), e))

    volume_result = dict()
    volume_result['index']=index
    volume_result['simulated_path'] = simulated_path
    volume_result['map_volume']=map_volume
    volume_result['map_path']=map_filename
    volume_result['pdb_path'] = pdb_filename
    volume_result['map_length'] = map_box[0]
    
    # Get simulated map volume
    if os.path.exists(simulated_path):
        # Simulated map contour must be a small number but not zero becouse eman sometimes generates a cube at 0 contour value
        sim_object = molecule.Molecule(simulated_path, recommendedContour=0.001, cutoffRatios=[1])
        # Get dictionary of volumes, choose element with key 1 (corresponding to 100% recommended contour level)
        pdb_volume = sim_object.getVolume()[1]
        # Simulated volume at contour level 0 sometimes gives a cube, so shuld use an epsilon instead.
        volume_result['pdb_volume']=pdb_volume
        volume_result['vol_capture'] = pdb_volume/map_volume

    else:
        volume_result['pdb_volume']=-1
        with open('error.txt', 'a') as error_file:
            error_file.write("Simulated map {} does not exist: {}\n".format(os.path.basename(simulated_path), e))

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
    result['contour']=0.001
    # Generate map
    try:
        #command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2pdb2mrc.py -A=1.0 -R=' + str(res) + ' --center '+  pdb_filepath + ' ' + simulated_path
        command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2pdb2mrc.py  --center '+  pdb_filepath + ' ' + simulated_path
        print(command)
        execute_command(command)
        map_object = molecule.Molecule(simulated_path, recommendedContour=0.001, cutoffRatios=[1])
        map_box = map_object.getGridSize()
        map_box_max_dim = np.max(map_box)
        if map_box_max_dim % 2 != 0:
            # Change map extension to mrc due to compatibility with EMAN
            map_box_max_dim = map_box_max_dim if map_box_max_dim % 2 == 0 else map_box_max_dim+1
            command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2proc3d.py '+ simulated_path + ' ' + simulated_path +' --clip='+str(map_box_max_dim)+','+str(map_box_max_dim)+','+str(map_box_max_dim)        
            print(command)
            execute_command(command)
            old_box = map_box
            map_object = molecule.Molecule(simulated_path, recommendedContour=0.001, cutoffRatios=[1])
            map_box = map_object.getGridSize()
        result['map_length'] = int(map_box[0])
        with open("generated.txt", 'a') as txt:
            txt.write("Map {} with size {} has been generated from {} \n".format(os.path.basename(simulated_path), map_box, os.path.basename(pdb_filepath)))
    except Exception as e:
        result['resolution']=-1
        result['map_length']=0
        with open('error.txt', 'a') as error_file:
            error_file.write("Error generating map {}: {}\n".format(os.path.basename(pdb_filepath), e))

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
    # OLD BINARY MASK
    #command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2proc3d.py ' + map_path + ' '+ mask_path +' --process threshold.binary:value=' + str(contourLvl)
    #command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2proc3d.py ' + map_path + ' '+ mask_path +' --process=mask.auto3d:nmaxseed=12:nshells=3:nshellsgauss=3:return_mask=1:threshold='+str(contourLvl)
    #print(command)
    #execute_command(command)
    
    try:
        simulated_object = molecule.Molecule(simulated_path, recommendedContour=0.001, cutoffRatios=[1])
        # Compute correlation and overlap after alignment
        map_object = molecule.Molecule(map_path, recommendedContour=contourLvl, cutoffRatios=[1])
        overlap_before = map_object.getOverlap(simulated_object)[1]
        corr_before = map_object.getCorrelation(simulated_object)[1]
        # Compute alignment
        #command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2proc3d.py  --align=rotate_translate_3d_tree --alignref='+simulated_path+' --multfile='+ mask_path +' '+map_path+' '+aligned_map_path
        command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2proc3d.py  '+map_path+' '+aligned_map_path
        print(command)
        execute_command(command)
        map_object = molecule.Molecule(aligned_map_path, recommendedContour=contourLvl, cutoffRatios=[1])
        overlap_after = map_object.getOverlap(simulated_object)[1]
        corr_after = map_object.getCorrelation(simulated_object)[1]
        corr_after = corr_after if not np.isnan(corr_after) else 0
    except RuntimeError as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Error aligning map {}: {}\n".format(os.path.basename(map_path), e))

    return {'index':index, 'aligned_path':aligned_map_path,'overlap_before': overlap_before, 'corr_before':corr_before, 'overlap_after':overlap_after, 'corr_after':corr_after}
     
# Generate simulated segments using EMAN
def generateSimulatedSegments(index, df, sim_model_path):
    if not os.path.exists(sim_model_path):
        os.makedirs(sim_model_path)
    pdb_filename = df.at[df.index[index], 'pdb_path']
    simulated_path = df.at[df.index[index], 'simulated_path']
    chain_id = df.at[df.index[index], 'chain_id']
    generated_segment_path = os.path.join(sim_model_path, os.path.basename(simulated_path.replace('.mrc','_'+chain_id+'.mrc')))
    map_box_length =df.at[df.index[index], 'map_length']
    # Get centered pdb
    pdb_filename = pdb_filename.replace('.ent','_centered.ent')
    pdb_filename_chain = pdb_filename.replace('.ent','_'+chain_id+'.ent')

    result = {}
    result['index'] = index
    try:
        #command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2procpdb.py --chains ' + chain_id + ' ' + pdb_filename + ' ' + pdb_filename_chain
        #print(command)
        #execute_command(command)
        # Generate map without specify box dimension
        command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2pdb2mrc.py -B '+str(int(map_box_length))+ '  '+ pdb_filename_chain + ' ' + generated_segment_path
        print(command)
        execute_command(command)
    except Exception as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Error generating segment map {} for pdb {}: {}\n".format(chain_id,os.path.basename(pdb_filename),e))
    # Get simulated map volume
    if os.path.exists(generated_segment_path):
        map_object = molecule.Molecule(generated_segment_path, recommendedContour=0.001, cutoffRatios=[1])        
        segment_volume = map_object.getVolume()[1]
        # Get subunit path back to main thread
        result['subunit_path'] = generated_segment_path
        result['segment_volume'] = segment_volume
    else: 
        result['subunit_path'] = ''
        result['segment_volume'] = -1
        with open('error.txt', 'a') as error_file:
            error_file.write("Error reading segment map {} from simulated {}: {}\n".format(chain_id,generated_segment_path,e)) 
    return result
    




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
def selectExperimentalDataset(result_dir, out_txt_path, percentil_boundaries_tuple=None):
    # Open result csv from volume calculation
    df = pd.read_csv('dataset_volume.csv')
    # Replace missing volumes (produced by e2pdb2mrc failure to generate map from pdb)
    df['pdb_volume'].replace('', np.nan, inplace=True)
    df['map_volume'].replace('', np.nan, inplace=True)
    df['vol_capture'].replace('', np.nan, inplace=True)
    df.dropna(inplace=True)

    # Plot hist
    plot_hist(df.vol_capture, '', os.path.join(result_dir,'data_capture_hist.eps'))   
    column_title_dict = {"map_volume":"volume exp", "pdb_volume":"vol sim", "vol_capture":"capture ratio"}
    plot_box(df, column_title_dict, os.path.join(result_dir, 'data_boxplot'))

    # Computing basic stats
    num_samples = len(df)
    avg = df['vol_capture'].mean()
    med = df['vol_capture'].median()
    std = df['vol_capture'].std()
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
    
    percentilA = df.vol_capture.quantile(int(percentil_boundaries_tuple[0])/100, interpolation = 'midpoint')
    percentilB = df.vol_capture.quantile(int(percentil_boundaries_tuple[1])/100,  interpolation = 'midpoint')

    # Print to file
    with open(out_txt_path, 'a') as out:
        out.write("Processing emdb data\n")
        out.write("Number of samples: {}\n".format(num_samples))
        out.write("Volume ratio avg: {}\n".format(avg))
        out.write("Volume ratio median: {}\n".format(med))
        out.write("Volume ratio std: {}\n".format(std))
        out.write("Printing vol capture distribution in percentiles\n")
        # Printing percentile values
        for i in range(0,100):
            print("Percentile {} value: {}\n".format(i, df.vol_capture.quantile(i/100, interpolation = 'midpoint')))
        out.write("Quartile Q1: {}\n".format(q1))
        out.write("Quartile Q3: {}\n".format(q3))
        out.write("IQR: {}\n".format(q3-q1))
        out.write("Lower boundary: {}\n".format(low_boundary))
        out.write("Upper boundary: {}\n".format(upper_boundary))
        out.write("Number of samples below lower boundary: {} and over upper boundary: {}\n".format(len( df[df.vol_capture < low_boundary] ), len( df[df.vol_capture >upper_boundary] )))

    if percentil_boundaries_tuple != None:
        selected_samples = df[ (df.vol_capture > percentilA) & (df.vol_capture < percentilB)]
        numMuestras = len( selected_samples)
        with open(out_txt_path, 'a') as out:
            out.write("Number of samples between percentil {} and {}: {}\n".format(percentil_boundaries_tuple[0],percentil_boundaries_tuple[1], numMuestras))
    else: 
        selected_samples = df[(df.vol_capture < low_boundary) | (df.vol_capture >upper_boundary) ]
        with open(out_txt_path, 'a') as out:
            out.write("Getting outliers following IQR...\n")
    # Save csv with selected samples
    selected_samples.to_csv('dataset_selected.csv', index = False)
    #Plot Number of protein subunits of EMDB maps
    plot_hist(df.subunit_count,'', os.path.join(result_dir,'data_subunitcount_hist.eps'))
    #Plot Protein volume values of EMDB maps
    plot_hist(df.map_volume, '', os.path.join(result_dir,'data_volume_exp_hist.eps'))
    # Plot Protein volume values of generated maps
    plot_hist(df.pdb_volume, '', os.path.join(result_dir,'data_volume_sim_hist.eps'))


def generateStatsFromSelectedData(df, result_dir):

    # Generate charts from selected data
    # Generate resolution histogram from EMDB data
    plot_hist(df.resolution,'', os.path.join(result_dir,'selected_resolution_hist.eps'))
    plot_hist(df.subunit_count,'', os.path.join(result_dir,'selected_subunitcount_hist.eps'))
    plot_hist(df.map_volume, '', os.path.join(result_dir,'selected_volume_exp_hist.eps'))
    plot_hist(df.pdb_volume, '', os.path.join(result_dir,'selected_volume_sim_hist.eps'))
    plot_hist(df.corr_before, '', os.path.join(result_dir,'selected_corr_before_hist.eps'))
    plot_hist(df.corr_after, '', os.path.join(result_dir,'selected_corr_after_hist.eps'))
    plot_hist(df.overlap_before, '', os.path.join(result_dir,'selected_overlap_before_hist.eps'))
    plot_hist(df.overlap_after, '', os.path.join(result_dir,'selected_overlap_after_hist.eps'))
 
    
def mergeMapAndChains(datasetSource, columns):
    # Result dataset
    dataset_result = pd.DataFrame(columns=datasetSource.columns.values.tolist())
    # Get Molecule segments from pdb lecture with lib 
    chain_df_to_merge = pd.DataFrame(columns=columns)
    for index, row in datasetSource.iterrows():
        entryid = row[columns[0]]
        filepath = row[columns[2]]
        centered_filepath = filepath.replace('.ent','_centered.ent')
        # Move pdb coords to center of mass
        command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2procpdb_mod.py ' + filepath + ' ' + centered_filepath + ' --center 1'
        execute_command(command)
        # Split centered pdb into chains
        pdb_split_and_save(entryid,centered_filepath)
        subunit_id_list, subunit_labels = pdb_get_chain_id_label_list(entryid, filepath)
        list_entry_name = [entryid for i in subunit_id_list]
        chain_df_to_merge = pd.DataFrame({columns[0]:list_entry_name,columns[3]:subunit_id_list, columns[1]:subunit_labels})
        dataset_df_to_merge = datasetSource[datasetSource[columns[0]]==entryid]
        cross_join = pd.merge(dataset_df_to_merge, chain_df_to_merge, on=columns[0], how='outer')
        dataset_result = dataset_result.append(cross_join, sort=False)
    # Reset indexes
    dataset_result.rename(columns = {columns[2]:'pdb_path', columns[4]:'simulated_path'}, inplace = True) 
    dataset_result.reset_index(inplace=True, drop=True)               
    return dataset_result

def parallelSimulationSegments(df, simulated_path, output_csv_name):
    # Get index list to schedule processess 
    index_list = df.index.tolist()
    print("Spawn procecess...")
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    with MPICommExecutor(comm, root=0, worker_size=size) as executor:
        if executor is not None:
            futures = []
            for i in index_list:
                futures.append(executor.submit(generateSimulatedSegments, i, df, simulated_path))
            wait(futures)
            for f in futures:
                try:
                    d = f.result()
                    print("received simulated segments dict: ",d)
                    res_index = d['index']
                    subunit_path = d['subunit_path']
                    segment_volume = d['segment_volume']
                    df.loc[res_index, 'subunit_path'] = subunit_path
                    df.loc[res_index, 'segment_volume'] = segment_volume

                except Exception as error:
                    print("Error calculating volume for simulated {}: {}".format(df.loc[res_index],error))
            df.to_csv(output_csv_name, index=False)


def main():
    global current_dir


    parser = argparse.ArgumentParser()
    parser.add_argument('--d', required=False, default=False, help='True if want to resync emdb data')
    parser.add_argument('--v', required=False, default=False, help='True if want to compare map volume with its simulated counterpart')
    parser.add_argument('--s', required=False, default=None, nargs=2, help='Please provide lower and upper experimental data percentil distribution to be included in final dataset')  
    parser.add_argument('--g', required=False, default=None, help='Please provide percent of simulated data to be included in final dataset')
    parser.add_argument('--f', required=False, default=False, help='True if want to align experimental data to simulated counterpart')
    parser.add_argument('--p', required=False, default=False, help='Generate segments')
    parser.add_argument('--header_dir', default='header', required=False,help='output directory to store emdb headers') 
    parser.add_argument('--models_dir', default='models', required=False, help='output directory to store emdb models') 
    parser.add_argument('--simulated_dir', default='simulated', required=False, help='output directory to store simulated maps') 
    parser.add_argument('--result_dir', default='output', required=False, help='output directory to store stats') 

    opt = parser.parse_args()
    current_dir = os.getcwd()
    header_path = os.path.join(current_dir, opt.header_dir)
    models_path = os.path.join(current_dir, opt.models_dir)
    simulated_path = os.path.join(current_dir, opt.simulated_dir)
    results_path = os.path.join(current_dir, opt.result_dir)
    output_txt_path = os.path.join(results_path, 'output.txt')
    pdb_path = os.path.join(models_path,'pdb')

    # Set seed 
    seed = 42 

    # Download and process data
    if int(opt.d):
        get_headers(header_path)
        generate_dataframe(header_path)
        processMetadata(output_txt_path)
        downloadModels(models_path, output_txt_path)
        downloadModelsPDB(pdb_path, models_path, output_txt_path)
        removeNonExistingModels(models_path, output_txt_path)
        processfiles(models_path, output_txt_path)
    #Calculate volume
    if int(opt.v):
        df = pd.read_csv('./dataset_metadata.csv')
        with open(output_txt_path, 'a') as out:
            out.write("Number of samples to process volume: {}\n".format( len(df.index)))
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
                wait(futures)
                for f in futures:
                    try:
                        d = f.result()
                        print("received vol dict: ",d)
                        res_index = d['index']
                        pdb_volume = d['pdb_volume']
                        map_volume = d['map_volume']
                        capture_ratio = d['vol_capture']
                        map_path = d['map_path']
                        pdb_path = d['pdb_path']
                        map_length = d['map_length']
                        sim_path = d['simulated_path']
                        df_volume.loc[res_index, 'simulated_path'] = sim_path
                        df_volume.loc[res_index, 'map_path'] = map_path
                        df_volume.loc[res_index, 'pdb_path'] = pdb_path
                        df_volume.loc[res_index, 'map_length'] = map_length
                        df_volume.loc[res_index, 'map_volume'] = map_volume
                        df_volume.loc[res_index, 'pdb_volume'] = pdb_volume
                        df_volume.loc[res_index, 'vol_capture'] = capture_ratio
                    except ValueError as error:
                        print("Error calculating volume for simulated {}: {}".format(df_volume.loc[i,'fitted_entries'],error))
                df_volume.to_csv('dataset_volume.csv', index=False)
    if opt.s != None:
        selectExperimentalDataset(results_path, output_txt_path, percentil_boundaries_tuple=opt.s)
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
                wait(futures)
                for f in futures:
                    try:
                        d = f.result()
                        print("received res dict: ",d)
                        res_index = d['index']
                        overlap_before = d['overlap_before']
                        corr_before = d['corr_before']
                        overlap_after = d['overlap_after']
                        corr_after = d['corr_after']
                        path = d['aligned_path']
                        df_selected.loc[res_index, 'aligned_path'] = path
                        df_selected.loc[res_index, 'overlap_before'] = overlap_before
                        df_selected.loc[res_index, 'corr_before'] = corr_before
                        df_selected.loc[res_index, 'overlap_after'] = overlap_after
                        df_selected.loc[res_index, 'corr_after'] = corr_after
                    except ValueError as error:
                        print("Error fitting map to simulated counterpart: {}".format(error))

                df_selected.to_csv('dataset_selected.csv', index=False)
        # Generate selected data stats
        generateStatsFromSelectedData(df_selected, results_path)

    # Generate simulated data to be included in final dataset
    if opt.g != None:
        # Open candidates to compute simulated map 
        df_synthetic = pd.read_csv('metadata_synthetic.csv')
        df_experimental = pd.read_csv('dataset_selected.csv')
        num_experimental = len(df_experimental)
        proportion = int(opt.g)/100
        num_samples_to_pick = int(round(proportion*num_experimental/1-proportion))
        with open(output_txt_path, 'a') as out:
            out.write("Number of candidate synthetic models: {}\n".format(len(df_synthetic.index)))
            out.write("A proportion of {} synthetic data in final dataset requires {} simulated and {} experimental samples.\n".format(proportion,num_samples_to_pick,num_experimental ))
        
        df_synthetic_selected = df_synthetic.sample(num_samples_to_pick, random_state=seed)
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
                wait(futures)
                for f in futures:
                    try:
                        d = f.result()
                        print("received res dict: ",d)
                        res_index = d['index']
                        res = d['resolution']
                        path = d['path']
                        map_length = d['map_length']
                        contour = d['contour']
                        df_synthetic_selected.loc[res_index, 'contourLevel'] = contour
                        df_synthetic_selected.loc[res_index, 'resolution'] = res
                        df_synthetic_selected.loc[res_index, 'map_path'] = path
                        df_synthetic_selected.loc[res_index, 'map_length'] = map_length
                    except Exception as error:
                        print("Error calculating volume for simulated {}: {}".format(i,error))
                df_synthetic_selected.to_csv('synthetic_selected.csv', index=False)
    
    
    if int(opt.p):
        # Merge data from experimental and synthetic sources, extract subunits from maps and save them in output directory
	
        with open(output_txt_path, 'a') as out:
            out.write("Generating segments with EMAN..\n")
        dataset_cryoem = pd.read_csv('dataset_selected.csv')
        dataset_synthetic = pd.read_csv('synthetic_selected.csv')

        subunits_path = os.path.join(simulated_path, 'subunits')

        # Create output directory for simulated segmetns
        if not(os.path.exists(subunits_path)):
            os.system('mkdir '+subunits_path)
        
        dataset_cryoem_segments = mergeMapAndChains(dataset_cryoem, ['fitted_entries','chain_label','pdb_path', 'chain_id','simulated_path'])
        with open(output_txt_path, 'a') as out:
            out.write("Experimental maps merged..\n")
        parallelSimulationSegments(dataset_cryoem_segments, simulated_path, 'dataset_exp_merged.csv')
        dataset_synthetic_merged = mergeMapAndChains(dataset_synthetic, ['entries','chain_label','path', 'chain_id', 'map_path'])
        with open(output_txt_path, 'a') as out:
            out.write("Simulated maps merged..\n")
        parallelSimulationSegments(dataset_synthetic_merged, simulated_path, 'dataset_sim_merged.csv')
                
        
     # Simulate give segment with EMAN and save file following pdbid_segment.mrc fale notation
      
        # Create Dataframe

        # Do segment annotation by multiplying the given mask the by int of each region
      
        # Check for voxel locality with exact and flexible localitty with a metric

        # Complete 
        # Save as files as .map 


         


       





       


       
        


if __name__ == '__main__':
    main()
