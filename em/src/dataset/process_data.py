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
from scipy.ndimage import zoom
import copy 
from miscellaneous import *
import mrcfile
from metrics import getCorrelation, getRelative_Masks_Overlap

# Sync headers to folder
def get_headers(header_path):
    if not os.path.exists(header_path):
        os.makedirs(header_path)
    print("Getting headers...")
    command = 'rsync -rlpt --ignore-existing -v -z --delete --port=33444 --include "emd-[0-9]+.xml" "rsync.rcsb.org::emdb/structures/EMD-*/header/" '+ header_path
    execute_command(command)


def clean(models_path):
    df = pd.read_csv('./dataset_valid.csv', dtype=str)

    pdb_files = glob(os.path.join(models_path,'*.ent'))
    pdb_valid = df['fitted_entries'].to_list()
    pdb_to_remove = [ pdb_id for pdb_id in pdb_files if os.path.basename(pdb_id)[3:-4] not in pdb_valid ]
    command = ' '.join(pdb_to_remove)
    try:
        execute_command("rm "+command)
    except Exception as e:
        print(e)
# Use function to map names
def restore_id(map_id):
    min_len = 4
    actual_len = len(map_id)
    if actual_len<min_len:
        diff = min_len-actual_len
        return diff*'0'+map_id
    else:
        return map_id

# Scan headers and create a dataframe with all samples
def generate_dataframe(header_path):
    # Get list of xml headers
    all_xml = glob(os.path.join(header_path,r"emd-*-v30.xml"))
    # Creates a dict with map path as key and map ID as value
    id_path_dict = {x:os.path.splitext(os.path.basename(x))[0][4:-4] for x in all_xml}
    # Creates a dataframe to store data
    df = pd.DataFrame(all_xml, columns=["path"])
    # Creates map id column
    df["id"] = df["path"].map(id_path_dict.get)
    
    df['fitted_entries']=None
    df['method']= None
    df['resolution']=None
    df['contourLevel']=None
    df['name']=None

    #iterate over all xml files
    for i,meta in enumerate(all_xml):
        #parse xml content
        with open(meta, 'rt') as f:
            tree = ElementTree.parse(f)
            resolution_node=None
            pdbe_fitted_entries = []
            for pdbe_entry in tree.findall('.//crossreferences/pdb_list/pdb_reference/pdb_id'):
                pdbe_fitted_entries.append(pdbe_entry.text)

            for child in tree.find('.//structure_determination_list/structure_determination'):
                if 'method' in child.tag:
                    processing_method_node = child
                if 'processing'  in child.tag:
                    for grandson in child:
                        if grandson.tag == 'final_reconstruction':
                            for ggrandson in grandson:
                                if ggrandson.tag == 'resolution':
                                    resolution_node = ggrandson
            contour_node = tree.find('.//map/contour_list/contour/level')
            name_node = tree.find('.//sample/name')
        # If has several entries, discard.
        if len(pdbe_fitted_entries)>0:
        #    if len(pdbe_fitted_entries)>1:
        #        print("Model {} has more than one fitted structures, discarting.".format(id_path_dict[meta]))
        #        continue
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
        if name_node is not None:
            df.loc[ df.id == id_path_dict[meta], ['name'] ] = name_node.text
    print("Processed {} files to create a csv with {} samples".format(len(all_xml),len(df)))
    df.to_csv('dataset_metadata.csv', index=False)


# Process initial metadata, drop records that don't satisfy criteria
def processMetadata(output_txt_path):
    df = pd.read_csv('./dataset_metadata.csv')
    total_num = len(df.index) 
    # Get number of null elements to be removed
    #non_valid= len(df.dropna().index)
    # Remove null values
    df_no_structure = df[df['fitted_entries'].isna()]
    df = df.dropna()
    # Number of samples out of resolution gap
    outside_gap = len(df[(df['resolution']<4.5) | (df['resolution']>10)])
    # Keep only samples within resolution gap
    df = df[(df.resolution>=4.5) & (df.resolution<=10)]

    # Print stats to output file
    with open(output_txt_path, 'w') as out:
        out.write("{} total entries in EMDB\n".format(total_num))
        out.write("{} invalid entries\n".format(len(df_no_structure)))
        out.write("{} entries out of resolution gap\n".format(outside_gap))
        out.write("This leads to {} candidates for dataset\n".format(len(df.index)))

    # Change dtype id to string and restore zeros before when needed
    df = df.astype({'id':str})

    # Use function to map names
    df['id'] = df['id'].map(lambda map_id: restore_id(map_id))
    df["pdb_rsync"] = df["fitted_entries"].map(lambda pdb_id: "rsync.rcsb.org::ftp_data/structures/all/pdb/pdb"+str(pdb_id)+".ent.gz")
    df['mmcif_rsync'] = df['fitted_entries'].map( lambda pdb_id:  "rsync.rcsb.org::ftp_data/structures/divided/mmCIF/"+str(pdb_id)[1:-1]+"/"+str(pdb_id)+".cif.gz")
    df["emdb_rsync"] = df["id"].map(lambda map_id: "rsync.rcsb.org::emdb/structures/EMD-"+str(map_id)+"/map/emd_"+str(map_id)+".map.gz")


    df.to_csv('dataset_metadata_selected.csv', index=False)

# Download candidate experimental models from database
def downloadModels(models_path, output_txt_path):
    if not os.path.exists(models_path):
        os.makedirs(models_path)
    df = pd.read_csv('./dataset_metadata_selected.csv', dtype=str)
    
    # Check existing models
    with open(output_txt_path, 'a') as out:
        out.write("Check existing files to download new models only..\n")
    pdb_files = [f for f in glob(os.path.join(models_path,'*.pdb'))]
    map_files = [f for f in glob(os.path.join(models_path,'*.mrc'))]

    pdb_filenames = [os.path.basename(f) for f in pdb_files]
    map_filenames = [os.path.basename(f) for f in map_files]

    with open(output_txt_path, 'a') as out:
        out.write("{} pdbs and {} maps found\n".format(len(pdb_files), len(map_files)))

    #Get entries with missing pdb or map
    df['map_found'] = df["id"].map(lambda map_id: True if map_id+'.mrc' in map_filenames else False)
    df['pdb_found'] = df["fitted_entries"].map(lambda pdb_id: True if pdb_id+'.pdb' in pdb_filenames else False)

    # Choose only non existing files to download
    df_maps = df[df['map_found'] == False] 
    df_pdbs = df[df['pdb_found'] == False] 
    with open(output_txt_path, 'a') as out:
        out.write("Downloading {} pdbs and {} maps \n".format(len(df_pdbs), len(df_maps)))
    emdb_ftp_list = df_maps["emdb_rsync"].tolist()
    pdb_ftp_list = df_pdbs["pdb_rsync"].tolist()
    mmcif_ftp_list = df_pdbs["mmcif_rsync"].tolist()
    emdb_id_list = df_maps["id"].tolist()
    pdb_id_list = df_pdbs["fitted_entries"].tolist()
     
    for uri,name in zip(emdb_ftp_list, emdb_id_list):
        command = 'rsync -rlpt --ignore-existing -v -z --delete --port=33444 '+ uri  +' '+ models_path+'/'
        try:
            execute_command(command) 
        except Exception as e:
            with open('error.txt', 'a') as error_file:
                error_file.write("Error excecuting command {}\n".format(command))
    
    for pdb_uri,mmcif_uri,name in zip(pdb_ftp_list, mmcif_ftp_list, pdb_id_list):
        command = 'rsync -rlpt --ignore-existing -v -z -L --delete --port=33444 '+  pdb_uri  +' '+ models_path+'/'
        try:
            execute_command(command)    
        except Exception as e:
            # try to download mmCIF
            print("Trying to download mmcif")
            command = 'rsync -rlpt --ignore-existing -v -z -L --delete --port=33444 '+  mmcif_uri  +' '+ models_path+'/'
            try:
                execute_command(command)
            except Exception as e:
                with open('error.txt', 'a') as error_file:
                    error_file.write("Error excecuting command {}\n".format(command))
    print("Extracting donwloaded models from emdb")
    command = 'find '+models_path+' -maxdepth 1 -name \'*.gz\' -exec gunzip -f -k -n {}  \;'
    try:
        execute_command(command)
    except Exception as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Error excecuting command {}\n".format(command))

    print("Rename map name  to map_id.mrc")
    command = 'rename \'s/emd_//\' models/*.map\;'
    try:
        execute_command(command)
    except Exception as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Error excecuting command {}\n".format(command))

    command = 'rename \'s/map/mrc/\' models/*.map\;'
    try:
        execute_command(command)
    except Exception as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Error excecuting command {}\n".format(command))

    
    command = 'rename \'s/pdb//\' models/*.ent\;'
    try:
        execute_command(command)
    except Exception as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Error excecuting command {}\n".format(command))

    command = 'rename \'s/ent/pdb/\' models/*.ent\;'
    try:
        execute_command(command)
    except Exception as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Error excecuting command {}\n".format(command))


    # Convert mmCIFF to pdb model 
    command = 'find '+models_path+' -maxdepth 1 -name \'*.cif\' -exec maxit -input {} -output {}.pdb -o 2   \;'
    try:
        execute_command(command)
    except Exception as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Error excecuting command {}\n".format(command))
    # Rename .cif.pdb to .pdb
    print("Rename cif.pdb to .pdb")
    command = 'rename \'s/.cif//\' models/*.cif.pdb '
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
    
    pdb_files = [f for f in glob(os.path.join(models_path,'*.pdb'))]

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
    df['pdb_found'] = df["fitted_entries"].map(lambda pdb_id: True if pdb_id+'.pdb' in pdb_filenames else False)

    df_map_not_found = df[df["map_found"]==False]
    df_pdb_not_found = df[df["pdb_found"]==False]
    # Remove entries that doesn't have neither map or pdb
    indexes_to_remove = df[ (df['map_found'] == False) | (df['pdb_found'] == False) ].index
    #indexes_to_remove = df[ df['pdb_found'] == False ].index
    df = df.drop(indexes_to_remove)
    
    df['map_path'] = df["id"].map(lambda map_id: map_id_path_dict[map_id.replace('-','_')+'.map'])
    df['pdb_path'] = df["fitted_entries"].map(lambda pdb_id: pdb_id_path_dict[pdb_id+'.pdb'])
    
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
        pdb_filename_chain = filepath.replace('_nohet.pdb','_'+chain._id+'.pdb')    
        io.set_structure(model[chain._id])
        io.save(pdb_filename_chain)
    print("{} splitted and saved".format(filepath)) 

# Define function to get number of chains in structure
def pdb_num_chain_mapper(chain, filepath):
    parser = PDBParser(PERMISSIVE = True, QUIET = True)
    io = PDBIO()
    num_chains = 0
    try:
        with open(filepath+'/'+chain+'.pdb', "r") as inputFile,open(filepath+'/'+chain+'_nohet.pdb',"w") as outFile:
            for line in inputFile:
                if not line.startswith("HETATM"):
                    outFile.write(line)
        pdb_obj = parser.get_structure(chain, filepath+'/'+chain+'_nohet.pdb')
        num_chains = len(list(pdb_obj.get_chains()))
        print("Structure {} chains num: {}\n".format(chain,num_chains))
    except Exception as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Error: Cant open pdb file {}: {}\n".format(chain, e))
    return num_chains

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
    df = pd.read_csv('./dataset_metadata_selected.csv', dtype=str)
    number_selected = len(df)
    # Check downloaded EM map to resize to 1A per voxel
    for i,row in df.iterrows():
        map_path = "models/"+row['id']+".mrc"
        if os.path.exists(map_path):
            with mrcfile.open(map_path, 'r+') as map_object:
                map_shape = list(map_object.data.shape)
                map_shape_in_angstroms = map_object.header['cella'].tolist()
                scaling_factor = [ y/x for x,y in zip(map_shape,map_shape_in_angstroms) ]
                scaling_factor_rounded = [round(x,1) for x in scaling_factor]
                print("Map with size {} and {} A/voxel, scale by {} to {}. ".format(row['id'], map_object.voxel_size.tolist(), scaling_factor_rounded, map_object.header[['nx','ny','nz']].tolist()))
                if ( (scaling_factor_rounded[0] != 1.0) | (scaling_factor_rounded[1] != 1.0) | (scaling_factor_rounded[2] != 1.0) ):
                    scaled_data = zoom(map_object.data, scaling_factor, order=1, grid_mode=True, mode='grid-constant')
                    new_shape = scaled_data.shape
                    #with open(out_txt_path, 'a') as out:
                    #    out.write("Map {} resized from {} to {} with scaling factor {}\n".format(row['id'], map_shape,new_shape, scaling_factor_rounded))
                    map_object.set_data(scaled_data)
                    new_map_shape = list(map_object.data.shape)
                    new_map_shape_in_angstroms = map_object.header['cella'].tolist()
                    new_voxel_size = [round(y/x,1) for x,y in zip(new_map_shape,new_map_shape_in_angstroms)]
                    assert new_voxel_size == [1.0, 1.0, 1.0], "Voxel size is not 1 for resized Map {}: {}".format(map_path,new_voxel_size)
                else: 
                    assert scaling_factor_rounded ==  [1.0, 1.0, 1.0], "Voxel size is not 1 for unresized Map {}: {}".format(map_path,scaling_factor_rounded)
        else:
            with open('error.txt', 'a') as error_file:
                error_file.write("Error: Map {} not found\n".format(map_path))



    df['subunit_count'] = df["fitted_entries"].map(lambda pdb_id: pdb_num_chain_mapper(pdb_id, models_path))
    # Get number of models with less than 2 subunits
    non_enough_subunit = len(df[df.subunit_count<2].index)
    # Remove entries with less than 2 subunits
    df = df[df.subunit_count>=2]

    df.to_csv('dataset_metadata_selected.csv', index=False)
     
    number_after_clean = len(df)
   
    with open(out_txt_path, 'a') as out:
        out.write("Removing {} of {} experimental entries with less than 2 subunits\n".format(non_enough_subunit, number_selected))
        out.write("{} experimental candidates for dataset\n".format(number_after_clean))
    '''
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
    '''
    

# Simulate map from pdb structure and calculate volume 
def simulateMapAndCompareVolume(index, df, sim_model_path):
    if not os.path.exists(sim_model_path):
        os.makedirs(sim_model_path)
    map_filename = 'models/'+df.at[df.index[index], 'id']+'.mrc'
    pdb_filename = 'models/'+df.at[df.index[index], 'fitted_entries']+'.pdb'
    res = float(df.at[df.index[index], 'resolution'])
    contourLevel = float(df.at[df.index[index], 'contourLevel'])
    try:
        # Get map volume
        map_object = mrcfile.open(map_filename, 'r', permissive=True)
        map_data = map_object.data
        map_shape = map_data.shape
    except Exception as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Failed to open map {}:{}\n".format(map_filename,e))
            volume_result['vol_overlap'] = None
            volume_result['vol_correlation'] = None
            return volume_result

    simulated_path = os.path.join(sim_model_path, 'sim_'+os.path.basename(map_filename))
    # Prepare command
    tmp_filename = "/tmp/{}_cmd.txt".format(os.path.basename(simulated_path))
    with open(tmp_filename, 'w') as tmp_file:
        tmp_file.write("1\n1\n1\n-{}\n1\n1\n1\n".format(str(res)))
    command =   "cat {}  |  pdb2vol {} {}".format(tmp_filename, pdb_filename, simulated_path)
    
    try:
        execute_command(command)
        #print(command)
    except Exception as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Error simulating map {}: {}\n".format(os.path.basename(simulated_path), e))

    volume_result = dict()
    volume_result['index']=index
    #volume_result['simulated_path'] = simulated_path
    volume_result['map_length'] = ','.join([str(map_shape[0]),str(map_shape[1]),str(map_shape[2])])
    # Get simulated map volume
    if os.path.exists(simulated_path):
        try:
            sim_map_object = mrcfile.open(simulated_path, 'r', permissive=True)
            sim_data = sim_map_object.data
            sim_shape = sim_data.shape
            # Open simulated padded object to store resized simulated map
            sim_padded_path = simulated_path.replace('.mrc','_resized.mrc')
            sim_padded_object = mrcfile.open(sim_padded_path, 'w+')
            # compute padding
            axis_len_odd = [ (m % 2 !=0) ^ (s % 2 != 0)   for m,s in zip(map_shape,sim_shape) ]
            padding_xyz = [ [(size_map-size_sim)//2,(size_map-size_sim)//2+1] if is_odd else [(size_map-size_sim)//2,(size_map-size_sim)//2] for size_map,size_sim,is_odd in zip(map_shape, sim_shape, axis_len_odd)]
            # Check for negative padding values and slice it
            if ( (padding_xyz[0][0]<0) | (padding_xyz[0][1]<0) | (padding_xyz[1][0]<0) | (padding_xyz[1][1]<0) | (padding_xyz[2][0]<0) | (padding_xyz[2][1]<0) ):
                for i,axis_padding in enumerate(padding_xyz):
                    new_start = 0
                    new_end = sim_data.shape[i]
                    if (axis_padding[0] <0) :
                        new_start = abs(axis_padding[0])
                        padding_xyz[i][0] = 0
                    if (axis_padding[1] < 0):
                        new_end = sim_shape[i] + axis_padding[1]
                        padding_xyz[i][1] = 0
                    if (i == 0):
                        sim_data = sim_data[new_start:new_end,:,:]
                    elif (i == 1):
                        sim_data = sim_data[:,new_start:new_end,:]
                    elif (i == 2):
                        sim_data = sim_data[:,:,new_start:new_end]

            sim_data_padded = np.pad(sim_data,padding_xyz)
            #print("Segment shape before {}".format(segment_array.shape))
            sim_padded_object.set_data(sim_data_padded)
            # Get dictionary of volumes, choose element with key 1 (corresponding to 100% recommended contour level)
            map_mask = map_data >= contourLevel
            sim_mask = sim_data_padded >= 0.1

            map_data_at_contour = map_data.copy()
            sim_data_at_contour = sim_data_padded.copy()
            map_data_at_contour[np.invert(map_mask)] = 0.0
            sim_data_at_contour[np.invert(sim_mask)] = 0.0

            overlap = getRelative_Masks_Overlap(map_mask,sim_mask)
            correlation = getCorrelation(map_data_at_contour, sim_data_at_contour)
            volume_result['vol_overlap'] = round(overlap, 2) if not np.isnan(overlap) else 0
            volume_result['vol_correlation'] = round(correlation, 2) if not np.isnan(correlation) else 0 
            print("{} Corr: {}, Overlap {}".format(os.path.basename(map_filename),volume_result['vol_correlation'],volume_result['vol_overlap'] ))

            # Close objects and release memory
            sim_map_object.close()
            sim_padded_object.close()
            del map_data_at_contour 
            del sim_data_at_contour

        except Exception as e:
            with open('error.txt', 'a') as error_file:
                error_file.write("Failed to get metrics for {} with size {}  to match {}:{}\n".format(os.path.basename(simulated_path),sim_shape,map_shape, e))
            volume_result['vol_overlap'] = None
            volume_result['vol_correlation'] = None
    else:
        with open('error.txt', 'a') as error_file:
            error_file.write("Simulated map {} does not exist:\n".format(os.path.basename(simulated_path)))
        volume_result['vol_overlap'] = None
        volume_result['vol_correlation'] = None

    map_object.close()

    return volume_result


# Fun:ction to generateSimulated samples from pdb
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

     
# Generate simulated segments using EMAN
def generateSimulatedSegments(index, df, sim_model_path):
    if not os.path.exists(sim_model_path):
        os.makedirs(sim_model_path)
    pdb_id = df.at[df.index[index], 'fitted_entries']
    map_id = df.at[df.index[index], 'id']
    res = df.at[df.index[index], 'resolution']
    pdb_filename = 'models/{}.pdb'.format(pdb_id)
    simulated_path = sim_model_path+ '/sim_{}.mrc'.format(map_id)
    chain_id = df.at[df.index[index], 'chain_id']
    generated_segment_path = os.path.join(sim_model_path, os.path.basename(simulated_path.replace('.mrc','_'+chain_id+'.mrc')))
    contour = float(df.at[df.index[index], 'contourLevel'])
    # Get pdb
    pdb_filename_chain = pdb_filename.replace('.pdb','_'+chain_id+'.pdb')
    map_filename = 'models/{}.mrc'.format(map_id)
    volume_result = {}
    volume_result['index'] = index
    tmp_filename = "/tmp/{}_cmd.txt".format(os.path.basename(generated_segment_path))
    with open(tmp_filename, 'w') as tmp_file:
        tmp_file.write("1\n1\n1\n-{}\n1\n1\n1\n".format(str(res)))
    command =   "cat {}  |  pdb2vol {} {}".format(tmp_filename, pdb_filename_chain, generated_segment_path)

    try:
        execute_command(command)
        #print(command)
    except Exception as e:
        with open('error.txt', 'a') as error_file:
            error_file.write("Error simulating map {}: {}\n".format(os.path.basename(generated_segment_path), e))

    # Get simulated map volume
    if os.path.exists(generated_segment_path):
        try:
            simulated_object = mrcfile.open(simulated_path, permissive=True)
            simulated_array = simulated_object.data
            simulated_shape = simulated_array.shape
            segment_object = mrcfile.open(generated_segment_path, 'r+', permissive=True)
            segment_array = segment_object.data
            segment_shape = segment_array.shape
            map_object = mrcfile.open(map_filename, permissive=True)
            map_array = map_object.data
            map_shape = map_array.shape
            sim_padded_path = simulated_path.replace('.mrc','_resized.mrc')
            sim_padded_object = mrcfile.open(sim_padded_path, 'r')
            sim_padded_shape = sim_padded_object.data.shape
            # compute padding ## NEED WORK BELOW
            segment_data_resized = np.zeros_like(simulated_array)
            # Data shape is in (z,y,x) so origin must reverse
            simulated_origin = list(simulated_object.header['origin'].tolist())
            segment_origin = list(segment_object.header['origin'].tolist())
            simulated_origin.reverse()
            segment_origin.reverse()
            offset = [int(abs(seg_orig-sim_orig))  for sim_orig,seg_orig in zip(simulated_origin,segment_origin)]
            len_x = int(offset[0]+segment_shape[0])
            len_y = int(offset[1]+segment_shape[1])
            len_z = int(offset[2]+segment_shape[2])
            print("{} Seg orig {}, shape {}\nSim orig {}, shape {}\nPadded shape {}\nTrying to assign segment shape {} starting from {} with len {} into array with shape {}\n".format(os.path.basename(generated_segment_path),segment_origin, segment_shape,simulated_origin, simulated_shape,sim_padded_shape, segment_shape,offset, (len_x,len_y,len_z), segment_data_resized.shape))
            
            segment_data_resized[offset[0]:len_x,offset[1]:len_y,offset[2]:len_z] = segment_array
            #print("Segment shape before {}".format(segment_array.shape))

            axis_len_odd = [ (m % 2 !=0) ^ (s % 2 != 0)   for m,s in zip(sim_padded_shape,simulated_shape) ]
            padding_xyz = [ [(size_padded-size_sim)//2,(size_padded-size_sim)//2+1] if is_odd else [(size_padded-size_sim)//2,(size_padded-size_sim)//2] for size_padded,size_sim,is_odd in zip(sim_padded_shape, simulated_shape, axis_len_odd)]
            # Check for negative padding values and slice it
            if ( (padding_xyz[0][0]<0) | (padding_xyz[0][1]<0) | (padding_xyz[1][0]<0) | (padding_xyz[1][1]<0) | (padding_xyz[2][0]<0) | (padding_xyz[2][1]<0) ):
                for i,axis_padding in enumerate(padding_xyz):
                    new_start = 0
                    new_end = segment_data_resized.shape[i]
                    if (axis_padding[0] <0) :
                        new_start = abs(axis_padding[0])
                        padding_xyz[i][0] = 0
                    if (axis_padding[1] < 0):
                        new_end = segment_data_resized.shape[i] + axis_padding[1]
                        padding_xyz[i][1] = 0
                    if (i == 0):
                        segment_data_resized = segment_data_resized[new_start:new_end,:,:]
                    elif (i == 1):
                        segment_data_resized = segment_data_resized[:,new_start:new_end,:]
                    elif (i == 2):
                        segment_data_resized = segment_data_resized[:,:,new_start:new_end]

            segment_data_padded = np.pad(segment_data_resized,padding_xyz)

            segment_object.set_data(segment_data_padded)
            segment_array = segment_object.data

            # Continuie after resizing
            map_mask = map_array >= contour
            segment_mask = segment_array >= 0.1

            map_data_at_contour = map_array.copy()
            sim_data_at_contour = segment_array.copy()
            map_data_at_contour[np.invert(map_mask)] = 0.0
            sim_data_at_contour[np.invert(segment_mask)] = 0.0

            overlap = getRelative_Masks_Overlap(segment_mask,map_mask)
            volume_result['segment_volume'] = np.sum(segment_mask)
            volume_result['segment_overlap'] = round(overlap, 2) if not np.isnan(overlap) else 0

            segment_object.close()
            map_object.close()
            simulated_object.close()
            sim_padded_object.close()
            del segment_data_padded
            del map_data_at_contour
            del sim_data_at_contour
        
        except Exception as e:
            with open('error.txt', 'a') as error_file:
                error_file.write("Seg orig {}, shape {}\nSim orig {}, shape {}\nPadded shape {}\nTrying to assign segment shape {} starting from {} with len {} into array with shape {}\nError with {}: {}\n".format(segment_origin, segment_shape,simulated_origin, simulated_shape,sim_padded_shape,
                    segment_shape,offset, (len_x,len_y,len_z), segment_data_resized.shape, os.path.basename(generated_segment_path),e))
            volume_result['segment_volume'] = -1
            volume_result['segment_overlap'] = None
    return volume_result
    




# Function to select experimental data acording to cut value for overlap
def selectExperimentalDataset(result_dir, out_txt_path, cut_overlap, cut_correlation):
    # Open result csv from volume calculation
    df = pd.read_csv('dataset_volume.csv', dtype=str)
    # Replace missing volumes (failure to generate map from pdb)
    df['vol_correlation'].replace('None', np.nan, inplace=True)
    df['vol_overlap'].replace('None', np.nan, inplace=True)
    df.dropna(inplace=True)

    df['vol_correlation'] = df['vol_correlation'].astype(float)
    df['vol_overlap'] = df['vol_overlap'].astype(float)


    # Plot hist
    #plot_hist(df[['vol_overlap','vol_correlation']], 'Distribution of maps overlap and correlation with their struture', os.path.join(result_dir,'volume_overlap_correlation.eps'))   

    selected_samples = df[ (df['vol_correlation']>=cut_correlation) & (df['vol_overlap']>=cut_overlap) ]
    numMuestras = len( selected_samples)
    with open(out_txt_path, 'a') as out:
        out.write("Number of samples over overlap cut of {} and correlation cut of {} : {}\n".format(cut_overlap, cut_correlation, numMuestras))
    
    selected_samples = selected_samples.drop(['vol_correlation','vol_overlap'], axis=1)
    # Save csv with selected samples
    selected_samples.to_csv('dataset_selected.csv', index = False)


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
        filepath = 'models/{}_nohet.pdb'.format(entryid)
        #centered_filepath = filepath.replace('.ent','_centered.ent')
        # Move pdb coords to center of mass
        #command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2procpdb_mod.py ' + filepath + ' ' + centered_filepath + ' --center 1'
        #execute_command(command)
        # Split centered pdb into chains
        #pdb_split_and_save(entryid,filepath)
        subunit_id_list, subunit_labels = pdb_get_chain_id_label_list(entryid, filepath)
        list_entry_name = [entryid for i in subunit_id_list]
        chain_df_to_merge = pd.DataFrame({columns[0]:list_entry_name,columns[2]:subunit_id_list, columns[1]:subunit_labels})
        dataset_df_to_merge = datasetSource[datasetSource[columns[0]]==entryid]
        cross_join = pd.merge(dataset_df_to_merge, chain_df_to_merge, on=columns[0], how='outer')
        dataset_result = dataset_result.append(cross_join, sort=False)
    # Reset indexes
    dataset_result.reset_index(inplace=True, drop=True)               
    return dataset_result

def parallelSimulationSegments(df, simulated_path, output_csv_name, comm, size, rank):
    # Get index list to schedule processess 
    index_list = df.index.tolist()
    print("Spawn procecess...")
    ''' 
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
                    segment_volume = d['segment_volume']
                    segment_overlap = d['segment_overlap']
                    df.loc[res_index, 'segment_volume'] = segment_volume
                    df.loc[res_index,'segment_overlap'] = segment_overlap

                except Exception as error:
                    print("Error calculating volume for simulated {}:{}".format(df.loc[res_index].id,error))
    '''                                                                                                                                                                                                                               
    if rank == 0:                                                                                                                                                                                                                     
                                                                                                                                                                                                                                          
        for i,row in df.iterrows():                                                                                                                                                                                                       
            d = generateSimulatedSegments(i, df, simulated_path)                                                                                                                                                                        
            res_index = d['index']                                                                                                                                                                                                        
            capture_ratio = d['segment_overlap']                                                                                                                                                                                              
            volume = d['segment_volume']                                                                                                                                                                                            
            df.loc[res_index, 'segment_overlap'] = str(capture_ratio)                                                                                                                                                                  
            df.loc[res_index, 'segment_volume'] = str(volume)                                                                                                                                                                
    
    df.to_csv(output_csv_name, index=False)


def main():
    global current_dir


    parser = argparse.ArgumentParser()
    parser.add_argument('--d', required=False, default=False, help='True if want to resync emdb data')
    parser.add_argument('--v', required=False, default=False, help='True if want to compare map volume with its simulated counterpart')
    parser.add_argument('--s', required=False, default=None, help='Please provide a overlap value cut to select maps to be included in final dataset', type=str)  
    parser.add_argument('--g', required=False, default=None, help='Please provide percent of simulated data to be included in final dataset')
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
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = MPI.COMM_WORLD.Get_rank()
    # Download and process data
    if int(opt.d):
        # only dowload data in rank 0
        if(rank == 0):
            get_headers(header_path)
            generate_dataframe(header_path)
            processMetadata(output_txt_path)
            downloadModels(models_path, output_txt_path)
            #downloadModelsPDB(pdb_path, models_path, output_txt_path)
            #removeNonExistingModels(models_path, output_txt_path)
            processfiles(models_path, output_txt_path)
    #Calculate volume
    if int(opt.v):
        df = pd.read_csv('./dataset_metadata_selected.csv', dtype=str)
        df_volume = df[['id', 'fitted_entries', 'resolution','contourLevel']].copy()
        # Get index list to schedule processess 
        index_list = df_volume.index.tolist()
        print("Spawn procecess...")
        with open(output_txt_path, 'a') as out:
            out.write("Number of samples to process volume: {}\n".format( len(index_list)))
            out.write("Number of avaliable threads {}\n".format(size))
         
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
                        capture_ratio = d['vol_overlap']
                        correlation = d['vol_correlation']
                        map_length = d['map_length']
                        df_volume.loc[res_index, 'map_length'] = map_length
                        df_volume.loc[res_index, 'vol_overlap'] = str(capture_ratio)
                        df_volume.loc[res_index, 'vol_correlation'] = str(correlation)
                    except ValueError as error:
                        print("Error calculating volume : {}".format(error))
        '''
        if rank == 0:
           
        for i,row in df.iterrows():
            d = simulateMapAndCompareVolume(i, df, simulated_path)
            res_index = d['index']
            capture_ratio = d['vol_overlap']
            correlation = d['vol_correlation']
            df_volume.loc[res_index, 'vol_overlap'] = str(capture_ratio)
            df_volume.loc[res_index, 'vol_correlation'] = str(correlation)
        '''    
        df_volume.to_csv('dataset_volume.csv', index=False)
    if opt.s != None:
        if rank == 0:
            cut_off_param = opt.s.split(',')            
            selectExperimentalDataset(results_path, output_txt_path, cut_overlap=float(cut_off_param[0]), cut_correlation=float(cut_off_param[1]))
        
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
        if rank == 0:
            with open(output_txt_path, 'a') as out:
                out.write("Generating segments...\n")
        dataset_cryoem = pd.read_csv('dataset_selected.csv', dtype=str)
        dataset_cryoem['contourLevel'] = dataset_cryoem['contourLevel'].astype(float)
        dataset_cryoem['resolution'] = dataset_cryoem['resolution'].astype(float)
        #dataset_synthetic = pd.read_csv('synthetic_selected.csv')

        #subunits_path = os.path.join(simulated_path, 'subunits')

        # Create output directory for simulated segmetns
        #if not(os.path.exists(subunits_path)):
        #    os.system('mkdir '+subunits_path)
        if rank == 0:        
        #    dataset_cryoem_segments = mergeMapAndChains(dataset_cryoem, ['fitted_entries','chain_label', 'chain_id'])
        #print(dataset_cryoem_segments)
        #    dataset_cryoem_segments.to_csv('dataset_tmp.csv', index=False)
            dataset_cryoem_segments = pd.read_csv('dataset_tmp.csv', dtype=str)
            #print(dataset_cryoem_segments)

#        with open(output_txt_path, 'a') as out:
#            out.write("Experimental maps merged..\n")
        #dataset_cryoem_segments = pd.read_csv('dataset_tmp.csv', dtype=str)
        parallelSimulationSegments(dataset_cryoem_segments, simulated_path, 'dataset_exp_merged.csv', comm, size, rank)
        #dataset_synthetic_merged = mergeMapAndChains(dataset_synthetic, ['entries','chain_label','path', 'chain_id', 'map_path'])
        #with open(output_txt_path, 'a') as out:
        #    out.write("Simulated maps merged..\n")
        #parallelSimulationSegments(dataset_synthetic_merged, simulated_path, 'dataset_sim_merged.csv')
                
        


         


       





       


       
        


if __name__ == '__main__':
    main()
