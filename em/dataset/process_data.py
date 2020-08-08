import pandas as pd
import os
from glob import glob
from xml.etree import ElementTree
import argparse
import sys
sys.path.append('..')
import molecule
from mpi4py import MPI
from mpi4py.futures import MPICommExecutor
from Bio.PDB import PDBParser
import subprocess


# Sync headers to folder
def get_headers(header_path):
    if not os.path.exists(header_path):
        os.makedirs(header_path)
    print("Getting v1.9 headers...")
    command = 'rsync -rlpt --ignore-existing -v -z --delete --port=33444 --include "emd-*-v19.xml" "rsync.rcsb.org::emdb/structures/EMD-*/header/" '+ header_path
    os.system(command)

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

# Download candidate models from database
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
        try:
            if os.system(command_emdb) != 0:
                raise Exception('Command "%s" does not exist' % command_emdb)
        except:
            print('Command "%s" does not work' % command_emdb)
            

    
    for uri,name in zip(pdb_ftp_list, pdb_id_list):
        command_pdb = 'rsync -rlpt --ignore-existing -v -z -L --delete --port=33444 '+  uri  +' '+ models_path+'/'
        try:
            if os.system(command_pdb) != 0:
                raise Exception('Command "%s" does not exist' % command_pdb)
        except:
            print('Command "%s" does not work' % command_pdb)
            
    try:
        command = 'gunzip -f ' +models_path+'/*.gz'
        if os.system(command) != 0:
            raise Exception('Command "%s" does not exist' % command)
    except:
        print('Command "%s" does not work' % command)
    

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
   
    print("{} entries with less than 2 subunits".format(non_enough_subunit))
    print("{} candidates for dataset".format(number_of_samples))


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
    map_volume = map_object.getVolume()[1]
    # Get voxel size in Angstroms to generate simulated map 
    voxel_size = map_object.getVoxelSize()[0] #using only one value, supose its the same for all axes
    if voxel_size == 0:
        print("Map {} has voxel volume of 0, header is: \n {}".format(os.path.basename(map_filename), map_object.emMap.rawHeader))
    # Get map bounding box
    map_box = map_object.getCellDim()

    simulated_filename = os.path.join(sim_model_path, 'sim_'+os.path.basename(map_filename).replace('.map','.mrc'))
    # Generate map
    try:
        command = '/work/mzumbado/EMAN2/bin/python /work/mzumbado/EMAN2/bin/e2pdb2mrc.py -A=' + str(voxel_size)+ ' -R=' + str(res) + ' -B='+ str(int(round(map_box[2]))) +','+ str(int(round(map_box[1])))+ ','+ str(int(round(map_box[0]))) + ' --center '+  pdb_filename + ' ' + simulated_filename
        print(command)
        if os.system(command) != 0:
            raise Exception('Command "%s" does not exist' % command)
    except:
        print('Command "%s" does not work' % command)

    volume_result = dict()
    volume_result['index']=index
    volume_result['map_volume']=map_volume
    
    # Get simulated map volume
    if os.path.exists(simulated_filename):
        map_object = molecule.Molecule(simulated_filename, recommendedContour=contourLevel, cutoffRatios=[1])
        # Get dictionary of volumes, choose element with key 1 (corresponding to 100% recommended contour level)
        pdb_volume = map_object.getVolume()[1]
        volume_result['pdb_volume']=pdb_volume
    else:
        volume_result['pdb_volume']=-1

    return volume_result
    



def main():
    global current_dir


    parser = argparse.ArgumentParser()
    parser.add_argument('--d', required=True, default=False, help='True if want to resync emdb data')
    parser.add_argument('--v', required=True, default=False, help='True if want to compare map volume with its simulated counterpart')
    parser.add_argument('--header_dir', default='header', required=False,help='output directory to store emdb headers') 
    parser.add_argument('--models_dir', default='models', required=False, help='output directory to store emdb models') 
    parser.add_argument('--simulated_dir', default='simulated', required=False, help='output directory to store simulated maps') 

    opt = parser.parse_args()
    current_dir = os.getcwd()
    header_path = os.path.join(current_dir, opt.header_dir)
    models_path = os.path.join(current_dir, opt.models_dir)
    simulated_path = os.path.join(current_dir, opt.simulated_dir)
    # Download and process data
    if int(opt.d):
        get_headers(header_path)
        generate_dataframe(header_path)
        processMetadata()
        downloadModels(models_path)
        removeNonExistingModels(models_path)
        processfiles(models_path)
    #Calculate volume
    elif int(opt.v):
        df = pd.read_csv('./dataset_metadata.csv')
        print("Number of samples to process volume: ", len(df.index))
        df_volume = df[['id', 'fitted_entries', 'resolution','contourLevel']].copy()
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
                        df_volume.loc[res_index, 'map_volume'] = map_volume
                        df_volume.loc[res_index, 'pdb_volume'] = pdb_volume
                    except ValueError as error:
                        print("Error calculating volume for simulated {}: {}".format(df_volume.loc[i,'fitted_entries'],error))
                df_volume.to_csv('dataset_volume.csv', index=False)

       


       
        


if __name__ == '__main__':
    main()
