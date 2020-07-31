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
import subprocess


# Get headers into folder
def get_headers(header_path):
    if not os.path.exists(header_path):
        os.makedirs(header_path)
    print("Getting v1.9 headers...")
    command = 'rsync -rlpt -v -z --delete --port=33444 --include "emd-*-v19.xml" "rsync.rcsb.org::emdb/structures/EMD-*/header/" '+ header_path
    os.system(command)


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
     
        if len(pdbe_fitted_entries)>0:
            if len(pdbe_fitted_entries)>1:
                print("Model {} has more than one fitted structures, using the first one.")
            df.loc[ df.id == id_path_dict[meta], ['fitted_entries'] ] =  pdbe_fitted_entries[0]
        if processing_method_node is not None:
            df.loc[ df.id == id_path_dict[meta], ['method'] ]  =  processing_method_node.text
        if resolution_node is not None:
            df.loc[ df.id == id_path_dict[meta], ['resolution'] ] =  resolution_node.text
        if contour_node is not None:
            df.loc[ df.id == id_path_dict[meta], ['contourLevel'] ] = contour_node.text 
        if i % 100:
            print("Processed {}".format(i))
    df.to_csv('emdb_metadata.csv', index=False)



def processMetadata():
    df = pd.read_csv('./emdb_metadata.csv')
    non_fitted= len(df.index) - len(df.dropna(subset=['fitted_entries']).index)
    non_res = len(df.index) - len(df.dropna(subset=['resolution']).index)
    outside_gap = len(df.index) - len(df[(df.resolution<4.5) & (df.resolution>10)].index)
    non_enough_subunit = len(df.index) - len(df[df.subunit_count<2].index)

    df = df.dropna(subset=['fitted_entries', 'resolution', 'subunit_count'])
    df = df[(df.resolution>=4.5) & (df.resolution<=10)]

    print("{} entries without fitted structure".format(non_fitted))
    print("{} entries without reported resolution".format(non_res))
    print("{} entries outside resolution gap".format(non_res))
    print("{} candidates for dataset".format(len(df.index)))

    df["pdb_rsync"] = df["fitted_entries"].map(lambda pdb_id: "rsync.rcsb.org::ftp_data/structures/all/pdb/pdb"+pdb_id+".ent.gz")
    df["emdb_rsync"] = df["id"].map(lambda map_id: "rsync.rcsb.org::emdb/structures/"+str.upper(map_id)+"/map/"+map_id.replace("-","_")+".map.gz")

    df.to_csv('emdb_cleaned.csv', index=False)


def downloadModels(models_path):
    if not os.path.exists(models_path):
        os.makedirs(models_path)
    df = pd.read_csv('./emdb_cleaned.csv')

    emdb_ftp_list = df["emdb_rsync"].tolist()
    pdb_ftp_list = df["pdb_rsync"].tolist()
    emdb_id_list = df["id"].tolist()
    pdb_id_list = df["fitted_entries"].tolist()

    
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
        command = 'gunzip ' +models_path+'/*.gz'
        if os.system(command) != 0:
            raise Exception('Command "%s" does not exist' % command)
    except:
        print('Command "%s" does not work' % command)
    
    print("Scanning existing models..")
    pdb_files = [os.path.basename(f) for f in glob(os.path.join(models_path,'*.ent'))]
    map_files = [os.path.basename(f) for f in glob(os.path.join(models_path,'*.map'))]
    #Get entries with missing pdb or map
    df['map_found'] = df["id"].map(lambda map_id: True if map_id.replace('-','_')+'.map' in map_files else False)
    df['pdb_found'] = df["fitted_entries"].map(lambda pdb_id: True if 'pdb'+pdb_id+'.ent' in pdb_files else False)
    
    df_map_not_found = df[df["map_found"]==False]
    df_pdb_not_found = df[df["pdb_found"]==False]

    df_map_not_found.to_csv('emdb_not_found.csv', index=False)
    df_pdb_not_found.to_csv('pdb_not_found.csv', index=False)

    indexes_to_remove = df[ (df['map_found'] == False) | (df['pdb_found'] == False) ].index
    df = df.drop(indexes_to_remove)
    df.to_csv('dataset_metadata.csv', index=False)

    number_of_samples = len( df.index )
   
    print("{} candidates for dataset".format(len(df.index)))
    print("{} maps not found".format(len(df_map_not_found.index)))
    print("{} pdbs not found".format(len(df_pdb_not_found.index)))


def processPDBfiles(models_path):
    df = pd.read_csv('./dataset_metadata.csv')



def simulateMapAndCompareVolume(index, df, sim_model_path):
    if not os.path.exists(sim_model_path):
        os.makedirs(sim_model_path)
    map_filename = df.at[df.index[index], 'map_file']
    pdb_filename = df.at[df.index[index], 'pdb_file']
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
    if int(opt.d):
        get_headers(header_path)
        generate_dataframe(header_path)
        processMetadata()
        downloadModels(models_path)
        processPDBfiles(models_path)
    elif int(opt.v):
        df = pd.read_csv('./dataset_metadata.csv')
        print("Number of samples to process volume: ", len(df.index))
        df_volume = df[['id', 'fitted_entries', 'resolution','contourLevel']].copy()
        #add filename columns for each downloaded file
        df_volume["pdb_file"] = df["fitted_entries"].map(lambda pdb_name: os.path.join(models_path,'pdb'+pdb_name+'.ent'))
        df_volume["map_file"] = df["id"].map(lambda map_name: os.path.join(models_path, map_name.replace("-","_")+'.map'))
        index_list = df_volume.index.tolist()
        print("Spawn procecess...")
        comm = MPI.COMM_WORLD
        size = comm.Get_size()
        with MPICommExecutor(comm, root=0, worker_size=size) as executor:
            if executor is not None:
                
                futures = []
                for i in index_list:
                    futures.append(executor.submit(simulateMapAndCompareVolume, i, df_volume, simulated_path))
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
                        print("error %s" % error)
                df_volume.drop(columns=['pdb_file', 'map_file'], inplace=True)
                df_volume.to_csv('dataset_volume.csv', index=False)

        


       
        


if __name__ == '__main__':
    main()
