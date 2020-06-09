#!EMAN2_py2/extlib/bin/python
import os
import argparse
import shutil

parser = argparse.ArgumentParser()
parser.add_argument("-if","--input_file", type=str, help="PDB file.")
parser.add_argument("-od","--output_dir", type=str, help="Out dir create directory with name of pdb.")
parser.add_argument("-r","--resolution", type=str, help="Map resolution.")
parser.add_argument("-c","--chains", help="Create mrc file for all chains in pdb", action="store_true")
parser.add_argument("-v","--verbose", help="increase output verbosity", action="store_true")

args = parser.parse_args()
name_of_pdb = os.path.splitext(args.input_file)[0]
out_file =  name_of_pdb + ".mrc"

dir = args.output_dir+"/"+name_of_pdb
if os.path.exists(dir):
    shutil.rmtree(dir)
os.makedirs(dir)
complete_file_path = dir+"/"+str(out_file)

stream = os.popen("EMAN2_py2/bin/e2pdb2mrc.py -R "+str(args.resolution)+" "+str(args.input_file)+" "+complete_file_path)
output = stream.read()
if args.verbose:
	print(output)

if args.chains:
    
    with open(args.input_file) as origin_file:
        actual_chain = ''
        for line in origin_file:
            if line[0:4]=="ATOM":
                if actual_chain == '':
                    actual_chain = line[21:22]
                elif actual_chain!= line[21:22]:
                    command = "EMAN2_py2/bin/e2pdb2mrc.py -R "+str(args.resolution)+" "+"--chains "+"'"+actual_chain+"' "+str(args.input_file)+" "+ dir+"/"+name_of_pdb+"_"+actual_chain+".mrc"
                    stream = os.popen(command)
                    output = stream.read()
                    if args.verbose:
                    	print(output)
                    actual_chain= line[21:22]
                    
        command = "EMAN2_py2/bin/e2pdb2mrc.py -R "+str(args.resolution)+" "+"--chains "+"'"+actual_chain+"' "+str(args.input_file)+" "+ dir+"/"+name_of_pdb+"_"+actual_chain+".mrc"
        stream = os.popen(command)
        output = stream.read()
        if args.verbose:
        	print(output)
print("Finish")