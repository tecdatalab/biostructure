#!/home/luis98/Documents/git-proyects/biostructure/Match_mrc_graph/pdb_mrc/EMAN2_py2/extlib/bin/python

# This example script will extract CTF parameters from bdb:e2ctf.parms into a CSV file readable by a spreadsheet

from EMAN2 import *

db=db_open_dict("bdb:e2ctf.parms",True)

out=file("ctf_parms.txt","w")
for k in db.keys():
	ctf=EMAN2Ctf()
	ctf.from_string(db[k][0])
	out.write( "%s,%1.3f,%1.1f\n"%(k,ctf.defocus,ctf.bfactor))
