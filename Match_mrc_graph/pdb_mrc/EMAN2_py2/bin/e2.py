#!/home/luis98/Documents/git-proyects/biostructure/Match_mrc_graph/pdb_mrc/EMAN2_py2/extlib/bin/python

import os

e2real=os.getenv("EMAN2DIR")+"/bin/e2_real.py"
os.execlp("ipython","ipython","-i",e2real)