import numpy as np
import mrcfile
import sys,os
from scipy import ndimage

def mrc_util(filename,sh):
	newsitus=filename.split(".")[0]+"_new.mrc"
	situs=filename
	sh=float(sh)
	if os.path.exists(situs):
	    with mrcfile.open(situs) as mrc:
	        nx,ny,nz,nxs,nys,nzs,mx,my,mz = mrc.header.nx,mrc.header.ny,mrc.header.nz,mrc.header.nxstart,mrc.header.nystart,mrc.header.nzstart,mrc.header.mx,mrc.header.my,mrc.header.mz
	        # nxs,nys,nzs = mrc.header.nxstart,mrc.header.nystart,mrc.header.nzstart
	        if nxs != 0 or nys !=0 or nzs !=0:
	        	print("nxs, nys and nzs are not symmetric... Skipping..")
	        	return 


	        orig=mrc.header.origin
	        data = np.array(mrc.data,dtype="float32")
	        # data=np.swapaxes(data,0,2)
	        # data_new=data
	        data_new=np.zeros((nx,ny,nz))
	        data_new = ndimage.interpolation.zoom(data, (sh,sh,sh),order=3)
	        # data_new=np.swapaxes(data_new,0,2)
	        mrc_new = mrcfile.new(newsitus,data=data_new, overwrite=True)
	        vsize=mrc_new.voxel_size
	        vsize.flags.writeable = True
	        vsize.x=1.0
	        vsize.y=1.0
	        vsize.z=1.0
	        mrc_new.voxel_size=vsize
	        mrc_new.update_header_from_data()   
	        mrc_new.header.nxstart=nxs*sh
	        mrc_new.header.nystart=nys*sh
	        mrc_new.header.nzstart=nzs*sh
	        mrc_new.header.origin=orig
	        mrc_new.update_header_stats()
	        #print(mrc_new.data.shape)
	        #print(mrc_new.voxel_size)
	        mrc.print_header()
	        mrc_new.print_header()
	        mrc_new.close()
	        del data
	        del data_new

def main():
	mrc_util(sys.argv[1],sys.argv[2])


main()
