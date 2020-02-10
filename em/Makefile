
ZERNIKE_WRAPPER_DIR=utils/pyzernike/ 
OUT=_zernike.so

.PHONY : clean

all : 
	cd $(ZERNIKE_WRAPPER_DIR) && make 
	
	
clean :
	cd $(ZERNIKE_WRAPPER_DIR) && make clean
	rm -vf utils/_zernike.so 