PYTHON_DIR = /home/manzumbado/miniconda3/envs/visual/include/python3.6m/ 
CC = g++
CFLAGS = -std=c++11 -fpic
OBJ = _zernike.so

.PHONY : clean

$(OBJ) : 
	swig -c++ -python  -o zernike_wrap.cpp zernike.i
	$(CC) -I$(PYTHON_DIR) $(CFLAGS) -c zernike_wrap.cpp -o zernike_wrap.o 
	$(CC) $(CFLAGS) -c zernike.c -o zernike.o 
	$(CC) zernike_wrap.o zernike.o -o $(OBJ) -shared -Wl,-soname,$(OBJ)
	rm -vf zernike_wrap.o zernike.o
	cp $(OBJ) ../../utils/
	
	
	
	
clean :
	rm -vf $(OBJ)zernike_wrap.o zernike.py