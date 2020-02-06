#include "zernike.h"

int maxOrder = 20;
bool VERBOSE = false;

void computeDescriptors(double* in_array, int x, int y, int z){
    Grid<double> g(in_array, x, y, z);
    ZernikeDescriptor<double> zd(g, maxOrder);
    for (auto &inv : zd)
    {
        printf("%f \n", inv/10);
    }
}