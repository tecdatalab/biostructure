#ifndef ZERNIKE_H
#define ZERNIKE_H

#include "../map2zernike/Grid.h"
#include "../map2zernike/ZernikeDescriptor.h"
#include "../map2zernike/util.h"
#include <stdexcept>


void computeDescriptors(double* in_array, int x, int y, int z, double** out_array, int* x_out);

#endif