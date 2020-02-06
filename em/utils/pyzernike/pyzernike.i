%module zernike

%{
    #define SWIG_FILE_WITH_INIT
    #include "zernike.h"


    void __call_at_begining()
    {
       Binomial<double>::computePascalsTriangle(60);
    }
%}

/* Include the NumPy typemaps library */
%include "numpy.i"


%init %{
    import_array();
    __call_at_begining();
%}


%apply( double* IN_ARRAY3, int DIM1, int DIM2, int DIM3) {(double * in_array, int x, int y, int z)};

%include  "zernike.h"
%rename (computeDescriptors) computeDescriptors_func;

%inline %{

    void computeDescriptors_func(double * myarray, int x, int y, int z){
        computeDescriptors(myarray, x, y, z);
    }

%}