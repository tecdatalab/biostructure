## EvoSeg:  An evolutionary automatic subunit segmentation

EvoSeg is a tool for automatic subunit segmentation of EM Maps. 

**[Tutorial (Jupyter Noteboook)](https://github.com/tecdatalab/biostructure/blob/master/em/evolutionary_segmentation/Evolutionary_Segmentation.ipynb)**

## Dependencies
 
 -  numpy>=1.19.5
 -  pandas>=1.1.5
 -  python>=3.6
 -  scikit-image>=0.17.2
 -  scikit-learn>=0.24.2
 -  scipy>=1.5.4
 -  mrcfile>=1.3.0

## Instalation
  
    pip install -i https://test.pypi.org/simple/ evoseg-tecdatalab

## Program command

### Required

- **--input_path** string 
    - EM Map in MRC Format with protein structure.
- **--output_dir** string
    - Output directory path to save output as npy array with numerical labels for each segment.
- **--level** float
    - Isosurface contour level
    
### Optional
- **--init_size** int
    - Size of initial population
- **--n_mates** int
    - Max number of matings per generation
- **--p_mates** float
    - Probability of combination 
- **--p_split** float
    - Probability of split mutation
- **--p_merge** float
    - Probability of merge mutation
- **--n_patience** int
    - Max number of generations without improvements before stopping.
- **--n_max** int
    - Max number of generations.
    
   
### Example
  
  python -m evoseg --input_path path_to_mrc_files/ --output_dir out/ --level 1.0
  
