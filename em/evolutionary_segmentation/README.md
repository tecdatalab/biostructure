## EvoSeg:  An evolutionary automatic subunit segmentation

EvoSeg is a tool for automatic subunit segmentation of EM Maps. 

**[Tutorial (Jupyter Noteboook)](https://github.com/tecdatalab/biostructure/blob/master/em/evolutionary_segmentation/Evolutionary_Segmentation.ipynb)**

**[Classifier Training (Jupyter Notebook)](https://github.com/tecdatalab/biostructure/blob/master/em/dataset/jupyter_notebook/Segmentation_classifier_train.ipynb)


## User Guide 

## Program command

### Required

- **--input_path** string 
    - EM Map in MRC Format with protein structure.
- **--output_dir** string
    - Output directory path to save output.
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
    
