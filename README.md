## About The Project
This code is the extension of the code written by Adrian Cierpka (https://github.com/aimat-lab/DFT_Task_Generation). It is designed to efficiently parallelise Density Functional Theory (DFT) calculations of molecules using different DFT "flavours" (basis set + functional combination) for applications in meta-learning. 

Added functionalities include:
- sampling the molecules equally from all of the data types in the dataset
- reading charges of each molecules from info.csv file and including them in DFT calculations
- generating the placeholder cat for only one flavour based on its index
- automatically appending calculations result to the rest of results for certain flavour
- excluding some flavours combinations (for example too computationally expensive ones) 

### Installation on HoreKa
1) Create a conda environment with packages included in the environment.yaml file:
```
conda env create -f environment.yaml
```

2) Install Turbomole 7.5.1 

Warning: for version 7.7.1 there is an issue with paralellisation - there are a few lines of code that create temporary directories for calculations, which is something that our code already takes care of, leading to interference that disrupts calculations. This issue might possibly also occur for other, newer Turbomole versions.

## Generating Placeholder Categories

In this context, placeholder category indicates a folder designated with one DFT flavour, containing the following files:
- sampled data in one xyz file, one molecule under another
- `task_info.json` file, containing the settings for the calculation (basis set and flavour) and the number of molecules

To execute the `create_placeholder_categories.py` file, adjust the paths in `scripts/test_script_cat` and use sbatch command to schedule it on the cluster.  
The last, optional argument `flavour_idx` allows for choosing only one flavour to be generated.

## Calculating Energies for Categories

This code reads in the data from placeholder categories and employs Turbomole to calculate energies for each category using appropriate DFT flavour.

The simplest way to run is to adjust the paths in `scripts/test_script_e` and use sbatch command to schedule it on the cluster. The code can be executed repeatedly, as it automatically detects if there are previously created results files and appends current results to them.

If more parallelization is required, each flavour can be run as a separate batch job. To achieve than, adjust the paths in `scripts/test_script_e_parallel` and use the `bash test_script_e_parallel` command to automatically schedule X jobs, where X = number of all flavours. 

## Postprocessing

Calculations that did not converge are still included in ther results files. For energies, nonconverges results are equal to 0, whereas for forces: [0, 0, 0]. It is assumed that the postprocessing is included in the ML part. 

## Roadmap

For now, the flavours are created by taking all possible combinations of basis sets and functionals included in `flavours.json` file and excluding those that are specified as "forbidden" in line 36 of `create_placeholder_categories.py`. However, an alternative approach providing more flexibility over the flavours choice could be importing them from a json file, calculated for example via using an example `create_flavours_file.py` script. It is not yet included in the code, though.


