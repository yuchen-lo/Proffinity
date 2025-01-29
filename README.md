# Proffinity

![image](https://github.com/user-attachments/assets/f9a9bb6f-e401-43d1-a70a-54fd70a7fd1c)

## Installation

To run Proffinity, clone this github repo and install the following packages using pip install:
````
pip install ipython ipywidgets numpy pandas plotly pycaret scikit_learn matplotlib seaborn voila
````
The Proffinity program is consisted of three modules, featurizev2.ipynb, ml_regressorv2.ipynb, and visualizationv2.ipynb, each is run using voila as user-interface. 

Inside the Proffinity directory,

````
voila featurizev2.ipynb
````

## Featurize module

The purpose of the featurize module is the extract machine learning features from binary complexes for the ML prediction module as well as the graph connectivity information for the visualzation module.

Inputs:
1. a list of protein_id and binding affinity data (optional), 1 row per PPI complex, as two-column csv file. If the binding affinity data is not available, the user will use 'nan' value for the entry. The input file should be named as
   input_"username".csv
3. complex structure data file in pdb format with filename that match the protein_id in the input list from 1. The structure should be inside the 'model' folder (note: do not change the folder name as the program will search for it).   

Outputs:
1. the extracted features will be saved to the file with name ppi_index_extract_"username".csv where "username" is defmined by the input csv file.
2. the graph connectivity file for the complex will be generated to the "raw_graphv2" folder (note: do not change the folder name as the program will search for it).

Module peration example: 
![image](https://github.com/user-attachments/assets/d02e995c-9f38-47fa-9532-065c6e0b687f)






