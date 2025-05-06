# Proffinity

![pic1](https://github.com/user-attachments/assets/af65266b-2f05-404b-b86e-57aa0b2943f3)


## Installation

1. Install [anaconda](https://www.anaconda.com/docs/getting-started/anaconda/install) with respect to each specific platform. Ensure that the _conda_ command is in path and can be executed directly from the command line.

2. Create a new conda working environment:
````
conda create -n proffinity python=3.9.18
````
3. Activate the proffinity conda environemnt:
````
conda activate proffinity
````
4. Install anaconda metapackage:
````
conda install anaconda
````
5. To run Proffinity, clone this github repo and install the following packages using pip install:

````
pip install ipython plotly scikit_learn seaborn voila ipywidgets==8.1.3 pycaret Bio
````
The Proffinity program is consisted of three modules, 
- Featurization Module
- ML Prediction Module
- Visualization Module

each is running through [Voila](https://voila.readthedocs.io/en/stable/) as graphic user interfaces (GUIs). 

## Featurization Module

The purpose of the featurize module is to extract machine learning features from binary complexes for the ML prediction module as well as the graph connectivity information for the visualzation module.

### Inputs:
1. a list of protein_id (_eg._ if the protein complex is "name.pdb" then the protein_id will be "name") and binding affinity data (optional), 1 row per PPI complex, as two-column csv file. If the binding affinity data is not available, the user will use "nan" value instead. The input file should be named as input_"inputname".csv
3. complex structure data file in pdb format with filename that match the protein_id in the input list from 1. The PDB files of complex structures should be inside the 'model' folder (note: do not change the folder name "model").   

### Outputs:
1. the extracted features will be saved to the file with name ppi_index_extract_"inputname".csv where "inputname" is the name of input csv file.<sub></sub>
2. the graph connectivity file for the complex will be generated to the "raw_graphv2" folder (_note_: do not change the folder name "raw_graphv2").

### To run:
Inside the Proffinity directory,

````
voila featurizev2.ipynb
````

Module operation steps: 
![image](https://github.com/user-attachments/assets/d02e995c-9f38-47fa-9532-065c6e0b687f)

The user has the opportunity to visualize and further filter the extracted feature when the feature extraction process is completed:
![image](https://github.com/user-attachments/assets/60490145-fda4-400c-b2a7-9b7804b69e5a)

## ML Prediction Module

The purpose of the ML prediction module is to predict binding affinity (K<sub>D</sub>) of the input complexes based on the extracted features from the featurization module.

We recommended running this module with a few complexes with known (K<sub>D</sub>) using different ML models to identify the best model that optimize the validation performance. A rule-of-thumb is to aim at obtaining R<sup>2</sup> and/or RMSE values better than the provided background performance based on the SKEMPIv2 datasets. If the user cannot identify an optimal model, we recommend generating cutomized models by incoporating user provided data using the ml_regressor_train.ipynb script. 

### Inputs:
1. a pre-trained extra tree regressor model in the .plk format. Please unzip the three model files start with "saved_model" namestem. The model has been pre-trained and can be selected directly from the GUI. (Note: We have provided a seperate ml_regressor_train.ipynb for user who want to re-train their cutomized model using their own data).
2. the extracted features from the featurization module (must be named in the format of ppi_index_extract_"inputrname".csv). The csv file can also be selected directly from the GUI.

### Output:

We divided the output data as validation set (K<sub>D</sub> data provided) and/or test set (K<sub>D</sub> data not provided).

1. for validation set:
   - correlation (R<sup>2</sup>) between predicted and experimental K<sub>D</sub> values.
   - bar graph compared experimental and predicted K<sub>D</sub> (RMSE).
   - a csv output file pred_ppi_index_extract_'inputname'_by_'modelname'_vset.csv    

2. for test set:
   - histogram distribution of predicted K<sub>D</sub> values.
   - a csv output file pred_ppi_index_extract_'inputname'_by_'modelname'_vset.csv

### To run:  
Inside the Proffinity directory,

````
voila ml_regressorv2.ipynb
````

Module operation steps:

1. select and load the desired ML model:
![image](https://github.com/user-attachments/assets/a9597fe2-c61e-4f99-ba0a-1902a2ade10b)

2. load input features:
![image](https://github.com/user-attachments/assets/8ab62342-6680-463f-8ce3-e3093985e6b7)

3. Analyze the prediction outputs:
![image](https://github.com/user-attachments/assets/b5b6556e-f26a-4480-8ca0-178b70bc1017)


## Optional: ML model retraining with customize dataset

### Input preparation:
1. the extracted features (ppi_index_extract_'username'.csv) from the featurize module.
2. open the "train" folder and determine the reference features (.csv) to be combined with user's customize feature.
3. remove the header from "username.csv" and combine user's features with one of the given reference features into one single .csv file.
4. rename the combined features file.

### Output:

The script will generate two outputs:

1. a customized ml regressor model file (.plk).
2. feature importance file (.csv)

### To run:
Inside the "Training" directory,

````
jupyter notebook ml_regresso_train.ipynb
````

## Visualization Module

The purpose of the visualization module is to identify key residues and interactions that contribute to complexes' binding activity (or inactivity).

### Inputs:
1. feature importance (FI) from a pre-trained extra tree regressor model. The FI is used as weight coefficent for subsequent scoring.
2. a target complex structure in the "model" folder and associated graph connectivity data in the "raw_graphv2" folder.

### Outputs:
Interative structure model that highlight key residues and interactions based on user defined criteria. 

### To run:  <sub></sub>
Inside the Proffinity directory,

````
voila visualizationv2.ipynb
````
Module operation steps:
1. load model FI and target structures:
![image](https://github.com/user-attachments/assets/657a28b5-e241-4bd8-ad12-bd598bc23c1b)


2. visualize key residues and interactions:
![image](https://github.com/user-attachments/assets/884ed8cc-8585-492b-bbfa-b6811ee75992)




