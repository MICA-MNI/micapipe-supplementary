# micapipe-supplementary

Code used to generate figures 3 and 4 of the paper "Micapipe: a standardized processing pipeline for multiscale imaging, from connectomes to gradients"

|               Script               |                                                                                                                                                     Description                                                                                                                                                     |
|:----------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------:|
| mics_gradients.py                  | This scripts loads all the MICs connectomes of the different modalities (GD, SC, MPC and FC), and generates:  a) mean matrix plot b) Scree plot of the c) 3D scatterplot of G1, G2 and G3 d) calculates and saves a MICs to MICs-mean gradients correlation (rho) CSV file in the format subject x gradient         |
| gradients_dataset.py               | Loads all connectomes by dataset and modality, and calculates the within group mean gradient correlation. Additionally it plots the group level G1-G3 on fs5 surface                                                                                                                                                |
| gradients_datasets-correlations.py | This scripts reads each subject connectomes by modality and dataset, then calculates a correlation matrix between datasets by modality .                                                                                                                                                                            |
| func_gradients.py                  | Here are stored different functions which: a) Loads connectomes by group and modality b) Plots scree plot and 3d scatterplot of the gradients c) Calculates a dataset-subject correlation table (subject x gradient) d) Calculates a dataset-dataset correlation matrix per gradient and plots them (G1,G2, G3, G9) |
