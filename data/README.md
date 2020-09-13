# Data folder

The content of the `data` subfolders are as follows
1. `UKB_summaries`: the mean bias-reduced intercept and age estimates from a real data set mass-univariate analysis (13,680 subjects from the UK Biobank) and the associated empirical lesion probability mask in MNI space (2mm^3 voxel size).
2. `simulated_data`: 1000 simulated lesion masks using our proposed simulation framework along with .dat file for the simulated data set (containing subject ID, age, lesion mask filename as columns) and empirical lesion probability mask for the simulated data set. `results/LesionMaskSimulation.Rmd` file used to simulate those lesion masks.
3. `temp`: temporary files used as part of running the mass-univariate generalised linear models (maximum likelihood estimates and mean bias-reduced estimates). 
