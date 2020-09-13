# Results folder
## Simulate lesion masks
File `LesionMaskSimulation.Rmd` contains the steps used to simulate 1000 lesion masks (simulated lesion masks saved in `data/simulated_data`). The Gaussian Random Field parameters were selected as in our manuscript [CITATION](CITATION) and the real data set used to obtain the intercept and age regression coefficients voxel-wise contained 13,680 UK Biobank individuals. 

## Run GLMs voxel-wise
File `main_RunGLM.Rmd` contains the steps used to run mass-univariate generalised linear models and to obtain either maximum likelihood estimates (method=1 and method_name="ML"), or mean bias-reduced estimates (method=2 and method_name="MeanBR"). The analysis is performed on the simulated 1000 steps as described above, but it can easily be adapted to run on your own data set. In such case, the directories pointing to the lesion masks and .dat file should be adjusted (section `Set paths` of the Rmd file) along with all the variables, such as number of categorical/continuous variables, number of cores to be used, etc. (section `Main variables`) of the Rmd file. 

The folders `ML` and `MeanBR` contain the GLM results, i.e. coeficient masks (e.g. `coef_age_nvars_1_method_2.nii.gz`) and standard error masks along with their standardised coefficient masks. 

## Figures
**TO DO** File `Plots.Rmd` contains an illustration of the spatial plots (axial slices) as in our manuscript. 
