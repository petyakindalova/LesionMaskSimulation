# Lesion mask simulation and analysis 
  
## Table of contents
   * [How to cite?](#how-to-cite)
   * [Contents overview](#contents-overview)
   * [Illustrative analysis](#illustrative-analysis)

## How to cite?

See our manuscript on bioRxiv [CITATION](CITATION).

# Contents overview

The repository containts code for (1) lesion mask simulation, and (2) voxel-wise analysis of binary lesion masks utilising binary-response models (generalised linear models) and obtaining maximum likelihood estimates (ML) and mean bias-reduced estimates (MeanBR). The second step ensures scalability to biobank-scale neuroimaging data sets.  

## Illustrative analysis
Our manuscript compares three regression approaches for modelling binary lesion masks including mass-univariate probit regression modelling with either maximum likelihood estimates, or mean bias-reduced estimates, and spatial Bayesian modelling, where the regression coefficients have a conditional autoregressive model prior to account for local spatial dependence.     

We design a novel simulation framework of artificial lesion maps to compare the three alternative lesion mapping methods. The age effect on lesion probability estimated from a reference data set (13,680 individuals from the UK Biobank) is used to simulate a realistic voxel-wise distribution of lesions across age. To mimic the real features of lesion masks, we suggest matching brain lesion summaries (total lesion volume, average lesion size and lesion count) across the reference data set and the simulated data sets. Thus, we allow for a fair comparison between the modelling approaches, under a realistic simulation setting.

### Simulate lesion masks
File `LesionMaskSimulation.Rmd` in folder `results/` contains the steps used to simulate 1000 lesion masks (simulated lesion masks saved in `data/simulated_data`). The Gaussian Random Field parameters were selected as in our manuscript [CITATION](CITATION) and the real data set used to obtain the intercept and age regression coefficients voxel-wise contained 13,680 UK Biobank individuals. 

### Run GLMs voxel-wise
File `main_RunGLM.Rmd` in folder `results/` contains the steps used to run mass-univariate generalised linear models and to obtain either maximum likelihood estimates (method=1 and method_name="ML"), or mean bias-reduced estimates (method=2 and method_name="MeanBR"). The analysis is performed on the simulated 1000 steps as described above, but it can easily be adapted to run on your own data set. In such case, the directories pointing to the lesion masks and .dat file should be adjusted (section `Set paths` of the Rmd file) along with all the variables, such as number of categorical/continuous variables, number of cores to be used, etc. (section `Main variables`) of the Rmd file. 

The folders `ML` and `MeanBR` contain the GLM results, i.e. coeficient masks (e.g. `coef_age_nvars_1_method_2.nii.gz`) and standard error masks along with their standardised coefficient masks. 

