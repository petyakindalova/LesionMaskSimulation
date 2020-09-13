---
title: "Run GLM and obtain meab bias-reduced (MeanBR) estimates"
author: "Petya Kindalova"
date: "September 2020"
output: 
  html_document:
    keep_md: true
    theme: united
    toc: yes
    toc_float: yes
---
# Set the ground for the simulations

## Install R Markdown 

```r
#install.packages("rmarkdown")
```



## Install required packages 

```r
#install.packages("oro.nifti")
#install.packages("pryr")
#install.packages("enrichwith")
#install.packages("lpSolveAPI")
#install.packages("brglm2")
#install.packages("detectseparation")
#install.packages("devtools")
#devtools::install_github("ikosmidis/waldi", force=T)
#install.packages("foreach")
#install.packages("doParallel")
#install.packages("iterators")
```

## Load libraries 

```r
library(oro.nifti)
library(pryr)
library(enrichwith)
library(lpSolveAPI)
library(brglm2)
library(detectseparation)
library(devtools) 
library(waldi)
library(foreach)
library(doParallel)
library(iterators)
```

## Set paths

```r
PATH_PROJ = "D:/Neuroimaging/Simulations/LesionMaskSimulation/" # Project path
#
PATH_DATA = file.path(PATH_PROJ, 'data', 'simulated_data') # Path where simulated lesion masks are stored
training_dataset = read.table(paste0(PATH_DATA, "/GLM_sample1000.dat"), header=T)
#
PATH_TEMP = file.path(PATH_PROJ, 'data', 'temp') # Path where temporary files will be stored
PATH_RESULTS = file.path(PATH_PROJ, 'results', 'MeanBR') # Path where GLM results are saved
#PATH_RESULTS = file.path(PATH_PROJ, 'results', 'ML') # change to ML if you are planning to obtain MLEs
#
PATH_SRC = file.path(PATH_PROJ, 'src') # Path where some help functions are stored
source(paste0(PATH_SRC, "/empir_prob_step1.R"))
source(paste0(PATH_SRC, "/lesion_matrix_step2.R"))
source(paste0(PATH_SRC, "/glm_step3.R"))
source(paste0(PATH_SRC, "/map_to_masks_step4.R"))
```

## Main variables 

```r
n_covs_cat = 0 # number of categorical covariates
n_covs_cont = 1 # number of continuous covariates
n_covs = n_covs_cat + n_covs_cont

empir_avail = 1 # flag: 0 (empirical probability is not available), 1 (empirical probability is available)

n_cores = 2 # number of cores for parallel GLM 
link_fn = "probit" # options are logit and probit
method = 2 # 1: ML, 2: MeanBR
method_name = "MeanBR" # "ML" or "MeanBR"

brain_mask = readNIfTI(paste0(PATH_PROJ, "/data/MNI152_T1_2mm_brain_mask.nii.gz")) #MNI 2mm brain mask
```

# Step 1: Calculate empirical lesion probability (if not available)

```r
if(!empir_avail) {
start_time = Sys.time()
temp = empir_prob(datafile = training_dataset,
                   imagedir = PATH_DATA,
                   outputdir = PATH_TEMP,
                   voxel_IDs = TRUE)
 
end_time = Sys.time()
end_time - start_time
 
print("time for empir prob and voxel IDs")
 
save(list = ls(all.names = TRUE), 
      file=paste0(tempdir, "/step1.RData"))
} 
```

# Step 2: Prepare binary lesion data for GLMs

```r
if(empir_avail) {
	empir_prob_mask = readNIfTI(paste0(PATH_DATA, "/empir_prob_mask.nii")) #read in empirical probability mask
	voxel_idx = as.matrix(which(empir_prob_mask!=0), col = 1)
    
  # add the first null voxel, i.e. lesion-free
  #voxel_idx = rbind(which(empir_prob_mask == 0)[1], voxel_idx)
  print(dim(voxel_idx))
    
  write.table(voxel_idx, file=paste0(PATH_TEMP,"/voxel_IDs.dat"), sep=" ", row.names = F, col.names = F)
  voxel_IDs = read.table(paste0(PATH_TEMP, "/voxel_IDs.dat"), header=F)
  
} else {
	voxel_IDs = read.table(paste0(PATH_TEMP, "/voxel_IDs.dat"), header=F)

}
```

```
## [1] 20300     1
```

```r
voxel_IDs = as.matrix(voxel_IDs, nrow = 1)

start_time = Sys.time()
lesions_subj = all_subj(datafile = training_dataset,
                         imagedir = PATH_DATA,
                         voxel_IDs = voxel_IDs)
 
end_time = Sys.time()
end_time - start_time
```

```
## Time difference of 2.548649 mins
```

```r
save(list = ls(all.names = TRUE),
     file = paste0(PATH_TEMP, "/step2.RData"))
```

# Step 3: Run GLMs

```r
load(paste0(PATH_TEMP, "/step2.RData"))

#our proposed option is to parallelize in subsets of voxels
subset_size = 1000
n_subsets = ceiling(nrow(lesions_subj)/subset_size)

for(i in 1:n_subsets){
  
  load(paste0(PATH_TEMP, "/step2.RData"))
  
  #determine indices
  if(i == n_subsets) {
    if(subset_size*(i-1)+1==nrow(lesions_subj)) {
	subset_idx = nrow(lesions_subj)
	} else {
	subset_idx = seq(subset_size*(i-1)+1, nrow(lesions_subj), by = 1)
	}
  } else {
    subset_idx = seq(subset_size*(i - 1)+1, subset_size*i , by = 1)
  }
  
  
  print(i) 
  lesions_subj_temp= lesions_subj[subset_idx,]
  print(dim(lesions_subj_temp))
  if(length(subset_idx)==1) {
  lesions_subj_temp = as.matrix(as.vector(lesions_subj[subset_idx,]), nrow=1)
  lesions_subj_temp = t(lesions_subj_temp)
  print(dim(lesions_subj_temp))
  }
  rm(lesions_subj)
  gc()
  #n_cores = 2
  #run GLM
  glm_results = fit_glm_fn(datafile = training_dataset,
                           lesionmat = lesions_subj_temp,
                           n_covs_cat = n_covs_cat, 
                           n_covs_cont = n_covs_cont,
                           GLMmethod = method,
                           link_fn = link_fn,
                           outputdir = PATH_TEMP,
                           subset = i, n_cores = n_cores)
  rm(glm_results)
  gc()
  
  print(i)
}
```

```
## [1] 1
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 1
## [1] 2
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 2
## [1] 3
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 3
## [1] 4
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 4
## [1] 5
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 5
## [1] 6
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 6
## [1] 7
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 7
## [1] 8
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 8
## [1] 9
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 9
## [1] 10
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 10
## [1] 11
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 11
## [1] 12
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 12
## [1] 13
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 13
## [1] 14
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 14
## [1] 15
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 15
## [1] 16
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 16
## [1] 17
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 17
## [1] 18
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 18
## [1] 19
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 19
## [1] 20
## [1] 1000 1000
## [1] "voxel_lesion ~ age"
## [1] 20
## [1] 21
## [1]  300 1000
## [1] "voxel_lesion ~ age"
## [1] 21
```

# Step 4: Map results for all voxels to brain locations and save as nifti images

```r
mapping = read.table(paste0(PATH_TEMP, "/voxel_IDs.dat"), header=F)
mapping = as.matrix(mapping, nrow = 1)
n_subsets = ceiling(length(mapping)/1000)
print("number of subsets of voxels")
```

```
## [1] "number of subsets of voxels"
```

```r
print(n_subsets)
```

```
## [1] 21
```

```r
filename_subset = paste0(PATH_TEMP, "/GLM_subset_",1, "_", method_name, "_Nvars_", n_covs,"_results.RData")
load(filename_subset)
n_coefs = length(output[[1]]$parameter) #obtain number of coefficients since it may differ from n_covs

Sys.time()
```

```
## [1] "2020-09-13 17:30:30 BST"
```

```r
mapping_fn_vol2(brain_mask = brain_mask, 
                tempdir = PATH_TEMP, 
                mapping = mapping, 
                n_subsets = n_subsets, 
                n_covs = n_covs, 
                n_coefs = n_coefs,
                GLMmethod = method, 
                outputdir = PATH_RESULTS)
```

```
## [1] 1
## [1] 2
## [1] 3
## [1] 4
## [1] 5
## [1] 6
## [1] 7
## [1] 8
## [1] 9
## [1] 10
## [1] 11
## [1] 12
## [1] 13
## [1] 14
## [1] 15
## [1] 16
## [1] 17
## [1] 18
## [1] 19
## [1] 20
## [1] 21
## [1] "% infinite/converging voxels out of mask voxels per variable"
## [1] 20300 20300
## [1] "% scores less than 1e-05"
## [1] 88 88
```

```r
Sys.time()
```

```
## [1] "2020-09-13 17:32:34 BST"
```

# Step 5: Illustration of plots - TO ADD?

