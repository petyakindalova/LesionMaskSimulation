##########################
#libraries

#install some packages to run bias-recuced GLMs in R
library(pryr)
library(oro.nifti, lib.loc = "/well/nichols/users/kindalov/FMRIB/local_libs/")
library(enrichwith, lib.loc = "/well/nichols/users/kindalov/FMRIB/local_libs/")
library(lpSolveAPI, lib.loc = "/well/nichols/users/kindalov/FMRIB/local_libs/")
library(brglm2, lib.loc = "/well/nichols/users/kindalov/FMRIB/local_libs/")
library(devtools) 
library(waldi, lib.loc ="/well/nichols/users/kindalov/FMRIB/local_libs/")

library(foreach)
library(doParallel)
library(iterators)
#########################
#functions
#setwd("/well/nichols/users/kindalov/FMRIB/CVR_application/code_Nov2019")
source("/well/nichols/users/kindalov/FMRIB/CVR_application/code_Nov2019/data_split_step0.R")
source("/well/nichols/users/kindalov/FMRIB/CVR_application/code_Nov2019/empir_prob_step1.R")
source("/well/nichols/users/kindalov/FMRIB/CVR_application/code_Nov2019/lesion_matrix_step2.R")
source("/well/nichols/users/kindalov/FMRIB/CVR_application/code_Nov2019/glm_step3.R")
source("/well/nichols/users/kindalov/FMRIB/CVR_application/code_Nov2019/map_to_masks_step4.R")

#-----------------------------------------------------------------------------------------------------------
#read in arguments from command line
#!/usr/bin/env Rscript
Nargs = 3
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=Nargs) {
  stop("Wrong number of arguments.")
}

n_covs_cat = as.numeric(args[1])
n_covs_cont = as.numeric(args[2])
n_covs = n_covs_cat + n_covs_cont

empir_avail = 0
model = args[3]

datadir = paste0("/well/nichols/users/kindalov/FMRIB/Simulations/data/repetitions_2Mar/seed_", model)
tempdir = paste0("/well/nichols/users/kindalov/FMRIB/Simulations/meanBR/repetitions_2Mar/sample1000/seed_", model,"/temp")
resultsdir = paste0("/well/nichols/users/kindalov/FMRIB/Simulations/meanBR/repetitions_2Mar/sample1000/seed_", model, "/results")
imagedir = datadir
#-----------------------------------------------------------------------------------------------------------
##########################
#
# STEP 0
#
#split the data
#print((datadir))
#GLM_model1 =  read.table("/well/nichols/users/kindalov/FMRIB/CVR_application/Analysis_Nov2019/GLM_base/GLM_base.dat", header=T)
GLM_model1 = read.table(paste0(datadir,"/GLM_sample1000.dat"), header=T)

data_split(datafile = GLM_model1,
            train_percent = 1,
            outputdir = tempdir,
            outputname = "nvs")
 
rm(GLM_model1)
# 
##########################
#
# STEP 1
#
#find empirical probability for training data
if(!empir_avail) {
training_dataset_ids = read.table(paste0(tempdir, "/nvs_training_ids.dat"), header=F)
training_dataset =  read.table(paste0(datadir,"/GLM_sample1000.dat"),
                                header=T)[training_dataset_ids$V1,]
 
start_time = Sys.time()
temp = empir_prob(datafile = training_dataset,
                   imagedir = imagedir,
                   outputdir = tempdir,
                   voxel_IDs = TRUE)
 
end_time = Sys.time()
end_time - start_time
 
print("time for empir prob and voxel IDs")
 
save(list = ls(all.names = TRUE), 
      file=paste0(tempdir, "/step1.RData"))
} 
# rm(training_dataset_ids, training_dataset, temp)
##########################
#
# STEP 2
#
#get binary lesion data for the GLMs in next step
training_dataset_ids = read.table(paste0(tempdir, "/nvs_training_ids.dat"),
                                   header=F)
training_dataset =  read.table(paste0(datadir, "/GLM_sample1000.dat"),
                               header=T)[training_dataset_ids$V1,]
 
if(empir_avail) {
	voxel_IDs = read.table(paste0(datadir, "/voxel_IDs.dat"), header=F)
} else {
	voxel_IDs = read.table(paste0(tempdir, "/voxel_IDs.dat"), header=F)

}
voxel_IDs = as.matrix(voxel_IDs, nrow = 1)

start_time = Sys.time()
lesions_subj = all_subj(datafile = training_dataset,
                         imagedir = imagedir,
                         voxel_IDs = voxel_IDs)
 
end_time = Sys.time()
end_time - start_time

save(list = ls(all.names = TRUE),
     file = paste0(tempdir, "/step2.RData"))
 
# ##########################
#
# STEP 3
#
#run GLMs for different models and different methods
load(paste0(tempdir, "/step2.RData"))
# 
# #takes 7 hours..
# start_time = Sys.time()
# glm_results = fit_glm_fn(datafile = training_dataset,
#                          lesionmat = lesions_subj,
#                          n_covs_binary = 0, 
#                          n_covs_cont = 2,
#                          GLMmethod = 1,
#                          link_fn = "logit",
#                          outputdir = "/data/greyplover/not-backed-up/oxwasp/oxwasp16/kindalov/FMRIB/CVR/code_Rpackage/temp/")
# end_time = Sys.time()
# end_time - start_time

#one option is to parallelize in subsets of voxels
subset_size = 1000
n_subsets = ceiling(nrow(lesions_subj)/subset_size)
for(i in 1:n_subsets){
  
  load(paste0(tempdir, "/step2.RData"))
  
  #get binary lesion data for the GLMs in next step
  training_dataset_ids = read.table(paste0(tempdir, "/nvs_training_ids.dat"),
                                    header=F)
  training_dataset =  read.table(paste0(datadir, "/GLM_sample1000.dat"),
                         header=T)[training_dataset_ids$V1,]
  
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
  
  #run ML
  #glm_results = fit_glm_fn(datafile = training_dataset,
  #                         lesionmat = lesions_subj_temp,
  #                         n_covs_cat = n_covs_cat, 
  #                         n_covs_cont = n_covs_cont,
  #                         GLMmethod = 1,
  #                         link_fn = "probit",
  #                         outputdir = tempdir,
  #                         subset = i, n_cores = 8)
  #rm(glm_results)
  #gc()
  
  #run meanBR
  glm_results = fit_glm_fn(datafile = training_dataset,
                           lesionmat = lesions_subj_temp,
                           n_covs_cat = n_covs_cat, 
                           n_covs_cont = n_covs_cont,
                           GLMmethod = 2,
                           link_fn = "probit",
                           outputdir = tempdir,
                           subset = i, n_cores = 8)
  rm(glm_results)
  gc()
  
  print(i)
}

##########################
#
# STEP 4
#
print("------------------------------")
print("step 4 starts")

#map results for all voxels to the brain locations and save nifti files
brain_mask = readNIfTI(paste0(imagedir, "/MNI152_T1_2mm_brain_mask.nii.gz"))
if(empir_avail) {
	mapping = read.table(paste0(datadir, "/voxel_IDs.dat"),
                     header=F)
} else {
	mapping = read.table(paste0(tempdir, "/voxel_IDs.dat"),
                     header=F)
}
mapping = as.matrix(mapping, nrow = 1)
n_subsets = ceiling(length(mapping)/1000)
print("number of subsets of voxels")
print(n_subsets)

#print("ML step 4 starts") 
#method="ML"
#filename_subset = paste0(tempdir, "/GLM_subset_",1, "_", method, "_Nvars_", n_covs,"_results.RData")
#load(filename_subset)
#n_coefs = length(output[[1]]$parameter)

#Sys.time()
#mapping_fn_vol2(brain_mask = brain_mask, 
#                tempdir = tempdir, 
#                mapping = mapping, 
#                n_subsets = n_subsets, 
#                n_covs = n_covs,
#                n_coefs = n_coefs,
#                GLMmethod = 1, 
#                outputdir = resultsdir)
#Sys.time()

print("------------------------------")
print("meanBR step 4 starts")
method="meanBR"
filename_subset = paste0(tempdir, "/GLM_subset_",1, "_", method, "_Nvars_", n_covs,"_results.RData")
load(filename_subset)
n_coefs = length(output[[1]]$parameter)

Sys.time()
mapping_fn_vol2(brain_mask = brain_mask, 
                tempdir = tempdir, 
                mapping = mapping, 
                n_subsets = n_subsets, 
                n_covs = n_covs, 
                n_coefs = n_coefs,
                GLMmethod = 2, 
                outputdir = resultsdir)
Sys.time()
