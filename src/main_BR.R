#-----------------------------
#Load libraries

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

#-----------------------------
#Load functions
source("/well/nichols/users/kindalov/FMRIB/Simulations/code_Sep2020/empir_prob_step1.R")
source("/well/nichols/users/kindalov/FMRIB/Simulations/code_Sep2020/lesion_matrix_step2.R")
source("/well/nichols/users/kindalov/FMRIB/Simulations/code_Sep2020/glm_step3.R")
source("/well/nichols/users/kindalov/FMRIB/Simulations/code_Sep2020/map_to_masks_step4.R")

#-----------------------------------------------------------------------------------------------------------
#read in arguments from command line
#!/usr/bin/env Rscript
Nargs = 7
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=Nargs) {
  stop("Wrong number of arguments.")
}

n_covs_cat = as.numeric(args[1])
n_covs_cont = as.numeric(args[2])
n_covs = n_covs_cat + n_covs_cont

n_cores = as.numeric(args[3])
empir_avail = as.numeric(args[4])

datadir = args[5]
imagedir = datadir

proj_path = args[6]
tempdir = file.path(proj_path, 'temp')
resultsdir = file.path(proj_path, 'results')

method = as.numeric(args[7])

#hard-coded variables - link function
link_fn = "probit" # options are logit and probit

if (method == 1) {
  method_name = "ML"
} else {
  if(method==2) {
    method_name = "meanBR"
  } else {
    stop("wrong method requested")
  }
}

brain_mask = readNIfTI(paste0(imagedir, "/MNI152_T1_2mm_brain_mask.nii")) #MNI 2mm brain mask

#--------------------------------------------
# Step 1: Calculate empirical lesion probability (if not available)
training_dataset = read.table(paste0(datadir, "/GLM_sample1000.dat"), header=T)

if(!empir_avail) {
  start_time = Sys.time()
  temp = empir_prob(datafile = training_dataset,
                    imagedir = datadir,
                    outputdir = tempdir,
                    voxel_IDs = TRUE)
  
  end_time = Sys.time()
  end_time - start_time
  
  print("Time for empir prob and voxel IDs")
  
  #save(list = ls(all.names = TRUE), 
  #      file=paste0(tempdir, "/step1.RData"))
}
##########################
# Step 2: Prepare binary lesion data for GLMs

if(empir_avail) {
  empir_prob_mask = readNIfTI(paste0(datadir, "/empir_prob_mask.nii")) #read in empirical probability mask
  voxel_idx = as.matrix(which(empir_prob_mask!=0), col = 1)
  
  # add the first null voxel, i.e. lesion-free
  #voxel_idx = rbind(which(empir_prob_mask == 0)[1], voxel_idx)
  print(dim(voxel_idx))
  
  write.table(voxel_idx, file=paste0(tempdir,"/voxel_IDs.dat"), sep=" ", row.names = F, col.names = F)
  voxel_IDs = read.table(paste0(tempdir, "/voxel_IDs.dat"), header=F)
  
} else {
  voxel_IDs = read.table(paste0(tempdir, "/voxel_IDs.dat"), header=F)
  
}
voxel_IDs = as.matrix(voxel_IDs, nrow = 1)

start_time = Sys.time()
lesions_subj = all_subj(datafile = training_dataset,
                        imagedir = datadir,
                        voxel_IDs = voxel_IDs)

end_time = Sys.time()
end_time - start_time

save(lesions_subj,
     file = paste0(tempdir, "/step2.RData"))

#-------------------------------------------------
# Step 3: Run GLMs
#
load(paste0(tempdir, "/step2.RData"))

#our proposed option is to parallelize in subsets of voxels
subset_size = 1000
n_subsets = ceiling(nrow(lesions_subj)/subset_size)

for(i in 1:n_subsets){
  
  load(paste0(tempdir, "/step2.RData"))
  
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

  #run GLM
  glm_results = fit_glm_fn(datafile = training_dataset,
                           lesionmat = lesions_subj_temp,
                           n_covs_cat = n_covs_cat, 
                           n_covs_cont = n_covs_cont,
                           GLMmethod = method,
                           link_fn = link_fn,
                           outputdir = tempdir,
                           subset = i, n_cores = n_cores)
  rm(glm_results)
  gc()
  
  print(i)
}

#--------------------------------------------------
# Step 4: Map results for all voxels to brain locations and save as nifti images

mapping = read.table(paste0(tempdir, "/voxel_IDs.dat"), header=F)
mapping = as.matrix(mapping, nrow = 1)
n_subsets = ceiling(length(mapping)/1000)
print("number of subsets of voxels")
print(n_subsets)

filename_subset = paste0(tempdir, "/GLM_subset_",1, "_", method_name, "_Nvars_", n_covs,"_results.RData")
load(filename_subset)
n_coefs = length(output[[1]]$parameter) #obtain number of coefficients since it may differ from n_covs

Sys.time()
mapping_fn(brain_mask = brain_mask, 
           tempdir = tempdir, 
           mapping = mapping, 
           n_subsets = n_subsets, 
           n_covs = n_covs, 
           n_coefs = n_coefs,
           GLMmethod = method, 
           outputdir = resultsdir)
Sys.time()
