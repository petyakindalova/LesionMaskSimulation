#'Get binary lesion data for 'lesion' voxels
#' @param datafile dataframe with information about each subject and last column lesion mask filename.
#'  
#' @param imagedir Path to lesion masks.
#' @param voxel_IDs List of voxel IDs.
#' @return Matrix of binary data of dimension number of voxels by number of subjects.
#'
all_subj = function(datafile, imagedir, voxel_IDs, n_cores = 8){
  
  ## Reading in the data
  #patient_characteristics = datafile
  n_subjects = nrow(datafile)
  n_voxels = length(voxel_IDs)
  
  each_subj = function(i) {
    
    patient_temp = readNIfTI(paste0(imagedir, "/", datafile$file_name[i]))
    
    subj_temp = rep(0, n_voxels)
    
    subj_temp = patient_temp[as.vector(voxel_IDs)]
    
    subj_temp
  }
  
  #registerDoParallel(cores=n_cores)
  #subj_all_output = foreach(k=1:n_subjects) %dopar% each_subj(k)
  subj_all_output = matrix(rep(0, n_voxels*n_subjects), ncol=n_subjects)
  
  for (k in 1:n_subjects) {
    subj_all_output[,k] = each_subj(k)
    
    # #to be removed later
    # if(k %% 100==0) {
    #  print(k)
    #  print(mem_used())
    # }
  }
  
  subj_all_output
}
