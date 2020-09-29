#'Calculate empirical lesion probability
#' @param datafile dataframe with information about each subject and last column lesion mask filename.
#'  
#' @param imagedir Path to lesion masks.
#' @param outputdir Path for saving the empirical probability map (and voxel_IDs if voxel_IDs=TRUE).
#' @param voxel_IDs Flag indicating whether voxel indices for 'lesion voxels' would be part of the output. Default is set to FALSE.
#' 
#' @return Empirical lesion probability map for the provided dataset is saved in the output directory.
#' @return List of voxel ids with non-zero lesion incidence. 
#' 
empir_prob = function(datafile, imagedir, outputdir, voxel_IDs = FALSE){
  
  n_subjects = nrow(datafile)
  
  #load the first image
  patient_temp = readNIfTI(paste0(imagedir, "/", datafile$file_name[1]))
  
  for (i in 2:n_subjects) {
    patient_temp = patient_temp + readNIfTI(paste0(imagedir, "/", datafile$file_name[i]))
  }
  empir_prob_mask = patient_temp / n_subjects
  
  writeNIfTI(empir_prob_mask, paste0(outputdir, "/empir_prob_mask"))
  
  ####
  
  if(voxel_IDs) {
    voxel_idx = as.matrix(which(empir_prob_mask!=0), col = 1)
    print(dim(voxel_idx))
  
    # add the first null voxel, i.e. lesion-free
    #voxel_idx = rbind(which(empir_prob_mask == 0)[1], voxel_idx)
    #print(dim(voxel_idx))
    
    write.table(voxel_idx, file=paste0(outputdir,"/voxel_IDs.dat"), sep=" ", row.names = F, col.names = F)
  }

  voxel_idx
}
