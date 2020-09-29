# map results for all voxels to the brain locations and save nifti files
#'@param brain_mask Brain mask file - make sure it has the same header as the binary lesion masks, e.g. FLOAT64  
#'@param tempdir Path with all temperary files saved in step 3
#'@param mapping A dataframe with voxel locations for 'lesion' voxels - saved as part of step 1
#'@param n_subsets Number of subsets used to split lesion matrix - used in step 3; hard-coded subset size of 1000!
#'@param n_covs Number of covariates in model
#'@param n_coefs Number of coefficients estimated
#'@param GLMmethod Number indicating method, i.e. 1 for ML and 2 for MeanBR
#'@param outputdir Path for saving the resulting nifti masks
#'
#'@return Saved nifti files for coefficients, standard errors, z-scores
#'
mapping_fn = function(brain_mask, tempdir, mapping, n_subsets = NULL, n_covs, n_coefs, GLMmethod, outputdir) {
  
  image_temp = brain_mask
  image_temp[image_temp!=0] = 0
  
  ## create the fits list from all subsets of unique configurations
  if (GLMmethod == 1) {
    method = "ML"
  } else {
    method = "meanBR"
  }
  
  fits_overall = list()
  
  setwd(tempdir)
  
  ## create a list of all results from GLMs
  if(is.null(n_subsets)) {
    
    filename_GLMresults = paste0("GLM_", method, "_Nvars_", n_covs,"_results.RData")
    load(filename_GLMresults)
    
    fits_overall = output
    
    # #outliers
    # temp_img_outliers = image_temp
    # temp_img_outliers[brain_mask!=0] = sum(abs(fits_overall[[1]]$stand_pearson_resid) > 2)
    # for (i in 1:length(fits_overall)) {
    #   temp_img_outliers[mapping[i]] = sum(abs(fits_overall[[i]]$stand_pearson_resid) > 2)
    #   fits_overall[[i]]$stand_pearson_resid = NULL
    #   gc()
    #   if(i %% 100 == 0) print (i)
    # }
    
    # temp_img_outliers[1:902629]  = as.nifti(temp_img_outliers, image_temp)
    # filename_img = paste0("outliercounts_nvars_", n_covs, "_method_", GLMmethod)
    # writeNIfTI(temp_img_outliers, filename_img)
    
    
    ## add +1 for the intercept term
    n_coefs = length(output[[1]]$parameter)
    
    fits_zscore = matrix(rep(0,length(fits_overall)*(n_coefs)), ncol=(n_coefs))
    fits_estimate = matrix(rep(0, length(fits_overall)*(n_coefs)), ncol=(n_coefs))
    fits_stderr = matrix(rep(0, length(fits_overall)*(n_coefs)), ncol=(n_coefs))
    
    #if(GLMmethod == 2) {
    #  fits_zscore_LA = matrix(rep(0, length(fits_overall)*(n_coefs)), ncol=(n_coefs))
    #}
    
    fits_converged = matrix(rep(0,length(fits_overall)*(n_coefs)), ncol=(n_coefs))
    fits_scoressmall = matrix(rep(0,length(fits_overall)*(n_coefs)), ncol=(n_coefs))
    
    for (i in 1:length(fits_overall)) {
      fits_zscore[i,] = fits_overall[[i]]$zscore
      fits_estimate[i,] = fits_overall[[i]]$estimate
      fits_stderr[i,] = fits_overall[[i]]$stderr
      
      #check convergence for ML
      if(GLMmethod == 1) {
        fits_converged[i,] = fits_overall[[i]]$infinite
      }
      if(GLMmethod == 2) {
        fits_scoressmall[i,] = fits_overall[[i]]$status_br
        fits_converged[i,] = fits_overall[[i]]$converged_br
        #fits_zscore_LA[i,] = fits_overall[[i]]$LA_zscore
      }
    }
    
    names_covs = output[[1]]$parameter
    names_covs
    
    setwd(outputdir)
    
    for(i in 1:(n_coefs)){
      temp_img_vars = image_temp
      #if zero-incidence voxel included
      #temp_img_vars[brain_mask!=0] = fits_zscore[1,i]
      temp_img_vars[mapping]  = as.nifti(fits_zscore[,i], image_temp)
      filename_img = paste0("zscore_", names_covs[i], "_nvars_", n_covs, "_method_", GLMmethod)
      writeNIfTI(temp_img_vars, filename_img)
      
      temp_img_vars = image_temp
      #temp_img_vars[brain_mask!=0] = fits_estimate[1,i]
      temp_img_vars[mapping]  = as.nifti(fits_estimate[,i], image_temp)
      filename_img = paste0("coef_",names_covs[i], "_nvars_", n_covs, "_method_", GLMmethod)
      writeNIfTI(temp_img_vars, filename_img)
      
      temp_img_vars = image_temp
      #temp_img_vars[brain_mask!=0] = fits_stderr[1,i]
      temp_img_vars[mapping]  = as.nifti(fits_stderr[,i], image_temp)
      filename_img = paste0("stderr_",names_covs[i], "_nvars_", n_covs, "_method_", GLMmethod)
      writeNIfTI(temp_img_vars, filename_img)
      
      
      #if(GLMmethod == 2) {
      #  temp_img_vars = image_temp
      #  temp_img_vars[brain_mask!=0] = fits_zscore_LA[1,i]
      #  temp_img_vars[mapping]  = as.nifti(fits_zscore_LA[,i], image_temp)
      #  filename_img = paste0("LA_",names_covs[i], "_nvars_", n_covs, "_method_", GLMmethod)
      #  writeNIfTI(temp_img_vars, filename_img)
      #}
    }
    
    
    print("% infinite (ML) /converging (meanBR) voxels out of mask voxels per variable")
    apply(fits_converged,2,sum)
    #/sum(image_temp!=0)
    
    if(GLMmethod==2){
      print("% scores less than 1e-05")
      apply(fits_scoressmall,2,sum)
      #/sum(image_temp!=0)
    }
    
    
  } else {
    
    ## add +1 for the intercept term
    fits_zscore = matrix(rep(0,length(mapping)*(n_coefs)), ncol=(n_coefs))
    fits_estimate = matrix(rep(0, length(mapping)*(n_coefs)), ncol=(n_coefs))
    fits_stderr = matrix(rep(0, length(mapping)*(n_coefs)), ncol=(n_coefs))
    
    fits_converged = matrix(rep(0,length(mapping)*(n_coefs)), ncol=(n_coefs))
    fits_scoressmall = matrix(rep(0,length(mapping)*(n_coefs)), ncol=(n_coefs))
    
    #if(GLMmethod == 2) {
    #  fits_zscore_LA = matrix(rep(0, length(mapping)*(n_coefs)), ncol=(n_coefs))
    #}
    
    # #outliers
    # temp_img_outliers = image_temp

    
    for (i in 1:n_subsets) {
      filename_subset = paste0("GLM_subset_",i, "_", method, "_Nvars_", n_covs,"_results.RData")
      load(filename_subset)
      
      # if(i == 1) {
      #   temp_img_outliers[brain_mask!=0] = sum(abs(output[[1]]$stud_pearson_resid) > 2)
      # }
      
      subset_size=1000
      #determine indices
      if(i == n_subsets) {
        subset_idx = seq(subset_size*(i-1)+1, length(mapping), by = 1)
        #print(subset_idx)
      } else {
        subset_idx = seq(subset_size*(i - 1)+1, subset_size*i , by = 1)
        #print(subset_idx)
      }
      
      for (j in 1:length(subset_idx)) {
        fits_zscore[subset_idx[j],] = output[[j]]$zscore
        fits_estimate[subset_idx[j],] = output[[j]]$estimate
        fits_stderr[subset_idx[j],] = output[[j]]$stderr
        
        #check convergence for ML
        if(GLMmethod == 1) {
          fits_converged[subset_idx[j],] = output[[j]]$infinite
        }
        if(GLMmethod == 2) {
          fits_scoressmall[subset_idx[j],] = output[[j]]$status_br
          fits_converged[subset_idx[j],] = output[[j]]$converged_br
          #fits_zscore_LA[subset_idx[j],] = output[[j]]$LA_zscore
        }
        
        # #outliers
        # temp_img_outliers[mapping[subset_idx[j]]] = sum(abs(output[[j]]$stud_pearson_resid) > 2)
      }
      
      print(i)
      
    }
    
    setwd(outputdir)
    
    # #outliers file
    # temp_img_outliers[1:902629]  = as.nifti(temp_img_outliers, image_temp)
    # filename_img = paste0("outliercounts_nvars_", n_covs, "_method_", GLMmethod)
    # writeNIfTI(temp_img_outliers, filename_img)
    
    names_covs = output[[1]]$parameter
    names_covs
    
    for(i in 1:(n_coefs)){
      temp_img_vars = image_temp
      #if zero-incidence voxel included
      #temp_img_vars[brain_mask!=0] = fits_zscore[1,i]
      temp_img_vars[mapping]  = as.nifti(fits_zscore[,i], image_temp)
      filename_img = paste0("zscore_", names_covs[i], "_nvars_", n_covs, "_method_", GLMmethod)
      writeNIfTI(temp_img_vars, filename_img)
      
      temp_img_vars = image_temp
      #temp_img_vars[brain_mask!=0] = fits_estimate[1,i]
      temp_img_vars[mapping]  = as.nifti(fits_estimate[,i], image_temp)
      filename_img = paste0("coef_",names_covs[i], "_nvars_", n_covs, "_method_", GLMmethod)
      writeNIfTI(temp_img_vars, filename_img)
      
      temp_img_vars = image_temp
      #temp_img_vars[brain_mask!=0] = fits_stderr[1,i]
      temp_img_vars[mapping]  = as.nifti(fits_stderr[,i], image_temp)
      filename_img = paste0("stderr_",names_covs[i], "_nvars_", n_covs, "_method_", GLMmethod)
      writeNIfTI(temp_img_vars, filename_img)
      
      
      #if(GLMmethod == 2) {
      #  temp_img_vars = image_temp
      #  temp_img_vars[brain_mask!=0] = fits_zscore_LA[1,i]
      #  temp_img_vars[mapping]  = as.nifti(fits_zscore_LA[,i], image_temp)
      #  filename_img = paste0("LA_",names_covs[i], "_nvars_", n_covs, "_method_", GLMmethod)
      #  writeNIfTI(temp_img_vars, filename_img)
      #}
    }
    
    
    print("% infinite (ML) /converging (meanBR) voxels out of mask voxels per variable")
    print(apply(fits_converged,2,sum))
    #/sum(image_temp!=0)
    
    if(GLMmethod==2){
      print("% scores less than 1e-05")
      print(apply(fits_scoressmall,2,sum))
      #/sum(image_temp!=0)
    }
    
  }
  

  
}
