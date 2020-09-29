# Perform GLM and save resulting coefficients, z-scores, predictions, etc. 
#'@param datafile dataframe with information about each subject and last column lesion mask filename.
#'@param lesionfile matrix resulting from step 2
#'@param n_covs_cont Number of continuous covariates
#'@param n_covs_cat Number of categorical covariates
#'@param GLMmethod Number indicating method, i.e. 1 for ML and 2 for meanBR
#'@param link_fn Link function to be used for performing the GLM, i.e. logit or probit
#'@param outputdir Path for saving the resulting lists of GLM outputs
#'@param subset subset ID, here used to split the number of regression; hard-coded subset size of 1000 regressions per subset
#'@param n_cores number of cores to be used to parallelise the GLM implementation
#'
#'@return List of results for each voxels
#'
fit_glm_fn = function(datafile, lesionmat, n_covs_cat, n_covs_cont, GLMmethod = 2, link_fn = "logit", outputdir, 
                      subset = NULL, n_cores = 8) {
  
  n_subjects = nrow(datafile)
  #we have subject id and file name as first and last column of dataframe, i.e. (-2) below
  n_covs = ncol(datafile) - 2 
  #check number of covariates
  if(n_covs != (n_covs_cat + n_covs_cont)) stop("Something is wrong with data table or number of covariate supplied!")
  
  
  ## Determine the model formula, to be used for each voxel-specific fit
  names_data = names(datafile)
  if (n_covs_cat > 0) {covs_cat_ind = seq(2, (n_covs_cat+1), by = 1)}
  if(n_covs_cont > 0) {covs_cont_ind = (n_covs_cat + 2):(ncol(datafile)-1)}
  factor_vars = ""
  cont_vars = ""
    

  if(n_covs_cat > 0) {
    for (i in 1: length(covs_cat_ind)) {
      if (i == 1) {
        factor_vars = paste0("factor(",names_data[covs_cat_ind[i]],")")
      } else {
        factor_vars = paste0(factor_vars, " + factor(", names_data[covs_cat_ind[i]], ")")
      }
    }
  }
  
  factor_vars
  
  if(n_covs_cont > 0) {
    for (i in 1:length(covs_cont_ind)) {
      if (i == 1) {
        cont_vars = names_data[covs_cont_ind[i]]
      } else {
        cont_vars = paste0(cont_vars, " + ", names_data[covs_cont_ind[i]])
      }
     #make sure the continuous covariates are demeaned
     datafile[,names_data[covs_cont_ind[i]]] = scale(datafile[,names_data[covs_cont_ind[i]]], scale=F)
    }
  }
  
  cont_vars
  
  if (factor_vars == "") {
    all_vars = paste0(cont_vars)
  } else if(cont_vars==""){
    all_vars = paste0(factor_vars)
  } else {
    all_vars = paste0(factor_vars," + ", cont_vars)
  }
  base_formula <- paste("voxel_lesion ~", all_vars)
  print(base_formula)
  
  ## GLM function to be used in foreach 
  if (GLMmethod == 1) {
    fit_glm_chosen = function(j) {
      
      current_data = data.frame(voxel_lesion = as.factor(as.vector(lesionmat[j,1:n_subjects])), datafile)
      
      fit_ml = glm(base_formula, family = binomial(link_fn), data = current_data,
                    maxit = 100)
      
      coefs_ml = coef(fit_ml)
      coefnames = names(coefs_ml)
      
      mod = glm(base_formula, family = binomial(link_fn), data = current_data,
                 method = "detect_separation", linear_program = "dual")
      is_inf = is.infinite(mod$beta)[coefnames]
      nvar = length(coefnames) 
      
      z_ml = coef(summary(fit_ml))[coefnames, "z value"]
      coef_ml = coef(summary(fit_ml))[coefnames, "Estimate"]
      stderr_ml = coef(summary(fit_ml))[coefnames, "Std. Error"]
      
      #uncomment if you want to explore outliers
      #fitted_ml = fit_ml$fitted.values
      #pearson_resid_ml = resid(fit_ml, type='pearson')
      #stud_pearson_resid_ml = resid(fit_ml, type='pearson')/sqrt(1 - hatvalues(fit_ml))

      
      ## Give z_ml value zero if infinite estimate
      z_ml[is_inf] = 0
      
      out <- list(parameter = coefnames,
                  zscore = z_ml,
                  estimate = coef_ml,
                  stderr = stderr_ml,
                  #fitted = fitted_ml,
		  #pearson_resid = pearson_resid_ml,
                  #stud_pearson_resid = stud_pearson_resid_ml,
                  infinite = is_inf,
                  voxel = j)
      
      out
      
    }
  } else {
    if(GLMmethod == 2) {
      fit_glm_chosen = function(j) {
        
        current_data = data.frame(voxel_lesion = as.factor(as.vector(lesionmat[j,1:n_subjects])), datafile)
        
        fit_br = glm(base_formula, family = binomial(link_fn), data = current_data,
                      method = "brglmFit", type = "AS_mean", maxit = 10000, epsilon = 1e-05, slowit=0.5)
        
        coefs_br = coef(fit_br)
        coefnames = names(coefs_br)
        nvar = length(coefnames) 
        
        z_br = coef(summary(fit_br))[coefnames, "z value"]
        coef_br = coef(summary(fit_br))[coefnames, "Estimate"]
        stderr_br = coef(summary(fit_br))[coefnames, "Std. Error"]
        
        #uncomment if you want to explore outliers
        #fitted_br = fit_br$fitted.values 
        #pearson_resid_br = resid(fit_br, type='pearson')
        #stud_pearson_resid_br = resid(fit_br, type='pearson')/sqrt(1 - hatvalues(fit_br))
        
        #uncomment if you want to explore location-adjusted z-scores
        #corz_br <- try(waldi(fit_br, null = 0, what = NULL, numerical=FALSE)[coefnames], silent = TRUE)
        #if (inherits(corz_br, "try-error")) {
        #  corz_br <- rep(NA, nvar)
        #}
        
        scores_br <- na.omit(fit_br$grad)
        status_br <- all(abs(scores_br) < 1e-05)
        converged_br <- fit_br$converged      
        
        out <- list(parameter = coefnames,
                    zscore = z_br,
                    estimate = coef_br,
                    stderr = stderr_br,
                    #LA_zscore = corz_br,
                    #fitted = fitted_br,
		                #pearson_resid = pearson_resid_br,	
                    #stud_pearson_resid = stud_pearson_resid_br,
                    status_br = status_br,
                    converged_br = converged_br,
                    voxel = j)
        
        out
        
      }
    } else {
      stop("Unknown GLM method suggested!")
    }
  }
  
  output = list()
  
  # if no subset id provided
  if (is.null(subset)) {
    for ( i in 1:nrow(lesionmat)) {
      output[[i]] = fit_glm_chosen(i)
      if(i%%100==0) {
        print(mem_used())
        print(i)
      }
    }
  } else {
    #if there is only 1 voxel in the subset
    if(nrow(lesionmat)==1) {
	
        #registerDoParallel(cores=1)
	      output = list()
        output[[1]] = fit_glm_chosen(1)
        
	} else {
	  
    registerDoParallel(cores=n_cores)
    
    output = foreach(i = 1:nrow(lesionmat), .packages = 'brglm2') %dopar% fit_glm_chosen(i)
	}
  }
 
  setwd(outputdir)
  
  rm(datafile)
  rm(lesionmat)
  
  if (GLMmethod == 1) {
    method = "ML"
  } else {
    method = "meanBR"
  }
  
  if(is.null(subset)) {
    filename_GLMresults = paste0("GLM_", method, "_Nvars_", n_covs,"_results.RData")
  } else {
    filename_GLMresults = paste0("GLM_subset_",subset, "_", method, "_Nvars_", n_covs,"_results.RData")
  }
  
  
  save(output, file = filename_GLMresults)
  
  output
}
