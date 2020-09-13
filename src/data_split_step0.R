#'Divide the full data set in training and test datasets 
#'
#' @param datafile dataframe with information about each subject and last column lesion mask filename.
#' @param train_percent Percent of training data out of the entire data set. Default 50%.
#' @param outputdir Path for saving the resulting two data sets.
#' @param outputname Core name of resulting two data sets.
#' 
#' @return Two new data sets - training and test - are saved in the output directory provided.
#' 
data_split = function(datafile, train_percent = 0.5, outputdir, outputname) {
  
  n_subjects = nrow(datafile)
  
  #decide on number of subjects in training data set
  n_train = ceiling(train_percent * n_subjects)
  
  #RANDOMLY allocate subjects to either train or test data set
  set.seed(8)
  train_idx = sample(1:n_subjects, n_train)
  
  #save training data set
  filename_temp = paste0(outputdir, "/", outputname, "_training_ids.dat")
  write.table(sort(train_idx), file=filename_temp, sep=" ", row.names = F, col.names = F)
  #save test data set
  filename_temp = paste0(outputdir, "/", outputname, "_test_ids.dat")
  write.table((1:n_subjects)[-train_idx], file=filename_temp, sep=" ", row.names = F, col.names = F)
  
}


