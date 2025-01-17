args <- commandArgs(TRUE)

Cross_Validation <- function(LabelsPath, col_Index = 1, OutputDir){
  "
  Cross_Validation
  Function returns train and test indices for 5 folds stratified across unique cell populations,
  also filter out cell populations with less than 10 cells.
  It return a 'CV_folds.RData' file which then used as input to classifiers wrappers.

  Parameters
  ----------
  LabelsPath : Cell population annotations file path (.csv).
  col_Index : column index (integer) defining which level of annotation to use,
  in case of multiple cell type annotations (default is 1)
  OutputDir : Output directory defining the path of the exported file.
  "

  Labels <- as.matrix(read.csv(LabelsPath))
  Labels <- as.vector(Labels[,col_Index])

  Removed_classes <- !(table(Labels) > 10)  # !(table(Labels) > 10)
  Cells_to_Keep <- !(is.element(Labels,names(Removed_classes)[Removed_classes]))
  Labels <- Labels[Cells_to_Keep]

  # Getting training and testing Folds
  dataset_index <- read.csv(sub('Labels', 'datasets', LabelsPath))
  dataset_index <- as.data.frame(dataset_index[Cells_to_Keep,])
  #dataset_index <- read.csv(paste(data_path, '/datasets.csv', sep=''))
  Train_Idx <- list()
  Test_Idx <- list()
  # unique value of dataset_index
  unival <- unique(dataset_index$dataset_index)
  n_folds <- length(unival)
  for (i in 1:n_folds){
      # train: exclude i dataset
      Train_Idx[[i]] <- which(dataset_index$dataset_index != unival[i])
      # test: include i dataset only
      Test_Idx[[i]] <- which(dataset_index$dataset_index == unival[i])
  }
  save(n_folds,Train_Idx,Test_Idx,col_Index,Cells_to_Keep,file = paste0(OutputDir, '/CV_folds.RData'))
}

Cross_Validation(args[1], as.numeric(args[2]), args[3])
