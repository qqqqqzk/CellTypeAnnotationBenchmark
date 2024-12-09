args <- commandArgs(TRUE)

run_scPred<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run scPred
  Wrapper script to run scPred on a benchmark dataset with 5-fold cross validation,
  outputs lists of true and predicted cell labels as csv files, as well as computation time.
  
  Parameters
  ----------
  DataPath : Data file path (.csv), cells-genes matrix with cell unique barcodes 
  as row names and gene names as column names.
  LabelsPath : Cell population annotations file path (.csv).
  CV_RDataPath : Cross validation RData file path (.RData), obtained from Cross_Validation.R function.
  OutputDir : Output directory defining the path of the exported file.
  GeneOrderPath : Gene order file path (.csv) obtained from feature selection, 
  defining the genes order for each cross validation fold, default is NULL.
  NumGenes : Number of genes used in case of feature selection (integer), default is NULL.
  "
  
  Data <- read.csv(DataPath,row.names = 1)
  Labels <- as.matrix(read.csv(LabelsPath))
  load(CV_RDataPath)
  Labels <- as.vector(Labels[,col_Index])
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                scPred                                     #
  #############################################################################
  library(scPred)
  library(Seurat)
  library(magrittr)
  library(harmony)  
  True_Labels_scPred <- list()
  Pred_Labels_scPred <- list()
  Training_Time_scPred <- list()
  Testing_Time_scPred <- list()
  #Training_Memory_scPred <- list()
  #Testing_Memory_scPred <- list()
  Data = t(as.matrix(Data))
  
  for (i in c(1:2)){
  #for (i in c(1:n_folds)){
    seuratobj <- CreateSeuratObject(counts = Data[,Train_Idx[[i]]])
    train_Labels <- as.data.frame(Labels[Train_Idx[[i]]], 
                                  row.names = colnames(seuratobj))
    colnames(train_Labels) <- c('Labels')
    seuratobj <- AddMetaData(seuratobj, metadata = train_Labels)
    
    test_seuratobj <- CreateSeuratObject(counts = Data[,Test_Idx[[i]]])
    test_Labels <- as.data.frame(Labels[Test_Idx[[i]]], 
                                 row.names = colnames(test_seuratobj))
    colnames(test_Labels) <- c('Labels')
    test_seuratobj <- AddMetaData(test_seuratobj, metadata = test_Labels)
    
    # scPred Training    
    #Rprof(paste0(OutputDir,"/Rprof.out"), memory.profiling=TRUE)
    start_time <- Sys.time()
      
    seuratobj <- seuratobj %>% 
      NormalizeData() %>% 
      FindVariableFeatures() %>% 
      ScaleData() %>% 
      RunPCA() %>% 
      RunUMAP(dims = 1:30)
    seuratobj <- getFeatureSpace(seuratobj, "Labels")
    seuratobj <- trainModel(seuratobj)
      
    end_time <- Sys.time()
    #Rprof(NULL)
    Training_Time_scPred[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    #print('memory usage:')
    #max_mem <- max(summaryRprof(paste0(OutputDir,"/Rprof.out"),memory="both")$by.total$mem.total)
    #print(max_mem)
    #Training_Memory_scPred[i] <- max_mem
    
    # scPred Prediction
    #Rprof(paste0(OutputDir,"/Rprof.out"), memory.profiling=TRUE)
    start_time <- Sys.time()
      
    test_seuratobj <- NormalizeData(test_seuratobj)
    test_seuratobj <- scPredict(test_seuratobj, seuratobj)
        
    end_time <- Sys.time()
    #Rprof(NULL)
    Testing_Time_scPred[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    print('memory usage:')
    #max_mem <- max(summaryRprof(paste0(OutputDir,"/Rprof.out"),memory="both")$by.total$mem.total)
    #print(max_mem)
    #Testing_Memory_scPred[i] <- max_mem
    
    True_Labels_scPred[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_scPred[i] <- list(getPredictions(scp)$predClass)
  }
  True_Labels_scPred <- as.vector(unlist(True_Labels_scPred))
  Pred_Labels_scPred <- as.vector(unlist(Pred_Labels_scPred))
  Training_Time_scPred <- as.vector(unlist(Training_Time_scPred))
  Testing_Time_scPred <- as.vector(unlist(Testing_Time_scPred))
  #Training_Memory_scPred <- as.vector(unlist(Training_Memory_scPred))
  #Testing_Memory_scPred <- as.vector(unlist(Testing_Memory_scPred))
  
  setwd(OutputDir)
  
  write.csv(True_Labels_scPred,'scPred_True_Labels.csv',row.names = FALSE)
  write.csv(Pred_Labels_scPred,'scPred_Pred_Labels.csv',row.names = FALSE)
  write.csv(Training_Time_scPred,'scPred_Training_Time.csv',row.names = FALSE)
  write.csv(Testing_Time_scPred,'scPred_Testing_Time.csv',row.names = FALSE)
  #write.csv(Training_Memory_scPred,'scPred_Training_Memory.csv',row.names = FALSE)
  #write.csv(Testing_Memory_scPred,'scPred_Testing_Memory.csv',row.names = FALSE)
}

if (args[6] == "0") {
  run_scPred(args[1], args[2], args[3], args[4])
} else {
  run_scPred(args[1], args[2], args[3], args[4], args[5], as.numeric(args[6]))
}
