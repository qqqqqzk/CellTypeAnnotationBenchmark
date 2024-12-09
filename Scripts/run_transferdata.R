args <- commandArgs(TRUE)

run_transferdata <- function(DataPath,LabelsPath,CV_RDataPath,OutputDir,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run transferdata
  Wrapper script to run scmap on a benchmark dataset with 5-fold cross validation,
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
  dataset_index <- read.csv(sub('Labels', 'datasets', LabelsPath))
  Data <- Data[Cells_to_Keep,]
  Labels <- Labels[Cells_to_Keep]
  dataset_index <- dataset_index[Cells_to_Keep,]
  Data <- t(as.matrix(Data))
  #if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
  #  GenesOrder = read.csv(GeneOrderPath)
  #}
  
  #############################################################################
  #                             transferdata                                  #
  #############################################################################
  library(Seurat)
  True_Labels_transferdata <- list()
  Pred_Labels_transferdata <- list()
  Training_Time_transferdata <- list()
  Testing_Time_transferdata <- list()
  Training_Memory_transferdata <- list()
  Testing_Memory_transferdata <- list()
  
  for (i in c(1:n_folds)){

    # Train set
    seuratobj <- CreateSeuratObject(counts = Data[,Train_Idx[[i]]])
    # add dataset_index
    train_dataset_index <- as.data.frame(dataset_index[Train_Idx[[i]]], 
                               row.names = colnames(seuratobj))
    colnames(train_dataset_index) <- c('dataset_index')
    seuratobj <- AddMetaData(seuratobj, metadata = train_dataset_index)

    # add labels
    train_Labels <- as.data.frame(Labels[Train_Idx[[i]]], 
                        row.names = colnames(seuratobj))
    colnames(train_Labels) <- c('Labels')
    seuratobj <- AddMetaData(seuratobj, metadata = train_Labels)
    
    # Test set
    test_seuratobj <- CreateSeuratObject(counts = Data[,Test_Idx[[i]]])
test_seuratobj
    # add dataset_index
    test_dataset_index <- as.data.frame(dataset_index[Test_Idx[[i]]], 
                             row.names = colnames(test_seuratobj))
    colnames(test_dataset_index) <- c('dataset_index')
    test_seuratobj <- AddMetaData(test_seuratobj, metadata = test_dataset_index)

    # add labels
    test_Labels <- as.data.frame(Labels[Test_Idx[[i]]], 
                      row.names = colnames(test_seuratobj))
    colnames(test_Labels) <- c('Labels')
    test_seuratobj <- AddMetaData(test_seuratobj, metadata = test_Labels)
      
    Rprof(paste0(OutputDir,"/Rprof.out"), memory.profiling=TRUE)
    start_time <- Sys.time()
    # Normalize + HVG
    seuratobj.list <- SplitObject(seuratobj, split.by = 'dataset_index')
    for (dataset_i in 1:length(seuratobj.list)) {
        seuratobj.list[[dataset_i]] <- NormalizeData(seuratobj.list[[dataset_i]], verbose = FALSE)
        seuratobj.list[[dataset_i]] <- FindVariableFeatures(seuratobj.list[[dataset_i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    }
      
    # Integration
    reference.list <- seuratobj.list[names(seuratobj.list)]
    seuratobj.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:30)
    seuratobj.integrated <- IntegrateData(anchorset = seuratobj.anchors, dims = 1:30)
    DefaultAssay(seuratobj.integrated) <- "integrated"

    # Dimension reduction for reference-query
    seuratobj.integrated <- ScaleData(seuratobj.integrated, verbose = FALSE)
    seuratobj.integrated <- RunPCA(seuratobj.integrated, npcs = 30, verbose = FALSE)
      
    end_time <- Sys.time()
    Rprof(NULL)
    Training_Time_transferdata[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    print('memory usage:')
    max_mem <- max(summaryRprof(paste0(OutputDir,"/Rprof.out"),memory="both")$by.total$mem.total)
    print(max_mem)
    Training_Memory_transferdata[i] <- max_mem
      
    # Prediction
    Rprof(paste0(OutputDir,"/Rprof.out"), memory.profiling=TRUE)
    start_time <- Sys.time()
    # Preprocessing
    test_seuratobj <- NormalizeData(test_seuratobj, verbose = FALSE)
    test_seuratobj <- FindVariableFeatures(test_seuratobj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

    # Integrating new dataset
    test.anchors <- FindTransferAnchors(reference = seuratobj.integrated, query = test_seuratobj, dims = 1:30, reference.reduction = "pca")
    
    # Annotation
    predictions <- TransferData(anchorset = test.anchors, refdata = seuratobj.integrated$Labels, dims = 1:30)
    #saveRDS(predictions, file='transferdata_prediction.rds')
    end_time <- Sys.time()
    Rprof(NULL)
    Testing_Time_transferdata[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    print('memory usage:')
    max_mem <- max(summaryRprof(paste0(OutputDir,"/Rprof.out"),memory="both")$by.total$mem.total)
    print(max_mem)
    Testing_Memory_transferdata[i] <- max_mem

    # Evaluation
    match <- test_seuratobj$Labels == predictions$predicted.id
    print(table(match))
      
    True_Labels_transferdata[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_transferdata[i] <- list(predictions$predicted.id)
  }
  
  True_Labels_transferdata <- as.vector(unlist(True_Labels_transferdata))
  Pred_Labels_transferdata <- as.vector(unlist(Pred_Labels_transferdata))
  Training_Time_transferdata <- as.vector(unlist(Training_Time_transferdata))
  Testing_Time_transferdata <- as.vector(unlist(Testing_Time_transferdata))
  Training_Memory_transferdata <- as.vector(unlist(Training_Memory_transferdata))
  Testing_Memory_transferdata <- as.vector(unlist(Testing_Memory_transferdata))
  
  write.csv(True_Labels_transferdata,paste0(OutputDir,'/transferdata_true.csv'),row.names = FALSE)
  write.csv(Pred_Labels_transferdata,paste0(OutputDir,'/transferdata_pred.csv'),row.names = FALSE)
  write.csv(Training_Time_transferdata,paste0(OutputDir,'/transferdata_training_Time.csv'),row.names = FALSE)
  write.csv(Testing_Time_transferdata,paste0(OutputDir,'/transferdata_test_Time.csv'),row.names = FALSE)
  write.csv(Training_Memory_transferdata,paste0(OutputDir,'/transferdata_training_Memory.csv'),row.names = FALSE)
  write.csv(Testing_Memory_transferdata,paste0(OutputDir,'/transferdata_test_Memory.csv'),row.names = FALSE)
}
if (args[6] == "0") {
  run_transferdata(args[1], args[2], args[3], args[4])
} else {
  run_transferdata(args[1], args[2], args[3], args[4], args[5], as.numeric(args[6]))
}
