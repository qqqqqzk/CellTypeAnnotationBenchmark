args <- commandArgs(TRUE)

run_transferdata <- function(DataPath,LabelsPath,CV_RDataPath,OutputDir,ExternalTestDataPath,ExternalTestLabelsPath,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run transferdata
  Wrapper script to run transferdata on a benchmark dataset with 5-fold cross validation,
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

  ExternalTestData <- read.csv(ExternalTestDataPath,row.names = 1)
  ExternalTestLabels <- as.matrix(read.csv(ExternalTestLabelsPath))
  ExternalTestData <- t(as.matrix(ExternalTestData))
  
  #############################################################################
  #                             transferdata                                  #
  #############################################################################
  library(Seurat)
  options(future.globals.maxSize = 8000 * 1024^2)  # SCTransform
  True_Labels_transferdata <- list()
  Pred_Labels_transferdata <- list()
  External_True_Labels_transferdata <- list()
  External_Pred_Labels_transferdata <- list()
  Prep_Time_transferdata <- list()
  Testing_Time_transferdata <- list()
  ExternalTesting_Time_transferdata <- list()
  Prep_Memory_transferdata <- list()
  Testing_Memory_transferdata <- list()
  ExternalTesting_Memory_transferdata <- list()
  
  for (i in c(1:n_folds)){
    print(paste0(i, ' start'))

    # Train set
    seuratobj <- CreateSeuratObject(counts = Data[,Train_Idx[[i]]])
    # 1. add cv dataset_index -> splitObject()
    train_dataset_index <- as.data.frame(dataset_index[Train_Idx[[i]]], 
                                         row.names = colnames(seuratobj))
    colnames(train_dataset_index) <- c('dataset_index')
    seuratobj <- AddMetaData(seuratobj, metadata = train_dataset_index)

    # 2. add labels
    train_Labels <- as.data.frame(Labels[Train_Idx[[i]]], 
                                  row.names = colnames(seuratobj))
    colnames(train_Labels) <- c('Labels')
    seuratobj <- AddMetaData(seuratobj, metadata = train_Labels)
    
    # Test set
    test_seuratobj <- CreateSeuratObject(counts = Data[,Test_Idx[[i]]])
    # 1. add dataset_index
    test_dataset_index <- as.data.frame(dataset_index[Test_Idx[[i]]], 
                                        row.names = colnames(test_seuratobj))
    colnames(test_dataset_index) <- c('dataset_index')
    test_seuratobj <- AddMetaData(test_seuratobj, metadata = test_dataset_index)

    # 2. add labels
    test_Labels <- as.data.frame(Labels[Test_Idx[[i]]], 
                                 row.names = colnames(test_seuratobj))
    colnames(test_Labels) <- c('Labels')
    test_seuratobj <- AddMetaData(test_seuratobj, metadata = test_Labels)

    # External test set
    external_test_seuratobj <- CreateSeuratObject(counts = ExternalTestData)
    external_test_Labels <- as.data.frame(ExternalTestLabels, 
                                          row.names = colnames(external_test_seuratobj))
    colnames(external_test_Labels) <- c('Labels')
    external_test_seuratobj <- AddMetaData(external_test_seuratobj, metadata = external_test_Labels)
      
    Rprof(paste0(OutputDir,"/Rprof.out"), memory.profiling=TRUE)
    prep_start_time <- Sys.time()
    # Normalize + HVG
    seuratobj.list <- SplitObject(seuratobj, split.by = 'dataset_index')

    # Filter out datasets with fewer than 500 cells
    filtered_seuratobj.list <- lapply(seuratobj.list, function(seurat_obj) {
                                            if (ncol(seurat_obj) >= 500) {
                                              return(seurat_obj)
                                            } else {
                                              return(NULL)
                                            }
                                          })
    # Remove NULL entries from the list
    seuratobj.list <- Filter(Negate(is.null), filtered_seuratobj.list)
    #print('seuratobj.list')
    print(length(seuratobj.list))

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
      
    # Testing Set Preprocessing
    test_seuratobj <- NormalizeData(test_seuratobj, verbose = FALSE)
    test_seuratobj <- FindVariableFeatures(test_seuratobj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

    # Integrating new dataset
    test.anchors <- FindTransferAnchors(reference = seuratobj.integrated, query = test_seuratobj, dims = 1:30, reference.reduction = "pca")

    # External Testing Set Preprocessing
    external_test_seuratobj <- NormalizeData(external_test_seuratobj, verbose = FALSE)
    external_test_seuratobj <- FindVariableFeatures(external_test_seuratobj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

    # Integrating new dataset
    external_test.anchors <- FindTransferAnchors(reference = seuratobj.integrated, query = external_test_seuratobj, dims = 1:30, reference.reduction = "pca")
    
    prep_end_time <- Sys.time()
    Rprof(NULL)
    Prep_Time_transferdata[i] <- as.numeric(difftime(prep_end_time,prep_start_time,units = 'secs'))
    print('memory usage:')
    max_mem <- max(summaryRprof(paste0(OutputDir,"/Rprof.out"),memory="both")$by.total$mem.total)
    print(max_mem)
    Prep_Memory_transferdata[i] <- max_mem

    # Testing Set Annotation
    Rprof(paste0(OutputDir,"/Rprof.out"), memory.profiling=TRUE)
    start_time <- Sys.time()
    predictions <- TransferData(anchorset = test.anchors, refdata = seuratobj.integrated$Labels, dims = 1:30)
    #saveRDS(predictions, file='transferdata_prediction.rds')
    end_time <- Sys.time()
    Rprof(NULL)
    Testing_Time_transferdata[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    print('memory usage:')
    max_mem <- max(summaryRprof(paste0(OutputDir,"/Rprof.out"),memory="both")$by.total$mem.total)
    print(max_mem)
    Testing_Memory_transferdata[i] <- max_mem

    True_Labels_transferdata[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_transferdata[i] <- list(predictions$predicted.id)

    # External Testing Set Annotation
    Rprof(paste0(OutputDir,"/Rprof.out"), memory.profiling=TRUE)
    start_time <- Sys.time()
    external_predictions <- TransferData(anchorset = external_test.anchors, refdata = seuratobj.integrated$Labels, dims = 1:30)
    end_time <- Sys.time()
    Rprof(NULL)
    ExternalTesting_Time_transferdata[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    print('memory usage:')
    max_mem <- max(summaryRprof(paste0(OutputDir,"/Rprof.out"),memory="both")$by.total$mem.total)
    print(max_mem)
    ExternalTesting_Memory_transferdata[i] <- max_mem

    # Evaluation
    match <- external_test_seuratobj$Labels == external_predictions$predicted.id
    print(table(match))
      
    External_True_Labels_transferdata[i] <- list(ExternalTestLabels)
    External_Pred_Labels_transferdata[i] <- list(external_predictions$predicted.id)
  }
  
  True_Labels_transferdata <- as.vector(unlist(True_Labels_transferdata))
  Pred_Labels_transferdata <- as.vector(unlist(Pred_Labels_transferdata))
  External_True_Labels_transferdata <- as.vector(unlist(External_True_Labels_transferdata))
  External_Pred_Labels_transferdata <- as.vector(unlist(External_Pred_Labels_transferdata))

  Prep_Time_transferdata <- as.vector(unlist(Prep_Time_transferdata))
  Testing_Time_transferdata <- as.vector(unlist(Testing_Time_transferdata))
  ExternalTesting_Time_transferdata <- as.vector(unlist(ExternalTesting_Time_transferdata))

  Prep_Memory_transferdata <- as.vector(unlist(Prep_Memory_transferdata))
  Testing_Memory_transferdata <- as.vector(unlist(Testing_Memory_transferdata))
  ExternalTesting_Memory_transferdata <- as.vector(unlist(ExternalTesting_Memory_transferdata))
  
  write.csv(True_Labels_transferdata,paste0(OutputDir,'/transferdata_true.csv'),row.names = FALSE)
  write.csv(Pred_Labels_transferdata,paste0(OutputDir,'/transferdata_pred.csv'),row.names = FALSE)
  write.csv(External_True_Labels_transferdata,paste0(OutputDir,'/transferdata_true_external.csv'),row.names = FALSE)
  write.csv(External_Pred_Labels_transferdata,paste0(OutputDir,'/transferdata_pred_external.csv'),row.names = FALSE)
  write.csv(Prep_Time_transferdata,paste0(OutputDir,'/transferdata_prep_Time.csv'),row.names = FALSE)
  write.csv(Testing_Time_transferdata,paste0(OutputDir,'/transferdata_test_Time.csv'),row.names = FALSE)
  write.csv(ExternalTesting_Time_transferdata,paste0(OutputDir,'/transferdata_external_test_Time.csv'),row.names = FALSE)
  write.csv(Prep_Memory_transferdata,paste0(OutputDir,'/transferdata_prep_Memory.csv'),row.names = FALSE)
  write.csv(Testing_Memory_transferdata,paste0(OutputDir,'/transferdata_test_Memory.csv'),row.names = FALSE)
  write.csv(ExternalTesting_Memory_transferdata,paste0(OutputDir,'/transferdata_external_test_Memory.csv'),row.names = FALSE)
}
if (args[8] == "0") {
  run_transferdata(args[1], args[2], args[3], args[4], args[5], args[6])
} else {
  run_transferdata(args[1], args[2], args[3], args[4], args[5], args[6], args[7], as.numeric(args[8]))
}
