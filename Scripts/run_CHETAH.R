args <- commandArgs(TRUE)

run_CHETAH<-function(DataPath,LabelsPath,CV_RDataPath,OutputDir,ExternalTestDataPath,ExternalTestLabelsPath,GeneOrderPath = NULL,NumGenes = NULL){
  "
  run CHETAH
  Wrapper script to run CHETAH on a benchmark dataset with 5-fold cross validation,
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

  ExternalTestData <- read.csv(ExternalTestDataPath,row.names = 1)
  ExternalTestLabels <- as.matrix(read.csv(ExternalTestLabelsPath))

  if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
    GenesOrder = read.csv(GeneOrderPath)
  }
  
  #############################################################################
  #                                CHETAH                                     #
  #############################################################################
  library(CHETAH)
  library(scuttle)
  library(Matrix)
  library(Rtsne)
  True_Labels_CHETAH <- list()
  Pred_Labels_CHETAH <- list()
  ExternalTestPred_Labels_CHETAH <- list()
  Total_Time_CHETAH <- list()
  Total_Memory_CHETAH <- list()
  Data = t(as.matrix(Data))
  ExternalTestData = t(as.matrix(ExternalTestData))
  
  for (i in c(1:n_folds)){
    print(paste0(i, ' start'))
    if(!is.null(GeneOrderPath) & !is.null (NumGenes)){
      sce <- SingleCellExperiment(list(counts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Train_Idx[[i]]]), 
                                  colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))      
      sce_test <- SingleCellExperiment(list(counts = Data[as.vector(GenesOrder[c(1:NumGenes),i])+1,Test_Idx[[i]]]), 
                                       colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      sce_test_external <- SingleCellExperiment(list(counts = ExternalTestData[as.vector(GenesOrder[c(1:NumGenes),i])+1,]), 
                                                colData = data.frame(cell_type1 = ExternalTestLabels))

      Rprof(paste0(OutputDir,"/Rprof.out"), memory.profiling=TRUE)
        
      prep_start_time <- Sys.time()
      sce <- logNormCounts(sce)
      counts_sce <- logcounts(sce)
      celltypes_sce <- sce$cell_type1

      sce_test <- logNormCounts(sce_test)
      pca_data <- prcomp(t(logcounts(sce_test)), rank=30)
      tsne_data <- Rtsne(pca_data$x[,1:30], pca = FALSE)
      reducedDims(sce_test) <- list(PCA=pca_data$x, TSNE=tsne_data$Y)
      counts_test <- Matrix(logcounts(sce_test), sparse = TRUE)
      tsne_test <- reducedDim(sce_test, "TSNE")[,1:2]

      sce_test_externl <- logNormCounts(sce_test_externl)
      pca_data_external <- prcomp(t(logcounts(sce_test_externl)), rank=30)
      tsne_data_external <- Rtsne(pca_data_external$x[,1:30], pca = FALSE)
      reducedDims(sce_test_externl) <- list(PCA=pca_data_external$x, TSNE=tsne_data_external$Y)
      counts_test_externl <- Matrix(logcounts(sce_test_externl), sparse = TRUE)
      tsne_test_externl <- reducedDim(sce_test_externl, "TSNE")[,1:2]
      
      ## For the reference we define a "counts" assay and "celltypes" metadata
      sce <- SingleCellExperiment(assays = list(counts = counts_sce),
                                     colData = DataFrame(celltypes = celltypes_sce))

      ## For the input we define a "counts" assay and "TSNE" reduced dimensions
      sce_test <- SingleCellExperiment(assays = list(counts = counts_test),
                                       reducedDims = SimpleList(TSNE = tsne_test))
        
      sce_test_external <- SingleCellExperiment(assays = list(counts = counts_test_external),
                                                reducedDims = SimpleList(TSNE = tsne_test_external))
      prep_end_time <- Sys.time()

      start_time <- Sys.time()
      sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce)
      sce_test_external <- CHETAHclassifier(input = sce_test_external, ref_cells = sce)
      end_time <- Sys.time()
      Rprof(NULL)
    }
    else{
      sce <- SingleCellExperiment(list(counts = Data[,Train_Idx[[i]]]), 
                                  colData = data.frame(cell_type1 = Labels[Train_Idx[[i]]]))      
      sce_test <- SingleCellExperiment(list(counts = Data[,Test_Idx[[i]]]), 
                                       colData = data.frame(cell_type1 = Labels[Test_Idx[[i]]]))
      sce_test_external <- SingleCellExperiment(list(counts = ExternalTestData), 
                                       colData = data.frame(cell_type1 = ExternalTestLabels))

      Rprof(paste0(OutputDir,"/Rprof.out"), memory.profiling=TRUE)
      prep_start_time <- Sys.time()
        
      sce <- logNormCounts(sce)
      counts_sce <- logcounts(sce)
      celltypes_sce <- sce$cell_type1

      sce_test <- logNormCounts(sce_test)
      pca_data <- prcomp(t(logcounts(sce_test)), rank=30)
      tsne_data <- Rtsne(pca_data$x[,1:30], pca = FALSE)
      reducedDims(sce_test) <- list(PCA=pca_data$x, TSNE=tsne_data$Y)
      counts_test <- Matrix(logcounts(sce_test), sparse = TRUE)
      tsne_test <- reducedDim(sce_test, "TSNE")[,1:2]

      sce_test_externl <- logNormCounts(sce_test_externl)
      pca_data_external <- prcomp(t(logcounts(sce_test_externl)), rank=30)
      tsne_data_external <- Rtsne(pca_data_external$x[,1:30], pca = FALSE)
      reducedDims(sce_test_externl) <- list(PCA=pca_data_external$x, TSNE=tsne_data_external$Y)
      counts_test_externl <- Matrix(logcounts(sce_test_externl), sparse = TRUE)
      tsne_test_externl <- reducedDim(sce_test_externl, "TSNE")[,1:2]
      
      ## For the reference we define a "counts" assay and "celltypes" metadata
      sce <- SingleCellExperiment(assays = list(counts = counts_sce),
                                     colData = DataFrame(celltypes = celltypes_sce))

      ## For the input we define a "counts" assay and "TSNE" reduced dimensions
      sce_test <- SingleCellExperiment(assays = list(counts = counts_test),
                                       reducedDims = SimpleList(TSNE = tsne_test))

      sce_test_external <- SingleCellExperiment(assays = list(counts = counts_test_external),
                                                reducedDims = SimpleList(TSNE = tsne_test_external))
      prep_end_time <- Sys.time()
        
      start_time <- Sys.time()
      sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce)      
      sce_test_external <- CHETAHclassifier(input = sce_test_external, ref_cells = sce)
      end_time <- Sys.time()
      Rprof(NULL)
    }
    
    Prep_Time_CHETAH[i] <- as.numeric(difftime(prep_end_time,prep_start_time,units = 'secs'))
    Total_Time_CHETAH[i] <- as.numeric(difftime(end_time,start_time,units = 'secs'))
    print('memory usage:')
    max_mem <- max(summaryRprof(paste0(OutputDir,"/Rprof.out"),memory="both")$by.total$mem.total)
    print(max_mem)
    Total_Memory_CHETAH[i] <- max_mem

    
    True_Labels_CHETAH[i] <- list(Labels[Test_Idx[[i]]])
    Pred_Labels_CHETAH[i] <- list(sce_test$celltype_CHETAH)
    ExternalTestPred_Labels_CHETAH[i] <- list(sce_test_external$celltype_CHETAH)
  }
  True_Labels_CHETAH <- as.vector(unlist(True_Labels_CHETAH))
  Pred_Labels_CHETAH <- as.vector(unlist(Pred_Labels_CHETAH))
  ExternalTestPred_Labels_CHETAH <- as.vector(unlist(ExternalTestPred_Labels_CHETAH))
  Total_Time_CHETAH <- as.vector(unlist(Total_Time_CHETAH))
  Total_Memory_CHETAH <- as.vector(unlist(Total_Memory_CHETAH))
  write.csv(True_Labels_CHETAH,paste0(OutputDir,'/CHETAH_true.csv'),row.names = FALSE)
  write.csv(Pred_Labels_CHETAH,paste0(OutputDir,'/CHETAH_pred.csv'),row.names = FALSE)
  write.csv(ExternalTestPred_Labels_CHETAH,paste0(OutputDir,'/CHETAH_pred_external.csv'),row.names = FALSE)
  write.csv(Total_Time_CHETAH,paste0(OutputDir,'/CHETAH_total_time.csv'),row.names = FALSE)
  write.csv(Prep_Time_CHETAH,paste0(OutputDir,'/CHETAH_prep_time.csv'),row.names = FALSE)
  write.csv(Total_Memory_CHETAH,paste0(OutputDir,'/CHETAH_total_memory.csv'),row.names = FALSE)
}

if (args[6] == "0") {
  run_CHETAH(args[1], args[2], args[3], args[4], args[5])
} else {
  run_CHETAH(args[1], args[2], args[3], args[4], args[5], args[6], as.numeric(args[7]))
}
