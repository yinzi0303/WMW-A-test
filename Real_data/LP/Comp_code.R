

# *** TMM normalization:

TMMnormalization <- function(countTable){
  ## TMM normalization based on edgeR package:
  require("edgeR")
  
  nf=calcNormFactors(countTable ,method= "TMM")
  nf= colSums(countTable)*nf
  scalingFactors = nf/mean(nf)
  countTableTMM <- t(t(countTable)/scalingFactors)
  
  return(countTableTMM)
  
}

# ******************************************************************************************* #

# *** Running ROTS:

Run_ROTS = function(rawData, gr1, gr2, B, K , FDR){
  
  require(ROTS)
  
  # rawData is the input raw count table
  # gr1 is the number of cells in condition 1
  # gr2 is the number of cells in condition 2
  # B is the number of bootstraps
  # k is the number of top ranked genes 
  # FDR is the fdr threshold for the final detections
  
  # First the genes with 0 value for all the cells were filtered out:
  
  filtered = apply(rawData,1, function(x) {if(all(x == 0)) return (FALSE) else return(TRUE)})
  FilteredData = rawData[filtered,]
  
  # TMM normalization over the raw count matrix:
  
  TMM_Filtered_Data = TMMnormalization(FilteredData)
  
  # Running ROTS over the filtered and normalized data:
  
  ROTS_filtered_TMM = ROTS(data = TMM_Filtered_Data, groups = c(rep(0,as.numeric(gr1)),rep(1,as.numeric(gr2))) , B = 1000, K = 6000 )
  
  # Rows/genes with FDR smaller than the desired threshold:
  
  # ROTS_results = summary(ROTS_filtered_TMM , fdr = 0.05)
  
  # save(ROTS_filtered_TMM , ROTS_results , file = "ROTS_results.RData")
  pVals=ROTS_filtered_TMM[["FDR"]]
  
  return(pVals)
  
}

# ******************************************************************************************* #


Run_DESeq = function(rawData , gs){
  
  require("DESeq2")
  
  dds <- DESeqDataSetFromMatrix(rawData, DataFrame(gs), ~ gs)
  ## factor levels were dropped which had no samples
  dds2 <- DESeq(dds)
  ## estimating size factors
  ## estimating dispersions
  ## gene-wise dispersion estimates
  ## mean-dispersion relationship
  ## final dispersion estimates
  ## fitting model and testing
  ## -- replacing outliers and refitting for 4 genes
  ## -- DESeq argument 'minReplicatesForReplace' = 7 
  ## -- original counts are preserved in counts(dds)
  ## estimating dispersions
  ## fitting model and testing
  res <-  results(dds2, contrast=c("gs","0","1"),alpha = sig_score)
  pVals <- res[,6]#padj
  pVals[which(is.na(pVals))]=1
  # names(pVals) <- rownames(res)
  # sigDE <- which(pVals < 0.05)
  
  return(pVals)
  
}



# ******************************************************************************************* #


Run_scde = function(rawData, gs){
  require(flexmix)
  require(scde)
  cnts <- apply(rawData,2,
                function(x) {
                  storage.mode(x) <- 'integer'
                  return(x)
                }
  )
  names(gs) <- 1:length(gs)
  colnames(cnts) <- 1:length(gs)
  o.ifm <- scde::scde.error.models(
    counts = cnts,
    groups = gs,
    n.cores = 1,
    threshold.segmentation = TRUE,
    save.crossfit.plots = FALSE,
    save.model.plots = FALSE,
    verbose = 0,
    min.size.entries = 2
  )
  priors <- scde::scde.expression.prior(
    models = o.ifm,
    counts = cnts,
    length.out = 400,
    show.plot = FALSE
  )
  resSCDE <- scde::scde.expression.difference(
    o.ifm,
    cnts,
    priors,
    groups = gs,
    n.randomizations = 100,
    n.cores = 1,
    verbose = 0
  )
  # Convert Z-scores into 2-tailed p-values
  pVals <- pnorm(abs(resSCDE$cZ), lower.tail = FALSE) * 2
  pVals <- p.adjust(pVals, method = "fdr")
  pVals[which(is.na(pVals))]=1
  
  return(pVals)
  # sigDE <- names(pVals)[pVals < 0.05]
  # SCDE_pVals=pVals
  # DE_Quality_rate(sigDE)
  # DE_Quality_AUC(pVals)
}

# ******************************************************************************************* #


Run_MAST = function(rawData, gs){
  require(MAST)
  log_counts <- log(rawData+1)/log(2)
  fData = data.frame(names=rownames(log_counts))
  rownames(fData) = rownames(log_counts);
  cData = data.frame(cond=gs)
  rownames(cData) = colnames(log_counts)
  obj <- FromMatrix(as.matrix(log_counts), cData, fData)
  colData(obj)$cngeneson <- scale(colSums(assay(obj)>0))
  cond <- factor(colData(obj)$cond)
  # Model expression as function of condition & number of detected genes
  zlmCond <- zlm(~cond + cngeneson, obj)
  summaryCond <- summary(zlmCond, doLRT="cond1")
  summaryDt <- summaryCond$datatable
  summaryDt <- as.data.frame(summaryDt)
  pVals <- unlist(summaryDt[summaryDt$component == "H",4]) # H = hurdle model
  names(pVals) <- unlist(summaryDt[summaryDt$component == "H",1])
  pVals <- p.adjust(pVals, method = "fdr")
  pVals[which(is.na(pVals))]=1
  
  pVals = as.matrix(pVals)
  pid <- match(row.names(rawData),row.names(pVals))
  pVals2 = pVals[pid,]
  pVals2 = as.numeric(pVals2)
  
  return(pVals2)

}



####monocle
Run_monocle = function(rawData, gr1,gr2 , gs ){
  
  require(monocle)
  pd <- data.frame(group=gs)
  rownames(pd) <- colnames(rawData)
  pd <- new("AnnotatedDataFrame", data = pd)
  fData <- data.frame(gene_short_name = row.names(rawData), row.names = row.names(rawData))
  fd <- new("AnnotatedDataFrame", data = fData)
  Obj <- newCellDataSet(as.matrix(rawData), phenoData=pd, featureData = fd,
                        expressionFamily=VGAM::negbinomial.size())
  Obj <- estimateSizeFactors(Obj)
  Obj <- estimateDispersions(Obj)
  group=gs
  res <- differentialGeneTest(Obj,fullModelFormulaStr="~group")
  pVals <- res[,3]
  names(pVals) <- rownames(res)
  pVals[which(is.na(pVals))]=1
  pVals <- p.adjust(pVals, method = "fdr")
  

  return(pVals)
}

# ******************************************************************************************* #
# edgeR

Run_edgeR = function(rawData, gr1,gr2 , gs ){
  
  require(edgeR)
  
  # assign samples to groups and set up design matrix
  # gs <- factor(sml)
  #group = gset$group
  # sel=c(1,2,3,4,5,93,94,95,96,97)
  # ex = ex_new[,sel]
  # sml <- c(rep(0,5),rep(1,5))
  # sml = as.character(sml)
  # gs <- factor(sml)
  
  
  dge <- DGEList(counts=rawData, norm.factors = rep(1, ncol(rawData)), group=gs)
  group_edgeR <- gs
  design <- model.matrix(~group_edgeR)
  #dge <- estimateDisp(dge)
  dge <- estimateDisp(dge, design = design)
  fit <- glmFit(dge, design)
  res <- glmLRT(fit)
  pVals <- res$table[,4]
  pVals[which(is.na(pVals))]=1
  # names(pVals) <- rownames(res$table)
  pVals <- p.adjust(pVals, method = "fdr")

  return(pVals)
  
  
}


# ******************************************************************************************* #

# *** Running Limma:

Run_Limma = function(rawData, gs){
  
  require(limma)
  require(edgeR)
  
  # filtered = apply(rawData,1, function(x) {if(all(x == 0)) return (FALSE) else return(TRUE)})
  # FilteredData = rawData[filtered,]
  
  
  # Samples' conditions:
  mType <- gs
  
  # Normalization factors from edgeR TMM method:
  
  # nf <- calcNormFactors(FilteredData)
  
  design <- model.matrix(~mType+0)
  rownames(design)=colnames(rawData)
  
  
  
  fit <- lmFit(rawData,design)
  
  groups <- make.names(c("mType0", "mType1"))
  cts <- paste(groups[1], groups[2], sep="-")
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  
  fit <- contrasts.fit(fit, cont.matrix)
  
  fit <- eBayes(fit)
  
  # Summary of the results:
  
  Limma_results <- topTable(fit,n=nrow(fit))
  
  pVals=Limma_results[,5]
  pVals[which(is.na(pVals))]=1
  pVals = as.matrix(pVals)
  rownames(pVals)=row.names(Limma_results)
  pid <- match(row.names(rawData),row.names(pVals))
  pVals2 = pVals[pid,]
  pVals2 = as.numeric(pVals2)
  pval <- p.adjust(pVals2, method = "fdr")
  
  return(pval)

}

#####################################
# Wilcox/Mann-Whitney-U Test

Run_wilx = function(rawData, gr1,gr2){
  
  pVals <- apply(rawData, 1, function(x) {
    wilcox.test(x[1:gr1], 
                x[(gr1+1):(gr1+gr2)])$p.value
  })
  # multiple testing correction
  pVals[which(is.na(pVals))]=1
  pVals <- p.adjust(pVals, method = "fdr")
  # sigDE <- which(pVals < 0.05)
  Wilcox_pVals=pVals

  
  return(pVals)
  
}


Run_ttest = function(rawData, gr1,gr2){
  
  V = NULL
  pVals = NULL
  n_sam = dim(rawData)[2]
  ave_sam = cumsum(1:n_sam)[n_sam]/n_sam
  # for (tj in 1:dim(rawData)[1]){
  #   V = rank(rawData[tj,])
  #   if( sum(V==ave_sam)==n_sam){pVals[tj] = 0} else {pVals[tj] = t.test(rawData[tj, 1:gr1], rawData[tj, (gr1+1):(gr1+gr2)])$p.value}
  #   
  # }
  
  
  for (tj in 1:dim(rawData)[1]){
    V = rank(rawData[tj,])
    pVals[tj] = tryCatch({
      t.test(rawData[tj, 1:gr1], rawData[tj, (gr1+1):(gr1+gr2)])$p.value
      
    },error = function(e){
      if (identical(sort(V[1:gr1]),sort(V[(gr1+1):(gr1+gr2)]))){1} else {0}
    })
  }
  
  
  
  # pVals <- apply(rawData, 1, function(x) {
  #   t.test(x[1:gr1], x[(gr1+1):(gr1+gr2)])$p.value
  # })
  # multiple testing correction
  pVals[which(is.na(pVals))]=1
  pVals <- p.adjust(pVals, method = "fdr")
  
  return(pVals)
  
}


Run_ttestR = function(rawData, gr1,gr2){
  
  V = NULL
  pVals = NULL
  n_sam = dim(ex)[2]
  ave_sam = cumsum(1:n_sam)[n_sam]/n_sam
  # for (tj in 1:dim(ex)[1]){
  #   V = rank(ex[tj,])
  #   if( sum(V==ave_sam)==n_sam){pVals[tj] = 0} else {pVals[tj] = t.test(V[1:gr1], V[(gr1+1):(gr1+gr2)])$p.value}
  #   
  # }
  
  
  for (tj in 1:dim(rawData)[1]){
    V = rank(rawData[tj,])
    pVals[tj] = tryCatch({
      t.test(V[1:gr1], V[(gr1+1):(gr1+gr2)])$p.value
      
    },error = function(e){
      if (identical(sort(V[1:gr1]),sort(V[(gr1+1):(gr1+gr2)]))){1} else {0}
    })
  }
  
  
  # pVals <- apply(rawData, 1, function(x) {
  #   
  #   t.test(rank(x)[1:gr1], rank(x)[(gr1+1):(gr1+gr2)])$p.value
  #   
  # })
  # multiple testing correction
  pVals[which(is.na(pVals))]=1
  pVals <- p.adjust(pVals, method = "fdr")
  return(pVals)
  
}


Run_welch = function(rawData, gr1,gr2){
  
  V = NULL
  pVals = NULL
  for (tj in 1:dim(rawData)[1]){
    V = rank(rawData[tj,])
    pVals[tj] = tryCatch({
      t.test(rawData[tj, 1:gr1], rawData[tj, (gr1+1):(gr1+gr2)],var.equal=T,alternative="two.sided")$p.value
      
    },error = function(e){
      if (identical(sort(V[1:gr1]),sort(V[1:gr1]))){1} else {0}
    })
  }
  
  # pVals <- apply(rawData, 1, function(x) {
  #   t.test(x[1:gr1], x[(gr1+1):(gr1+gr2)],var.equal=T,alternative="two.sided")$p.value
  # })
  # multiple testing correction
  pVals[which(is.na(pVals))]=1
  pVals <- p.adjust(pVals, method = "fdr")
  # sigDE <- which(pVals < 0.05)
  Wilcox_pVals=pVals

  return(pVals)
  
}







Execute_DEsingle <- function(object, DEsingle.parallel = FALSE){
  
  object_DEsingle <- as.matrix(SummarizedExperiment::assay(object, "counts"))
  results_DEsingle <- DEsingle::DEsingle(counts = object_DEsingle,
                                         group = factor(object$label),
                                         parallel = DEsingle.parallel)
  result_DEsingle_DE <- list(gene_names = row.names(results_DEsingle), pvalue = results_DEsingle$pvalue,
                             FDR = results_DEsingle$pvalue.adj.FDR)
  
  pval=result_DEsingle_DE[["FDR"]]# results_DEsingle$pvalue.adj.FDR
  pval = as.data.frame(pval)
  rownames(pval)=rownames(results_DEsingle)
  genes = rownames(object@assays@data@listData[["counts"]])
  a=match(rownames(pval), genes)
  Dep = data.frame(pval=matrix(1, nrow = length(genes),ncol = 1))
  
  Dep[which(is.na(Dep))]=1

  return(Dep)
}


# perform zingeR.edgeR
Execute_zingeR.edgeR <- function(object, maxit.EM = 200){
  
  object_zingerR.edgeR <- as.matrix(SummarizedExperiment::assay(object, "counts"))
  
  dge <- edgeR::DGEList(round(object_zingerR.edgeR), group = factor(object$label))
  dge <- edgeR::calcNormFactors(dge)
  groups <- object$label
  design <- stats::model.matrix(~groups)
  weights <- zeroWeightsLS(counts=dge$counts, design=design, maxit = maxit.EM, normalization="TMM")
  dge$weights <- weights
  
  dge = edgeR::estimateDisp(dge, design)
  fit = edgeR::glmFit(dge, design)
  lrt = glmWeightedF(fit, coef=2, independentFiltering = TRUE)
  
  result_zingerR.edgeR <- list(gene_names = rownames(lrt$table),
                               pvalue = lrt$table$PValue,
                               FDR = lrt$table$padjFilter)
  Dep = result_zingerR.edgeR[["pvalue"]]

  Dep[which(is.na(Dep))]=1
  pval <- p.adjust(Dep, method = "fdr")
  return(pval)
}


data_process <- function(Data, group, norm.form = "CPM",  is.normalized = FALSE){
  options(warn = -1)
  
  if(is.normalized){
    
    normcounts <- Data
    gene_df <- data.frame(Gene = rownames(Data))
    cell_df <- data.frame(cell = colnames(Data))
    # # pd <- new("AnnotatedDataFrame", data = cell_df)
    # # fd <- new("AnnotatedDataFrame", data = gene_df)
    # # transfer <- new("CellDataSet", exprs = as.matrix(Data))
    # cds <- monocle::newCellDataSet(cellData = Data, expressionFamily = tobit())
    
    
    exp.matrix <- as(as.matrix(Data), 'sparseMatrix')
    feature_ann <- data.frame(gene_id=rownames(exp.matrix),gene_short_name=rownames(exp.matrix))
    rownames(feature_ann) <-row.names(Data)
    fd <- new("AnnotatedDataFrame", data = feature_ann)
    sample_ann <- data.frame(group=group)
    rownames(sample_ann) <- colnames(Data)
    pd <- new("AnnotatedDataFrame", data = sample_ann)
    cds <- newCellDataSet(exp.matrix, phenoData =pd, featureData =fd, expressionFamily=VGAM::negbinomial.size())
    
    # 
    counts_relative <- monocle::relative2abs(cds)
    counts_relative <- floor(counts_relative)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts_relative, normcounts = normcounts))
    gene_df <- data.frame(Gene = rownames(sce))
    cell_df <- data.frame(label = group, cell = colnames(sce))
    rownames(gene_df) <- gene_df$Gene
    rownames(cell_df) <- cell_df$cell
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = cds, normcounts = normcounts),
                                                      colData = cell_df,
                                                      rowData = gene_df)
    # sce <- scater::calculateQCMetrics(sce)
  } else
  {
    normcounts <- normalized(Data, method = norm.form)
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = Data, normcounts = normcounts))
    gene_df <- data.frame(Gene = rownames(sce))
    cell_df <- data.frame(label = group, cell = colnames(sce))
    rownames(gene_df) <- gene_df$Gene
    rownames(cell_df) <- cell_df$cell
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = Data, normcounts = normcounts),
                                                      colData = cell_df,
                                                      rowData = gene_df)
    # sce <- scater::calculateQCMetrics(sce)
  }
  
  return(sce)
}



normalized <- function(counts.matrix, method = "CPM"){
  if(method == "TMM"){
    norm_factor <- edgeR::calcNormFactors(counts.matrix, method = method)
    norm.item <- t(t(counts.matrix)/norm_factor)
    return(norm.item)
  }
  if(method == "RLE"){
    geomeans <- exp(rowMeans(log(counts.matrix)))
    SF <-function(cnts){
      stats::median((cnts/geomeans)[is.finite(geomeans) & geomeans >0])
    }
    norm_factor <- apply(counts.matrix, 2, SF)
    norm.item <- t(t(counts.matrix)/norm_factor)
    return(norm.item)
  }
  if(method == "CPM"){
    # sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts.matrix))
    # gene_df <- DataFrame(Gene = rownames(sce))
    # cell_df <- DataFrame(cell = colnames(sce))
    # rownames(gene_df) <- gene_df$Gene
    # rownames(cell_df) <- cell_df$cell
    # sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts.matrix),
    #                                                   colData = cell_df,
    #                                                   rowData = gene_df)
    # norm.item <- scater::calculateCPM(sce)
    norm_factor <- colSums(counts.matrix)
    norm.item <- t(t(counts.matrix)/norm_factor) * 1e6
    return(norm.item)
  }
  # if(method == "TPM"){
  #   sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts.matrix))
  #   gene_df <- DataFrame(Gene = rownames(sce))
  #   cell_df <- DataFrame(cell = colnames(sce))
  #   rownames(gene_df) <- gene_df$Gene
  #   rownames(cell_df) <- cell_df$cell
  #   sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts.matrix),
  #                                                     colData = cell_df,
  #                                                     rowData = gene_df)
  #   norm.item <- scater::calculateTPM(sce, exprs_values = "counts")
  #   return(norm.item)
  # }
}


Run_Agg <- function(rawdata, metas){
  
  cell_id2 = as.factor(as.vector(colnames(rawdata)) )  #get the cell id from the rawdata
  gene_id2 = as.factor(rownames(rawdata))   #get the gene id from the rawdata
  
  fData = data.frame(primerid = gene_id2)
  cData = data.frame(wellKey = cell_id2)
  
  
  
  
  sca2=SingleCellExperiment(list(counts=rawdata),
                            colData=data.frame(wellKey=cell_id2),
                            rowData=data.frame(primerid=gene_id2))
  
  # diagnosis = as.factor(metas$diagnosis)
  # diagnosis2=matrix("1",ncol=1,nrow=length(cell_id2))
  # diagnosis2[which(diagnosis == "Control")] = "0"
  # diagnosis= as.factor(diagnosis)
  
  # tmp = table(metas$individual)
  # inds = NULL
  # 
  # inds = ref[rep(1, tmp[[1]])]
  # for (i in 2:length(tmp)){
  #   inds = c(inds, ref[rep(i, tmp[[i]])])
  # }
  
  
  colData(sca2)$diagnosis = as.factor(metas$diagnosis)
  colData(sca2)$ind = as.factor(metas$individual)#as.factor(inds)#
  colData(sca2)$celltype = as.factor(metas$cluster)
  
  # colData(sca2)$diagnosis =as.factor(meta$diagnosis)[1:1314]
  # colData(sca2)$ind = as.factor(meta$individual)[1:10]
  # colData(sca2)$cluster = as.factor(meta$cluster)[1:1314]
  
  colData(sca2)
  # aggregateBioVar(scExp = sca2,subjectVar = "ind", cellVar = "celltype")
  
  agrBV_pval0=NA
  agrBV_pval1=NA
  
  cell_dds=NA
  subj_dds=NA
  
  
  #agreBV_pval0: DEseq2.
  # cell_dds =tryCatch(DESeqDataSetFromMatrix(countData = assay(sca2, "counts"),colData = colData(sca2),design = ~ diagnosis), error = function(e) {NA} )
  # cell_dds = tryCatch(DESeq(cell_dds), error = function(e) {NA} )
  # agrBV_pval0 = tryCatch(results(cell_dds)$pvalue, error = function(e) {NA} )
  
  #agreBV_pval1: do aggregate and then DEseq2.
  sca_argBV1 = tryCatch(aggregateBioVar(scExp = sca2,subjectVar = "ind", cellVar = "celltype"), error = function(e) {NA} )
  # subj_dds =tryCatch(DESeqDataSetFromMatrix(countData = assay(sca_argBV$type1, "counts"),
  #                                           colData = colData(sca_argBV$type1),design = ~ diagnosis),
  #                    error = function(e) {NA} )
  
  subj_dds =tryCatch(DESeqDataSetFromMatrix(
    countData = assay(sca_argBV1[[2]], "counts"),
    colData = colData(sca_argBV1[[2]]),
    design = ~ diagnosis),
    error = function(e) {NA} )
  subj_dds = tryCatch(DESeq(subj_dds), error = function(e) {NA} )
  
  agrBV_pval = results(subj_dds)@listData[["pvalue"]]
  pval[which(is.na(agrBV_pval))]=1
  pval <- p.adjust(pval, method = "fdr")
  
  return(pval)
}

###########################################################################################

Run_muscat <- function(rawdata, metas){
  
  
  
  
  #construct the format
  cell_id = colnames(rawdata)   #get the cell id from the rawdata
  gene_id = rownames(rawdata)   #get the gene id from the rawdata
  
  fData = data.frame(primerid = gene_id)
  cData = data.frame(wellKey = cell_id)
  
  
  sca2=SingleCellExperiment(list(counts=rawdata), 
                            colData=data.frame(wellKey=cell_id), 
                            rowData=data.frame(primerid=gene_id))
  
  metas$cluster = as.factor(metas$cluster)
  metas$individual = as.factor(metas$individual)
  metas$diagnosis = as.factor(metas$diagnosis)
  colData(sca2)$cluster_id = as.factor(metas$cluster)
  colData(sca2)$sample_id = as.factor(metas$individual)
  colData(sca2)$group_id = as.factor(metas$diagnosis)
  
  n_cells=table(metas$individual)
  n_cells=as.numeric(n_cells[match(unique(metas$individual),names(n_cells))])
  group_id=as.factor(metas$diagnosis[match(unique(metas$individual),metas$individual)])
  sample_id=as.factor(unique(metas$individual))
  
  metadata(sca2)$experiment_info=data.frame(sample_id=sample_id,group_id=group_id,n_cells=n_cells)
  
  # test
  pb2 = aggregateData(sca2,assay = "counts", fun = "sum",by = c("cluster_id", "sample_id"))
  res1 = pbDS(pb2, verbose = FALSE, method = c("edgeR"), filter = c("both"),min_cells = 0)
  # res1 = pbDS(pb1, verbose = FALSE, method = c("edgeR"), filter = c("both"),design = design,min_cells = 1)#pbDS(pb1, verbose = FALSE, min_cells = 1)
  tbl1 = res1$table[[1]]
  pb_pval1=tryCatch(tbl1[[1]]$p_val, error = function(e) {NA} )
  names(pb_pval1)=tryCatch(tbl1[[1]]$gene, error = function(e) {NA} )
  pb_pval1 = as.data.frame(pb_pval1)
  pval = matrix(data=1,nrow=length(gene_id),ncol=1)
  rownames(pval) = gene_id
  pval[match(rownames(pb_pval1), gene_id)] = pb_pval1$pb_pval1
  
  pval[which(is.na(pval))]=1
  pval <- p.adjust(pval, method = "fdr")
  
  return(pval)
}

