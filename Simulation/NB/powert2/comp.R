
rm(list= ls())
gc()


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

# *** Running DESeq:

# Run_DESeq = function(rawData , gr1 , gr2){
# 
#   require("DESeq2")
#   require("DESeq")
# 
#   # First the genes with 0 value for all the cells were filtered out:
# 
#   filtered = apply(rawData,1, function(x) {if(all(x == 0)) return (FALSE) else return(TRUE)})
# 
#   FilteredData = rawData[filtered,]
# 
#   # Defining the sample's description and conditions:
# 
#   dataPack = data.frame(row.names= colnames(FilteredData), condition = c(rep("GR1",gr1) ,rep("GR2",gr2)))
# 
#   conds <- factor (c(rep("GR1",gr1) ,rep("GR2",gr2)))
# 
#   # Initializing a CountDataSet (DESeq data structure for count matrix):
# 
#   cds <- newCountDataSet(FilteredData, conds)
# 
#   # DESeq normaliztion:
# 
#   cds <- estimateSizeFactors (cds)
#   print (sizeFactors(cds))
#   print (head(counts(cds, normalized = TRUE)))
# 
#   # Variance estimation:
# 
#   cds <- estimateDispersions (cds)
# 
#   print(str(fitInfo(cds)))
# 
#   # Negative binomial test:
# 
#   DESeq_results = nbinomTest(cds, "GR1", "GR2")
# 
#   save(DESeq_results, file ="DESeq_results.RData")
# 
# }


Run_DESeq = function(rawData , gs){
  
  require("DESeq2")
  
  dds <- DESeqDataSetFromMatrix(rawData, DataFrame(gs), ~ gs)
  ## factor levels were dropped which had no samples
  dds2 <- DESeq(dds,sfType = "poscounts")
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
  res <-  results(dds2, contrast=c("gs","0","1"),alpha = 0.05)
  pVals <- res@listData[["pvalue"]]#padj
  pVals[which(is.na(pVals))]=1
  # names(pVals) <- rownames(res)
  # sigDE <- which(pVals < 0.05)
  
  return(pVals)
  
}

# ******************************************************************************************* #

# *** Running Limma:

Run_Limma = function(rawData, gs){
  
  require(limma)
  require(edgeR)
  
  # filtered = apply(rawData,1, function(x) {if(all(x == 0)) return (FALSE) else return(TRUE)})
  # FilteredData = rawData[filtered,]
  
  
  mType <- gs
  
  # Normalization factors from edgeR TMM method:
  
  # nf <- calcNormFactors(FilteredData)
  
  design <- model.matrix(~mType+0)
  rownames(design)=colnames(rawData)
  dge <- DGEList(counts=rawData)
  dge <- calcNormFactors(dge)
  v <- voom( dge, design, normalize="quantile") #
  fit <- lmFit(v, design)
  
  groups <- make.names(c("mType0", "mType1"))
  cts <- paste(groups[1], groups[2], sep="-")
  cont.matrix <- makeContrasts(contrasts=cts, levels=design)
  
  fit <- contrasts.fit(fit, cont.matrix)
  
  fit <- eBayes(fit)
  
  # Summary of the results:
  
  Limma_results <- topTable(fit,n=nrow(fit))
  
  pVals=Limma_results[,4]
  pVals[which(is.na(pVals))]=1
  pVals = as.matrix(pVals)
  rownames(pVals)=row.names(Limma_results)
  pid <- match(row.names(rawData),row.names(pVals))
  pVals2 = pVals[pid,]
  pVals2 = as.numeric(pVals2)
  
  return(pVals2)
  # save(Limma_results , file = "Limma_results.RData")
  
  
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
  # pVals <- p.adjust(pVals, method = "fdr")
  pVals[which(is.na(pVals))]=1
  
  return(pVals)
  # sigDE <- names(pVals)[pVals < 0.05]
  # SCDE_pVals=pVals
  # DE_Quality_rate(sigDE)
  # DE_Quality_AUC(pVals)
}

# ******************************************************************************************* #

# *** Running MAST:


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
  # pVals <- p.adjust(pVals, method = "fdr")
  pVals[which(is.na(pVals))]=1
  
  pVals = as.matrix(pVals)
  pid <- match(row.names(rawData),row.names(pVals))
  pVals2 = pVals[pid,]
  pVals2 = as.numeric(pVals2)
  
  return(pVals2)
  
  
  # sigDE <- names(pVals)[pVals < 0.05]
  # MAST_pVals=pVals
  # DE_Quality_rate(sigDE)
  # DE_Quality_AUC(pVals)
  # return(pVals)
}


# ******************************************************************************************* #
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
  
  
  
  # sigDE <- which(pVals < 0.05)
  # monocle_pVals=pVals
  # DE_Quality_rate(sigDE)
  # DE_Quality_AUC(pVals)
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
  # pVals <- p.adjust(pVals, method = "fdr")
  # sigDE <- which(pVals < 0.05)
  # edgeR_pVals=pVals
  # DE_Quality_rate(sigDE)
  ## 0.8692022 0.3948121
  # DE_Quality_AUC(pVals)
  return(pVals)
  
  
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
  # pVals <- p.adjust(pVals, method = "fdr")
  
  # sigDE <- which(pVals < 0.05)
  Wilcox_pVals=pVals
  #DE_Quality_rate(sigDE)
  ## 0.3816953 0.08816031
  #DE_Quality_AUC(pVals) 
  ## [1] 0.8320326
  
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
      1
    })
  }
  
  
  
  # pVals <- apply(rawData, 1, function(x) {
  #   t.test(x[1:gr1], x[(gr1+1):(gr1+gr2)])$p.value
  # })
  # multiple testing correction
  pVals[which(is.na(pVals))]=1
  # pVals <- p.adjust(pVals, method = "fdr")
  # sigDE <- which(pVals < 0.05)
  Wilcox_pVals=pVals
  #DE_Quality_rate(sigDE)
  ## 0.3816953 0.08816031
  #DE_Quality_AUC(pVals) 
  ## [1] 0.8320326
  
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
      1
    })
  }
  
  
  # pVals <- apply(rawData, 1, function(x) {
  #   
  #   t.test(rank(x)[1:gr1], rank(x)[(gr1+1):(gr1+gr2)])$p.value
  #   
  # })
  # multiple testing correction
  pVals[which(is.na(pVals))]=1
  # pVals <- p.adjust(pVals, method = "fdr")
  # sigDE <- which(pVals < 0.05)
  # Wilcox_pVals=pVals
  #DE_Quality_rate(sigDE)
  ## 0.3816953 0.08816031
  #DE_Quality_AUC(pVals) 
  ## [1] 0.8320326
  
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
      1
    })
  }
  
  # pVals <- apply(rawData, 1, function(x) {
  #   t.test(x[1:gr1], x[(gr1+1):(gr1+gr2)],var.equal=T,alternative="two.sided")$p.value
  # })
  # multiple testing correction
  pVals[which(is.na(pVals))]=1
  # pVals <- p.adjust(pVals, method = "fdr")
  # sigDE <- which(pVals < 0.05)
  Wilcox_pVals=pVals
  #DE_Quality_rate(sigDE)
  ## 0.3816953 0.08816031
  #DE_Quality_AUC(pVals) 
  ## [1] 0.8320326
  
  return(pVals)
  
}


wresult_Limma = function(data){
  
  DX = data$X
  DY = data$Y
  Times = data$Times
  
  result=list()
  
  result$power = list()

  result$power$Limma = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  
  result$p = list()
  
  for (d in 1:length(DX)){
    
    result$p$Limma[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    
    for (i in 1:length(DX[[1]])){
      X = DX[[d]][[i]]
      Y = DY[[d]][[i]]
      m = dim(X)[2]
      n = dim(Y)[2]
      N = m+n
      
      sml <- c(rep(0,m),rep(1,n))
      sml = as.character(sml)
      gs <- factor(sml)
      ex = cbind(X,Y)
      rownames(ex) = 1:Times
      colnames(ex) = 1:N
      ex = as.data.frame(ex)
      
      ####################################################################################
      Res_Limma = Run_Limma(rawData = ex, gs )
      
      result$p$Limma[[d]][,i] = Res_Limma
      
      result$power$Limma[d, i] = length(which(result$p$Limma[[d]][,i] <= 0.05)) / Times
      
    }
  }
  
  
  return(result)
  
}




wresult_comp = function(data){
  
  DX = data$X
  DY = data$Y
  Times = data$Times
  
  result=list()
  
  result$power = list()
  
  result$power$monocle = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  result$power$MAST = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  result$power$scde = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  
  result$power$DESeq = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  result$power$Limma = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  result$power$edgeR = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  
  result$power$ttest = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  result$power$ttestR = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  result$power$welch = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  result$power$wilx = matrix(nrow = length(DX), ncol = length(DX[[1]]))

  
  result$p = list()
  
  for (d in 1:length(DX)){
    
    result$p$monocle[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    result$p$MAST[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    result$p$scde[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    
    result$p$DESeq[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    result$p$Limma[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    result$p$edgeR[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    
    result$p$ttest[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    result$p$ttestR[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    result$p$welch[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
    result$p$wilx[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))

    
    for (i in 1:length(DX[[1]])){
      X = DX[[d]][[i]]
      Y = DY[[d]][[i]]
      m = dim(X)[2]
      n = dim(Y)[2]
      N = m+n

      sml <- c(rep(0,m),rep(1,n))
      sml = as.character(sml)
      gs <- factor(sml)
      ex = cbind(X,Y)
      rownames(ex) = 1:Times
      colnames(ex) = 1:N
      ex = as.data.frame(ex)
      
      ####################################################################################
      # Res_ROTS = Run_ROTS(rawData = ex, gr1= m, gr2=n, B=1000, K=3000 , FDR = sig_score )
      Res_DESeq = Run_DESeq(rawData = ex, gs)
      Res_Limma = Run_Limma(rawData = ex, gs)#Run_Limma(rawData = ex, gr1= m, gr2=n )
      # Res_MAST = Run_MAST(rawData = ex, gs)
      # Res_monocle = Run_monocle(rawData = ex, gr1= m, gr2=n, gs )
      Res_edgeR = Run_edgeR(rawData = ex, gr1= m, gr2=n, gs )
      Res_wilx = Run_wilx(rawData = ex, gr1= m, gr2=n)
      # Res_scde = Run_scde(rawData = ex, gs)
      Res_ttest = Run_ttest(rawData = ex, gr1= m, gr2=n)
      Res_ttestR = Run_ttestR(rawData = ex, gr1= m, gr2=n)
      Res_welch = Run_welch(rawData = ex, gr1= m, gr2=n)
      
      # Res_MAST = tryCatch({
      #   Run_MAST(rawData = ex, gs)
      # 
      # },error = function(e){
      #   matrix(nrow = Times, ncol = 1)
      # })
      # 
      # Res_monocle = tryCatch({
      #   Run_monocle(rawData = ex, gr1= m, gr2=n, gs)
      # 
      # },error = function(e){
      #   matrix(nrow = Times, ncol = 1)
      # })
      # 
      # Res_scde = tryCatch({
      #   Run_scde(rawData = ex, gs)
      # 
      # },error = function(e){
      #   matrix(nrow = Times, ncol = 1)
      # })
      
      
      # Res$ROTS[[d]][,i] = Res_ROTS
      
      # result$p$MAST[[d]][,i] = Res_MAST
      # result$p$monocle[[d]][,i] = Res_monocle
      # result$p$scde[[d]][,i] = Res_scde
      
      result$p$edgeR[[d]][,i] = Res_edgeR
      result$p$DESeq[[d]][,i] = Res_DESeq
      result$p$Limma[[d]][,i] = Res_Limma
      
      result$p$ttest[[d]][,i] = Res_ttest
      result$p$ttestR[[d]][,i] = Res_ttestR
      result$p$welch[[d]][,i] = Res_welch
      result$p$wilx[[d]][,i] = Res_wilx
      
      
      # if (is.na(Res_scde[1])){
      #   result$power$scde[d, i] = NA
      # }else{
      #   result$power$scde[d, i] = length(which(result$p$scde[[d]][,i] <= 0.05)) / Times
      # }
      # 
      # if (is.na(Res_MAST[1])){
      #   result$power$MAST[d, i] = NA
      # }else{
      #   result$power$MAST[d, i] = length(which(result$p$MAST[[d]][,i] <= 0.05)) / Times
      # }
      # 
      # if (is.na(Res_monocle[1])){
      #   result$power$monocle[d, i] = NA
      # }else{
      #   result$power$monocle[d, i] = length(which(result$p$monocle[[d]][,i] <= 0.05)) / Times
      # }
        
      
      
      result$power$DESeq[d, i] = length(which(result$p$DESeq[[d]][,i] <= 0.05)) / Times
      result$power$Limma[d, i] = length(which(result$p$Limma[[d]][,i] <= 0.05)) / Times
      result$power$edgeR[d, i] = length(which(result$p$edgeR[[d]][,i] <= 0.05)) / Times
      
      result$power$ttest[d, i] = length(which(result$p$ttest[[d]][,i] <= 0.05)) / Times
      result$power$ttestR[d, i] = length(which(result$p$ttestR[[d]][,i] <= 0.05)) / Times
      result$power$welch[d, i] = length(which(result$p$welch[[d]][,i] <= 0.05)) / Times
      result$power$wilx[d, i] = length(which(result$p$wilx[[d]][,i] <= 0.05)) / Times
      
      
      }
  }
  
  
  return(result)
  
}


# ******************************************************************************************* #

# *** Read in the raw data:
setwd("C:/Users/dell/Desktop/wmwa/WMWA-R/NB")

data = readRDS("data_nb1_eq.rds")
re_NB1_eq_comp = wresult(data)
saveRDS(re_NB1_eq_comp, "result_NB1_eq_comp.rds")

data = readRDS("data_nb1_ueq.rds")
re_NB1_ueq_comp = wresult(data)
saveRDS(re_NB1_ueq_comp, "result_NB1_ueq_comp.rds")

data = readRDS("data_nb2_eq.rds")
re_NB2_eq_comp = wresult_comp(data)
saveRDS(re_NB2_eq_comp, "result_NB2_eq_comp.rds")

data = readRDS("data_nb2_ueq.rds")
re_NB2_ueq_comp = wresult(data)
saveRDS(re_NB2_ueq_comp, "result_NB2_ueq_comp.rds")


##############################################################################################

# data = readRDS("data_nb1_eq2.rds")
# re_NB1_eq2_comp = wresult(data)
# saveRDS(re_NB1_eq2_comp, "result_NB1_eq2_comp.rds")

# data = readRDS("data_nb1_ueq2.rds")
# re_NB1_ueq2_comp = wresult(data)
# saveRDS(re_NB1_ueq2_comp, "result_NB1_ueq2_comp.rds")
# 
# data = readRDS("data_nb2_eq2.rds")
# re_NB2_eq2_comp = wresult(data)
# saveRDS(re_NB2_eq2_comp, "result_NB2_eq2_comp.rds")
# 
# data = readRDS("data_nb2_ueq2.rds")
# re_NB2_ueq2_comp = wresult(data)
# saveRDS(re_NB2_ueq2_comp, "result_NB2_ueq2_comp.rds")
  
##############################################################################################

data = readRDS("data_nb1_eqt.rds")
re_NB1_eq = wresult_comp(data)
saveRDS(re_NB1_eq, "result_NB1_eqt_comp.rds")


data = readRDS("data_nb1_ueqt.rds")
re_NB1_ueq = wresult_comp(data)
saveRDS(re_NB1_ueq, "result_NB1_ueqt_comp.rds")


data = readRDS("data_nb2_eqt.rds")
re_NB2_eq = wresult_comp(data)
saveRDS(re_NB2_eq, "result_NB2_eqt_comp.rds")


data = readRDS("data_nb2_ueqt.rds")
re_NB2_ueq = wresult_comp(data)
saveRDS(re_NB2_ueq, "result_NB2_ueqt_comp.rds")



####################################################################################################
####################################################################################################
data = readRDS("data_nb1_eq.rds")
re_NB1_eq_comp = readRDS("result_NB1_eq_comp.rds")
relimma = wresult_Limma(data)
re_NB1_eq_comp[["power"]][["Limma"]] = relimma[["power"]][["Limma"]]
re_NB1_eq_comp[["p"]][["Limma"]] = relimma[["p"]][["Limma"]]
saveRDS(re_NB1_eq_comp, "result_NB1_eq_comp.rds")

data = readRDS("data_nb1_ueq.rds")
re_NB1_ueq_comp = readRDS("result_NB1_ueq_comp.rds")
relimma = wresult_Limma(data)
re_NB1_ueq_comp[["power"]][["Limma"]] = relimma[["power"]][["Limma"]]
re_NB1_ueq_comp[["p"]][["Limma"]] = relimma[["p"]][["Limma"]]
saveRDS(re_NB1_ueq_comp, "result_NB1_ueq_comp.rds")

data = readRDS("data_nb2_eq.rds")
re_NB2_eq_comp = readRDS("result_NB1_ueq_comp.rds")
relimma = wresult_Limma(data)
re_NB2_eq_comp[["power"]][["Limma"]] = relimma[["power"]][["Limma"]]
re_NB2_eq_comp[["p"]][["Limma"]] = relimma[["p"]][["Limma"]]
saveRDS(re_NB2_eq_comp, "result_NB2_eq_comp.rds")

data = readRDS("data_nb2_ueq.rds")
re_NB2_ueq_comp = readRDS("result_NB2_ueq_comp.rds")
relimma = wresult_Limma(data)
re_NB2_ueq_comp[["power"]][["Limma"]] = relimma[["power"]][["Limma"]]
re_NB2_ueq_comp[["p"]][["Limma"]] = relimma[["p"]][["Limma"]]
saveRDS(re_NB2_ueq_comp, "result_NB2_ueq_comp.rds")


##############################################################################################

# 
# 
# ##############################################################################################

data = readRDS("data_nb1_eqt.rds")
re_NB1_eq = readRDS("result_NB1_eqt_comp.rds")
relimma = wresult_Limma(data)
re_NB1_eq[["power"]][["Limma"]] = relimma[["power"]][["Limma"]]
re_NB1_eq[["p"]][["Limma"]] = relimma[["p"]][["Limma"]]
saveRDS(re_NB1_eq, "result_NB1_eqt_comp.rds")


data = readRDS("data_nb1_ueqt.rds")
re_NB1_ueq = readRDS("result_NB1_ueqt_comp.rds")
relimma = wresult_Limma(data)
re_NB1_ueq[["power"]][["Limma"]] = relimma[["power"]][["Limma"]]
re_NB1_ueq[["p"]][["Limma"]] = relimma[["p"]][["Limma"]]
saveRDS(re_NB1_ueq, "result_NB1_ueqt_comp.rds")


data = readRDS("data_nb2_eqt.rds")
re_NB2_eq = readRDS("result_NB2_eqt_comp.rds")
relimma = wresult_Limma(data)
re_NB2_eq[["power"]][["Limma"]] = relimma[["power"]][["Limma"]]
re_NB2_eq[["p"]][["Limma"]] = relimma[["p"]][["Limma"]]
saveRDS(re_NB2_eq, "result_NB2_eqt_comp.rds")


data = readRDS("data_nb2_ueqt.rds")
re_NB2_ueq = readRDS("result_NB2_ueqt_comp.rds")
relimma = wresult_Limma(data)
re_NB2_ueq[["power"]][["Limma"]] = relimma[["power"]][["Limma"]]
re_NB2_ueq[["p"]][["Limma"]] = relimma[["p"]][["Limma"]]
saveRDS(re_NB2_ueq, "result_NB2_ueqt_comp.rds")


