

rm(list= ls())
gc()


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
  
  pVals=Limma_results[,4]
  pVals[which(is.na(pVals))]=1
  # pVals = as.matrix(pVals)
  # rownames(pVals)=row.names(Limma_results)
  # pid <- match(row.names(rawData),row.names(pVals))
  # pVals2 = pVals[pid,]
  pVals = as.numeric(pVals)
  
  return(pVals)
}

# ******************************************************************************************* #




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
  
 
  result$power$Limma = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  result$power$ttest = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  result$power$ttestR = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  result$power$welch = matrix(nrow = length(DX), ncol = length(DX[[1]]))
  result$power$wilx = matrix(nrow = length(DX), ncol = length(DX[[1]]))

  
  result$p = list()
  
  for (d in 1:length(DX)){
    
    
    result$p$Limma[[d]] = matrix(nrow = Times, ncol = length(DX[[1]]))
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
      
      Res_Limma = Run_Limma(rawData = ex, gs)#Run_Limma(rawData = ex, gr1= m, gr2=n )
      Res_wilx = Run_wilx(rawData = ex, gr1= m, gr2=n)
      Res_ttest = Run_ttest(rawData = ex, gr1= m, gr2=n)
      Res_ttestR = Run_ttestR(rawData = ex, gr1= m, gr2=n)
      Res_welch = Run_welch(rawData = ex, gr1= m, gr2=n)
      
      

      result$p$Limma[[d]][,i] = Res_Limma
      result$p$ttest[[d]][,i] = Res_ttest
      result$p$ttestR[[d]][,i] = Res_ttestR
      result$p$welch[[d]][,i] = Res_welch
      result$p$wilx[[d]][,i] = Res_wilx
      
      
      
      
      result$power$Limma[d, i] = length(which(result$p$Limma[[d]][,i] <= 0.05)) / Times
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

data = readRDS("Simulation/gamma/powert/data_gamma_eq.rds")
re_gamma_eq_comp = wresult_comp(data)
saveRDS(re_gamma_eq_comp, "Simulation/gamma/powert/result_gamma_eq_comp.rds")

data = readRDS("Simulation/gamma/powert/data_gamma_ueq.rds")
re_gamma_ueq_comp = wresult_comp(data)
saveRDS(re_gamma_ueq_comp, "Simulation/gamma/powert/result_gamma_ueq_comp.rds")

