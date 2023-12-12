

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




# ******************************************************************************************* #

# *** Read in the raw data:
data = readMat('Simulation/gamma/comb/rate/data_rate.mat')

DX = list()
DY = list()
dim1 = dim(data$X)[1]
dim2 = dim(data$X)[2]
dim3 = dim(data$X)[3]
for (d in 1:dim1){
  DX[[d]] = list()
  DY[[d]] = list()
  for (i in 1:dim2){
    X_tmp = matrix(nrow = dim3, ncol = length(data[["X"]][[dim1*(i-1)+d]][[1]]))
    Y_tmp = matrix(nrow = dim3, ncol = length(data[["Y"]][[dim1*(i-1)+d]][[1]]))
    for (j in 1:dim3){
      X_tmp[j,] = data[["X"]][[dim1*dim2*(j-1)+dim1*(i-1)+d]][[1]]
      Y_tmp[j,] = data[["Y"]][[dim1*dim2*(j-1)+dim1*(i-1)+d]][[1]]
    }
    DX[[d]][[i]] = X_tmp
    DY[[d]][[i]] = Y_tmp
  }
}

Data_trans = list(X=DX, Y=DY, M=data[["n.case"]], N=data[["n.control"]], Times=dim3)


res = wresult_Limma(Data_trans)
re_limma = res[["power"]][["Limma"]]
saveRDS(re_limma, "Simulation/gamma/powert/rate_limma.rds")



