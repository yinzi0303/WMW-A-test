DE_Quality_AUC <- function(golden, pVals) {
  pVals <- pVals[which(pVals >= 0.05) %in% golden$DE | 
                   which(pVals >= 0.05) %in% golden$notDE]
  truth <- rep(1, times = length(pVals));
  truth[which(pVals >= 0.05) %in% golden$DE] = 0;
  pred <- ROCR::prediction(pVals, truth)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  # ROCR::plot(perf)
  aucObj <- ROCR::performance(pred, "auc")
  return(aucObj@y.values[[1]])
}
DE_Quality_rate1 <- function(golden, sigDE, unsigDE) {
  tp <- sum(golden$DE %in% sigDE)
  fp <- sum(golden$notDE %in% sigDE)
  tn <- sum(golden$notDE %in% unsigDE)
  fn <- sum(golden$DE %in% unsigDE)
  
  tpr <- tp/(tp + fn)
  fpr <- fp/(fp + tn)
  tnr = tn/(fp+tn)
  fnr = fn/(tp+fn)
  
  acc <- (tp + tn)/(tp + fn+fp + tn)
  
  return(acc)
}

DE_Quality_rate2 <- function(golden, sigDE, unsigDE) {
  tp <- sum(golden$DE %in% sigDE)
  fp <- sum(golden$notDE %in% sigDE)
  tn <- sum(golden$notDE %in% unsigDE)
  fn <- sum(golden$DE %in% unsigDE)
  
  precesion <- tp/(tp + fp)

  return(precesion)
}

DE_Quality_rate3 <- function(golden, sigDE, unsigDE) {
  tp <- sum(golden$DE %in% sigDE)
  fp <- sum(golden$notDE %in% sigDE)
  tn <- sum(golden$notDE %in% unsigDE)
  fn <- sum(golden$DE %in% unsigDE)
  
  tpr <- tp/(tp + fn)
  fpr <- fp/(fp + tn)
  
  # tpp <- sum(golden$DE %in% sigDE)/length(golden$DE)
  # tnp <- sum(golden$notDE %in% unsigDE)/length(golden$notDE)
  return(tpr)
}

DE_Quality_rate4 <- function(golden, sigDE, unsigDE) {
  tp <- sum(golden$DE %in% sigDE)
  fp <- sum(golden$notDE %in% sigDE)
  tn <- sum(golden$notDE %in% unsigDE)
  fn <- sum(golden$DE %in% unsigDE)
  
  tpr <- tp/(tp + fn)
  fpr <- fp/(fp + tn)
  
  # tpp <- sum(golden$DE %in% sigDE)/length(golden$DE)
  # tnp <- sum(golden$notDE %in% unsigDE)/length(golden$notDE)
  return(fpr)
}



#########################################################################################



dataset = readRDS("Real_data/MESC/data_count_MESC.rds")
ex_new = dataset[['expr']]
label = dataset$label

casesamples = which(label==1)
contsamples = which(label==2)



result = readRDS("Real_data/MESC/compresult.rds")

ntimes = 100

golden = NULL
deg = read.csv("Real_data/MESC/deg.csv",header = F)
de = match(deg[,1], rownames(data))
sigde = na.omit(de)
golden$DE = sigde
golden$notDE = setdiff(1:18370, sigde)


ntest = 8
sig_score = 0.01
ngenes = 1:length(result[["Res"]][["welch"]][[1]][[1]])


result$rate1 = list()
result$rate1$monocle = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$scde = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$MAST = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$DESeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$Limma = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$edgeR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$ttest = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$ttestR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$welch = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$wilx = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$zingeR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$DEsingle = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$ZIAQ = matrix(data = NA, nrow = ntimes, ncol = ntest)


result$rate2 = list()
result$rate2$monocle = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$scde = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$MAST = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$DESeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$Limma = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$edgeR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$ttest = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$ttestR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$welch = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$wilx = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$zingeR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$DEsingle = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$ZIAQ = matrix(data = NA, nrow = ntimes, ncol = ntest)




result$rate3 = list()
result$rate3$monocle = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$scde = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$MAST = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$DESeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$Limma = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$edgeR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$ttest = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$ttestR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$welch = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$wilx = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$zingeR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$DEsingle = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$ZIAQ = matrix(data = NA, nrow = ntimes, ncol = ntest)


result$rate4 = list()
result$rate4$monocle = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$scde = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$MAST = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$DESeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$Limma = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$edgeR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$ttest = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$ttestR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$welch = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$wilx = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$zingeR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$DEsingle = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$ZIAQ = matrix(data = NA, nrow = ntimes, ncol = ntest)



result$auc = list()
result$auc$monocle = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$auc$scde = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$auc$MAST = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$auc$DESeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$auc$Limma = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$auc$edgeR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$auc$ttest = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$auc$ttestR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$auc$welch = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$auc$wilx = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$auc$zingeR = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$auc$DEsingle = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$auc$ZIAQ = matrix(data = NA, nrow = ntimes, ncol = ntest)


ngenes = 1:dim(data)[1]

for (t in 1:ntimes){


  result$DE_gene$DESeq = list()
  result$DE_gene$Limma = list()
  result$DE_gene$MAST = list()
  result$DE_gene$monocle = list()
  result$DE_gene$edgeR = list()
  result$DE_gene$wilx = list()
  result$DE_gene$scde = list()
  result$DE_gene$ttest = list()
  result$DE_gene$ttestR = list()
  result$DE_gene$welch = list()
  result$DE_gene$zingeR = list()
  result$DE_gene$DEsingle = list()
  result$DE_gene$ZIAQ = list()
  
  

  result$UNDE_gene$DESeq = list()
  result$UNDE_gene$Limma = list()
  result$UNDE_gene$MAST = list()
  result$UNDE_gene$monocle = list()
  result$UNDE_gene$edgeR = list()
  result$UNDE_gene$wilx = list()
  result$UNDE_gene$scde = list()
  result$UNDE_gene$ttest = list()
  result$UNDE_gene$ttestR = list()
  result$UNDE_gene$welch = list()
  result$UNDE_gene$zingeR = list()
  result$UNDE_gene$DEsingle = list()
  result$UNDE_gene$ZIAQ = list()




  for (j in 1:ntest){


    result$DE_gene$DESeq <- which(result$Res$DESeq[[t]][[j]] < sig_score)
    result$UNDE_gene$DESeq = setdiff(ngenes,result$DE_gene$DESeq)
    result$rate1$DESeq[t,j] = DE_Quality_rate1(golden, result$DE_gene$DESeq, result$UNDE_gene$DESeq)
    result$rate2$DESeq[t,j] = DE_Quality_rate2(golden, result$DE_gene$DESeq, result$UNDE_gene$DESeq)
    result$rate3$DESeq[t,j] = DE_Quality_rate3(golden, result$DE_gene$DESeq, result$UNDE_gene$DESeq)
    result$rate4$DESeq[t,j] = DE_Quality_rate4(golden, result$DE_gene$DESeq, result$UNDE_gene$DESeq)
    result$auc$DESeq[t,j] = DE_Quality_AUC(golden, result$Res$DESeq[[t]][[j]])
    
    pVals = result$Res$Limma[[t]][[j]]
    pVals[which(is.na(pVals))]=1
    result$Res$Limma[[t]][[j]] = pVals 
    result$DE_gene$Limma <- which(result$Res$Limma[[t]][[j]] < sig_score)
    result$UNDE_gene$Limma = setdiff(ngenes,result$DE_gene$Limma)
    result$rate1$Limma[t,j] = DE_Quality_rate1(golden, result$DE_gene$Limma, result$UNDE_gene$Limma)
    result$rate2$Limma[t,j] = DE_Quality_rate2(golden, result$DE_gene$Limma, result$UNDE_gene$Limma)
    result$rate3$Limma[t,j] = DE_Quality_rate3(golden, result$DE_gene$Limma, result$UNDE_gene$Limma)
    result$rate4$Limma[t,j] = DE_Quality_rate4(golden, result$DE_gene$Limma, result$UNDE_gene$Limma)
    result$auc$Limma[t, j] = DE_Quality_AUC(golden, result$Res$Limma[[t]][[j]])

    result$DE_gene$MAST <- which(result$Res$MAST[[t]][[j]] < sig_score)
    result$UNDE_gene$MAST = setdiff(ngenes,result$DE_gene$MAST)
    result$rate1$MAST[t,j] = DE_Quality_rate1(golden, result$DE_gene$MAST, result$UNDE_gene$MAST)
    result$rate2$MAST[t,j] = DE_Quality_rate2(golden, result$DE_gene$MAST, result$UNDE_gene$MAST)
    result$rate3$MAST[t,j] = DE_Quality_rate3(golden, result$DE_gene$MAST, result$UNDE_gene$MAST)
    result$rate4$MAST[t,j] = DE_Quality_rate4(golden, result$DE_gene$MAST, result$UNDE_gene$MAST)
    result$auc$MAST[t, j] = DE_Quality_AUC(golden, result$Res$MAST[[t]][[j]])

    result$DE_gene$monocle <- which(result$Res$monocle[[t]][[j]] < sig_score)
    result$UNDE_gene$monocle = setdiff(ngenes,result$DE_gene$monocle)
    result$rate1$monocle[t,j] = DE_Quality_rate1(golden, result$DE_gene$monocle, result$UNDE_gene$monocle)
    result$rate2$monocle[t,j] = DE_Quality_rate2(golden, result$DE_gene$monocle, result$UNDE_gene$monocle)
    result$rate3$monocle[t,j] = DE_Quality_rate3(golden, result$DE_gene$monocle, result$UNDE_gene$monocle)
    result$rate4$monocle[t,j] = DE_Quality_rate4(golden, result$DE_gene$monocle, result$UNDE_gene$monocle)
    result$auc$monocle[t, j] = DE_Quality_AUC(golden, result$Res$monocle[[t]][[j]])

    result$DE_gene$edgeR <- which(result$Res$edgeR[[t]][[j]] < sig_score)
    result$UNDE_gene$edgeR = setdiff(ngenes,result$DE_gene$edgeR)
    result$rate1$edgeR[t,j] = DE_Quality_rate1(golden, result$DE_gene$edgeR, result$UNDE_gene$edgeR)
    result$rate2$edgeR[t,j] = DE_Quality_rate2(golden, result$DE_gene$edgeR, result$UNDE_gene$edgeR)
    result$rate3$edgeR[t,j] = DE_Quality_rate3(golden, result$DE_gene$edgeR, result$UNDE_gene$edgeR)
    result$rate4$edgeR[t,j] = DE_Quality_rate4(golden, result$DE_gene$edgeR, result$UNDE_gene$edgeR)
    result$auc$edgeR[t, j] = DE_Quality_AUC(golden, result$Res$edgeR[[t]][[j]])

    result$DE_gene$wilx <- which(result$Res$wilx[[t]][[j]] < sig_score)
    result$UNDE_gene$wilx = setdiff(ngenes,result$DE_gene$wilx)
    result$rate1$wilx[t,j] = DE_Quality_rate1(golden, result$DE_gene$wilx, result$UNDE_gene$wilx)
    result$rate2$wilx[t,j] = DE_Quality_rate2(golden, result$DE_gene$wilx, result$UNDE_gene$wilx)
    result$rate3$wilx[t,j] = DE_Quality_rate3(golden, result$DE_gene$wilx, result$UNDE_gene$wilx)
    result$rate4$wilx[t,j] = DE_Quality_rate4(golden, result$DE_gene$wilx, result$UNDE_gene$wilx)
    result$auc$wilx[t, j] = DE_Quality_AUC(golden, result$Res$wilx[[t]][[j]])

    result$DE_gene$scde <- which(result$Res$scde[[t]][[j]] < sig_score)
    result$UNDE_gene$scde = setdiff(ngenes,result$DE_gene$scde)
    result$rate1$scde[t,j] = DE_Quality_rate1(golden, result$DE_gene$scde, result$UNDE_gene$scde)
    result$rate2$scde[t,j] = DE_Quality_rate2(golden, result$DE_gene$scde, result$UNDE_gene$scde)
    result$rate3$scde[t,j] = DE_Quality_rate3(golden, result$DE_gene$scde, result$UNDE_gene$scde)
    result$rate4$scde[t,j] = DE_Quality_rate4(golden, result$DE_gene$scde, result$UNDE_gene$scde)
    result$auc$scde[t, j] = DE_Quality_AUC(golden, result$Res$scde[[t]][[j]])


    result$DE_gene$ttest <- which(result$Res$ttest[[t]][[j]] < sig_score)
    result$UNDE_gene$ttest = setdiff(ngenes,result$DE_gene$ttest)
    result$rate1$ttest[t,j] = DE_Quality_rate1(golden, result$DE_gene$ttest, result$UNDE_gene$ttest)
    result$rate2$ttest[t,j] = DE_Quality_rate2(golden, result$DE_gene$ttest, result$UNDE_gene$ttest)
    result$rate3$ttest[t,j] = DE_Quality_rate3(golden, result$DE_gene$ttest, result$UNDE_gene$ttest)
    result$rate4$ttest[t,j] = DE_Quality_rate4(golden, result$DE_gene$ttest, result$UNDE_gene$ttest)
    result$auc$ttest[t, j] = DE_Quality_AUC(golden, result$Res$ttest[[t]][[j]])


    result$DE_gene$ttestR <- which(result$Res$ttestR[[t]][[j]] < sig_score)
    result$UNDE_gene$ttestR = setdiff(ngenes,result$DE_gene$ttestR)
    result$rate1$ttestR[t,j] = DE_Quality_rate1(golden, result$DE_gene$ttestR, result$UNDE_gene$ttestR)
    result$rate2$ttestR[t,j] = DE_Quality_rate2(golden, result$DE_gene$ttestR, result$UNDE_gene$ttestR)
    result$rate3$ttestR[t,j] = DE_Quality_rate3(golden, result$DE_gene$ttestR, result$UNDE_gene$ttestR)
    result$rate4$ttestR[t,j] = DE_Quality_rate4(golden, result$DE_gene$ttestR, result$UNDE_gene$ttestR)
    result$auc$ttestR[t, j] = DE_Quality_AUC(golden, result$Res$ttestR[[t]][[j]])


    result$DE_gene$welch <- which(result$Res$welch[[t]][[j]] < sig_score)
    result$UNDE_gene$welch = setdiff(ngenes,result$DE_gene$welch)
    result$rate1$welch[t,j] = DE_Quality_rate1(golden, result$DE_gene$welch, result$UNDE_gene$welch)
    result$rate2$welch[t,j] = DE_Quality_rate2(golden, result$DE_gene$welch, result$UNDE_gene$welch)
    result$rate3$welch[t,j] = DE_Quality_rate3(golden, result$DE_gene$welch, result$UNDE_gene$welch)
    result$rate4$welch[t,j] = DE_Quality_rate4(golden, result$DE_gene$welch, result$UNDE_gene$welch)
    result$auc$welch[t, j] = DE_Quality_AUC(golden, result$Res$welch[[t]][[j]])
    
    
    
    result$DE_gene$DEsingle <- which(result$Res$DEsingle[[t]][[j]] < sig_score)
    result$UNDE_gene$DEsingle = setdiff(ngenes,result$DE_gene$DEsingle)
    result$rate1$DEsingle[t,j] = DE_Quality_rate1(golden, result$DE_gene$DEsingle, result$UNDE_gene$DEsingle)
    result$rate2$DEsingle[t,j] = DE_Quality_rate2(golden, result$DE_gene$DEsingle, result$UNDE_gene$DEsingle)
    result$rate3$DEsingle[t,j] = DE_Quality_rate3(golden, result$DE_gene$DEsingle, result$UNDE_gene$DEsingle)
    result$rate4$DEsingle[t,j] = DE_Quality_rate4(golden, result$DE_gene$DEsingle, result$UNDE_gene$DEsingle)
    result$auc$DEsingle[t,j] = DE_Quality_AUC(golden, result$Res$DEsingle[[t]][[j]])
    
    
    
    result$DE_gene$zingeR <- which(result$Res$zingeR[[t]][[j]] < sig_score)
    result$UNDE_gene$zingeR = setdiff(ngenes,result$DE_gene$zingeR)
    result$rate1$zingeR[t,j] = DE_Quality_rate1(golden, result$DE_gene$zingeR, result$UNDE_gene$zingeR)
    result$rate2$zingeR[t,j] = DE_Quality_rate2(golden, result$DE_gene$zingeR, result$UNDE_gene$zingeR)
    result$rate3$zingeR[t,j] = DE_Quality_rate3(golden, result$DE_gene$zingeR, result$UNDE_gene$zingeR)
    result$rate4$zingeR[t,j] = DE_Quality_rate4(golden, result$DE_gene$zingeR, result$UNDE_gene$zingeR)
    result$auc$zingeR[t, j] = DE_Quality_AUC(golden, result$Res$zingeR[[t]][[j]])
    
    
    result$DE_gene$ZIAQ <- which(result$Res$ZIAQ[[t]][[j]] < sig_score)
    result$UNDE_gene$ZIAQ = setdiff(ngenes,result$DE_gene$ZIAQ)
    result$rate1$ZIAQ[t,j] = DE_Quality_rate1(golden, result$DE_gene$ZIAQ, result$UNDE_gene$ZIAQ)
    result$rate2$ZIAQ[t,j] = DE_Quality_rate2(golden, result$DE_gene$ZIAQ, result$UNDE_gene$ZIAQ)
    result$rate3$ZIAQ[t,j] = DE_Quality_rate3(golden, result$DE_gene$ZIAQ, result$UNDE_gene$ZIAQ)
    result$rate4$ZIAQ[t,j] = DE_Quality_rate4(golden, result$DE_gene$ZIAQ, result$UNDE_gene$ZIAQ)
    result$auc$ZIAQ[t, j] = DE_Quality_AUC(golden, result$Res$ZIAQ[[t]][[j]])
    
    

  }

}

saveRDS(result,file="Real_data/MESC/compresult_eva.rds")

