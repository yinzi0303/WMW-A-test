
DE_Quality_AUC <- function(pVals) {
  pVals <- pVals[which(pVals >= 0.05) %in% golden$DE | 
                   which(pVals >= 0.05) %in% golden$notDE]
  truth <- rep(1, times = length(pVals));
  truth[which(pVals >= 0.05) %in% golden$DE] = 0;
  pred <- ROCR::prediction(pVals, truth)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  ROCR::plot(perf)
  aucObj <- ROCR::performance(pred, "auc")
  return(aucObj@y.values[[1]])
}
DE_Quality_rate1 <- function(sigDE, unsigDE) {
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

DE_Quality_rate2 <- function(sigDE, unsigDE) {
  tp <- sum(golden$DE %in% sigDE)
  fp <- sum(golden$notDE %in% sigDE)
  tn <- sum(golden$notDE %in% unsigDE)
  fn <- sum(golden$DE %in% unsigDE)
  
  precesion <- tp/(tp + fp)
  
  return(precesion)
}

DE_Quality_rate3 <- function(sigDE, unsigDE) {
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

DE_Quality_rate4 <- function(sigDE, unsigDE) {
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




myresult=readRDS("Real_data/MESC/result/WMWAresult.rds") 
result = readRDS("Real_data/MESC/compresult_eva.rds")
degene_comp = readRDS("Real_data/MESC/degene_comp.rds")
golden = NULL
deg = read.csv("Real_data/MESC/deg.csv",header = F)
de = match(deg[,1], rownames(data))
sigde = na.omit(de)
golden$DE = sigde
golden$notDE = setdiff(1:18370, sigde)


for (i in 1:length(myresult[["Res"]])){
  for (j in 1:length(myresult[["Res"]][[1]])){
    for (k in 1:6){
      pval = myresult[["Res"]][[i]][[j]][[k]]
      myresult[["Res"]][[i]][[j]][[k]]<- p.adjust(pval, method = "fdr")
    }
    
  }
}

result$Res$WMWAN = myresult$Res$WMWANL
result$Res$WMWANM = myresult$Res$WMWANLM
result$Res$WMWANB = myresult$Res$WMWANB
result$Res$WMWANBM = myresult$Res$WMWANBM
result$Res$WMWA_count = myresult$Res$WMWA_count
result$Res$WMWA_norm = myresult$Res$WMWA_norm


ntest = 6
sig_score = 0.01
ngenes = 1:length(result[["Res"]][[1]][[1]][[1]])
ntime = length(myresult[["Res"]][[1]])


result$rate1$WMWAN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate1$WMWANM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate1$WMWANB = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate1$WMWANBM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate1$WMWA_count = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate1$WMWA_norm = matrix(data = NA, nrow = ntime, ncol = ntest)


result$rate2$WMWAN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate2$WMWANM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate2$WMWANB = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate2$WMWANBM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate2$WMWA_count = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate2$WMWA_norm = matrix(data = NA, nrow = ntime, ncol = ntest)


result$rate3$WMWAN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate3$WMWANM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate3$WMWANB = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate3$WMWANBM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate3$WMWA_count = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate3$WMWA_norm = matrix(data = NA, nrow = ntime, ncol = ntest)


result$rate4$WMWAN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate4$WMWANM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate4$WMWANB = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate4$WMWANBM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate4$WMWA_count = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate4$WMWA_norm = matrix(data = NA, nrow = ntime, ncol = ntest)


result$auc$WMWAN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$auc$WMWANM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$auc$WMWANB = matrix(data = NA, nrow = ntime, ncol = ntest)
result$auc$WMWANBM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$auc$WMWA_count = matrix(data = NA, nrow = ntime, ncol = ntest)
result$auc$WMWA_norm = matrix(data = NA, nrow = ntime, ncol = ntest)




for (t in 1:ntime){
  
  result$DE_gene$WMWAN = list()
  result$DE_gene$WMWANM = list()
  result$DE_gene$WMWANB = list()
  result$DE_gene$WMWANBM = list()
  result$DE_gene$WMWA_count = list()
  result$DE_gene$WMWA_norm = list()
  
  
  result$UNDE_gene$WMWAN = list()
  result$UNDE_gene$WMWANM = list()
  result$UNDE_gene$WMWANB = list()
  result$UNDE_gene$WMWANBM = list()
  result$UNDE_gene$WMWA_count = list()
  result$UNDE_gene$WMWA_norm = list()
  
  
  
  for (j in 1:ntest){
    
    
    result$DE_gene$WMWAN <- which(result$Res$WMWAN[[t]][[j]] < sig_score)
    result$UNDE_gene$WMWAN = setdiff(ngenes,result$DE_gene$WMWAN)
    result$rate1$WMWAN[t,j] = DE_Quality_rate1(golden, result$DE_gene$WMWAN, result$UNDE_gene$WMWAN)
    result$rate2$WMWAN[t,j] = DE_Quality_rate2(golden, result$DE_gene$WMWAN, result$UNDE_gene$WMWAN)
    result$rate3$WMWAN[t,j] = DE_Quality_rate3(golden, result$DE_gene$WMWAN, result$UNDE_gene$WMWAN)
    result$rate4$WMWAN[t,j] = DE_Quality_rate4(golden, result$DE_gene$WMWAN, result$UNDE_gene$WMWAN)
    result$auc$WMWAN[t,j] = DE_Quality_AUC(golden, result$Res$WMWAN[[t]][[j]])
    
    result$DE_gene$WMWANM <- which(result$Res$WMWANM[[t]][[j]] < sig_score)
    result$UNDE_gene$WMWANM = setdiff(ngenes,result$DE_gene$WMWANM)
    result$rate1$WMWANM[t,j] = DE_Quality_rate1(golden, result$DE_gene$WMWANM, result$UNDE_gene$WMWANM)
    result$rate2$WMWANM[t,j] = DE_Quality_rate2(golden, result$DE_gene$WMWANM, result$UNDE_gene$WMWANM)
    result$rate3$WMWANM[t,j] = DE_Quality_rate3(golden, result$DE_gene$WMWANM, result$UNDE_gene$WMWANM)
    result$rate4$WMWANM[t,j] = DE_Quality_rate4(golden, result$DE_gene$WMWANM, result$UNDE_gene$WMWANM)
    result$auc$WMWANM[t,j] = DE_Quality_AUC(golden, result$Res$WMWANM[[t]][[j]])
    
    
    result$DE_gene$WMWANB <- which(result$Res$WMWANB[[t]][[j]] < sig_score)
    result$UNDE_gene$WMWANB = setdiff(ngenes,result$DE_gene$WMWANB)
    result$rate1$WMWANB[t,j] = DE_Quality_rate1(golden, result$DE_gene$WMWANB, result$UNDE_gene$WMWANB)
    result$rate2$WMWANB[t,j] = DE_Quality_rate2(golden, result$DE_gene$WMWANB, result$UNDE_gene$WMWANB)
    result$rate3$WMWANB[t,j] = DE_Quality_rate3(golden, result$DE_gene$WMWANB, result$UNDE_gene$WMWANB)
    result$rate4$WMWANB[t,j] = DE_Quality_rate4(golden, result$DE_gene$WMWANB, result$UNDE_gene$WMWANB)
    result$auc$WMWANB[t,j] = DE_Quality_AUC(golden, result$Res$WMWANB[[t]][[j]])
    
    result$DE_gene$WMWANBM <- which(result$Res$WMWANBM[[t]][[j]] < sig_score)
    result$UNDE_gene$WMWANBM = setdiff(ngenes,result$DE_gene$WMWANBM)
    result$rate1$WMWANBM[t,j] = DE_Quality_rate1(golden, result$DE_gene$WMWANBM, result$UNDE_gene$WMWANBM)
    result$rate2$WMWANBM[t,j] = DE_Quality_rate2(golden, result$DE_gene$WMWANBM, result$UNDE_gene$WMWANBM)
    result$rate3$WMWANBM[t,j] = DE_Quality_rate3(golden, result$DE_gene$WMWANBM, result$UNDE_gene$WMWANBM)
    result$rate4$WMWANBM[t,j] = DE_Quality_rate4(golden, result$DE_gene$WMWANBM, result$UNDE_gene$WMWANBM)
    result$auc$WMWANBM[t,j] = DE_Quality_AUC(golden, result$Res$WMWANBM[[t]][[j]])
    

    result$DE_gene$WMWA_count <- which(result$Res$WMWA_count[[t]][[j]] < sig_score)
    result$UNDE_gene$WMWA_count = setdiff(ngenes,result$DE_gene$WMWA_count)
    result$rate1$WMWA_count[t,j] = DE_Quality_rate1(golden, result$DE_gene$WMWA_count, result$UNDE_gene$WMWA_count)
    result$rate2$WMWA_count[t,j] = DE_Quality_rate2(golden, result$DE_gene$WMWA_count, result$UNDE_gene$WMWA_count)
    result$rate3$WMWA_count[t,j] = DE_Quality_rate3(golden, result$DE_gene$WMWA_count, result$UNDE_gene$WMWA_count)
    result$rate4$WMWA_count[t,j] = DE_Quality_rate4(golden, result$DE_gene$WMWA_count, result$UNDE_gene$WMWA_count)
    result$auc$WMWA_count[t,j] = DE_Quality_AUC(golden, result$Res$WMWA_count[[t]][[j]])
    
    
    result$DE_gene$WMWA_norm <- which(result$Res$WMWA_norm[[t]][[j]] < sig_score)
    result$UNDE_gene$WMWA_norm = setdiff(ngenes,result$DE_gene$WMWA_norm)
    result$rate1$WMWA_norm[t,j] = DE_Quality_rate1(golden, result$DE_gene$WMWA_norm, result$UNDE_gene$WMWA_norm)
    result$rate2$WMWA_norm[t,j] = DE_Quality_rate2(golden, result$DE_gene$WMWA_norm, result$UNDE_gene$WMWA_norm)
    result$rate3$WMWA_norm[t,j] = DE_Quality_rate3(golden, result$DE_gene$WMWA_norm, result$UNDE_gene$WMWA_norm)
    result$rate4$WMWA_norm[t,j] = DE_Quality_rate4(golden, result$DE_gene$WMWA_norm, result$UNDE_gene$WMWA_norm)
    result$auc$WMWA_norm[t,j] = DE_Quality_AUC(golden, result$Res$WMWA_norm[[t]][[j]])
    
    
    
    
  }
  
}

comb = list(golden=golden, degene=degene, result=result)
saveRDS(comb,file="Real_data/MESC/result_comb.rds")
