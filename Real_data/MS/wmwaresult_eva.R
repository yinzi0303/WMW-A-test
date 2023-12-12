
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



myresult=readRDS("Real_data/MS/WMWAresult.rds") 
result = readRDS("Real_data/MS/compresult_eva.rds")
degene_comp = readRDS("Real_data/MS/degene_comp.rds")
golden = degene_comp[["golden"]]
degene = degene_comp[["degene"]]

for (i in 1:length(myresult[["Res"]])){
  for (j in 1:length(myresult[["Res"]][[1]])){
    for (k in 1:8){
      pval = myresult[["Res"]][[i]][[j]][[k]]
      myresult[["Res"]][[i]][[j]][[k]]<- p.adjust(pval, method = "fdr")
    }
    
  }
}

result$Res$WMWAN = myresult$Res$WMWAN
result$Res$WMWANM = myresult$Res$WMWANM
result$Res$WMWALN = myresult$Res$WMWALN
result$Res$WMWALNM = myresult$Res$WMWALNM
result$Res$WMWANB = myresult$Res$WMWANB
result$Res$WMWANBM = myresult$Res$WMWANBM
result$Res$WMWAreal = myresult$Res$WMWAreal


ntest = 8
sig_score = 0.01
ngenes = 1:length(result[["Res"]][[1]][[1]][[1]])
ntime = length(myresult[["Res"]][["WMWAN"]])


result$rate1$WMWAN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate1$WMWANM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate1$WMWALN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate1$WMWALNM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate1$WMWANB = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate1$WMWANBM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate1$WMWAreal = matrix(data = NA, nrow = ntime, ncol = ntest)


result$rate2$WMWAN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate2$WMWANM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate2$WMWALN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate2$WMWALNM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate2$WMWANB = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate2$WMWANBM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate2$WMWAreal = matrix(data = NA, nrow = ntime, ncol = ntest)


result$rate3$WMWAN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate3$WMWANM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate3$WMWALN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate3$WMWALNM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate3$WMWANB = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate3$WMWANBM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate3$WMWAreal = matrix(data = NA, nrow = ntime, ncol = ntest)


result$rate4$WMWAN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate4$WMWANM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate4$WMWALN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate4$WMWALNM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate4$WMWANB = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate4$WMWANBM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$rate4$WMWAreal = matrix(data = NA, nrow = ntime, ncol = ntest)



result$auc$WMWAN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$auc$WMWANM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$auc$WMWALN = matrix(data = NA, nrow = ntime, ncol = ntest)
result$auc$WMWALNM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$auc$WMWANB = matrix(data = NA, nrow = ntime, ncol = ntest)
result$auc$WMWANBM = matrix(data = NA, nrow = ntime, ncol = ntest)
result$auc$WMWAreal = matrix(data = NA, nrow = ntime, ncol = ntest)




for (t in 1:ntime){
  
  result$DE_gene$WMWAN = list()
  result$DE_gene$WMWANM = list()
  result$DE_gene$WMWALN = list()
  result$DE_gene$WMWALNM = list()
  result$DE_gene$WMWANB = list()
  result$DE_gene$WMWANBM = list()
  result$DE_gene$WMWAreal = list()
  
  
  result$UNDE_gene$WMWAN = list()
  result$UNDE_gene$WMWANM = list()
  result$UNDE_gene$WMWALN = list()
  result$UNDE_gene$WMWALNM = list()
  result$UNDE_gene$WMWANB = list()
  result$UNDE_gene$WMWANBM = list()
  result$UNDE_gene$WMWAreal = list()
  
  
  
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
    
    
    result$DE_gene$WMWALN <- which(result$Res$WMWALN[[t]][[j]] < sig_score)
    result$UNDE_gene$WMWALN = setdiff(ngenes,result$DE_gene$WMWALN)
    result$rate1$WMWALN[t,j] = DE_Quality_rate1(golden, result$DE_gene$WMWALN, result$UNDE_gene$WMWALN)
    result$rate2$WMWALN[t,j] = DE_Quality_rate2(golden, result$DE_gene$WMWALN, result$UNDE_gene$WMWALN)
    result$rate3$WMWALN[t,j] = DE_Quality_rate3(golden, result$DE_gene$WMWALN, result$UNDE_gene$WMWALN)
    result$rate4$WMWALN[t,j] = DE_Quality_rate4(golden, result$DE_gene$WMWALN, result$UNDE_gene$WMWALN)
    result$auc$WMWALN[t,j] = DE_Quality_AUC(golden, result$Res$WMWALN[[t]][[j]])
    
    result$DE_gene$WMWALNM <- which(result$Res$WMWALNM[[t]][[j]] < sig_score)
    result$UNDE_gene$WMWALNM = setdiff(ngenes,result$DE_gene$WMWALNM)
    result$rate1$WMWALNM[t,j] = DE_Quality_rate1(golden, result$DE_gene$WMWALNM, result$UNDE_gene$WMWALNM)
    result$rate2$WMWALNM[t,j] = DE_Quality_rate2(golden, result$DE_gene$WMWALNM, result$UNDE_gene$WMWALNM)
    result$rate3$WMWALNM[t,j] = DE_Quality_rate3(golden, result$DE_gene$WMWALNM, result$UNDE_gene$WMWALNM)
    result$rate4$WMWALNM[t,j] = DE_Quality_rate4(golden, result$DE_gene$WMWALNM, result$UNDE_gene$WMWALNM)
    result$auc$WMWALNM[t,j] = DE_Quality_AUC(golden, result$Res$WMWALNM[[t]][[j]])
    
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
    
    result$DE_gene$WMWAreal <- which(result$Res$WMWAreal[[t]][[j]] < sig_score)
    result$UNDE_gene$WMWAreal = setdiff(ngenes,result$DE_gene$WMWAreal)
    result$rate1$WMWAreal[t,j] = DE_Quality_rate1(golden, result$DE_gene$WMWAreal, result$UNDE_gene$WMWAreal)
    result$rate2$WMWAreal[t,j] = DE_Quality_rate2(golden, result$DE_gene$WMWAreal, result$UNDE_gene$WMWAreal)
    result$rate3$WMWAreal[t,j] = DE_Quality_rate3(golden, result$DE_gene$WMWAreal, result$UNDE_gene$WMWAreal)
    result$rate4$WMWAreal[t,j] = DE_Quality_rate4(golden, result$DE_gene$WMWAreal, result$UNDE_gene$WMWAreal)
    result$auc$WMWAreal[t,j] = DE_Quality_AUC(golden, result$Res$WMWAreal[[t]][[j]])
    
    
  }
  
}

comb = list(golden=golden, degene=degene, result=result)
saveRDS(comb,file="result_comb2.rds")
