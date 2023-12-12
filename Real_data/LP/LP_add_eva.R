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



ntest = 6
sig_score = 0.01
# result<-readRDS("Real_data/LP/compresult.rds")

comb<-readRDS("Real_data/LP/result_comb.rds") 
result = comb$result

ngenes = 1:length(result[["Res"]][[1]][[1]][[1]])
degene = readRDS("Real_data/LP/degene_comp.rds")
golden = degene[["golden"]]

ntimes = 100

result$rate1$DEGnext = matrix(data = NA, nrow = ntimes, ncol = ntest)

result$rate2$DEGnext = matrix(data = NA, nrow = ntimes, ncol = ntest)

result$rate3$DEGnext = matrix(data = NA, nrow = ntimes, ncol = ntest)

result$rate4$DEGnext = matrix(data = NA, nrow = ntimes, ncol = ntest)

library(reticulate) 
np<-import("numpy") #datareading 
path = '/Users/guoyin/Desktop/github/'
mat_LP2<-np$load(paste0(path, "DEGnext_code/results/LP/LP_deg3.npy"), allow_pickle=TRUE) 
result[["Res"]][["DEGnext"]] = mat_LP2[[1]]

for (t in 1:ntimes){
  
  result$DE_gene$DEGnext = list()

  result$UNDE_gene$DEGnext = list()
  
  
  
  
  for (j in 1:ntest){
    

    result$DE_gene$DEGnext <- which(result$Res$DEGnext[[as.character(t-1)]][[as.character(j-1)]]==1)
    result$UNDE_gene$DEGnext = setdiff(ngenes,result$DE_gene$DEGnext)
    result$rate1$DEGnext[t,j] = DE_Quality_rate1(golden, result$DE_gene$DEGnext, result$UNDE_gene$DEGnext)
    result$rate2$DEGnext[t,j] = DE_Quality_rate2(golden, result$DE_gene$DEGnext, result$UNDE_gene$DEGnext)
    result$rate3$DEGnext[t,j] = DE_Quality_rate3(golden, result$DE_gene$DEGnext, result$UNDE_gene$DEGnext)
    result$rate4$DEGnext[t,j] = DE_Quality_rate4(golden, result$DE_gene$DEGnext, result$UNDE_gene$DEGnext)
    # result$auc$DEGnext[t, j] = DE_Quality_AUC(golden, result$Res$DEGnext[[as.character(t-1)]][[as.character(j-1)]])
    
  }
  
}

comb[["result"]] = result
saveRDS(comb,file="Real_data/LP/result_comb.rds")
# saveRDS(result,file="Real_data/LP/compresult_eva.rds")

