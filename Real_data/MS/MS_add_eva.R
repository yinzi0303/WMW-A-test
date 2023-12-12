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



ntimes = 100
ntest = 8
sig_score = 0.01
result_comb<-readRDS("Real_data/MS/result_comb2.rds") 
result = result_comb[["result"]]
degene = result_comb[["degene"]]
ngenes = 1:length(result[["Res"]][["welch"]][[1]][[1]])

deg = list()
undeg = list()
for (d in 1:15){
  deg[[d]] = which(degene[[d]]<=sig_score)
  undeg[[d]] = setdiff(ngenes, deg[[d]])
}

d = 1:15
d = d[-c(2,3,4,5,11,12)]
de_union = NULL
de_inter = deg[[1]]
for (e in d){
  de_union = union(de_union, deg[[e]])
  de_inter = intersect(de_inter, deg[[e]])
}
de_not_inter = setdiff(ngenes, de_inter)
golden = list(notDE=de_not_inter, DE=de_inter)



# golden = matrix(data=NA, nrow = length(result[["Res"]][["welch"]][[1]][[1]]), ncol = 1)
# golden[de_not_inter] = 0
# golden[de_inter] = 1
# write.csv(golden, 'tung/golden.csv')


result$rate1$consexpression = matrix(data = NA, nrow = ntimes, ncol = ntest)
# result$rate1$consexpression2 = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$BaySeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$NOISeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate1$DEGnext = matrix(data = NA, nrow = ntimes, ncol = ntest)


result$rate2$consexpression = matrix(data = NA, nrow = ntimes, ncol = ntest)
# result$rate2$consexpression2 = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$BaySeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$NOISeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate2$DEGnext = matrix(data = NA, nrow = ntimes, ncol = ntest)


result$rate3$consexpression = matrix(data = NA, nrow = ntimes, ncol = ntest)
# result$rate3$consexpression2 = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$BaySeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$NOISeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate3$DEGnext = matrix(data = NA, nrow = ntimes, ncol = ntest)


result$rate4$consexpression = matrix(data = NA, nrow = ntimes, ncol = ntest)
# result$rate4$consexpression2 = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$BaySeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$NOISeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
result$rate4$DEGnext = matrix(data = NA, nrow = ntimes, ncol = ntest)

# result$auc$consexpression = matrix(data = NA, nrow = ntimes, ncol = ntest)
# result$auc$consexpression2 = matrix(data = NA, nrow = ntimes, ncol = ntest)
# result$auc$BaySeq = matrix(data = NA, nrow = ntimes, ncol = ntest)
# result$auc$NOISeq = matrix(data = NA, nrow = ntimes, ncol = ntest)



library(reticulate) 
np<-import("numpy") #datareading 
path = '/Users/guoyin/Desktop/github/'
mat_MS<-np$load(paste0(path, "consexpression-master/results/MS/MS_deg.npy"), allow_pickle=TRUE) 
result[["Res"]][["consexpression"]] = mat_MS[[1]][["Consexp"]]
# result[["Res"]][["consexpression2"]] = mat_MS[[1]][["Consexp_2"]]
result[["Res"]][["BaySeq"]] = mat_MS[[1]][["BaySeq"]]
result[["Res"]][["NOISeq"]] = mat_MS[[1]][["NOISeq"]]

mat_MS2<-np$load(paste0(path, "DEGnext_code/results/MS/MS_deg3.npy"), allow_pickle=TRUE) 
result[["Res"]][["DEGnext"]] = mat_MS2[[1]]

for (t in 1:ntimes){
  
  
  # result$DE_gene$consexpression2 = list()
  result$DE_gene$BaySeq = list()
  result$DE_gene$consexpression = list()
  result$DE_gene$NOISeq = list()
  result$DE_gene$DEGnext = list()
  
  
  
  # result$UNDE_gene$consexpression2 = list()
  result$UNDE_gene$BaySeq = list()
  result$UNDE_gene$consexpression = list()
  result$UNDE_gene$NOISeq = list()
  result$UNDE_gene$DEGnext = list()
  
  
  
  
  for (j in 1:ntest){
    
    # result$DE_gene$consexpression2 <- result$Res$consexpression2[[as.character(t-1)]][[as.character(j-1)]]
    # result$UNDE_gene$consexpression2 = setdiff(ngenes,result$DE_gene$consexpression2)
    # result$rate1$consexpression2[t,j] = DE_Quality_rate1(golden, result$DE_gene$consexpression2, result$UNDE_gene$consexpression2)
    # result$rate2$consexpression2[t,j] = DE_Quality_rate2(golden, result$DE_gene$consexpression2, result$UNDE_gene$consexpression2)
    # result$rate3$consexpression2[t,j] = DE_Quality_rate3(golden, result$DE_gene$consexpression2, result$UNDE_gene$consexpression2)
    # result$rate4$consexpression2[t,j] = DE_Quality_rate4(golden, result$DE_gene$consexpression2, result$UNDE_gene$consexpression2)
    # result$auc$consexpression2[t,j] = DE_Quality_AUC(golden, result$Res$consexpression2[[as.character(t-1)]][[as.character(j-1)]])
    

    result$DE_gene$BaySeq <- result$Res$BaySeq[[as.character(t-1)]][[as.character(j-1)]]
    result$UNDE_gene$BaySeq = setdiff(ngenes,result$DE_gene$BaySeq)
    result$rate1$BaySeq[t,j] = DE_Quality_rate1(golden, result$DE_gene$BaySeq, result$UNDE_gene$BaySeq)
    result$rate2$BaySeq[t,j] = DE_Quality_rate2(golden, result$DE_gene$BaySeq, result$UNDE_gene$BaySeq)
    result$rate3$BaySeq[t,j] = DE_Quality_rate3(golden, result$DE_gene$BaySeq, result$UNDE_gene$BaySeq)
    result$rate4$BaySeq[t,j] = DE_Quality_rate4(golden, result$DE_gene$BaySeq, result$UNDE_gene$BaySeq)
    # result$auc$BaySeq[t, j] = DE_Quality_AUC(golden, result$Res$BaySeq[[as.character(t-1)]][[as.character(j-1)]])
    

    result$DE_gene$consexpression <- result$Res$consexpression[[as.character(t-1)]][[as.character(j-1)]]
    result$UNDE_gene$consexpression = setdiff(ngenes,result$DE_gene$consexpression)
    result$rate1$consexpression[t,j] = DE_Quality_rate1(golden, result$DE_gene$consexpression, result$UNDE_gene$consexpression)
    result$rate2$consexpression[t,j] = DE_Quality_rate2(golden, result$DE_gene$consexpression, result$UNDE_gene$consexpression)
    result$rate3$consexpression[t,j] = DE_Quality_rate3(golden, result$DE_gene$consexpression, result$UNDE_gene$consexpression)
    result$rate4$consexpression[t,j] = DE_Quality_rate4(golden, result$DE_gene$consexpression, result$UNDE_gene$consexpression)
    # result$auc$consexpression[t, j] = DE_Quality_AUC(golden, result$Res$consexpression[[as.character(t-1)]][[as.character(j-1)]])
    
    result$DE_gene$NOISeq <- result$Res$NOISeq[[as.character(t-1)]][[as.character(j-1)]]
    result$UNDE_gene$NOISeq = setdiff(ngenes,result$DE_gene$NOISeq)
    result$rate1$NOISeq[t,j] = DE_Quality_rate1(golden, result$DE_gene$NOISeq, result$UNDE_gene$NOISeq)
    result$rate2$NOISeq[t,j] = DE_Quality_rate2(golden, result$DE_gene$NOISeq, result$UNDE_gene$NOISeq)
    result$rate3$NOISeq[t,j] = DE_Quality_rate3(golden, result$DE_gene$NOISeq, result$UNDE_gene$NOISeq)
    result$rate4$NOISeq[t,j] = DE_Quality_rate4(golden, result$DE_gene$NOISeq, result$UNDE_gene$NOISeq)
    # result$auc$NOISeq[t, j] = DE_Quality_AUC(golden, result$Res$NOISeq[[as.character(t-1)]][[as.character(j-1)]])
    
    
    result$DE_gene$DEGnext <- which(result$Res$DEGnext[[as.character(t-1)]][[as.character(j-1)]]==1)
    result$UNDE_gene$DEGnext = setdiff(ngenes,result$DE_gene$DEGnext)
    result$rate1$DEGnext[t,j] = DE_Quality_rate1(golden, result$DE_gene$DEGnext, result$UNDE_gene$DEGnext)
    result$rate2$DEGnext[t,j] = DE_Quality_rate2(golden, result$DE_gene$DEGnext, result$UNDE_gene$DEGnext)
    result$rate3$DEGnext[t,j] = DE_Quality_rate3(golden, result$DE_gene$DEGnext, result$UNDE_gene$DEGnext)
    result$rate4$DEGnext[t,j] = DE_Quality_rate4(golden, result$DE_gene$DEGnext, result$UNDE_gene$DEGnext)
    # result$auc$DEGnext[t, j] = DE_Quality_AUC(golden, result$Res$DEGnext[[as.character(t-1)]][[as.character(j-1)]])
    
  }
  
}

result_comb[["result"]] = result
saveRDS(result_comb,file="Real_data/MS/result_comb2.rds")

