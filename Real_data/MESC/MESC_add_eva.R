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



# result = readRDS("Real_data/MESC/compresult.rds")
result<-readRDS("Real_data/MESC/result_comb.rds") 
ntimes = 100

golden = NULL
deg = read.csv("Real_data/MESC/deg.csv",header = F)
de = match(deg[,1], rownames(ex_new))
sigde = na.omit(de)
golden$DE = sigde
golden$notDE = setdiff(1:18370, sigde)


ntest = 8
sig_score = 0.01
ngenes = 1:length(result[["Res"]][["welch"]][[1]][[1]])



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
mat_MESC<-np$load(paste0(path, "consexpression-master/results/MESC/MESC_deg.npy"), allow_pickle=TRUE) 
result[["Res"]][["consexpression"]] = mat_MESC[[1]][["Consexp"]]
# result[["Res"]][["consexpression2"]] = mat_MESC[[1]][["Consexp_2"]]
result[["Res"]][["BaySeq"]] = mat_MESC[[1]][["BaySeq"]]
result[["Res"]][["NOISeq"]] = mat_MESC[[1]][["NOISeq"]]

mat_MESC2<-np$load(paste0(path, "DEGnext_code/results/MESC/MESC_deg3.npy"), allow_pickle=TRUE) 
result[["Res"]][["DEGnext"]] = mat_MESC2[[1]]

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



saveRDS(result,file="Real_data/MESC/result_comb.rds")

# saveRDS(result,file="Real_data/MESC/compresult_eva.rds")

