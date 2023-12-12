

rm(list= ls())
gc()





source('Real_data/MS/Comp_code.R')
meta_new = readRDS("Real_data/MS/meta.rds")


# del_cell =NULL
# for (i in 1:dim(ex_new)[2]){
#   if (length(which(ex_new[,i]>0))<=10){
#     del_cell = c(del_cell,i)
#   }
# }
# del_gene =NULL
# for (i in 1:dim(ex_new)[1]){
#   if (length(which(ex_new[i,]>0))<3){
#     del_gene = c(del_gene,i)
#   }
# }
# genes = setdiff(1:dim(ex_new)[1],del_gene)
# ex_new = as.matrix(ex_new[genes,])

ex_new = readRDS("ex_new.rds")

num = 100



sig_score = 0.01
result = list()
result$Res=list()



for (t in 1:num){
  
  
  result$Res$DESeq[[t]] = list()
  result$Res$Limma[[t]] = list()
  result$Res$MAST[[t]] = list()
  result$Res$monocle[[t]] = list()
  result$Res$edgeR[[t]] = list()
  result$Res$wilx[[t]] = list()
  result$Res$scde[[t]] = list()
  result$Res$ttest[[t]] = list()
  result$Res$ttestR[[t]] = list()
  result$Res$welch[[t]] = list()
  result$Res$Agg[[t]] = list()
  result$Res$muscat[[t]] = list()
  result$Res$DEsingle[[t]] = list()
  result$Res$zingeR[[t]] = list()
  result$Res$ZIAQ[[t]] = list()
  
  
  
  n1=c(5,5,5,5,10,15,20,25)
  n2=c(5,10,15,20,10,15,20,25)
  #length(n1)
  for (j in 1:length(n1)){
    sel_1 = sample(1:623,n1[j])#sample(1:dim(casesamples)[2],n1[j])
    sel_2 = 623+sample(1:691,n2[j])#dim(casesamples)[2]+sample(1:dim(contsamples)[2],n2[j])
    sel = c(sel_1,sel_2)
    ex = ex_new[,sel]
    sml <- c(rep(0,n1[j]),rep(1,n2[j]))
    gs <- factor(as.character(sml))
    colDat = data.frame(condition = sml)
    me = meta_new[sel,]
    inds = table(me$individual)
    mc = order(match(me$individual, names(inds)))
    ex = ex[,mc]
    me = me[mc,]
    
    
    ####################################################################################
    Res_DESeq = Run_DESeq(rawData = ex, gs)
    Res_Limma = Run_Limma(rawData = ex, gr1= n1[j], gr2=n2[j] )
    Res_MAST = Run_MAST(rawData = ex, gs)
    Res_monocle = Run_monocle(rawData = ex, gr1= n1[j], gr2=n2[j], gs )
    Res_edgeR = Run_edgeR(rawData = ex, gr1= n1[j], gr2=n2[j], gs )
    Res_wilx = Run_wilx(rawData = ex, gr1= n1[j], gr2=n2[j])
    Res_scde = Run_scde(rawData = ex, gs)
    Res_ttest = Run_ttest(rawData = ex, gr1= n1[j], gr2=n2[j])
    Res_ttestR = Run_ttestR(rawData = ex, gr1= n1[j], gr2=n2[j])
    Res_welch = Run_welch(rawData = ex, gr1= n1[j], gr2=n2[j])
    Res_Agg = tryCatch(Run_Agg(rawdata = ex,metas = me), error = function(e) {NA})
    Res_muscat = tryCatch(Run_muscat(rawdata = ex, metas = me), error = function(e) {NA})
    
    Data <- data_process(Data = ex, group = sml, is.normalized = F)
    Res_zingeR <- Execute_zingeR.edgeR(Data,  maxit.EM = 200)
    Res_DEsingle <- Execute_DEsingle(object = Data, DEsingle.parallel = F)
    
    res = ziaq(ex, colDat, formula = ~ condition,
               group = 'condition', probs = c(0.25, 0.5, 0.75),
               log_i = TRUE, parallel = FALSE, no.core = 1)
    Res_ZIAQ = as.matrix(res[["pvalue"]])

   
    
    result$Res$DESeq[[t]][[j]] = Res_DESeq
    result$Res$Limma[[t]][[j]] = Res_Limma
    result$Res$MAST[[t]][[j]] = Res_MAST
    result$Res$monocle[[t]][[j]] = Res_monocle
    result$Res$edgeR[[t]][[j]] = Res_edgeR
    result$Res$wilx[[t]][[j]] = Res_wilx
    result$Res$scde[[t]][[j]] = Res_scde
    result$Res$ttest[[t]][[j]] = Res_ttest
    result$Res$ttestR[[t]][[j]] = Res_ttestR
    result$Res$welch[[t]][[j]] = Res_welch
    result$Res$Agg[[t]][[j]] = Res_Agg
    result$Res$muscat[[t]][[j]] = Res_muscat
    result$Res$DEsingle[[t]][[j]] = Res_DEsingle    
    result$Res$zingeR[[t]][[j]] = Res_zingeR
    result$Res$ZIAQ[[t]][[j]] = Res_ZIAQ
    
  }
  

  saveRDS(result,"Real_data/MS/compresult.rds")
  
}


