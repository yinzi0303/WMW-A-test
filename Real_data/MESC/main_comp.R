
source('Real_data/MESC/Comp_code.R')

# ******************************************************************************************* #

# *** Read in the raw data:
dataset = readRDS("Real_data/MESC/data_count_MESC.rds")
ex_new = dataset[['expr']]
label = dataset$label

casesamples = which(label==1)
contsamples = which(label==2)

sig_score = 0.01
result = list()

n1=c(5,5,5,5,10,15,20,25)
n2=c(5,10,15,20,10,15,20,25)

for (t in 1:num){
  
  result$Res$monocle[[t]] = list()
  result$Res$scde[[t]] = list()
  result$Res$MAST[[t]] = list()
  result$Res$DESeq[[t]] = list()
  result$Res$Limma[[t]] = list()
  result$Res$edgeR[[t]] = list()
  result$Res$ttest[[t]] = list()
  result$Res$ttestR[[t]] = list()
  result$Res$welch[[t]] = list()
  result$Res$wilx[[t]] = list()
  result$Res$DEsingle[[t]] = list()
  result$Res$zingeR[[t]] = list()
  result$Res$ZIAQ[[t]] = list()
  

  #length(n1)
  for (j in 1:length(n1)){
    sel_1 = sample(1:length(casesamples),n1[j])
    sel_2 = length(casesamples)+sample(1:length(contsamples),n2[j])
    sel = c(sel_1,sel_2)
    ex = ex_new[,sel]
    sml <- c(rep(0,n1[j]),rep(1,n2[j]))
    sml = as.character(sml)
    gs <- factor(sml)
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
    result$Res$DEsingle[[t]][[j]] = Res_DEsingle    
    result$Res$zingeR[[t]][[j]] = Res_zingeR
    result$Res$ZIAQ[[t]][[j]] = Res_ZIAQ
    
  }
  

  saveRDS(result,"Real_data/MESC/compresult.rds")
  
}
    

