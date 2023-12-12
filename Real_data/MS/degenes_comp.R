# ******************************************************************************************* #
#                                                                                             #
# ***                   Running DE analysis over Islam et al. dataset                     *** #
#                                                                                             #
# ******************************************************************************************* #

# *** The input file (Islam.txt) includes the count table.
# *** The count table includes 7284 genes and 92 cells.

# *** Set the path:

rm(list= ls())
gc()

source('Real_data/MS/Comp_code.R')

# ******************************************************************************************* #

# *** Read in the raw data:
# IN_VIP = readRDS("Real_data/MS/data_IN_VIP_new.rds")
# data = IN_VIP[["data"]]
# control = c("25G_P3C1", "5546_BA9", "6_C5")
# MS = c("371G_A3D3","7_D5","MS200_A303")
# wsel = c(control, MS)
# a=IN_VIP[["meta"]][["sample"]] %in% control
# b=IN_VIP[["meta"]][["sample"]] %in% MS
# casesamples = IN_VIP[["data"]][,a]
# contsamples = IN_VIP[["data"]][,b]
# counts = cbind(casesamples, contsamples)
# del_cell =NULL
# for (i in 1:dim(counts)[2]){
#   if (length(which(counts[,i]>0))<=10){
#     del_cell = c(del_cell,i)
#   }
# }
# del_gene =NULL
# for (i in 1:dim(counts)[1]){
#   if (length(which(counts[i,]>0))<3){
#     del_gene = c(del_gene,i)
#   }
# }
# # write.csv(del_gene,"sel5.csv")
# genes = setdiff(1:dim(counts)[1],del_gene)
# ex_new = counts[genes,]
# ex_new = as.matrix(counts[genes,])

# 
num = 100
Res = list()
# 


ex_new = readRDS("Real_data/MS/ex_new.rds")
meta_new = readRDS("Real_data/MS/meta.rds")

sig_score = 0.01



result=list()
result$Res=list()


result$Res$DESeq = list()
result$Res$Limma = list()
result$Res$MAST = list()
result$Res$monocle = list()
result$Res$edgeR = list()
result$Res$wilx = list()
result$Res$scde = list()
result$Res$ttest = list()
result$Res$ttestR = list()
result$Res$welch = list()
result$Res$Agg = list()
result$Res$muscat = list()
result$Res$DEsingle = list()
result$Res$zingeR = list()
result$Res$ZIAQ = list()





sel_1 = 1:623
sel_2 = 623+(1:691)#dim(casesamples)[2]+sample(1:dim(contsamples)[2],n2[j])
sel = c(sel_1,sel_2)
ex = ex_new
sml <- c(rep(0,623),rep(1,691))
sml = as.character(sml)
gs <- factor(sml)
me = meta_new[sel,]
inds = table(me$individual)
mc = order(match(me$individual, names(inds)))
ex = ex[,mc]
me = me[mc,]

n1 = 623
n2 = 691

####################################################################################
# Res_ROTS = Run_ROTS(rawData = ex, gr1= n1, gr2=n2, B=1000, K=3000 , FDR = sig_score )
Res_DESeq = Run_DESeq(rawData = ex, gs)
Res_Limma = Run_Limma(rawData = ex, gr1= n1, gr2=n2 )
Res_MAST = Run_MAST(rawData = ex, gs)
Res_monocle = Run_monocle(rawData = ex, gr1= n1, gr2=n2, gs )
Res_edgeR = Run_edgeR(rawData = ex, gr1= n1, gr2=n2, gs )
Res_wilx = Run_wilx(rawData = ex, gr1= n1, gr2=n2)
Res_scde = Run_scde(rawData = ex, gs)
Res_ttest = Run_ttest(rawData = ex, gr1= n1, gr2=n2)
Res_ttestR = Run_ttestR(rawData = ex, gr1= n1, gr2=n2)
Res_welch = Run_welch(rawData = ex, gr1= n1, gr2=n2)

Res_Agg = tryCatch(Run_Agg(rawdata = ex,metas = me), error = function(e) {NA})
Res_muscat = tryCatch(Run_muscat(rawdata = ex, metas = me), error = function(e) {NA})

Data <- data_process(Data = ex, group = sml, is.normalized = F)
Res_zingeR <- Execute_zingeR.edgeR(Data,  maxit.EM = 200)
Res_DEsingle <- Execute_DEsingle(object = Data, DEsingle.parallel = F)

res = ziaq(ex, colDat, formula = ~ condition,
           group = 'condition', probs = c(0.25, 0.5, 0.75),
           log_i = TRUE, parallel = FALSE, no.core = 1)
Res_ZIAQ = as.matrix(res[["pvalue"]])





result$Res$DESeq = Res_DESeq
result$Res$Limma = Res_Limma
result$Res$MAST = Res_MAST
result$Res$monocle = Res_monocle
result$Res$edgeR = Res_edgeR
result$Res$wilx = Res_wilx
result$Res$scde = Res_scde
result$Res$ttest = Res_ttest
result$Res$ttestR = Res_ttestR
result$Res$welch = Res_welch 
result$Res$Agg = Res_Agg
result$Res$muscat = Res_muscat
result$Res$DEsingle = Res_DEsingle    
result$Res$zingeR = Res_zingeR
result$Res$ZIAQ = Res_ZIAQ
  

ngenes = length(Res_Limma)
deg = list()
undeg = list()
for (d in 1:15){
  deg[[d]] = which(result[["Res"]][[d]]<=sig_score)
  undeg[[d]] = setdiff(ngenes, deg[[d]])
}
rem_id = NULL
for (i in 1:15){
  if (length(deg[[i]]) < 2000 | length(deg[[i]]) > 10000)  
    rem_id = c(rem_id, i)
}
d = 1:15; 
d = d[-rem_id]
de_union = NULL
de_inter = deg[[1]]
for (e in d){
  de_union = union(de_union, deg[[e]])
  de_inter = intersect(de_inter, deg[[e]])
}
de_not_inter = setdiff(ngenes, de_inter)
golden = list(notDE=de_not_inter, DE=de_inter)

degene_comp = list(degene = result$Res, golden = golden)

saveRDS(degene_comp,"Real_data/MS/degene_comp.rds")  
  
  
  


