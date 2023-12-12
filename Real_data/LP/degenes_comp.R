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

source('Real_data/LP/Comp_code.R')


# ******************************************************************************************* #

# *** Read in the raw data:


dataset = readRDS('Real_data/LP/GSE54456.rds') 

data = dataset[["expr"]]
label = dataset$label



num = 100
Res = list()


sig_score = 0.01
result = list()



casesamples = which(label==1)
contsamples = which(label==2)


# result$Res$ROTS = list()
# result$Res$DESeq = list()
result$Res$Limma = list()
# result$Res$MAST = list()
result$Res$monocle = list()
# result$Res$edgeR = list()
result$Res$wilx = list()
# result$Res$scde = list()
result$Res$ttest = list()
result$Res$ttestR = list()
result$Res$welch = list()




n1 = length(casesamples)
n2 = length(contsamples)
sel = c(casesamples,contsamples)
ex = data[,sel]

sml <- c(rep(0,n1),rep(1,n2))
sml = as.character(sml)
gs <- factor(sml)


####################################################################################
# Res_ROTS = Run_ROTS(rawData = ex, gr1= n1, gr2=n2, B=1000, K=3000 , FDR = sig_score )
# Res_DESeq = Run_DESeq(rawData = ex, gs)
Res_Limma = Run_Limma(rawData = ex, gs )
# Res_MAST = Run_MAST(rawData = ex, gs)
Res_monocle = Run_monocle(rawData = ex, gr1= n1, gr2=n2, gs )
# Res_edgeR = Run_edgeR(rawData = ex, gr1= n1, gr2=n2, gs )
Res_wilx = Run_wilx(rawData = ex, gr1= n1, gr2=n2)
# Res_scde = Run_scde(rawData = ex, gs)
Res_ttest = Run_ttest(rawData = ex, gr1= n1, gr2=n2)
Res_ttestR = Run_ttestR(rawData = ex, gr1= n1, gr2=n2)
Res_welch = Run_welch(rawData = ex, gr1= n1, gr2=n2)

individual = colnames(ex)
cluster = c(rep("spike",n1[j]+n2[j]))
diagnosis = c(rep("healthy",n1[j]),rep("disease",n2[j]))
metas = data.frame(individual=individual, diagnosis=diagnosis, cluster=cluster)
Data <- data_process(Data = ex, group = sml, is.normalized = T)

Res_zingeR <- Execute_zingeR.edgeR(Data,  maxit.EM = 200)
Res_DEsingle <- Execute_DEsingle(object = Data, DEsingle.parallel = F)
# Res_Agg = tryCatch(Run_Agg(rawdata = ex,metas = metas), error = function(e) {NA})
Res_muscat = tryCatch(Run_muscat(rawdata = ex, metas = metas), error = function(e) {NA})


# result$Res$ROTS = NULL
# result$Res$DESeq = Res_DESeq
result$Res$Limma = Res_Limma
# result$Res$MAST = Res_MAST
result$Res$monocle = Res_monocle
# result$Res$edgeR = Res_edgeR
result$Res$wilx = Res_wilx
# result$Res$scde = NULL
result$Res$ttest = Res_ttest
result$Res$ttestR = Res_ttestR
result$Res$welch = Res_welch
result$Res$zingeR = Res_zingeR
result$Res$DEsingle = Res_DEsingle
result$Res$muscat = Res_muscat

ngenes = length(Res_Limma)
deg = list()
undeg = list()
for (d in 1:9){
  deg[[d]] = which(result[["Res"]][[d]]<=sig_score)
  undeg[[d]] = setdiff(ngenes, deg[[d]])
}
d = 1:9 
de_union = NULL
de_inter = deg[[1]]
for (e in d){
  de_inter = intersect(de_inter, deg[[e]])
}
de_not_inter = setdiff(ngenes, de_inter)
golden = list(notDE=de_not_inter, DE=de_inter)

degene_comp = list(degene = result$Res, golden = golden)

saveRDS(degene_comp,"Real_data/LP/degene_comp.rds")  
  
 
