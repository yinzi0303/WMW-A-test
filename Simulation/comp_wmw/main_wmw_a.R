

rm(list= ls())
gc()


source("Simulation/comp_wmw/WMWA.R")
source('Simulation/comp_wmw/z_generate.R')


comp_wmw <- function(x, y){
  ex_xy = c(x, y)
  # ex_xy = apply(ex_xy, 2, as.numeric)
  
  result = list()
  result$Res$WMWANL = list()
  result$Res$WMWANLM = list()
  result$Res$WMWALN = list()
  result$Res$WMWALNM = list()
  result$Res$WMWANB = list()
  result$Res$WMWANBM = list()
  
  result$Res$WMW = list()
  
  k = 5
  m = length(x)
  n = length(y)
  num_l = length(ex_xy) * k
  
  ####################################################################################
  
  ex_z1 = z_generateNL(comb_xy = ex_xy, k=k, num_l=num_l)
  ex_z1[which(is.na(ex_z1))]=0
  ex = c(ex_xy, ex_z1)
  N = length(ex)
  
  Res_WMWANL = wmwa(ex, m=m, n=n, N=N, nperm=10000)
  result$Res$WMWANL = Res_WMWANL
  
  ####################################################################################
  
  ex_z2 = z_generateLN(comb_xy = ex_xy, k=k, num_l=num_l)
  ex_z2[which(is.na(ex_z2))]=0
  ex = c(ex_xy, ex_z2)
  N = length(ex)
  
  Res_WMWALN = wmwa(ex, m=m, n=n, N=N, nperm=10000)
  result$Res$WMWALN = Res_WMWALN
  
  ####################################################################################
  
  ex_z3 = z_generateNLM(comb_xy = ex_xy, k=k, num_l=num_l)
  ex_z3[which(is.na(ex_z3))]=0
  ex = c(ex_xy, ex_z3)
  N = length(ex)
  
  Res_WMWANLM = wmwa(ex, m=m, n=n, N=N, nperm=10000)
  result$Res$WMWANLM = Res_WMWANLM
  
  ####################################################################################
  
  ex_z4 = z_generateLNM(comb_xy = ex_xy, k=k, num_l=num_l)
  ex_z4[which(is.na(ex_z4))]=0
  ex = c(ex_xy, ex_z4)
  N = length(ex)
  
  Res_WMWALNM = wmwa(ex, m=m, n=n, N=N, nperm=10000)
  result$Res$WMWALNM = Res_WMWALNM
  
  
  
  ####################################################################################
  
  ex_z5 = z_generateNB(comb_xy = ex_xy, k=k, num_l=num_l)
  ex_z5[which(is.na(ex_z5))]=0
  ex = c(ex_xy, ex_z5)
  N = length(ex)
  
  Res_WMWANB = wmwa(ex, m=m, n=n, N=N, nperm=10000)
  result$Res$WMWANB = Res_WMWANB
  
  ####################################################################################
  
  ex_z6 = z_generateNBM(comb_xy = ex_xy, k=k, num_l=num_l)
  ex_z6[which(is.na(ex_z6))]=0
  ex = c(ex_xy, ex_z6)
  N = length(ex)
  
  Res_WMWANBM = wmwa(ex, m=m, n=n, N=N, nperm=10000)
  result$Res$WMWANBM = Res_WMWANBM
  
  Res_WMW = wilcox.test(x, y)$p.value
  result$Res$WMW = Res_WMW
  
  return(result)
}




# ******************************************************************************************* #

# *** Read in the raw data:

x = c(1,2)
y = c(7,8)
res_sim1 = comp_wmw(x, y)

x = c(1:4)
y = c(3:6)
res_sim2 = comp_wmw(x, y)

x = c(1:4,5)
y = c(3:6,5)
res_sim3 = comp_wmw(x, y)

x = c(1:4,7)
y = c(6:9,3)
res_sim4 = comp_wmw(x, y)



saveRDS(res_sim1, 'Simulation/comp_wmw/sim1.rds')
saveRDS(res_sim2, 'Simulation/comp_wmw/sim2.rds')
saveRDS(res_sim3, 'Simulation/comp_wmw/sim3.rds')
saveRDS(res_sim4, 'Simulation/comp_wmw/sim4.rds')



x = c(1,2)
y = c(7,8)
res_sim1 = wilcox.test(x,y)$p.value

x = c(1:4)
y = c(3:6)
res_sim2 = wilcox.test(x,y)$p.value

x = c(1:4,5)
y = c(3:6,5)
res_sim3 = wilcox.test(x,y)$p.value

x = c(1:4,7)
y = c(6:9,3)
res_sim4 = wilcox.test(x,y)$p.value

res_sim_wmw = list(res_sim1=res_sim1,res_sim2=res_sim2,
                   res_sim3=res_sim3,res_sim4=res_sim4)

saveRDS(res_sim_wmw, 'Simulation/comp_wmw/sim_wmw.rds')

####################################################

res_sim1 = readRDS('Simulation/comp_wmw/sim1.rds')
res_sim2 = readRDS('Simulation/comp_wmw/sim2.rds')
res_sim3 = readRDS('Simulation/comp_wmw/sim3.rds')
res_sim4 = readRDS('Simulation/comp_wmw/sim4.rds')
res_sim_wmw = readRDS('Simulation/comp_wmw/sim_wmw.rds')

res_wmw = unlist(res_sim_wmw)


combined_matrix <- do.call(rbind, c(res_sim1, res_sim2, res_sim3, res_sim4))
print(combined_matrix)
write.csv(combined_matrix, 'Simulation/comp_wmw/combined_matrix.csv')




# 创建示例向量
vector1 <- c(1, 2, 3, 4, 5)
vector2 <- c(2, 3, 3, 4, 5)
vector <- unique(c(vector1, vector2))
a=as.data.frame(table(vector1))
b=as.data.frame(table(vector2))

sa = setdiff(as.factor(vector), a$vector1)
sb = setdiff(as.factor(vector), b$vector2)

if(length(sa)>0){
  new_row <- data.frame(vector1 = sa, Freq = 0)
  an <- rbind(a, new_row)
}else{an = a}

if(length(sb)>0){
  new_row <- data.frame(vector2 = sb, Freq = 0)
  bn <- rbind(b, new_row)
}else{bn=b}

barplot(an, col=rgb(0, 0, 1, alpha=0.5), xlab="Value")
par(new=TRUE)
barplot(bn, col=rgb(1, 0, 0, alpha=0.5), xlab="Value")

hist(vector1, col=rgb(0, 0, 1, alpha=0.5), xlab="Value", ylab="Frequency", main="Combined Histogram", xlim=c(min(vector1, vector2), max(vector1, vector2)))
par(new=TRUE)  # 允许在同一绘图区域上叠加绘图
hist(vector2, col=rgb(1, 0, 0, alpha=0.5), xlab="", ylab="", main="Combined Histogram", axes=FALSE)

# 添加图例
legend("topright", legend=c("Vector 1", "Vector 2"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))




# 创建示例向量
vector1 <- c(1, 2, 3, 4, 5)
vector2 <- c(2, 3, 3, 4, 5)

# 创建直方图
hist(vector1, col=rgb(0, 0, 1, alpha=0.5), xlab="Value", ylab="Frequency", main="Combined Histogram", xlim=c(min(vector1, vector2), max(vector1, vector2)))
par(new=TRUE)  # 允许在同一绘图区域上叠加绘图
hist(vector2, col=rgb(1, 0, 0, alpha=0.5), xlab="", ylab="", main="Combined Histogram", axes=FALSE)

# 添加图例
legend("topright", legend=c("Vector 1", "Vector 2"), fill=c(rgb(0, 0, 1, alpha=0.5), rgb(1, 0, 0, alpha=0.5)))




carrots <- data.frame(length = rnorm(100000, 6, 2))
cukes <- data.frame(length = rnorm(50000, 7, 2.5))
carrots$veg <- 'carrot'
cukes$veg <- 'cuke'
vegLengths <- rbind(carrots, cukes)
ggplot(vegLengths, aes(length, fill = veg)) + geom_density(alpha = 0.2)
ggplot(vegLengths, aes(length, fill = veg)) +
  geom_histogram(alpha = 0.5, aes(y = ..density..), position = 'identity')

x = c(1:4,7)
y = c(6:9,3)
x = data.frame(value=x)
y = data.frame(value=y)
x$veg <- 'x'
y$veg <- 'y'
z <- rbind(x,y)
ggplot(z, aes(value, fill = veg)) + geom_density(alpha = 0.2)
ggplot(z, aes(value, fill = veg)) +
  geom_histogram(alpha = 0.5, position = 'identity', binwidth=0.5) #aes(y = ..density..),


