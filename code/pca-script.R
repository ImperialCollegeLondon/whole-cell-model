library(ggplot2)
library(reshape2)
library(gridExtra)
# library(cowplot)
library(R.utils)
library(ade4)
library(diffusionMap)



the_good<- read.csv("the-good.csv", row.names = 1)
the_good$ID<- "g"
the_bad<- read.csv("the-bad.csv", row.names = 1)
the_bad$ID<-"b"
the_ugly<- read.csv("the-ugly.csv", row.names = 1)
the_ugly$ID<- "b" #can be its own factor if 'u' or can be collapsed into bad if 'b'

all_data<- rbind(the_good, the_bad, the_ugly)
log_data<- log10(all_data[,1:6])
log_data$ID<- all_data$ID

cat_data<-log_data #changing 0.1 to low, 10 to med and 1000 to high
# cat_data[,1:6]<- as.numeric(cat_data[,1:6])
cat_data[cat_data == '-1']<- 'A'
cat_data[cat_data == '1']<- 'B'
cat_data[cat_data == '3']<- 'C'
cat_data$k_cat_AA <- as.factor(cat_data$k_cat_AA)
cat_data$k_a_NH4 <- as.factor(cat_data$k_a_NH4)
cat_data$k_NH4 <- as.factor(cat_data$k_NH4)
cat_data$k_a_AA <- as.factor(cat_data$k_a_AA)
cat_data$k_NH4_AA <- as.factor(cat_data$k_NH4_AA)
cat_data$k_ribo_a <- as.factor(cat_data$k_ribo_a)
cat_data$ID <- as.factor(cat_data$ID)
#change the characters to factors


# melted_data<- melt(sorted_data, id.var="timestep")
#can change data to long format if needed



mca2= dudi.acm(cat_data[,1:6], scannf = F, nf = 6)

a = rep("green",sum(all_data$ID=="g"))
b = rep("red",sum(all_data$ID=="b"))
colr2 = c(a,b)
eig = 100*mca1$eig/sum(mca1$eig)
x1 = paste("PC1",round(eig[1],2),"%")
x2 = paste("PC2",round(eig[2],2),"%")

plot(mca2$li[,1], mca2$li[,2], col = colr2,xlab = x1, ylab = x2)
plot(mca2$li[,2], mca2$li[,3], col = colr2,xlab = x1, ylab = x2)




# pca1= dudi.pca(df = log_data[,1:6], center = T, scale = T)
# eig = 100*pca1$eig/sum(pca1$eig)
# x1 = paste("PC1",round(eig[1],2),"%")
# x2 = paste("PC2",round(eig[2],2),"%")
# 
# a = rep("green",sum(all_data$ID=="g"))
# b = rep("red",sum(all_data$ID=="b"))
# c = rep("blue", sum(all_data$ID =="u"))
# colr2 = c(a,b,c)
# 
# plot(pca1$li[,1], pca1$li[,2], xlab = x1, ylab = x2, col = colr2, pch = 16, cex = 0.8)
# plot(pca1$li[,1], pca1$li[,3], xlab = x1, ylab = x2, col = colr2, pch = 16, cex = 0.8)
# plot(pca1$li[,3], pca1$li[,2], xlab = x1, ylab = x2, col = colr2, pch = 16, cex = 0.8)
#######################################################################################
##trying a diffusion map for dimensionality reduction
# 
# D = dist(log_data, method = "maximum")
# # D=dist(data_AG140_RCC)
# dmap = diffuse(D,maxdim  = 6)
# plot.dmap(dmap,col=colr2)

#######################################################################################
#trying a PCA method for categorical data
library(FactoMineR)
library(factoextra)

# mca1<- MCA(cat_data[,1:6])#, quali.sup = cat_data$ID)
mca1<- MCA(cat_data)
plot(mca1, col.ind = cat_data$ID)
fviz_screeplot(mca1, addlabels = T )


var <- get_mca_var(mca1)
fviz_mca_var(mca1, repel = T)#repel means the labels wont overlap but this takes a while to plot 
# fviz_mca_biplot(mca1, col.ind = cat_data$ID, axes = c(1,2)) #the default plot the two most informative dimensions
# fviz_mca_biplot(mca1, col.ind = cat_data$ID, axes = c(1,3))
# fviz_mca_biplot(mca1, col.ind = cat_data$ID, axes = c(1,4))
# fviz_mca_biplot(mca1, col.ind = cat_data$ID, axes = c(1,5))
fviz_mca_biplot(mca1, col.ind = cat_data$ID)#, axes = c(2,3))
fviz_mca_ind(mca1, label = "none", habillage = cat_data$ID, addEllipses = T,  axes = c(1,2))#ellipse.type = "confidence",)

# fviz_mca_biplot(mca1, col.ind = cat_data$ID, axes = c(2,4))
# fviz_mca_biplot(mca1, col.ind = cat_data$ID, axes = c(2,5))
# fviz_mca_biplot(mca1, col.ind = cat_data$ID, axes = c(3,4))
# fviz_mca_biplot(mca1, col.ind = cat_data$ID, axes = c(3,5))
# fviz_mca_biplot(mca1, col.ind = cat_data$ID, axes = c(4,5))

fviz_mca_var(mca1, col.var = "cos2", gradient.cols = c("red", "green") , axes = c(2,5))
fviz_cos2(mca1, choice = "var", axes = 2:3)#shows how well represented each factor is by the two dimenstions indicated. if there are
#very long tail with very low values then those are not well represented by the two dimensions so a higher dimesional solution is needed

# round(var$contrib,2)
fviz_contrib(mca1, choice = "var", axes = 1:5)

ind<- get_mca_ind(mca1)
# head(ind$contrib)
fviz_contrib(mca1, choice = "ind", axes = 1:5, top = 20)

result_des<- dimdesc(mca1, axes = c(1,2))
result_des[[2]]
#######################################################################################
library(ExPosition)

mca3 <- epMCA(cat_data[,1:6], graph = FALSE, correction = "bg")
plot(mca3$Plotting.Data$fj.pch)



