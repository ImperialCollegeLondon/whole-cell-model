# rm(list=ls())
# graphics.off()

library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(R.utils)

#pick a few points rom green cloud and plot all odes manually. try to find param combinations where reasonable end points. 
#use thresholds. try to increase k_cat_AA as wel as increasing k_a to keep ratio roughly the same


important_data <- read.csv("param-sweep2-results.csv")
important_data<- important_data[-1]
# important_data$log_AA <- log(important_data$AminoAcid)+0.1
# important_data$log_atp <- log(important_data$ATP)+0.1

non_neg_data <- subset(important_data, important_data$AminoAcid >= 0 & important_data$ATP >= 0)
non_neg_data$log_AA<- log(non_neg_data$AminoAcid)
non_neg_data$log_ATP<- log(non_neg_data$ATP)
non_neg_data$product<- as.numeric(as.character(non_neg_data[,3]))* as.numeric(as.character(non_neg_data[,4]))* as.numeric(as.character(non_neg_data[,5]))* as.numeric(as.character(non_neg_data[,6]))* as.numeric(as.character(non_neg_data[,7]))* as.numeric(as.character(non_neg_data[,8]))* as.numeric(as.character(non_neg_data[,9]))
non_neg_data$product_col<- log(non_neg_data$product)
non_neg_data$product_col<- non_neg_data$product_col - min(non_neg_data$product_col)

plt <- ggplot(data = non_neg_data, aes(x = k_ribo_AA_a, y = k_ribo_a_AA, colour = AminoAcid)) +geom_point() +facet_grid(rows = non_neg_data$k_ribo_a)
plt
# plot(x = important_data$k_ribo_AA_a, y = important_data$k_ribo_a_AA, col = important_data$AminoAcid)
plot(x = non_neg_data$AminoAcid, y = non_neg_data$ATP)
plot(x = non_neg_data$log_AA, y = non_neg_data$log_ATP, col = non_neg_data$k_ribo_a_AA, main = "colored by k_NH4_AA")
plot(x = non_neg_data$log_AA, y = non_neg_data$log_ATP, col = non_neg_data$k_a, main = "colored by k_a")



df1<- subset(non_neg_data, non_neg_data$k_a== '0.1')
df2<- subset(non_neg_data, non_neg_data$k_a== '10.0')
df3<- subset(non_neg_data, non_neg_data$k_a== '1000.0')
plot(x = df1$log_AA, y = df1$log_ATP, main = "k_a = 0.1")
plot(x = df2$log_AA, y = df2$log_ATP, main = "k_a = 10")
plot(x = df3$log_AA, y = df3$log_ATP, main = "k_a = 1000", col = df3$k_ribo_AA_a)


hist(x = non_neg_data$ATP)
hist(x = non_neg_data$AminoAcid)

boxplot(non_neg_data$log_AA ~ non_neg_data$log_ATP)



# require(graphics); require(grDevices)
# x  <- as.matrix(mtcars)
# rc <- rainbow(nrow(x), start = 0, end = .3)
# cc <- rainbow(ncol(x), start = 0, end = .3)
# hv <- heatmap(x, col = cm.colors(256), scale = "column",
#               RowSideColors = rc, ColSideColors = cc, margins = c(5,10),
#               xlab = "specification variables", ylab =  "Car Models",
#               main = "heatmap(<Mtcars data>, ..., scale = \"column\")")
# utils::str(hv) # the two re-ordering index vectors


