rm(list=ls())
graphics.off()

library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)
library(R.utils)


important_data <- read.csv("param-sweep-results.csv")
important_data<- important_data[-1]
# important_data$log_AA <- log(important_data$AminoAcid)+0.1
# important_data$log_atp <- log(important_data$ATP)+0.1

non_neg_data <- subset(important_data, important_data$AminoAcid >= 0 & important_data$ATP >= 0)
non_neg_data$log_AA<- log(non_neg_data$AminoAcid)
non_neg_data$log_ATP<- log(non_neg_data$ATP)

plt <- ggplot(data = non_neg_data, aes(x = k_ribo_AA_a, y = k_ribo_a_AA, colour = AminoAcid)) +geom_point() +facet_grid(rows = non_neg_data$k_ribo_a)
plt
# plot(x = important_data$k_ribo_AA_a, y = important_data$k_ribo_a_AA, col = important_data$AminoAcid)
plot(x = non_neg_data$AminoAcid, y = non_neg_data$ATP)
plot(x = non_neg_data$log_AA, y = non_neg_data$log_ATP)
hist(x = non_neg_data$log_ATP)
hist(x = non_neg_data$log_AA)



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


