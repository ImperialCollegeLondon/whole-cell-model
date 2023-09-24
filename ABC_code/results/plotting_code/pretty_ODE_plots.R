# rm(list=ls())
# graphics.off()

library(ggplot2)
library(reshape2)
# library(gridExtra)
# library(cowplot)
# library(R.utils)

#pick a few points rom green cloud and plot all odes manually. try to find param combinations where reasonable end points. 
#use thresholds. try to increase k_cat_AA as wel as increasing k_a to keep ratio roughly the same


# important_data <- read.csv("param-sweep5-results.csv")
# raw_data <- read.csv("../data/0.1/param-sweep-0.1-0.1-0.1-0.1-1000.0-0.1.csv")#oscillation
# raw_data <- read.csv("../data/param-sweep-10.0-10.0-10.0-10.0-10.0-10.0.csv")
# raw_data<- read.csv("../data/latest-output.csv")
# raw_data<- read.csv("multilevel-output.csv")
scenarios <- list("norm_C_N","C_run_out","N_run_out","C_N_run_out","C_N_low")


for (scenario in scenarios){
  filepath <- paste("../ATP_adjusted/",scenario,"_index-10_OG_model.csv", sep = "")
  raw_data<- read.csv(filepath)
  raw_data$"timestep" = seq.int(nrow(raw_data))
  melted_data = melt(raw_data, id.var = "timestep")
  plt1 = ggplot(data = melted_data, aes(x = timestep, y = log10(value+0.1)))
  plt1 = plt1 + labs(y="Log Molecules", x = "Time", title = paste(scenario,"no AA")) + theme_bw() +theme(legend.position = "none", plot.title = element_text(hjust =0.5))
  plt1 = plt1 +annotate("rect",xmin = -Inf,xmax = 10000, ymin = -Inf, ymax = Inf,fill = "lightgray",alpha = 0.5)
  plt1 = plt1 + annotate("rect",xmin = 20000,xmax = Inf,ymin = -Inf, ymax = Inf,fill = "lightgrey",alpha = 0.5)
  plt1 = plt1 +geom_line(aes(color = variable), size = 2)+ facet_wrap(~ variable, ncol = 4)
  savefile1 = paste(scenario,"_index-10_OG_model.png" ,sep = "")
  ggsave(filename =  savefile1)
  # plt1
  
  
  filepath <- paste("../ATP_adjusted/",scenario,"_index-10_AA_model.csv", sep = "")
  raw_data<- read.csv(filepath)
  raw_data$"timestep" = seq.int(nrow(raw_data))
  melted_data = melt(raw_data, id.var = "timestep")
  plt2 = ggplot(data = melted_data, aes(x = timestep, y = log10(value+0.1))) + geom_rect(aes(xmin = 0,xmax = 10000,ymin = -Inf, ymax = Inf),fill = "lightgrey",alpha = 0.2)
  plt2 = plt2+ geom_rect(aes(xmin = 20000,xmax = 30000,ymin = -Inf, ymax = Inf),fill = "lightgrey",alpha = 0.2 ) +theme_bw() 
  plt2 = plt2 +theme(legend.position = "none", plot.title = element_text(hjust =0.5))+ labs(y="Log Molecules", x = "Time", title = paste(scenario,"with AA"))+geom_line(aes(color = variable), size = 2)+ facet_wrap(~ variable, ncol = 4)
  # plt2
  savefile2 = paste(scenario,"_index-10_AA_model.png" ,sep = "")
  ggsave(filename =  savefile2)
  
  
  # filepath <- paste("../ATP_adjusted/",scenario,"index-10_AA_model.csv", sep = "")
  # raw_data<- read.csv(filepath)
  # raw_data$"timestep" = seq.int(nrow(raw_data))
  # melted_data = melt(raw_data, id.var = "timestep")
  # plt2 = ggplot(data = melted_data, aes(x = timestep, y = log10(value+0.1))) +theme_bw() +theme(legend.position = "none")+ labs(y="Log Molecules", x = "Time", title = "" )+geom_line(aes(color = variable), size = 2)+ facet_wrap(~ variable, ncol = 4)
  # plt2

  
}


# important_data$un_ov<- 
# important_data$log_AA <- log(important_data$AminoAcid)+0.1
# important_data$log_atp <- log(important_data$ATP)+0.1
# plot(log10(sorted_again$`22AA`), log10(sorted_again$`21ATP`))#, col = sorted_again$timestep)
# plot(log10(important_data$AminoAcid), log10(important_data$ATP))
par(mfrow= c(2,1))
# times<- c(1,1000)
plot(log10(sorted_again$`22AA`), x = sorted_again$timestep)#, xlim = times, main = "params all 1.0")
plot(log10(sorted_again$`21ATP`),x = sorted_again$timestep)#, xlim = times)
plot(log10(sorted_again$`20num cells`),x = sorted_again$timestep)#, xlim = times)
# plot((sorted_again$`22AA`), (sorted_again$`21ATP`))

thresholded<- subset(important_data, important_data$AminoAcid >= 1 & important_data$ATP >= 1)
plot(log10(thresholded$AminoAcid), log10(thresholded$ATP), col = ceiling(log10(thresholded$k_a_NH4 +1)), main = "coloured by k_a_NH4")
legend(x = 7, y = 7, legend = c('0.1','10','1000'), pch = 1,col = unique(ceiling(log10(thresholded$k_a_NH4 +1))))

selected<- subset(important_data, log10(important_data$ATP)>6 & log10(important_data$AminoAcid)>5)
plot(log10(selected$AminoAcid), log10(selected$ATP))
max(selected$AminoAcid)

# under<- subset(important_data, important_data$AminoAcid < 1000 | important_data$ATP < 1000)
# plot(log(under$AminoAcid), log(under$ATP))
# 
# ordered<- thresholded[order(thresholded$AminoAcid,thresholded$ATP),]
# hist(log(ordered$k_ribo_a))
# 
# good_points<- subset(under, log10(under$ATP)>5.5)
# plot(log(good_points$AminoAcid), log(good_points$ATP))


moles<- 0.0015 #moles of ATP in an e coli
molecules<- moles * 6 *10^23

# selected<- subset(good_points, log(good_points$AminoAcid)> 5.24)
# plot(log(selected$AminoAcid), log(selected$ATP))

getmode(selected$k_ribo_a)
hist(log(selected$k_ribo_a))

getmode<- function(v){
  uniqv<- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
  
}

# hists<- ggplot(under, aes(log(under$k_a_NH4) )) + geom_histogram(alpha = 0.2)
# hists

non_neg_data <- subset(important_data, important_data$AminoAcid >= 0 & important_data$ATP >= 0)
non_neg_data$log_AA<- log(non_neg_data$AminoAcid)
non_neg_data$log_ATP<- log(non_neg_data$ATP)
non_neg_data$product<- as.numeric(as.character(non_neg_data[,3]))* as.numeric(as.character(non_neg_data[,4]))* as.numeric(as.character(non_neg_data[,5]))* as.numeric(as.character(non_neg_data[,6]))* as.numeric(as.character(non_neg_data[,7]))* as.numeric(as.character(non_neg_data[,8]))* as.numeric(as.character(non_neg_data[,9]))
non_neg_data$product_col<- log(non_neg_data$product)
non_neg_data$product_col<- non_neg_data$product_col - min(non_neg_data$product_col)

# plt <- ggplot(data = non_neg_data, aes(x = k_ribo_AA_a, y = , colour = AminoAcid)) +geom_point() +facet_grid(rows = non_neg_data$k_ribo_a)
# plt
# plot(x = important_data$k_ribo_AA_a, y = important_data$k_ribo_a_AA, col = important_data$AminoAcid)
plot(x = non_neg_data$AminoAcid, y = non_neg_data$ATP)
plot(x = non_neg_data$log_AA, y = non_neg_data$log_ATP)#, col = non_neg_data$k_ribo_AA_a, main = "colored by k_a_NH4")
plot(x = non_neg_data$log_AA, y = non_neg_data$log_ATP, col = non_neg_data$k_cat_AA ,main = "colored by k_cat_AA")



df1<- subset(thresholded, thresholded$k_a_NH4== '0.1')
df2<- subset(thresholded, thresholded$k_a_NH4== '10')
df3<- subset(thresholded, thresholded$k_a_NH4== '1000')
plot(x = log10(df1$AminoAcid), y = log10(df1$ATP), main = "k_a_NH4 = 0.1")
plot(x = log10(df2$AminoAcid), y = log10(df2$ATP), main = "k_a_NH4 = 10")
plot(x = log10(df3$AminoAcid), y = log10(df3$ATP), main = "k_a_NH4 = 1000")


hist(x = log10(non_neg_data$k_cat_AA))
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


