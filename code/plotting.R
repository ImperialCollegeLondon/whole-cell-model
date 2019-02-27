rm(list=ls())
graphics.off()

library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)

raw_data = read.csv("../data/AA-results-k_cat_AA-1.0k_a_NH41.0k_a11.0k_NH41.0k_a_AA1.0k_NH4_AA1.0.csv")
raw_data = raw_data[,-c(1)]
# test[ , order(names(test))]
colnames(raw_data)<- as.numeric(substring(names(raw_data),6))
sorted_data= raw_data[,order(names(raw_data))]

colnames(sorted_data)= c("S_external 1","ribo mrna comp 2","metab enzyme 3","housekpng mrna comp 4","trans mrna comp 5","transporter prot 6","metab mrna comp 7","trans mrna 8","metab mrna 9","housekepng prot 10","si11","housekpng mrna 12","ribo mrna 13","free ribo 14","NH4 int 15","nit mrna 16","nit mrna comp 17","nitrogenase18","cumulative NH4 19","num cells 20","ATP 21","AA 22","AA prot 23","AA mRNA 24","AA mrna comp 25", "timestep")

melted_data<- melt(sorted_data, id.var="timestep")
melted_data['log_molecules']= log(melted_data['value']+0.1)

plt = ggplot(data = melted_data, aes(x = timestep, y = log_molecules))+geom_line(aes(color = variable))
plt2=ggplot(data = melted_data, aes(x = timestep, y = log_molecules))+geom_line(aes(color = variable))
# plt+ labs(x = "time", y = "title")

plt
