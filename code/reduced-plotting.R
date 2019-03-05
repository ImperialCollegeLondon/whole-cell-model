rm(list=ls())
graphics.off()

library(ggplot2)
library(reshape2)
library(gridExtra)
library(cowplot)

#Plotted and saved model k_cat_AA: 0.01 k_a_NH4: 0.01k_a: 0.01k_NH4: 0.01k_a_AA: 0.01k_NH4_AA: 0.01k_ribo_a: 100.0k_ribo_AA: 0.01k_ribo_a_AA: 0.01k_ribo_AA_a: 100.0

# raw_data= read.csv("../data/high-AA-results-k_AA_0.01k_a_NH40.01k_a0.01k_NH40.01k_a_AA0.01k_NH4_AA0.01k_ribo_a0.01k_ribo_AA0.01k_ribo_a_AA0.01k_ribo_AA_a1.0.csv")
# raw_data = read.csv("../data/high-AA-results-k_AA_0.01k_a_NH40.01k_a0.01k_NH40.01k_a_AA0.01k_NH4_AA0.01k_ribo_a100.0k_ribo_AA100.0k_ribo_a_AA100.0k_ribo_AA_a100.0.csv")
# raw_data= read.csv("../data/high-AA-results-k_AA_0.01k_a_NH40.01k_a0.01k_NH40.01k_a_AA0.01k_NH4_AA0.01k_ribo_a0.01k_ribo_AA0.01k_ribo_a_AA0.01k_ribo_AA_a10000.0.csv")
raw_data= read.csv("reduced-model-output.csv")


raw_data = raw_data[,-c(1)]
sorted_data= raw_data[,order(names(raw_data))]
colnames(sorted_data)<- substring(names(sorted_data),6)

sorted_again<-sorted_data[,c("tamp","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25")]
colnames(sorted_again)= c("timestep","S_external 1","ribo mrna comp 2","metab enzyme 3","housekpng mrna comp 4","trans mrna comp 5","transporter prot 6","metab mrna comp 7","trans mrna 8","metab mrna 9","housekepng prot 10","si11","housekpng mrna 12","ribo mrna 13","free ribo 14","NH4 int 15","nit mrna 16","nit mrna comp 17","nitrogenase18","cumulative NH4 19","num cells 20","ATP 21","AA 22","AA prot 23","AA mRNA 24","AA mrna comp 25")

melted_data<- melt(sorted_again, id.var="timestep")
melted_data['log_molecules']= log(melted_data['value']+0.1)

plt = ggplot(data = melted_data, aes(x = timestep, y = log(value)))+geom_line(aes(color = variable))
plt2=ggplot(data = melted_data, aes(x = timestep, y = log_molecules))+geom_line(aes(color = variable))
plt
plt2




par(mfcol=c(4,1))
params_data = read.csv("../data/high-AA-params-file-0.01-0.01-0.01-0.01-0.01-0.01-0.01-0.01-0.01-1.0.csv")
AA_3= raw_data[raw_data$timestamp>20000,]
plot(y= AA_3$value22, x=AA_3$timestamp, ylab = "num of AA (k_ribo_AA_a:1.0)")
plot(params_data$new_AA, ylab = "num of AA produced")
plot(params_data$AA_prot, ylab ="AA proteins")
plot(params_data$AA_req, ylab = "AA requirement")