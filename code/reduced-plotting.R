# rm(list=ls())
# graphics.off()

library(ggplot2)
library(reshape2)
# library(gridExtra)
# library(cowplot)

#Plotted and saved model k_cat_AA: 0.01 k_a_NH4: 0.01k_a: 0.01k_NH4: 0.01k_a_AA: 0.01k_NH4_AA: 0.01k_ribo_a: 100.0k_ribo_AA: 0.01k_ribo_a_AA: 0.01k_ribo_AA_a: 100.0

# raw_data= read.csv("../data/high-AA-results-k_AA_0.01k_a_NH40.01k_a0.01k_NH40.01k_a_AA0.01k_NH4_AA0.01k_ribo_a0.01k_ribo_AA0.01k_ribo_a_AA0.01k_ribo_AA_a1.0.csv")
# raw_data = read.csv("../data/high-AA-results-k_AA_0.01k_a_NH40.01k_a0.01k_NH40.01k_a_AA0.01k_NH4_AA0.01k_ribo_a100.0k_ribo_AA100.0k_ribo_a_AA100.0k_ribo_AA_a100.0.csv")
# raw_data= read.csv("../data/high-AA-results-k_AA_0.01k_a_NH40.01k_a0.01k_NH40.01k_a_AA0.01k_NH4_AA0.01k_ribo_a0.01k_ribo_AA0.01k_ribo_a_AA0.01k_ribo_AA_a10000.0.csv")
raw_data= read.csv("../data/0.1/param-sweep-0.1-0.1-1000.0-0.1-0.1-0.1-1000.0-1000.0-10.0.csv")
# raw_data= read.csv("testfile.csv")
# raw_data<- read.csv("reduced-model-output.csv")

raw_data = raw_data[,-c(1)]
sorted_data= raw_data[,order(names(raw_data))]
colnames(sorted_data)<- substring(names(sorted_data),6)

sorted_again<-sorted_data[,c("tamp","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25")]
colnames(sorted_again)= c("timestep","1S_external","2ribo mrna comp","3metab enzyme","4housekpng mrna comp","5trans mrna comp","6transporter prot","7metab mrna comp","8trans mrna","9metab mrna","10housekepng prot","11si","12housekpng mrna","13ribo mrna","14free ribo","15NH4 int","16nit mrna","17nit mrna comp","18nitrogenase","19cumulative NH4","20num cells","21ATP","22AA","23AA prot","24AA mRNA","25AA mrna comp")

melted_data<- melt(sorted_again, id.var="timestep")
melted_data['log_molecules']= log(melted_data['value']+0.1)

# plt = ggplot(data = melted_data, aes(x = timestep, y = log(value)))+geom_line(aes(color = variable))
# plt2=ggplot(data = melted_data, aes(x = timestep, y = log_molecules))+geom_line(aes(color = variable))
# plt
# plt2
############################
#change this to plot a different single ODE
plot(type = "l",y=sorted_again$`22AA`,x= sorted_again$timestep)#, xlim = c(9900,10400), ylab = "num of cells", xlab = "timestep")#, main = "NH4 constant at 1000. gamma modded: k_ribo_a = 0.01 k_ribo_AA= 0.01 k_ribo_a_AA = 0.01 k_ribo_AA_a = 1000.\n AA_prot_params: k_AA= 1000.0 k_cat_AA= 1000.0 k_a_NH4 = 1000.0 k_a = 1000.0 k_NH4 = 1000.0 k_a_AA = 0.01 k_NH4_AA= 0.01")
###########################
tail((sorted_again$`22AA`))

plt3 = ggplot(data = melted_data, aes(x = timestep, y = log(value+0.1))) +theme(legend.position = "none")+geom_point(aes(color = variable))+ facet_wrap(~ variable, ncol = 5)
plt3
# min(sorted_again$`2ribo mrna comp`)


# gam_data = read.csv("../data/gamma-file-1000.0-1000.0-1000.0-1000.0-0.01-0.01-0.01-0.01-0.01-1000.0.csv", header = F)
# gam_data["timestep"]<-seq(1,nrow(gam_data))
# colnames(gam_data)<- c("gamma","AA","atp","timestep")
# melted_gam<- melt(gam_data, id.var = "timestep")
# plot(gam_data$atp)
# plot(gam_data$gamma)
# plot(gam_data$AA)
# plt4<- ggplot(data = melted_gam, aes(y=log(value), x = timestep)) +geom_line(aes(color = variable))
# plt4



# par(mfcol=c(4,1))
# params_data = read.csv("../data/high-AA-params-file-0.01-0.01-0.01-0.01-0.01-0.01-0.01-0.01-0.01-1.0.csv")
# #AA_3= raw_data[raw_data$timestamp>20000,]
# plot(y= raw_data$value22, x=raw_data$timestamp, ylab = "num of AA (k_ribo_AA_a:1.0)")
# plot(params_data$new_AA, ylab = "num of AA produced")
# plot(params_data$AA_prot, ylab ="AA proteins")
# plot(params_data$AA_req, ylab = "AA requirement")