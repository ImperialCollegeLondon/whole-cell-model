rm(list=ls())
graphics.off()

library(ggplot2)
library(reshape2)
library(gridExtra)
# library(cowplot)
library(R.utils)


# name_fun<-function(file_name){assign("line_num",countLines(file_name))}


important_data <- data.frame(matrix(data = NA, nrow=0, ncol = 12))
e=NULL
# colnames(important_data)<- c("AminoAcid","ATP","k_cat_AA","k_a_NH4","k_a","k_NH4","k_a_AA","k_NH4_AA","k_ribo_a","k_ribo_AA","k_ribo_a_AA","k_ribo_AA_a")


#######################################################################
# "../data/param-sweep-0.1-0.1-0.1-1.0-10.0-100.0-0.1-10.0-1.0-1.0.csv"
#this is the latest file reached. start processing from here onward to
#save time
#######################################################################
  
for (i1 in c('0.1','10.0','1000.0')){
  for (i2 in c('0.1','10.0')){#,'1000.0')){
    for (i3 in c('0.1','10.0','1000.0')){
      for (i4 in c('1000.0')){ #k_aa fixed as 1000
        for (i5 in c('0.1','10.0','1000.0')){
          for (i6 in c('0.1','10.0','1000.0')){
            for (i7 in c('0.1','10.0','1000.0')){
              for (i8 in c('1.0')){#,'10.0','1000.0')){
                for (i9 in c('1.0')){#,'10.0','1000.0')){

closeAllConnections()
filename<- paste("../data/",i1,"/param-sweep-",i1,"-",i2,"-",i3,"-",i4,"-",i5,"-",i6,"-",i7,"-",i8,"-",i9,".csv", sep = "")

e<-tryCatch({expr=assign("line_num",countLines(filename))},
            error = function(e){},
            warning= function(w){})



# e<-tryCatch({expr=assign("line_num",countLines(filename))},
#             error = function(e) e,
#             warning= function(w) w)
            # ,finally = next)

if (typeof(e[1]) != "integer" ){
  # print("file doesnt exist")
  next
}else{
  
# print("reading in data")
raw_data<- read.csv(filename,header = F, skip = (line_num-1))
# print("data read in")

current_vals<- data.frame(raw_data$V14,raw_data$V23,i1,i2,i3,i4,i5,i6,i7,i8,i9)
important_data <- rbind(important_data,current_vals, stringsAsFactors= FALSE)
next
}

# raw_data<- read.csv("../data/param-sweep-0.1-0.1-0.1-0.1-0.1-0.1-0.1-0.1-10.0-100.0.csv")#, skip = (line_num-5))
# raw_data = raw_data[,-c(1)]
# 
# 
# unsorted_colnames<- names(raw_data)
# shortened_colnames<-substring(unsorted_colnames,6) 

# sorted_data= raw_data[,order(names(raw_data))]

# raw_data<- read.csv("../data/param-sweep-0.1-0.1-0.1-0.1-0.1-0.1-0.1-0.1-10.0-100.0.csv",header = F, skip = (line_num-1))
# raw_data = raw_data[,-c(1)]
# colnames(raw_data)<- shortened_colnames
# sorted_data<-raw_data[,c("tamp","1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25")]
# colnames(sorted_data)= c("timestep","1S_external","2ribo mrna comp","3metab enzyme","4housekpng mrna comp","5trans mrna comp","6transporter prot","7metab mrna comp","8trans mrna","9metab mrna","10housekepng prot","11si","12housekpng mrna","13ribo mrna","14free ribo","15NH4 int","16nit mrna","17nit mrna comp","18nitrogenase","19cumulative NH4","20num cells","21ATP","22AA","23AA prot","24AA mRNA","25AA mrna comp")

# melted_data<- melt(sorted_again, id.var="timestep")
# melted_data['log_molecules']= log(melted_data['value']+0.1)

# param_values<- 

# important_data <- data.frame(raw_data$V14,raw_data$V23)
# colnames(important_data)<- c("AminoAcid","ATP")

                  }}}}}}}}}
colnames(important_data)<- c("AminoAcid","ATP","k_cat_AA","k_a_NH4","k_a","k_NH4","k_a_AA","k_NH4_AA","k_ribo_a","k_ribo_a_AA","k_ribo_AA_a")

print("all done!")

#need to append to file not overwrite it!!!#
write.csv(important_data, file = "param-sweep2-results.csv")

