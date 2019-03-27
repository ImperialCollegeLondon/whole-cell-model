# rm(list=ls())
# graphics.off()

library(ggplot2)
library(reshape2)
library(gridExtra)
# library(cowplot)
library(R.utils)



good_combos<- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
bad_combos<-as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
oss_combos<-as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))

#initalise the two data frames to store the good/bad param combinations as well as those that produce oscillations


e=NULL
for (i1 in c('0.1','10.0','1000.0')){
  for (i2 in c('0.1','10.0','1000.0')){
    for (i3 in c('0.1','10.0','1000.0')){
      for (i4 in c('0.1','10.0','1000.0')){
        for (i5 in c('0.1','10.0','1000.0')){
          for (i6 in c('0.1','10.0','1000.0')){
            
            
            closeAllConnections()
            filename<- paste("../data/",i1,"/param-sweep-",i1,"-",i2,"-",i3,"-",i4,"-",i5,"-",i6,".csv", sep = "")
            # print(paste("trying file", filename))
            e<-tryCatch({expr=assign("line_num",countLines(filename))},
                        error = function(e){},
                        warning= function(w){})
        
            
            if (typeof(e[1]) != "integer" ){
              # print("file doesnt exist")
              current_vals<- data.frame(i1,i2,i3,i4,i5,i6)
              bad_combos <- rbind(bad_combos,current_vals, stringsAsFactors= FALSE)
              print("real bad one!")
              next
            }else{
              if (line_num> 25000){#if there are more rows than timepoints that means this parameter combination causes oscillations
                current_vals<- data.frame(i1,i2,i3,i4,i5,i6)
                oss_combos <- rbind(oss_combos,current_vals, stringsAsFactors= FALSE)
                print("adding to oss")
                next
              }
              
              
              # raw_data<- read.csv(filename,header = F, skip = (line_num-1))
              # reads in only the last line of data to use less memory
              # non_zero<- raw_data[raw_data != 0]
              
              
              full_data<- read.csv(filename) #read in all rows of data to implement check for minimum boundary on all species
              full_data = full_data[-c(1),-c(1)]#remove the row number and the initial values since they have lots of zeros

              sorted_data= full_data[,order(names(full_data))]
              colnames(sorted_data)<- substring(names(sorted_data),6)
              sorted_data<-sorted_data[,c("tamp","1","2","3","4","5","6","7","8","9","10","11","12","13","21","22","23","24","25")]
              colnames(sorted_data)= c("timestep","1S_external","2ribo mrna comp","3metab enzyme","4housekpng mrna comp","5trans mrna comp","6transporter prot","7metab mrna comp","8trans mrna","9metab mrna","10housekepng prot","11si","12housekpng mrna","13ribo mrna","21ATP","22AA","23AA prot","24AA mRNA","25AA mrna comp")
              #sorts and names all the columns properly
              
              end_min<- min(tail(sorted_data))
              print(paste("current min val is", end_min))
              if (end_min < 5){
                #add parameter combination to exclusion table
                current_vals<- data.frame(i1,i2,i3,i4,i5,i6)
                bad_combos <- rbind(bad_combos,current_vals, stringsAsFactors= FALSE)
                print("adding to bad")
    
              }else{
                #add parameter combination to good table
                current_vals<- data.frame(i1,i2,i3,i4,i5,i6)
                good_combos <- rbind(good_combos,current_vals, stringsAsFactors= FALSE)
                print("adding to good")
              }#need to check for oscillations somehow
              #####################################################################################
              #####################################################################################
              
              # melted_data<- melt(sorted_data, id.var="timestep")
              # melted_data['log_molecules']= log10(melted_data['value']+0.1) #puts data into long format for ggplot
              # plt3 = ggplot(data = melted_data, aes(x = timestep, y = log10(value+0.1))) +theme(legend.position = "none")+geom_point(aes(color = variable))+ facet_wrap(~ variable, ncol = 5)
              # plt3
              
              # min_val<- min(sorted_again)
              
              # current_vals<- data.frame(raw_data$V14,raw_data$V23,i1,i2,i3,i4,i5,i6)
              # important_data <- rbind(important_data,current_vals, stringsAsFactors= FALSE)
              next
            }
            
            
          }}}}}}

print("all done!")

colnames(good_combos)<- c("k_cat_AA","k_a_NH4","k_NH4","k_a_AA","k_NH4_AA","k_ribo_a")
colnames(bad_combos)<- c("k_cat_AA","k_a_NH4","k_NH4","k_a_AA","k_NH4_AA","k_ribo_a")
colnames(oss_combos)<- c("k_cat_AA","k_a_NH4","k_NH4","k_a_AA","k_NH4_AA","k_ribo_a")
#need to append to file not overwrite it!!!#
# write.csv(important_data, file = "param-sweep6-results.csv")
write.csv(good_combos, file = "the-good.csv")
write.csv(bad_combos, file = "the-bad.csv")
write.csv(oss_combos, file = "the-ugly.csv")




