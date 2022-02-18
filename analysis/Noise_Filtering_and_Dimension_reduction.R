library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(data.table)
library(matrixStats)

#FUNCTION DEFINITION
chi.square.test <- function(filter1){
  
  deg.f = ncol(filter1)-1
  
  #Calculating upper and lower limit for critical values for the t-test 
  chilower <- qchisq(0.01/2,df=deg.f)
  chiupper <- qchisq(1-0.01/2,df=deg.f)
  
  #adding a variance column at the end of every row
  filter1$variance = rowVars(as.matrix(filter1))
  
  #Calculating median variance
  median_variance = median(filter1$variance)
  
  filter1 <- transform(filter1, t_stat = (( deg.f*variance / median_variance) ) )
  
  #filtering values not lying in either tails 
  filter2 <- filter(filter1, t_stat > chiupper | t_stat < chilower)
  filter2 <- subset(filter2, select = -c(variance,t_stat) )
  
  return(filter2)
}


filename = "/projectnb/bf528/users/im_not_dead_yet/project_1/edata.csv"

#reading in the file
data.table <- read.csv(filename)
# converting first column to rownames
data <- data.frame(data.table[,-1], row.names=data.table[,1])

#Filter 1
filter_20<-data[rowSums(data > log2(15)) >= (0.2*ncol(data)), ]

#Filter 2
chi_square_result <- chi.square.test(filter_20)

#Filter 3
coeff_variation<-subset(chi_square_result, apply(chi_square_result, 1, function(x) sd(x)/mean(x)) > 0.186)

#DELIVERABLES 

#1
write.csv(coeff_variation, "Noise_filtered_data.csv")

#2
print(paste0('The number of genes that pass all of these thresholds: ', nrow(coeff_variation)))
