#### R-SCRIPT: Multi-class VAR estimation of Multi-Store sales application ####
rm(list=ls())
library("tictoc")
#### Set working directory to Directory "Code" ####
# In R: File --> Change dir... 
# In RStudio: Session --> Set Working Directory --> Choose Directory
# Alternatively, specify the path in the setwd function and execute the following line
# setwd("...../Code") # Specify the path on .....

#### Check if all necessary packages are installed ####
checkpackage<-function(U){
  if((U %in% rownames(installed.packages()))==F){
    install.packages(U)
  }
}
packagelist<-list("igraph","corrplot","ggplot2","scales","Matrix","MASS","magic","rARPACK","JGL","plyr","MTS")
lapply(packagelist,checkpackage)
# Load packages
suppressMessages(suppressWarnings(packages <- lapply(packagelist, FUN = function(x) {
  library(x, character.only = TRUE)
})))

#### SOURCE FUNCTIONS ####
source('Functions/MultiClass_VAR.R')
source('Functions/BootSE.R')
source('Functions/FunctionsVisualTools.R')


#### IMPORT DATA ####
new=load("Data/Data.RData") # Array of dimension NxJxK=76x15x15
y = get(new)
namedir = "./data_R_formulationS/"

density = c(1,5)
TTT = 100
nnn = 20
KKK=5


dd=1
itr = 1

tmp = read.csv(file = paste(namedir,'K',KKK,'_data_',density[dd],'percent_',itr,'.csv',sep=""),header=FALSE)
Data = as.matrix(tmp)
dim(Data)<-c(TTT,nnn,KKK)
colnames(Data) <- NULL


#### Multi-class VAR FIT ####
tic("total fitting time")
FIT<-MultiClass_VAR(Data=Data, P=1,lambda1_min=0.1, lambda1_max=50, lambda1_steps=5,
                                   lambda2_min=0.1, lambda2_max=2, lambda2_steps=5)

toc()
write.csv(FIT$Beta_array,paste(namedir,"result_K",KKK,"_",density[dd],"percent_",itr,".csv",sep=""))
