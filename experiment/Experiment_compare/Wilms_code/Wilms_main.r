#### R-SCRIPT: Multi-class VAR estimation of Multi-Store sales application ####
rm(list=ls())

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
load("Data/Data.RData") # Array of dimension NxJxK=76x15x15

#### Multi-class VAR FIT ####
FIT<-MultiClass_VAR(Data=Data, P=1)
