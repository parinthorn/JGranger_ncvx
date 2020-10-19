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
# To load the object: load("FIT.RData")

#### Residual Bootstrap to obtain standard errors ####
set.seed(20170516)
bootFIT<-bootSE(FITOBJECT=FIT,Nboot=1000,parallell=T)  
# Important Note: Make sure that the sourcing of the functions under #### SOURCE FUNCTIONS #### works 
# as provided in this R-script, since this sourcing is also done within the "bootSE" function.
# To load the object: load("bootFIT.RData")

#### FIGURES ####
##### PRELIMINARIES ####
J <- dim(Data)[2] 
K <- dim(Data)[3] 
G <- 5

#### CLUSTER PLOTS ####
catnames <- c("BER", "BJC", "RFJ", "FRJ", "SDR")
CLUSTER.PLOT(i=which(catnames=="BER"), j=which(catnames=="RFJ"), bootFIT=bootFIT, type="PRICE", G=G, K=K, J=J,
             main="Beer Prices on Refrigerated Juices Sales") # Plot BER Prices on RFJ Sales  
CLUSTER.PLOT(i=which(catnames=="BER"), j=which(catnames=="FRJ"), bootFIT=bootFIT, type="PROMO", G=G, K=K, J=J,
             main="Beer Promotion on Frozen Juices Sales") # Plot BER Promotion on FRJ Sales  
CLUSTER.PLOT(i=which(catnames=="BER"), j=which(catnames=="BJC"), bootFIT=bootFIT, type="SLS", G=G, K=K, J=J,
             main="Beer Sales on Bottled Juices Sales") # Plot BER Sales on BJC Sales 


#### NETWORK PLOTS ####
Kstores <- c(4, 9, 1, 10, 2, 6, 8, 12, 13, 14, 15, 3, 5, 7, 11) # Store Numbers Dominick Data
NETWORK.PLOT(bootFIT=bootFIT, type="PRICE", Kclasses=Kstores, G=G, K=K, J=J) # Prices on Sales


#### SIMILARITY PLOTS ####
SIMILARITY.PLOT(bootFIT=bootFIT, type="PRICE", Kclasses=Kstores, G=G, K=K, J=J, title="(a) Prices on Sales") # Prices on Sales
SIMILARITY.PLOT(bootFIT=bootFIT, type="PROMO", Kclasses=Kstores, G=G, K=K, J=J, title="(b) Promotion on Sales") # Promotion on Sales
SIMILARITY.PLOT(bootFIT=bootFIT, type="SLS", Kclasses=Kstores, G=G, K=K, J=J,title="(c) Sales on Sales") # Sales on Sales


#### Extensions: Use of Main Functions ####

#### Multi-class VAR with adaptive group lasso penalty ####
groupindex<-c()
i.index<-seq(from=1,to=(G*J*K-(G-1)),by=G)
for(i in i.index){
  groupindex<-c(groupindex,rep(seq(from=i,length=G,by=1),J/G)) # each group consists of J/G=3 time series: price, promo and sales
}
FITGroupLasso<-MultiClass_VAR(Data,P=1,type="GrLasso",group=groupindex) 
bootFITGroupLasso<-bootSE(FITOBJECT=FITGroupLasso,Nboot=1000,parallell=T,type="GrLasso",group=groupindex)  

#### Multi-class VAR incorporating additional clustering information ####
load("Data/Distances.RData") # Import additional clustering information
FITDemo<-MultiClass_VAR(Data,P=1,lambda_weights=T, Clusterinfo=Distances) 
bootFITDemo<-bootSE(FITOBJECT=FITDemo,Nboot=1000,parallell=T,lambda_weights=T, Clusterinfo=Distances)

#### Multi-class VAR allowing for cross-error correlation ####
FITOmega<-MultiClass_VAR_generalOmega(Data,P=1) 
bootFITOmega<-bootSE_GeneralOmega(FITOBJECT=FITOmega,Nboot=1000,parallell=T)