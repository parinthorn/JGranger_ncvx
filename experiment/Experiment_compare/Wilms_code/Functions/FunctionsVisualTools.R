##########################################################################
#### FunctionsVisualTools.R: Functions to make the three visual tools ####
##########################################################################

CLUSTER.PLOT<-function(i, j, bootFIT, main= " ", G, K, J, type=c("PRICE", "PROMO","SLS")){
  #### Function to obtain Cluster Plot ####
  # INPUT
  # i,j: effect of time series i on time series j
  # bootFIT: Bootstrap FIT
  # G: number of product categories
  # K: number of classes
  # J: number of time series
  # type: c("PRICE", "PROMO", "SLS"): effects of prices on sales (type="PRICE"); promotion on sales (type="PROMO"), or sales on sales (type="SLS")
  
  # OUTPUT
  # CLUSTER PLOT estimated effect of time series i on time series j
  
  plotSE<-function(y=0.5,x,SE,vlength=0.05,col="black"){
    # AUX Function to plot 95% confidence interval based on point estimate x and standard error SE
    up<-x + 1.96*SE
    low<-x - 1.96*SE
    segments(low,y , up, y, col=col)
    segments(low, y-vlength , low, y+vlength, col=col)
    segments(up, y-vlength , up, y+vlength, col=col)
  }
  
  # Collect estimated effects : Prices on Sales, Promo on Sales and Sales on Sales #
  beta_classes <- matrix(bootFIT$Beta, J^2, K)
  betasd_classes <- matrix(apply(bootFIT$BetaBoot,2,sd), J^2, K)
  
  beta_PRICE <- beta_PROMO <- beta_SLS <-array(NA,c(G,G,K)) 
  for( i.class in 1:K){
    beta_perclass <- matrix(beta_classes[,i.class], ncol=J, nrow=J)
    beta_SLS[,,i.class] <- matrix(beta_perclass[(1:G), (1:G)], G, G)  
    beta_PROMO[,,i.class] <- matrix(beta_perclass[((G+1):(G*2)), (1:G)], G, G) 
    beta_PRICE[,,i.class] <- matrix(beta_perclass[((G*2+1):(G*3)), (1:G)], G, G) 
  }
  
  # Collect standard errors : Prices on Sales, Promo on Sales and Sales on Sales #
  betasd_PRICE <- betasd_PROMO <- betasd_SLS <-array(NA,c(G,G,K)) 
  for( i.class in 1:K){
    betasd_perclass <- matrix(betasd_classes[,i.class], ncol=J, nrow=J)
    betasd_SLS[,,i.class] <- matrix(betasd_perclass[(1:G), (1:G)], G, G)  
    betasd_PROMO[,,i.class] <- matrix(betasd_perclass[((G+1):(G*2)), (1:G)], G, G) 
    betasd_PRICE[,,i.class] <- matrix(betasd_perclass[((G*2+1):(G*3)), (1:G)], G, G) 
  }
  
  
  catnames <- c("BER", "BJC", "RFJ", "FRJ", "SDR")
  orderstores <- c(3, 5, 12, 1, 13, 6, 14, 7, 2, 4, 15, 8, 9, 10, 11) # Sort stores along price tier type 
  listnames <- list(catnames, catnames, paste(orderstores))
  dimnames(beta_PRICE)<-dimnames(beta_PROMO)<-dimnames(beta_SLS)<-listnames
  
  
  if(type=="PRICE"){
    beta <- beta_PRICE
    betasd <- betasd_PRICE
  }
  
  if(type=="PROMO"){
    beta <- beta_PROMO
    betasd <- betasd_PROMO
  }
  
  if(type=="SLS"){
    beta <- beta_SLS
    betasd <- betasd_SLS
  }  
  
  
  #### CLUSTER PLOT ####
  
  unique<-unique(sort(round(beta[i,j,],2))) # Extract unique estimated coefficients
  cluster.position<-cluster.label<-c()
  for(i.cluster in 1:length(unique)){
    cluster.elements<-as.numeric(names(beta[i,j,])[is.element(round(beta[i,j,],2),unique[i.cluster])])
    cluster.count<-length(cluster.elements)
    cluster.position<-c(cluster.position,1:cluster.count)
    cluster.label<-c(cluster.label,sort(cluster.elements))
  }
  par(mfrow=c(1,1),mar=c(5,2,3,0.5))
  cluster.position<-1:15
  # Plot point estimates
  plot(main=main, sort(round(beta[i,j,],2)),cluster.position,bty="n",yaxt="n",ylab="",xlab="Estimated effect",pch=16,xlim=c(min(unique)-0.2,max(unique)+0.2),ylim=c(0,max(cluster.position)+1),cex.main=1.7,cex.lab=1.4)

  
  # Plot confidence interval
  xpos <- sort(round(beta[i,j,],2))
  sd <- betasd[i,j,]
  dimnames(betasd)[[3]] <- dimnames(beta)[[3]]
  names(sd) <- dimnames(beta)[[3]]
  
  for(i.sd in 1:length(cluster.position)){
    sd_idx<-which(names(sd)==names(xpos)[i.sd])
    sd_i<- sd[sd_idx]
    
    # Significant negative effect 
    if ((xpos[i.sd] + 1.96*sd_i)<0){
      plotSE(x=xpos[i.sd],y=cluster.position[i.sd],SE=sd_i, col="orangered")
      text(min(unique)-0.175,cluster.position[i.sd],labels=cluster.label[i.sd],pos=4,cex=1.2, col="orangered") 
      points(sort(round(beta[i,j,],2))[i.sd],cluster.position[i.sd], pch=15, col="orangered", cex=1.5)
    }
    
    # Non-significant effect 
    if (((xpos[i.sd] + 1.96*sd_i)>0) & ((xpos[i.sd] - 1.96*sd_i)<0)){
      plotSE(x=xpos[i.sd],y=cluster.position[i.sd],SE=sd_i, col="darkgray")
      text(min(unique)-0.175,cluster.position[i.sd],labels=cluster.label[i.sd],pos=4,cex=1.2, col="darkgray")
      points(sort(round(beta[i,j,],2))[i.sd],cluster.position[i.sd], pch=16, col="darkgray", cex=1.5)
    }
    
    # Significant positive effect
    if (((xpos[i.sd] - 1.96*sd_i)>0)){
      plotSE(x=xpos[i.sd],y=cluster.position[i.sd],SE=sd_i, col="blue")
      text(min(unique)-0.175,cluster.position[i.sd],labels=cluster.label[i.sd],pos=4,cex=1.2, col="blue") 
      points(sort(round(beta[i,j,],2))[i.sd],cluster.position[i.sd], pch=17, col="blue", cex=1.5)
    }
  }
  
}

NETWORK.PLOT<-function(bootFIT, G, K, J, type, Kclasses, plot.row=3, plot.column=5, edge.tick=0.18, main.label="Store"){
  #### Function to obtain Network plot ####
  
  # INPUT
  # bootFIT: Bootstrap FIT
  # G: number of product categories
  # K: number of classes
  # J: number of time series
  # type: c("PRICE", "PROMO", "SLS"): effects of prices on sales (type="PRICE"); promotion on sales (type="PROMO"), or sales on sales (type="SLS")
  # Kclasses: labels of class numbers 
  # plot.row: number of network rows in plot
  # plot.column: number of network columns in plot
  
  # OUTPUT
  # NETWORK PLOT OF PRODUCT CATEGORIES: K networks, each consisting of J nodes
  
  # Collect significant parameter estimates
  beta_classes <- matrix(bootFIT$BetaSign, J^2, K)
  beta_PRICE <- beta_PROMO <- beta_SLS <-array(NA,c(G,G,K)) 
  for( i.class in 1:K){
    beta_perclass <- matrix(beta_classes[,i.class], ncol=J, nrow=J)
    beta_SLS[,,i.class] <- matrix(beta_perclass[(1:G), (1:G)], G, G)  
    beta_PROMO[,,i.class] <- matrix(beta_perclass[((G+1):(G*2)), (1:G)], G, G) 
    beta_PRICE[,,i.class] <- matrix(beta_perclass[((G*2+1):(G*3)), (1:G)], G, G) 
  }
  orderstores <- c(3,5,12,1,13,6,14,7,2,4,15,8,9,10,11) # Sort stores along price tier type 
  listnames <- list(catnames, catnames, paste(orderstores))
  dimnames(beta_PRICE) <- dimnames(beta_PROMO) <- dimnames(beta_SLS) <- listnames
  
  if(type=="PRICE"){
    beta <- beta_PRICE
  }
  
  if(type=="PROMO"){
    beta <- beta_PROMO
  }
  
  if(type=="SLS"){
    beta <- beta_SLS
  } 
  
  radian.rescale <- function(x, start=0, direction=1) {
    c.rotate <- function(x) (x + start) %% (2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  lab.locs <- radian.rescale(x=1:J, direction=-1, start=0) 
  
  par(mfrow=c(plot.row, plot.column), mar=c(2, 2, 1.5, 1.5) + 0.1) #mar: c(bottom, left, top, right)
  for( i.class in Kclasses){
    GRAPH<-graph.adjacency(adjmatrix=beta[,,i.class], mode="directed", diag=FALSE,weighted=T)
    if( i.class==Kclasses[1]){ #Fix locations for first store
      FIX_LOCATION<-layout.circle(GRAPH) 
    }
    checkzero<-try(abs(E(GRAPH)$weight),silent=T)
    if(class(checkzero)=="try-error"){
      E(GRAPH)$width <- 1
      E(GRAPH)$color<-ifelse(1, "orangered", "blue")
    }else{
      E(GRAPH)$width <- abs(E(GRAPH)$weight)*edge.tick*10^2 
      E(GRAPH)$color<-ifelse(E(GRAPH)$weight<0, "orangered", "blue")
    }
    
    set.seed(1)
    plot(margin=rep(10^-10,4),GRAPH, layout=FIX_LOCATION, vertex.size=2,vertex.label.dist=0.8, vertex.color="black", edge.arrow.size=0.5,vertex.label.cex=1.5,vertex.shape="circle",vertex.label.color="black")
    title(paste(main.label,which(i.class==Kclasses)),cex.main=1.2)
  }
}

SIMILARITY.PLOT<-function(bootFIT, G, K, J, type, Kclasses, main.label="Class", title=" "){
  #### Function to obtain Similarity Matrices ####
  
  # INPUT
  # bootFIT: Bootstrap FIT
  # G: number of product categories
  # K: number of classes
  # J: number of time series
  # type: c("PRICE", "PROMO", "SLS"): effects of prices on sales (type="PRICE"); promotion on sales (type="PROMO"), or sales on sales (type="SLS")
  # Kclasses: labels of class numbers 
  
  # OUTPUT
  # SIMILARITY MATRIX (KxK) OF SHARED EFFECTS BETWEEN CLASS i and CLASS j
  
  # Collect significant parameter estimates
  beta_classes <- matrix(bootFIT$BetaSign, J^2, K)
  beta_PRICE <- beta_PROMO <- beta_SLS <-array(NA,c(G,G,K)) 
  for( i.class in 1:K){
    beta_perclass <- matrix(beta_classes[,i.class], ncol=J, nrow=J)
    beta_SLS[,,i.class] <- matrix(beta_perclass[(1:G), (1:G)], G, G)  
    beta_PROMO[,,i.class] <- matrix(beta_perclass[((G+1):(G*2)), (1:G)], G, G) 
    beta_PRICE[,,i.class] <- matrix(beta_perclass[((G*2+1):(G*3)), (1:G)], G, G) 
  }
  orderstores <- c(3,5,12,1,13,6,14,7,2,4,15,8,9,10,11) # Sort stores along price tier type 
  listnames <- list(catnames, catnames, paste(orderstores))
  dimnames(beta_PRICE) <- dimnames(beta_PROMO) <- dimnames(beta_SLS) <- listnames
  
  if(type=="PRICE"){
    beta <- beta_PRICE
  }
  
  if(type=="PROMO"){
    beta <- beta_PROMO
  }
  
  if(type=="SLS"){
    beta <- beta_SLS
  } 
  
  SHARED_NONZERO<-matrix(NA,ncol=K,nrow=K)
  colnames(SHARED_NONZERO)<-paste("Class",Kstores)
  rownames(SHARED_NONZERO)<-paste("Class",Kstores)
  
  TPR <- function(A,B){
    # This is a function that compares the structure of two matrixes A and B                    #
    # It outputs the number of entries that A and B have in common that are different from zero #
    # A and B need to have the same number of rows and columns                                  #
    out <- sum(A!=0&B!=0)/sum(A!=0)
  }
  
  for(i in Kclasses){
    for(j in Kclasses){
      SHARED_NONZERO[which(Kclasses==i),which(Kclasses==j)]<-TPR(beta[,,i],beta[,,j])
    }
  }
  SHARED_NONZERO[which(is.na(SHARED_NONZERO))]<-0
  par(mfrow=c(1,1))
  colnames(SHARED_NONZERO)<-rownames(SHARED_NONZERO)<-paste(main.label,1:K)
  corrplot(tl.col = "black",SHARED_NONZERO, method = "circle",cl.lim=c(0,1),order="hclust", hclust.method="ward.D",main=title,mar=c(5, 4, 4, 2) + 0.1)
  SHARED_NONZERO
  }