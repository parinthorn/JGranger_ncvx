source("0_library_declare.r")
# put data in format same as an input for mat.setup
# variable D, p, K
run<-1
Thresh2 <- 0.02 
criter <- "BIC"                                   # selection criterion for group lasso (first stage), might be "BIC", "AIC", "AICc"
print(c("criter",criter))
criter.second <- "BIC"                            # selection criterion for sparse lasso (second stage), might be "BIC", "AIC", "AICc"
print(c("criter.second:",criter.second))
max.iter <- 30                                    # maximum number of iterations for the two-stage estimation algorithm
print(c("max.iter:",max.iter))
eps <- 0.0001                                     # stopping criterion value for two-stage estimation algorithm
print(c("eps:",eps))
intercept <- FALSE
standardize <- TRUE
D <-1
p <- 20
t <- 100
K <- 50
density = c(1,5)
realz <- 20
Final.Est <- make.list(make.list(matrix(0,D*p,p),K),1)
Final.Comm.Est <- make.list(make.list(matrix(0,D*p,p),K),1)
Final.Ind.Est <-  make.list(make.list(matrix(0,D*p,p),K),1)
for (itr in c(40,49,57,65,73,81,89,97)){
  for (dd in 1:1){
    namedir <- './data_R_formulationD/'
    tmp <- read.csv(file = paste(namedir,'K',K,'_data_',density[dd],'percent_',itr,'.csv',sep=""),header=FALSE)
    DATA <-as.matrix(tmp)
    colnames(DATA) <- NULL
    
    
    
    
    sigma2 <- rep(0,K)
    sds <- apply(DATA,1,function(x) sd(x))
    for (k in 1:K) sigma2[k] <- OLS.tseries(DATA[(k-1)*p + 1:p,],D=D)$sigma2
    
    
    
    Group.Est <- make.list(matrix(0,D*p,p),K)
    Sep.Est.Second <-  make.list(matrix(0,D*p,p),K)
    Group.Final <- make.list(matrix(0,D*p,p),K)
    
    #A.list <- list()
    #
    #for (d in 1:D){
    #  A.full <- list()
    #  for (k in 1:K) A.full[[k]] <- A.true$A.true[[k]][[d]]
    #  A.list[[d]] <- block.diag(A.full)
    #}
    
    Sigma <- list()
    for(i in 1:K) Sigma[[i]] <- diag(1,p)               # simply setting the error covariance to be identity matrix
    
    Sigma.full <- block.diag(Sigma)
    
    #############################
    ### FULL PROBLEM SETUP   ####
    #############################
    M.setup <- mat.setup(DATA,t,K,p,D=D)
    C.list <- M.setup$C
    X.list <- M.setup$X
    
    ### Vector of group number assignments for group lasso
    group <- c(1:(D*p))
    if (K>1){
      for (j in 2:K) group <- c(group,1:(D*p))
    }
    
    ## Initializing vectors to contain estimates during algorithm iterations
    grouplasso.result_before <- make.list(numeric(D*K*p),p)
    seplasso.result_before <- make.list(numeric(D*K*p),p)
    seplasso.result <- make.list(numeric(D*K*p),p)
    grouplasso.result <- make.list(numeric(D*K*p),p)
    gl.zeros <- list()
    
    ## Going column-by-column, j - column index.
    for (j in 1:p){
      cat("\n")
      print(paste("j:",j,sep=""))
      # cat("\n")
      
      
      it <- 0
      flag <- 0
      
      ### First 5 iterations to tune up the regularization parameters.
      ### (or less than 5, if regularization parameter values stop changing => we can just fix them)
      while ((it<5) & (flag==0)){
        it <- it+1
        #print(paste("Tune Iter:",it,sep=""))
        
        ###############################################################
        ## FIRST STAGE: Group lasso estimation of common component ####
        ###############################################################
        
        ### Initializing response vector and data matrix for standard regression problems for the first stage
        Y <- numeric(K*(t-D))
        Xbeta <- X.list %*% seplasso.result[[j]]
        for (k in 1:K){
          Y[(k-1)*(t-D) + (1:(t-D))] <- C.list[1:(t-D),j + (k-1)*p] - Xbeta[(k-1)*(t-D) + (1:(t-D))]
        }
        
        ### Doing group lasso optimization
        D.sigma <- sqrt(diag(c(sapply(sigma2, function(x) return(rep(x,(t-D)))))))
        
        r <- grpreg(solve(D.sigma) %*% X.list,
                    solve(D.sigma) %*% Y,
                    group=group,
                    penalty="grLasso",
                    family="gaussian",
                    intercept=intercept,
                    warn=FALSE)
        # lambda=lambda.group)
        
        lambda_G.path <- r$lambda
        est <- r$beta[-1,]
        grouplasso.result.df <- r$df
        
        ### Tuning parameter selection
        if (criter == "AIC")
          group.est.out <- AIC(as.matrix(est),X.list,as.matrix(Y),p=p,K=K,lambda.path=lambda_G.path,df.path=grouplasso.result.df)
        if (criter == "BIC")
          group.est.out <- BIC(as.matrix(est),X.list,as.matrix(Y),p=p,K=K,df.path=grouplasso.result.df,lambda.path=lambda_G.path)
        if (criter == "AICc")
          group.est.out <- AICc(as.matrix(est),X.list,as.matrix(Y),p=p,K=K,df.path=grouplasso.result.df,lambda.path=lambda_G.path)
        
        # grouplasso.result[[j]] <- sparsify(group.est.out$Est, Thresh2)
        grouplasso.result[[j]] <- group.est.out$Est
        lambda.group <- group.est.out$lambda1
        
        
        ####################################################################
        ## SECOND STAGE: Sparse lasso estimation of individual component ###
        ####################################################################
        
        ### Initializing response vector and data matrix for standard regression problem for the second stage
        gl.zeros[[j]] <- (grouplasso.result[[j]] == 0)
        Xbeta <- X.list %*% grouplasso.result[[j]]
        for (k in 1:K){
          Y[(k-1)*(t-D) + (1:(t-D))] <- C.list[1:(t-D),j + (k-1)*p] - Xbeta[(k-1)*(t-D) + (1:(t-D))]
        }
        X.zero <- X.list[,gl.zeros[[j]]]
        
        ### Doing sparse lasso optimization
        r <- glmnet(solve(D.sigma) %*% X.zero,
                    solve(D.sigma) %*% Y,
                    family="gaussian",
                    standardize=standardize,
                    intercept=intercept)
        
        lambda_SPARS.path <- r$lambda
        est <- r$beta
        sep.df <- r$df
        
        ### Tuning parameter selection
        if (criter.second == "AIC")
          sep.est.out <- AIC(as.matrix(est),X.zero,as.matrix(Y),p=p,K=K,df.path=sep.df,lambda.path=lambda_SPARS.path)
        if (criter.second == "BIC")
          sep.est.out <- BIC(as.matrix(est),X.zero,as.matrix(Y),p=p,K=K,df.path=sep.df,lambda.path=lambda_SPARS.path)
        if (criter.second == "AICc")
          sep.est.out <- AICc(as.matrix(est),X.zero,as.matrix(Y),p=p,K=K,df.path=sep.df,lambda.path=lambda_SPARS.path)
        
        lambda.sparse <- sep.est.out$lambda1
        seplasso.result[[j]] <- rep(0,D*p*K)
        # seplasso.result[[j]][gl.zeros[[j]]] <- sparsify(sep.est.out$Est, Thresh2)
        seplasso.result[[j]][gl.zeros[[j]]] <- sep.est.out$Est
        
        ### Recording estimates
        for (k in 1:K){
          Group.Est[[k]][j,] <- grouplasso.result[[j]][(k-1)*p+(D-1)*(p) + (1:p)]
          Sep.Est.Second[[k]][j,] <- seplasso.result[[j]][(k-1)*p+(D-1)*(p) + (1:p)]
        }
        
        if (it>1){
          #####
          ## STOPPING CRITERION for the tuning parameter initialization steps (if it hasn't reached 5 iterations yet)
          #####
          diff_magnitudes <- sum((c(grouplasso.result[[j]],seplasso.result[[j]])-c(grouplasso.result_before[[j]],seplasso.result_before[[j]]))^2)
          
          ## Printing out current difference in optimized values for consecutive iterations.
          print(paste("Tune-up:", diff_magnitudes))
          
          if (diff_magnitudes < eps){
            flag <- 1
            #n.iter[run,j] <- it
          }
        }
        ### keeping track of estimates from previous iterations
        grouplasso.result_before <- grouplasso.result
        seplasso.result_before <- seplasso.result
      }
      
      # cat("\n")
      
      # print("lambdas selected")
      # print(paste("Group lambda:", lambda.group))
      # print(paste("Sparse lambda:", lambda.sparse))
      
      # cat("\n")
      
      
      #####
      ## Now, FIX the values of lambda.group and lambda.sparse.
      ## Iterate until convergence for fixed lambda.group and lambda.sparse
      #####
      
      it <- 0
      flag <- 0
      
      while ((it<max.iter) & (flag==0)){
        it <- it+1
        #  print(paste("Iter:",it,sep=""))
        
        ###############################################################
        ## FIRST STAGE: Group lasso estimation of common component ####
        ###############################################################
        
        ### Initializing response vector and data matrix for standard regression problems for the first stage
        Y <- numeric(K*(t-D))
        Xbeta <- X.list %*% seplasso.result[[j]]
        for (k in 1:K){
          Y[(k-1)*(t-D) + (1:(t-D))] <- C.list[1:(t-D),j + (k-1)*p] - Xbeta[(k-1)*(t-D) + (1:(t-D))]
        }
        
        ### Doing group lasso optimization
        D.sigma <- sqrt(diag(c(sapply(sigma2, function(x) return(rep(x,(t-D)))))))
        
        r <- grpreg(solve(D.sigma) %*% X.list,
                    solve(D.sigma) %*% Y,
                    group=group,
                    penalty="grLasso",
                    family="gaussian",
                    intercept=intercept,
                    warn=FALSE,
                    lambda=lambda.group)
        
        #lambda_G.path <- r$lambda
        est <- r$beta[-1]
        grouplasso.result.df <- r$df
        
        #grouplasso.result[[j]] <- sparsify(est, Thresh2)
        grouplasso.result[[j]] <- est
        # lambda.group <- group.est.out$lambda1
        
        
        ####################################################################
        ## SECOND STAGE: Sparse lasso estimation of individual component ###
        ####################################################################
        
        ### Initializing response vector and data matrix for standard regression problem for the second stage
        gl.zeros[[j]] <- (grouplasso.result[[j]] == 0)
        Xbeta <- X.list %*% grouplasso.result[[j]]
        for (k in 1:K){
          Y[(k-1)*(t-D) + (1:(t-D))] <- C.list[1:(t-D),j + (k-1)*p] - Xbeta[(k-1)*(t-D) + (1:(t-D))]
        }
        X.zero <- X.list[,gl.zeros[[j]]]
        
        ### Doing sparse lasso optimization
        r <- glmnet(solve(D.sigma) %*% X.zero,
                    solve(D.sigma) %*% Y,
                    family="gaussian",
                    standardize=standardize,
                    intercept=intercept,
                    lambda=lambda.sparse)
        
        #lambda_SPARS.path <- r$lambda
        est <- as.matrix(r$beta)
        sep.df <- r$df
        
        #lambda.sparse <- sep.est.out$lambda1
        seplasso.result[[j]] <- rep(0,D*p*K)
        # seplasso.result[[j]][gl.zeros[[j]]] <- sparsify(sep.est.out$Est, Thresh2)
        seplasso.result[[j]][gl.zeros[[j]]] <- est
        
        ### Recording estimates
        for (k in 1:K){
          Group.Est[[k]][j,] <- sparsify(grouplasso.result[[j]][(k-1)*p+(D-1)*(p) + (1:p)], Thresh2)
          Sep.Est.Second[[k]][j,] <- sparsify(seplasso.result[[j]][(k-1)*p+(D-1)*(p) + (1:p)], Thresh2)
        }
        
        if (it>1){
          #####
          ## STOPPING CRITERION for the two-stage estimation algorithm
          #####
          diff_magnitudes <- sum((c(grouplasso.result[[j]],seplasso.result[[j]])-c(grouplasso.result_before[[j]],seplasso.result_before[[j]]))^2)
          
          ## Printing out current difference in optimized values for consecutive iterations.
          print(paste("Convergence:", diff_magnitudes))
          
          if (diff_magnitudes < eps){
            flag <- 1
            #n.iter[run,j] <- it
          }
        }
        ### keeping track of estimates from previous iterations
        grouplasso.result_before <- grouplasso.result
        seplasso.result_before <- seplasso.result
      }
      
      #if (it == max.iter) n.iter[run,j] <- max.iter
      # print(it)
    }
    
    
    Group.Est.Common <- list()
    Group.Est.Common <- Group.Est
    
    for(k in 1:K){
      Group.Final[[k]] <- Group.Est[[k]] + Sep.Est.Second[[k]]
      Final.Comm.Est[[run]][[k]] <- Group.Est[[k]]
      Final.Ind.Est[[run]][[k]] <- Sep.Est.Second[[k]]
      Group.Final[[k]] <- sparsify(Group.Final[[k]],Thresh2)
      Final.Est[[run]][[k]] <- Group.Final[[k]]
    }
    #saveRDS(Final.Est,paste(namedir,"/",itr,"_","Final_Est.rds",sep=""))
    #saveRDS(Final.Comm.Est,paste(namedir,"/",itr,"_","Common_Est.rds",sep=""))
    #saveRDS(Final.Ind.Est,paste(namedir,"/",itr,"_","Ind_Est.rds",sep=""))
    write.csv(Final.Est,paste(namedir,"K",K,"_Final_Est_",density[dd],"percent_",itr,".csv",sep=""))
    write.csv(Final.Comm.Est,paste(namedir,"K",K,"Common_Est_",density[dd],"percent_",itr,".csv",sep=""))
    write.csv(Final.Ind.Est,paste(namedir,"K",K,"Ind_Est_",density[dd],"percent_",itr,".csv",sep=""))
  }
}
