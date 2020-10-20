######################################################################
# MultiClass_VAR.R:  Functions for sparse multi-class VAR estimation #
######################################################################

MultiClass_VAR<-function(Data=NULL, P,
                         Y=NULL, X=NULL, X0=NULL, XX_prod=NULL, XY_prod=NULL, J=NULL, K=NULL, n=NULL,
                         calculate.C=TRUE, C=NULL, CNorm=NULL, b_init=NULL, mu=1e-4,
                         maxit.both=10, maxit.beta=10, maxit.omega=10, tol.both=10^-5, tol.beta=0.01, tol.omega=0.01,
                         lambda1_min=0.1, lambda1_max=25, lambda1_steps=5,
                         lambda2_min=0.1, lambda2_max=2, lambda2_steps=5,
                         gamma1_min=0.1, gamma1_max=1, gamma1_steps=3,
                         gamma2_min=0.1, gamma2_max=1, gamma2_steps=3,
                         lambda1_OPT=NULL, lambda2_OPT=NULL, gamma1_OPT=NULL, gamma2_OPT=NULL,
                         lambdaRIDGE_min=0.1, lambdaRIDGE_max=10, lambdaRIDGE_steps=5, lambdaRIDGE_OPT=NULL,
                         criterion="BIC", type="AdLasso", gamma_ridge=NULL, group=NULL,
                         lambda_weights=FALSE, Clusterinfo=NULL){

  #### Function for the Multi-class VAR ####

  #########
  # INPUT #
  #########
  # Data: a N*J*K array, where:
          #  - N is the times series length
          #  - J is the number of time series
          #  - K is the number of classes
  # P: order of VAR
  # Y : vector of dimension NJK x 1 of responses (built with 'stack_fun.R'). Default is NULL.
  # X : matrix of dimension NJK x J�PK of predictors (built with 'stack_fun.R'). Default is NULL.
  # X0 : array of dimension n x JP x K of predictors (built with 'stack_fun.R'). Default is NULL.
  # XX_prod : cross-product matrix of dimension J�PK x J�PK. Default is NULL.
  # XY_prod : cross-product matrix of dimension J�PK x 1. Default is NULL.
  # J : number of time series in each class.
  # K : number of classes (K>1).
  # n : time series length.
  # calculate.C : logical for calculating the C matrix.
  # C : C matrix (built with 'Cmatrix.R').
  # CNorm : norm of C (built with 'Cmatrix.R').
  # b_init : vector of dimension J�PK x 1 with the initial values of the autoregressive coefficients.
  # mu : smoothness parameter. Default 1e-4.
  # maxit.both : maximum iterations of the Multiclass VAR
  # maxit.beta : maximum iterations  SPG algorithm
  # maxit.omega : maximum iterations JGL algorithm
  # tol.both : tolerance  Multiclass VAR
  # tol.beta : tolerance SPG algorithm
  # tol.omega : tolerance JGL algorithm
  # lambda1_min : minimum value of the regularization parameter on Beta (lasso penalty, lambda2 in paper!)
  # lambda1_max : maximum value of the regularization parameter on Beta (lasso penalty)
  # lambda1_steps : number of steps in the grid lambda1
  # lambda2_min : minimum value of the regularization parameter on Beta (fusion penalty, lambda1 in paper!)
  # lambda2_max : maximum value of the regularization parameter on Beta (fusion penalty)
  # lambda2_steps : number of steps in the grid lambda2
  # gamma1_min : minimum value of the regularization parameter on Omega (lasso penalty, gamma2 in paper!)
  # gamma1_max : maximum value of the regularization parameter on Omega (lasso penalty)
  # gamma1_steps : number of steps in the grid gamma1
  # gamma2_min : minimum value of the regularization parameter on Omega (fusion penalty, gamma1 in paper!)
  # gamma2_max : maximum value of the regularization parameter on Omega (fusion penalty)
  # gamma2_steps : number of steps in the grid gamma2
  # lambda1_OPT: optimal regularization parameter on B (lasso penalty). Default NULL.
  # lambda2_OPT: optimal regularization parameter on B (fusion penalty). Default NULL.
  # gamma1_OPT: optimal regularization parameter on Omega (lasso penalty). Default NULL.
  # gamma2_OPT: optimal regularization parameter on Omega (fusion penalty). Default NULL.
  # lambdaRIDGE_min : minimum value of the regularization parameter on B for the ridge estimator (only if type="AdLasso").
  # lambdaRIDGE_max : maximum value of the regularization parameter on B for the ridge estimator (only if type="AdLasso").
  # lambdaRIDGE_steps : number of steps in the grid lambda_ridge
  # lambdaRIDGE_OPT : array (1,1,K) containing the optimal regularization parameter on B for the ridge estimator
  # criterion: "BIC" or "AICc", selection criterion used for the selection of the regularization paramaters. Default is BIC.
  # type: "Lasso", "AdLasso", "GrLasso" for Lasso, Adapative Lasso and Group Lasso. Default is Adapative Lasso.
  # gamma_ridge: exponent for the weight for the Adpative Lasso. Default is NULL.
  # group : grouping structure of the Group Lasso. Vector of length KPJ^2. Default is NULL.
  # lambda_weights : logical. If TRUE, then do the weighting of fusion parameter based on Clusterinfo. Default is FALSE.
  # Clusterinfo : K x K matrix of additional clustering information. K is the number of classes. Default is NULL.

  ##########
  # OUTPUT #
  ##########
  # Beta_new : vector of dimension J�PK x 1 with the autoregressive coefficients. Given B^{(k)}_{p,ij}, with i referring to the dependent variable and j to the independent variable, the structure of Beta_new=[B^{(1)}_{1,11}, B^{(1)}_{1,12}, ..., B^{(1)}_{1,1J}, B^{(1)}_{2,11}, B^{(1)}_{2,12}, ..., B^{(1)}_{P,1J}, B^{(1)}_{1,21}, ...,B^{(1)}_{P,JJ}, B^{(2)}_{1,11}, B^{(2)}_{1,12}, ..., B^{(K)}_{P,JJ}]'
  # Beta_array : array (P*J,J,K) of estimated coeffficients. [,,k] contains the B for class k. In [,,k], the first [1:J,1:J] contains the B at lag p=1 in this order c(row1: B_11,.., B_1J; row2: B_21, ..., B_2J; ...; rowJ: B_J1, ..., B_JJ). In [,,k], the obs [1:J, (J(P-1)+1):JP] contain the at lag p=P in the order explained before.
  # CFIT: list containing the C matrix and the CNorm built with the function "buildCmatrix"
  # J : number of time series in each class.
  # K : number of classes.
  # P : VAR order.
  # n : time series length.
  # Y_array : array (n,J,K) of original responses. [,,k] contains the responses for class k. In [,,k], column j contains the responses of variable j from 1:n.
  # X_array : array (n,P*J,K) of original inputs. [,,k] contains the inputs for class k. In [,,k], the first [1:n,1:J] contains the inputs at lag p=1: column j, contains the inputs from (P-p+1):(T-p). In [,,k], the [1:n, (J(P-1)+1):JP] contains the J inputs at lag P.
  # E_array : array (n,J,K) of errors. [,,k] contains the errors for class k. In [,,k], column j contains the error of variable j from 1:n.
  # omega_new_list : list of J x J inverse error covariance matrices
  # critvalue: value of information criterion
  # lambda1_opt : regularization parameter on the autoregressive coeffcients
  # lambda2_opt : regularization parameter on the difference between autoregressive coeffcients
  # gamma1_opt : regularization parameter on the elements of inverse error covariance matrix
  # gamma2_opt : regularization parameter on the difference between elements of inverse error covariance matrix
  # lambda_RIDGE : regularization paramater for RIDGE
  # beta_RIDGE_arr : array of ridge estimates

  #### Additional auxiliary functions #####
  # Spectral Decomposition
  spectral_dec<-function(Omega){
    decomp<-eigen(Omega, symmetric = TRUE)
    CC <- decomp$vectors
    Lambda <- diag(sqrt(decomp$values))
    P <- CC %*% Lambda %*% t(CC)
  }

  # Opposite of is.null
  is.not.null <- function(x) ! is.null(x)

  # Opposite of is.array
  is.not.array <- function(x) ! is.array(x)

  # Log-likelihood of the Multi-class VAR
  LogLik_function <- function (X_list, Y_stack_list,  Beta_new_list, omega_tilde_list, omega_new_list, n, J){
    AA <- (Y_stack_list - X_list %*% Beta_new_list)
    ll<- ((t(AA) %*% omega_tilde_list %*% AA) - (n*J*log(det(omega_new_list))))
  }
  # Kronecker product
  transpose_function<-function(U,J){
    kronecker(diag(1,J),t(U)%*%U)
  }


  #### START CODE ####

  # Preliminary steps #
  # Check Data
  if (is.null(Data) & (is.null(Y) | is.null(X) | is.null(X0) | is.null(J) | is.null(K) | is.null(n)) ){
    stop("Either give as input Data or give as input X, Y, X0, J, K, n.")
  }

  # Generate the inputs
  if (is.not.null(Data)){
    INPUTS<-Inputs_MultiClass_VAR(Data, P=P,calculate.C=calculate.C)
    Y <- INPUTS$Y
    X <- INPUTS$X
    X0 <- INPUTS$X0
    XX_prod <- INPUTS$XX
    XY_prod <- INPUTS$XY

    # Dimensions
    Ndata<-dim(Data)[1]
    J<-dim(Data)[2]
    K<-dim(Data)[3]
    n<-Ndata-P
  }

  if ((calculate.C==F & (is.null(C) | is.null(CNorm)))  ){
    stop("Either set calculate.C=TRUE or give the input C and CNorm.")
  }

  if (calculate.C==TRUE){
    C <- INPUTS$C
    CNorm <- INPUTS$CNorm
    CFIT<-list("C"=INPUTS$C,"CNorm"=INPUTS$CNorm)
  }

  if (is.not.null(C) & is.not.null(CNorm)){
    CFIT<-list("C"=C,"CNorm"=CNorm)
  }

  # Additional inputs
  if(is.null(b_init)){
    b_init<-array(0,c(ncol(X),1))
  }
  if(is.null(XX_prod)){
    X0_list<-alply(X0,3,.dims=TRUE)
    XtX<-lapply(X0_list,transpose_function,J=J)
    XX_prod<-bdiag(XtX)
  }

  if(is.null(XY_prod)){
    XY_prod<-t(X)%*%Y
  }

  if (type=="GrLasso" & is.null(group)){
    stop("Set the grouping structure.")
  }

  if (is.not.null(lambdaRIDGE_OPT) & is.not.array(lambdaRIDGE_OPT)){
    stop("lambdaRIDGE_OPT must be an array c(1,1,K).")
  }

  if (lambda_weights==T & is.null(Clusterinfo)){
    stop("Insert a value for Clusterinfo.")
  }


  # Check the optimal regularization paramters
  if ((is.not.null(lambda1_OPT) | is.not.null(lambda2_OPT) | is.not.null(gamma1_OPT) | is.not.null(gamma2_OPT))
      & (is.not.null(lambda1_steps) | is.not.null(lambda2_steps) | is.not.null(gamma1_steps) | is.not.null(gamma2_steps)) ){
    stop("Either set the optimal values of lambda&gamma or set the values for the grid search.")
  }
  if ((is.null(lambda1_OPT) & is.null(lambda2_OPT) & is.null(gamma1_OPT) & is.null(gamma2_OPT))
      & (is.null(lambda1_steps) & is.null(lambda2_steps) & is.null(gamma1_steps) & is.null(gamma2_steps))
      & (is.null(lambda1_min) & is.null(lambda1_max) & is.null(gamma1_min) & is.null(gamma1_max))
      & (is.null(lambda2_min) & is.null(lambda2_max) & is.null(gamma2_min) & is.null(gamma2_max)) ){
    stop("Either set the optimal values of lambda&gamma or set the values for the grid search.")
  }
  if ( ((type=="AdLasso") | (type=="GrLasso")) & (is.null(lambdaRIDGE_OPT) & is.null(lambdaRIDGE_min) & is.null(lambdaRIDGE_max) & is.null(lambdaRIDGE_steps)) ){
    stop("Either set the optimal values of lambdaRIDGE or set the values for the grid search.")
  }
  if ( is.not.null(lambdaRIDGE_OPT) & (is.not.null(lambdaRIDGE_min) | is.not.null(lambdaRIDGE_max) & is.not.null(lambdaRIDGE_steps) )){
    stop("Either set the optimal values of lambdaRIDGE or set the values for the grid search.")
  }

  if (is.null(lambda1_OPT) & is.not.null(lambda2_OPT)){
    stop("Set the optimal value of lambda1.")
  }
  if (is.not.null(lambda1_OPT) & is.null(lambda2_OPT)){
    stop("Set the optimal value of lambda2.")
  }
  if (is.null(gamma1_OPT) & is.not.null(gamma2_OPT)){
    stop("Set the optimal value of gamma1.")
  }
  if (is.not.null(gamma1_OPT) & is.null(gamma2_OPT)){
    stop("Set the optimal value of gamma2.")
  }


  # Check the grid of the regularization paramters
  if ((is.null(lambda1_min) | is.null(lambda1_max) | is.null(lambda1_steps) ) & (is.not.null(lambda1_min) | is.not.null(lambda1_max) | is.not.null(lambda1_steps))){
    stop("Set all the values of the lambda1 grid.")
  }
  if ((is.null(lambda2_min) | is.null(lambda2_max) | is.null(lambda2_steps) ) & (is.not.null(lambda2_min) | is.not.null(lambda2_max) | is.not.null(lambda2_steps))){
    stop("Set all the values of the lambda2 grid.")
  }
  if ((is.null(gamma1_min) | is.null(gamma1_max) | is.null(gamma1_steps) ) & (is.not.null(gamma1_min) | is.not.null(gamma1_max) | is.not.null(gamma1_steps))){
    stop("Set all the values of the gamma1 grid.")
  }
  if ((is.null(gamma2_min) | is.null(gamma2_max) | is.null(gamma2_steps) ) & (is.not.null(gamma2_min) | is.not.null(gamma2_max) | is.not.null(gamma2_steps))){
    stop("Set all the values of the gamma2 grid.")
  }

  # Set the grid for the regularization parameters
  if (is.not.null(lambda1_min) & is.not.null(lambda1_max) & is.not.null(lambda1_steps) & is.not.null(lambda2_min) & is.not.null(lambda2_max) & is.not.null(lambda2_steps) ){
    if(type=="GrLasso"){
      lambda1_max<-lambda1_max/lambda1_max
    }
    lambda1<-seq(from=lambda1_min, to=lambda1_max, length=lambda1_steps)
    lambda2<-seq(from=lambda2_min, to=lambda2_max, length=lambda2_steps)
    l1<-list(lambda1_l=sort(rep(lambda1,length(lambda2))), lambda2_l=rep(lambda2,length(lambda1)))
  }

  if (is.not.null(gamma1_min) & is.not.null(gamma1_max) & is.not.null(gamma1_steps) & is.not.null(gamma2_min) & is.not.null(gamma2_max) & is.not.null(gamma2_steps)){
    gamma1<-seq(from=gamma1_min,to=gamma1_max,length=gamma1_steps)
    gamma2<-seq(from=gamma2_min,to=gamma2_max,length=gamma2_steps)
    g1<-list(gamma1_l=sort(rep(gamma1,length(gamma2))), gamma2_l=rep(gamma2,length(gamma1)))
  }


  # Data preliminaries
  Ynew <- Y; YData<-Y;rm(Y)
  Xnew <- X; XData<-X;rm(X)
  XX_prodnew<-XX_prod;rm(XX_prod)
  XY_prodnew<-XY_prod;rm(XY_prod)
  maxeig <- eigs(XX_prodnew, k=1, opts = list(retvec = FALSE))

  Obj=rep(1,maxit.both)
  Obj_Conv=rep(1,maxit.both)
  iter<-2

  while (Obj_Conv[iter]>tol.both & iter<maxit.both) {

    #########
    ## SPG ##
    #########

    # Estimation of the autoregressive coeffcients using SPG algorithm

      #########
      # Lasso #
      #########

      if (type=="Lasso"){
        if (is.null(lambda1_OPT) &  is.null(lambda2_OPT)){
          # Use BIC as criterion
          if (criterion == "BIC"){
            morearg=list(Y=Ynew, X=Xnew, XX_prod=XX_prodnew,XY_prod=XY_prodnew, J=J, P=P,  C=C, CNorm=CNorm,
                         maxiter=maxit.beta, b_init=b_init, mu=mu, tol=tol.beta, n=n,maxeig=maxeig,
                         lambda2weights_SPG=lambda_weights, Clusterinfo_SPG=Clusterinfo,K=K)
            if (iter==2){
              BIC_values<-mapply(BICfunction,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
            } else{
              BIC_values<-mapply(BICfunction,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
              BIC_values<-c(do.call("cbind",BIC_values))
              BIC_values<-c(matrix(BIC_values[[1]]))
            }

            which.min(BIC_values)
            lambda1_opt<-l1$lambda1_l[which.min(BIC_values)]
            lambda2_opt<-l1$lambda2_l[which.min(BIC_values)]
          }

          # Use AICc as criterion
          if (criterion == "AICc"){
            morearg=list(Y=Ynew, X=Xnew, XX_prod=XX_prodnew,XY_prod=XY_prodnew,C=C, CNorm=CNorm,
                         maxiter=maxit.beta, b_init=b_init, mu=mu, tol=tol.beta, n=n,maxeig=maxeig, J=J, K=K, P=P,
                         lambda2weights_SPG=lambda_weights, Clusterinfo_SPG=Clusterinfo)
            if (iter==2){
              AICc_values<-mapply(AICcfunction,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
            } else{
              AICc_values<-mapply(AICcfunction,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
              AICc_values<-c(do.call("cbind",AICc_values))
              AICc_values<-c(matrix(AICc_values[[1]]))
            }

            which.min(AICc_values)
            lambda1_opt<-l1$lambda1_l[which.min(AICc_values)]
            lambda2_opt<-l1$lambda2_l[which.min(AICc_values)]
          }
        }

        if (is.not.null(lambda1_OPT) &  is.not.null(lambda2_OPT)){
          lambda1_opt<-lambda1_OPT
          lambda2_opt<-lambda2_OPT
        }


        #  SPG estimation
        spg<-SPG(Y=Ynew, X=Xnew, XX_prod=XX_prodnew,XY_prod=XY_prodnew, J=J, P=P,
                 lambda1=lambda1_opt, lambda2=lambda2_opt, C=C, CNorm=CNorm,
                 maxiter=maxit.beta, b_init=b_init, mu=mu, tol=tol.beta,maxeig=maxeig, type=type,
                 lambda2weights_SPG=lambda_weights, Clusterinfo_SPG=Clusterinfo)

      } # End of type=="Lasso"

      ##################
      # Adaptive Lasso #
      ##################

      if (type=="AdLasso"){

        # RIDGE estimator for each class --> weight for the Adaptive Lasso
        Y_arr <- array(Ynew, c((n*J),1,K))    # Array nJx1xK: In [,,k], you have the first [1:n] obs of series j=1, [n+1:2n] obs of series j=2...
        Y_list<-alply(Y_arr,3,.dims=TRUE);rm(Y_arr)     # List of K matrices of dimension nJx1
        X0_list<-alply(X0,3,.dims=TRUE)       # List of nxJ matrices
        X_list<-lapply(X0_list,function(U){kronecker(diag(J),U)});rm(X0_list) # Knonecker diag(J) and X0_list

        # Ridge estimator with the optimal regularization parameters
        if (is.not.null(lambdaRIDGE_OPT)){
          lambdaRIDGE_OPT_list<-alply(lambdaRIDGE_OPT,3,.dims=TRUE)       # List of Optimal RIDGE regularization parameters
          mr<-list(J=J, P=P, n=n,  criterion=criterion)
          LIST_RIDGE<-mapply(RIDGE, Y=Y_list, X=X_list, l_ridge_opt=lambdaRIDGE_OPT_list, MoreArgs = mr)
        }
        # Ridge estimator with the regularization parameters selected via BIC or AICc
        if (is.not.null(lambdaRIDGE_min) & is.not.null(lambdaRIDGE_max) & is.not.null(lambdaRIDGE_steps)){
          mr<-list(J=J, P=P, n=n, l_ridge_min=lambdaRIDGE_min, l_ridge_max=lambdaRIDGE_max, l_ridge_steps=lambdaRIDGE_steps,
                   criterion=criterion)
          LIST_RIDGE<-mapply(RIDGE, Y=Y_list, X=X_list, MoreArgs = mr)
        }
        # Ridge output
        beta_RIDGE_arr<-array(NA, c((P*J^2), 1, K))
        lambda_RIDGE<-array(NA, c(1,1,K))
        for (ik in 1:K){
          beta_RIDGE_arr[,,ik]<-LIST_RIDGE[[2*ik-1]]
          lambda_RIDGE[,,ik]<-LIST_RIDGE[[2*ik]]
        }
        beta_RIDGE_stacked<-matrix(stack(data.frame(beta_RIDGE_arr))[,1], nrow = K*P*J^2, ncol=1)
        if (is.null(gamma_ridge)){
          gamma_ridge=1
        }
        weight_RIDGE<-(abs(beta_RIDGE_stacked))^gamma_ridge # Stacked coefficients. The first PJ^2 elements belong to k=1, ...

        if (is.null(lambda1_OPT)  &  is.null(lambda2_OPT)){
          # Select via BIC the regularization paramaters for the Adaptive Lasso
          if (criterion == "BIC"){
            morearg=list(Y=Ynew, X=Xnew, XX_prod=XX_prodnew,XY_prod=XY_prodnew,C=C, CNorm=CNorm, J=J, P=P,
                         maxiter=maxit.beta, b_init=b_init, mu=mu, tol=tol.beta, n=n,maxeig=maxeig, weight=weight_RIDGE,
                         lambda2weights_SPG=lambda_weights, Clusterinfo_SPG=Clusterinfo,K=K)
            if (iter==2){
              BIC_values<-mapply(BICfunction_AdLasso,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
            } else{
              BIC_values<-mapply(BICfunction_AdLasso,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
              BIC_values<-c(do.call("cbind",BIC_values))
              BIC_values<-c(matrix(BIC_values[[1]]))
            }

            which.min(BIC_values)
            lambda1_opt<-l1$lambda1_l[which.min(BIC_values)]
            lambda2_opt<-l1$lambda2_l[which.min(BIC_values)]
          }

          # Select via AICc the regularization paramaters for the Adaptive Lasso
          if (criterion == "AICc"){
            morearg=list(Y=Ynew, X=Xnew, XX_prod=XX_prodnew,XY_prod=XY_prodnew,C=C, CNorm=CNorm,
                         maxiter=maxit.beta, b_init=b_init, mu=mu, tol=tol.beta, n=n,maxeig=maxeig, J=J, K=K, P=P,
                         weight=weight_RIDGE, lambda2weights_SPG=lambda_weights, Clusterinfo_SPG=Clusterinfo)
            if (iter==2){
              AICc_values<-mapply(AICcfunction_AdLasso,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
            } else{
              AICc_values<-mapply(AICcfunction_AdLasso,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
              AICc_values<-c(do.call("cbind",AICc_values))
              AICc_values<-c(matrix(AICc_values[[1]]))
            }

            which.min(AICc_values)
            lambda1_opt<-l1$lambda1_l[which.min(AICc_values)]
            lambda2_opt<-l1$lambda2_l[which.min(AICc_values)]
          }
        }

        if (is.not.null(lambda1_OPT) &  is.not.null(lambda2_OPT)){
          lambda1_opt<-lambda1_OPT
          lambda2_opt<-lambda2_OPT
        }

        # SPG estimation for the Adaptive Lasso
        spg<-SPG(Y=Ynew, X=Xnew, XX_prod=XX_prodnew,XY_prod=XY_prodnew, J=J, P=P,
                 lambda1=lambda1_opt, lambda2=lambda2_opt, C=C, CNorm=CNorm,
                 maxiter=maxit.beta, b_init=b_init, mu=mu, tol=tol.beta,maxeig=maxeig, type=type, weight=weight_RIDGE,
                 lambda2weights_SPG=lambda_weights, Clusterinfo_SPG=Clusterinfo)


      } # End of type=="AdLasso"

      ###############
      # Group Lasso #
      ###############

      if (type=="GrLasso"){
        # RIDGE estimator for each class --> weight for the Group Lasso
        Y_arr <- array(Ynew, c((n*J),1,K))    # Array nJx1xK: In [,,k], you have the first [1:n] obs of series j=1,  [n+1:2n] obs of series j=2...
        Y_list<-alply(Y_arr,3,.dims=TRUE)     # List of K matrices of dimension nJx1
        X0_list<-alply(X0,3,.dims=TRUE)       # List of nxJP matrices
        X_list<-lapply(X0_list,function(U){kronecker(diag(J),U)}) # Knonecker diag(J) and X0_list

        # Ridge estimator with the optimal regularization parameters
        if (is.not.null(lambdaRIDGE_OPT)){
          lambdaRIDGE_OPT_list<-alply(lambdaRIDGE_OPT,3,.dims=TRUE)       # List of Optimal RIDGE regularization parameters
          mr<-list(J=J, P=P, n=n,  criterion=criterion)
          LIST_RIDGE<-mapply(RIDGE, Y=Y_list, X=X_list, l_ridge_opt=lambdaRIDGE_OPT_list, MoreArgs = mr)
        }
        # Ridge estimator with the regularization parameters selected via BIC or AICc
        if (is.not.null(lambdaRIDGE_min) & is.not.null(lambdaRIDGE_max) & is.not.null(lambdaRIDGE_steps)){
          mr<-list(J=J, P=P, n=n, l_ridge_min=lambdaRIDGE_min, l_ridge_max=lambdaRIDGE_max, l_ridge_steps=lambdaRIDGE_steps,
                   criterion=criterion)
          LIST_RIDGE<-mapply(RIDGE, Y=Y_list, X=X_list, MoreArgs = mr)
        }
        # Ridge output
        beta_RIDGE_arr<-array(NA, c((P*J^2), 1, K))
        lambda_RIDGE<-array(NA, c(1,1,K))
        for (ik in 1:K){
          beta_RIDGE_arr[,,ik]<-LIST_RIDGE[[2*ik-1]]
          lambda_RIDGE[,,ik]<-LIST_RIDGE[[2*ik]]
        }
        beta_RIDGE_stacked<-matrix(stack(data.frame(beta_RIDGE_arr))[,1], nrow = K*P*J^2, ncol=1)
        if (is.null(gamma_ridge)){
          gamma_ridge=1
        }
        beta_RIDGE<-beta_RIDGE_stacked # Stacked coefficients. The first PJ^2 elements belong to k=1, ...


        if (is.null(lambda1_OPT)  &  is.null(lambda2_OPT)){
          # Select via BIC the regularization paramaters for the Group Lasso
          if (criterion == "BIC"){
            morearg=list(Y=Ynew, X=Xnew, XX_prod=XX_prodnew,XY_prod=XY_prodnew,C=C, CNorm=CNorm,J=J, P=P,
                         maxiter=maxit.beta, b_init=b_init, mu=mu, tol=tol.beta, n=n,maxeig=maxeig,
                         weight=beta_RIDGE, group=group, lambda2weights_SPG=lambda_weights,Clusterinfo_SPG=Clusterinfo,K=K)
            if (iter==2){
              BIC_values<-mapply(BICfunction_GrLasso,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
            } else{
              BIC_values<-mapply(BICfunction_GrLasso,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
              BIC_values<-c(do.call("cbind",BIC_values))
              BIC_values<-c(matrix(BIC_values[[1]]))
            }

            which.min(BIC_values)
            lambda1_opt<-l1$lambda1_l[which.min(BIC_values)]
            lambda2_opt<-l1$lambda2_l[which.min(BIC_values)]
          }

          # Select via AICc the regularization paramaters for the Group Lasso
          if (criterion == "AICc"){
            morearg=list(Y=Ynew, X=Xnew, XX_prod=XX_prodnew,XY_prod=XY_prodnew,C=C, CNorm=CNorm,
                         maxiter=maxit.beta, b_init=b_init, mu=mu, tol=tol.beta, n=n,maxeig=maxeig, P=P, J=J, K=K,
                         weight=beta_RIDGE, group=group, lambda2weights_SPG=lambda_weights, Clusterinfo_SPG=Clusterinfo)
            if (iter==2){
              AICc_values<-mapply(AICcfunction_GrLasso,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
            } else{
              AICc_values<-mapply(AICcfunction_GrLasso,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
              AICc_values<-c(do.call("cbind",AICc_values))
              AICc_values<-c(matrix(AICc_values[[1]]))
            }

            which.min(AICc_values)
            lambda1_opt<-l1$lambda1_l[which.min(AICc_values)]
            lambda2_opt<-l1$lambda2_l[which.min(AICc_values)]
          }
        }

        if (is.not.null(lambda1_OPT) &  is.not.null(lambda2_OPT)){
          lambda1_opt<-lambda1_OPT
          lambda2_opt<-lambda2_OPT
        }

        # SPG estimation for the Group Lasso

        spg<-SPG(Y=Ynew, X=Xnew, XX_prod=XX_prodnew,XY_prod=XY_prodnew, J=J, P=P,
                 lambda1=lambda1_opt, lambda2=lambda2_opt, C=C, CNorm=CNorm,
                 maxiter=maxit.beta, b_init=b_init, mu=mu, tol=tol.beta,maxeig=maxeig,
                 type=type, weight=beta_RIDGE, group=group, lambda2weights_SPG=lambda_weights, Clusterinfo_SPG=Clusterinfo)


      } # End of type=="GrLasso"


    # Coefficients from the SPG
    Beta_new<-spg$beta                          # autoregressive coeffcients
    Beta_new_arr<-array(Beta_new, c(P*(J^2),1,K))
    Beta_new_list<-alply(Beta_new_arr,3,.dims=TRUE)
    iter_spg<-spg$it                            # iteration
    rm(spg);rm(Ynew);rm(Xnew)

    # Residuals
    e = (YData-XData%*%Beta_new)                        # stacked from the SPG
    e_array<-array(e, c(n,J,K))
    e_array_list<-alply(e_array,3,.dims=TRUE)   # create a list from the array
    S<-array(apply(e_array,3,cov),c(J,J,K))     # VAR-COV matrices
    S_list<-alply(S,3,.dims=TRUE)               # create a list from the array


    #########
    ## JGL ##
    #########

    # Estimation of the inverse error covariance matrix using JGL package

    if (is.null(gamma1_OPT)  &  is.null(gamma2_OPT)){
      # Use BIC as criterion
      if (criterion == "BIC"){
        morearg=list(Y=e_array_list, S_list=S_list, penalty="fused", n=n, J=J, maxit = maxit.omega, tol = tol.omega)

        BIC_values_JGL<-mapply(BICfunction_JGL,g1$gamma1_l, g1$gamma2_l, MoreArgs = morearg)

        which.min(BIC_values_JGL)
        gamma1_opt<-g1$gamma1_l[which.min(BIC_values_JGL)]
        gamma2_opt<-g1$gamma2_l[which.min(BIC_values_JGL)]
      }

      # Use the AICc as criterion
      if (criterion == "AICc"){
        morearg=list(Y=e_array_list, S_list=S_list, penalty="fused", n=n, J=J, K=K, maxit = maxit.omega, tol = tol.omega)

        AICc_values_JGL<-mapply(AICcfunction_JGL,g1$gamma1_l, g1$gamma2_l, MoreArgs = morearg)

        which.min(AICc_values_JGL)
        gamma1_opt<-g1$gamma1_l[which.min(AICc_values_JGL)]
        gamma2_opt<-g1$gamma2_l[which.min(AICc_values_JGL)]
      }
    }

    if (is.not.null(gamma1_OPT) &  is.not.null(gamma2_OPT)){
      gamma1_opt<-gamma1_OPT
      gamma2_opt<-gamma2_OPT
    }

    # JGLestimation
    jgl<-JGL(Y=e_array_list, penalty="fused", lambda1=gamma1_opt, lambda2=gamma2_opt,
             maxiter = maxit.omega, tol =tol.omega, return.whole.theta = TRUE )

    omega_new_list<-jgl$theta   # list of precision matrices

    # Spectral Decomposition on Omega
    P_list<-lapply(omega_new_list, spectral_dec)
    Ptilde_list<-lapply(P_list, kronecker, diag(n))

    # Define the new elements
    Y_arr <- array(YData, c((n*J),1,K))
    Y_list <- alply(Y_arr,3,.dims=TRUE)
    Y_new_list <- lapply(1:K,function(i,m,v){v[[i]]%*%m[[i]]},m=Y_list,v=Ptilde_list)
    Y_new_data <- as.data.frame(Y_new_list)
    Y_new1 <- stack(Y_new_data)
    Y_new <- data.matrix(Y_new1[,1])

    X0_list<-alply(X0,3,.dims=TRUE)
    X_list<-lapply(X0_list,function(U){kronecker(diag(J),U)})

    X_new_list <- lapply(1:K,function(i,m,v){v[[i]]%*%m[[i]]},m=X_list,v=Ptilde_list)
    X_new <- bdiag(X_new_list)


    # Check  Convergence  #
    # Objective function
    Y_stack_arr <- array(YData, c((n*J),1,K))
    Y_stack_list <- alply(Y_stack_arr, 3,.dims=TRUE)
    omega_tilde_list <- lapply(omega_new_list, kronecker, diag(n))
    Morearg = list(n=n, J=J)
    LL <- mapply(LogLik_function, X_list=X_list, Y_stack_list=Y_stack_list, Beta_new_list=Beta_new_list,
                 omega_tilde_list=omega_tilde_list, omega_new_list=omega_new_list,MoreArgs = Morearg)
    Obj[iter] <- sum(LL)
    Obj_Conv[iter+1]<-(abs(Obj[iter]-Obj[iter-1])/Obj[iter-1])

    # Create the new INPUTS
    iter=iter+1
    Ynew <- Y_new
    Xnew <- X_new

    XtX<-lapply(X_new_list,function(U){t(U)%*%U})
    XX_prodnew<-bdiag(XtX)
    XY_prodnew<-t(Xnew)%*%Ynew
    maxeig <- eigs(XX_prodnew, k=1, opts = list(retvec = FALSE))

  } # end while loop

  # OUTPUT #

  Beta_array<-aperm(array(Beta_new, c(P*J,J,K)), c(2,1,3))  # Array (P*J,J,K) of estimated coeffficients. [,,k] contains the B for class k. In [,,k], the first [1:J,1:J] contains the B at lag p=1 in this order c(row1: B_11,.., B_1J; row2: B_21, ..., B_2J; ...; rowJ: B_J1, ..., B_JJ). In [,,k], the obs [1:J, (J(P-1)+1):JP] contain the at lag p=P in the order explained before.
  Y_array<-array(YData, c(n,J,K))                           # Array (n,J,K) of original responses. [,,k] contains the responses for class k. In [,,k], column j contains the responses of variable j from 1:n.
  X_array<-X0                                               # Array (n,P*J,K) of original inputs. [,,k] contains the inputs for class k. In [,,k], the first [1:n,1:J] contains the inputs at lag p=1: column j, contains the inputs from (P-p+1):(T-p). In [,,k], the [1:n, (J(P-1)+1):JP] contains the J inputs at lag P.
  e_arr = (YData-XData%*%Beta_new)
  E_array<-array(e_arr, c(n,J,K))                           # Array (n,J,K) of errors. [,,k] contains the errors for class k. In [,,k], column j contains the error of variable j from 1:n.

  if(criterion=="BIC"){
    # Compute the BIC
    Log_lik<-(1/(n*J*K))*t((Ynew-Xnew%*%Beta_new))%*%(Ynew-Xnew%*%Beta_new)
    df<-length(which(Beta_new!=0))
    critvalue<-(-2*Log_lik+log(n*J*K)*df)
  }

  if(criterion=="AICc"){
    # Compute AICc
    Log_lik<-(1/(n*J*K))*t((Ynew-Xnew%*%Beta_new))%*%(Ynew-Xnew%*%Beta_new)
    df<-length(which(Beta_new!=0))
    AIC_value<-(-2*Log_lik+2*df)
    correction<-2*df*(df+1)/(n*J*K-df-1)
    critvalue<-AIC_value+correction
  }

  if (type=="Lasso"){
    lambda_RIDGE=NA
    beta_RIDGE=NA
    beta_RIDGE_arr=NA
  }



MultiClass_VAR <- list(Beta_new=Beta_new, Beta_array=Beta_array, CFIT=CFIT,
                       J=J, K=K, P=P, n=n,
                       Y_array=Y_array, X_array=X_array, E_array=E_array,
                       omega_new_list=omega_new_list, critvalue=critvalue,
                       lambda1_opt=lambda1_opt, lambda2_opt=lambda2_opt,
                       gamma1_opt=gamma1_opt, gamma2_opt=gamma2_opt,
                       lambda_RIDGE=lambda_RIDGE, beta_RIDGE_arr=beta_RIDGE_arr)
}

Inputs_MultiClass_VAR<-function(Data, P,calculate.C=F){

  #### Function to genrate the inputs of the function MultiClass_VAR ####

  ##########
  # INPUTS #
  ##########
  # Data: a N*J*K array, where:
  #  - N is the times series length
  #  - J is the number of time series
  #  - K is the number of classes
  # P: order of VAR
  # calculate.C : logical for calculating the C matrix. Default is FALSE.

  ###########
  # OUTPUTS #
  ###########
  # Y: NJKx1 vector of responses
  # X: NJKxP(J^2) matrix of inputs
  # X0: NxJPxK array of inputs
  # XX: product X'X
  # XY: product X'Y
  # C: C matrix for SPG algorithm
  # CNorm: norm of C for the SPG algorithm


  ## -- Dimensions -- ##
  Ndata<-dim(Data)[1]
  J<-dim(Data)[2]
  K<-dim(Data)[3]
  ndata<-Ndata-P

  ## -- Stack Data -- ##
  # Y
  y<-apply(Data, 3, stack_y, burnin=NULL, P=P, N=Ndata) # Matrix ncol=K. In [,k=1], you have the first [1:n] obs of series j=1, [n+1:2n] obs of series j=2...
  Y_stacked<-stack(data.frame(y))[,1]                   # The first [1:nJ] obs about k=1 (of which the first [1:n] obs of series j=1, [n+1:2n] obs of series j=2...); the [nJ+1:2nJ] obs about k=2 (of which the first [1:n] obs of series j=1, [n+1:2n] obs of series j=2...);
  Y_out = data.matrix(Y_stacked)

  # X
  Xbig<-apply(Data, 3, stack_Xbig, burnin=NULL, P=P, N=Ndata,J=J) # Column 1 is k=1, column K is k=K
  X_arr<-array(Xbig, dim=c(ndata*J,(J^2)*P, K))                   # Array. [,,1] is the block diagonal matrix for k=1 ...
  X_List<-alply(X_arr, 3, .dims=TRUE)
  X_stacked<-do.call(adiag, X_List)                               # Matrix KJn x KPJ^2 witha block diagonal structure. First block is k=1, last block is k-K.
  X_out = data.matrix(X_stacked)

  # X0
  X0s<-apply(Data, 3, X0function, burnin=NULL, P=P, N=Ndata)
  X0<-array(X0s, dim=c(ndata,J*P,K))                              # Array. [,,k] is about class k: the first J columns are at lag p=1, the last J columns are at lag p=P. The rows go from P-p to N-p.

  # XX
  # Transpose
  transpose_function<-function(U,J){
    kronecker(diag(1,J),t(U)%*%U)
  }
  X0_list<-alply(X0,3,.dims=TRUE)
  XtX<-lapply(X0_list,transpose_function,J=J)
  XX_prod<-bdiag(XtX)

  # XY
  XY_prod<-t(X_out)%*%Y_out

  # C
  if(calculate.C==T){
    CFIT<-buildCmatrix(J=J, K=K, P=P)
    CNorm<-CFIT$CNorm
    C<-CFIT$C
  }else{
    CNorm<-NULL
    C<-NULL
  }

  ## -- Output -- ##
  Inputs_MultiClass_VAR<-list( "Y"=Y_out, "X"=X_out, "X0"=X0, "XX"=XX_prod, "XY"=XY_prod, "C"=C, "CNorm"=CNorm)
}

buildCmatrix<-function(J,K,P){

  ##### Function to build the C matrix ####

  #########
  # INPUT #
  #########
  # J : number of time series per class
  # K : number of classes
  # P : order of VAR

  ##########
  # OUTPUT #
  ##########
  # C : C matrix
  # CNorm : norm of the C matrix

  c1<-diag(1,P*(J)^2,P*(J)^2)
  c2<-diag(-1,P*(J)^2,P*(J)^2)
  Cmatrix<-c()
  for(k in 1:(K-1)){
    Cmatrixnew<-matrix(0,nrow=(K-k)*P*J^2,ncol=(K*P*J^2))

    for( i in 1:(K-k)){
      Cmatrixnew[ ((i-1)*P*J^2+1):(i*P*J^2), ((i-1)*P*J^2+1):(i*P*J^2)      ]<-c1
      Cmatrixnew[ ((i-1)*P*J^2+1):(i*P*J^2), (((i+k)-1)*P*J^2+1):((i+k)*P*J^2)     ]<-c2
    }

    Cmatrix<-rbind(Cmatrix,Cmatrixnew)

  }

  C_sparse<-Matrix(Cmatrix, sparse=TRUE)
  C_ct<-Conj(t(Cmatrix))
  C_sparse_ct<-Matrix(C_ct,sparse=TRUE)
  rm(C_ct)
  CX_sparse<-C_sparse_ct%*%C_sparse
  C_sparse_eig<-eigen(CX_sparse)
  rm(CX_sparse)
  CNorm_sparse<-sqrt(max(C_sparse_eig$values))
  rm(C_sparse_eig)
  CNorm<-CNorm_sparse

  out<-list(CNorm=CNorm,C=Cmatrix)
}

stack_y <- function(A,burnin,P,N){

  #### Functions to stack the data for the Multi-class VAR: Stacking Y ####

  #########
  # INPUT #
  #########
  # A : data
  # burnin : burnin to remove at the beginnnig. If not present, set to NULL
  # P : order of the VAR
  # N : time series length

  ##########
  # OUTPUT #
  ##########
  # a stacked vector Y

  if(is.null(burnin)){
    burnin = 1
  }
  a<-data.frame(A)
  a<-a[(burnin+P):N,]
  a_stacked<-stack(a) # Put one column below the other
  return(a_stacked[,1]) #
}

stack_Xbig <- function(A,burnin,P,N, J){

  #### Functions to stack the data for the Multi-class VAR: Stacking X ####

  #########
  # INPUT #
  #########
  # A : data
  # burnin : burnin to remove at the beginnnig. If not present, set to NULL
  # P : order of the VAR
  # N : time series length
  # J : number of time series

  ##########
  # OUTPUT #
  ##########
  # Xbig_mat : a block matrix X

  if(is.null(burnin)){
    burnin = 1
  }
  X0_arr<-array(NA, c((N-P-burnin+1), ncol(A), P)) # X0 for class k in array form
  for (ip in P:1){
    X0_arr[,,ip]<-A[(burnin+ip-1):(N-P+ip-1),] # [,,1]: elements at lag p=1, J columns and N rows going from P-p to N-p. [,,P] elements at lag p=P
  }
  X0_list<-alply(X0_arr, 3,.dims=TRUE) # Convert the array into a list
  X0<-do.call(cbind,X0_list) # Matrix; the first J columns are at lag p=1, the last J columns are at lag p=P. The rows go from P-p to N-p.
  Xbig_mat<-kronecker(diag(J), X0)
  return(Xbig_mat)
}

X0function <- function(A,burnin,P,N){
  #### Fuction to create X0 ####

  #########
  # INPUT #
  #########
  # A : data
  # burnin : burnin to remove at the beginnnig. If not present, set to NULL
  # P : order of the VAR
  # N : time series length

  ##########
  # OUTPUT #
  ##########
  # X0: a list of X matrix for all classes

  if(is.null(burnin)){
    burnin = 1
  }
  X0_arr<-array(NA, c((N-P-burnin+1), ncol(A), P)) # X0 for class k in array form
  for (ip in P:1){
    X0_arr[,,ip]<-A[(burnin+ip-1):(N-P+ip-1),] # [,,1]: elements at lag p=1, J columns and N rows going from P-p to N-p. [,,P] elements at lag p=P
  }
  X0_list<-alply(X0_arr, 3,.dims=TRUE) # Convert the array into a list
  X0<-do.call(cbind,X0_list) # Matrix; the first J columns are at lag p=1, the last J columns are at lag p=P. The rows go from P-p to N-p.

  return(X0)
}


"SPG" = function(Y, X, XX_prod, XY_prod, P=P, J=J, lambda1, lambda2, C, CNorm=NULL, maxiter,
                 b_init, mu, tol, maxeig, type, weight=NULL, group=NULL, lambda2weights_SPG=FALSE,
                 Clusterinfo_SPG=NULL){

  ##### Function for the estimation of the autoregressive coefficients using the SPG algorithm ####

  #########
  # INPUT #
  #########
  # Y : vector of dimension NJK x 1 of dependent variables
  # X : matrix of dimension NJK x J�PK of independent variables
  # XX_prod : matrix of dimension J�PK x J�PK
  # XY_prod : matrix of dimension J�PK x 1
  # J : number of variables in each class
  # P : order of VAR
  # lambda1 = regularization parameter on the autoregressive coeffcients (lasso)
  # lambda2 = regularization parameter on the difference between autoregressive coeffcients (fusion)
  # C : C matrix (built with 'Cmatrix.R')
  # CNorm : norm of C (built with 'Cmatrix.R')
  # maxiter ; maxiter SPG
  # b_init : vector of dimension J�PK x 1 with the initial values of the autoregressive coefficients
  # mu : smoothness parameter
  # tol : tolerance for SPG
  # maxeig : first eigenvalue of the matrix XX
  # type : type of penalty not included in the smooth function ("Lasso", "AdLasso", "GrLasso")
  # weight: weight for the Adaptive Lasso or the Group Lasso. Default is NULL.
  # group : grouping structure of the Group Lasso. Default is NULL.
  # lambda_weights_SPG : logical. If TRUE, then do the weighting based on Clusterinfo. Default is FALSE.
  # Clusterinfo_SPG : K x K matrix of additional clustering information. K is the number of classes. Default is NULL.

  ##########
  # OUTPUT #
  ##########
  # beta : vector of dimension J�PK x 1 with the autoregressive coefficients
  # obj : value of the objective function
  # it : iteration reached by SPG

  obj=rep(1,maxiter)
  obj_con=rep(1,maxiter)
  density=array(0,c(maxiter,1))
  beta=b_init
  W=beta
  theta=1

  a=nrow(Y); a
  b=ncol(X); b

  XX=XX_prod
  XY=XY_prod

  if (lambda2weights_SPG==F){
    C=lambda2*C
    L=maxeig$values+lambda2^2*CNorm/mu # Lipschitz-constant
  }


  if (lambda2weights_SPG==T){
    K<-nrow(Clusterinfo_SPG)
    weights_lambda2<-lambda2weights(X=Clusterinfo_SPG, K=K, P=P, J=J) # Vector of length (K-1)*J*J*P*K/2
    Cweight<-C/weights_lambda2
    C=lambda2*Cweight
    L=maxeig$values+(lambda2/min(weights_lambda2))^2*CNorm/(mu)  # Lipschitz-constant
  }

  # While loop
  it<-2
  while (obj_con[it]>tol & it<maxiter) {

    W1=(C%*%W)/mu

    # Hard-Threshold #
    ht<-hard_threshold(W1, 1)
    A<-hard_threshold(W1, 1)$A
    grad=XX%*%W-XY+t(C)%*%A   # Gradient as in eq.3.10 of Chen et al. (2012)

    # Step of size (1/L) in the direction of the Gradient
    V=W-(1/L)*grad # Dim PKJ^2 x 1: the first P*J^2 obs belong to k=1 in this order B_11, B_21, B_12, B_22 (i.e. stacked by columns)


    # Proximal Operator #
    ## -- Lasso -- ##
    if (type=="Lasso"){
      L1<-lambda1/L
      st<-soft_threshold(V,L1)
      beta_new<-st$B
      density[it,]<-st$nonzero/b
    }

    ## -- Adaptive Lasso -- ##
    if (type=="AdLasso"){
      L1<-lambda1/L
      L1_AdLasso<-L1%/%(weight)
      st<-soft_threshold_AdLasso(W=V,lambda1=L1_AdLasso)
      beta_new<-st$B
      density[it,]<-st$nonzero/b
    }

    ## -- Group Lasso -- ##
    if (type=="GrLasso"){
      L1<-lambda1/L
      st<-soft_threshold_GrLasso(W=V, b.ridge=weight, lambda1=L1, group=group)
      beta_new<-st$B
      if (any(is.na(beta_new))){
        stop("Problem in the thresholding function of the Group Lasso: beta_new is NA. Check whether group has length KPJ^2.")
      }
    }


    # Update #
    theta_new<-(sqrt(theta^4+4*theta^2)-theta^2)/2
    W=beta_new+(1-theta)/theta*theta_new*(beta_new - beta)

    obj[it]<-(sum((Y-X%*%beta_new)^2))/2+sum(abs(C%*%beta_new))+lambda1*sum(abs(beta_new))

    obj_con[it+1]<-(abs(obj[it]-obj[it-1])/obj[it-1])

    beta=beta_new
    theta=theta_new
    it=it+1

  } # end while loop

  SPG <- list(beta=beta, obj=obj, it=it)
}


"hard_threshold" = function(W, lambda1){
  #### Hard thresholding SPG algorithm ####

  #########
  # INPUT #
  #########
  # W : autoregressive coefficients
  # lambda1 : threshold

  ##########
  # OUTPUT #
  ##########
  # A : autoregressive coefficients

  m=nrow(W)
  n=ncol(W)
  len=m*n

  s<-matrix(NA,len,1)

  for (i in 1:len){
    if (W[i]>lambda1){
      s[i]=lambda1
    }
    else {
      if (W[i]<(-lambda1)){
        s[i]=-lambda1
      }
      else {
        s[i]=W[i]
      }
    }
  }


  hard_threshold <- list(A=s)
}

"soft_threshold" = function(W, lambda1){
  #### Soft thresholding in SPG algorithm for Lasso ####

  #########
  # INPUT #
  #########
  # W : autoregressive coefficients
  # lambda1 : threshold

  ##########
  # OUTPUT #
  ##########
  # B : autoregressive coefficients
  # nz : non-zero values

  m=nrow(W)
  n=ncol(W)
  len=m*n
  s<-matrix(NA,len,1)

  for (i in 1:len){
    if (W[i]>lambda1){
      s[i]=W[i]-lambda1
      nz=1
    }
    else {
      if (W[i]<(-lambda1)){
        s[i]=W[i]+lambda1
        nz=1
      }
      else {
        s[i]=0
        nz=0
      }
    }
  }

  soft_threshold <- list(B=s, nonzero=nz)

}

"soft_threshold_AdLasso" = function(W, lambda1){

  #### Soft thresholding in SPG algorithm for the Adaptive Lasso ####

  #########
  # INPUT #
  #########
  # W : autoregressive coefficients
  # lambda1 : threshold already weighted

  ##########
  # OUTPUT #
  ##########
  # B : autoregressive coefficients
  # nz : non-zero values

  m=nrow(W)
  n=ncol(W)
  len=m*n
  s<-matrix(NA,len,1)

  for (i in 1:len){
    if (W[i]>lambda1[i]){
      s[i]=W[i]-lambda1[i]
      nz=1
    }
    else {
      if (W[i]<(-lambda1[i])){
        s[i]=W[i]+lambda1[i]
        nz=1
      }
      else {
        s[i]=0
        nz=0
      }
    }
  }

  soft_threshold_AdLasso <- list(B=s, nonzero=nz)

}

"soft_threshold_GrLasso" = function(W, b.ridge, lambda1, group){

  #### Soft threshold in SPG algorithm for the Group Lasso ####

  #########
  # INPUT #
  #########
  # W : autoregressive coefficients
  # b.ridge : Ridge estimates of the autoregressive coefficients
  # lambda1 : threshold altready weighted
  # group : vector indicating the grouping structure

  ##########
  # OUTPUT #
  ##########
  # B : autoregressive coefficients


  m=nrow(W)
  n=ncol(W)
  len=m*n
  G<-max(group) # maximum number of groups

  s<-matrix(NA,len,1)

  for (i.g in 1:G){
    var.index<-which(group==i.g)    # variables belonging to group i.g
    bg.ridge<-b.ridge[var.index]
    bg.norm<-sqrt(sum(bg.ridge^2))
    p.var<-length(var.index)        # number of elements in group i.g
    lambdag<-lambda1*sqrt(p.var)/bg.norm
    Wg<-W[var.index]
    Wgnorm<-sqrt(sum(Wg^2))

    arg<-(1-lambdag/Wgnorm)
    if(arg<0){
      s[var.index]<-0
    }else{
      s[var.index]<-arg*Wg
    }
  }

  soft_threshold_AdLasso <- list(B=s)

}


"BICfunction" = function(Y, X, XX_prod, XY_prod, J, K, P, lambda1, lambda2,
                         C, CNorm, maxiter, b_init, mu, tol, n, maxeig, lambda2weights_SPG=NULL, Clusterinfo_SPG=NULL){

  #### Function for the choice of the regularization parameters based on BIC ####

  #########
  # INPUT #
  #########
  # Y : vector of dimension NJK x 1 of dependent variables
  # X : matrix of dimension NJK x J�PK of independent variables
  # XX_prod : matrix of dimension J�PK x J�PK
  # XY_prod : matrix of dimension J�PK x 1
  # J : number of variables in each class
  # K : number of classes
  # P : lag order of VAR
  # lambda1 : regularization parameter for B (lasso)
  # lambda2 : regularization parameter for B (fusion)
  # C : C matrix (built with 'Cmatrix.R')
  # CNorm : norm of C (built with 'Cmatrix.R')
  # maxiter : maxiter SPG
  # b_init : vector of dimension J�PK x 1 with the initial values of the autoregressive coefficients (if no default, take all zeros)
  # mu : smoothness parameter (if no default, take 1e-4)
  # tol : tolerance for SPG
  # n : time series length
  # maxeig : first eigenvalue of the matrix XX
  # lambda2weights_SPG : logical. If TRUE, then do the weighting based on Clusterinfo. Default is FALSE.
  # Clusterinfo_SPG : K x K matrix of additional clustering information. K is the number of classes. Default is NULL.

  ##########
  # OUTPUT #
  ##########
  # BIC_values: BIC values

  # Estimate the autoregressive coefficients using SPG algorithm
  FIT<-SPG(Y=Y, X=X, XX_prod=XX_prod,XY_prod=XY_prod, J=J, P=P, lambda1=lambda1, lambda2=lambda2, C=C, CNorm=CNorm,
           maxiter=maxiter, b_init=b_init, mu=mu, tol=tol,maxeig=maxeig, type="Lasso", lambda2weights_SPG=lambda2weights_SPG, Clusterinfo_SPG=Clusterinfo_SPG)

  # Compute the BIC
  Log_lik<-(1/(n*J*K))*t((Y-X%*%FIT$beta))%*%(Y-X%*%FIT$beta)
  df<-length(which(FIT$beta!=0))
  BIC_value<-(-2*Log_lik+log(n*J*K)*df)

  return(BIC_value)

}

"BICfunction_JGL" = function(Y, S_list, penalty, lambda1, lambda2, n, J, maxit, tol){

  #### BIC for the inverse error covariance matrix with JGL ####

  #########
  # INPUT #
  #########
  # Y : list of data matrices for JGL
  # S_list : list of inverse error covariance matrices
  # penalty : penalty of JGL (choose "fused")
  # lambda1 : regularization parameter on the elements of inverse error covariance matrix (lasso)
  # lambda2 : regularization parameter on the difference between elements of inverse error covariance matrix (fusion)
  # n : time series length
  # J : number of time series per class
  # maxit : maximum iteration for JGL
  # tol : tolerance for JGL

  ##########
  # OUTPUT #
  ##########
  # BIC_values: BIC values

  # Fit the JGL
  FIT<-JGL(Y=Y, penalty=penalty, lambda1=lambda1, lambda2=lambda2, maxiter = maxit, tol = tol, return.whole.theta = TRUE)

  more_arg=list(n=n, J=J)

  # Compute the BIC
  BIC_function<-function(A, B, J, n){
    nk<-n*J
    M<-length(which(B!=0))
    E<-(M-J)/2
    BIC<-(nk*(sum(diag(A%*%B)))-nk*log(det(B))+log(nk)*E)
  }

  BIC_value<-mapply(BIC_function, A=S_list, B=FIT$theta, MoreArgs = more_arg)
  BIC_values<-sum(BIC_value)
  BIC_values

  return(BIC_values)
}

"BICfunction_RIDGE" = function(l_ridge_BIC, Y_BIC, X_BIC, J_BIC, P_BIC, n_BIC){

  #### BIC for the autoregressive coeffcients with RIDGE  ####

  #########
  # INPUT #
  #########
  # l_ridge_BIC : grid of regularization paramters for the Ridge estimator
  # Y_BIC : vector of responses of length nJ
  # X_BIC : matrix of inputs of dimension nJ x PJ^2
  # J_BIC : number of time series
  # P_BIC : order of VAR
  # n_BIC : time series length


  ##########
  # OUTPUT #
  ##########
  # BIC_values: BIC values

  # Ridge estimator
  Ridge_Function<-function(RX, RY, RJ, RP, Rl_ridge){
    B_r<- solve(t(RX)%*%RX + Rl_ridge*diag(RP*RJ^2)) %*% (t(RX)%*%RY)
  }

  # Estimate the autoregressive coefficients using SPG algorithm
  B_r<-Ridge_Function(RY=Y_BIC, RX=X_BIC, RJ=J_BIC, RP=P_BIC, Rl_ridge=l_ridge_BIC)

  # Compute the BIC
  Log_lik<-(1/n_BIC)*t((Y_BIC-X_BIC%*%B_r))%*%(Y_BIC-X_BIC%*%B_r)
  df<-length(which(B_r!=0))/J_BIC
  BIC_value<-(-2*Log_lik+log(n_BIC)*df)

  return(BIC_value)

}

"BICfunction_AdLasso" = function(Y, X, XX_prod, XY_prod, J, K, P, lambda1, lambda2, C, CNorm, maxiter, b_init,
                                 mu, tol, n, maxeig, weight, lambda2weights_SPG=NULL, Clusterinfo_SPG=NULL){

  ##### BIC for for the autoregressive coeffcients with Adaptive Lasso ####


  #########
  # INPUT #
  #########
  # Y : vector of dimension NJK x 1 of dependent variables
  # X : matrix of dimension NJK x J�PK of independent variables
  # XX_prod : matrix of dimension J�PK x J�PK
  # XY_prod : matrix of dimension J�PK x 1
  # J : number of time series in each class
  # K : number of classes
  # P : order of VAR
  # lambda1 : regularization parameter for the autoregressive coeffcients (lasso)
  # lambda2 : regularization parameter for the autoregressive coeffcients (fusion)
  # C : C matrix (built with 'Cmatrix.R')
  # CNorm : norm of C (built with 'Cmatrix.R')
  # maxiter : maxiter SPG
  # b_init : vector of dimension J�PK x 1 with the initial values of the autoregressive coefficients (if no default, take all zeros)
  # mu : smoothness parameter (if no default, take 1e-4)
  # tol : tolerance for SPG
  # n : time series length
  # maxeig : first eigenvalue of the matrix XX
  # weight: weight for the Adaptive Lasso.
  # lambda2weights_SPG : logical. If TRUE, then do the weighting based on Clusterinfo. Default is FALSE.
  # Clusterinfo_SPG : K x K matrix of additional clustering information. K is the number of classes. Default is NULL.

  ##########
  # OUTPUT #
  ##########
  # BIC_values: BIC values

  # Estimate the autoregressive coefficients using SPG algorithm
  FIT<-SPG(Y=Y, X=X, XX_prod=XX_prod,XY_prod=XY_prod, J=J, P=P, lambda1=lambda1, lambda2=lambda2, C=C, CNorm=CNorm,
           maxiter=maxiter, b_init=b_init, mu=mu, tol=tol,maxeig=maxeig, type="AdLasso", weight=weight,
           lambda2weights_SPG=lambda2weights_SPG, Clusterinfo_SPG=Clusterinfo_SPG)

  # Compute the BIC
  Log_lik<-(1/(n*J*K))*t((Y-X%*%FIT$beta))%*%(Y-X%*%FIT$beta)
  df<-length(which(FIT$beta!=0))
  BIC_value<-(-2*Log_lik+log(n*J*K)*df)

  return(BIC_value)

}

"BICfunction_GrLasso" = function(Y, X, XX_prod, XY_prod, J, K, P, lambda1, lambda2, C, CNorm, maxiter, b_init,
                                 mu, tol, n, maxeig, weight, group, lambda2weights_SPG=NULL, Clusterinfo_SPG=NULL){

  ##### BIC for for the autoregressive coeffcients with Group Lasso ####

  #########
  # INPUT #
  #########
  # Y : vector of dimension NJK x 1 of dependent variables
  # X : matrix of dimension NJK x J�PK of independent variables
  # XX_prod : matrix of dimension J�PK x J�PK
  # XY_prod : matrix of dimension J�PK x 1
  # J : number of time series in each class
  # K : number of classes
  # P : lag order of VAR
  # lambda1 : regularization parameter for the autoregressive coeffcients (lasso)
  # lambda2 : regularization parameter for the autoregressive coeffcients (fusion)
  # C : C matrix (built with 'Cmatrix.R')
  # CNorm : norm of C (built with 'Cmatrix.R')
  # maxiter : maxiter SPG
  # b_init : vector of dimension J�PK x 1 with the initial values of the autoregressive coefficients (if no default, take all zeros)
  # mu : smoothness parameter (if no default, take 1e-4)
  # tol : tolerance for SPG
  # n : time series length
  # maxeig : first eigenvalue of the matrix XX
  # weight: weight for the Adaptive Lasso.
  # group : grouping structure of the Group Lasso
  # lambda2weights_SPG : logical. If TRUE, then do the weighting based on Clusterinfo. Default is FALSE.
  # Clusterinfo_SPG : K x K matrix of additional clustering information. K is the number of classes. Default is NULL.

  ##########
  # OUTPUT #
  ##########
  # BIC_values: BIC values

  # Estimate the autoregressive coefficients using SPG algorithm
  FIT<-SPG(Y=Y, X=X, XX_prod=XX_prod,XY_prod=XY_prod, J=J, P=P, lambda1=lambda1, lambda2=lambda2, C=C, CNorm=CNorm,
           maxiter=maxiter, b_init=b_init, mu=mu, tol=tol,maxeig=maxeig, type="GrLasso", weight=weight, group=group,
           lambda2weights_SPG=lambda2weights_SPG, Clusterinfo_SPG=Clusterinfo_SPG)

  # Compute the BIC
  Log_lik<-(1/(n*J*K))*t((Y-X%*%FIT$beta))%*%(Y-X%*%FIT$beta)
  df<-length(which(FIT$beta!=0))
  BIC_value<-(-2*Log_lik+log(n*J*K)*df)

  return(BIC_value)

}

"AICcfunction" = function(Y, X, XX_prod, XY_prod, lambda1, lambda2, C, CNorm, maxiter,
                          b_init, mu, tol, n, maxeig, J, K, P, lambda2weights_SPG=FALSE, Clusterinfo_SPG=NULL){

  ##### Functions for the selection of regularizaion parameters based on AICc. For the autoregressive coeffcients with Lasso ####

  #########
  # INPUT #
  #########
  # Y : vector of dimension NJK x 1 of dependent variables
  # X : matrix of dimension NJK x J�PK of independent variables
  # XX_prod : matrix of dimension J�PK x J�PK
  # XY_prod : matrix of dimension J�PK x 1
  # lambda1 : regularization parameter for the autoregressive coeffcients (lasso)
  # lambda2 : regularization parameter for the autoregressive coeffcients (fusion)
  # C : C matrix (built with 'Cmatrix.R')
  # CNorm : norm of C (built with 'Cmatrix.R')
  # maxiter : maxiter SPG
  # b_init : vector of dimension J�PK x 1 with the initial values of the autoregressive coefficients (if no default, take all zeros)
  # mu : smoothness parameter (if no default, take 1e-4)
  # tol : tolerance for SPG
  # n : time series length
  # maxeig : first eigenvalue of the matrix XX
  # J : number of time series in each class
  # K : number of classes
  # P : order of VAR
  # lambda2weights_SPG : logical. If TRUE, then do the weighting based on Clusterinfo. Default is FALSE.
  # Clusterinfo_SPG : K x K matrix of additional clustering information. K is the number of classes. Default is NULL.

  ##########
  # OUTPUT #
  ##########
  # AICc_values: AICc values

  FIT<-SPG(Y=Y, X=X, XX_prod=XX_prod,XY_prod=XY_prod, J=J, P=P, lambda1=lambda1, lambda2=lambda2, C=C, CNorm=CNorm,
           maxiter=maxiter, b_init=b_init, mu=mu, tol=tol,maxeig=maxeig, type="Lasso", lambda2weights_SPG=lambda2weights_SPG, Clusterinfo_SPG=Clusterinfo_SPG)

  Log_lik<-(1/(n*J*K))*t((Y-X%*%FIT$beta))%*%(Y-X%*%FIT$beta)
  df<-length(which(FIT$beta!=0))
  AIC_value<-(-2*Log_lik+2*df)
  correction<-2*df*(df+1)/(n*J*K-df-1)
  AICc_value<-AIC_value+correction
  return(AICc_value)
}

"AICcfunction_JGL" = function(Y, S_list, penalty, lambda1, lambda2, n, J, K, maxit, tol){

  ##### AICc for the inverse error covariance matrix with JGL ####

  FIT<-JGL(Y=Y, penalty=penalty, lambda1=lambda1, lambda2=lambda2, maxiter = maxit, tol = tol, return.whole.theta = T)

  #########
  # INPUT #
  #########
  # Y : list of data matrices for JGL
  # S_list : list of inverse error covariance matrices
  # penalty : penalty of JGL (choose "fused")
  # lambda1 : regularization parameter on the elements of inverse error covariance matrix (lasso)
  # lambda2 : regularization parameter on the difference between elements of inverse error covariance matrix (fusion)
  # n : time series length
  # J : number of time series per class
  # maxit : maximum iteration for JGL
  # tol : tolerance for JGL

  ##########
  # OUTPUT #
  ##########
  # AICc_values: AICc values


  more_arg=list(n=n, J=J, K=K)

  AICc_function<-function(A, B, J, K, n){
    nk<-n*J
    M<-length(which(B!=0))
    E<-(M-J)/2
    AIC<-(nk*(sum(diag(A%*%B)))-nk*log(det(B))+log(2)*E)
    correction<-2*E*(E+1)/(nk-E-1)
    AICc<-AIC+correction
  }

  AICc_value<-mapply(AICc_function, A=S_list, B=FIT$theta, MoreArgs = more_arg)
  AICc_values<-sum(AICc_value)
  AICc_values

  return(AICc_values)
}

"AICcfunction_RIDGE" = function(l_ridge_AICc, Y_AICc, X_AICc, J_AICc, P_AICc, n_AICc){

  #### AICc for the autoregressive coeffcients with RIDGE ####

  #########
  # INPUT #
  #########
  # l_ridge_AICc : grid of regularization paramters for the Ridge estimator
  # Y_AICc : vector of responses of length nJ
  # X_AICc : smatrix of inputs of dimension nJ x PJ^2
  # J_AICc : number of time series
  # P_AICc : order of VAR
  # n_AICc : time series length

  ##########
  # OUTPUT #
  ##########
  # AICc_values: AICc values

  # Ridge estimator
  Ridge_Function<-function(RX, RY, RJ, RP, Rl_ridge){
    B_r<- solve(t(RX)%*%RX + Rl_ridge*diag(RP*RJ^2)) %*% (t(RX)%*%RY)
  }

  # Estimate the autoregressive coefficients using SPG algorithm
  B_r<-Ridge_Function(RY=Y_AICc, RX=X_AICc, RJ=J_AICc, RP=P_AICc, Rl_ridge=l_ridge_AICc)

  # Compute the AICc
  Log_lik<-(1/n_AICc)*t((Y_AICc-X_AICc%*%B_r))%*%(Y_AICc-X_AICc%*%B_r)
  df<-length(which(B_r!=0))/J_AICc
  AIC_value<-(-2*Log_lik+2*df)
  correction<-2*df*(df+1)/(n_AICc-df-1)
  AICc_value<-AIC_value+correction
  return(AICc_value)

}

"AICcfunction_AdLasso" = function(Y, X, XX_prod, XY_prod, lambda1, lambda2, C, CNorm, maxiter,
                                  b_init, mu, tol, n,maxeig, J, K, P, weight, lambda2weights_SPG=NULL, Clusterinfo_SPG=NULL){

  #### AICc for for the autoregressive coeffcients with Adaptive Lasso ####

  #########
  # INPUT #
  #########
  # Y : vector of dimension NJK x 1 of dependent variables
  # X : matrix of dimension NJK x J�PK of independent variables
  # XX_prod : matrix of dimension J�PK x J�PK
  # XY_prod : matrix of dimension J�PK x 1
  # lambda1 : regularization parameter for the autoregressive coeffcients (lasso)
  # lambda2 : regularization parameter for the autoregressive coeffcients (fusion)
  # C : C matrix (built with 'Cmatrix.R')
  # CNorm : norm of C (built with 'Cmatrix.R')
  # maxiter : maxiter SPG
  # b_init : vector of dimension J�PK x 1 with the initial values of the autoregressive coefficients (if no default, take all zeros)
  # mu : smoothness parameter (if no default, take 1e-4)
  # tol : tolerance for SPG
  # n : time series length
  # maxeig : first eigenvalue of the matrix XX
  # J : number of time series in each class
  # K : number of classes
  # P : order of the VAR
  # weight: weight for the Adaptive Lasso.
  # lambda2weights_SPG : logical. If TRUE, then do the weighting based on Clusterinfo. Default is FALSE.
  # Clusterinfo_SPG : K x K matrix of additional clustering information. K is the number of classes. Default is NULL.

  ##########
  # OUTPUT #
  ##########
  # AICc_values: AICc values

  FIT<-SPG(Y=Y, X=X, XX_prod=XX_prod,XY_prod=XY_prod, J=J, P=P, lambda1=lambda1, lambda2=lambda2, C=C, CNorm=CNorm,
           maxiter=maxiter, b_init=b_init, mu=mu, tol=tol,maxeig=maxeig, type="AdLasso", weight=weight,
           lambda2weights_SPG=lambda2weights_SPG, Clusterinfo_SPG=Clusterinfo_SPG)

  Log_lik<-(1/(n*K*J))*t((Y-X%*%FIT$beta))%*%(Y-X%*%FIT$beta)
  df<-length(which(FIT$beta!=0))
  AIC_value<-(-2*Log_lik+2*df)
  correction<-2*df*(df+1)/(n*K*J-df-1)
  AICc_value<-AIC_value+correction
  return(AICc_value)
}

"AICcfunction_GrLasso" = function(Y, X, XX_prod, XY_prod, lambda1, lambda2, C, CNorm, maxiter,
                                  b_init, mu, tol, n,maxeig, J, K, P, weight, group,lambda2weights_SPG=NULL, Clusterinfo_SPG=NULL){

  #### AICc for for the autoregressive coeffcients with Group Lasso ####

  #########
  # INPUT #
  #########
  # Y : vector of dimension NJK x 1 of dependent variables
  # X : matrix of dimension NJK x J�PK of independent variables
  # XX_prod : matrix of dimension J�PK x J�PK
  # XY_prod : matrix of dimension J�PK x 1
  # lambda1 : regularization parameter for the autoregressive coeffcients (lasso)
  # lambda2 : regularization parameter for the autoregressive coeffcients (fusion)
  # C : C matrix (built with 'Cmatrix.R')
  # CNorm : norm of C (built with 'Cmatrix.R')
  # maxiter : maxiter SPG
  # b_init : vector of dimension J�PK x 1 with the initial values of the autoregressive coefficients (if no default, take all zeros)
  # mu : smoothness parameter (if no default, take 1e-4)
  # tol : tolerance for SPG
  # n : time series length
  # maxeig : first eigenvalue of the matrix XX
  # J : number of time series in each class
  # K : number of classes
  # P : order of the VAR
  # weight : Ridge estimator for the Group Lasso
  # group : grouping structure of the Group Lasso
  # lambda2weights_SPG : logical. If TRUE, then do the weighting based on Clusterinfo. Default is FALSE.
  # Clusterinfo_SPG : K x K matrix of additional clustering information. K is the number of classes. Default is NULL.

  ##########
  # OUTPUT #
  ##########
  # AICc_values: AICc values

  FIT<-SPG(Y=Y, X=X, XX_prod=XX_prod,XY_prod=XY_prod, J=J, P=P, lambda1=lambda1, lambda2=lambda2, C=C, CNorm=CNorm,
           maxiter=maxiter, b_init=b_init, mu=mu, tol=tol,maxeig=maxeig, type="GrLasso", weight=weight, group=group,
           lambda2weights_SPG=lambda2weights_SPG, Clusterinfo_SPG=Clusterinfo_SPG)

  Log_lik<-(1/(n*K*J))*t((Y-X%*%FIT$beta))%*%(Y-X%*%FIT$beta)
  df<-length(which(FIT$beta!=0))
  AIC_value<-(-2*Log_lik+2*df)
  correction<-2*df*(df+1)/(n*K*J-df-1)
  AICc_value<-AIC_value+correction
  return(AICc_value)
}


RIDGE<-function(Y, X, J, P, n, l_ridge_min=NULL, l_ridge_max=NULL, l_ridge_steps=NULL, l_ridge_opt=NULL, criterion){

  #### function to compute the Ridge estimator ####

  #########
  # INPUT #
  #########
  # Y : vector of responses of length nJ
  # X : matrix of inputs of dimension nJ x PJ^2
  # J : number of time series
  # P : order of VAR
  # n : time series length
  # l_ridge_min : minimum value of the regularization parameter on the autoregressive coeffcients in the RIDGE estimator. Default is NULL.
  # l_ridge_max : maximum value of the regularization parameter on the autoregressive coeffcients in the RIDGE estimator. Default is NULL.
  # l_ridge_steps : number of steps in the grid search for the regularization parameter on the autoregressive coeffcients in the RIDGE estimator. Default is NULL.
  # l_ridge_opt : optimal regularization parameter on the autoregressive coeffcients in the RIDGE estimator. Default is NULL.
  # criterion: "BIC" or "AICc", selection criterion used for the selction of the regularization paramaters.

  ##########
  # OUTPUT #
  ##########
  # B_ridge :  vector of estimated RIDGE coefficients of length PJ^2
  # l_ridge_opt : optimal regularization parameter


  # Opposite of null
  is.not.null <- function(x) ! is.null(x)

  if (is.not.null(l_ridge_min) & is.not.null(l_ridge_max) & is.not.null(l_ridge_steps) ){
    lr<-seq(from=l_ridge_min, to=l_ridge_max, length=l_ridge_steps)
    grid<-alply(lr, 1,.dims=TRUE)

    if (criterion=="BIC"){
      BIC_values_RIDGE<-lapply(X=grid, FUN=BICfunction_RIDGE, Y_BIC=Y, X_BIC=X, J_BIC=J, P_BIC=P,  n_BIC=n)
      l_ridge_OPT<-grid[[which.min(BIC_values_RIDGE)]]

    }
    if (criterion=="AICc"){
      AICc_values_RIDGE<-lapply(X=grid, FUN=AICcfunction_RIDGE, Y_AICc=Y, X_AICc=X, J_AICc=J, P_AICc=P,  n_AICc=n)
      l_ridge_OPT<-grid[[which.min(AICc_values_RIDGE)]]
    }
  }


  if (is.not.null(l_ridge_opt)){
    l_ridge_OPT<-l_ridge_opt
  }

  # Ridge estimator
  B_ridge<- solve(t(X)%*%X + l_ridge_OPT*diag(P*J^2)) %*% (t(X)%*%Y)

  RIDGE<-list(B_ridge=B_ridge, l_ridge_opt=l_ridge_OPT)

}

lambda2weights<-function(X,K,P,J){

  ##### Function to compute adaptive fusion regularization parameters, needed when incorporating additional clustering information ####

  #########
  # INPUT #
  #########
  # X : KxK dissimilarity matrix
  # K : number of classes (K>1)
  # P : order of VAR
  # J : number of time series in each class

  #########
  # OUTPUT #
  #########
  # lambda2_weight weights for lambda2


  DISTMAX<-matrix(apply(X,1,max),nrow=K,ncol=1)
  DISTscale<-matrix(NA,ncol=K,nrow=K) #Maximal distance set to one, all other values below one --> lower distance --> higher similarity
  for(i in 1:K){
    DISTscale[i,]<-X[i,]/DISTMAX[i,1]
  }

  lambda2_weight<-c()
  for(i.kout in 1:(K-1)){
    for(i.kin in 1:(K-i.kout)){
      newdist<-DISTscale[i.kin,i.kin+i.kout]
      lambda2_weight<-c(lambda2_weight,rep(newdist,J*J*P))
    }
  }
  return(lambda2_weight)
}

MultiClass_VAR_generalOmega<-function(Data=NULL, P,
                                      Y=NULL, X=NULL, X0=NULL, XX_prod=NULL, XY_prod=NULL, J=NULL, K=NULL, n=NULL,
                                      calculate.C=TRUE, C=NULL, CNorm=NULL, b_init=NULL, mu=1e-4,
                                      maxit.both=10, maxit.beta=10, maxit.omega=10, tol.both=10^-5, tol.beta=0.01, tol.omega=0.01,
                                      lambda1_min=0.1, lambda1_max=25, lambda1_steps=5,
                                      lambda2_min=0.1, lambda2_max=2, lambda2_steps=5,
                                      gamma1_min=0.1, gamma1_max=1, gamma1_steps=3,
                                      gamma2_min=0.1, gamma2_max=1, gamma2_steps=3,
                                      lambda1_OPT=NULL, lambda2_OPT=NULL, gamma1_OPT=NULL, gamma2_OPT=NULL,
                                      lambdaRIDGE_min=0.1, lambdaRIDGE_max=10, lambdaRIDGE_steps=5, lambdaRIDGE_OPT=NULL,
                                      criterion="BIC", type="AdLasso", gamma_ridge=NULL, group=NULL,
                                      lambda_weights=FALSE, Clusterinfo=NULL){

  ####  Function for the Multi-class VAR estimation allowing for cross-error correlations ####

  #########
  # INPUT #
  #########
  # Data: a N*J*K array, where:
  #  - N is the times series length
  #  - J is the number of time series
  #  - K is the number of classes
  # P: order of VAR
  # Y : vector of dimension NJK x 1 of responses (built with 'stack_fun.R'). Default is NULL.
  # X : matrix of dimension NJK x J�PK of predictors (built with 'stack_fun.R'). Default is NULL.
  # X0 : array of dimension n x JP x K of predictors (built with 'stack_fun.R'). Default is NULL.
  # XX_prod : cross-product matrix of dimension J�PK x J�PK. Default is NULL.
  # XY_prod : cross-product matrix of dimension J�PK x 1. Default is NULL.
  # J : number of time series in each class.
  # K : number of classes (K>1).
  # n : time series length.
  # calculate.C : logical for calculating the C matrix.
  # C : C matrix (built with 'Cmatrix.R').
  # CNorm : norm of C (built with 'Cmatrix.R').
  # b_init : vector of dimension J�PK x 1 with the initial values of the autoregressive coefficients.
  # mu : smoothness parameter. Default 1e-4.
  # maxit.both : maximum iterations of the Multiclass VAR
  # maxit.beta : maximum iterations  SPG algorithm
  # maxit.omega : maximum iterations JGL algorithm
  # tol.both : tolerance  Multiclass VAR
  # tol.beta : tolerance SPG algorithm
  # tol.omega : tolerance JGL algorithm
  # lambda1_min : minimum value of the regularization parameter on Beta (lasso penalty)
  # lambda1_max : maximum value of the regularization parameter on Beta (lasso penalty)
  # lambda1_steps : number of steps in the grid lambda1
  # lambda2_min : minimum value of the regularization parameter on Beta (fusion penalty)
  # lambda2_max : maximum value of the regularization parameter on Beta (fusion penalty)
  # lambda2_steps : number of steps in the grid lambda2
  # gamma1_min : minimum value of the regularization parameter on Omega (lasso penalty)
  # gamma1_max : maximum value of the regularization parameter on Omega (lasso penalty)
  # gamma1_steps : number of steps in the grid gamma1
  # gamma2_min : minimum value of the regularization parameter on Omega (fusion penalty)
  # gamma2_max : maximum value of the regularization parameter on Omega (fusion penalty)
  # gamma2_steps : number of steps in the grid gamma2
  # lambda1_OPT: optimal regularization parameter on B (lasso penalty). Default NULL.
  # lambda2_OPT: optimal regularization parameter on B (fusion penalty). Default NULL.
  # gamma1_OPT: optimal regularization parameter on Omega (lasso penalty). Default NULL.
  # gamma2_OPT: optimal regularization parameter on Omega (fusion penalty). Default NULL.
  # lambdaRIDGE_min : minimum value of the regularization parameter on B for the ridge estimator (only if type="AdLasso").
  # lambdaRIDGE_max : maximum value of the regularization parameter on B for the ridge estimator (only if type="AdLasso").
  # lambdaRIDGE_steps : number of steps in the grid lambda_ridge
  # lambdaRIDGE_OPT : array (1,1,K) containing the optimal regularization parameter on B for the ridge estimator
  # criterion: "BIC" or "AICc", selection criterion used for the selection of the regularization paramaters. Default is BIC.
  # type: "AdLasso"
  # gamma_ridge: exponent for the weight for the Adpative Lasso. Default is NULL.
  # group : grouping structure of the Group Lasso. Vector of length KPJ^2. Default is NULL.
  # lambda2_weights : logical. If TRUE, then do the weighting of fusion parameter based on X_weight. Default is FALSE.
  # Clusterinfo : K x K matrix of additional clustering information. K is the number of classes. Default is NULL.

  ##########
  # OUTPUT #
  ##########
  # Beta_new : vector of dimension J�PK x 1 with the autoregressive coefficients. Given B^{(k)}_{p,ij}, with i referring to the dependent variable and j to the independent variable, the structure of Beta_new=[B^{(1)}_{1,11}, B^{(1)}_{1,12}, ..., B^{(1)}_{1,1J}, B^{(1)}_{2,11}, B^{(1)}_{2,12}, ..., B^{(1)}_{P,1J}, B^{(1)}_{1,21}, ...,B^{(1)}_{P,JJ}, B^{(2)}_{1,11}, B^{(2)}_{1,12}, ..., B^{(K)}_{P,JJ}]'
  # Beta_array : array (P*J,J,K) of estimated coeffficients. [,,k] contains the B for class k. In [,,k], the first [1:J,1:J] contains the B at lag p=1 in this order c(row1: B_11,.., B_1J; row2: B_21, ..., B_2J; ...; rowJ: B_J1, ..., B_JJ). In [,,k], the obs [1:J, (J(P-1)+1):JP] contain the at lag p=P in the order explained before.
  # CFIT: list containing the C matrix and the CNorm built with the function "buildCmatrix"
  # J : number of time series in each class.
  # K : number of classes.
  # P : VAR order.
  # n : time series length.
  # Y_array : array (n,J,K) of original responses. [,,k] contains the responses for class k. In [,,k], column j contains the responses of variable j from 1:n.
  # X_array : array (n,P*J,K) of original inputs. [,,k] contains the inputs for class k. In [,,k], the first [1:n,1:J] contains the inputs at lag p=1: column j, contains the inputs from (P-p+1):(T-p). In [,,k], the [1:n, (J(P-1)+1):JP] contains the J inputs at lag P.
  # E_array : array (n,J,K) of errors. [,,k] contains the errors for class k. In [,,k], column j contains the error of variable j from 1:n.
  # omegahat_matrix : KJxKJ matrix of estimated error correlations
  # lambda1_opt : regularization parameter on the autoregressive coeffcients
  # lambda2_opt : regularization parameter on the difference between autoregressive coeffcients
  # gamma1_opt : regularization parameter on the elements of inverse error covariance matrix
  # gamma2_opt : regularization parameter on the difference between elements of inverse error covariance matrix
  # lambda_RIDGE : regularization paramater for RIDGE


  # Additional functions #
  # Spectral Decomposition
  spectral_dec<-function(Omega){
    decomp<-eigen(Omega, symmetric = TRUE)
    CC <- decomp$vectors
    Lambda <- diag(sqrt(decomp$values))
    P <- CC %*% Lambda %*% t(CC)
  }

  # Opposite of is.null
  is.not.null <- function(x) ! is.null(x)

  # Opposite of is.array
  is.not.array <- function(x) ! is.array(x)

  # Log-likelihood of the Multi-class VAR
  LogLik_function <- function (X_list, Y_stack_list,  Beta_new_list, omega_tilde_list, omega_new_list, n, J){
    AA <- (Y_stack_list - X_list %*% Beta_new_list)
    ll<- ((t(AA) %*% omega_tilde_list %*% AA) - (n*J*log(det(omega_new_list))))
  }
  # Kronecker product
  transpose_function<-function(U,J){
    kronecker(diag(1,J),t(U)%*%U)
  }

  BIC_ADMMgeneral<-function(e_list, lambda1, lambda2, n, J, maxit, tol,K){

    # Fit the JGL
    FIT<-ADMM.GeneralOmega(e_list=e_list,rho=1,lambda1=lambda1,lambda2=lambda2,eps=2*tol,max.iter=maxit)

    BIC_Omegafull<-function(e_list,J,K,omegahat,n){
      S<-matrix(NA,ncol=J*K,nrow=J*K)
      for(i.row in 1:(K-1)){
        for(j.col in (i.row+1):K){
          S[((i.row-1)*J+1):(i.row*J),((j.col-1)*J+1):(j.col*J)]<-cov(e_list[[i.row]],e_list[[j.col]])
        }
      }
      for(i.row in 1:K){
        S[((i.row-1)*J+1):(i.row*J),((i.row-1)*J+1):(i.row*J)]<-cov(e_list[[i.row]],e_list[[i.row]])
      }
      S[lower.tri(S)]<-t(S)[lower.tri(S)] # dimension JK x JK

      omegadf<-length(which(diag(omegahat)!=0))
      omegacross<-omegahat
      diag(omegacross)<-0
      omegadf<-omegadf+length(which(omegacross!=0))/2

      BIC<- sum(diag(S%*%omegahat)) - log(det(omegahat))+(log(n*K*(K+1)/2)/(n*K*(K+1)/2))*omegadf
      return(BIC)
    }

    BICvalue<-BIC_Omegafull(e_list=e_list,omegahat=FIT$Omegahat,n=n,J=J,K=K)
    return(BICvalue)
  }

  AICc_ADMMgeneral<-function(e_list, lambda1, lambda2, n, J, maxit, tol,K){

    # Fit the JGL
    FIT<-ADMM.GeneralOmega(e_list=e_list,rho=1,lambda1=lambda1,lambda2=lambda2,eps=2*tol,max.iter=maxit)

    AICc_Omegafull<-function(e_list,J,K,omegahat,n){
      S<-matrix(NA,ncol=J*K,nrow=J*K)
      for(i.row in 1:(K-1)){
        for(j.col in (i.row+1):K){
          S[((i.row-1)*J+1):(i.row*J),((j.col-1)*J+1):(j.col*J)]<-cov(e_list[[i.row]],e_list[[j.col]])
        }
      }
      for(i.row in 1:K){
        S[((i.row-1)*J+1):(i.row*J),((i.row-1)*J+1):(i.row*J)]<-cov(e_list[[i.row]],e_list[[i.row]])
      }
      S[lower.tri(S)]<-t(S)[lower.tri(S)] # dimension JK x JK

      omegadf<-0
      omegacross<-omegahat
      diag(omegacross)<-0
      omegadf<-omegadf+length(which(omegacross!=0))/2

      AIC<- sum(diag(S%*%omegahat)) - log(det(omegahat))+(2/(n*K*(K+1)/2))*omegadf
      correction<- 2*omegadf*(omegadf+1)/(n*K*(K+1)/2-omegadf-1)
      AICc<-AIC + correction
      return(AICc)
    }

    AICcvalue<-AICc_Omegafull(e_list=e_list,omegahat=FIT$Omegahat,n=n,J=J,K=K)
    return(AICcvalue)
  }

  LogLik_objective<-function(e_list,omegahat_matrix,J,K){
    S<-matrix(NA,ncol=J*K,nrow=J*K)
    for(i.row in 1:(K-1)){
      for(j.col in (i.row+1):K){
        S[((i.row-1)*J+1):(i.row*J),((j.col-1)*J+1):(j.col*J)]<-cov(e_list[[i.row]],e_list[[j.col]])
      }
    }
    for(i.row in 1:K){
      S[((i.row-1)*J+1):(i.row*J),((i.row-1)*J+1):(i.row*J)]<-cov(e_list[[i.row]],e_list[[i.row]])
    }
    S[lower.tri(S)]<-t(S)[lower.tri(S)] # dimension JK x JK

    nll<- sum(diag(S%*%omegahat_matrix))-log(det(omegahat_matrix))
    return(nll)
  }


  # Preliminary steps #
  # Check Data
  if (is.null(Data) & (is.null(Y) | is.null(X) | is.null(X0) | is.null(J) | is.null(K) | is.null(n)) ){
    stop("Either give as input Data or give as input X, Y, X0, J, K, n.")
  }

  # Generate the inputs
  if (is.not.null(Data)){
    INPUTS<-Inputs_MultiClass_VAR(Data, P=P,calculate.C=calculate.C)
    Y <- INPUTS$Y
    X <- INPUTS$X
    X0 <- INPUTS$X0
    XX_prod <- INPUTS$XX
    XY_prod <- INPUTS$XY

    # Dimensions
    Ndata<-dim(Data)[1]
    J<-dim(Data)[2]
    K<-dim(Data)[3]
    n<-Ndata-P
  }

  if ((calculate.C==F & (is.null(C) | is.null(CNorm)))  ){
    stop("Either set calculate.C=TRUE or give the input C and CNorm.")
  }

  if (calculate.C==TRUE){
    C <- INPUTS$C
    CNorm <- INPUTS$CNorm
    CFIT<-list("C"=INPUTS$C,"CNorm"=INPUTS$CNorm)
  }

  if (is.not.null(C) & is.not.null(CNorm)){
    CFIT<-list("C"=C,"CNorm"=CNorm)
  }

  # Additional inputs
  if(is.null(b_init)){
    b_init<-array(0,c(ncol(X),1))
  }

  if(is.null(XX_prod)){
    X0_list<-alply(X0,3,.dims=TRUE)
    XtX<-lapply(X0_list,transpose_function,J=J)
    XX_prod<-bdiag(XtX)
  }

  if(is.null(XY_prod)){
    XY_prod<-t(X)%*%Y
  }


  if (is.not.null(lambdaRIDGE_OPT) & is.not.array(lambdaRIDGE_OPT)){
    stop("lambdaRIDGE_OPT must be an array c(1,1,K).")
  }

  if (lambda_weights==T & is.null(Clusterinfo)){
    stop("Insert a value for Clusterinfo.")
  }


  # Check the optimal regularization paramters
  if ((is.not.null(lambda1_OPT) | is.not.null(lambda2_OPT) | is.not.null(gamma1_OPT) | is.not.null(gamma2_OPT))
      & (is.not.null(lambda1_steps) | is.not.null(lambda2_steps) | is.not.null(gamma1_steps) | is.not.null(gamma2_steps)) ){
    stop("Either set the optimal values of lambda&gamma or set the values for the grid search.")
  }
  if ((is.null(lambda1_OPT) & is.null(lambda2_OPT) & is.null(gamma1_OPT) & is.null(gamma2_OPT))
      & (is.null(lambda1_steps) & is.null(lambda2_steps) & is.null(gamma1_steps) & is.null(gamma2_steps))
      & (is.null(lambda1_min) & is.null(lambda1_max) & is.null(gamma1_min) & is.null(gamma1_max))
      & (is.null(lambda2_min) & is.null(lambda2_max) & is.null(gamma2_min) & is.null(gamma2_max)) ){
    stop("Either set the optimal values of lambda&gamma or set the values for the grid search.")
  }

  if ( is.not.null(lambdaRIDGE_OPT) & (is.not.null(lambdaRIDGE_min) | is.not.null(lambdaRIDGE_max) & is.not.null(lambdaRIDGE_steps) )){
    stop("Either set the optimal values of lambdaRIDGE or set the values for the grid search.")
  }

  if (is.null(lambda1_OPT) & is.not.null(lambda2_OPT)){
    stop("Set the optimal value of lambda1.")
  }
  if (is.not.null(lambda1_OPT) & is.null(lambda2_OPT)){
    stop("Set the optimal value of lambda2.")
  }
  if (is.null(gamma1_OPT) & is.not.null(gamma2_OPT)){
    stop("Set the optimal value of gamma1.")
  }
  if (is.not.null(gamma1_OPT) & is.null(gamma2_OPT)){
    stop("Set the optimal value of gamma2.")
  }


  # Check the grid of the regularization paramters
  if ((is.null(lambda1_min) | is.null(lambda1_max) | is.null(lambda1_steps) ) & (is.not.null(lambda1_min) | is.not.null(lambda1_max) | is.not.null(lambda1_steps))){
    stop("Set all the values of the lambda1 grid.")
  }
  if ((is.null(lambda2_min) | is.null(lambda2_max) | is.null(lambda2_steps) ) & (is.not.null(lambda2_min) | is.not.null(lambda2_max) | is.not.null(lambda2_steps))){
    stop("Set all the values of the lambda2 grid.")
  }
  if ((is.null(gamma1_min) | is.null(gamma1_max) | is.null(gamma1_steps) ) & (is.not.null(gamma1_min) | is.not.null(gamma1_max) | is.not.null(gamma1_steps))){
    stop("Set all the values of the gamma1 grid.")
  }
  if ((is.null(gamma2_min) | is.null(gamma2_max) | is.null(gamma2_steps) ) & (is.not.null(gamma2_min) | is.not.null(gamma2_max) | is.not.null(gamma2_steps))){
    stop("Set all the values of the gamma2 grid.")
  }



  # Set the grid for the regularization parameters
  if (is.not.null(lambda1_min) & is.not.null(lambda1_max) & is.not.null(lambda1_steps) & is.not.null(lambda2_min) & is.not.null(lambda2_max) & is.not.null(lambda2_steps) ){
    lambda1<-seq(from=lambda1_min, to=lambda1_max, length=lambda1_steps)
    lambda2<-seq(from=lambda2_min, to=lambda2_max, length=lambda2_steps)
    l1<-list(lambda1_l=sort(rep(lambda1,length(lambda2))), lambda2_l=rep(lambda2,length(lambda1)))
  }

  if (is.not.null(gamma1_min) & is.not.null(gamma1_max) & is.not.null(gamma1_steps) & is.not.null(gamma2_min) & is.not.null(gamma2_max) & is.not.null(gamma2_steps)){
    gamma1<-seq(from=gamma1_min,to=gamma1_max,length=gamma1_steps)
    gamma2<-seq(from=gamma2_min,to=gamma2_max,length=gamma2_steps)
    g1<-list(gamma1_l=sort(rep(gamma1,length(gamma2))), gamma2_l=rep(gamma2,length(gamma1)))
  }



  #### START CODE #####

  # Data preliminaries
  Ynew <- Y; YData<-Y;rm(Y)
  Xnew <- X; XData<-X;rm(X)
  XX_prodnew<-XX_prod;rm(XX_prod)
  XY_prodnew<-XY_prod;rm(XY_prod)
  maxeig <- eigs(XX_prodnew, k=1, opts = list(retvec = FALSE))

  Obj=rep(1,maxit.both)
  Obj_Conv=rep(1,maxit.both)
  iter<-2


  while (Obj_Conv[iter]>tol.both & iter<maxit.both) {

    #########
    ## SPG ##
    #########

    # Estimation of the autoregressive coeffcients using SPG algorithm

    ##################
    # Adaptive Lasso #
    ##################

    # RIDGE estimator for each class --> weight for the Adaptive Lasso
    Y_arr <- array(Ynew, c((n*J),1,K))    # Array nJx1xK: In [,,k], you have the first [1:n] obs of series j=1, [n+1:2n] obs of series j=2...
    Y_list<-alply(Y_arr,3,.dims=TRUE);rm(Y_arr)     # List of K matrices of dimension nJx1
    X0_list<-alply(X0,3,.dims=TRUE)       # List of nxJ matrices
    X_list<-lapply(X0_list,function(U){kronecker(diag(J),U)});rm(X0_list) # Knonecker diag(J) and X0_list

    # Ridge estimator with the OPTIMAL regularization parameters
    if (is.not.null(lambdaRIDGE_OPT)){
      lambdaRIDGE_OPT_list<-alply(lambdaRIDGE_OPT,3,.dims=TRUE)       # List of Optimal RIDGE regularization parameters
      mr<-list(J=J, P=P, n=n,  criterion=criterion)
      LIST_RIDGE<-mapply(RIDGE, Y=Y_list, X=X_list, l_ridge_opt=lambdaRIDGE_OPT_list, MoreArgs = mr)
    }
    # Ridge estimator with the regularization parameters selected via BIC or AICc
    if (is.not.null(lambdaRIDGE_min) & is.not.null(lambdaRIDGE_max) & is.not.null(lambdaRIDGE_steps)){
      mr<-list(J=J, P=P, n=n, l_ridge_min=lambdaRIDGE_min, l_ridge_max=lambdaRIDGE_max, l_ridge_steps=lambdaRIDGE_steps,
               criterion=criterion)
      LIST_RIDGE<-mapply(RIDGE, Y=Y_list, X=X_list, MoreArgs = mr)
    }
    # Ridge output
    beta_RIDGE_arr<-array(NA, c((P*J^2), 1, K))
    lambda_RIDGE<-array(NA, c(1,1,K))
    for (ik in 1:K){
      beta_RIDGE_arr[,,ik]<-LIST_RIDGE[[2*ik-1]]
      lambda_RIDGE[,,ik]<-LIST_RIDGE[[2*ik]]
    }
    beta_RIDGE_stacked<-matrix(stack(data.frame(beta_RIDGE_arr))[,1], nrow = K*P*J^2, ncol=1);rm(beta_RIDGE_arr)
    if (is.null(gamma_ridge)){
      gamma_ridge=1
    }
    weight_RIDGE<-(abs(beta_RIDGE_stacked))^gamma_ridge # Stacked coefficients. The first PJ^2 elements belong to k=1, ...

    if (is.null(lambda1_OPT)  &  is.null(lambda2_OPT)){
      # Select via BIC the regularization paramaters for the Adaptive Lasso
      if (criterion == "BIC"){
        morearg=list(Y=Ynew, X=Xnew, XX_prod=XX_prodnew,XY_prod=XY_prodnew,C=C, CNorm=CNorm, J=J, P=P,
                     maxiter=maxit.beta, b_init=b_init, mu=mu, tol=tol.beta, n=n,maxeig=maxeig, weight=weight_RIDGE,
                     lambda2weights_SPG=lambda_weights, Clusterinfo_SPG=Clusterinfo,K=K)
        if (iter==2){
          BIC_values<-mapply(BICfunction_AdLasso,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
        } else{
          BIC_values<-mapply(BICfunction_AdLasso,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)

          if(class(BIC_values)=="list"){
            BIC_values<-c(do.call("cbind",BIC_values))
            BIC_values<-c(matrix(BIC_values[[1]]))
          }

        }

        which.min(BIC_values)
        lambda1_opt<-l1$lambda1_l[which.min(BIC_values)]
        lambda2_opt<-l1$lambda2_l[which.min(BIC_values)]
      }

      # Select via AICc the regularization paramaters for the Adaptive Lasso
      if (criterion == "AICc"){
        morearg=list(Y=Ynew, X=Xnew, XX_prod=XX_prodnew,XY_prod=XY_prodnew,C=C, CNorm=CNorm,
                     maxiter=maxit.beta, b_init=b_init, mu=mu, tol=tol.beta, n=n,maxeig=maxeig, J=J, K=K, P=P,
                     weight=weight_RIDGE, lambda2weights_SPG=lambda_weights, Clusterinfo_SPG=Clusterinfo)
        if (iter==2){
          AICc_values<-mapply(AICcfunction_AdLasso,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
        } else{
          AICc_values<-mapply(AICcfunction_AdLasso,l1$lambda1_l, l1$lambda2_l, MoreArgs = morearg)
          if(class(AICc_values)=="list"){
            AICc_values<-c(do.call("cbind",AICc_values))
            AICc_values<-c(matrix(AICc_values[[1]]))
          }

        }

        which.min(AICc_values)
        lambda1_opt<-l1$lambda1_l[which.min(AICc_values)]
        lambda2_opt<-l1$lambda2_l[which.min(AICc_values)]
      }
    }

    if (is.not.null(lambda1_OPT) &  is.not.null(lambda2_OPT)){
      lambda1_opt<-lambda1_OPT
      lambda2_opt<-lambda2_OPT
    }

    # SPG estimation for the Adaptive Lasso
    spg<-SPG(Y=Ynew, X=Xnew, XX_prod=XX_prodnew,XY_prod=XY_prodnew, J=J, P=P,
             lambda1=lambda1_opt, lambda2=lambda2_opt, C=C, CNorm=CNorm,
             maxiter=maxit.beta, b_init=b_init, mu=mu, tol=tol.beta,maxeig=maxeig, type="AdLasso", weight=weight_RIDGE,
             lambda2weights_SPG=lambda_weights, Clusterinfo_SPG=Clusterinfo)





    # Coefficients from the SPG
    Beta_new<-spg$beta                          # autoregressive coeffcients
    Beta_new_arr<-array(Beta_new, c(P*(J^2),1,K))
    Beta_new_list<-alply(Beta_new_arr,3,.dims=TRUE)
    iter_spg<-spg$it                            # iteration
    rm(spg);rm(Ynew);rm(Xnew)

    # Residuals
    e = (YData-XData%*%Beta_new)                        # stacked from the SPG
    e_array<-array(e, c(n,J,K))
    e_array_list<-alply(e_array,3,.dims=TRUE)   # create a list from the array
    S<-array(apply(e_array,3,cov),c(J,J,K))     # VAR-COV matrices
    S_list<-alply(S,3,.dims=TRUE)               # create a list from the array


    #########
    ## JGL ##
    #########

    # Estimation of the inverse error covariance matrix using JGL package

    if (is.null(gamma1_OPT)  &  is.null(gamma2_OPT)){
      # Use BIC as criterion
      if (criterion == "BIC"){
        morearg=list(e_list=e_array_list, n=n, J=J, K=K, maxit = maxit.omega, tol = tol.omega)
        BIC_values_JGL<-mapply(BIC_ADMMgeneral,g1$gamma1_l, g1$gamma2_l, MoreArgs = morearg)
        gamma1_opt<-g1$gamma1_l[which.min(BIC_values_JGL)]
        gamma2_opt<-g1$gamma2_l[which.min(BIC_values_JGL)]
      }

      # Use the AICc criterion
      if (criterion == "AICc"){
        morearg=list(e_list=e_array_list, n=n, J=J, K=K, maxit = maxit.omega, tol = tol.omega)
        AICc_values_JGL<-mapply(AICc_ADMMgeneral,g1$gamma1_l, g1$gamma2_l, MoreArgs = morearg)
        gamma1_opt<-g1$gamma1_l[which.min(AICc_values_JGL)]
        gamma2_opt<-g1$gamma2_l[which.min(AICc_values_JGL)]
      }
    }

    if (is.not.null(gamma1_OPT) &  is.not.null(gamma2_OPT)){
      gamma1_opt<-gamma1_OPT
      gamma2_opt<-gamma2_OPT
    }

    # General Omega estimation
    jgl<-ADMM.GeneralOmega(e_list=e_array_list,rho=1,lambda1=gamma1_opt,lambda2=gamma2_opt,eps=2*tol.omega,max.iter=maxit.omega)
    omegahat_matrix<-jgl$Omegahat;rm(jgl)

    # Spectral Decomposition on Omega using Th.4.2.12 Horn-Johnson (1991)
    eigen_omega<-eigen(omegahat_matrix, symmetric = TRUE)
    e_omega<-eigen_omega$values
    e_omega_tilde<-Diagonal(length(rep(e_omega, each=n)),rep(e_omega, each=n)) # eigenvalues of the kronecker product in matrix form
    v_omega<-eigen_omega$vectors
    rm(e_omega);rm(eigen_omega)
    v_identity<-Matrix(eigen(Diagonal(n))$vector)# eigenvectors of I_{n}
    v_omega_tilde<-kronecker(v_omega,v_identity) # eigenvectors of the kronecker product

    # Define the new matrices
    Y_new<-matrix(v_omega_tilde%*%sqrt(e_omega_tilde)%*%t(v_omega_tilde)%*%YData,ncol=1)
    X_new<-matrix(v_omega_tilde%*%sqrt(e_omega_tilde)%*%t(v_omega_tilde)%*%XData,ncol=ncol(XData))
    rm(v_omega_tilde);rm(v_omega);rm(v_identity)

    X_new_arr<-array(NA, c(n*J, P*J^2, K))
    X_new_arr[,,1]<-X_new[1:(J*n),1:(P*J^2)]
    for (ik in 2:K){
      X_new_arr[,,ik]<-X_new[((ik-1)*J*n+1):(ik*n*J), ((ik-1)*P*J^2+1):(ik*P*J^2)]
    }
    X_new_list<-alply(X_new_arr, 3,.dims=TRUE)


    #    Check Convergence   #

    # Objective function
    Obj[iter] <- LogLik_objective(e_list=e_array_list,omegahat_matrix=omegahat_matrix,J=J,K=K)
    Obj_Conv[iter+1]<-(abs(Obj[iter]-Obj[iter-1])/Obj[iter-1])


    # Create the new INPUTS
    iter=iter+1
    Ynew <- Y_new
    Xnew <- X_new
    rm(Y_new);rm(X_new)

    XtX<-lapply(X_new_list,function(U){t(U)%*%U})
    XX_prodnew<-bdiag(XtX)
    XY_prodnew<-t(Xnew)%*%Ynew
    maxeig <- eigs(XX_prodnew, k=1, opts = list(retvec = FALSE))
    rm(XtX)

  } # end while loop


  ##########
  # OUTPUT #
  ##########

  Beta_array<-aperm(array(Beta_new, c(P*J,J,K)), c(2,1,3))  # Array (P*J,J,K) of estimated coeffficients. [,,k] contains the B for class k. In [,,k], the first [1:J,1:J] contains the B at lag p=1 in this order c(row1: B_11,.., B_1J; row2: B_21, ..., B_2J; ...; rowJ: B_J1, ..., B_JJ). In [,,k], the obs [1:J, (J(P-1)+1):JP] contain the at lag p=P in the order explained before.
  Y_array<-array(YData, c(n,J,K))                           # Array (n,J,K) of original responses. [,,k] contains the responses for class k. In [,,k], column j contains the responses of variable j from 1:n.
  X_array<-X0                                               # Array (n,P*J,K) of original inputs. [,,k] contains the inputs for class k. In [,,k], the first [1:n,1:J] contains the inputs at lag p=1: column j, contains the inputs from (P-p+1):(T-p). In [,,k], the [1:n, (J(P-1)+1):JP] contains the J inputs at lag P.
  e_arr = (YData-XData%*%Beta_new)
  E_array<-array(e_arr, c(n,J,K))                           # Array (n,J,K) of errors. [,,k] contains the errors for class k. In [,,k], column j contains the error of variable j from 1:n.


  MultiClass_VAR_generalOmega <- list(Beta_new=Beta_new, Beta_array=Beta_array, CFIT=CFIT,
                                      J=J, K=K, P=P, n=n,
                                      Y_array=Y_array, X_array=X_array, E_array=E_array,
                                      omegahat_matrix=omegahat_matrix,
                                      lambda1_opt=lambda1_opt, lambda2_opt=lambda2_opt,
                                      gamma1_opt=gamma1_opt, gamma2_opt=gamma2_opt,
                                      lambda_RIDGE=lambda_RIDGE)
}

ADMM.GeneralOmega<-function(e_list,rho=1,lambda1=0.1,lambda2=0.1,eps=2*10^-2,max.iter=50){
  #### Alternating Direction Method of Multipliers (ADMM) for the MultiClass General Omega ####

  #########
  # INPUT #
  #########
  # e_LIST: list of residuals of the multi-class VAR model
  # rho : step size parameter. Default is 1.
  # lambda1 : lasso penalty
  # lambda2 : fusion penalty
  # eps : convergence tolerance level
  # max.iter : maximum number of iterations

  ##########
  # OUTPUT #
  ##########
  # Omegahat : estimated Omega


  # Preliminaries
  J<-dim(e_list[[1]])[2]
  K<-length(e_list)

  # Center data matrices
  for (k in 1:K) {
    for (j in 1:J) {
      e_list[[k]][, j] = e_list[[k]][, j] - mean(e_list[[k]][, j])
    }
  }

  # Sample covariance matrix
  S<-matrix(NA,ncol=J*K,nrow=J*K)
  for(i.row in 1:(K-1)){
    for(j.col in (i.row+1):K){
      S[((i.row-1)*J+1):(i.row*J),((j.col-1)*J+1):(j.col*J)]<-cov(e_list[[i.row]],e_list[[j.col]])
    }
  }
  for(i.row in 1:K){
    S[((i.row-1)*J+1):(i.row*J),((i.row-1)*J+1):(i.row*J)]<-cov(e_list[[i.row]],e_list[[i.row]])
  }
  S[lower.tri(S)]<-t(S)[lower.tri(S)] # dimension JK x JK


  # # INIT
  Zold<-Uold<-Omegaold<-matrix(0,ncol=ncol(S),nrow=nrow(S))
  diag(Omegaold)<-1
  diff<-eps*10
  it<-1
  while( (diff>eps) & (it<max.iter)  ){

    # STEP 1: COMPUTE SINGULAR VALUE DECOMPOSITION
    svd.decomp<-svd(S-rho*Zold+rho*Uold)
    Djj<-(1/(2*rho))*(-svd.decomp$d + (svd.decomp$d^2 + 4*rho)^(1/2))
    Omeganew<-svd.decomp$v%*%diag(Djj)%*%t(svd.decomp$v)


    # STEP 2.1 LASSO on OFF-DIAGONAL ELEMENTS
    A.off<-Omeganew+Uold
    Znew<-soft(A.off,lam=lambda1/rho,penalize.diagonal=F)

    # STEP 2.2: FUSED LASSO ON MAIN DIAGONAL
    Omegalist=list()
    Alist = list()
    Ulist=list()
    for(k in 1:K){
      Omegalist[[k]]<-Omeganew[((k-1)*J+1):(k*J),((k-1)*J+1):(k*J)]
      Ulist[[k]]<-Uold[((k-1)*J+1):(k*J),((k-1)*J+1):(k*J)]
      Alist[[k]]<-Omegalist[[k]] + Ulist[[k]]
    }
    Z = flsa.general(Alist, L=rho, lam1=lambda1, lam2=lambda2, penalize.diagonal = TRUE)

    for(k in 1:K){
      Znew[((k-1)*J+1):(k*J),((k-1)*J+1):(k*J)]<-Z[[k]]
    }


    # STEP 3: update dual variable U
    Unew <-Uold + Omeganew - Znew

    # CHECK CONVERGENCE
    diff<-max(abs(Omegaold-Omeganew))
    # Update
    Zold<-Znew
    Uold<-Unew
    Omegaold<-Omeganew
    it<-it+1
  }

  out<-list("Omegahat"=Znew)

}



flsa.general<-function (A, L, lam1, lam2, penalize.diagonal)
{
  #### JGL PACKAGE: Function to compute the Fused Lasso on the elements on the main diagonal ####
  trueA = A
  if (is.matrix(A[[1]])) {
    p = dim(A[[1]])[1]
  }
  if (is.vector(A[[1]])) {
    p = length(A[[1]])
  }
  K = length(A)
  X = list()
  for (k in 1:K) {
    X[[k]] = A[[1]] * NA
  }
  if (is.matrix(A[[1]])) {
    fusions = array(FALSE, dim = c(K, K, p, p))
  }
  if (is.vector(A[[1]])) {
    fusions = array(FALSE, dim = c(K, K, p, 1))
  }
  newc = list()
  for (k in 1:K) {
    others = setdiff(1:K, k)
    others.smaller.k = 1:(k - 1)
    newc[[k]] = A[[1]] * 0
    for (o in others) {
      newc[[k]] = newc[[k]] + (A[[o]] - A[[k]] < -1e-04) -
        (A[[o]] - A[[k]] > 1e-04)
    }
  }
  for (iter in 1:(K - 1)) {
    ordermats = list()
    for (k in 1:K) {
      others = setdiff(1:K, k)
      others.smaller.k = 1:(k - 1)
      ordermats[[k]] = A[[1]] * 0
      for (o in others) {
        ordermats[[k]] = ordermats[[k]] + (A[[k]] - A[[o]] >
                                             1e-04)
      }
      if (k > 1) {
        for (o in others.smaller.k) {
          ordermats[[k]] = ordermats[[k]] + (abs(A[[o]] -
                                                   A[[k]]) < 1e-04)
        }
      }
      ordermats[[k]] = ordermats[[k]] + 1
    }
    betas.g = list()
    for (k in 1:K) {
      betas.g[[k]] = A[[k]] - lam2/L * newc[[k]]
    }
    new.ordermats = list()
    for (k in 1:K) {
      others = setdiff(1:K, k)
      others.smaller.k = 1:(k - 1)
      new.ordermats[[k]] = A[[1]] * 0
      for (o in others) {
        new.ordermats[[k]] = new.ordermats[[k]] + (betas.g[[k]] -
                                                     betas.g[[o]] > 1e-04)
      }
      if (k > 1) {
        for (o in others.smaller.k) {
          new.ordermats[[k]] = new.ordermats[[k]] + (abs(betas.g[[o]] -
                                                           betas.g[[k]]) < 1e-04)
        }
      }
      new.ordermats[[k]] = new.ordermats[[k]] + 1
    }
    for (k in 1:K) {
      for (kp in 1:K) {
        fusions[k, kp, , ] = fusions[k, kp, , ] + ((ordermats[[k]] -
                                                      1 == ordermats[[kp]]) & (new.ordermats[[k]] <
                                                                                 new.ordermats[[kp]])) + ((ordermats[[k]] +
                                                                                                             1 == ordermats[[kp]]) & (new.ordermats[[k]] >
                                                                                                                                        new.ordermats[[kp]])) + (abs(A[[k]] - A[[kp]]) <
                                                                                                                                                                   1e-04)
        fusions = (fusions > 0) * 1
      }
    }
    for (k in 1:K) {
      for (kp in 1:K) {
        others = setdiff(1:K, c(k, kp))
        for (o in others) {
          bothfused = fusions[k, o, , ] & fusions[kp,
                                                  o, , ]
          fusions[k, kp, , ] = fusions[k, kp, , ] | bothfused
        }
      }
    }
    for (k in 1:K) {
      others = setdiff(1:K, k)
      fusemean = trueA[[k]]
      denom = A[[1]] * 0 + 1
      for (o in others) {
        fusemean = fusemean + fusions[k, o, , ] * trueA[[o]]
        denom = denom + fusions[k, o, , ]
      }
      A[[k]] = fusemean/denom
    }
    newc = list()
    for (k in 1:K) {
      others = setdiff(1:K, k)
      others.smaller.k = 1:(k - 1)
      newc[[k]] = A[[1]] * 0
      for (o in others) {
        newc[[k]] = newc[[k]] + (A[[o]] - A[[k]] < -1e-04) -
          (A[[o]] - A[[k]] > 1e-04)
      }
    }
  }
  for (k in 1:K) {
    betas.g[[k]] = A[[k]] - lam2/L * newc[[k]]
  }
  for (k in 1:K) {
    X[[k]] = soft(betas.g[[k]], lam = lam1/L, penalize.diagonal = penalize.diagonal)
  }
  return(X)
}

soft<-function (a, lam, penalize.diagonal)
{ # Function for soft thresholding on the off-diagonal
  out <- sign(a) * pmax(0, abs(a) - lam)
  if (!penalize.diagonal)
    diag(out) <- diag(a)
  return(out)
}
