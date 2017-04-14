#### PAIRED BOOTSTRAP in linear regression
#### 2017 - 4-12:  code works on data analysis examples and a multivariate data import
####  --retrieval of response in the anova.bootstrap will not work generally
####  --depends on formula.tools package for checks

###no hypothesis testing (no bootstrap under null right now...)
###Need to check response and covariates more thoroughly at the front...
###add data=   option





library(formula.tools)
paired.boot <- function(formula, B=1000, seed=NULL){

  ##Check that a formula is input
  if(inherits(formula, "formula")==FALSE){
    stop("The input model must be a formula.")
  }

  ##Check that the response variable is in the correct format
  #print(lhs.vars(formula)); flush.console();
  #if(length(lhs.vars(formula)) >1){
  #  stop("Response must be a vector (not multivariate).")
  #}
  #else{
    resp <- eval(lhs(formula)) ##get the response variable

    if(is.matrix(resp)!=TRUE && is.atomic(resp)!=TRUE){
       stop("Response must be a vector or matrix.")
    }
    else if((dim(resp)[1]==0 || dim(resp)[2]==0) && length(resp)==0){
       stop("Response must have entries.")
    }
    else if(mode(resp)!="numeric"){
       stop("Response must be of type numeric.")
    }
    else if(anyNA(resp)==TRUE){
       stop("Response must not have any missing values.")
    }
  #}

 
  ##THIS IS NOT COMPLETE
  ##Check that the covariate variables are in the correct format
  if(length(rhs.vars(formula))<1){
    stop("Linear model must have at least 1 predictor variable.")
  }
  else{
    p <- length(rhs.vars(formula)) #number of covariates, does not include intercept -- does not account for matrices...

  }




  ##Check that the number of bootstrap samples is appropriate
  if(is.matrix(resp)==FALSE){ n <- length(resp) }
  else if(dim(resp)[1]!=1){ n <- dim(resp)[1] }
  else{ n <- dim(resp)[2] }  

  if(mode(B)!="numeric"){
    stop("Number of bootstrap samples, B, must be of type numeric.")
  }
  else if(is.atomic(B)!=TRUE){
    stop("Number of bootstrap samples, B, must be a constant.")
  }
  else if( B < n){
    stop("Number of bootstrap samples is recommended to be more than the data length.")
  }


  ##Check that the seed is numeric, or set a seed if unspecified
  if(is.null(seed)==TRUE){
    seed <- sample(seq(1,100000000), size=1)
  }
  else{
    if(mode(seed)!="numeric"){
      stop("The seed must be of type numeric.")
    }
    else if(is.atomic(seed)!=TRUE){
      stop("The seed must be a constant.")
    }
  }
  set.seed(seed)


  ##Remove any variables that are not used below
  remove(p, resp)


  #######################################################
  ## Least Squares Fit
  #######################################################
  obsDataregFit <- lm(formula)                        #fit the linear model specified in formula input
  resp <- obsDataregFit$model[,1]                     #get the response variable
  estParam <- matrix(obsDataregFit$coef, ncol=1)      #keep the param. estimates in a vector
  obsDataResid <- as.vector(residuals(obsDataregFit)) #keep the original residuals
  ParamNames <- names(obsDataregFit$coefficients)     #keep the coefficient name/association
  rownames(estParam) <- ParamNames                    #name the rows for the parameters so we know what they are

  ##Check there were not missing values in the covariates
  if(length(resp) != n){
    stop("There is at least one missing value within a predictor variable.  Missing values are not allowed.")
  }

  modelMat <- model.matrix(obsDataregFit)                   #model matrix (X)
  hatMat <- solve(t(modelMat) %*% modelMat) %*% t(modelMat) #projection matrix (X^TX)^-1 X^T



  ######################################################
  ## Bootstrap
  ######################################################
  ##Objects to keep Bootstrap Observations
  bootEstParam <- matrix(NA, nrow=B, ncol=dim(estParam)[1])  #bootstrap param. estimates
  colnames(bootEstParam) <- ParamNames

  boot.env <- new.env()                            #code scoping variable so we can sub the response in 'lm'

  for(i in 1:B){
    #Sample the indices of data frame
    index.sel <- sample(1:n, size=n, replace=TRUE)
    while(length(unique(index.sel)) <= dim(estParam)[1]){
      index.sel <- sample(1:n, size=n, replace=TRUE)
    }

    bootResp <- matrix(resp[index.sel], ncol=1)  #bootstrap response
    bootPred <- modelMat[index.sel, ]            #bootstrap predictors
    bootEstParam[i,]<- as.vector(solve(t(bootPred) %*% bootPred) %*% t(bootPred) %*% bootResp) #boot param est
  }


  #####################################################
  ## Returns
  #####################################################
  structure(invisible(list(bootEstParam=bootEstParam, 
                           origEstParam=estParam, seed=seed)))

}