###Currently this is only set up for one-way ANOVA
###only implements residual bootstrap
###pretty slow since a matrix is being inverted twice in each bootstrap iteration -- update to directly implement formula
###put in option to set a seed
###Need to check it picks up the response appropriately every time.... seems to work on examples where dataset has response last naturally...
###add data=   option

library(formula.tools)
anova.boot <- function(formula, B=1000, type="residual", seed=NULL){

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



  ##Check that the type of bootstrap is implementable
  if(type!="residual"){
    stop("Only residual bootstrap is allowed for type.")
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
    # Change this to a warning instead
    # stop("Number of bootstrap samples is recommended to be more than the data length.")
  }


  ##Check the the seed is numeric, or set a seed if unspecified
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




  ##################################################
  ## Least Squares Fit
  ##  -- update to use Type II SS...
  ##################################################
  obsDataregFit <- aov(formula)                       #fit the linear model specified in formula input
  resp <- obsDataregFit$model[,1]                     #get the response variable
  estParam <- matrix(obsDataregFit$coef, ncol=1)      #keep the param. estimates in a vector
  obsDataResid <- as.vector(residuals(obsDataregFit)) #keep the original residuals
  ParamNames <- names(obsDataregFit$coefficients)     #keep the coefficient name/association
  rownames(estParam) <- ParamNames                    #name the rows for the parameters so we know what they are

  obsFStats <- summary(obsDataregFit)[[1]][,4]        #get the F test stats from orig fit
  obsFStats <- obsFStats[-length(obsFStats)]          #remove NA next to MSE

  modelMat <- model.matrix(obsDataregFit)                    #model matrix (X)
  quadMat <- solve(t(modelMat) %*% modelMat) %*% t(modelMat) #projection matrix (X^TX)^-1 X^T



  ###################################################
  ## Bootstrap
  ###################################################
  ##Objects to keep Bootstrap Observations
  bootEstParam <- matrix(NA, nrow=B, ncol=dim(estParam)[1])  #bootstrap param. estimates
  colnames(bootEstParam) <- ParamNames

  bootFStats <- matrix(NA, nrow=B, ncol=length(obsFStats))   #bootstrap F statistics
  boot.env <- new.env()                            #code scoping variable so we can sub the response in 'lm'

  for(i in 1:B){
    bootResid <- matrix(sample(obsDataResid, replace=TRUE), ncol=1)  #bootstrap residuals

    #Bootstrap under H0 to get F test statistics
    boot.env$bootResp <- mean(resp)+ bootResid       #bootstrap response under H0
    newFormula <- update.formula(formula, boot.env$bootResp~.)
    bootregFit <- lm(newFormula, boot.env)
    holdFstat <- anova(bootregFit)[[4]]
    bootFStats[i,]<- holdFstat[-length(holdFstat)]

    #bootstrap to get param estimates
    bootEstParam[i,]<- as.vector( estParam + quadMat %*% bootResid )  #bootstrap parameter estimates
  }

  structure(invisible(list(bootEstParam=bootEstParam, bootFStats=bootFStats,
                           origEstParam=estParam, origFStats=obsFStats)))

}





