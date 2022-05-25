#Sandwich Variance Function
library(MASS)
library(mvtnorm)
library(survival)
library(sandwich)
library(xtable)
library(dplyr)
library(numDeriv)

####Functions in this code:
############################
# SandwichRegCal.lm() - computes the sandwich when the stage 2 model is a linear model (returns a list)
# SandwichRegCal.glm() - computes the sandwich when the stage 2 model is a generalized linear model (returns a list)
# SandwichRegCal.coxph() - computes the sandwich when the stage 2 model is a Cox PH model (returns a list)
# print.SandwichRegCal() - print function for SandwichRegCal
# coef.SandwichRegCal() - returns vector of coefficients from stage 1 and stage 2 models 
# vcov.SandwichRegCal() - returns sandwich variance matrix from stacked estimating equation approach 

####Things to develop in this code:
################ 1 ######################
#Adopt the trick for hidden functions inside a package 
# example is survey:::estfuns.lm, survey:::estfuns.glm, and survey:::estfuns.coxph --
# idea is that when we call a function, like potentially SandwichRegCal(), it will call an internal function
#   in order to know whether to use the .lm, .glm, or .coxph hidden function
################ 2 ######################
#Need to add more error checks
# Specific one I can think of - the third error check below asks the user to...
#   "Please run the stage 1 model and stage 2 model using the same data set and survey design"
#      ...But this isn't entirely what we want The stage 2 model could be fit with a design created using subset function.
#          I think we need to go with idea that we might have 2 design objects, and they have a row level ID used to join them. 
################ 3 ######################
#Cox function - make this work for more complex settings like time varying covariate
# Have commented



# Function SandwichRegCal.lm() - computes the sandwich when the stage 2 model is a linear model (returns a list)
# inputs: stage 1 model object, stage 2 model object, name of xstar, name of xhat
##################################################################################################################

SandwichRegCal.lm<-function(stage1.model,stage2.model,xstar="xstar",xhat="xhat"){
  
  if(!inherits(stage1.model,"svyglm")){
    stop("Please run the stage 1 model using model-fitting functions from the survey package")
  }
  
  if(!inherits(stage2.model,"svyglm")){
    stop("Please run the stage 2 model using model-fitting functions from the survey package")
  }
  
  if(!stage1.model$call$design==stage2.model$call$design){
    print("Please run the stage 1 model and stage 2 model using the same data set and survey design")
  }
  
  if(!(xstar %in% all.vars(formula(stage1.model)[[3]]))){
    stop(paste(xstar, " is not in the stage 1 model"))
  }
  
  if(!(xhat %in% all.vars(formula(stage2.model)[[3]]))){
    stop(paste(xhat, " is not in the stage 2 model"))
  }
  
  #rewrite stage 1 formula so that the error-prone covariate comes first 
  vars.stage1<-  all.vars(formula(stage1.model)[[3]])
  vars.stage1.all<- c(xstar,setdiff(vars.stage1, xstar))

  formula.stage1<-as.formula( paste(formula(stage1.model)[2], "~", paste(vars.stage1.all, collapse = "+")))

  #rewrite stage 2 formula so that the xhat comes first 
  vars.stage2<-  all.vars(formula(stage2.model)[[3]])
  vars.stage2.all<- c(xhat,setdiff(vars.stage2, xhat))
  
  formula.stage2<-as.formula( paste(formula(stage2.model)[2], "~", paste(vars.stage2.all, collapse = "+")))
  
  #Update the models with new variable order 
  stage1.model.new<-update(stage1.model,formula.stage1)
  stage2.model.new<-update(stage2.model,formula.stage2)
  
  #fit naive stage 2 model - will use its model matrix 
  stage2.model.naive<-update(stage2.model,update(formula.stage1,paste(formula.stage2[[2]],"~.")))

  #we will have j unknown alphas and k unknown betas - assign these for constructing the sandwich 
  j_dim<-length(coef(stage1.model.new))
  k_dim<-length(coef(stage2.model.new))  
  
  #Need to obtain the number of subjects in main study (N) fit in stage 2 model
  #  and number of subjects in subset (n_sub) from stage 1 model
  n_sub<-nobs(stage1.model.new)
  N<-nobs(stage2.model.new)
  
  #Save variable being used for calibration on LHS of equation (xstarstar)
  xstarstar<-as.data.frame(stage2.model.new$data)[, paste0(formula.stage1[2])]

  #Now we can begin constructing meat of sandwich -- use existing functions
  #Create the estfun for the stage 1 model by first constructing a N by j_dim matrix of zeroes
  #  and insert estfun contributions for those who are in stage 1 model
  #  this way, all contributions for those NOT in stage 1 model are zero
  # with estimating equation contributions for those in the stage 1 model
  estfun.stage1<-matrix(0,nrow=N,ncol=j_dim)
  is.calibration <- !is.na(xstarstar)
  
  estfuns.lmNEW<-function (model) {
    Y=model$y
    XMat<-model.matrix(model)
    MU<-as.vector(XMat%*%coef(model))
    Resids<-(Y - MU)
    
    myestfun<-Resids * XMat
    return(myestfun)
  }
  
  #Don't divide by weights because survey:::estfuns.lm already gives something unweighted
  #estfun.stage1[is.calibration,] <- survey:::estfuns.lm(stage1.model.new) #/weights(stage1.model.new$survey.design,"sampling")
  estfun.stage1[is.calibration,] <- estfuns.lmNEW(stage1.model.new) 
  
  #estfun for stage 2 is straightforward -- just pull from estfun 
  estfun.stage2<-estfuns.lmNEW(stage2.model.new) #/weights(stage2.model.new$survey.design,"sampling")
  
  #Combine estfuns for both stage 1 and stage 2 models 
  estfun.all<-cbind(estfun.stage1,estfun.stage2)

 ###Now let's obtain the pieces for matrix A 
  #Easiest part: upper-right quadrant (all elements are zero), j times k
  A.upperright<- matrix(0,nrow=j_dim,ncol=k_dim)
  
  #Here is the upper left computed using the Hessian from stage 1 
  A.upperleft<- -solve(stage1.model.new$naive.cov/mean(stage1.model.new$prior.weights))/N
  
  #Bottom-right quadrant is just Hessian for stage 2 model
  A.bottomright<- -solve(stage2.model.new$naive.cov/mean(stage2.model.new$prior.weights))/N
  
  #Bottom-left quadrant is a little trickier depends on stage 2 model.
  #For now, let's assume GLM for stage 2. 

 #Save alphas - coefficients from stage 1 model
  alphas.stage1<-coef(stage1.model.new)

 #Now, for stage 2 model, we write an estfun function as a function of alphas, holding beta constant 
  stage2.alphas.estfuns.lm<-function (alphas) {
    MyX_Naive<-model.matrix(stage2.model.naive)
    MyX_RC<-model.matrix(stage2.model.new)
    y.j=stage2.model.new$y
    myxhat<-as.vector(MyX_Naive%*%alphas)
    xmat.j <- as.matrix(cbind(MyX_RC[,1],myxhat,MyX_RC[,c(3:ncol(MyX_RC))]))
    mu.j<-as.vector(stage2.model.new$family$linkinv(xmat.j%*%coef(stage2.model.new)))
    r.j<-(y.j - mu.j)/(stage2.model.new$family$variance(mu.j))

    myestfun<-r.j * xmat.j*stage2.model.new$prior.weights
    return(myestfun)
  }
  
  #Now, use Jacobian function from numderiv to obtain derivatives of stage 2 estimating equation wrt alpha
  #Then split it up, because it computes large matrix with each person's contributions for each parameter
  JacobianAll.stage2<-jacobian(func=stage2.alphas.estfuns.lm, x=alphas.stage1)
  JacobianList.stage2<-lapply(split(JacobianAll.stage2,rep(c(1:N),each=1)),matrix,nrow=k_dim)
  
  #Add these matrices for subjects 1...N then divide by N
  add <- function(x) Reduce("+", x)
  
  A.bottomleft<-add(JacobianList.stage2)*mean(stage2.model.new$prior.weights)/N
  

  #Now we can combine all four quadrants of our A matrix and obtain AAll
  AAll<-rbind(cbind(A.upperleft,A.upperright),cbind(A.bottomleft,A.bottomright))
  
  #Invert this matrix; needed for sandwich 
  A.inv<-solve(AAll)
  
  #Compute influence functions
  infl<- as.matrix(estfun.all)%*%t(A.inv)/N
   
  #Compute the sandwich 2 ways for SRS -- directly in sandwich from, or using vcov(svytotal())
  sandmat<-vcov(svytotal(infl,stage2.model$survey.design))
  
  #Combine all estimated coefficients from stage 1 and stage 2
  coef.all<-c(coef(stage1.model.new),coef(stage2.model.new))
  
  
  #Give names
  colnames(sandmat)<-rownames(sandmat)<-names(coef.all)<-
    c(paste0("stage1.", names(coef(stage1.model.new))),paste0("stage2.", names(coef(stage2.model.new))))

  
  #Return list with stage 1 and stage 2 model calls, formulas, coefficients,and sandwich variances 
  list(stage1.call = stage1.model.new$call,stage2.call = stage2.model.new$call,
       stage1.formula = stage1.model.new$formula,stage2.formula = stage2.model.new$formula,
       coef = coef.all, sand.var=sandmat)
  
  
}

# Function SandwichRegCal.glm() - computes the sandwich when the stage 2 model is a generalized linear model (returns a list)
# inputs: stage 1 model object, stage 2 model object, name of xstar, name of xhat
##############################################################################################################################



SandwichRegCal.glm<-function(stage1.model,stage2.model,xstar="xstar",xhat="xhat"){
  
  if(!inherits(stage1.model,"svyglm")){
    stop("Please run the stage 1 model using model-fitting functions from the survey package")
  }
  
  if(!inherits(stage2.model,"svyglm")){
    stop("Please run the stage 2 model using model-fitting functions from the survey package")
  }
  
  if(!stage1.model$call$design==stage2.model$call$design){
    print("Please run the stage 1 model and stage 2 model using the same data set and survey design")
  }
  
  if(!(xstar %in% all.vars(formula(stage1.model)[[3]]))){
    stop(paste(xstar, " is not in the stage 1 model"))
  }
  
  if(!(xhat %in% all.vars(formula(stage2.model)[[3]]))){
    stop(paste(xhat, " is not in the stage 2 model"))
  }
  
  #rewrite stage 1 formula so that the error-prone covariate comes first 
  vars.stage1<-  all.vars(formula(stage1.model)[[3]])
  vars.stage1.all<- c(xstar,setdiff(vars.stage1, xstar))
  
  formula.stage1<-as.formula( paste(formula(stage1.model)[2], "~", paste(vars.stage1.all, collapse = "+")))
  
  #rewrite stage 2 formula so that the xhat comes first 
  vars.stage2<-  all.vars(formula(stage2.model)[[3]])
  vars.stage2.all<- c(xhat,setdiff(vars.stage2, xhat))
  
  formula.stage2<-as.formula( paste(formula(stage2.model)[2], "~", paste(vars.stage2.all, collapse = "+")))
  
  #Update the models with new variable order 
  stage1.model.new<-update(stage1.model,formula.stage1)
  stage2.model.new<-update(stage2.model,formula.stage2)
  
  #fit naive stage 2 model - will use its model matrix 
  stage2.model.naive<-update(stage2.model,update(formula.stage1,paste(formula.stage2[[2]],"~.")))
  
  #we will have j unknown alphas and k unknown betas - assign these for constructing the sandwich 
  j_dim<-length(coef(stage1.model.new))
  k_dim<-length(coef(stage2.model.new))  
  
  #Need to obtain the number of subjects in main study (N) fit in stage 2 model
  #  and number of subjects in subset (n_sub) from stage 1 model
  n_sub<-nobs(stage1.model.new)
  N<-nobs(stage2.model.new)
  
  #Save variable being used for calibration on LHS of equation (xstarstar)
  xstarstar<-as.data.frame(stage2.model.new$data)[, paste0(formula.stage1[2])]
  
  #Now we can begin constructing meat of sandwich -- use existing functions
  #Create the estfun for the stage 1 model by first constructing a N by j_dim matrix of zeroes
  #  and insert estfun contributions for those who are in stage 1 model
  #  this way, all contributions for those NOT in stage 1 model are zero
  # with estimating equation contributions for those in the stage 1 model
  estfun.stage1<-matrix(0,nrow=N,ncol=j_dim)
  is.calibration <- !is.na(xstarstar)
  
  estfuns.lmNEW<-function (model) {
    Y=model$y
    XMat<-model.matrix(model)
    MU<-as.vector(XMat%*%coef(model))
    Resids<-(Y - MU)
    
    myestfun<-Resids * XMat
    return(myestfun)
  }
  
  #Don't divide by weights because survey:::estfuns.lm already gives something unweighted
  #  estfun.stage1[is.calibration,] <- survey:::estfuns.lm(stage1.model.new) #/weights(stage1.model.new$survey.design,"sampling")
  estfun.stage1[is.calibration,] <- estfuns.lmNEW(stage1.model.new) 
  
  #estfun for stage 2 is straightforward -- just pull from estfun 
  #estfun.stage2<-(survey:::estfuns.glm(stage2.model.new)/weights(stage2.model.new$survey.design,"sampling"))
  estfun.stage2<-survey:::estfuns.glm(stage2.model.new)/(stage2.model.new$prior.weights)

  #Combine estfuns for both stage 1 and stage 2 models 
  estfun.all<-cbind(estfun.stage1,estfun.stage2)
  
  #Function to add matrices
  add <- function(x) Reduce("+", x)
  
  ###Now let's obtain the pieces for matrix A 
  #Easiest part: upper-right quadrant (all elements are zero), j times k
  A.upperright<- matrix(0,nrow=j_dim,ncol=k_dim)
  
  #Save alphas and betas - coefficients from stage 1 and stage 2 models
  alphas.stage1<-coef(stage1.model.new)
  
  #Here is the upper left computed using the Hessian from stage 1 
  A.upperleft<- -solve(stage1.model.new$naive.cov/mean(stage1.model.new$prior.weights))/N
  
  #Bottom-right quadrant is just Hessian for stage 2 model
  A.bottomright<- -solve(stage2.model.new$naive.cov/mean(stage2.model.new$prior.weights))/N
  
  #Bottom-left quadrant is a little trickier depends on stage 2 model.
  #For now, let's assume GLM for stage 2. 
  
  #Now, for stage 2 model, we write an estfun function as a function of alphas, holding beta constant 
  stage2.alphas.estfuns.glm<-function (alphas) {
    MyX_Naive<-model.matrix(stage2.model.naive)
    MyX_RC<-model.matrix(stage2.model.new)
    y.j=stage2.model.new$y
    myxhat<-as.vector(MyX_Naive%*%alphas)
    xmat.j <- as.matrix(cbind(MyX_RC[,1],myxhat,MyX_RC[,c(3:ncol(MyX_RC))]))
    mu.j<-as.vector(stage2.model.new$family$linkinv(xmat.j%*%coef(stage2.model.new)))
    r.j<-(y.j - mu.j)/(stage2.model.new$family$variance(mu.j))
    workingweights.j <- as.vector(((stage2.model.new$family$variance(mu.j))^2)/stage2.model.new$family$variance(mu.j))
    
    myestfun<-r.j * workingweights.j * xmat.j*stage2.model.new$prior.weights
    
    return(myestfun)
  }
  
  
  #Now, use Jacobian function from numderiv to obtain derivatives of stage 2 estimating equation wrt alpha
  #Then split it up, because it computes large matrix with each person's contributions for each parameter
  JacobianAll.stage2<-jacobian(func=stage2.alphas.estfuns.glm, x=alphas.stage1)
  JacobianList.stage2<-lapply(split(JacobianAll.stage2,rep(c(1:N),each=1)),matrix,nrow=k_dim)

  #Add these matrices for subjects 1...N then divide by N
  A.bottomleft<-add(JacobianList.stage2)*mean(stage2.model.new$prior.weights)/N
  
  #Now we can combine all four quadrants of our A matrix and obtain AAll
  AAll<-rbind(cbind(A.upperleft,A.upperright),cbind(A.bottomleft,A.bottomright))
  
  #Invert this matrix; needed for sandwich 
  A.inv<-solve(AAll)
  
  #Compute influence functions
  infl<- as.matrix(estfun.all)%*%t(A.inv)/N
  
  #Compute the sandwich using vcov(svytotal())
  sandmat<-vcov(svytotal(infl,stage2.model$survey.design))
  
  
  #Combine all estimated coefficients from stage 1 and stage 2
  coef.all<-c(coef(stage1.model.new),coef(stage2.model.new))
  
  
  #Give names
  colnames(sandmat)<-rownames(sandmat)<-names(coef.all)<-
    c(paste0("stage1.", names(coef(stage1.model.new))),paste0("stage2.", names(coef(stage2.model.new))))
  
  
  #Return list with stage 1 and stage 2 model calls, formulas, coefficients,and sandwich variances 
  list(stage1.call = stage1.model.new$call,stage2.call = stage2.model.new$call,
       stage1.formula = stage1.model.new$formula,stage2.formula = stage2.model.new$formula,
       coef = coef.all, sand.var=sandmat, estfun=estfun.all)
  
}



# Function SandwichRegCal.coxph() - computes the sandwich when the stage 2 model is a Cox PH model (returns a list)
# inputs: stage 1 model object, stage 2 model object, name of xstar, name of xhat
###################################################################################################################


SandwichRegCal.coxph<-function(stage1.model,stage2.model,xstar="xstar",xhat="xhat"){
  
  if(!inherits(stage1.model,"svyglm")){
    print("Please run the stage 1 model using model-fitting functions from the survey package")
  }
  
  if(!inherits(stage2.model,"svycoxph")){
    print("Please run the stage 2 model using model-fitting functions from the survey package")
  }
  
  if(!stage1.model$call$design==stage2.model$call$design){
    print("Please run the stage 1 model and stage 2 model using the same data set and survey design")
  }
  
  if(!(xstar %in% all.vars(formula(stage1.model)[[3]]))){
    stop(paste(xstar, " is not in the stage 1 model"))
  }
  
  if(!(xhat %in% all.vars(formula(stage2.model)[[3]]))){
    stop(paste(xhat, " is not in the stage 2 model"))
  }
  
  #rewrite stage 1 formula so that the error-prone covariate comes first 
  vars.stage1<-  all.vars(formula(stage1.model)[[3]])
  vars.stage1.all<- c(xstar,setdiff(vars.stage1, xstar))
  
  formula.stage1<-as.formula( paste(formula(stage1.model)[2], "~", paste(vars.stage1.all, collapse = "+")))
  
  #Find out if strata in stage 2 model
  stage2formula_char <- as.character(formula(stage2.model))[3]
  
  #Split up all character elements of formula so we can see if strata
  stage2formula_char_split<-strsplit(stage2formula_char, " ")[[1]]
  
  #Save part of formula with strata and use in writing of stage 2 model
  strataforstage2<-stage2formula_char_split[str_detect(stage2formula_char_split, "strata")]
  
  #rewrite stage 2 formula so that the xhat comes first 
  datavariables_stage2<-survival:::coxph.getdata(stage2.model, y=TRUE, x=TRUE, stratax=TRUE)
  datavariables_stage2_x<-datavariables_stage2$x
  vars.stage2<-colnames(datavariables_stage2_x)
  vars.stage2.all<- c(xhat,setdiff(vars.stage2, xhat))
  
  if(identical(strataforstage2,character(0))){
    formula.stage2<-as.formula( paste(formula(stage2.model)[2], "~", paste(vars.stage2.all, collapse = "+")))
  } else {
    formula.stage2<-as.formula( paste(formula(stage2.model)[2], "~",paste(strataforstage2)," +", paste(vars.stage2.all, collapse = "+")))
  }
  
  
  
  #Update the models with new variable order 
  stage1.model.new<-update(stage1.model,formula.stage1)
  stage2.model.new<-update(stage2.model,formula.stage2)
  
  #fit naive stage 2 model - will use its model matrix 
  #stage2.model.naive<-update(stage2.model,update(formula.stage1,paste(formula.stage2[[2]],"~.")))
  
  #we will have j unknown alphas and k unknown betas - assign these for constructing the sandwich 
  j_dim<-length(coef(stage1.model.new))
  k_dim<-length(coef(stage2.model.new))  
  
  #Need to obtain the number of subjects in main study (N) fit in stage 2 model
  #  and number of subjects in subset (n_sub) from stage 1 model
  n_sub<-nobs(stage1.model.new)
  N<-stage2.model.new$n
  
  
  #Save variable being used for calibration on LHS of equation (xstarstar)
  xstarstar<-stage2.model.new$survey.design$variables[, paste0(formula.stage1[2])]
  
  #Now we can begin constructing meat of sandwich -- use existing functions
  #Create the estfun for the stage 1 model by first constructing a N by j_dim matrix of zeroes
  #  and insert estfun contributions for those who are in stage 1 model
  #  this way, all contributions for those NOT in stage 1 model are zero
  # with estimating equation contributions for those in the stage 1 model
  estfun.stage1<-matrix(0,nrow=N,ncol=j_dim)
  is.calibration <- !is.na(xstarstar)
  estfun.stage1[is.calibration,] <- survey:::estfuns.lm(stage1.model.new)
  
  #estfun for stage 2 is straightfroward -- just pull from estfun 
  estfun.stage2<-as.matrix(estfun(stage2.model.new))
  
  #Combine estfuns for both stage 1 and stage 2 models 
  estfun.all<-cbind(estfun.stage1,estfun.stage2)
  
  #Function to add matrices
  add <- function(x) Reduce("+", x)
  
  ###Now let's obtain the pieces for matrix A 
  #Easiest part: upper-right quadrant (all elements are zero), j times k
  A.upperright<- matrix(0,nrow=j_dim,ncol=k_dim)
  
  #Here is the upper left computed using the Hessian from stage 1 
  A.upperleft<- -solve(stage1.model.new$naive.cov)/N
  
  #Bottom-right quadrant is just Hessian for stage 2 model
  A.bottomright<- -solve(stage2.model.new$inv.info)/N
  
  #Bottom-left quadrant is a little trickier depends on stage 2 model.
  #For now, let's assume GLM for stage 2. 
  
  #Save alphas - coefficients from stage 1 model
  alphas.stage1<-coef(stage1.model.new)
  
  #Now, for stage 2 model, we write a function for obtaining cox score residuals
  #   as a function of alphas, holding beta constant 
  stage2.alphas.estfuns.coxph<-function(alphas){
    type<-"score"
    
    otype <- type
    n <- length(stage2.model.new$residuals)
    rr <- stage2.model.new$residuals
    y <- stage2.model.new$y
    x <- stage2.model.new[['x']]  # avoid matching stage2.model.new$xlevels
    
    vv <- drop(stage2.model.new$naive.var)
    
    if (is.null(vv)) vv <- drop(stage2.model.new$var)
    weights <- stage2.model.new$weights
    if (is.null(weights)) weights <- rep(1,n)
    strat <- stage2.model.new$strata
    method <- stage2.model.new$method
    
    # I need Y, and perhaps the X matrix (and strata)
    Terms <- stage2.model.new$terms
    if (!inherits(Terms, 'terms'))
      stop("invalid terms component of object")
    strats <- attr(Terms, "specials")$strata
    
    
    
    if (is.null(y)  ||  (is.null(x) && type!= 'deviance')) {
      temp <- survival:::coxph.getdata(stage2.model.new, y=TRUE, x=TRUE, stratax=TRUE)
      y <- temp$y
      x <- temp$x
      
      myvarscox<-colnames(temp$x)[-1]
      
      zvars<-x[,(myvarscox)]
      
      svydata<-stage2.model.new$survey.design$variables
      xmatstage1<-as.matrix(cbind(rep(1,N),svydata[,all.vars(formula(stage1.model.new)[[3]])]))
      myxhat<-as.vector(xmatstage1%*%alphas)
      xmatstage2 <- as.matrix(cbind(myxhat,zvars))
      
      if (length(strats)) strat <- temp$strata
    }
    ny <- ncol(y)
    status <- y[,ny,drop=TRUE]
    
    
    
    nstrat <- as.numeric(strat)
    nvar <- ncol(x)
    if (is.null(strat)) {
      ord <- order(y[,ny-1], -status)
      newstrat <- rep(0,n)
    }else {
      ord <- order(nstrat, y[,ny-1], -status)
      newstrat <- c(diff(as.numeric(nstrat[ord]))!=0 ,1)
    }
    newstrat[n] <- 1
    
    # sort the data
    x <- x[ord,]
    xmatstage2 <- xmatstage2[ord,]
    
    y <- y[ord,]
    linpred<-exp(xmatstage2%*%coef(stage2.model.new))
    mymeans<-apply(xmatstage2,2,mean)
    expmean<-as.vector(exp(mymeans%*%coef(stage2.model.new)))
    score<-as.vector(linpred/expmean)
    
    
    
    #Score
    storage.mode(y) <- storage.mode(x) <- storage.mode(xmatstage2) <- "double"
    storage.mode(newstrat) <- "integer"
    storage.mode(score) <- storage.mode(weights) <- "double"
    #if (ny==2) {
    resid <-  .Call(survival:::Ccoxscore2, 
                    y, 
                    xmatstage2, 
                    newstrat,
                    score,
                    weights[ord],
                    as.integer(method=='efron'))
    
    summary(resid)
    #} else {
    #  resid<- .Call(survival:::Cagscore2,
    #                y, 
    #                x, 
    #                newstrat,
    #                score,
    #                weights[ord],
    #                as.integer(method=='efron'))
    #}
    
    if (nvar >1) {
      rr <- matrix(0, n, nvar)
      rr[ord,] <- resid
      dimnames(rr) <- list(names(stage2.model.new$residuals), 
                           names(stage2.model.new$coefficients))
    }else {
      rr[ord] <- resid
    }
    
    return(rr)
  }
  
  
  JacobianAll.stage2<-jacobian(func=stage2.alphas.estfuns.coxph, x=alphas.stage1)
  JacobianList.stage2<-lapply(split(JacobianAll.stage2,rep(c(1:N),each=1)),matrix,nrow=k_dim)
  A.bottomleft<-add(JacobianList.stage2)/N
  
  #Now we can combine all four quadrants of our A matrix and obtain AAll
  AAll<-rbind(cbind(A.upperleft,A.upperright),cbind(A.bottomleft,A.bottomright))
  
  #Invert this matrix; needed for sandwich 
  A.inv<-solve(AAll)
  
  #Compute influence functions
  infl<- as.matrix(estfun.all)%*%t(A.inv)/N
  
  #Compute the sandwich 2 ways for SRS -- directly in sandwich from, or using vcov(svytotal())
  sand2<-vcov(svytotal(infl,stage2.model$survey.design))
  
  #Combine all estimated coefficients from stage 1 and stage 2
  coef.all<-c(coef(stage1.model.new),coef(stage2.model.new))
  
  
  #Give names
  colnames(sand2)<-rownames(sand2)<-names(coef.all)<-
    c(paste0("stage1.", names(coef(stage1.model.new))),paste0("stage2.", names(coef(stage2.model.new))))
  
  #look at results
  
  list(stage1.call = stage1.model.new$call,stage2.call = stage2.model.new$call,
       stage1.formula = stage1.model.new$formula,stage2.formula = stage2.model.new$formula,
       coef = coef.all, sand.var=sand2)
  
  
}


# Function print.SandwichRegCal() - print function for SandwichRegCal
# inputs: object from SandwichRegCal.lm, SandwichRegCal.glm, or SandwichRegCal.coxph
######################################################################################

print.SandwichRegCal <- function(object, digits = max(3L, getOption("digits") - 3L), ...)
{
  
  cat("\nStage 1 Model Call:  ",
      paste(deparse(object$stage1.call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("\nStage 2 Model Call:  ",
      paste(deparse(object$stage2.call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  
  cat("Coefficients:  ")
  cat("\n")
  
  print.default(format(object$coef, digits = digits),
                print.gap = 2, quote = FALSE)
  
  cat("\nSandwich Variance: ")
  cat("\n")
  print(object$sand.var)
  
}

# Function coef.SandwichRegCal() - returns vector of coefficients from stage 1 and stage 2 models 
# inputs: object from SandwichRegCal.lm, SandwichRegCal.glm, or SandwichRegCal.coxph
##################################################################################################

coef.SandwichRegCal<-function(object) {
  beta<-object$sand.var
  beta
}

# Function vcov.SandwichRegCal() - returns sandwich variance matrix from stacked estimating equation approach 
# inputs: object from SandwichRegCal.lm, SandwichRegCal.glm, or SandwichRegCal.coxph
###############################################################################################################

vcov.SandwichRegCal<-function(object) {
  v<-object$cov.unscaled
  v
}





