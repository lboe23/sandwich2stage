#' @title Implement sandwich variance estimation in 2-stage regression.
#' @description This function provides an estimate of the sandwich variance from the 2-stage regression approach of regression calibration, which is obtained by stacking the stage 1 and stage 2 model estimating equations.
#' @details This function can be used for a linear stage 1 (calibration) model and when the stage 2 model is a linear, GLM, or Cox regression model. This function requires the following arguments: the fitted stage 1 and 2 model objects (fit using functions from the \code{survey} package, e.g. \code{svyglm()}), the names of the error-prone and estimated exposure variables in the data, and the names of the subject ID variables in the data sets used to fit the stage 1 and 2 models.
#' @param stage1.model the fitted stage 1 (calibration) model object fit using \code{svyglm()}
#' @param stage2.model the fitted stage 2 (outcome) model object fit using \code{svyglm()} or \code{svycoxph()}
#' @param xstar the name of the error-prone variable in the main study data
#' @param xhat the name of the estimated exposure variable from the stage 1 model
#' @param Stage1ID the name of the subject-level ID variable in the data used to fit the stage 1 model
#' @param Stage2ID the name of the subject-level ID variable in the data used to fit the stage 2 model
#' @import survey survival stats sandwich numDeriv stringr
#' @useDynLib sandwich2stage
#' @examples
#' #Begin by loading our package and then loading the data sandwichdata_SRS.
#' library(sandwich2stage)
#' data("sandwichdata_SRS")
#' #Next, we will create a survey design object (e.g. simple random sample).
#' sampdesign <- survey::svydesign(id=~1, data=sandwichdata_SRS)
#' #We may then fit the stage 1 and stage 2 models, saving the estimated nuisance parameters from the
#' #stage 1 model and using them to obtain an estimated of the unknown exposure (xhat), which will be
#' #used in the stage 2 model. We will then fit out stage 2 outcome model.
#' stage1.model<-survey::svyglm(xstarstar~xstar+z,design=sampdesign,family=gaussian(),subset=v==1)
#' alphas.stage1<-coef(stage1.model)
#' sampdesign <- update(sampdesign,xhat =predict(stage1.model,newdata=sampdesign$variables) )
#' stage2.model<-  survey::svyglm(y ~ xhat+z,design=sampdesign,family=binomial())
#' #Finally, we can obtain an estimate of the sandwich variance using our 2-step regression procedure.
#' sandwich.object<-sandwich2stage(stage1.model,stage2.model,
#' xstar="xstar",xhat="xhat",
#' Stage1ID="ID",Stage2ID="ID")
#' sandwichvar<-vcov(sandwich.object)
#' @references Binder, D. A. (1983). On the variances of asymptotically normal estimators from complex surveys. International Statistical Review/Revue Internationale de Statistique, pages 279–292.
#' @references Boos, D. and Stefanski, L. (2013). Essential Statistical Inference: Theory and Methods. Springer, New York, NY.
#' @references Lumley, T. (2011). Complex surveys: a guide to analysis using R, volume 565. John Wiley & Sons.
#' @references Lumley, T. and Scott, A. (2017). Fitting regression models to survey data. Statistical Science, pages 265–278.
#' @return returns a list containing the stage 1 and stage 2 model calls and formulas, the vector of estimated coefficients from the stage 1 and stage 2 models, the estimated sandwich variance matrix, and the matrix of influence functions.
#' @export
sandwich2stage <- function(stage1.model,stage2.model,xstar="xstar",xhat="xhat",Stage1ID="ID",Stage2ID="ID")
{
  UseMethod("sandwich2stage",stage2.model)
}

#'@export
sandwich2stage.lm<-function(stage1.model,stage2.model,xstar="xstar",xhat="xhat",Stage1ID="ID",Stage2ID="ID"){

  if(!inherits(stage1.model,"svyglm")){
    stop("Please run the stage 1 model using model-fitting functions from the survey package")
  }

  if(!inherits(stage2.model,"svyglm")){
    stop("Please run the stage 2 model using model-fitting functions from the survey package")
  }


  if(!(xstar %in% all.vars(formula(stage1.model)[[3]]))){
    stop(paste(xstar, " is not in the stage 1 model"))
  }

  if(!(xhat %in% all.vars(formula(stage2.model)[[3]]))){
    stop(paste(xhat, " is not in the stage 2 model"))
  }

  if(!is.null(stage2.model$na.action)) warning(length(stage2.model$na.action)," observations were missing in stage 2")

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

  #Get IDs for stage 1 and stage 2 participants
  allIDs.Stage1<-stage1.model.new$survey.design$variables[,(paste(Stage1ID))]
  allIDs.Stage2<-stage2.model.new$survey.design$variables[,(paste(Stage2ID))]

  fromstage2.instage1<-allIDs.Stage2 %in% allIDs.Stage1 #this is T/F for whether or not people from stage 2 in stage 1
  fromstage1.instage2<-allIDs.Stage1 %in% allIDs.Stage2  #this is T/F for whether or not people from stage 1 in stage 2

  #Error check for whether or not subset is nested
  if(any(fromstage1.instage2==FALSE)){
    stop("The subset used to fit the stage 1 model must be nested in the main study data.")
  }

  estfuns.lmNEW<-function (model) {
    Y=model$y
    XMat<-model.matrix(model)
    MU<-as.vector(XMat%*%coef(model))
    Resids<-(Y - MU)

    myestfun<-Resids * XMat
    return(myestfun)
  }

  #Don't divide by weights because estfuns.lmNEW already gives something unweighted
  estfun.stage1[fromstage2.instage1,] <- estfuns.lmNEW(stage1.model.new)[fromstage1.instage2,]

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
  sandwichresults<-list(stage1.call = stage1.model.new$call,stage2.call = stage2.model.new$call,
             stage1.formula = stage1.model.new$formula,stage2.formula = stage2.model.new$formula,
             betas = coef.all, sand.var=sandmat,infl.func=infl)
  class(sandwichresults)<-"SandwichObject"
  sandwichresults
}


# Function sandwich2stage.glm() - computes the sandwich when the stage 2 model is a generalized linear model (returns a list)
# inputs: stage 1 model object, stage 2 model object, name of xstar, name of xhat
##############################################################################################################################
#'@export
sandwich2stage.svyglm<-function(stage1.model,stage2.model,xstar="xstar",xhat="xhat",Stage1ID="ID",Stage2ID="ID"){

  if(!inherits(stage1.model,"svyglm")){
    stop("Please run the stage 1 model using model-fitting functions from the survey package")
  }

  if(!inherits(stage2.model,"svyglm")){
    stop("Please run the stage 2 model using model-fitting functions from the survey package")
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
  #xstarstar<-as.data.frame(stage2.model.new$data)[, paste0(formula.stage1[2])]

  #Now we can begin constructing meat of sandwich -- use existing functions
  #Create the estfun for the stage 1 model by first constructing a N by j_dim matrix of zeroes
  #  and insert estfun contributions for those who are in stage 1 model
  #  this way, all contributions for those NOT in stage 1 model are zero
  # with estimating equation contributions for those in the stage 1 model
  estfun.stage1<-matrix(0,nrow=N,ncol=j_dim)

  #Get IDs for stage 1 and stage 2 participants
  allIDs.Stage1<-stage1.model.new$survey.design$variables[,paste(Stage1ID)]
  allIDs.Stage2<-stage2.model.new$survey.design$variables[,paste(Stage2ID)]

  fromstage2.instage1<-allIDs.Stage2 %in% allIDs.Stage1 #this is T/F for whether or not people from stage 2 in stage 1
  fromstage1.instage2<-allIDs.Stage1 %in% allIDs.Stage2  #this is T/F for whether or not people from stage 1 in stage 2

  #Error check for whether or not subset is nested
  if(any(fromstage1.instage2==FALSE)==TRUE){
    print("The subset used to fit the stage 1 model must be nested in the main study data.")
  }

  estfuns.lmNEW<-function(model) {
    Y=model$y
    XMat<-model.matrix(model)
    MU<-as.vector(XMat%*%coef(model))
    Resids<-(Y - MU)

    myestfun<-Resids * XMat
    return(myestfun)
  }

  #Don't divide by weights because estfuns.lmNEW already gives something unweighted
  estfun.stage1[fromstage2.instage1,] <- estfuns.lmNEW(stage1.model.new)[fromstage1.instage2,]

  estfuns.glmNEW<-function(model)
  {
    xmat <- model.matrix(model)
    residuals(model, "working") * model$weights * xmat
  }

  #estfun for stage 2 is straightforward -- just pull from estfun
  estfun.stage2<-estfuns.glmNEW(stage2.model.new)/(stage2.model.new$prior.weights)

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
  sandwichresults<-list(stage1.call = stage1.model.new$call,stage2.call = stage2.model.new$call,
             stage1.formula = stage1.model.new$formula,stage2.formula = stage2.model.new$formula,
             betas = coef.all, sand.var=sandmat,infl.func=infl)
  class(sandwichresults)<-"SandwichObject"
  sandwichresults
}




#'@export
sandwich2stage.svycoxph<-function(stage1.model,stage2.model,xstar="xstar",xhat="xhat",Stage1ID="ID",Stage2ID="ID"){

  if(!inherits(stage1.model,"svyglm")){
    print("Please run the stage 1 model using model-fitting functions from the survey package")
  }

  if(!inherits(stage2.model,"svycoxph")){
    print("Please run the stage 2 model using model-fitting functions from the survey package")
  }

  if(!(xstar %in% all.vars(formula(stage1.model)[[3]]))){
    stop(paste(xstar, " is not in the stage 1 model"))
  }

  if(!(xhat %in% all.vars(formula(stage2.model)[[3]]))){
    stop(paste(xhat, " is not in the stage 2 model"))
  }


  stacker <- function(cmap, smap, istate, X, Y, strata, states, dropzero=TRUE) {
    from.state <- as.numeric(sub(":.*$", "", colnames(cmap)))
    to.state   <- as.numeric(sub("^.*:", "", colnames(cmap)))

    # just in case cmap has columns I don't need (I don't think this can
    #  happen
    check <- match(from.state, istate, nomatch=0)
    if (any(check==0)){
      # I think that this is impossible
      warning("extra column in cmap, this is a bug")  # debugging line
      # browser()
      cmap <- cmap[,check>0]
      from.state <- from.state[check>0]
      to.state <- to.state[check>0]
    }

    # Don't create X and Y matrices for transitions with no covariates, for
    #  coxph calls.  But I need them for survfit.coxph.
    zerocol <- apply(cmap==0, 2, all)
    if (dropzero && any(zerocol)) {
      cmap <- cmap[,!zerocol, drop=FALSE]
      smap <- smap[,!zerocol, drop=FALSE]
      smap[,] <- match(smap, sort(unique(c(smap)))) # relabel as 1, 2,...
      from.state <- from.state[!zerocol]
      to.state <- to.state[!zerocol]
    }

    endpoint <- c(0, match(attr(Y, "states"), states))
    endpoint <- endpoint[ 1 + Y[,ncol(Y)]]  # endpoint of each row, 0=censor

    # Jan 2021: changed from looping once per strata to once per transition.
    #  Essentially, a block of data for each unique column of cmap.  If two
    #  of those columns have the same starting state, it makes me nervous
    #  (statistically), but forge onward and sort the issues out in the
    #  fits.
    # Pass 1 to find the total data set size
    nblock <- ncol(cmap)
    n.perblock <- integer(nblock)
    for (i in 1:nblock) {
      n.perblock[i] <- sum(istate == from.state[i]) # can participate
    }

    # The constructed X matrix has a block of rows for each column of cmap
    n2 <- sum(n.perblock)  # number of rows in new data
    newX <- matrix(0, nrow=n2, ncol=max(cmap))
    k <- 0
    rindex <- integer(n2)   # original row for each new row of data
    newstat <- integer(n2)  # new status
    Xcols   <- ncol(X)      # number of columns in X
    for (i in 1:nblock) {
      subject <- which(istate == from.state[i]) # data rows in strata
      nr <- k + seq(along.with =subject)  # rows in the newX for this strata
      rindex[nr] <- subject
      nc <- cmap[,i]
      if (any(nc > Xcols)) { # constructed PH variables
        newX[nr, nc[nc>Xcols] ] <- 1
        nc <- nc[1:Xcols]
      }
      newX[nr, nc[nc>0]] <- X[subject, which(nc>0)] # row of cmap= col of X

      event.that.counts <- (endpoint[subject] == to.state[i])
      newstat[nr] <- ifelse(event.that.counts, 1L, 0L)
      k <- max(nr)
    }

    # which transition each row of newX represents
    transition <- rep(1:nblock, n.perblock)

    # remove any rows where X is missing
    #  these arise when a variable is used only for some transitions
    #  the row of data needs to be tossed for the given ones, but will be
    #    okay for other transitions
    keep <- !apply(is.na(newX), 1, any)
    if (!all(keep)) {
      newX <- newX[keep,, drop=FALSE]
      rindex <- rindex[keep]
      newstat <- newstat[keep]
      transition <- transition[keep]
    }

    if (ncol(Y) ==2) newY <- Surv(Y[rindex,1], newstat)
    else newY <- Surv(Y[rindex,1], Y[rindex,2], newstat)

    # newstrat will be an integer vector.
    newstrat <- smap[1, transition]   # start with strata induced by multi-state
    # then add any strata from the users strata() terms
    if (is.matrix(strata)){
      # this is the most complex case.
      maxstrat <- apply(strata, 2, max)  # max in each colum of strata
      mult <- cumprod(c(1, maxstrat))
      temp <- max(mult) * newstrat
      for (i in 1:ncol(strata)) {
        k <- smap[i+1, transition]
        temp <- temp + ifelse(k ==0, 0L, strata[i, rindex]* temp[i] -1L)
      }
      newstrat <- match(temp, sort(unique(temp)))
    }
    else if (length(strata) > 0) {
      # strata will be an integer vector with elements of 1, 2 etc
      mult <- max(strata)
      temp <- mult * newstrat + ifelse(smap[2,transition]==0, 0L, strata[rindex] -1L)
      newstrat <- match(temp, sort(unique(temp)))
    }

    # give variable names to the new data  (some names get used more than once)
    #    vname <- rep("", ncol(newX))
    #    vname[cmap[cmap>0]] <- colnames(X)[row(cmap)[cmap>0]]
    first <- match(sort(unique(cmap[cmap>0])), cmap) #first instance of each value
    vname <- rownames(cmap)[row(cmap)[first]]
    colnames(newX) <- vname
    list(X=newX, Y=newY, strata=as.integer(newstrat),
         transition= as.integer(transition), rindex=rindex)
  }

  survcheck2 <- function(y, id, istate=NULL, istate0="(s0)") {
    n <- length(id)
    ny <- ncol(y)
    # the next few line are a debug for my code; survcheck2 is not visible
    #  to users so only survival can call it directly
    if (!is.Surv(y) || is.null(attr(y, "states")) ||
        any(y[,ncol(y)] > length(attr(y, "states"))))
      stop("survcheck2 called with an invalid y argument")
    to.names <- c(attr(y, "states"), "(censored)")

    if (length(istate)==0) {
      inull<- TRUE
      cstate <- factor(rep(istate0, n))
    }
    else {
      if (length(istate) !=n) stop ("wrong length for istate")
      if (is.factor(istate)) cstate <- istate[, drop=TRUE] #drop unused levels
      else cstate <- as.factor(istate)
      inull <- FALSE
    }

    ystate <- attr(y, "states")
    # The vector of all state names is put in a nice printing order:
    #   initial states that are not destination states, then
    #   the destination states.  This keeps destinations in the order the
    #   user chose, while still putting initial states first.
    index <- match(levels(cstate), ystate, nomatch=0)
    states <- c(levels(cstate)[index==0], ystate)
    cstate2 <- factor(cstate, states)

    # Calculate the counts per id for each state, e.g., 10 subjects had
    #  3 visits to state 2, etc.
    # Count the censors, so that each subject gets a row in the table,
    #  but then toss that column
    tab1 <- table(id, factor(y[,ncol(y)], 0:length(ystate)))[,-1, drop=FALSE]
    tab1 <- cbind(tab1, rowSums(tab1))
    tab1.levels <- sort(unique(c(tab1)))  #unique counts
    if (length(tab1.levels) ==1) {
      # In this special case the table command does not give a matrix
      #  A data set with no events falls here, for instance
      events <- matrix(tab1.levels, nrow=1, ncol= (1 + length(ystate)))
    }
    else events <- apply(tab1, 2, function(x) table(factor(x, tab1.levels)))
    dimnames(events) = list("count"= tab1.levels,
                            "state"= c(ystate, "(any)"))
    # remove columns with no visits
    novisit <- colSums(events[-1,, drop=FALSE]) ==0
    if (any(novisit)) events <- events[,!novisit]

    # Use a C routine to create 3 variables: a: an index of whether this is
    #   the first (1) or last(2) observation for a subject, 3=both, 0=neither,
    #  b. current state, and
    #  c. sign of (start of this interval - end of prior one)
    # start by making stat2 = status re-indexed to the full set of states
    ny <- ncol(y)
    sindx <- match(ystate, states)
    stat2 <- ifelse(y[,ny]==0, 0L, sindx[pmax(1L, y[,ny])])
    id2 <- match(id, unique(id))  # we need unique integers
    if (ncol(y)==2) {
      index <- order(id, y[,1])
      check <- .Call("Cmulticheck", rep(0., n), y[,1], stat2, id2,
                     as.integer(cstate2), index- 1L)
    } else {
      index <- order(id, y[,2], y[,1])
      check <- .Call("Cmulticheck", y[,1], y[,2], stat2, id2,
                     as.integer(cstate2), index- 1L)
    }

    if (inull && ny> 2) {
      # if there was no istate entered in, use the constructed one from
      # the check routine
      # if ny=2 then every row starts at time 0
      cstate2 <-factor(check$cstate, seq(along.with=states), states)
    }

    # create the transtions table
    # if someone has an intermediate visit, i.e., (0,10, 0)(10,20,1), don't
    #  report the false 'censoring' in the transitions table
    # make it compact by removing any cols that are all 0, and rows of
    #  states that never occur (sometimes the starting state is a factor
    #  with unused levels)
    keep <- (stat2 !=0 | check$dupid > 1)  # not censored or last obs of this id
    transitions <- table(from=cstate2[keep],
                         to= factor(stat2[keep], c(seq(along.with=states), 0),
                                    c(states, "(censored)")),
                         useNA="ifany")
    nr <- nrow(transitions)
    never <- (rowSums(transitions) + colSums(transitions[,1:nr]))==0
    transitions <- transitions[!never, colSums(transitions)>0, drop = FALSE]

    # now continue with error checks
    # A censoring hole in the middle, such as happens with survSplit,
    #  uses "last state carried forward" in Cmultistate, which also
    #  sets the "gap" to 0 for the first obs of a subject
    mismatch <- (as.numeric(cstate2) != check$cstate)

    # gap = 0   (0, 10], (10, 15]
    # gap = 1   (0, 10], (12, 15]  # a hole in the time
    # gap = -1  (0, 10], (9, 15]   # overlapping times
    flag <- c(overlap= sum(check$gap < 0),
              gap =    sum(check$gap > 0 & !mismatch),
              jump =   sum(check$gap > 0 & mismatch),
              teleport = sum(check$gap==0 & mismatch & check$dupid%%2 ==0))

    rval <- list(states=states, transitions=transitions,
                 events= t(events), flag=flag,
                 istate= factor(check$cstate, seq(along.with=states), states))

    # add error details, if necessary
    if (flag["overlap"] > 0) {
      j <- which(check$gap < 0)
      rval$overlap <- list(row=j, id= unique(id[j]))
    }
    if (flag["gap"] > 0) {
      j <- which(check$gap > 0 & !mismatch)
      rval$gap <- list(row=j, id= unique(id[j]))
    }
    if (flag["jump"] > 0) {
      j <- which(check$gap > 0 & mismatch)
      rval$jump <- list(row=j, id= unique(id[j]))
    }
    if (flag["teleport"] > 0) {
      j <- which(check$gap==0 & mismatch)
      rval$teleport <- list(row=j, id= unique(id[j]))
    }

    rval
  }




  model.matrix.coxph<-function (object, data = NULL, contrast.arg = object$contrasts,
                                ...)
  {
    if (is.null(data) && !is.null(object[["x"]]))
      return(object[["x"]])
    Terms <- delete.response(object$terms)
    if (is.null(data))
      mf <- stats::model.frame(object)
    else {
      if (is.null(attr(data, "terms")))
        mf <- stats::model.frame(Terms, data, xlev = object$xlevels)
      else mf <- data
    }
    cluster <- attr(Terms, "specials")$cluster
    if (length(cluster)) {
      temp <- untangle.specials(Terms, "cluster")
      dropterms <- temp$terms
    }
    else dropterms <- NULL
    strats <- attr(Terms, "specials")$strata
    hasinteractions <- FALSE
    if (length(strats)) {
      stemp <- untangle.specials(Terms, "strata", 1)
      if (length(stemp$vars) == 1)
        strata.keep <- mf[[stemp$vars]]
      else strata.keep <- strata(mf[, stemp$vars], shortlabel = TRUE)
      istrat <- as.integer(strata.keep)
      for (i in stemp$vars) {
        if (any(attr(Terms, "order")[attr(Terms, "factors")[i,
        ] > 0] > 1))
          hasinteractions <- TRUE
      }
      if (!hasinteractions)
        dropterms <- c(dropterms, stemp$terms)
    }
    else istrat <- NULL
    if (length(dropterms)) {
      Terms2 <- Terms[-dropterms]
      X <- model.matrix(Terms2, mf, constrasts.arg = contrast.arg)
      temp <- attr(X, "assign")
      shift <- sort(dropterms)
      for (i in seq(along.with = shift)) temp <- temp + 1 *
        (shift[i] <= temp)
      attr(X, "assign") <- temp
    }
    else X <- model.matrix(Terms, mf, contrasts.arg = contrast.arg)
    Xatt <- attributes(X)
    if (hasinteractions)
      adrop <- c(0, untangle.specials(Terms, "strata")$terms)
    else adrop <- 0
    xdrop <- Xatt$assign %in% adrop
    X <- X[, !xdrop, drop = FALSE]
    attr(X, "assign") <- Xatt$assign[!xdrop]
    attr(X, "contrasts") <- Xatt$contrasts
    X
  }

  Coxgetdata<-function (fit, y = TRUE, x = TRUE, stratax = TRUE, weights = TRUE,
                        offset = FALSE)
  {
    ty <- fit[["y"]]
    tx <- fit[["x"]]
    twt <- fit[["weights"]]
    toff <- fit[["offset"]]
    if (!is.null(ty))
      n <- nrow(ty)
    else if (!is.null(tx))
      n <- nrow(tx)
    else n <- NULL
    coxms <- inherits(fit, "coxphms")
    Terms <- fit$terms
    if (!inherits(Terms, "terms"))
      stop("invalid terms component of fit")
    if (!is.null(n)) {
      if (is.null(fit$call$weights))
        twt <- rep(1, n)
      if (is.null(attr(terms(fit), "offset")))
        toff <- rep(0, n)
    }
    strat <- fit$strata
    strats <- attr(Terms, "specials")$strata
    if (length(strats) == 0 && length(strat) == 0 & !coxms)
      stratax <- FALSE
    if ((y && is.null(ty)) || (x && is.null(tx)) || (weights &&
                                                     is.null(twt)) || (stratax && is.null(strat)) || (offset &&
                                                                                                      is.null(toff))) {
      mf <- stats::model.frame(fit)
      n <- nrow(mf)
      if (weights) {
        twt <- model.extract(mf, "weights")
        if (is.null(twt))
          twt <- rep(1, n)
      }
      if (offset) {
        toff <- model.extract(mf, "offset")
        if (is.null(toff))
          toff <- rep(0, n)
      }
      if (inherits(fit, "coxphms") && ((y && is.null(ty)) ||
                                       ((x | stratax) && is.null(tx)))) {
        id <- model.extract(mf, "id")
        istate <- model.extract(mf, "istate")
        ty <- model.response(mf)
        if (is.null(fit$timefix) || fit$timefix)
          ty <- aeqSurv(ty)
        check <- survcheck2(ty, id, istate)
        tx <- model.matrix.coxph(fit, data = mf)
        if (length(strats)) {
          temp <- untangle.specials(Terms, "strata", 1)
          strat <- as.integer(strata(mf[temp$vars], shortlabel = T))
        }
        else strat <- NULL
        xstack <- stacker(fit$cmap, fit$stratum_map, as.integer(check$istate),
                          tx, ty, strat, check$states)
        tx <- xstack$X
        ty <- xstack$Y
        strat <- xstack$strata
        stratax <- TRUE
        if (offset)
          toff <- toff[xstack$rindex]
        if (weights)
          twt <- twt[xstack$rindex]
        ismiss <- is.nan(ty) | apply(is.na(tx), 1, any)
        if (offset)
          ismiss <- ismiss | is.nan(toff)
        if (weights)
          ismiss <- ismiss | is.nan(twt)
        if (any(ismiss)) {
          if (offset)
            toff <- toff[!ismiss]
          if (weights)
            twt <- twt[!ismiss]
          if (y)
            ty <- ty[!ismiss]
          if (x)
            tx <- tx[!ismiss, , drop = FALSE]
          if (stratax)
            strat <- strat[!ismiss]
        }
      }
      else {
        if (y && is.null(ty)) {
          ty <- model.extract(mf, "response")
          if (is.null(fit$timefix) || fit$timefix)
            ty <- aeqSurv(ty)
        }
        if ((x || stratax) && is.null(tx)) {
          if (stratax) {
            temp <- untangle.specials(Terms, "strata",
                                      1)
            strat <- strata(mf[temp$vars], shortlabel = T)
          }
          tx <- model.matrix.coxph(fit, data = mf)
        }
      }
    }
    temp <- list()
    if (y)
      temp$y <- ty
    if (x)
      temp$x <- tx
    if (stratax)
      temp$strata <- strat
    if (offset)
      temp$offset <- toff
    if (weights)
      temp$weights <- twt
    temp
  }

  #
  # The multi-state routines need to check their input data
  #  y = survival object
  #  id = subject identifier
  #  istate = starting state for each row, this can be missing.
  # The routine creates a proper current-state vector accounting for censoring
  #  (which should be called cstate for 'current state', but istate is retained)
  #  If a subject started in state A for instance, and had observations
  #  of (0, 10, B), (10, 15, censor), (15, 20, C) the current state is
  #  (A, B, B).  It checks that against the input 'istate' if it is present
  #  to generate checks.
  # Multiple other checks as well
  #


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
  datavariables_stage2<-Coxgetdata(stage2.model, y=TRUE, x=TRUE, stratax=TRUE)
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

  #Get IDs for stage 1 and stage 2 participants
  allIDs.Stage1<-stage1.model.new$survey.design$variables[,paste(Stage1ID)]
  allIDs.Stage2<-stage2.model.new$survey.design$variables[,paste(Stage2ID)]

  fromstage2.instage1<-allIDs.Stage2 %in% allIDs.Stage1 #this is T/F for whether or not people from stage 2 in stage 1
  fromstage1.instage2<-allIDs.Stage1 %in% allIDs.Stage2  #this is T/F for whether or not people from stage 1 in stage 2

  #Error check for whether or not subset is nested
  if(any(fromstage1.instage2==FALSE)==TRUE){
    print("The subset used to fit the stage 1 model must be nested in the main study data.")
  }

  estfuns.lmNEW<-function (model) {
    Y=model$y
    XMat<-model.matrix(model)
    MU<-as.vector(XMat%*%coef(model))
    Resids<-(Y - MU)

    myestfun<-Resids * XMat
    return(myestfun)
  }

  #Don't divide by weights because estfuns.lmNEW already gives something unweighted
  estfun.stage1[fromstage2.instage1,] <- estfuns.lmNEW(stage1.model.new)[fromstage1.instage2,]

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
      temp <- Coxgetdata(stage2.model.new, y=TRUE, x=TRUE, stratax=TRUE)
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
    resid <-  .Call("coxscore2",
                    y,
                    xmatstage2,
                    newstrat,
                    score,
                    weights[ord],
                    as.integer(method=='efron'))


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

  #Return list with stage 1 and stage 2 model calls, formulas, coefficients,and sandwich variances
  sandwichresults<-list(stage1.call = stage1.model.new$call,stage2.call = stage2.model.new$call,
             stage1.formula = stage1.model.new$formula,stage2.formula = stage2.model.new$formula,
             betas = coef.all, sand.var=sand2,infl.func=infl)
  class(sandwichresults)<-"SandwichObject"
  sandwichresults
}


#' @rdname print.SandwichObject
#' @export
print.SandwichObject <- function(x,digits = max(3L, getOption("digits") - 3L), ...)
{

  cat("\nStage 1 Model Call:  ",
      paste(deparse(x$stage1.call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("\nStage 2 Model Call:  ",
      paste(deparse(x$stage2.call), sep = "\n", collapse = "\n"), "\n\n", sep = "")

  cat("Coefficients:  ")
  cat("\n")

  print.default(format(x$betas, digits = digits),
                print.gap = 2, quote = FALSE)

  cat("\nSandwich Variance: ")
  cat("\n")
  print(x$sand.var)

}

#' @rdname coef.SandwichObject
#' @export
coef.SandwichObject<-function(object,...) {
  beta<-object$betas
  beta

}

#' @rdname vcov.SandwichObject
#' @export
vcov.SandwichObject<-function(object,...) {
  v<-object$sand.var
  v

}



