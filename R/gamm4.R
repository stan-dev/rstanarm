# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
# Copyright (C) 2009, 2010, 2011, 2012, 2013 Simon N. Wood
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

## Reparameterization trick as Wood (2004,2006). 
## fooling lmer using Fabian Scheipl's trick (now adapted for new lme4).
## Hopefully unnecessary now but check that mcgv::gamm.setup amounts to the same thing
gamm4.setup<-function(G)
## set up the model matrix, penalty matrices and auxilliary information about the smoothing bases
## needed for a gamm4 fit.
## There is an implicit assumption that any rank deficient penalty does not penalize 
## the constant term in a basis. 
## 1. Calls gam.setup, as for a gam to produce object G suitable for estimating a gam.
## 2. Works through smooth list, G$smooth, modifying so that... 
##    i) Smooths are reparameterized to have a sequence of (portion of) identity matrix
##       penalties. 
##    ii) 'random' list is accumulated containing random effect model matrices for terms.     
##    iii) Sparse version of full model matrix in original parameterization is also accumulated
##    iv) Various indices are created for moving between the parameterizations.
{ 

  if (!is.null(G$L)) 
    stop("stan_gamm4 can not handle linked smoothing parameters (probably from use of `id' or adaptive smooths)")
  # now perform re-parameterization...

  first.f.para <- G$nsdf+1 
  first.r.para <- 1
 
  random <- list()
  
  if (G$nsdf > 0) ind <- 1:G$nsdf else ind <- rep(0,0)  
  X <- G$X[,ind,drop = FALSE] # accumulate fixed effects into here

  xlab <- rep("", 0)
  
  G$Xf <- as(X,"dgCMatrix") ## sparse version of full matrix, treating smooths as fixed

  first.para <- G$nsdf + 1

  used.names <- names(data) ## keep track of all variable names already used

  if (G$m) for (i in 1:G$m) { ## work through the smooths
    sm <- G$smooth[[i]]
    sm$X <- G$X[,sm$first.para:sm$last.para, drop = FALSE]
    ## convert smooth to random effect and fixed effects
    ## Can this ::: be avoided?
    rasm <- mgcv:::smooth2random(sm,used.names,type=2)
    used.names <- c(used.names,names(rasm$rand))    

    sm$fixed <- rasm$fixed

    ## deal with creation of sparse full model matrix  
    if (!is.null(sm$fac)) { 
      flev <- levels(sm$fac) ## grouping factor for smooth
      n.lev <- length(flev)
      for (k in 1:n.lev) {
        G$Xf <- cbind2(G$Xf,as(sm$X * as.numeric(sm$fac==flev[k]),"dgCMatrix"))
      }
    } else { 
      n.lev <- 1
      G$Xf <- cbind2(G$Xf,as(sm$X,"dgCMatrix"))
    }

    ## now append random effects to main list
    n.para <- 0 ## count random coefficients
    rinc <- rind <- rep(0,0)
    if (!sm$fixed) {
      for (k in 1:length(rasm$rand)) n.para <- n.para + ncol(rasm$rand[[k]])
      sm$lmer.name <- names(rasm$rand)
      random <- c(random,rasm$rand)
      sm$trans.D <- rasm$trans.D
      sm$trans.U <- rasm$trans.U ## matrix mapping fit coefs back to original
    }

    ## ensure stored first and last para relate to G$Xf in expanded version

    sm$last.para <- first.para + ncol(rasm$Xf) + n.para - 1
    sm$first.para <- first.para
    first.para <- sm$last.para + 1    

    if (ncol(rasm$Xf)) {
      Xfnames <- rep("",ncol(rasm$Xf)) 
      k <- length(xlab) + 1
      for (j in 1:ncol(rasm$Xf)) {
        xlab[k] <- Xfnames[j] <-
          new.name(paste(sm$label,"Fx",j,sep=""),xlab)
        k <- k + 1
      } 
      colnames(rasm$Xf) <- Xfnames
    }

    X <- cbind(X, rasm$Xf) # add fixed model matrix to overall fixed X
   
    sm$first.f.para <- first.f.para
    first.f.para <- first.f.para + ncol(rasm$Xf)
    sm$last.f.para <- first.f.para - 1 ## note less than sm$first.f.para => no fixed

    ## store indices of random parameters in smooth specific array
    sm$rind <- rasm$rind #- 1 + first.r.para
    sm$rinc <- rasm$rinc 

    sm$pen.ind <- rasm$pen.ind ## pen.ind==i TRUE for coef penalized by ith penalty

    sm$n.para <- n.para
 
    sm$X <- NULL ## delete model matrix
  
    G$smooth[[i]] <- sm  ## replace smooth object with extended version 
  }
 
  G$random <- random ## named list of random effect matrices
  G$X <- X  ## fixed effects model matrix

  return(G)
} ## end of gamm4.setup

# Routine to fit a GAMM to some data. Fixed and smooth terms are defined in the formula, but the wiggly 
# parts of the smooth terms are treated as random effects. The onesided formula random defines additional 
# random terms.
stan_gamm4 <- function(formula, random = NULL, family = gaussian(), data = list(), 
                       weights = NULL, subset = NULL, na.action, knots = NULL, 
                       drop.unused.levels = TRUE, ..., 
                       prior = normal(), prior_intercept = normal(),
                       prior_ops = prior_options(),
                       prior_covariance = decov(), prior_PD = FALSE, 
                       algorithm = c("sampling", "meanfield", "fullrank"), 
                       adapt_delta = NULL, QR = FALSE) {
  if (!is.null(random)) {
    if (!inherits(random,"formula")) stop("stan_gamm4 requires `random' to be a formula")
    random.vars <- all.vars(random)
    random_list <- lme4::findbars(random)
    names(random_list) <- names(lme4::mkReTrms(random_list, 
                                               data, drop.unused.levels)$cnms)
    for (i in seq_along(random_list)) { 
      random_list[[i]] <- as.formula(paste("~", random_list[i]))
    }
  } 
  else random_list <- random.vars <- NULL

  # create model frame.....
  gp <- mgcv::interpret.gam(formula) # interpret the formula 
  # call mgcv::gamm() with zero iterations to set everything (?) up
  # calling gamm4::gamm4() is too slow because we cannot do zero iterations
  settings <- list(maxIter = 0, msMaxIter = 0, niterEM = 0, msMaxEval = 0, 
                   gradHess = FALSE, apVar = FALSE, opt = "optim", 
                   optimMethod = "BFGS", natural = FALSE)
  # G <- mgcv::gamm(formula, random_list, correlation = NULL, family, 
  #                 data, weights, subset, 
  #                 na.action, knots, control = settings, niterPQL = 0,
  #                 drop.unused.levels = drop.unused.levels)$gam
  call <- match.call(expand.dots = TRUE)
  mc <- match.call(expand.dots = FALSE)
  mc[[1]] <- quote(mgcv::gamm)
  mc$... <- mc$prior <- mc$prior_intercept <- mc$prior_ops <- 
    mc$prior_covariance <- mc$prior_PD <- mc$algorithm <- mc$adapt_delta <- NULL
  mc$control <- settings
  mc$random <- random_list
  G <- eval(mc, parent.frame(1L))$gam
  mf <- G$model
  gmf <- eval(mf, parent.frame()) # the model frame now contains all the data, for the gam part only 
  gam.terms <- attr(gmf,"terms") # terms object for `gam' part of fit -- need this for prediction to work properly
  X <- model.matrix(gam.terms, data)[,-1,drop = FALSE]
  if (length(random.vars)) {
    mc$formula <- as.formula(paste(paste(deparse(gp$fake.formula,
            backtick = TRUE), collapse = ""), "+", paste(random.vars,
            collapse = "+")))
    mf <- eval(mc, parent.frame())
  }
  else mf <- gmf
  rm(gmf)

  Terms <- attr(mf,"terms")
  
  ## summarize the *raw* input variables
  ## note can't use get_all_vars here -- buggy with matrices
  vars <- all.vars(gp$fake.formula[-2]) ## drop response here
  inp <- parse(text = paste("list(", paste(vars, collapse = ","),")"))
  dl <- eval(inp, data, parent.frame())
  names(dl) <- vars ## list of all variables needed
  var.summary <- G$var.summary # FIXME: verify this is adequate

  ## lmer offset handling work around...
  ## variables not in mf raw -- can cause lmer problem
  mvars <- vars[!(vars %in% names(mf))]
  if (length(mvars) > 0) ## append raw versions to mf
    for (i in 1:length(mvars)) mf[[mvars[i]]] <- dl[[mvars[i]]]

  rm(dl) ## save space 

  # pmf$formula <- gp$pf
  # pmf <- eval(pmf, parent.frame()) # pmf contains all data for non-smooth part 
  # pTerms <- attr(pmf,"terms")

  linear <- FALSE
  pTerms <- G$pTerms
  
  n.sr <- length(G$random) # number of random smooths (i.e. s(...,fx=FALSE,...) terms)

  if (is.null(random) && n.sr == 0) 
    stop("stan_gamm4 models must have at least 1 smooth with unknown smoothing parameter or at least one other random effect")

  offset.name <- attr(mf,"names")[attr(attr(mf,"terms"),"offset")]

  yname <- mgcv::new.name("y",names(mf))
  eval(parse(text=paste("mf$",yname,"<-G$y",sep=""))) # what does this do?
  Xname <- new.name("X",names(mf))
  eval(parse(text=paste("mf$",Xname,"<-G$X",sep=""))) # what does this do?
    
  lme4.formula <- paste(yname,"~",Xname,"-1")
  if (length(offset.name)) 
    lme4.formula <- paste(lme4.formula,"+",offset.name) 

  ## Basic trick is to call (g)lFormula to set up model, with simple i.i.d. dummy random effects for the 
  ## penalized component of each smooth. This results in columns of Z being produced for these dummy's,
  ## which can be over-written with the right thing. NOTE: that lambdat could also be modified, I think!!

  ## Add the random effect dummy variables for the smooth
  r.name <- names(G$random) 
  if (n.sr) for (i in 1:n.sr) { # adding the constructed variables to the model frame avoiding name duplication
    mf[[r.name[i]]] <- factor(rep(1:ncol(G$random[[i]]), length=nrow(G$random[[i]])))
    lme4.formula <- paste(lme4.formula,"+ (1|",r.name[i],")")
  }
  
  if (!is.null(random)) ## append the regular random effects
    lme4.formula <- paste(lme4.formula, "+" , substring(deparse(random), first = 2))
  
  lme4.formula <- as.formula(lme4.formula)


  b <- glFormula(lme4.formula, data = data, family = family, weights = G$w,
                 control = glmerControl(check.nlev.gtreq.5 = "ignore",
                                        check.nlev.gtr.1 = "stop",
                                        check.nobs.vs.rankZ = "ignore",
                                        check.nobs.vs.nlev = "ignore",
                                        check.nobs.vs.nRE = "ignore"))
  y <- b$fr[, as.character(b$formula[2L])]
  if (is.matrix(y) && ncol(y) == 1L) y <- as.vector(y)
  offset <- model.offset(glmod$fr) %ORifNULL% double(0)
  weights <- validate_weights(weights)
  if (is.null(prior))
    prior <- list()
  if (is.null(prior_intercept)) 
    prior_intercept <- list()
  if (!length(prior_ops)) 
    prior_ops <- list(scaled = FALSE, prior_scale_for_dispersion = Inf)
  group <- b$reTrms
  group$decov <- prior_covariance
  algorithm <- match.arg(algorithm)
  
  if (n.sr) { ## Fabian Scheipl's trick of overwriting dummy slots revised for new structure
     tn <- names(b$reTrms$cnms) ## names associated with columns of Z (same order as Gp)
     ind <- 1:length(tn)
     sn <- names(G$random) ## names of smooth random components
     for (i in 1:n.sr) { ## loop through random effect smooths
       k <- ind[sn[i]==tn] ## which term should contain G$random[[i]] 
       ii <- (b$reTrms$Gp[k]+1):b$reTrms$Gp[k+1]
       b$reTrms$Zt[ii,] <- as(t(G$random[[i]]),"dgCMatrix")
       b$reTrms$cnms[[k]] <- attr(G$random[[i]],"s.label") 
     }
  }
  rm(b)
  
  ## now do the actual fitting...
  ret <- list()
  stanfit <- stan_glm.fit(x = X, y = y, weights = weights,
                          offset = offset, family = family,
                          prior = prior, prior_intercept = prior_intercept,
                          prior_ops = prior_ops, prior_PD = prior_PD, 
                          algorithm = algorithm, adapt_delta = adapt_delta,
                          group = group, QR = QR, ...)

  Z <- pad_reTrms(Z = t(group$Zt), cnms = group$cnms, 
                  flist = group$flist)$Z
  colnames(Z) <- b_names(names(stanfit), value = TRUE)
  fit <- nlist(stanfit, family, formula, offset, weights, 
               x = if (getRversion() < "3.2.0") cBind(X, Z) else cbind2(X, Z), 
               y = y, data, call, terms = NULL, model = NULL, 
               prior.info = get_prior_info(call, formals()),
               na.action, contrasts, algorithm, glmod)
  out <- stanreg(fit)
  class(out) <- c(class(out), "lmerMod")
  
  ## now fake a gam object, replace / add elements to this as needed
  object <- G
  object$y <- y
  pvars <- all.vars(delete.response(object$terms))
  object$pred.formula <- if (length(pvars) > 0) reformulate(pvars) else NULL

  ## to unpack coefficients look at names(ret$lme$flist), ret$lme@Zt, ranef(), fixef()
 
    ## let the GAM coefficients in the original parameterization be beta,
    ## and let them be beta' in the fitting parameterization. 
    ## Then beta = B beta'. B and B^{-1} can be efficiently accumulated
    ## and are useful for stable computation of the covariance matrix
    ## etc... 
  
    B <- diag(1, ncol(G$Xf), ncol(G$Xf))
    Xfp <- G$Xf
    ## Transform  parameters back to the original space....
    bf <- as.numeric(fixef(out)) ## the fixed effects
    br <- ranef(out) ## a named list
    if (G$nsdf) p <- bf[1:G$nsdf] else p <- array(0,0) ## fixed parametric component
    if (G$m > 0) for (i in 1:G$m) {
      fx <- G$smooth[[i]]$fixed 
      first <- G$smooth[[i]]$first.f.para; last <- G$smooth[[i]]$last.f.para
      if (first <= last) beta <- bf[first:last] else beta <- array(0,0)
      if (fx) b <- beta else { ## not fixed so need to undo transform of random effects etc. 
        b <- rep(0,0)
        for (k in 1:length(G$smooth[[i]]$lmer.name)) ## collect all coefs associated with this smooth
          b <- c(b,as.numeric(br[[G$smooth[[i]]$lmer.name[k]]][[1]]))     
        b <- b[G$smooth[[i]]$rind] ## make sure coefs are in order expected by smooth
        b <- c(b,beta) 
        b <- G$smooth[[i]]$trans.D*b
        if (!is.null(G$smooth[[i]]$trans.U)) b <- G$smooth[[i]]$trans.U%*%b ## transform back to original 
      }
      p <- c(p,b)
     
      ## now fill in B...
      ind <- G$smooth[[i]]$first.para:G$smooth[[i]]$last.para
      if (!fx) { 
         D <- G$smooth[[i]]$trans.D
         if (is.null(G$smooth[[i]]$trans.U)) B[ind,ind] <- Diagonal(length(D),D) else
         B[ind,ind] <- t(D*t(G$smooth[[i]]$trans.U))
      }
      ## and finally transform G$Xf into fitting parameterization...
      Xfp[,ind] <- G$Xf[,ind,drop = FALSE] %*% B[ind, ind, drop = FALSE]

    }
 
    object$coefficients <- p

    ## need to drop smooths from Zt and then
    ## form Z'phiZ + I \sigma^2

    vr <- VarCorr(out) ## list of ranef variance components in the same order as Gp
    
    scale <- as.numeric(attr(vr,"sc"))^2 ## get the scale parameter
    if (!is.finite(scale) || scale == 1) { ## NOTE: better test???
      scale <- 1
      object$scale.estimated <- FALSE
    } 
    else object$scale.estimated <- TRUE
    
    sp <- rep(-1, n.sr)

    Zt <- Matrix(0, 0, ncol(group$Zt))
    if (n.sr == 0) sn <- NULL ## names by which smooths are known in mer
    rn <- names(vr)
    ind <- rep(0, 0) ## index the non-smooth random effects among the random effects
    for (i in 1:length(vr)) {
      if (is.null(sn) || !(rn[i] %in% sn)) { ## append non smooth r.e.s to Zt
        Gp <- group$Gp ## group index ends
        ind <- c(ind,(Gp[i]+1):Gp[i+1])
      } 
      else if (!is.null(sn)) { ## extract smoothing parameters for smooth r.e.s
        k <- (1:n.sr)[rn[i] == sn] ## where in original smooth ordering is current smoothing param
        if (as.numeric(vr[[i]] > 0)) sp[k] <- scale / as.numeric(vr[[i]]) else 
          sp[k] <- 1e10
      }
    }

    if (length(ind)) { ## extract columns corresponding to non-smooth r.e.s 
      Zt <- group$Zt[ind,] ## extracting random effects model matrix
      root.phi <- group$Lambdat[ind,ind] ## and corresponding sqrt of cov matrix (phi)
    }

    object$prior.weights <- G$w
    if (linear) {
      object$weights <- object$prior.weights 
      V <- Matrix::Diagonal(n = length(object$weights), x =scale / object$weights) 
    } 
    else {
     # mu <- getME(ret$mer,"mu")
     # eta <- family$linkfun(mu)
      # object$weights <- ret$mer@resp$sqrtWrkWt()^2
      object$weights <- object$prior.weights
      ## object$prior.weights*family$mu.eta(eta)^2/family$variance(mu)
      V <- Matrix::Diagonal(x = scale / object$weights)
      #V <- Diagonal(x=scale*family$variance(mu)/object$prior.weights)
    }

  
    if (nrow(Zt) > 0) 
      V <- V + crossprod(root.phi %*% Zt) * scale ## data or pseudodata cov matrix, treating smooths as fixed now

    ## NOTE: Cholesky probably better in the following - then pivoting 
    ##       automatic when solving....

    R <- Matrix::chol(V, pivot = TRUE); piv <- attr(R,"pivot") 

    G$Xf <- as(G$Xf,"dgCMatrix")
    Xfp <- as(Xfp,"dgCMatrix")
    
    if (is.null(piv)) {
      WX <- as(solve(t(R), Xfp),"matrix")    ## V^{-.5}Xp -- fit parameterization
      XVX <- as(solve(t(R), G$Xf),"matrix")  ## same in original parameterization 
    } 
    else {
      WX <- as(solve(t(R), Xfp[piv,]),"matrix")    ## V^{-.5}Xp -- fit parameterization
      XVX <- as(solve(t(R), G$Xf[piv,]),"matrix")  ## same in original parameterization
    }
    qrz <- qr(XVX, LAPACK = TRUE)
    object$R <- qr.R(qrz); object$R[,qrz$pivot] <- object$R

    XVX <- crossprod(object$R) ## X'V^{-1}X original parameterization

    object$sp <- sp
    
    colx <- ncol(G$Xf)
    Sp <- matrix(0, colx, colx) # penalty matrix - fit param
    first <- G$nsdf + 1
    k <- 1
    if (G$m > 0) for (i in 1:G$m) { # Accumulate the total penalty matrix
      if (!object$smooth[[i]]$fixed) {
        ii <- object$smooth[[i]]$first.para:object$smooth[[i]]$last.para ## index this smooth's params
        for (j in 1:length(object$smooth[[i]]$S)) { ## work through penalty list
          ind <- ii[object$smooth[[i]]$pen.ind == j] ## index of currently penalized
          diag(Sp)[ind] <-  sqrt(object$sp[k]) ## diagonal penalty
          k <- k + 1
        }									              
      }
      first <- last + 1 
    }
   
    ## Alternative cov matrix calculation. Basic
    ## idea is that cov matrix is computed stably in
    ## fitting parameterization, and then transformed to
    ## original parameterization. 
    qrx <- qr(rbind(WX, Sp / sqrt(scale)), LAPACK = TRUE)
    Ri <- backsolve(qr.R(qrx), diag(ncol(WX)))
    ind <- qrx$pivot;ind[ind] <- 1:length(ind) ## qrx$pivot
    Ri <- Ri[ind,] ## unpivoted square root of cov matrix in fitting parameterization Ri Ri' = cov
    Vb <- B %*% Ri; Vb <- tcrossprod(Vb)

    object$edf <- rowSums(Vb * t(XVX))
   
    object$df.residual <- length(object$y) - sum(object$edf)

    object$sig2 <- scale
    if (linear) { object$method <- "lmer.REML"
    } else { object$method <- "glmer.ML"}

    object$Vp <- as(Vb, "matrix")
  
    object$Ve <- as(Vb %*% XVX %*% Vb, "matrix")
   
    class(object) <- "gam"
   
    ## Restore original smooth list, if it was split to deal with t2 terms...
    if (!is.null(G$original.smooth)) {
      object$smooth <- G$smooth <- G$original.smooth
    }

    ## If prediction parameterization differs from fit parameterization, transform now...
    ## (important for t2 smooths, where fit constraint is not good for component wise 
    ##  prediction s.e.s)

    if (!is.null(G$P)) {
      object$coefficients <- G$P %*% object$coefficients
      object$Vp <- G$P %*% object$Vp %*% t(G$P) 
      object$Ve <- G$P %*% object$Ve %*% t(G$P) 
    }

    object$linear.predictors <- predict.gam(object,type="link")
    object$fitted.values <- object$family$linkinv(object$linear.predictors)
    
    object$residuals <- residuals(ret$mer) 

    if (G$nsdf>0) term.names <- colnames(G$X)[1:G$nsdf] else term.names <- array("", 0)
    n.smooth<-length(G$smooth) 
    if (n.smooth)
    for (i in 1:n.smooth) {
      k <-1
      for (j in object$smooth[[i]]$first.para:object$smooth[[i]]$last.para) {
        term.names[j]<-paste(object$smooth[[i]]$label,".",as.character(k),sep="")
        k <- k + 1
      }
    }
    names(object$coefficients) <- term.names  # note - won't work on matrices!!
    names(object$edf) <- term.names
    names(object$sp) <- names(G$sp)

    object$gcv.ubre <- if (isREML(ret$mer)) REMLcrit(ret$mer) else deviance(ret$mer)

    if (!is.null(G$Xcentre)) object$Xcentre <- G$Xcentre ## any column centering applied to smooths

    out$gam <- object
    return(out)
} ## end of stan_gamm4

