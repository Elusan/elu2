#' Fit ELU2 Model
#'
#' Fits the customized ELU2 stock assessment model using TMB.
#'
#' @param inp List of input variables as output by check.inp.
#' @param verbose Should detailed outputs be provided (default: TRUE).
#' @param dbg Debugging option. Will print out runtime information useful for debugging if set to 1. Will print even more if set to 2.
#' @return A result report containing estimates of model parameters, random effects (biomass and fishing mortality), reference points (Fmsy, Bmsy, MSY) including uncertainties given as standard deviations.
#' @export
#' @examples
#' data(pol)
#' rep <- fit.spict(pol$albacore)
#' Bmsy <- get.par('logBmsy', rep, exp=TRUE)
#' summary(rep)
#' plot(rep)
#' @import TMB
fit.elu2 <- function(inp, verbose=TRUE, dbg=0){
  rep <- NULL
  # Check input list
  inp <- check.inp(inp, verbose = verbose)
  datin <- make.datin(inp, dbg)
  pl <- inp$parlist
  tic <- Sys.time()
  # Cycle through phases
  for (i in 1:inp$nphases){
    if (inp$nphases > 1) cat(paste('Estimating - phase', i, '\n'))
    # Create TMB object
    obj <- make.obj(datin, pl, inp, phase=i)
    if (dbg<1){
      # Do estimation
      if (inp$optimiser == 'nlminb'){
        opt <- try(nlminb(obj$par, obj$fn, obj$gr, control=inp$optimiser.control))
        if (class(opt)!='try-error'){
          pl <- obj$env$parList(opt$par)
        }
      }
      if (inp$optimiser == 'optim'){
        if (inp$optim.method == 'SANN'){
          tmpgr <- obj$gr
          obj$gr <- NULL # SANN doesn't use gradient
        }
        opt <- try(optim(par=obj$par, fn=obj$fn, gr=obj$gr, method=inp$optim.method,
                         control=inp$optimiser.control))
        if (inp$optim.method == 'SANN' & class(opt) != 'try-error'){
          obj$gr <- tmpgr # Restore gradient function
          cat('SANN optimisation done, switching to BFGS to refine...\n')
          opt <- try(optim(par=obj$par, fn=obj$fn, gr=obj$gr, method='BFGS',
                           control=inp$optimiser.control))
          #obj$fn(opt$par)
          #obj$gr(opt$par)
        }
        if (class(opt) != 'try-error'){
          opt$objective <- opt$value
          pl <- obj$env$parList(opt$par)
        }

      }
    }
  }
  if (dbg<1){
    optfailflag <- class(opt)=='try-error'
    sdfailflag <- FALSE
    if (optfailflag){ # Optimisation failed
      cat('obj$par:\n')
      print(obj$par)
      cat('obj$fn:\n')
      print(obj$fn())
      cat('obj$gr:\n')
      print(obj$gr())
      stop('Could not fit model. Error msg:', opt)
    } else {
      if (inp$do.sd.report){
        # Calculate SD report
        rep <- try(TMB::sdreport(obj,
                                 getJointPrecision=inp$getJointPrecision,
                                 bias.correct=inp$bias.correct,
                                 bias.correct.control=inp$bias.correct.control,
                                 getReportCovariance = inp$getReportCovariance))
        sdfailflag <- class(rep) == 'try-error'
        if (sdfailflag){
          warning('Could not calculate sdreport.\n')
          rep <- NULL
        }
      }
      if (is.null(rep)){ # If sdreport failed or was not calculated
        rep <- list()
        if (sdfailflag){
          rep <- list()
          rep$sderr <- 1
          rep$par.fixed <- opt$par
          rep$cov.fixed <- matrix(NA, length(opt$par), length(opt$par))
        }
      }
      rep$inp <- inp
      rep$obj <- obj
      rep$opt <- opt
      rep$opt$gr <- rep$obj$gr(rep$opt$par)
      rep$pl <- obj$env$parList(opt$par)
      obj$fn()
      rep$Cp <- obj$report()$Cp
      rep$report <- obj$report()
      if (!sdfailflag & inp$reportall){
        #  - Calculate Prager's statistics -
        #rep <- calc.prager.stats(rep)
        # - Built-in OSAR -
        if (!inp$osar.method == 'none'){
          reposar <- try(calc.osa.resid(rep))
          if (class(reposar) != 'try-error'){
            rep <- reposar
          }
        }
      }
    }
  }
  toc <- Sys.time()
  if (!is.null(rep)){
    rep$computing.time <- as.numeric(toc - tic)
    class(rep) <- "spictcls"
  }
  return(rep)
}



calc.prager.stats <- function(rep){
  if (!'stats' %in% names(rep)){
    rep$stats <- list()
  }
  K <- get.par('logK', rep, exp=TRUE)[2]
  Bests <- get.par('logB', rep, exp=TRUE)[rep$inp$indest, 2]
  Bmsy <- get.par('logBmsy', rep, exp=TRUE)[2]
  if (!any(is.na(Bests)) & !is.na(Bmsy)){
    Bdiff <- Bmsy - Bests
    # Prager's nearness
    if (any(diff(sign(Bdiff))!=0)){
      rep$stats$nearness <- 1
    } else {
      rep$stats$nearness <- 1 - min(abs(Bdiff))/Bmsy
    }
    # Prager's coverage
    rep$stats$coverage <- min(c(2, (min(c(K, max(Bests))) - min(Bests))/Bmsy))
  }
  return(rep)
}

#' @name make.datin
#' @title Create data list used as input to TMB::MakeADFun.
#' @param inp List of input variables as output by check.inp.
#' @param dbg Debugging option. Will print out runtime information useful for debugging if set to 1.
#' @return List to be used as data input to TMB::MakeADFun.
#' @export
make.datin <- function(inp, dbg=0){
  datin <- list(reportall=as.numeric(inp$reportall),
                reportRel=as.numeric(inp$reportRel),
                dt=inp$dt,
                dtpredcinds=inp$dtpredcinds,
                dtpredcnsteps=inp$dtpredcnsteps,
                dtprediind=inp$dtprediind,
                dtpredeinds=inp$dtpredeinds,
                dtpredensteps=inp$dtpredensteps,
                indlastobs=inp$indlastobs,
                obssrt=inp$obssrt,
                stdevfacc=inp$stdevfacC,
                stdevfaci=inp$stdevfacIin,
                stdevface=inp$stdevfacE,
                isc=inp$isc,
                isi=inp$isi,
                ise=inp$ise,
                nobsC=inp$nobsC,
                nobsI=sum(inp$nobsI),
                nobsE=inp$nobsE,
                ic=inp$ic,
                nc=inp$nc,
                ie=inp$ie,
                ne=inp$ne,
                ii=inp$iiin,
                iq=inp$iqin,
                isdi=inp$isdiin,
                ir=inp$ir,
                isdf=inp$isdf,
                logmcov=inp$logmcovariatein,
                seasons=inp$seasons,
                seasonindex=inp$seasonindex,
                nseasons=inp$nseasons,
                seasonindex2=inp$seasonindex2,
                splinemat=inp$splinemat,
                splinematfine=inp$splinematfine,
                omega=inp$omega,
                seasontype=inp$seasontype,
                efforttype=inp$efforttype,
                timevaryinggrowth=as.numeric(inp$timevaryinggrowth),
                logmcovflag=as.numeric(inp$logmcovflag),
                ffacvec=inp$ffacvec,
                fconvec=inp$fconvec,
                indpred=inp$indpred,
                robflagc=inp$robflagc,
                robflagi=inp$robflagi,
                robflage=inp$robflage,
                stochmsy=ifelse(inp$msytype=='s', 1, 0),
                stabilise=inp$stabilise,
                MSYregime=inp$MSYregime,
                iuse=as.numeric(inp$iuse),
                residFlag=as.numeric(inp$residFlag),

                priorn=inp$priors$logn,
                priorngamma=inp$priors$logngamma,
                priorr=inp$priors$logr,
                priorK=inp$priors$logK,
                priorm=inp$priors$logm,
                priormu=inp$priors$mu,
                priorq=inp$matrixpriors$logq,
                priorqf=inp$priors$logqf,
                priorbkfrac=inp$priors$logbkfrac,
                priorsdb=inp$priors$logsdb,
                priorsdm=inp$priors$logsdm,
                priorsdf=inp$priors$logsdf,
                priorsdi=inp$matrixpriors$logsdi,
                priorsde=inp$priors$logsde,
                priorsdc=inp$priors$logsdc,
                prioralpha=inp$priors$logalpha,
                priorbeta=inp$priors$logbeta,
                priorpsi=inp$priors$logpsi,
                priorB=inp$priors$logB,
                priorF=inp$priors$logF,
                priorBBmsy=inp$priors$logBBmsy,
                priorFFmsy=inp$priors$logFFmsy,
                priorBmsyB0=inp$priors$BmsyB0,

                simple=inp$simple,
                reportmode=inp$reportmode,
                simRandomEffects=inp$sim.random.effects,
                dbg=dbg)
  return(datin)
}

#' @name make.obj
#' @title Create TMB obj using TMB::MakeADFun and squelch screen printing.
#' @param datin Data list.
#' @param pl Parameter list.
#' @param inp List of input variables as output by check.inp.
#' @param phase Estimation phase, integer.
#' @return List to be used as data input to TMB.
#' @export
#' @import TMB
make.obj <- function(datin, pl, inp, phase=1){
  obj <- TMB::MakeADFun(data=datin, parameters=pl, random=inp$RE, DLL=inp$scriptname,
                        hessian=TRUE, map=inp$map[[phase]])
  TMB:::config(trace.optimize=0, DLL=inp$scriptname)
  verbose <- FALSE
  obj$env$tracemgc <- verbose
  obj$env$inner.control$trace <- verbose
  obj$env$silent <- ! verbose
  obj$fn(obj$par)
  return(obj)
}
