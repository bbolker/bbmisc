## gradients and utility functions for SEIR with Erlang E and I periods
## version 1.1
## version 1.2: BUG FIX in erlangseirmod, various cleanup; numerical integrators
##   and discretizers for serial interval
## version 1.3: serial interval stuff continued
## version 1.4: serial interval stuff, discrete simulation
##    (discrsim tested, gillesp is not)
##    wvec_to_wmat
##  1.4.2: tweak to wvec/wmat
##  1.4.3: much work on death-death interval, aggregated discrete simulations
##  1.4.4: DE: more precise incidence calculation

##  1.5.0: BMB 2008/03/12: various tweaks, working on stoch simulators
##                         added "unbounded" flag
##  1.5.1: BMB 2008/03/13: changed discsim to operate (always) on absolute
##        rates, rather than per capita: ratefun.pcerlSEIR now obsolete
##  1.5.2: BMB 2008/03/13: added dt (time step) to discsim; added
##                       recovery output (sim. death)
##  1.6: BMB 2008/03/14: various fixes, added ddesolve
##  1.7: BMB 2008/05/06: discsim fix (nt[i-1]); added est.lambda
##                         added times=times to discsim call
##                         calc stats for stoch models as well
##  1.8: BMB 2008/05/14: more discsim fixes
##  1.9: BMB 2010/07/08: fixed retval/match.arg bugs; moved to PBSddesolve;
## ?? BMB 2018/11/23: merge with gillesp.R

## TO DO:
##   check results vs. WL2007 expressions for lambda/R relationship
##   documentation/packaging
##   change interface to put everything
##     (nE,nI,unbounded etc.) into params
##  "ring buffer" version of discsim?
##  "queue" version of gillespie?

if (!require(deSolve)) cat("please install the deSolve package\n")
source("gillesp.R")

## e.g:
## set parameters and output times
## default params: R0=1.8
params1 = c(mu=0, N=1, beta=0.36, gamma=1/5, sigma=1/3)

## times  = 0:200

use.betafun <- FALSE

est.lambda <- function(x,I0,N,times,minSprop=0.98) {
    suppressWarnings(lambda <- try(coef(lm(log(I)~times,
                                       subset=I>I0 & S>minSprop*N,
                                       data=x))[2]))
    if (class(lambda)=="try-error") NA else lambda
}

## est lambda for discrete data with only time-series of I available
est.lambda_2 <- function(x,min=5,max=20,minval=0.001) {
    times <- seq(along=x)
    suppressWarnings(lambda <- coef(lm(log(x+minval)~times,
                                       subset=times>min & times<max))[2])
    lambda
}

## RUN AN ERLANG SIMULATION
erlang.sim = function(nI,nE,params, I0=1e-6,times,
  retval=c("incid","raw","recovery","stats","all"),
  type=c("ODE","dstoch","gillespie"),
  aggreg=FALSE,
  use.betafun=FALSE,
  period=7,
  unbounded=FALSE,
  ...) {
    type <- match.arg(type)
    retval <- match.arg(retval)
    ## set initial values for all subcompartments from given initial state
    EIvec = c(E=rep(0,ifelse(nE==Inf,1,nE)),
      I=c(I0,rep(0,ifelse(nI==Inf,0,nI-1))))
    start = params["N"]*c(S=1-I0,EIvec,R=0)
    ## append shape parameters to params vector (being sure to name them)
    params = c(params,nI=nI,nE=nE,unbounded=unbounded)
    param_names <<- names(params) ### UGH! HACK!
    y_names <<- names(start)
    ## make param names usable, e.g., beta = params["beta"]
    with(as.list(params), {
    if (type=="ODE") {
        ## solve the ODEs for times
        ## add state variables for incidence/recovery
        start = c(start,cuminf=0,cumonset=0,cumrecovery=0)
        y_names <<- names(start)
        ## (approx. = death)
        ## browser()
        use_dde = (nE==Inf)
        if (use_dde) {
            if (!require(PBSddesolve)) stop("please install the PBSddesolve package")
            if (nI!=Inf) stop("can only do discrete lags for both E and I")
            gfun = ddeseirmod
            times = with(as.list(params),c(-1/sigma-0.0001,times))
            ## start["E"] <- start["I"]
            ## start["I"] <- 0 
            intfun = dde
        } else {
            intfun = lsoda
            if (nE==0) { gfun = erlangsirmod
                     } else gfun = erlangseirmod
        }
        x = as.data.frame(intfun(y=start, func=gfun, times=times, parms=params))
        if (use_dde) {
            x = x[-1,]
            times = times[-1]
        }
        ## sum all E and I compartments into combined E, I values:
        ## I think this works generically even if
        ##  nI==1 or nE==1 but am not sure
        ## browser()
        x = condense(x,nI=nI,nE=nE) 
        finalsize = x$R[nrow(x)]
        if (is.na(finalsize)) finalsize = x$R[nrow(x)-1]  ## minor kluge for Gillespie
        ## approx. estimate of lambda (?)
        lambda <- est.lambda(x,I0,N,times)
        inc = c(diff(x$cumonset),0); ## diff is one shorter, but need same length
        recovery = c(diff(x$cumrecovery),0); ## diff is one shorter, but need same length
        if (aggreg) inc = aggreg(inc,period)
        stats <- c(finalsize,lambda)
    } else if (type=="dstoch") { ## discrete stochastic sim
        start = round(start)
        if (sum(start[1+nE+(1:nI)])==0) {
            stop("zero infectives: rounding problem?")
        }
        names(start) = statenames.eSEIR(nE,nI)
        x = discsim(start,ratefun.eSEIR,
          transmat.eSEIR(nE,nI,unbounded),
          pars=params,savetrans=TRUE,times=times,...)
        x = condense(x,nI=nI,nE=nE)
        inc = x$onset[-nrow(x)]  ## drop final NA
        recovery = x$recovery[-nrow(x)]  ## drop final NA
        if (aggreg) {
            inc = aggreg(inc,period)
            recovery = aggreg(recovery,period)
        }
        if (nrow(x)<20) {
            finalsize <- NA
            warning("not computing final size with < 20 points")
        } else finalsize = mean(x$R[(nrow(x)-20):nrow(x)])
        lambda <- est.lambda(x,I0,N,times)
        stats = c(finalsize,lambda)
    } else if (type=="gillespie") {
        start = round(start)
        if (sum(start[1+nE+(1:nI)])==0) {
            stop("zero infectives: rounding problem?")
        }
        names(start) = statenames.eSEIR(nE,nI)
        x = gillesp(start,ratefun.eSEIR,
          transmat.eSEIR(nE,nI),
          pars=params,times=times,...)
        ## browser()
        inc = x$onset[-nrow(x)]  ## drop final NA
        recovery = x$recovery[-nrow(x)]  ## drop final NA
        x = condense(x,nI,nE)
        if (aggreg) {
            inc = aggreg(inc,period)
            recovery = aggreg(recovery,period)
        }
        if (nrow(x)<20) {
            finalsize <- NA
            warning("not computing final size with < 20 points")
        } else finalsize = mean(x$R[(nrow(x)-20):nrow(x)])
        lambda <- est.lambda(x,I0,N,times)
        stats = c(finalsize,lambda)
    }
    else stop("unknown simulation type")
    return(switch(retval,
                  incid=inc,
                  recovery=recovery,
                  stats=stats,
                  raw=x,
                  all=list(incid=inc,stats=stats,x=x)))
  }) ## end with
}

## GRADIENT OF ERLANG SEIR MODEL
## with() is similar to attach(), i.e., makes the names of the params
## available.
erlangseirmod = function(t, x, params) {
  with(c(as.list(params),as.list(x)), {
    S = x[1]
    E = x[2:(nE+1)]
    I = x[(nE+2):(nE+nI+1)]
    R = x[nE+nI+2]
    sigma2=sigma*nE ## transition rate through E boxes
    gamma2=gamma*nI ## transition rate through I boxes
    #### HACK TO CREATE BEHAVIOURAL CHANGES: really need beta as vector
    if(use.betafun) {
      infrate = betafun(t,beta)*S*sum(I)/N
    } else {
      infrate = beta*S*sum(I)/N
    }
    dS = mu * (N  - S)
    if (!unbounded) dS = dS - infrate
    ## handle exposed
    if (nE>1) {
      dE = sigma2*(c(0,E[1:(nE-1)])-E)
    } else dE=-E*sigma
    dE[1] = dE[1] + infrate
    dE = dE - mu*E
    ## handle infected
    if (nI>1) {
      dI = gamma2*(c(0,I[1:(nI-1)])-I)
    } else dI=-I*gamma
    onset = sigma2*E[nE]
    dI[1] = dI[1]+onset
    dI = dI - mu*I
    recovery = gamma2*I[nI]
    dR = recovery - mu * R
    other = c(incid=infrate,death=gamma2*I[nI])
    ## dummy variable to integrate incidence so differencing C
    ## (cumulative incidence) will yield incidence between time steps
    res=c(dS, dE, dI, dR, infrate, onset, recovery)
    list(res,other)
  })
}

## GRADIENT OF ERLANG SEIR MODEL
## with() is similar to attach(), i.e., makes the names of the params
## available.
ddeseirmod = function(t, x, params) {
    names(x) <- y_names ### UGH! HACK!
    names(params) <- param_names ### UGH! HACK!
    with(c(as.list(params),as.list(x)), {
        if(use.betafun) {
            infrate = betafun(t,beta)*S*I/N
        } else {
            infrate = beta*S*I/N
        }
        infrate = betafun(t,beta)*S*I/N
        ggrad <- function(t,m) {
            z <- if (t<m) rep(0,length(x)) else pastgradient(t-m)
            names(z) <- y_names
            z
        }
        onset = ggrad(t,1/sigma)["cuminf"]
        recovery = ggrad(t,1/gamma)["cumonset"]
        dS = mu * (N  - S)
        if (!unbounded) dS = dS - infrate
        dE=infrate-mu*E-onset
        dI=onset-recovery-mu*I
        dR = recovery-mu*R
        res=c(dS, dE, dI, dR, infrate, onset, recovery)
        list(res,NULL)
    })
}

## make beta a function of t
## FIX HACK: this is not the best approach...
##           esp the method of shutting off the function...
betafun = function(t,beta) {
    if (TRUE) {
        fac = ifelse (t < 20, 0.5,
          ifelse( t < 40, 1,
                 ifelse( t < 60, 0.5, 0.2 )
                 )
          )
    } else {
        ##junling's example:
        fac = ifelse (t < 40, 1, 0.8 )
    }
    return(if (use.betafun)fac*beta else (beta))
}

## gradients for SIR with Erlang I periods
erlangsirmod = function(t, x, params) {
  with(c(as.list(params)), { ## don't include as.list(x)
    S = x[1]
    I = x[2:(nI+1)]
    R = x[nI+2]
    gamma2=gamma*nI
    #### HACK TO CREATE BEHAVIOURAL CHANGES: really need beta as vector
    if(use.betafun) {
        infrate = betafun(t,beta)*S*sum(I)/N
    } else {
        infrate = beta*S*sum(I)/N
    }
    dS = mu * (N  - S)
    if (!unbounded) dS = dS - infrate
    ## handle infected
    if (nI>1) {
      dI = gamma2*(c(0,I[1:(nI-1)])-I)
    } else dI=-I*gamma
    dI[1] = dI[1]+ infrate
    dI = dI - mu*I
    recovery = gamma2*I[nI]
    onset = infrate ## equal for SIR
    dR = recovery - mu * R
    ## dummy variables to integrate rates
    res=c(dS, dI, dR, infrate, onset, recovery)
    list(res)
  })
}

## collapse an n-period output into a single E and I variable.
## drop=FALSE keeps the output as a matrix even if it has only
## one column, which is important for applying rowSums to it
## elsewhere (so cases n=1 work)
condense <- function(x,nI,nE) {
    ## don't really need to condense for delay
    if (nE==Inf) nE=1
    if (nI==Inf) nI=1
    if (nE>0) {
    data.frame(time=x[,1],
               S=x[,2],
               E=rowSums(x[, 2+(1:nE) ,drop=FALSE]),
               I=rowSums(x[, ((nE+2)+(1:nI)) ,drop=FALSE]),
               R=x[,nE+nI+3],
               x[,-(1:(nE+nI+3))])
  } else {
      data.frame(time=x[,1],
                 S=x[,2],
                 I=rowSums(x[, 2+(1:nI), drop=FALSE]),
                 R=x[,nE+nI+3],
                 x[,-(1:(nE+nI+3))])
  }
}


## compute (approximate) incidence
## correct incidence is obtained after integrating the equations
## via dummy incidence variable, C[t]-C[t-1].
incidence <- function(x,beta,N) {
  ###beta*x$S*x$I/N
  inc = betafun(x$time,beta)*x$S*x$I/N
  #cat("incidence: dumping\n")
  #print(beta)
  #print(x)
  #print(N)
  #print(inc)
  #stop("incidence: finished dump\n")
  return(inc)
}

aggreg = function(x,period=7) {
  uneven = length(x) %% period
  if (uneven != 0) {  ## pad with zeros
    x = c(x,rep(0,period-uneven))
  }
  ## re-form into a matrix (by column)
  ## and take column sums to aggregate
  x = colSums(matrix(x,nrow=period))
  x
}

## OBSOLETE
## pull incidence out of discrete daily simulation
## incidence_d = function(x,aggreg=FALSE,period=7) {
##   inc = x$onset[-nrow(x)] ## drop final NA 
## }

## absolute rates for individual transitions in the
##  Erlang-SEIR model (for Gillespie or discsim
ratefun.eSEIR <- function(X,pars,time,percap=FALSE)  {
  vals <- c(as.list(pars),as.list(X))   ## attach state and pars as lists
  rates <- with(vals, {
    E = X[1+dseq(nE)]
    I = X[nE+1+dseq(nI)]
    infection=beta*S*sum(I)/N
    sigma2=sigma*nE ## transition rate through E boxes
    gamma2=gamma*nI ## transition rate through I boxes
    Etrans= if (nE>1) sigma2*E[1:(nE-1)] else numeric(0)
    onset = if (nE>0) sigma2*E[nE]       else numeric(0)
    Itrans= if (nI>1) gamma2*I[1:(nI-1)] else numeric(0)
    recovery = gamma2*I[nI]
    birth=mu*N
    death = mu*X
    res=c(infection,Etrans,onset,Itrans,recovery,birth,death)
    names(res) = transnames.eSEIR(nE,nI)
    res
  })
}

## per capita rates for individual transitions in the
##  Erlang-SEIR model (except birth, which is absolute)
## ratefun.pcerlSEIR <- function(X,pars,time,percap=FALSE)  {
##   vals <- c(as.list(pars),as.list(X))   ## attach state and pars as lists
##   rates <- with(vals, {
##     E = X[1+dseq(nE)]
##     I = X[nE+1+dseq(nI)]
##     infection=beta*sum(I)/N
##     sigma2=sigma*nE ## transition rate through E boxes
##     gamma2=gamma*nI ## transition rate through I boxes
##     Etrans= if (nE>1) rep(sigma2,nE-1) else numeric(0)
##     onset = if (nE>0) sigma2           else numeric(0)
##     Itrans= if (nI>1) rep(gamma2,nI-1) else numeric(0)
##     recovery = gamma2
##     birth=mu*N
##     death = rep(mu,length(X))
##     res=c(infection,Etrans,onset,Itrans,recovery,birth,death)
##     names(res) = transnames.eSEIR(nE,nI)
##     res
##   })
## }


## directed sequence: returns numeric(0) for decreasing
##   sequences
dseq <- function(a,b) {
  if (missing(b)) { b=a; a=1 }
  if (b<a) numeric(0) else a:b
}

## paste function  that returns empty
##  string if any arguments are length zero
dpaste <- function(...,sep=" ") {
  L = list(...)
  if (min(sapply(L,length))==0) {
    character(0)
  } else paste(...,sep=sep)
}
  
statenames.eSEIR <- function(nE,nI) {
  c("S",dpaste("E",dseq(nE),sep=""),
    dpaste("I",1:nI,sep=""),"R")
}
transnames.eSEIR <- function(nE,nI) {
  c("infection",
    dpaste("Etrans",dseq(nE-1),dseq(nE-1)+1,sep="."),
    if (nE>0) "onset" else character(0),
    dpaste("Itrans",dseq(nI-1),dseq(nI-1)+1,sep="."),
    "recovery",
    "birth",
    paste("death",statenames.eSEIR(nE,nI),sep="."))
}

transmat.eSEIR <- function(nE,nI,unbounded=FALSE) {
    if (nE>0) {
        Etransmat <- -diag(nE)[-nE,,drop=FALSE]
        Etransmat[col(Etransmat)==row(Etransmat)+1] <- 1
        Etransmat <- cbind(rep(0,nE-1),Etransmat,matrix(0,nrow=nE-1,ncol=nI+1))
    }
    else Etransmat = matrix(nrow=0,ncol=nE+nI+2)
    if (nI>0) {
        Itransmat <- -diag(nI)[-nI,,drop=FALSE]
        Itransmat[col(Itransmat)==row(Itransmat)+1] <- 1
        Itransmat <- cbind(matrix(0,nrow=nI-1,ncol=nE+1),Itransmat,rep(0,nI-1))
    }
    else Itransmat = matrix(nrow=0,ncol=nE+nI+2)
    Evals <- if (nE==0) numeric(0) else rep(0,nE-1)
    inftrans <- if (unbounded) c(0,1) else c(-1,1)
    trans.eSEIR <- rbind(c(inftrans,Evals,rep(0,nI),0), ## infection
                         Etransmat,  ## E-trans
                         ## onset: redundant -1/+1 to denote S/R at beg,end
                         if (nE>0) c(rep(0,1+nE-1),-1,1,rep(0,nI-1+1)) else {
                             matrix(nrow=0,ncol=nE+nI+2)
                         },
                         Itransmat,  ## I-trans
                         c(rep(0,1+nE+nI-1),-1,1), ## recovery
                         c(1,rep(0,nE+nI+1)), ## birth
                         -diag(nE+nI+2))  ## death
    dimnames(trans.eSEIR) <- list(transnames.eSEIR(nE,nI),
                                  statenames.eSEIR(nE,nI))
    trans.eSEIR
}



#####################
require(emdbook)
finalsize = function(R0) {
   1+1/R0*lambertW(-R0*exp(-R0))
}

## SERIAL INTERVAL CALCULATORS

## calculate the probability that the serial interval is in [tstart,tend]
## given latent distribution function and infectious upper cumulative distribution
## VERY VERY VERY SLOW!!!
serial_intprob_gen = function(tstart,tend,latdistrib,infuppercum) {
  a1 = adapt(2,c(tstart,0),c(tend,tend),
    f=function(v) {
      t=v[1]
      tau=v[2]
      ifelse(tau>t,0,
             latdistrib(tau)*infuppercum(t-tau))
    })
    a1$value
}

## calculate the probability density of the serial interval at t
## given latent distribution function and infectious upper cumulative distribution
## 1D numeric convolution:
##  prob(still infectious at time t) = 
##   int(tau<t, p(latent=tau)*p(infectious period >= t-tau))
serial_dens_gen = function(t,latdistrib,infuppercum) {
  a1 = integrate(
    f=function(tau) {
      latdistrib(tau)*infuppercum(t-tau)
    },
    lower=0,upper=t)
    a1$value
}
#
## calculate the probability density of the serial interval at t
## given mean and shape parameters for gamma-distributed latent and infectious distributions
serial_dens_gamma = function(t,lat.mean,lat.shape,inf.mean,inf.shape) {
  latdistrib=function(t) dgamma(t,shape=lat.shape,scale=lat.mean/lat.shape)
  infuppercum=function(t) pgamma(t,shape=inf.shape,scale=inf.mean/inf.shape,lower.tail=FALSE)
  serial_dens_gen(t,latdistrib,infuppercum)
}

## calculate the probability density of the DEATH-serial interval at t
## given mean and shape parameters for gamma-distributed latent,
##   death, and infectious distributions
serial_dens_deathgamma = function(t,lat.mean,lat.shape,inf.mean,inf.shape) {
  latdistrib=function(t) dgamma(t,shape=lat.shape,scale=lat.mean/lat.shape)
  infuppercum=function(t) pgamma(t,shape=inf.shape,scale=inf.mean/inf.shape,lower.tail=FALSE)
  deathdistrib=function(t) dgamma(t,shape=death.shape,scale=death.mean/death.shape)
  s1 = serial_dens_gen(t,latdistrib,infuppercum)
}


## calculate the probability that the serial interval is in [tstart,tend]
## given mean and shape parameters for gamma-distributed latent and infectious distributions
serial_intprob_gamma = function(tstart,tend,lat.mean,lat.shape,inf.mean,inf.shape) {
  latdistrib=function(t) dgamma(t,shape=lat.shape,scale=lat.mean/lat.shape)
  infuppercum=function(t) pgamma(t,shape=inf.shape,scale=inf.mean/inf.shape,lower.tail=FALSE)
  serial_intprob_gen(tstart,tend,latdistrib,infuppercum)
}

## SIR model: calculate discretized serial interval, which is just the
## normalized upper cumulative tail of the appropriate gamma distribution
serial_intcurve_gamma_SIR = function(breaks,inf.mean,inf.shape) {
  pvals = pgamma(breaks,shape=inf.shape,scale=inf.mean/inf.shape,lower.tail=FALSE)
  ser = diff(pvals) ## calculate slice areas
  ser/sum(ser) ## normalize
}

## calculate discretized serial interval for
## specified breakpoints; then normalize
serial_intcurve_gamma = function(breaks,...,verbose=FALSE) {
  ## code calling slow numerical integrator  -- don't
  ##     res = numeric(length(breaks)-1)
  ##     for (i in 1:(length(breaks)-1)) {
  ##       if (verbose) cat(i,breaks[i],breaks[i+1],"\n")
  ##       res[i] = serial_intprob_gamma(breaks[i],breaks[i+1],...)
  ##       if (verbose) cat("    ",res[i],"\n")
  ##     }
  s1 = serial_denscurve_gamma(breaks,...)
  ## trapezoid rule
  res = (s1[1:(length(s1)-1)]+s1[2:length(s1)])/2*diff(breaks)
  ## normalize
  res/sum(res)
}

## construct a matrix substituting values of x[ind]
## if ind<=0, substitute 0
posind = function(x,ind) {
  ind[ind<=0] <- NA
  y <- x[ind]
  y[is.na(y)] <- 0
  matrix(y,nrow=nrow(ind),ncol=ncol(ind))
}

## calculate the death-to-death interval for gamma-distributed latent, infectious,
##   and death periods
serial_denscurve_gamma_death = function(tvec,lat.mean,lat.shape,inf.mean,inf.shape,
  death.mean,death.shape,plot.it=FALSE,n=256) {
  ## figure out how to do this on basis of extreme values of inf, lat, death?
  tvec2 = seq(-3*death.mean,3*inf.mean+2*lat.mean+6*death.mean,length=n)
  dt = 1/diff(tvec[1:2])
  if (lat.mean>0) {
    latdistrib=function(t) dgamma(t,shape=lat.shape,scale=lat.mean/lat.shape)
    latvec = latdistrib(tvec2)
  }
  infuppercum=function(t) pgamma(t,shape=inf.shape,scale=inf.mean/inf.shape,lower.tail=FALSE) 
  infuvec = infuppercum(tvec2)
  infuvec = infuvec/(sum(infuvec)*dt)  ## normalize
  deathdistrib=function(t) dgamma(t,shape=death.shape,scale=death.mean/death.shape)
  deathvec = deathdistrib(tvec2)
  ## want to work with probability DENSITIES: sum(x)*dt should equal 1
  if (lat.mean>0) {  ## serial distribution: convolve latvec with upper-inf-cum
    ilvec = convolve(infuvec,rev(latvec),type="o")*dt
  } else ilvec=infuvec
  ## now compute death-death-serial convolution
  method = "brute"
  if (method=="brute") {
    v = numeric(length(tvec2))
    jkmat = outer(1:length(tvec2),1:length(tvec2),"-")  ## (j-k) -- integers
    deathmat = outer(deathvec,deathvec,"*")             ## D(j)*D(k)
    for (i in 1:length(tvec2)) {
      Smat = posind(ilvec,i+jkmat)                      ## S(i+j-k) [or 0 if i+j-k <= 0]
      v[i] = sum(Smat*deathmat)*dt^2                    ## integral/sum over (d1,d2) of triple product
    }
  }
  if (plot.it) {
    plot(tvec2,v,xlab="t",ylab="Death-to-death density")
  }
  v2 = approx(tvec2,v,tvec)$y
  sum(v2)/v2
}

  ## ild1vec = convolve(ilvec,rev(deathvec),type="o")
  ## ilvec = ilvec*diff(tvec2)[1]/sum(ilvec)

if (FALSE) {
s1 = serial_denscurve_gamma_death(0:100,lat.mean=2,lat.shape=100,inf.mean=2,inf.shape=100,
  death.mean=2,death.shape=100,n=100,plot=TRUE)
}

## MODIFY: take w vector and create a w matrix
wvec_to_wmat = function(tvec, wvec, times, digits=5) {
  ##  roundoff error would mess up match, so we round
  vdiff = length(times)-length(wvec)
  if (vdiff>0) wvec = c(wvec,rep(0,vdiff))
  ## cat("transposed w: don't forget to change your code if necessary\n")
  dist = t(round(outer(times,times,"-"),digits))
  ints = unique(sort(c(dist)))
  ints = ints[ints>=0]  ## drop negative times
  wmat = matrix(wvec[match(dist,ints)],nrow=length(times))
  ## all negative vals in dist become NAs because they don't match anything.
  ## will generalize to deaths data with negative serial intervals easily.
  wmat[is.na(wmat)] = 0
  wmat
}

## note that the tvec need not involve evenly spaced times.
w_numeric_gamma =  function( tvec, inf.mean, lat.mean, inf.shape, lat.shape, digits=5 ) {
  ## round to 5 digits
  dist = t(round(outer(tvec,tvec,"-"),digits))
  ## grab only that times at which we need to evaluate w(t):
  ints = unique(sort(c(dist)))
  ints = ints[ints>=0]  ## drop negative times
  ## construct breaks: first element (probably 0), followed
  ## by breaks halfway between the other elements, followed
  ## by a last point that puts the last element in an even interval
  breaks = c(ints[1],(ints[-length(ints)]+ints[-1])/2)
  ## 2nd entry above is list of averages of adjacent time pts, except 1st
  ## add an additional end time that is the length of the last time interval:
  breaks = c(breaks,breaks[length(breaks)]+(breaks[length(breaks)]-breaks[length(breaks)-1]))
  if (lat.mean==0 | lat.shape==0) {
    wvec = serial_intcurve_gamma_SIR(breaks,inf.mean, inf.shape)
  } else {
    wvec = serial_intcurve_gamma(breaks,lat.mean,lat.shape,inf.mean, inf.shape)
  }
  wvec
}

## calculate serial interval density for
## particular values
serial_denscurve_gamma = function(breaks,...) {
  res = sapply(breaks,
         serial_dens_gamma,...)
  res/sum(res)
}

## convolution of gamma(a,b) gamma(c,d)

## b^(-a) d^(-c) Exp[-y/d] y^(a+c-1) Hypergeometric1F1  (a,a+c,(-1/b+1/d)*y)

