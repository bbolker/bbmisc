## see also GillespieSSA package on CRAN:
##  description at http://www.jstatsoft.org/v25/i12/paper

## TO DO:
##  extended error checking by checking (negative) values of trans
##   vs rate gradient??
##  documentation, examples?
##  use built-in progress bar?
##  implement SSA modes?
## compilation? NIMBLE ???

## FIXME:
##  

##' generic Gillespie driver
##' @param start starting value of state vectors (length ns)
##' @param ratefun function to return vector of rates (length nr): arguments
##' @param trans transition matrix (nr x ns)
##' @param pars named parameter vector
##' @param times time vector for output
##' @param progress show progress bar?
##' @param debug print debugging output?
gillesp <- function(start,ratefun,trans,pars,
                    times=0:20,progress=FALSE,debug=FALSE) {
  if (!is.matrix(trans)) trans <- matrix(trans)
  ## check dimensions etc.
  t0 <- times[1]                 ## set time to starting time
  ntimes <- length(times)
  nstate <- length(start)
  statenames <- names(start)
  transnames <- rownames(trans)
  X <- start                     ## set state to starting state
  rates <- ratefun(X,pars,t0) ## calculate current rates
  ntrans <- length(rates)
  ##
  if (nrow(trans)!=ntrans)
    stop("mismatch in number of rates between rate function and trans")
  if (!is.null(transnames) && !is.null(names(rates))) {
    if (length(setdiff(transnames,names(rates)))>0 ||
        length(setdiff(names(rates),transnames))>0)
      stop("rate name mismatch between rate function and trans")
    if (!all(transnames==names(rates)))
      stop("rate name order mismatch between rate function and trans")
  }
  if (ncol(trans)!=nstate)
    stop("mismatch in number of states between start and trans")
  if (!is.null(statenames) && !is.null(colnames(trans))) {
    if (length(setdiff(statenames,colnames(trans)))>0 ||
        length(setdiff(colnames(trans),statenames))>0)
      stop("state name mismatch between start and trans")
    if (!all(statenames==colnames(trans)))
      stop("rate name order mismatch between start and trans")
  }
  ## should use built-in progress bar code instead ...
  if (progress) pctr <- 0
  res <- matrix(nrow=length(times),ncol=nstate,  ## matrix for results
                dimnames=list(times,names(start)))
  transmat <- matrix(nrow=length(times),ncol=ntrans,
                     dimnames=list(times,names(rates)))
  for (ctr in 1:(ntimes-1)) {     ## loop over reporting times
    transctr <- rep(0,ntrans)
    res[ctr,] <- X                ## record current state
    while (t0<times[ctr+1]) {
      rates <- ratefun(X,pars,t0) ## calculate current rates
      if (any(rates<0)) {
        msg <-  paste("(rate(s)",
                   paste(which(rates<0),collapse=","))
        if (!is.null(rownames(trans))) {
          msg <- paste(msg,"=",paste(rownames(trans)[rates<0],collapse=","))
        }
        msg <- paste(msg,                   "=",
                   paste(rates[rates<0],collapse=","),")")
        stop("negative rates computed ",msg)
      }           
      if (all(rates==0)) break    ## extinction?
      totrate <- sum(rates)       
      elapsed <- rexp(1,totrate)  ## sample elapsed time
      which.trans <- sample((1:ntrans)[rates>0],
                            size=1,prob=rates[rates>0])
                                        # pick transition
      transctr[which.trans] <- transctr[which.trans]+1
      t0 <- t0+elapsed            ## update time
      X <- X+trans[which.trans,]  ## add transition values to current state
    } ## while loop
    transmat[ctr,] <- transctr
    ## progress bar -- 60 #
    while (progress && ctr/length(times)>pctr/60) {
      cat("#"); pctr <- pctr+1
    }
  }
  if (progress) cat("\n")
  as.data.frame(cbind(times,res,transmat))
}


##  FIXME: also check for all rates==0?
##' generic discrete simulation driver
##' @inheritParams gillesp
##' @param dt time step
##' @param savetrans ??
discsim <- function(start,ratefun,trans,pars,times=0:200,
                    dt=1,savetrans=FALSE) {
    ## number of time steps in each period
    nt <- diff(times)/dt
    if (!all(abs(round(nt)-nt) < 1e-3))
      stop("dt not commensurate with time steps")
    nt <- round(nt)
    ntrans = nrow(trans)
    res = matrix(nrow=length(times),ncol=length(start))
    colnames(res) <- names(start)
    if (savetrans) savetransmat <- matrix(0,nrow=length(times)-1,
                                          ncol=nrow(trans))
    res[1,] = start
    binom = apply(trans==-1,1,any)
    bpos = apply(trans[binom,]<0,1,which)
    delta = numeric(ntrans)
    for (i in 2:length(times)) {
        tmpres = res[i-1,]
        for (j in 1:(nt[i-1])) {
            rates = ratefun(tmpres,pars,times[i-1]+(j-1)*dt)
            ## ttype = rowSums(trans)
            ## ttype = -1: loss (binomial)
            ##       =  0: transfer (binomial)
            ##       =  1: gain (Poisson)
            bN = tmpres[bpos]
            ## how should binomial rates be calculated, really?
            ## if we have an overall rate RATE, we want the
            ## expected number of transitions to be RATE per
            ## unit time, so the probability should be RATE/bN
            ## -- but also need to let this saturate towards
            ## 1.0 appropriately
            bprob = 1-exp(-(rates[binom]/bN*dt))
            delta[binom] = rbinom(sum(binom),
                   prob=ifelse(bN==0,0,bprob),
                   size=bN)
            delta[!binom] = rpois(sum(!binom),rates[!binom]*dt)
            ## hack: need to keep tmpres as a vector
            tmpres = drop(tmpres+ t(trans) %*% delta) ## flipped order
            if (savetrans) savetransmat[i-1,] <- savetransmat[i-1,]+delta
        }
        res[i,] <- tmpres
        if (i<nrow(res) && sum(res[i,c(-1,-ncol(res))])==0) {
            ## extinct
            ## warning(paste("extinction at time",times[i]))
            res[(i+1):nrow(res),] = rep(res[i,],each=nrow(res)-i)
            if(savetrans) savetransmat[i:nrow(savetransmat),] =
              rep(0,(nrow(savetransmat)-i+1)*ncol(savetransmat))
            break
        }
    }
    ## FIXME
    ## colnames(res) <- statenames.eSEIR(pars["nE"],pars["nI"])
    ## if (savetrans) colnames(transmat) = transnames.eSEIR(pars["nE"],pars["nI"])
    d <- data.frame(times,res)
    if (savetrans) {
      d <- data.frame(d,rbind(savetransmat,rep(NA,ncol(trans))))
    }
    d
}  

##############

##' state transitions etc. for logistic function
glogistfun <- function(state,params,time) {
  with(c(as.list(state),as.list(params)),
       c(birth=f*N,
         death=(mu+alpha*N)*N))
}
##' @rdname glogistfun
gtrans <- matrix(c(1,-1),ncol=1,
                 dimnames=list(c("birth","death"),"N"))

##' state transitions etc. for SIR model
sirfun <- function(state,params,time) {
  with(c(as.list(state),as.list(params)),
       c(infection=beta*S*I/N,
         recovery =gamma*I))
}

sirtrans <- matrix(c(-1,1,0,
                     0,-1,1),byrow=TRUE,nrow=2,
                   dimnames=list(c("infection","recovery"),
                     c("S","I","R")))

##' SIR with birth/death
bdsirfun <- function(state,params,time) {
  with(c(as.list(state),as.list(params)),
       c(infection=beta*S*I/(S+I+R),
         recovery=gamma*I,
         birth=mu*(S+I+R),
         Sdeath=mu*S,
         Ideath=mu*I,
         Rdeath=mu*R))
}

bdsirtrans <- matrix(c(-1, 1, 0,
                      0,-1, 1,
                      1, 0, 0,
                     -1, 0, 0,
                      0, -1,0,
                      0,  0,-1),
                     byrow=TRUE,ncol=3,
                     dimnames=list(c("infection","recovery",
                       "birth",
                       "Sdeath",
                       "Ideath",
                       "Rdeath"),
                     c("S","I","R")))


## EXAMPLES
if (FALSE) {
  ppars <- c(beta=2,gamma=1,mu=0.1)
  t1 <- system.time(g1 <-   gillesp(start=c(S=400,I=10,R=590),
                                  pars=ppars,
                                  time=1:1000,
                                  ratefun=bdsirfun,
                                  trans=bdsirtrans,
                                  progress=TRUE))
  matplot(g1[c("S","I","R")],type="l",lty=1,las=1,
          ylab="density",bty="l",log="y",col=c(1,2,4))
  with(as.list(ppars),
       abline(h=gamma/beta*1000,col=4,lty=2))

  ## will take about 20 minutes to run
  zz <- replicate(20,
                  gillesp(start=c(S=400,I=10,R=590),
                          pars=c(beta=2,gamma=1,mu=0.1),
                          time=1:1000,
                          ratefun=bdsirfun,
                          trans=bdsirtrans),simplify=FALSE)
  matplot(zz[[1]][c("S","I","R")],type="n",las=1,
          ylab="density",bty="l",log="y",ylim=c(1,2000))
  invisible(lapply(zz,
                   function(x) {
                     matlines(x[c("S","I","R")],
                              lty=1,col=c(1,2,4))}))
  
}
