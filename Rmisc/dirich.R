library(gtools)
#' @param n dimensionality
#' @param N total number of values
#' @param nsim number of simulations (-> n*nsim values)
simfun <- function(n,N,nsim=1000,scale=TRUE,verbose=TRUE) {
    if (!missing(N)) {
        nsim <- round(N/n)
    }
    if (verbose) cat("n:",n,"\n")
    r <- gtools::rdirichlet(nsim,alpha=rep(1,n))
    H <- rowSums(r^2)
    if (!scale) return(H)
    return(sqrt(n/2)*(H*(n/2)-1))
}

dd <- simfun(n=10,N=1e6)

sumfun <- function(dd, npt=50, xr=10,
                   brk=seq(-xr,xr,length=npt+1),
                   trim=TRUE) {
    bad_ind <- dd<=(-xr) | dd>=(xr)
    bad <- dd[bad_ind]
    if (length(bad)>0) {
        cat("bad values:\n")
        print(bad)
        if (trim) dd <- dd[!bad_ind]
    }
    h <- hist(dd, breaks=brk, plot=FALSE)
    d <- density(dd,from=-xr,to=xr,n=npt)
    data.frame(mids=h$mids,
               histdens=h$density,
               densx=d$x,
               dens=d$y)
}

ff <- function(n,N=5e6) sumfun(simfun(n, N=N))
## compute mean homozygosity; wasteful but ...
ff2 <- function(n,N=1e6) data.frame(m=mean(simfun(n, N=N, scale=FALSE)))

nvec <- unique(round(exp(seq(log(2),log(1000),length=100))))
fn <- "dirich_data.rds"
Hfn <- "dirich_meanH.rds"
if (!file.exists(fn)) {
    set.seed(101)
    res <- plyr::ldply(setNames(nvec, nvec), .fun=ff,.progress="text", .id="n")
    res$n <- as.numeric(as.character(res$n))
    saveRDS(res,file=fn)
    res2 <- plyr::ldply(setNames(nvec, nvec), .fun=ff2,.progress="text", .id="n")
    res2$n <- as.numeric(as.character(res2$n))
    saveRDS(res2,file=Hfn)
} else {
    res <- readRDS(fn)
    res2 <- readRDS(Hfn)
}

## basic plot
pfun <- function(dd,...) {
    hist(dd,col="gray",breaks=100,
         freq=FALSE,main="",...)
    lines(density(dd),col="red")
}

pfun(dd)

library(animation)

## fancy plot
pfun2 <- function(nval,dd=subset(res,n==nval),
                  ymax=0.5,xr=10) {
    dx <- diff(dd$mids)[1]
    brks <- c(dd$mids[1]-dx/2,dd$mids+dx/2)
    h <- list(breaks=brks,mids=dd$mids,density=dd$histdens,
              xname="",equidist=TRUE)
    class(h) <- "histogram"
    par(las=1)
    plot(h,freq=FALSE,col="gray",ylim=c(0,ymax),xlim=c(-xr,xr),
         main=sprintf("Dirichlet homozygosity with n=%d",nval))
    lines(dd$densx,dd$dens,col="red")
}

pfun2(2,ymax=2)
par(las=1,bty="l")
for (i in unique(res$n)) {
    pfun2(i,ymax=1)
    box()
}

## clean up
system("rm -Rf js css images dirich_homozyg.html")
system("ssh ms.mcmaster.ca 'cd public_html/misc; rm -Rf js css images dirich_homozyg.html'")

library(animation)
opts <- ani.options(interval=0.2)
saveHTML({
    for (i in unique(res$n)) pfun2(i,ymax=1)
},
img.name="dirich_homozyg",
htmlfile="dirich_homozyg.html")

system("scp -r js css images dirich_homozyg.html ms.mcmaster.ca:~/public_html/misc")

## TODO:
##  - shiny/sliders?
##  - add Gaussian curve?
##  - Q-Q version?
##  - compute moments etc.?



library(ggplot2); theme_set(theme_bw())
dd2 <- transform(data.frame(n=res2$n),
                 ninv=1/n,
                 ninv2=2/n)
ggplot(res2,aes(n,m))+geom_point()+
    scale_x_log10()+scale_y_log10()+
    geom_line(data=dd2,aes(y=ninv),colour="blue")+
    geom_line(data=dd2,aes(y=ninv2),colour="red")

