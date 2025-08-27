library(lme4)
library(Matrix) ## for image()
data("sleepstudy", package = "lme4")
set.seed(101)
grp <- c("A", "B", "C")
ns <- length(levels(dd$Subject))
## make three copies of sleepstudy, tagged as different groups
dd <- lapply(1:3,
             function(i) data.frame(group = grp[i],
                                    sleepstudy)) |>
  ## add a bit of noise so cov matrices are slightly different
  transform(Reaction = rnorm(nrow(sleepstudy), Reaction, sd = 0.25),
            grp = factor(grp))

## ?modular
lf <- lFormula(Reaction ~ Days + (Days|grp:Subject), data = dd)
tt <- lf$reTrms

## update components appropriately
tt$theta <- rep(tt$theta, 3)
tt$lower <- rep(tt$lower, 3)
tt$Lind <- unlist(lapply(list(1:3, 4:6, 7:9), rep, times = ns))

## illustrate what this looks like (by substituting in Lind values)
M <- tt$Lambdat
M@x <- tt$Lind
image(M, useAbs=FALSE, useRaster = TRUE)

lf$reTrms <- tt
devfun <- do.call(mkLmerDevfun, lf)
opt <- optimizeLmer(devfun)
m <- mkMerMod(environment(devfun), opt, lf$reTrms, fr = lf$fr)
print(m)

## works except that VarCorr isn't constructed correctly
M <- tt$Lambdat[1:6, 1:6]
M@x <- m@optinfo$val
V <- crossprod(M) * sigma(m)^2
print(V, digits = 3)
sqrt(diag(V))
sapply(1:3,
       function(i) {
         inds <- (2*i-1):(2*i)
         V0 <- V[inds, inds]
         cov2cor(V0)[2,1]
       })
## not sure why blocks 2 and 3 are so close to singular
## try simulating with known (sensible) true values?
