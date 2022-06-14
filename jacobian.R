## https://twitter.com/BjornstadOttar/status/1535773594252353537
Jacobian <- function (.vars, ...) {
    ## capture the elements
    vf <- substitute(list(...))[-1L]
    ## capture
    vars <- sapply(substitute(.vars), deparse)
    if (length(vars)>1) vars <- vars[-1L]
    ## differentiatiate each of the components
    sapply(
        vars,
        \(var) sapply(vf,D,name=var)
   ) -> jac
   ## set up dimensions for the output
   dd <- c(length(vf),length(vars))
   dim(jac) <- NULL
   ## figure out names
   dn <- list(
      ifelse(
        nzchar(names(vf)),
        names(vf),
        sapply(vf,deparse)
    ),
    vars
   )
  ## the hard part
  fun <- function (...) {
    ## figure out the correct environment in which to evaluate everything  
    e <- c(as.list(sys.frame(sys.nframe())),...)
    ## evaluate each element of the jacobian
    J <- vapply(jac,eval,numeric(1L),envir=e)
    ## reset the dimensions and names  
    dim(J) <- dd
    dimnames(J) <- dn
    J
  }
  ## !!
  formals(fun) <- eval(
    parse(text=paste0("alist(",paste0(c(vars,"..."),"=",collapse=","),")"))
  )
  ## return function as the result
  fun
}

f <- Jacobian(c(x,y),sin(x),cos(x),atan(y/x),tan(x+y))
f(y=2,x=3)

## cOde

library(cOde)
jacobianSymb(c(A="A*B", B="A+B"))

library(calculus)

## also: numDeriv, Deriv, madness, TMB (pytorch thing), (Nash autodiff thing?)

## symbolic and numerical?
## user-extended derivs table?
## autodiff?


