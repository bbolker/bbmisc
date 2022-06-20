## run via
## export R_MAX_NUM_DLLS=1000
## R CMD BATCH find_identity.R
## see run_identity
## https://developer.r-project.org/Blog/public/2018/03/23/maximum-number-of-dlls/
is_identity <- function(f) {
    identical(deparse(body(f)), names(formals(f))[1])
}

## check
stopifnot(is_identity(identity))
stopifnot(is_identity(force))
stopifnot(!is_identity(lm))

## https://stackoverflow.com/questions/8696158/find-all-functions-including-private-in-a-package
get_identity_funs <- function(i) {
    cat(i,"\n")
    ns <- try(asNamespace(i),silent=TRUE)
    if (inherits(ns, "try-error")) return(NA)
    fun_nms <- c(lsf.str(ns))
    if (length(fun_nms)==0) return(character(0))
    funs <- lapply(fun_nms, get, envir=ns)
    return(fun_nms[sapply(funs,is_identity)])
}
stopifnot(identical(get_identity_funs("stats"),
                    c("as.dendrogram.dendrogram", "formula.formula",
                      "logLik.logLik", 
                      "na.pass", "offset", "terms.terms")))

## https://regex101.com/r/rK7qX5/1
## patterns of the form "junk.junk" or "as.junk.junk"
exclude_self_methods <- function(aa) {
    aa[!grepl("^(as\\.)?(\\w+)\\.\\2$",aa)]
}

pkgs <- c("base",getOption("defaultPackages"))
base_funs <- setNames(lapply(pkgs, get_identity_funs), pkgs)
exclude_self_methods(unlist(base_funs))

all_pkgs <- rownames(installed.packages())
all_funs <- setNames(lapply(all_pkgs, get_identity_funs),all_pkgs)
save.image("identity.rda")
              
L <- load("identity.rda")
length(all_pkgs) ## 
aa <- na.omit(unlist(all_funs))
length(aa) ## 225
aa2 <- exclude_self_methods(aa)
length(aa2)  ## 24

