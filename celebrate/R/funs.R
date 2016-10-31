## https://www.freesound.org/people/primordiality/sounds/78823/
## https://www.freesound.org/people/kirbydx/sounds/175409/

.onLoad <- function(libname,pkgname) {
  options(celebrate.sounds=list(
    success=system.file("sounds",
      "78822__primordiality__fanfare-1.wav",
      package="celebrate"),
    failure=system.file("sounds",
                        "175409__kirbydx__wah-wah-sad-trombone.wav",
                        package="celebrate")))
    load_sounds()
}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage("WARNING: this package is augmenting the printout ",
                          "of statistical model summaries with audio output. ",
                          "Prepare to be annoyed and/or embarrassed ..."
                          )
}


#' Load sound files
#'
#' @param sound_files list of names of WAV files
#' @importFrom audio load.wave
#' @export
load_sounds <- function(sound_files=getOption("celebrate.sounds")) { 
  sound_list  <- lapply(sound_files,
                        audio::load.wave)
  assign("sound_list",sound_list,environment(print.summary.lm))
}

##' @importFrom stats coef
##' @importFrom audio play
##' @export
print.summary.lm <- function(x,...) {
    ## CRAN will never be happy about this, but then again
    ##   this package will never go on CRAN.  Don't know of
    ##   a way to do this without ::: ...
    stats:::print.summary.lm(x,...)
    cc <- coef(x)
    cc <- cc[rownames(cc)!="(Intercept)",,drop=FALSE]
    pvals <- cc[,"Pr(>|t|)"]
    if (any(pvals<0.05)) {
        audio::play(sound_list$success)
    } else {
        audio::play(sound_list$failure)
    }
    return(invisible(x))
}

