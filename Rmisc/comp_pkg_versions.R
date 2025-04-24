## compare  sessionInfo() information from two sources

library(tidyverse)

## for example
comp_txt <- readLines("outputs/tidyMRE.Rout")
new_txt <- readLines("outputs/tidyMRE_BMB.Rout")

## locate all lines in text that look like print-outs of
## vectors (i.e. lines starting with [number])
## split into chunks by vector
split_vecs <- function(txt, pat = "^ *\\[[0-9]+") {
    LL <- grep(pat, txt)
    ## find contiguous sequences
    split_indices <- c(0,cumsum(diff(LL) > 1))
    ss <- split(LL, split_indices)
    vals <- lapply(ss, \(i) txt[i])
    names(vals) <- vapply(ss, \(s) txt[s[1]-1], character(1))
    return(vals)
}

## extract a tibble {package_name, version_number}
## from output of split_vecs()
##' @param txt character vector representing the output
##' of `sessionInfo()`
##' (e.g. via `capture.output(sessionInfo())` or from a `.Rout` file)
##' @return a tibble with columns `pkg` (package name), `version` (version number)
get_versions <- function(txt) {
    pkgs <- (split_vecs(txt)[c("other attached packages:",
                               "loaded via a namespace (and not attached):")]
        |> unlist()
        |> unname()
        |> gsub(pattern = "\\[[0-9]+\\]", replacement = "")
    )
    pkgs2 <- scan(textConnection(pkgs), what = character(1), quiet = TRUE)
    vers_df <- (tibble(vals = pkgs2)
        |> tidyr::separate(vals, sep = "_",
                           into=c("pkg", "version"))
    )
    return(vers_df)
}

## run get_versions
## compare ...
(list(old = comp_txt, new = new_txt)
    |> map_dfr(get_versions, .id = "w")
    |> pivot_wider(names_from = w, values_from = version)
    |> filter(is.na(old) | is.na(new) | old != new)
)

