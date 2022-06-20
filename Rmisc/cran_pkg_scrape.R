pkg_info <- function(pkg = "numDeriv") {
    pkg_url <- paste0("https://CRAN.R-project.org/package=", pkg)
    r <- readLines(pkg_url)
    latest_rel <- (
        r[grep("Published", r)[1] + 1]
        |> gsub(pattern = "</?td>", replacement = "")
    )
    arch_url <- paste0("https://cran.r-project.org/src/contrib/Archive/", pkg)
    r <- readLines(arch_url)
    re_date <- "[0-9]{4}-[0-9]{2}-[0-9]{2}"
    arch_dates <- (grep(sprintf("[^_]%s[^.]", re_date), r, value = TRUE)
        |> stringr::str_extract(re_date)
        |> as.Date()
    )
    n_rel <- length(arch_dates) + 1
    n_rdep <- length(tools::package_dependencies(pkg, reverse = TRUE)[[pkg]])
    return(tibble::lst(n_rdep, n_rel, first_rel = arch_dates[1], latest_rel))
}
    
