#' Look up mathematicians by name
#'
#' Queries the Mathematics Genealogy Project's quick search (a simple
#' substring match against full names) and returns every match.
#'
#' @param name search string (e.g. a surname, or part of one)
#' @param timeout request timeout in seconds
#' @param retries number of retries on request failure (not on empty results)
#' @param quiet if FALSE, warn/message about failures and empty results
#' @return a data.frame with columns id, surname, given_name, school, year
#'   (zero rows if nothing matched or the request failed)
#' @export
mathgen_lookup <- function(name, timeout = 30, retries = 1, quiet = FALSE) {
    stopifnot(is.character(name), length(name) == 1, nzchar(name))
    resp <- mathgen_fetch("quickSearch.php", body = list(searchTerms = name),
                           method = "POST", timeout = timeout, retries = retries,
                           quiet = quiet)
    mathgen_parse_results(resp, label = sQuote(name), quiet = quiet)
}

#' Search the Math Genealogy Project by individual field
#'
#' Wraps the site's advanced search form
#' (https://genealogy.math.ndsu.nodak.edu/search.php), which lets each of
#' name, school, year, thesis title, country, and MSC subject code be
#' filtered independently. As with the quick search, each field does a
#' substring match rather than requiring an exact match. Leave a field as
#' "" to not filter on it; at least one field must be non-empty.
#'
#' @param given_name,other_names,family_name name fields
#' @param school degree-granting institution
#' @param year year the degree was granted (as it appears on the site, e.g. "1932")
#' @param thesis substring of the thesis title
#' @param country country of the school
#' @param msc two-digit top-level Mathematics Subject Classification code, e.g. "11" (Number theory)
#' @param chrono if TRUE, sort results chronologically instead of by name
#' @param timeout request timeout in seconds
#' @param retries number of retries on request failure (not on empty results)
#' @param quiet if FALSE, warn/message about failures and empty results
#' @return a data.frame with columns id, surname, given_name, school, year
#'   (zero rows if nothing matched or the request failed)
#' @export
mathgen_advanced_search <- function(given_name = "", other_names = "", family_name = "",
                                     school = "", year = "", thesis = "",
                                     country = "", msc = "", chrono = FALSE,
                                     timeout = 30, retries = 1, quiet = FALSE) {
    fields <- list(given_name = given_name, other_names = other_names,
                    family_name = family_name, school = school, year = year,
                    thesis = thesis, country = country, msc = msc)
    is_set <- vapply(fields, nzchar, logical(1))
    if (!any(is_set)) {
        stop("at least one of the search fields must be non-empty")
    }
    body <- c(fields, list(chrono = if (isTRUE(chrono)) "1" else "0", submit = "Submit"))
    resp <- mathgen_fetch("query-prep.php", body = body, method = "POST",
                           timeout = timeout, retries = retries, quiet = quiet)
    label <- paste(sprintf("%s=%s", names(fields)[is_set], sQuote(unlist(fields[is_set]))),
                    collapse = ", ")
    mathgen_parse_results(resp, label = label, quiet = quiet)
}
