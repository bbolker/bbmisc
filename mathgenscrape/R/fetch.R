mathgen_base_url <- "https://genealogy.math.ndsu.nodak.edu"

## Low-level fetch with a timeout and retries. The site is prone to long
## hangs (60s+, no response at all) on some queries -- in particular ones
## that turn up few or no matches -- rather than an HTTP error, so a
## finite timeout is essential; without one a bad search can block forever.
##
## A failed attempt that a subsequent retry recovers from is only ever
## `message()`d (informational, suppressible with suppressMessages(), and
## not something the caller's result was actually affected by); it is
## never `warning()`d here. Only a *total* failure (every attempt
## exhausted, NULL returned) is worth a formal warning, and that is left
## to the caller, which already has context-specific wording -- otherwise
## a single flaky request could raise one warning per retry plus another
## from the caller for the same underlying failure.
mathgen_fetch <- function(path, body = NULL, method = c("POST", "GET"),
                           timeout = 30, retries = 1, pause = 2,
                           quiet = FALSE) {
    method <- match.arg(method)
    url <- paste0(mathgen_base_url, "/", path)
    for (attempt in seq_len(retries + 1)) {
        resp <- tryCatch({
            if (method == "POST") {
                httr::POST(url, body = body, encode = "form", httr::timeout(timeout))
            } else {
                httr::GET(url, query = body, httr::timeout(timeout))
            }
        }, error = function(e) e)
        if (inherits(resp, "error")) {
            problem <- sprintf("failed: %s", conditionMessage(resp))
        } else if (httr::status_code(resp) >= 400) {
            problem <- sprintf("HTTP %d", httr::status_code(resp))
        } else {
            return(resp)
        }
        if (!quiet) {
            message(sprintf("attempt %d/%d fetching %s: %s%s", attempt, retries + 1, url,
                             problem, if (attempt <= retries) "; retrying" else ""))
        }
        if (attempt <= retries) Sys.sleep(pause)
    }
    NULL
}
