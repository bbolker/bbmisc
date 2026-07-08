## Split "Surname, Given Names" into surname/given_name character vectors.
## Returns a data.frame (rather than a matrix) so that single-row inputs
## don't trip the well-known R gotcha where `m[, "col"]` on a one-row
## matrix takes its `names()` from the column label instead of dropping
## them -- that spurious name then leaks into the final data.frame's
## row.names.
mathgen_split_name <- function(full_name) {
    parts <- strsplit(full_name, ",\\s*")
    surname <- unname(vapply(parts, `[`, character(1), 1))
    given <- unname(vapply(parts, function(x) if (length(x) > 1) x[2] else NA_character_,
                            character(1)))
    data.frame(surname = surname, given_name = given, stringsAsFactors = FALSE)
}

mathgen_empty_result <- function() {
    data.frame(id = integer(0), surname = character(0), given_name = character(0),
               school = character(0), year = integer(0), stringsAsFactors = FALSE)
}

## Extract "school, year" of a degree from the id.php profile div that
## holds it (e.g. "Ph.D. Helsingin yliopisto 1932"); the school name is
## pulled from its own inner span rather than split out of the degree text
## since the degree label itself ("Ph.D.", "Dr. rer. nat.", "Theol. Dr.",
## sometimes blank) is not a fixed vocabulary. Shared by node parsing and
## by the search-results single-match fallback below.
##
## The xpath below deliberately excludes divs that themselves contain
## *another* matching div: every ancestor of the true target div (up to
## and including #paddingWrapper) also "contains a span styled 006633" as
## a descendant, so a plain `[1]` selects the outermost such ancestor --
## #paddingWrapper itself -- rather than the small div actually meant.
## When that target div's year is genuinely blank (e.g. Avicenna, id
## 298616), the regex below then matched the first stray 4-digit run
## anywhere on the whole page, which turned out to be the leading digits
## of the "N descendants" count -- the source of the "2450"/"2445" garbage
## years that showed up in a few dozen entries.
mathgen_extract_school_year <- function(main) {
    school <- NA_character_
    year <- NA_integer_
    deg_div <- rvest::html_element(
        main, xpath = paste0(".//div[.//span[contains(@style,'006633')] and ",
                              "not(.//div[.//span[contains(@style,'006633')]])][1]"))
    if (!is.na(deg_div)) {
        school <- trimws(rvest::html_text2(rvest::html_element(
            deg_div, xpath = ".//span[contains(@style, '006633')]")))
        deg_text <- rvest::html_text2(deg_div)
        yr <- regmatches(deg_text, regexpr("[0-9]{4}", deg_text))
        if (length(yr)) year <- as.integer(yr)
    }
    list(school = school, year = year)
}

## Split a "Given Middle Names Surname" string (as used in <h2> headers and
## advisor-paragraph links, no comma) into surname/given_name by treating
## the last word as the surname. This is only a heuristic -- it breaks on
## multi-word surnames (e.g. "Yemeli Tido") -- but there's no delimiter to
## do better, and it is only used as a fallback (see below); results.php
## and quickSearch.php list pages instead give "Surname, Given" directly,
## which mathgen_split_name handles exactly.
mathgen_split_name_heuristic <- function(full_name) {
    words <- strsplit(trimws(full_name), "\\s+")[[1]]
    if (length(words) < 2) return(c(surname = full_name, given_name = NA_character_))
    c(surname = words[length(words)],
      given_name = paste(words[-length(words)], collapse = " "))
}

## Shared response -> data.frame parsing for both quickSearch.php (used by
## mathgen_lookup) and query-prep.php/results.php (used by
## mathgen_advanced_search) -- both are served from the same site template.
mathgen_parse_results <- function(resp, label, quiet) {
    if (is.null(resp)) {
        if (!quiet) {
            warning(sprintf(paste(
                "search for %s timed out or failed; the Math Genealogy Project",
                "is known to hang (rather than error) on some queries -- try",
                "different search terms or increase `timeout`"), label))
        }
        return(mathgen_empty_result())
    }

    doc <- httr::content(resp, as = "parsed", type = "text/html", encoding = "UTF-8")
    main <- rvest::html_element(doc, "#mainContent")
    if (is.na(main)) {
        ## Not a normal site page at all -- e.g. the backend has returned a
        ## raw error such as "MDB2 Error: connect failed". This happens
        ## occasionally when the site's database is overloaded/down; it is
        ## a server-side failure, not "zero matches".
        if (!quiet) {
            warning(sprintf("unexpected response for %s, possible server-side error: %s",
                             label, trimws(substr(rvest::html_text(doc), 1, 200))))
        }
        return(mathgen_empty_result())
    }
    ## When a search matches exactly one person, the site sometimes
    ## redirects straight to that person's own id.php profile page instead
    ## of showing a one-row results list -- confirmed by observing the
    ## same query return each form on different requests. A profile page
    ## always has an <h2> name heading, and *may or may not* have a table
    ## of its own (its "Students:" table, present whenever that person has
    ## any) -- so this has to be checked before assuming any table found
    ## on the page is a results-list table. Getting the order wrong here
    ## previously broke on anyone whose own students table has a header
    ## row (e.g. Victor Guillemin, id 26899): the header row has no link,
    ## so it was one row short of the link count and tripped the
    ## length(ids) == nrow(tbl) check below.
    h2 <- rvest::html_element(main, "h2")
    if (!is.na(h2)) {
        id_text <- rvest::html_text(doc)
        own_id <- suppressWarnings(as.integer(
            sub(".*MGP ID of\\s*([0-9]+).*", "\\1", id_text)))
        name_parts <- mathgen_split_name_heuristic(trimws(rvest::html_text2(h2)))
        sy <- mathgen_extract_school_year(main)
        return(data.frame(id = own_id, surname = name_parts["surname"],
                           given_name = name_parts["given_name"],
                           school = sy$school, year = sy$year,
                           stringsAsFactors = FALSE, row.names = NULL))
    }

    ## A zero-record search still renders an empty `<table>` element (no
    ## `<tr>`s at all) rather than omitting the table entirely, so an
    ## `is.na(tbl_node)` check alone doesn't catch it -- html_table() would
    ## then return a 0x0 tibble and blow up downstream (e.g. `tbl$X1` is
    ## NULL, which strsplit() rejects with "non-character argument").
    tbl_node <- rvest::html_element(main, "table")
    if (is.na(tbl_node) || length(rvest::html_elements(tbl_node, "tr")) == 0) {
        if (!quiet) {
            msg <- trimws(rvest::html_text(rvest::html_element(main, "p")))
            message(sprintf("no matches found for %s (%s)", label, msg))
        }
        return(mathgen_empty_result())
    }

    links <- rvest::html_elements(tbl_node, "a")
    ids <- as.integer(sub(".*[?&]id=(\\d+).*", "\\1", rvest::html_attr(links, "href")))
    tbl <- rvest::html_table(tbl_node, header = FALSE)
    stopifnot(length(ids) == nrow(tbl))

    ## the results page reports how many records it found; if that disagrees
    ## with the number of rows actually parsed, results may be paginated
    ## (not observed in practice, but the site has not been checked with
    ## large enough result sets to be sure it never happens)
    reported_text <- rvest::html_text(rvest::html_element(main, "p"))
    n_reported <- suppressWarnings(
        as.integer(sub(".*found ([0-9]+) records.*", "\\1", reported_text)))
    if (!quiet && !is.na(n_reported) && n_reported != nrow(tbl)) {
        warning(sprintf(
            "results page for %s reports %d records but only %d were parsed (possible pagination)",
            label, n_reported, nrow(tbl)))
    }

    name_cols <- mathgen_split_name(tbl$X1)
    school <- tbl$X2
    school[!nzchar(school)] <- NA_character_

    data.frame(id = unname(ids), surname = name_cols$surname,
               given_name = name_cols$given_name, school = unname(school),
               year = unname(suppressWarnings(as.integer(tbl$X3))),
               stringsAsFactors = FALSE, row.names = NULL)
}
