## Extract id+name pairs from the <a> links inside a set of <p> nodes,
## deduplicated by id. Shared by the "Advisor" and "Mentor" cases below,
## which have identical markup ("Label: <a href=id.php?id=...>Name</a>").
mathgen_extract_person_links <- function(p_nodes) {
    links <- rvest::html_elements(p_nodes, "a")
    out <- data.frame(
        id = as.integer(sub(".*[?&]id=(\\d+).*", "\\1", rvest::html_attr(links, "href"))),
        name = rvest::html_text(links), stringsAsFactors = FALSE)
    out[!duplicated(out$id), , drop = FALSE]
}

## Parse a single id.php profile page into its own id/name/school/year plus
## the ids+names of its advisor(s) and student(s). Returns NULL if `doc`
## isn't a normal profile page (nonexistent id, or a server-side error page).
##
## Some mathematicians have multiple degree entries on one page (e.g. an
## earlier degree from a different school), each with its own "Advisor"
## paragraph -- Leibniz has three. All such paragraphs are unioned to get
## the full advisor set; only the first degree entry is used for this
## node's own school/year, since the return schema has room for only one.
mathgen_parse_node <- function(doc, id, include_mentors = TRUE) {
    main <- rvest::html_element(doc, "#mainContent")
    if (is.na(main)) return(NULL)

    name <- trimws(rvest::html_text2(rvest::html_element(main, "h2")))
    if (is.na(name) || !nzchar(name)) return(NULL)

    sy <- mathgen_extract_school_year(main)
    school <- sy$school
    year <- sy$year

    adv_p <- rvest::html_elements(main, xpath = ".//p[starts-with(normalize-space(.), 'Advisor')]")
    advisors <- mathgen_extract_person_links(adv_p)
    if (nrow(advisors) == 0 && include_mentors) {
        ## Some mathematicians have no formal advisor/supervisor recorded
        ## but do have a "Mentor:" paragraph in the same position -- e.g.
        ## Andrew Mattei Gleason (id 13307), mentored by George Mackey but
        ## with no PhD advisor of his own. Treat that mentor as the parent
        ## node in the advisor graph, but only as a fallback: a genuine
        ## Advisor entry always takes precedence over a Mentor one, and
        ## this fallback can be disabled entirely via `include_mentors`.
        mentor_p <- rvest::html_elements(main, xpath = ".//p[starts-with(normalize-space(.), 'Mentor')]")
        advisors <- mathgen_extract_person_links(mentor_p)
    }

    tbl_node <- rvest::html_element(main, "table")
    if (is.na(tbl_node)) {
        students <- data.frame(id = integer(0), name = character(0),
                                school = character(0), year = integer(0),
                                stringsAsFactors = FALSE)
    } else {
        links <- rvest::html_elements(tbl_node, "a")
        stu_ids <- as.integer(sub(".*[?&]id=(\\d+).*", "\\1", rvest::html_attr(links, "href")))
        tbl <- rvest::html_table(tbl_node, header = TRUE)
        stopifnot(length(stu_ids) == nrow(tbl))
        stu_school <- tbl$School
        stu_school[!nzchar(stu_school)] <- NA_character_
        students <- data.frame(id = unname(stu_ids), name = tbl$Name, school = stu_school,
                                year = unname(suppressWarnings(as.integer(tbl$Year))),
                                stringsAsFactors = FALSE, row.names = NULL)
    }

    list(id = id, name = name, school = school, year = year,
         advisors = advisors, students = students)
}
