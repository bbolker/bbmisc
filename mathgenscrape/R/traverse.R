#' Recursively traverse the advisor or student graph
#'
#' BFS over id.php profile pages, following either the advisor links
#' (`direction = "advisors"`, i.e. walking backward toward ancestors) or
#' the student links (`direction = "students"`, walking forward toward
#' descendants). Regardless of direction, edges are always returned as
#' advisor_id -> student_id, so the two directions produce graphs that
#' compose consistently. [get_advisors()] and [get_students()] are thin
#' wrappers around this function and are usually more convenient to call
#' directly.
#'
#' Nodes can have more than one advisor (e.g. Felix Klein was supervised by
#' both Rudolf Lipschitz and Julius Plücker), so both branches are always
#' followed; a visited/queued set prevents re-fetching or infinite loops
#' where lineages reconverge on a shared ancestor/descendant.
#'
#' Descendant trees in particular can be enormous (some 19th-century
#' mathematicians have 100,000+ recorded descendants), so traversal is
#' capped by `max_nodes` and/or `max_depth`; hitting either cap yields a
#' partial graph plus a warning rather than an unbounded crawl.
#'
#' @param start_id starting Math Genealogy Project id
#' @param direction "advisors" to walk backward toward ancestors, or
#'   "students" to walk forward toward descendants
#' @param include_mentors if TRUE (the default), a person with no formal
#'   advisor listed but a "Mentor" recorded instead (e.g. Andrew Mattei
#'   Gleason, id 13307, mentored by George Mackey) has that mentor treated
#'   as their advisor. Only affects `direction = "advisors"`; set to FALSE
#'   to only ever follow genuine Advisor relationships.
#' @param max_nodes maximum number of profile pages to fetch before
#'   stopping (returns a partial graph with a warning if reached)
#' @param max_depth maximum number of generations to recurse (Inf for
#'   unlimited, subject to `max_nodes`)
#' @param delay seconds to pause between page fetches (politeness/rate-limiting)
#' @param timeout request timeout in seconds, passed to the underlying fetch
#' @param retries number of retries per page fetch on failure
#' @param quiet if FALSE, warn/message about failures, truncation, and progress
#' @return a "mathgen_df" object: a data.frame of edges (`from`/`to` ids,
#'   oriented advisor -> student) that also behaves as a plain data.frame
#'   (print, subset, `write.csv()`, ...). The vertex table (`id`/`name`/
#'   `school`/`year`/`visited`) and traversal metadata are attached as
#'   attributes `"vertices"`, `"start_id"`, and `"direction"`; see
#'   [plot.mathgen_df()] for a ready-made plot method.
#' @export
mathgen_traverse <- function(start_id, direction = c("advisors", "students"),
                              include_mentors = TRUE,
                              max_nodes = 500, max_depth = Inf, delay = 0.2,
                              timeout = 30, retries = 1, quiet = FALSE) {
    direction <- match.arg(direction)
    start_id <- as.integer(start_id)

    visited <- new.env(parent = emptyenv())
    queued <- new.env(parent = emptyenv())
    vertices <- list()
    edges <- vector("list", 0)

    add_vertex <- function(id, name, school = NA_character_, year = NA_integer_, seen = FALSE) {
        key <- as.character(id)
        if (is.null(vertices[[key]]) || seen) {
            vertices[[key]] <<- list(id = id, name = name, school = school,
                                      year = year, visited = seen)
        }
    }

    queue <- list(list(id = start_id, depth = 0L))
    assign(as.character(start_id), TRUE, envir = queued)
    truncated <- FALSE
    n_fetched <- 0L

    while (length(queue) > 0) {
        current <- queue[[1]]
        queue <- queue[-1]
        id <- current$id
        key <- as.character(id)
        if (!is.null(visited[[key]])) next
        if (n_fetched >= max_nodes) {
            truncated <- TRUE
            break
        }

        if (n_fetched > 0) Sys.sleep(delay)
        resp <- mathgen_fetch("id.php", body = list(id = id), method = "GET",
                               timeout = timeout, retries = retries, quiet = quiet)
        n_fetched <- n_fetched + 1L
        assign(key, TRUE, envir = visited)
        if (!quiet && n_fetched %% 25 == 0) {
            message(sprintf("...visited %d/%d nodes", n_fetched, max_nodes))
        }

        if (is.null(resp)) {
            if (!quiet) warning(sprintf("failed to fetch id %d; skipping its subtree", id))
            next
        }
        doc <- httr::content(resp, as = "parsed", type = "text/html", encoding = "UTF-8")
        node <- mathgen_parse_node(doc, id, include_mentors = include_mentors)
        if (is.null(node)) {
            if (!quiet) warning(sprintf("could not parse profile page for id %d (nonexistent id, or a server-side error page)", id))
            next
        }
        add_vertex(id, node$name, node$school, node$year, seen = TRUE)

        neighbors <- if (direction == "advisors") node$advisors else node$students
        for (i in seq_len(nrow(neighbors))) {
            nb_id <- neighbors$id[i]
            nb_school <- if (!is.null(neighbors$school)) neighbors$school[i] else NA_character_
            nb_year <- if (!is.null(neighbors$year)) neighbors$year[i] else NA_integer_
            add_vertex(nb_id, neighbors$name[i], nb_school, nb_year)
            edges[[length(edges) + 1]] <- if (direction == "advisors") {
                c(from = nb_id, to = id)
            } else {
                c(from = id, to = nb_id)
            }
            nb_key <- as.character(nb_id)
            if (current$depth + 1 <= max_depth && is.null(visited[[nb_key]]) && is.null(queued[[nb_key]])) {
                queue[[length(queue) + 1]] <- list(id = nb_id, depth = current$depth + 1L)
                assign(nb_key, TRUE, envir = queued)
            }
        }
    }

    if (truncated && !quiet) {
        warning(sprintf(paste(
            "traversal stopped after visiting max_nodes = %d nodes;",
            "the result is a partial graph. Increase `max_nodes` (or narrow",
            "`max_depth`) for a complete traversal."), max_nodes))
    }

    edges_df <- if (length(edges) > 0) {
        m <- do.call(rbind, edges)
        unique(data.frame(from = m[, "from"], to = m[, "to"], row.names = NULL))
    } else {
        data.frame(from = integer(0), to = integer(0))
    }

    vertices_df <- if (length(vertices) > 0) {
        do.call(rbind, lapply(vertices, function(v) {
            data.frame(id = v$id, name = v$name, school = v$school, year = v$year,
                       visited = v$visited, stringsAsFactors = FALSE)
        }))
    } else {
        data.frame(id = integer(0), name = character(0), school = character(0),
                   year = integer(0), visited = logical(0))
    }
    rownames(vertices_df) <- NULL

    new_mathgen_df(edges_df, vertices = vertices_df, start_id = start_id, direction = direction)
}

## Construct a "mathgen_df" object: the edges data.frame itself, classed
## and carrying the vertex table plus traversal metadata as attributes, so
## it still behaves as a plain data.frame (print, subset, write.csv, ...)
## while also supporting plot.mathgen_df().
new_mathgen_df <- function(edges, vertices, start_id, direction) {
    attr(edges, "vertices") <- vertices
    attr(edges, "start_id") <- start_id
    attr(edges, "direction") <- direction
    class(edges) <- c("mathgen_df", class(edges))
    edges
}

## Resolve a get_advisors()/get_students() call's starting id: use `id`
## directly if supplied, otherwise resolve `name` via mathgen_lookup() and
## require it to identify exactly one person.
##
## `id` is always numeric on this site, so a character `id` that isn't a
## bare number (e.g. get_advisors("Catalin Zara"), where the string lands
## on `id` because it's the first positional argument) almost certainly
## means the caller meant `name`. Treated as such here rather than letting
## it silently coerce to NA and fail deep inside mathgen_traverse()'s
## fetch loop with a confusing "NAs introduced by coercion" warning.
mathgen_resolve_id <- function(id, name, timeout, retries, quiet) {
    if (!is.null(id)) {
        if (is.character(id) && is.na(suppressWarnings(as.integer(id)))) {
            if (!quiet) {
                message(sprintf("`id` = %s is not numeric; treating it as `name` instead",
                                 sQuote(id)))
            }
            name <- id
        } else {
            return(id)
        }
    }
    if (is.null(name)) {
        stop("must supply either `id` or `name`")
    }
    lookup <- mathgen_lookup(name, timeout = timeout, retries = retries, quiet = quiet)
    if (nrow(lookup) == 0) {
        stop(sprintf("no matches found for name = %s", sQuote(name)))
    }
    if (nrow(lookup) > 1) {
        choices <- sprintf("%d (%s %s)", lookup$id, lookup$given_name, lookup$surname)
        stop(sprintf("multiple matches found for name = %s; supply `id` directly, one of:\n%s",
                      sQuote(name), paste(choices, collapse = "\n")))
    }
    lookup$id[1]
}

#' Trace a mathematician's advisors, recursively
#'
#' Starting from a Math Genealogy Project id, walks backward through advisor
#' links to build the ancestor tree (advisors, their advisors, and so on).
#' Nodes may have more than one advisor -- e.g. Felix Klein was supervised
#' by both Rudolf Lipschitz and Julius Plücker -- in which case every
#' branch is followed. A thin wrapper around [mathgen_traverse()].
#'
#' @param id starting Math Genealogy Project id. If `NULL`, `name` is used
#'   instead to resolve one via [mathgen_lookup()]. A non-numeric string
#'   passed positionally as `id` (e.g. `get_advisors("Catalin Zara")`) is
#'   treated as `name` rather than failing, since ids on this site are
#'   always numeric.
#' @param name a name to look up via [mathgen_lookup()] when `id` is not
#'   supplied. Ignored if `id` is given. Fails if the lookup has zero hits,
#'   or more than one (the error message lists the matching ids so you can
#'   retry with `id` directly).
#' @param include_mentors if TRUE (the default), a person with no formal
#'   advisor listed but a "Mentor" recorded instead (e.g. Andrew Mattei
#'   Gleason, id 13307, mentored by George Mackey) has that mentor treated
#'   as their advisor; set to FALSE to only follow genuine Advisor
#'   relationships.
#' @param max_nodes maximum number of profile pages to fetch before
#'   stopping (returns a partial graph with a warning if reached)
#' @param max_depth maximum number of generations to recurse (Inf for
#'   unlimited, subject to `max_nodes`)
#' @param delay seconds to pause between page fetches (politeness/rate-limiting)
#' @param timeout,retries,quiet passed to the underlying fetch; see [mathgen_traverse()]
#' @return a "mathgen_df" object; see [mathgen_traverse()] for details and
#'   [plot.mathgen_df()] for a ready-made plot method.
#' @export
get_advisors <- function(id = NULL, name = NULL, include_mentors = TRUE,
                          max_nodes = 500, max_depth = Inf, delay = 0.2,
                          timeout = 30, retries = 1, quiet = FALSE) {
    id <- mathgen_resolve_id(id, name, timeout, retries, quiet)
    mathgen_traverse(id, direction = "advisors", include_mentors = include_mentors,
                      max_nodes = max_nodes, max_depth = max_depth, delay = delay,
                      timeout = timeout, retries = retries, quiet = quiet)
}

#' Trace a mathematician's students, recursively
#'
#' Starting from a Math Genealogy Project id, walks forward through student
#' links to build the descendant tree (students, their students, and so
#' on). Descendant trees can be very large -- some historical
#' mathematicians have well over 100,000 recorded descendants -- so
#' `max_nodes`/`max_depth` matter much more here than for [get_advisors()];
#' consider setting them explicitly rather than relying on the defaults for
#' anyone with a long academic lineage. A thin wrapper around
#' [mathgen_traverse()].
#'
#' @inheritParams get_advisors
#' @return a "mathgen_df" object; see [mathgen_traverse()] for details and
#'   [plot.mathgen_df()] for a ready-made plot method.
#' @export
get_students <- function(id = NULL, name = NULL, max_nodes = 500, max_depth = Inf, delay = 0.2,
                          timeout = 30, retries = 1, quiet = FALSE) {
    id <- mathgen_resolve_id(id, name, timeout, retries, quiet)
    mathgen_traverse(id, direction = "students", max_nodes = max_nodes,
                      max_depth = max_depth, delay = delay, timeout = timeout,
                      retries = retries, quiet = quiet)
}
