## BFS distance (in edges) from `root`, following the tree's own
## advisor/student direction -- i.e. "generations away from root", the
## same notion of depth mathgen_traverse() itself used while fetching.
## Returns a named integer vector, id (as character) -> depth.
mathgen_bfs_depth <- function(edges, root, direction) {
    root <- as.character(root)
    if (identical(direction, "students")) {
        src <- as.character(edges$from); dst <- as.character(edges$to)
    } else {
        src <- as.character(edges$to); dst <- as.character(edges$from)
    }
    nbrs <- split(dst, src)

    depth <- new.env(parent = emptyenv())
    assign(root, 0L, envir = depth)
    queue <- list(root)
    while (length(queue) > 0) {
        cur <- queue[[1]]
        queue <- queue[-1]
        d <- get(cur, envir = depth)
        for (nb in nbrs[[cur]]) {
            if (!exists(nb, envir = depth, inherits = FALSE)) {
                assign(nb, d + 1L, envir = depth)
                queue[[length(queue) + 1]] <- nb
            }
        }
    }
    ids <- ls(depth)
    stats::setNames(vapply(ids, get, integer(1), envir = depth), ids)
}

## Restrict an edges/vertices pair to vertices within `maxdepth`
## generations of `root` and with year >= minyear (vertices with no
## recorded year are always kept on the year criterion alone, since it's
## unknown whether they'd violate the cutoff). `root` itself is always
## kept. Lets a large, already-fetched tree be plotted as a smaller slice
## without re-fetching.
##
## Keeping every NA-year vertex regardless of `minyear` means one can
## survive the year filter while the ancestor connecting it back to root
## does not (its own year is known and too old) -- severing its only edge
## and leaving it disconnected. A disconnected vertex isn't meaningfully
## part of "the tree rooted at `root`" for plotting purposes, and
## layout_as_tree() has nowhere principled to put it: it was observed
## dumping every such orphan at the same depth as root's immediate
## advisor, reading as a large, spurious fan-out at generation 1 that
## wasn't really there. So a final connectivity pass drops anything left
## unreachable from `root` once the year/depth filters have been applied.
mathgen_filter_tree <- function(x, root, maxdepth, minyear) {
    vertices <- attr(x, "vertices")
    edges <- as.data.frame(x)
    direction <- attr(x, "direction")
    keep_ids <- vertices$id

    if (is.finite(maxdepth)) {
        depths <- mathgen_bfs_depth(edges, root, direction)
        reachable <- as.integer(names(depths))[depths <= maxdepth]
        keep_ids <- intersect(keep_ids, reachable)
    }

    if (is.finite(minyear)) {
        yr <- vertices$year[match(keep_ids, vertices$id)]
        keep_ids <- keep_ids[is.na(yr) | yr >= minyear]
    }

    keep_ids <- union(as.integer(root), keep_ids)
    edges <- edges[edges$from %in% keep_ids & edges$to %in% keep_ids, , drop = FALSE]

    ## final connectivity pass: drop anything no longer reachable from
    ## root over the (possibly now sparser) surviving edge set
    reachable <- as.integer(names(mathgen_bfs_depth(edges, root, direction)))
    keep_ids <- union(as.integer(root), intersect(keep_ids, reachable))

    vertices <- vertices[vertices$id %in% keep_ids, , drop = FALSE]
    edges <- edges[edges$from %in% keep_ids & edges$to %in% keep_ids, , drop = FALSE]

    list(edges = edges, vertices = vertices)
}

## `graph_from_data_frame()` treats a vertices data frame's own `name`
## column (the person's name) as the special vertex-naming attribute if
## one is present, which would silently discard the numeric MGP id instead
## of using it to match edges. Renaming it to `label` keeps the id as the
## vertex name while the display name stays available for plotted labels.
mathgen_plot_igraph <- function(edges, vertices, root, direction, main, ...) {
    if (!requireNamespace("igraph", quietly = TRUE)) {
        stop("method = \"igraph\" requires the igraph package")
    }
    vertices_for_graph <- vertices
    names(vertices_for_graph)[names(vertices_for_graph) == "name"] <- "label"
    g <- igraph::graph_from_data_frame(edges, directed = TRUE, vertices = vertices_for_graph)

    root_idx <- which(igraph::V(g)$name == as.character(root))
    mode <- if (identical(direction, "students")) "out" else "in"

    igraph::plot.igraph(g, vertex.label = igraph::V(g)$label, vertex.size = 8,
                         vertex.label.cex = 0.7, vertex.label.dist = 1,
                         vertex.label.degree = -pi / 2, edge.arrow.size = 0.4,
                         layout = igraph::layout_as_tree(g, root = root_idx, mode = mode),
                         main = main, ...)
    invisible(g)
}

## Renders via the system `dot` binary (not DiagrammeR/Rgraphviz -- see
## mathgen_extract_person_links()'s neighbor, the note in ?plot.mathgen_df:
## DiagrammeR::grViz() unconditionally replaces every apostrophe with a
## double quote before rendering, which corrupts real names like
## "d'Alembert" into invalid DOT syntax) and opens the resulting SVG in a
## browser, since the point of this method is a crisp, native-resolution,
## scrollable/zoomable layout -- shrinking it into a small plot pane would
## defeat that.
mathgen_plot_dot <- function(edges, vertices, root, direction, main) {
    if (!nzchar(Sys.which("dot"))) {
        stop("method = \"graphviz-dot\" requires a system Graphviz install ",
             "providing the \"dot\" command (e.g. `apt install graphviz` or ",
             "`brew install graphviz`)")
    }
    esc <- function(s) gsub("\"", "\\\"", s, fixed = TRUE)

    lab <- ifelse(is.na(vertices$year), vertices$name,
                  sprintf("%s (%d)", vertices$name, vertices$year))
    node_lines <- sprintf("  \"%d\" [label=\"%s\"];", vertices$id, esc(lab))
    edge_lines <- sprintf("  \"%d\" -> \"%d\";", edges$from, edges$to)
    ## edges are always advisor -> student; for an advisor tree, root is
    ## the sink (the person we started from), so it belongs at the top
    ## with edges flowing up into it (rankdir=BT); for a student tree,
    ## root is the source, so it belongs at the top flowing down (TB).
    rankdir <- if (identical(direction, "students")) "TB" else "BT"

    dot <- paste0(
        "digraph mathgen {\n",
        sprintf(paste("  graph [rankdir=%s, bgcolor=\"transparent\", splines=spline,",
                       "nodesep=0.22, ranksep=0.45, label=\"%s\", labelloc=t, fontsize=18,",
                       "fontname=\"Helvetica\"];\n"), rankdir, esc(main)),
        "  node [shape=box, style=\"rounded,filled\", fillcolor=\"#E8A33D\", color=\"#8a5a12\", fontname=\"Helvetica\", fontsize=11, margin=0.09];\n",
        "  edge [color=\"#8a8578\", arrowsize=0.6];\n",
        paste(node_lines, collapse = "\n"), "\n",
        paste(edge_lines, collapse = "\n"), "\n",
        "}\n")

    dot_file <- tempfile(fileext = ".dot")
    svg_file <- tempfile(fileext = ".svg")
    writeLines(dot, dot_file)
    status <- system2("dot", c("-Tsvg", shQuote(dot_file), "-o", shQuote(svg_file)))
    if (!identical(status, 0L) || !file.exists(svg_file)) {
        stop("the \"dot\" command failed to render the graph")
    }
    utils::browseURL(svg_file)
    invisible(svg_file)
}

## Returns the visNetwork htmlwidget (not invisibly -- letting it
## autoprint is how htmlwidgets are meant to be displayed, in the RStudio
## viewer, a browser, or an R Markdown document).
mathgen_plot_visnetwork <- function(edges, vertices, root, direction, main) {
    if (!requireNamespace("visNetwork", quietly = TRUE)) {
        stop("method = \"visNetwork\" requires the visNetwork package")
    }
    nodes <- data.frame(
        id = vertices$id,
        label = ifelse(is.na(vertices$year), vertices$name,
                        sprintf("%s\n(%d)", vertices$name, vertices$year)),
        title = vertices$name, stringsAsFactors = FALSE)
    edges_vn <- data.frame(from = edges$from, to = edges$to)
    ## mirrors the rankdir logic in mathgen_plot_dot(): root at top either
    ## way, "DU" for an advisor tree (root is the sink) and "UD" for a
    ## student tree (root is the source).
    layout_direction <- if (identical(direction, "students")) "UD" else "DU"

    visNetwork::visNetwork(nodes, edges_vn, main = main) |>
        visNetwork::visEdges(arrows = "to", color = list(color = "#8a8578", opacity = 0.6),
                              smooth = list(type = "cubicBezier", roundness = 0.5)) |>
        visNetwork::visNodes(shape = "box",
                              color = list(background = "#E8A33D", border = "#8a5a12",
                                           highlight = "#f2b955"),
                              font = list(color = "#2a1c05", size = 15), margin = 10) |>
        visNetwork::visHierarchicalLayout(direction = layout_direction, sortMethod = "directed",
                                           levelSeparation = 110, nodeSpacing = 140) |>
        visNetwork::visInteraction(dragNodes = TRUE, dragView = TRUE, zoomView = TRUE,
                                    navigationButtons = TRUE) |>
        visNetwork::visPhysics(enabled = FALSE)
}

#' Plot a mathgen_df advisor/student graph
#'
#' Plots the tree returned by [mathgen_traverse()], [get_advisors()], or
#' [get_students()], rooted at the traversal's starting id. Three
#' rendering backends are available:
#'
#' - `"igraph"` (default): a quick base-graphics plot via
#'   `igraph::plot.igraph()`. No extra system dependencies, but its tree
#'   layout crowds badly past a few dozen nodes.
#' - `"graphviz-dot"`: calls the system `dot` command directly (*not* via
#'   the DiagrammeR/Rgraphviz packages) for a proper layered DAG layout --
#'   by far the prettiest of the three, and true vector SVG. It has no
#'   pan/zoom of its own, so the result is opened in a browser at native
#'   resolution, where scrolling and the browser's own zoom stay crisp.
#'   Requires a system Graphviz install (the `dot` binary). DiagrammeR's
#'   `grViz()` is avoided deliberately: it unconditionally replaces every
#'   apostrophe with a double quote before rendering, which silently
#'   corrupts real names like "d'Alembert" into invalid DOT syntax.
#' - `"visNetwork"`: an interactive htmlwidget (vis.js) with native
#'   pan/zoom/drag. CRAN-only install, no system dependency. Its
#'   auto-layout is a bit busier than dot's, but it's the only one of the
#'   three that's actually explorable.
#'
#' `maxdepth` and `minyear` subset an already-fetched tree for plotting,
#' without re-fetching: `maxdepth` keeps only vertices within that many
#' generations of `root` (following the same advisor/student direction the
#' tree was traversed in), and `minyear` drops vertices whose recorded
#' year is earlier than it (vertices with no recorded year are always
#' kept on the year criterion alone, since it isn't known whether they'd
#' violate the cutoff). Whatever remains is then pruned to only what's
#' still connected to `root`: keeping every NA-year vertex regardless of
#' `minyear` can otherwise leave one stranded once the (dated, and so
#' excluded) ancestor connecting it back to `root` is cut, and a
#' disconnected vertex has no principled place in a tree layout --
#' `layout_as_tree()` was observed dumping every such orphan at the same
#' depth as `root`'s immediate advisor, reading as a large, spurious
#' fan-out one generation down that wasn't really there.
#'
#' @param x a "mathgen_df" object
#' @param method one of "igraph" (default), "graphviz-dot", or "visNetwork"
#' @param root id of the vertex to root the tree layout at; defaults to
#'   the id the traversal was started from (`attr(x, "start_id")`)
#' @param maxdepth maximum number of generations from `root` to include
#'   (Inf, the default, includes everything)
#' @param minyear earliest year to include (-Inf, the default, includes
#'   everything); vertices with no recorded year are kept unless doing so
#'   would leave them disconnected from `root` (see Details)
#' @param main plot title; defaults to "Advisor tree: <name>" or "Student
#'   tree: <name>" depending on the traversal direction
#' @param ... additional arguments passed on to `igraph::plot.igraph()`;
#'   only used when `method = "igraph"`
#' @return depends on `method`: the `igraph` graph object, invisibly, for
#'   `"igraph"`; the path to the rendered SVG file, invisibly, for
#'   `"graphviz-dot"` (which is also opened in a browser); the
#'   `visNetwork` htmlwidget (visibly, so it displays) for `"visNetwork"`.
#'   These three are not interchangeable -- a base-graphics side effect, a
#'   file path, and an auto-printing htmlwidget -- so code that scripts
#'   against the return value of `plot()` (rather than just calling it for
#'   its display side effect) will need to branch on `method`.
#' @export
plot.mathgen_df <- function(x, method = c("igraph", "graphviz-dot", "visNetwork"),
                             root = attr(x, "start_id"), maxdepth = Inf, minyear = -Inf,
                             main = NULL, ...) {
    method <- match.arg(method)
    direction <- attr(x, "direction")

    filtered <- mathgen_filter_tree(x, root = root, maxdepth = maxdepth, minyear = minyear)
    edges <- filtered$edges
    vertices <- filtered$vertices

    if (is.null(main)) {
        label <- if (identical(direction, "students")) "Student tree" else "Advisor tree"
        root_name <- vertices$name[vertices$id == root]
        main <- if (length(root_name) > 0) sprintf("%s: %s", label, root_name[1]) else label
    }

    switch(method,
        "igraph" = mathgen_plot_igraph(edges, vertices, root, direction, main, ...),
        "graphviz-dot" = mathgen_plot_dot(edges, vertices, root, direction, main),
        "visNetwork" = mathgen_plot_visnetwork(edges, vertices, root, direction, main))
}
