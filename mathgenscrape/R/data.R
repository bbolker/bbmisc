#' Ethan Bolker's advisor tree
#'
#' An example "mathgen_df" object: the recursive advisor tree for Ethan
#' Bolker (Math Genealogy Project id 94376), as returned by
#' `get_advisors(id = 94376)`. Included for convenience in examples,
#' `plot()`/`print()` demonstrations, and vignettes, so they don't depend
#' on a live, sometimes-slow network call to the Math Genealogy Project.
#'
#' @format A "mathgen_df" object (see [mathgen_traverse()]): a data.frame
#'   of 144 edges (`from`/`to` ids, oriented advisor -> student), with the
#'   131-row vertex table and traversal metadata attached as attributes
#'   `"vertices"`, `"start_id"` (94376), and `"direction"` ("advisors").
#' @source <https://genealogy.math.ndsu.nodak.edu/id.php?id=94376>,
#'   fetched with `get_advisors(id = 94376)`.
#' @examples
#' print(ethan_bolker)
#' plot(ethan_bolker)
"ethan_bolker"

#' Catalin Zara's advisor tree
#'
#' An example "mathgen_df" object: the recursive advisor tree for Catalin
#' Zara (Math Genealogy Project id 47185), as returned by
#' `get_advisors(id = 47185)`. Larger and more densely interconnected than
#' [ethan_bolker], with several vertices having 2-3 recorded advisors, so
#' it's a useful stress test for layout/rendering (e.g. it shows how
#' `method = "graphviz-dot"` and `method = "igraph"` can place the same
#' vertex at visually different "generations" when its shortest and
#' longest path back from the root differ -- not a data error, just what
#' happens when a DAG this convergent is forced into a layered layout).
#' Included for convenience in examples, `plot()`/`print()`
#' demonstrations, and vignettes, so they don't depend on a live,
#' sometimes-slow network call to the Math Genealogy Project.
#'
#' @format A "mathgen_df" object (see [mathgen_traverse()]): a data.frame
#'   of 299 edges (`from`/`to` ids, oriented advisor -> student), with the
#'   243-row vertex table and traversal metadata attached as attributes
#'   `"vertices"`, `"start_id"` (47185), and `"direction"` ("advisors").
#' @source <https://genealogy.math.ndsu.nodak.edu/id.php?id=47185>,
#'   fetched with `get_advisors(id = 47185)`.
#' @examples
#' print(catalin_zara)
#' plot(catalin_zara)
"catalin_zara"
