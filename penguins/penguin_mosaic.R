library(magick)
library(httr2)

species <- c(
  "Aptenodytes forsteri",    "Aptenodytes patagonicus",
  "Eudyptes chrysocome",     "Eudyptes chrysolophus",
  "Eudyptes moseleyi",       "Eudyptes pachyrhynchus",
  "Eudyptes robustus",       "Eudyptes schlegeli",
  "Eudyptes sclateri",       "Eudyptula minor",
  "Megadyptes antipodes",    "Pygoscelis adeliae",
  "Pygoscelis antarcticus",  "Pygoscelis papua",
  "Spheniscus demersus",     "Spheniscus humboldti",
  "Spheniscus magellanicus", "Spheniscus mendiculus"
)

## Query Wikipedia API for the lead image thumbnail of a species page.
## Wikipedia requires a descriptive User-Agent or returns 429.
USERAGENT <- "PenguinMosaic/1.0 (R/magick educational script; bbolker@gmail.com)"

get_wiki_thumb <- function(sp) {
  title <- URLencode(gsub(" ", "_", sp), reserved = FALSE)
  url   <- paste0("https://en.wikipedia.org/api/rest_v1/page/summary/", title)
  r <- tryCatch(
    request(url) |>
      req_headers(`User-Agent` = USERAGENT) |>
      req_perform() |>
      resp_body_json(),
    error = function(e) NULL
  )
  src <- r$thumbnail$source
  if (is.null(src)) NA_character_ else src
}

## Resolve a Wikipedia File: page URL to a direct thumbnail URL via imageinfo API
get_wiki_file_thumb <- function(file_page_url, width = 400) {
  title <- sub(".*wiki/", "", file_page_url)  # e.g. "File:Foo.jpg" (URL-encoded)
  api_url <- paste0(
    "https://en.wikipedia.org/w/api.php?action=query",
    "&titles=", title,
    "&prop=imageinfo&iiprop=url&iiurlwidth=", width,
    "&format=json"
  )
  r <- tryCatch(
    request(api_url) |>
      req_headers(`User-Agent` = USERAGENT) |>
      req_perform() |>
      resp_body_json(),
    error = function(e) NULL
  )
  if (is.null(r)) return(NA_character_)
  src <- r$query$pages[[1]]$imageinfo[[1]]$thumburl
  if (is.null(src)) NA_character_ else src
}

## Species whose Wikipedia lead image is wrong (e.g. a range map).
## Values are Wikipedia File: page URLs pointing to a better photo.
img_overrides <- list(
  "Pygoscelis papua" = "https://en.wikipedia.org/wiki/File:Brown_Bluff-2016-Tabarin_Peninsula%E2%80%93Gentoo_penguin_(Pygoscelis_papua)_03.jpg"
)

## Cache downloaded images to avoid re-fetching on each run
img_dir <- "penguin_images"
dir.create(img_dir, showWarnings = FALSE)

TARGET_W <- 300
TARGET_H <- 380
LABEL_H  <- 36  # pixels reserved at bottom for species name

imgs <- lapply(species, function(sp) {
  fname <- file.path(img_dir, paste0(gsub(" ", "_", sp), ".jpg"))
  override <- img_overrides[[sp]]
  if (!is.null(override)) {
    ## Force re-download if a manual override is set, in case wrong image is cached
    unlink(fname)
  }
  if (!file.exists(fname)) {
    Sys.sleep(0.5)  # respect Wikipedia rate limits
    url <- if (!is.null(override)) {
      message("Using override image for: ", sp)
      get_wiki_file_thumb(override)
    } else {
      get_wiki_thumb(sp)
    }
    if (is.na(url)) {
      message("No Wikipedia image found for: ", sp)
      return(NULL)
    }
    message("Downloading: ", sp)
    request(url) |>
      req_headers(`User-Agent` = USERAGENT) |>
      req_perform() |>
      resp_body_raw() |>
      writeBin(fname)
  }
  tryCatch({
    img <- image_read(fname)
    ## Fit within target dimensions preserving aspect ratio, then pad to exact size
    img <- image_resize(img, sprintf("%dx%d", TARGET_W, TARGET_H - LABEL_H))
    img <- image_extent(img, sprintf("%dx%d", TARGET_W, TARGET_H - LABEL_H),
                        gravity = "Center", color = "white")
    ## Add a white strip at the bottom for the label
    strip <- image_blank(TARGET_W, LABEL_H, color = "white")
    img <- image_append(c(img, strip), stack = TRUE)
    ## Abbreviated italic-style label: "A. forsteri"
    shortname <- sub("^([A-Z])[a-z]+ ", "\\1. ", sp)
    image_annotate(img, shortname,
                   size = 16, font = "DejaVu-Sans-Oblique",
                   gravity = "South", color = "black",
                   location = "+0+6")
  }, error = function(e) {
    message("Could not process image for: ", sp, " — ", conditionMessage(e))
    NULL
  })
})
names(imgs) <- species

valid <- Filter(Negate(is.null), imgs)
n <- length(valid)
message(n, " of ", length(species), " images loaded")

## Arrange in a roughly square grid (6 columns × 3 rows for 18 species)
ncols <- ceiling(sqrt(n))
tile <- image_montage(
  image_join(valid),
  geometry = sprintf("%dx%d+4+4", TARGET_W, TARGET_H),
  tile      = sprintf("%dx", ncols),
  bg        = "white"
)

image_write(tile, "penguin_mosaic.png", format = "png")
message("Written to penguin_mosaic.png")
