library(ape)
library(magick)

tree    <- read.tree("tt_submission/penguin.nwk")
img_dir <- "penguin_images"
n_tips  <- Ntip(tree)

## 3 species (chrysocome, antipodes, papua) are absent from the OTL tree
missing <- setdiff(
  gsub(" ", "_", list.files(img_dir, pattern = "\\.jpg$", full.names = FALSE)),
  paste0(tree$tip.label, ".jpg")
)
if (length(missing))
  message("Images present but not in tree: ", paste(missing, collapse = ", "))

## ── output device ────────────────────────────────────────────────────────────
png("penguin_phylo_image.png", width = 1400, height = 1400, res = 130)
par(mar = c(8, 8, 8, 8))  # wide margins so tip images aren't clipped

## ── draw radial tree ─────────────────────────────────────────────────────────
plot(tree, type = "radial", show.tip.label = FALSE,
     edge.width = 2, edge.color = "grey40", open.angle = 10)

## ── recover node coordinates ─────────────────────────────────────────────────
pp    <- get("last_plot.phylo",
             envir = get(".PlotPhyloEnv", envir = asNamespace("ape")))
tip_x <- pp$xx[1:n_tips]
tip_y <- pp$yy[1:n_tips]

## root is node n_tips+1 and sits at the radial centre
root_x <- pp$xx[n_tips + 1]
root_y <- pp$yy[n_tips + 1]

## ── image sizing in plot coordinates ─────────────────────────────────────────
usr      <- par("usr")
img_half <- min(usr[2] - usr[1], usr[4] - usr[3]) * 0.06

## ── place each image at tip ──────────────────────────────────────────────────
for (i in seq_len(n_tips)) {
  sp    <- tree$tip.label[i]
  fname <- file.path(img_dir, paste0(sp, ".jpg"))
  if (!file.exists(fname)) {
    message("No image cached for: ", sp)
    next
  }

  ## outward direction from root → tip; push image centre beyond the tip point
  angle  <- atan2(tip_y[i] - root_y, tip_x[i] - root_x)
  push   <- img_half * 1.2
  img_cx <- tip_x[i] + cos(angle) * push
  img_cy <- tip_y[i] + sin(angle) * push

  img <- image_read(fname) |>
    image_resize("180x180") |>
    image_extent("180x180", gravity = "Center", color = "white")

  rasterImage(as.raster(img),
              xleft   = img_cx - img_half,
              ybottom = img_cy - img_half,
              xright  = img_cx + img_half,
              ytop    = img_cy + img_half,
              interpolate = TRUE,
              xpd = NA)  # allow drawing outside plot region
}

dev.off()
message("Written to penguin_phylo_image.png")
