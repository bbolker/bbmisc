library(ape)
library(ggtree)
library(ggimage)
library(magick)
library(dplyr)

tree     <- read.tree("tt_submission/penguin.nwk")
img_dir  <- normalizePath("penguin_images")
img_name <- "gg_penguin_phylo.png"

## ── 1. Create scaled copies with equal diagonal ──────────────────────────────
## Target diagonal in pixels; all images are rescaled to this.
DIAG_PX <- 400L

tip_labels <- tree$tip.label

for (sp in tip_labels) {
  src  <- file.path(img_dir, paste0(sp, ".jpg"))
  dest <- file.path(img_dir, paste0(sp, "_scaled.png"))
  if (file.exists(dest)) next
  img  <- image_read(src)
  info <- image_info(img)
  diag <- sqrt(info$width^2 + info$height^2)
  new_w <- round(info$width  * DIAG_PX / diag)
  new_h <- round(info$height * DIAG_PX / diag)
  image_resize(img, sprintf("%dx%d!", new_w, new_h)) |>
    image_extent(sprintf("%dx%d", DIAG_PX, DIAG_PX),   # pad to square so
                 gravity = "Center", color = "white") |> # geom_image sizes consistently
    image_strip() |>   # remove ICC profiles that confuse ggimage/magick
    image_write(dest, format = "png")
  message(sprintf("Scaled %s: %dx%d → %dx%d (diag %.0f → %d)",
                  sp, info$width, info$height, new_w, new_h, diag, DIAG_PX))
}

## ── 2. Build plot ─────────────────────────────────────────────────────────────
p0 <- ggtree(tree, layout = "circular", branch.length = "none",
             color = "grey40", linewidth = 0.8)

tip_data <- p0$data |>
  filter(isTip) |>
  mutate(
    img_path  = file.path(img_dir, paste0(label, "_scaled.png")),
    x_img     = x + 1.2,
    ## abbreviated binomial: "Spheniscus_humboldti" → "S. humboldti"
    sp_label  = sub("^([A-Z])[a-z]+_", "\\1. ", label) |> gsub("_", " ", x = _),
    x_lbl     = x + 3.0,
    ## left-half tips get hjust=1 so text extends leftward from anchor
    hjust_lbl = ifelse(angle > 90 & angle < 270, 1, 0)
  )

p <- p0 +
  geom_image(data  = tip_data,
             aes(x = x_img, y = y, image = img_path),
             size  = 0.105, angle = 0, by = "height") +
  geom_text(data  = tip_data,
            aes(x = x_lbl, y = y, label = sp_label, hjust = hjust_lbl),
            angle = 0, size = 2.5, fontface = "italic") +
  xlim(c(0, 11)) +
  theme(plot.margin = margin(1, 1, 1, 1, "mm"))

ggsave(img_name, p, width = 12, height = 12, dpi = 150, bg = "white")

## trim surrounding whitespace (ggplot2/coord_polar always leaves some padding)
image_read(img_name) |>
  image_trim() |>
  image_border(color = "white", geometry = "20x20") |>  # restore a thin margin
  image_write(img_name, format = "png")

message("Written to ", img_name)
