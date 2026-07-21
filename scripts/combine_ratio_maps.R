#!/usr/bin/env Rscript
# Combine se_ratio_map.png (left) and rl_ratio_map.png (right) into a single
# side-by-side figure saved as figures/ratio_maps.png.
#
# Uses ImageMagick's convert (no extra R packages needed).
# Output: 6.5 in wide at 300 DPI = 1950 px per panel, full text width for
# 1-inch margin documents.

panel_w <- 1950  # pixels per panel at 300 DPI for 6.5 in / 2

# Resize each panel to target width, preserving aspect ratio, then append
cmd <- sprintf(
  "convert figures/se_ratio_map.png figures/rl_ratio_map.png -resize %dx +append -density 300 figures/ratio_maps.png",
  panel_w
)
cat("Running:", cmd, "\n")
status <- system(cmd)
if (status != 0) stop("ImageMagick convert failed")

info <- file.info("figures/ratio_maps.png")
cat(sprintf("Saved figures/ratio_maps.png (%.0f KB)\n", info$size / 1024))
