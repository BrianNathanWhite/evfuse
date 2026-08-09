# Shared figure theme: base-R look for ggplot figures.
# White background, full black panel box, no grid, black axis text,
# plain facet strips. Matches the base-graphics figures
# (stage1_qq_plots, xi_bootstrap_qq) produced with default par().
theme_bw_nogrid <- function(base_size = 11, ...) {
  theme_bw(base_size = base_size, ...) +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"),
          strip.background = element_blank())
}
