wes_palette <- function(name = "Zissou1", n = 4, type = c("discrete", "continuous"), ...) {
  type <- match.arg(type)
  if (is.null(n) || !is.numeric(n) || n <= 0) n <- 4
  n <- as.integer(n)

  # First try wesanderson if available
  if (requireNamespace("wesanderson", quietly = TRUE)) {
    tryCatch({
      return(wesanderson::wes_palette(name = name, n = n, type = type))
    }, error = function(e) {
      # If wesanderson fails, fall back to RColorBrewer or default
      message("Note: wesanderson palette '", name, "' not found, using fallback.")
    })
  }

  # Fallback to RColorBrewer for common palette names
  if (requireNamespace("RColorBrewer", quietly = TRUE)) {
    # Map common palette names to RColorBrewer palettes
    rcb_name <- switch(name,
      "Dark2" = "Dark2",
      "Set2" = "Set2",
      "Set1" = "Set1",
      "Set3" = "Set3",
      "Pastel1" = "Pastel1",
      "Pastel2" = "Pastel2",
      "Accent" = "Accent",
      "Paired" = "Paired",
      "Zissou1" = "Set2",  # approximate
      "GrandBudapest1" = "Set2",
      "Moonrise1" = "Set2",
      "Royal1" = "Set2",
      "Darjeeling1" = "Set2",
      "Chevalier1" = "Set2",
      "FantasticFox1" = "Set2",
      "Cavalcanti1" = "Set2",
      "Rushmore1" = "Set2",
      "BottleRocket1" = "Set2",
      "AsteroidCity1" = "Set2",
      "Set2"  # default fallback
    )
    if (type == "discrete") {
      max_colors <- RColorBrewer::brewer.pal.info[rcb_name, "maxcolors"]
      if (n <= max_colors) {
        return(RColorBrewer::brewer.pal(n = n, name = rcb_name))
      } else {
        # Generate colors using grDevices
        return(grDevices::hcl.colors(n, "Dark 3"))
      }
    } else {
      # Continuous palette - use viridis
      return(grDevices::hcl.colors(n, "Viridis"))
    }
  }

  # Final fallback
  if (type == "discrete") {
    return(grDevices::hcl.colors(n, "Dark 3"))
  }
  grDevices::hcl.colors(n, "Viridis")
}
