#' @export
colormaps <- function(paper = "pallidum") {

  if (paper == "pallidum") {
    type <- c("hfd" = "#1B9E77",
                          "hfd-p" = "#D95F02",
                          "lfd" = "#7570B3",
                          "lfd-b" = "#E7298A")

    area_type <- c("gpe.hfd" = "#1B9E77",
                               "gpe.hfd-p" = "#D95F02",
                               "gpe.lfd" = "#7570B3",
                               "gpe.lfd-b" = "#E7298A",
                               "gpi.hfd" = "#555555")

    id <- c("flocky" = "red", "tess" = "black")
  } else {

  }

  list(Type = type, AreaType = area_type, Id = id)
}

#' @export
percentile_p <- function(x, h0) {
  half.pval <- mean(x > h0) + 0.5*mean(x == 0.5)
  2*min(c(half.pval, 1 - half.pval))
}
