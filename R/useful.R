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
                   "gpi.hfd" = "#555555",
                   "stn.pos" = "#1B9E77",
                   "stn.neg" = "#D95F02",
                   "stn.polypos" = "#7570B3",
                   "stn.polyneg" = "#E7298A",
                   "stn.NULL" = "#555555"
    )

    id <- c("flocky" = "red", "tess" = "black", "chanel" = "pink")

    condition <- c("go_con" = "#228833", "go" = "#4477AA", "nogo" = "#EE6677")
  } else {

  }

  list(Type = type, AreaType = area_type, Id = id)
}

#' @export
percentile_p <- function(x, h0) {
  half.pval <- mean(x > h0) + 0.5*mean(x == 0.5)
  2*min(c(half.pval, 1 - half.pval))
}

#' @export
find_peaks <- function(x, y, w=1, ...) {
  #https://stats.stackexchange.com/questions/36309/how-do-i-find-peaks-in-a-dataset
  #https://rpubs.com/mengxu/peak_detection
  require(zoo)
  n <- length(y)
  y.smooth <- loess(y ~ x, ...)$fitted
  y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align="center")
  delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]
  i.max <- which(delta <= 0) + w
  list(x=x[i.max], i=i.max, y.hat=y.smooth)
}
