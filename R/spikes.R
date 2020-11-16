#' @export
read_spike <- function(fname) {
  library(magrittr)
  library(dplyr)
  library(tibble)

  dat <- R.matlab::readMat(fname)

  ## Neuron info
  varNames    <- names(dat$neuron.info[,,1])
  datList     <- dat$neuron.info
  temp = (t(matrix(unlist(datList), nrow = length(varNames))))
  colnames(temp) <- varNames
  neuron_info = as_tibble(temp)
  neuron_info %<>% mutate(tStart = as.numeric(tStart),
                          tEnd = as.numeric(tEnd),
                          channel = as.integer(channel),
                          unit = as.integer(unit),
                          gap = as.integer(gap)
                          ) %>%
    relocate(tStart, .before = tEnd)

  session = temp = unlist(dat$session[, , 1])

  spk = matrix(dat$times, ncol = length(neuron_info$name))
  colnames(spk) <- neuron_info$name

  f <- function(x) as.vector(unlist(x))
  spk <- apply(spk, c(1,2), f)
  spktimes <- as_tibble(spk) %>%
    add_column(counter_total_trials = 1:nrow(spk), .before = 1)

  quality <- dat$quality
  colnames(quality) <- neuron_info$name

  out <- list(
    fname_ephys = session["eFname"],
    fname_eventide = session["bFname"],
    fname_vicon = session["vFname"],
    trigger = as.logical(session["trigger"]),
    trigger_timestamps = as.vector(dat$event.timestamps),
    artifact = as.logical(session["artifact"]),
    target = session["target"],
    grid_x = session["grid.x"],
    grid_y = session["grid.y"],
    tip_depth = as.numeric(session["depth"]),
    n_trials = as.numeric(session["n.trials"]),
    neuron_info = neuron_info,
    spktimes = spktimes,
    quality = quality
  )

  class(out) <- "GNGspikedata"

  return(out)

}
