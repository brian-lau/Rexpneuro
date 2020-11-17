#' @export
read_matching_spkdata <- function(obj,
                          basedir = getwd()
) {
  # Incoming eventide files to match to corresponding spike data
  fnames_eventide <- obj$fname
  base_name <- purrr::map_chr(stringr::str_split(fnames_eventide, ".txt"), 1)

  # Available spike data
  fnames <- list.files(path = basedir, pattern =
                         glob2rx(paste0("*GNG", "*.mat")))
  available_name <- purrr::map_chr(stringr::str_split(fnames, ".mat"), 1)

  # Which eventide sessions having corresponding spike data?
  ind <- which(base_name %in% available_name)

  fnames_spk <- paste0(base_name[ind], ".mat")

  spkdat <- purrr::map(fnames_spk, read_spike)

  out <- list(
    evi = obj,
    spk = spkdat,
  )

  return(out)
}
