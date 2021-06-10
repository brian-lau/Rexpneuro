#' @export
epcs.test <- function(data.cells, z.adjust.method=c('none','hommel')){
  # https://rpubs.com/seriousstats/epcs_test
  #https://seriousstats.wordpress.com/2019/09/05/chi-square-and-the-egon-pearson-correction/
  #https://sites.google.com/a/lakeheadu.ca/bweaver/Home/statistics/notes/chisqr_assumptions
  ucs.test <- suppressWarnings(stats::chisq.test(data.cells, simulate.p.value=FALSE, correct = FALSE))
  N <- sum(data.cells) ; nrows <- dim(data.cells)[1] ; ncols <- dim(data.cells)[2]
  corrected.stat <- ucs.test$stat[[1]] * (N-1)/N
  pval <- pchisq(corrected.stat, ucs.test$par, lower.tail = FALSE)
  p.resids <- ucs.test$resid ; as.resids <- ucs.test$stdres
  pr.pv <- pchisq(ucs.test$resid^2,1, lower.tail=F)
  pr.apv <- matrix(as.matrix(p.adjust(pchisq(p.resids^2, 1, lower.tail=F), z.adjust.method[1])), nrows, ncols, byrow=F)
  asr.apv <- matrix(as.matrix(p.adjust(pchisq(as.resids^2,1, lower.tail=F), z.adjust.method[2])), nrows, ncols, byrow=F)
  output.list <- list(
    'Uncorrected Pearson Chi-Square'= ucs.test$stat,
    'Egon Pearson Chi-Square' = corrected.stat,
    'df'=prod(dim(data.cells)-1), 'p'=pval,
    'Smallest expected value (should be greater than 1)' = min(ucs.test$expected),
    'Raw data' = data.cells,
    'Pearson (standardized) residuals (z)' = p.resids,
    'two-sided p (for z)' = pr.pv, 'adjusted p (for z)' = pr.apv,
    'Chi-square contribution per cell' = p.resids^2,
    'Adjusted standardized residuals' = as.resids,
    'adjusted p (for ASRs)' = asr.apv,
    'Pearson / ASR p adjustment' =z.adjust.method)
  return(output.list)
}
