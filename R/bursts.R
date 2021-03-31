## maxinterval.R --- maxinterval burst detection (from Neuroexplorer).
## Author: Stephen Eglen
## Copyright: GPL
## Fri 23 Feb 2007

#' @export
mi.find.bursts <- function(spikes,
                           beg.isi=0.006,
                           end.isi=0.009,
                           min.ibi=0.020,
                           min.durn=0.005,
                           min.spikes=3,
                           debug=FALSE) {

  ## For one spike train, find the burst using max interval method.
  ## e.g.
  ## find.bursts(s$spikes[[5]])
  ## init.
  ## params currently in MI.PAR
  ##

  ## TODO: all our burst analysis routines should use the same
  ## value to indiciate "no bursts" found.
  ##no.bursts = NA;                       #value to return if no bursts found.
  no.bursts = matrix(nrow=0,ncol=1)     #emtpy value nrow()=length() = 0.

  # par = mi.par
  # beg.isi =    par$beg.isi
  # end.isi =    par$end.isi
  # min.ibi =    par$min.ibi
  # min.durn =   par$min.durn
  # min.spikes = par$min.spikes

  nspikes = length(spikes)

  ## Create a temp array for the storage of the bursts.  Assume that
  ## it will not be longer than Nspikes/2 since we need at least two
  ## spikes to be in a burst.

  max.bursts <- floor(nspikes/2)
  bursts <- matrix(NA, nrow=max.bursts, ncol=3)
  colnames(bursts) = c("beg", "end", "IBI")
  burst <- 0                            #current burst number

  ## Phase 1 -- burst detection.  Here a burst is defined as starting
  ## when two consecutive spikes have an ISI less than BEG.ISI apart.
  ## The end of the burst is given when two spikes have an ISI greater
  ## than END.ISI.

  ## Find ISIs closer than beg.isi, and end with end.isi.


  ## LAST.END is the time of the last spike in the previous burst.
  ## This is used to calculate the IBI.
  ## For the first burst, this is no previous IBI
  last.end = NA;                        #for first burst, there is no IBI.

  n = 2
  in.burst = FALSE

  while ( n <= nspikes) {

    next.isi = spikes[n] - spikes[n-1]
    if (in.burst) {
      if (next.isi > end.isi) {
        ## end of burst
        end = n-1; in.burst = FALSE


        ibi =  spikes[beg] - last.end; last.end = spikes[end]
        res = c(beg, end, ibi)
        burst = burst + 1
        if (burst > max.bursts) {
          print("too many bursts!!!")
          browser()
        }
        bursts[burst,] <- res
      }
    } else {
      ## not yet in burst.
      if (next.isi < beg.isi) {
        ## Found the start of a new burst.
        beg = n-1; in.burst = TRUE
      }
    }
    n = n+1
  }

  ## At the end of the burst, check if we were in a burst when the
  ## train finished.
  if (in.burst) {
    end = nspikes
    ibi =  spikes[beg] - last.end
    res = c(beg, end, ibi)
    burst = burst + 1
    if (burst > max.bursts) {
      print("too many bursts!!!")
      browser()
    }
    bursts[burst,] <- res
  }

  ## Check if any bursts were found.
  if (burst > 0 ) {
    ## truncate to right length, as bursts will typically be very long.
    bursts = bursts[1:burst,,drop=FALSE]
  } else {
    ## no bursts were found, so return an empty structure.
    return(no.bursts)
  }

  if (debug) {
    print("End of phase1\n")
    print(bursts)
  }


  ## Phase 2 -- merging of bursts.  Here we see if any pair of bursts
  ## have an IBI less than MIN.IBI; if so, we then merge the bursts.
  ## We specifically need to check when say three bursts are merged
  ## into one.


  ibis = bursts[,"IBI"]
  merge.bursts = which(ibis < min.ibi)

  if (any(merge.bursts)) {
    ## Merge bursts efficiently.  Work backwards through the list, and
    ## then delete the merged lines afterwards.  This works when we
    ## have say 3+ consecutive bursts that merge into one.

    for (burst in rev(merge.bursts)) {
      bursts[burst-1, "end"] = bursts[burst, "end"]
      bursts[burst, "end"] = NA         #not needed, but helpful.
    }
    bursts = bursts[-merge.bursts,,drop=FALSE] #delete the unwanted info.
  }

  if (debug) {
    print("End of phase 2\n")
    print(bursts)
  }


  ## Phase 3 -- remove small bursts: less than min duration (MIN.DURN), or
  ## having too few spikes (less than MIN.SPIKES).
  ## In this phase we have the possibility of deleting all spikes.

  ## LEN = number of spikes in a burst.
  ## DURN = duration of burst.
  len = bursts[,"end"] - bursts[,"beg"] + 1
  durn = spikes[bursts[,"end"]] - spikes[bursts[,"beg"]]
  bursts = cbind(bursts, len, durn)

  rejects = which ( (durn < min.durn) | ( len < min.spikes) )

  if (any(rejects)) {
    bursts = bursts[-rejects,,drop=FALSE]
  }

  if (nrow(bursts) == 0) {
    ## All the bursts were removed during phase 3.
    bursts = no.bursts
  } else {
    ## Compute mean ISIS
    len = bursts[,"end"] - bursts[,"beg"] + 1
    durn = spikes[bursts[,"end"]] - spikes[bursts[,"beg"]]
    mean.isis = durn/(len-1)

    ## Recompute IBI (only needed if phase 3 deleted some cells).
    if (nrow(bursts)>1) {
      ibi2 = c(NA, calc.ibi(spikes, bursts))
    } else {
      ibi2 = NA
    }
    bursts[,"IBI"] = ibi2

    SI = rep(1, length(mean.isis ))
    bursts = cbind(bursts, mean.isis, SI)
  }

  ## End -- return burst structure.
  bursts

}


# Wrapper for PS.method
burst_poisson_surprise <- function(x) {
  burst = map(.x=x$times, .f=~PS.method(.x))

  # Deal with no spikes
  trial_dur = map(.x=x$times, .f=~ifelse(length(.x)>1, max(.x)-min(.x), NA))
  total_spks = map_dbl(.x=x$times, .f=~length(.x))

  rate <- map2_dbl(.x=burst, .y=trial_dur, .f=~ifelse(nrow(.x)>0,
                                                      nrow(.x)/.y,
                                                      0))
  dur <- map2_dbl(.x=burst, .y=trial_dur, .f=~ifelse(nrow(.x)>0,
                                                     mean(.x$durn),
                                                     NA))
  IBI <- map2_dbl(.x=burst, .y=trial_dur, .f=~ifelse(nrow(.x)>1,
                                                     mean(.x$IBI, na.rm = T),
                                                     NA))
  isi <- map2_dbl(.x=burst, .y=trial_dur, .f=~ifelse(nrow(.x)>0,
                                                     mean(.x$mean.isis),
                                                     NA))
  SI <- map2_dbl(.x=burst, .y=trial_dur, .f=~ifelse(nrow(.x)>0,
                                                    mean(.x$SI),
                                                    NA))

  spks_in_burst <- map_dbl(.x=burst, .f=~ifelse(nrow(.x)>0,
                                                sum(.x$len),
                                                0))
  total_burst_dur <- map_dbl(.x=burst, .f=~ifelse(nrow(.x)>0,
                                                sum(.x$durn),
                                                0))

  frac_spks_in_burst <- spks_in_burst/total_spks

  fr <- total_spks/unlist(trial_dur)
  fr_in_burst <- spks_in_burst/total_burst_dur
  fr_out_burst <- (total_spks - spks_in_burst)/(unlist(trial_dur) - total_burst_dur)

  return(data.frame(burst_rate = mean(rate, na.rm = T),
                    burst_dur = mean(dur, na.rm = T),
                    burst_IBI = mean(IBI, na.rm = T),
                    burst_isi = mean(isi, na.rm = T),
                    burst_SI = mean(SI, na.rm = T),
                    burst_frac_spks_in_burst = mean(frac_spks_in_burst, na.rm = T),
                    fr = mean(fr, na.rm = T),
                    fr_in_burst = mean(fr_in_burst, na.rm = T),
                    fr_out_burst = mean(fr_out_burst, na.rm = T)))
}

# Modified from https://github.com/ellesec/burstanalysis/blob/master/Burst_detection_methods/PS_method.R
# Returns empty data.frame for no bursts
PS.method<-function(spike.train, si.thresh=5) {
  si.thresh<-ifelse(is.null(si.thresh), 5, si.thresh)
  burst <- si.find.bursts.thresh(spike.train)
  if (is.null(dim(burst))){
    #result<-NA
    result <- data.frame(beg=numeric(), end=numeric(), IBI=numeric(), len=numeric(), durn=numeric(), mean.isis=numeric(), SI=numeric())
  } else {
    burst.rem<-which(burst[,"SI"]<si.thresh)
    if (length(burst.rem)) {
      burst<-burst[-burst.rem,]
    }
    if (length(dim(burst))<1) {
      burst<-data.frame(beg=burst[1], len=burst[2], SI=burst[3], durn=burst[4], mean.isis=burst[5])
    } else if (dim(burst)[1]==0){
      result <- data.frame(beg=numeric(), end=numeric(), IBI=numeric(), len=numeric(), durn=numeric(), mean.isis=numeric(), SI=numeric())
      return(result)
    }
    beg<-burst[,"beg"]
    len<-burst[,"len"]
    N.burst<-length(beg)
    end<-beg+len-1
    IBI<-c(NA, spike.train[beg[-1]]-spike.train[end[-N.burst]])
    result<-data.frame(beg=beg, end=end, IBI=IBI, len=len, durn=burst[,"durn"], mean.isis=burst[,"mean.isis"], SI=burst[,"SI"])
    rownames(result)<-NULL
  }
  result
}

si.find.bursts.thresh<- function (spikes, debug = FALSE)
{
  nspikes = length(spikes)
  mean.isi = mean(diff(spikes))
  threshold = mean.isi/2
  n = 1
  max.bursts <- floor(nspikes/3)
  bursts <- matrix(NA, nrow = max.bursts, ncol = burst.info.len)
  burst <- 0
  while (n < nspikes - 2) {
    if (debug)
      print(n)
    if (((spikes[n + 1] - spikes[n]) < threshold) && ((spikes[n +
                                                              2] - spikes[n + 1]) < threshold)) {
      res <- si.find.burst.thresh2(n, spikes, nspikes, mean.isi,
                                   burst.isi.max, debug)
      if (is.na(res[1])) {
        n <- n + 1
      }
      else {
        burst <- burst + 1
        if (burst > max.bursts) {
          print("too many bursts")
          browser()
        }
        bursts[burst, ] <- res
        n <- res[1] + res[2]
        names(n) <- NULL
      }
    }
    else {
      n = n + 1
    }
  }
  if (burst > 0) {
    res <- bursts[1:burst, , drop = FALSE]
    colnames(res) <- burst.info
  }
  else {
    res <- NA
  }
  res
}



si.find.burst.thresh2<-function(n, spikes, nspikes, mean.isi, threshold=NULL,
                                debug=FALSE) {
  ## Find a burst starting at spike N.
  ## Include a better phase 1.


  ## Determine ISI threshold.
  if (is.null(threshold))
    isi.thresh = 2 * mean.isi
  else
    isi.thresh = threshold

  if (debug)
    cat(sprintf("** find.burst %d\n", n))

  i=3  ## First three spikes are in burst.
  s = surprise(n, i, spikes, nspikes, mean.isi)

  ## Phase 1 - add spikes to the train.
  phase1 = TRUE
  ##browser()

  ## in Phase1, check that we still have spikes to add to the train.
  while( phase1 ) {

    ##printf("phase 1 s %f\n", s);

    i.cur = i;

    ## CHECK controls how many spikes we can look ahead until SI is maximised.
    ## This is normally 10, but will be less at the end of the train.
    check = min(10, nspikes-(i+n-1))

    looking = TRUE; okay = FALSE;
    while (looking) {

      if (check==0) {
        ## no more spikes left to check.
        looking=FALSE;
        break;
      }
      check=check-1; i=i+1
      s.new = surprise(n, i, spikes, nspikes, mean.isi)
      if (debug)
        printf("s.new %f s %f n %d i %d check %d\n", s.new, s, n, i, check)

      if (s.new > s) {
        okay=TRUE; looking=FALSE;
      } else {
        ## See if we should keep adding spikes?
        if ( (spikes[i] - spikes[i-1]) > isi.thresh ) {
          looking = FALSE;
        }

      }
    }
    ## No longer checking, see if we found an improvement.
    if (okay) {
      if (s > s.new) {
        ## This should not happen.
        printf("before s %f s.new %f\n", s, s.new)
        browser()
      }
      s = s.new
    } else {
      ## Could not add more spikes onto the end of the train.
      phase1 = FALSE
      i = i.cur
    }
  }


  ## start deleting spikes from the start of the burst.
  phase2 = TRUE
  while(phase2) {
    if (i==3) {
      ## minimum length of a burst must be 3.
      phase2=FALSE
    } else {
      s.new = surprise(n+1, i-1, spikes, nspikes, mean.isi)
      if (debug)
        cat(sprintf("phase 2: n %d i %d s.new %.4f\n", n, i, s.new))
      if (s.new > s) {
        if (debug)
          print("in phase 2 acceptance\n")
        n = n+1; i = i-1
        s = s.new
      } else {
        ## removing front spike did not improve SI.
        phase2 = FALSE
      }
    }
  }


  ## End of burst detection; accumulate result.


  ## compute the ISIs, and then the mean ISI.

  ## Fencepost issue: I is the number of spikes in the burst, so if
  ## the first spike is N, the last spike is at N+I-1, not N+I.
  isis = diff(spikes[n+(0:(i-1))])
  mean.isis = mean(isis)

  durn = spikes[n+i-1] - spikes[n]
  res <- c(n=n, i=i, s=s, durn=durn, mean.isis=mean.isis)

  ##browser()
  res

}
