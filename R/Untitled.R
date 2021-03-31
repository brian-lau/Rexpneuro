
basedir <- '~/ownCloud/behaviordata/data_analyses_Brian'

obj <- read_eventide(name = c('flocky','tess'), include_tracker = F, basedir = basedir) %>%
  read_matched_spike_data(basedir = basedir)

neuron_info <- obj$info %>%
  unnest(cols = neuron_info, names_repair = 'universal')
# %>%
#   mutate(rel_depth = depth - intralaminar_depth, .after=depth)

# clean database
# Exclude outside of pallidum
# Exclude poor isolation
# min_trials
neuron_info %<>% filter(depth >= gpe_depth) %>%
  filter(isi_lt_1_5 < 2)
  #filter(IsoS_f >= 0.8)

# Create unique label for neurons
neuron_info %<>%
  mutate(uname = fct_cross(as_factor(id), as_factor(session), name), .after = name) %>%
  arrange(id, session, uname) %>%
  select(-version, -fname_eventide, -fname_ephys, -trigger, -artifact, -probe_id, -filename,
         -channel, -channel_robust_sd, -exclude_times, -stab_slope, -neg, -pos, -neg_peak, -pos_peak,
         -neg_lead)


pre_cue <- spks_in_window(obj, align = "cue_onset_time", t_start = -0.5, t_end = 0,
                           binwidth = 0.5)


fano_factor <- function(x) {
  n <- map_dbl(x,length)
  var(n) / mean(n)
}

cv_isi <- function(x) {
  isi <- unlist(map(x, diff))
  sd(isi)/mean(isi)
}

isi_mode <- function(x) {
  isi <- unlist(map(x, diff))
  modeest::mlv(isi, method = "parzen", kernel = "gaussian")
}

isi_mode <- function(x) {
  isi <- unlist(map(x, diff))
  modeest::mlv(isi, method = "parzen", kernel = "gaussian")
}

burst_index <- function(x) {
  isi <- unlist(map(x, diff))
  sum(isi<0.010)/sum(isi<.5)
}


stats <- pre_cue %>% arrange(session, uname, counter_total_trials) %>%
  group_by(session, uname) %>%
  summarise(ff = fano_factor(times), cv = cv_isi(times2), isi_mode = isi_mode(times), bi = burst_index(times))

tic()
  stats2 <- pre_cue %>% arrange(session, uname, counter_total_trials) %>%
    group_by(session, uname) %>%
    group_modify(~burst_poisson_surprise(.x), .keep=TRUE)
toc()

# stats2 <- pre_cue[1:10000,] %>% arrange(session, uname, counter_total_trials) %>%
#   group_by(session, uname) %>%
#   group_modify(~burst_poisson_surprise(.x), .keep=TRUE)
#
# stats_pre_cue <- pre_cue %>% arrange(session, uname, counter_total_trials) %>%
#   group_by(session, uname) %>%
#   summarise(ff = fano_factor(times2), cv = cv_isi(times), isi_mode = isi_mode(times))
#
# temp = neuron_info %>%
#   select(id, session, name, uname, depth, area, halfpeak_dur, peak_to_trough_dur, fr, isi_mode, cv, cv2, lv, lvr, pause_rate) %>%
#   left_join(stats, by = c("uname"))

stats2 %<>% mutate(burst_dur = ifelse(is.nan(burst_dur), 0, burst_dur))
stats2 %<>% mutate(burst_IBI = ifelse(is.nan(burst_IBI), 15, burst_IBI))
stats2 %<>% mutate(burst_isi = ifelse(is.nan(burst_isi), 0.1, burst_isi))
stats2 %<>% mutate(burst_SI = ifelse(is.nan(burst_SI), 5, burst_SI))
stats2 %<>% mutate(fr_in_burst = ifelse(is.nan(fr_in_burst), fr, fr_in_burst))

temp = neuron_info %>%
  select(id, session, name, uname, rel_depth, area, halfpeak_dur, peak_to_trough_dur, afc, isi_mode, psp, fr, cv, cv2, lv, lvr, pause_rate, mean_interpause_interval,pause_fraction) %>%
  left_join(stats2, by = c("uname")) %>%
  left_join(stats, by = c("uname")) %>%
  mutate(halfpeak_dur = ifelse(halfpeak_dur>0.0003, 0.0003, halfpeak_dur)) %>%
  mutate(burst_rate = sqrt(ifelse(burst_rate==0, 0.01, burst_rate))) %>%
  mutate(fr_out_burst = sqrt(fr_out_burst))
#mutate(isi_mode.x = log(isi_mode.x+0.001))

df = temp %>%
  select(area, rel_depth, halfpeak_dur, peak_to_trough_dur, burst_frac_spks_in_burst, fr_out_burst)
# df = temp %>%
#   select(area, rel_depth, halfpeak_dur, peak_to_trough_dur, lv, pause_fraction, burst_rate, burst_frac_spks_in_burst, fr_out_burst)

M = as.matrix(df %>% filter(area=="gpe") %>% select(-area, -rel_depth))
rownames(M) <- temp %>% filter(area=="gpe") %>% pull(uname)
heatmaply(M, dendrogram = c("both"), scale = c("column"), k_row = 10, dist_method = "euclidean", hclust_method = "ward.D2", fontsize_row = 4)

heatmaply(M, dendrogram = c("row"), scale = c("column"), k_row = 5, hclust_method = NA,
          row_side_colors = temp %>% select(id,area))
