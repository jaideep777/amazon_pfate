library(tidyverse)
rm(list=ls())


input_dir  = here::here("input_data")
output_dir = here::here("pfate_output_mip")

read_pfate_outputs = function(input_dir, output_dir, expt_dir){
  wd_back = getwd()
  setwd(paste0(output_dir,"/",expt_dir))

  l = list(
    # seeds1 = read.delim("seeds.csv", header=F, col.names = paste0("V", 1:(n_species+2)))
    # Zp = read.csv("z_star.csv", header=F, col.names = paste0("V", 1:50)),
    # BA1 = read.csv("basal_area.csv", header=F, col.names = paste0("V", 1:(n_species+2)))
    # co = read.csv("canopy_openness.csv", header=F, col.names = paste0("V", 1:50)),
    # lai_v = read.csv("lai_profile.csv", header=F, col.names = paste0("V", 1:27)),
    traits = read.csv("traits.csv"),
    dat_d = readr::read_csv("D_PFATE.csv"),
    # dat$YEAR = decimal_date(as_date(dat$YEAR, format = "%Y-%m-%d %H:%M:%S GMT (doy = %j)"))
    dat2 = read.csv("Y_PFATE.csv"),
    dat3 = read.csv("Y_mean_PFATE.csv"),
    # dist = readr::read_csv("size_distributions.csv", col_names = F),
    x = exp(seq(log(0.01), log(10), length.out=100))
  )

  # l$dist = l$dist[,-ncol(l$dist)]
  # names(l$dist)[1:2] = c("YEAR", "SPP")
  # names(l$Zp)[1] = c("YEAR")
  # names(l$co)[1] = c("YEAR")
  # names(l$lai_v)[1] = c("YEAR")

  l$dat = l$dat_d %>%
    mutate(YEAR = as.integer(YEAR)) %>%
    group_by(YEAR) %>%
    summarize_all(mean)

  l$traits = l$traits %>% filter(!grepl("probe", .$SPP))

  n_species = l$dat2 %>% filter(!grepl("probe", .$PID)) %>% pull(PID) %>% unique() %>% length()
  n_year = length(unique(l$dat2$YEAR))

  setwd(wd_back)

  l
}

cat_outputs = function(list1, list2){
  keys <- unique(c(names(list1), names(list2)))
  l <- lapply(setNames(keys, keys), function(x) {rbind(list1[[x]], list2[[x]])})
  l$x = list1$x
  l
}

subsample = function(l, interval = 10){
  keys <- unique(c(names(l)))
  keys <- keys[keys != "x"]
  keys <- keys[keys != "dat_d"]

  years = unique(l$dat2$YEAR)

  lsub = lapply(setNames(keys, keys), function(x) {l[[x]] = l[[x]] %>% filter(as.integer(YEAR) %% interval == 0)})
  lsub$x = l$x
  lsub$dat_d = l$dat_d
  lsub
}

slice_time = function(l, ymin, ymax){
  keys <- unique(c(names(l)))
  keys <- keys[keys != "x"]

  years = unique(l$dat2$YEAR)

  lsub = lapply(setNames(keys, keys), function(x) {l[[x]] = l[[x]] %>% filter(YEAR < ymax & YEAR > ymin)})
  lsub$x = l$x
  lsub
}


dirs = list(
  sitename = c("AmzMIP"),
  scenario = c("AMB", "ELE"),
  div = c("evol", "ld")) |>
  cross_df() |>
  filter(!(div=="ld" & scenario=="AMB")) |>
  mutate(dirs = paste0(
    sitename, "_HIST_",
    scenario, "_",
    div, "_",
    case_match(div,
      "evol"~"20ky",
      "ld"~"1ky"),
    "_c2_rs0.04"
  )) |>
  distinct() 

## READ DATA 

data = dirs %>%
  mutate(raw_eco = purrr::map(.x = dirs, .f = function(x){
                  read_pfate_outputs(input_dir, output_dir, x)
    }))

## Process outputs

data_proc = data %>%
  mutate(D = purrr::map(.x = raw_eco, .f = function(x){
    x$dat_d %>%
      mutate(DoY = lubridate::yday(lubridate::date_decimal(YEAR))) %>% 
      mutate(ET = TRANS) %>% 
      select(YEAR, IYEAR, DoY, GPP, NPP, RAU, MORT, GS, ET, VCMAX) %>% 
      mutate(
        GPP = GPP * 1e3, # convert kgC m-2 d-1 --> gC m-2 d-1
        NPP = NPP * 1e3, # convert kgC m-2 d-1 --> gC m-2 d-1
        RAU = RAU * 1e3, # convert kgC m-2 d-1 --> gC m-2 d-1
        MORT = MORT * 1e3, # convert kgC m-2 d-1 --> gC m-2 d-1
        GS = GS * 1e3, # convert mol m-2 s-1 --> mmol m-2 s-1
        ET = ET * 1 # already in mm day-1
      ) %>% 
      left_join(x$dat3 %>%
                  mutate(IYEAR = floor(YEAR)) %>%
                  select(IYEAR, CL, CW, CCR, CFR, CR, LAI) %>% 
                  mutate(across(CL:CR,
                                ~.x*1e3 # convert kgC m-2 --> gC m-2
                                )) 
                ) %>% 
      select(-YEAR) %>% 
      rename(YEAR = IYEAR)
  })) %>%
  mutate(Y = purrr::map(.x = raw_eco, .f = function(x){
    x$dat2 %>%
      left_join(x$traits) %>% 
      filter(!grepl("probe", PID)) %>% 
      mutate(YEAR = floor(YEAR)) %>% 
      mutate(OC = -9999,
             MO = -9999, 
             SLA = 1/LMA) %>% 
      select(YEAR, PID, DE, OC, PH, MH=HMAT, CA, BA, TB, WD, MO, SLA, P50=P50X) # All densities (DE, CA, BA, TB) are per m2
  })) 


## Visualize results

cairo_pdf(here::here("figures/mip_amb_vs_ele.pdf"), width=11.3, height=7)
bind_rows(
  data_proc %>% 
    select(scenario, div, D) %>% 
    unnest(D) %>%
    filter(div == "evol") %>% 
    filter(YEAR %in% seq(-19000, 20000, by=51)) %>% 
    pivot_longer(-(scenario:DoY)) %>% 
    mutate(YEAR=YEAR+DoY/365),
  
  data_proc %>% 
    select(scenario, div, Y) %>% 
    unnest(Y) %>%
    filter(div == "evol") %>% 
    filter(YEAR %in% seq(-19000, 20000, by=51)) %>% 
    pivot_longer(-(scenario:PID))
) %>% 
  ggplot(aes(x=YEAR, y=value, col=scenario)) +
  geom_line(alpha=0.15) + 
  geom_smooth(alpha=0.5, method="loess", span=0.05, se=F) + 
  annotate(geom="rect", xmax=2020, xmin=-Inf, ymin=-Inf, ymax=Inf, fill="grey", alpha=0.2)+
  facet_wrap(~name, scales="free_y")+
  scale_x_continuous(n.breaks = 3)+
  theme_bw()
dev.off()

cairo_pdf(here::here("figures/mip_hd_vs_ld.pdf"), width=11.3, height=7)
bind_rows(
  data_proc %>% 
    select(scenario, div, D) %>% 
    unnest(D) %>%
    group_by(div, scenario) %>% 
    filter(YEAR > max(YEAR)-100) %>% 
    ungroup() %>% 
    pivot_longer(-(scenario:DoY)) %>% 
    group_by(name, scenario, div) %>% 
    summarize(value = mean(value)),
  
  data_proc %>% 
    select(scenario, div, Y) %>% 
    unnest(Y) %>%
    group_by(div, scenario) %>% 
    filter(YEAR > max(YEAR)-100) %>% 
    ungroup() %>% 
    pivot_longer(-(scenario:PID)) %>% 
    group_by(name, scenario, div, PID) %>% 
    summarize(value = mean(value))
) %>% 
  ggplot() +
  geom_boxplot(aes(x=paste(scenario, div), y=value, col=paste(scenario,div))) + 
  facet_wrap(~name, scales="free_y")+
  scale_y_continuous(expand = expansion(mult = 0.5))+
  theme_bw()
dev.off()

## Filter and write outputs

data_proc %>%
  mutate(D = purrr::map(.x = D, .f = function(x){
    x %>% filter(YEAR > max(YEAR)-100)
    }
    )) %>%
  mutate(Y = purrr::map(.x = Y, .f = function(x){
    x %>% filter(YEAR > max(YEAR)-100)
  }
  )) %>%
  pivot_longer(D:Y, names_to ="dy") %>% 
  mutate(filename = paste0(
    sitename, "_", 
    dy, "_PFATE_",
    scenario, "_",
    case_match(div,
               "evol"~"HD",
               "ld"~"LD"
    ),
    ".csv"
  )) %>% 
  mutate(file = here::here("mip_data_submitted", filename)) %>%
  group_by(file) %>%
  do(a = readr::write_csv(as.data.frame(.$value), file = .$file))


# 
# ## Write long-term outputs
# data_proc %>%
#   select(filename_prefix, processed_eco, processed_ind) %>%
#   pivot_longer(processed_eco:processed_ind, names_sep = "_", names_to = c("blurb", "level")) %>%
#   mutate(filename = paste0(filename_prefix, "_", level, ".csv")) %>%
#   mutate(file = here::here("data_submitted_v3", "long_term", filename)) %>%
#   group_by(file) %>%
#   do(a = write_csv(as.data.frame(.$value), file = .$file))
# 
# 

## Traits for LD run

data_proc %>% 
  select(scenario, div, Y) %>% 
  unnest(Y) %>% 
  filter(scenario == "AMB") %>% 
  filter(div == "evol") %>% 
  filter(YEAR > 1900 & YEAR < 2000) %>% 
  select(YEAR, PID, MH, WD, P50) %>% 
  pivot_longer(-c(YEAR, PID)) %>% 
  group_by(PID, name) %>% 
  summarize(value = mean(value))
