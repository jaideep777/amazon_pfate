rm(list=ls())
library(tidyverse)

source(here::here("Rscripts/definitions.R"))

input_dir = here::here("input_data/")
output_dir = data_path

expt_dir = "AmzMortality_AMB_rs0.5e-3m_inf0.0020_K_leaf0.5e-16_zeta0.2nspecies2" 

read_scan = function(expt_dir, ymin=-500, ymax=0){
  # expt_dir = paste0("AmzMIP_HIST_",co2,"_evol_20ky_4")
  wd_back = getwd()
  setwd(paste0(output_dir,"/",expt_dir))

  dat3 = readr::read_csv("Y_mean_PFATE.csv") |> 
    dplyr::mutate(YEAR = as.integer(YEAR))

  df = dat3 |> 
    mutate(AGB = CL + CW,
            BGB = CCR + CFR,
            RootShoot = CFR/CL) |>
    filter(YEAR >= ymin & YEAR <= ymax) |> 
    colMeans()
  
  setwd(wd_back)  
  
  df    
}
  

## Final runs (minf scan)
df_zeta <- tibble(zeta = seq(0.2, 0.3, by = 0.05)) %>%
  mutate(expt_dir = sprintf("AmzMortality_AMB_rs0.5e-3m_inf0.0020_K_leaf0.5e-16_zeta%gnspecies2", zeta)) |>
  mutate(data = purrr::map(expt_dir, ~as_tibble_row(read_scan(.x, ymin = -500, ymax = 0)))) %>%
  tidyr::unnest(cols = data)

df_zeta1 = df_zeta |> 
  pivot_longer(-(zeta:YEAR)) |> 
  select(-expt_dir) |>
  pivot_wider(names_from=zeta, values_from=value) |> 
  rename(zeta_0.2 = `0.2`) |> 
  pivot_longer(`0.25`:`0.3`, names_to = "zeta_new") |> 
  mutate(abs_change = value - zeta_0.2) |> 
  mutate(pc_change = (value/zeta_0.2-1)*100) |> 
  mutate(beta = log(value/zeta_0.2)/log(614/369)) |> 
  filter(name %in% c("AGB", "BGB", "RootShoot")) |> 
  filter(zeta_new == 0.25) 

df_zeta1 |> write.csv(here::here("summarized_outputs/zeta_change_betas.csv"), row.names = F)
