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

  dat = readr::read_csv("D_PFATE.csv") |> 
    dplyr::mutate(YEAR = as.integer(YEAR)) |> 
    dplyr::group_by(YEAR) |> 
    dplyr::summarize_all(mean)

  dat3 = readr::read_csv("Y_mean_PFATE.csv") |> 
    dplyr::mutate(YEAR = as.integer(YEAR))

  dat2 = readr::read_csv("Y_PFATE.csv") |> 
    dplyr::mutate(YEAR = as.integer(YEAR)) |> 
    filter(!grepl("probe", PID))

  traits = readr::read_csv("traits.csv") |>
    rename(PID = SPP) |>
    dplyr::mutate(YEAR = as.integer(YEAR)) |>
    filter(!grepl("probe", PID))

  traits_cwm <- traits |> 
    select(YEAR, PID, WD, P50X, HMAT) |> 
    left_join(
      dat2 |> select(YEAR, PID, BA)
    ) |> 
    group_by(YEAR) |> 
    summarize(
      WD = weighted.mean(WD, BA, na.rm=TRUE),
      P50X = weighted.mean(P50X, BA, na.rm=TRUE),
      HMAT = weighted.mean(HMAT, BA, na.rm=TRUE)
    )

  df = dat |>
    select(YEAR, GPP, MORT) |>
    left_join(
      dat3 |> 
        mutate(AGB = CL + CW) |>
        select(YEAR, AGB) 
    ) |>
    left_join(
      traits_cwm
    ) |> 
    filter(YEAR >= ymin & YEAR <= ymax) |> 
    colMeans()
  
  setwd(wd_back)  
  
  df    
}
  

## Final runs (minf scan)
df_minf <- tibble(m_inf = seq(0.0016, 0.0028, by = 0.0002)) %>%
  mutate(expt_dir = sprintf("AmzMortality_AMB_rs0.5e-3m_inf%0.4f_K_leaf0.5e-16_zeta0.2nspecies2", m_inf)) |>
  mutate(data = purrr::map(expt_dir, ~as_tibble_row(read_scan(.x, ymin = -500, ymax = 0)))) %>%
  tidyr::unnest(cols = data)

df_minf %>%
  write.csv(here::here("summarized_outputs/scan_minf_summary.csv"), row.names = F)

## ------------

df_minf = read.csv(here::here("summarized_outputs/scan_minf_summary.csv"))

p1 <- df_minf %>%
  select(MINF=m_inf, AGB, MORT, WD) %>%
  mutate(MORT = MORT*365.2425) %>%  # Convert kg m-2 day-1 --> kg m-2 yr-1
  pivot_longer(-MINF) %>% 
  mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)])) %>%
  ggplot(aes(x=MINF, y=value)) +
  geom_point(size=1)+
  geom_smooth(method = "lm", se=F, linewidth=0.5)+
  facet_wrap(~name, scales="free_y", 
             strip.position = "left") +
  labs(x = labels2["MINF"], y="")+
  geom_label(data = . %>% count(name) %>% mutate(label = letters[1:n()]),
              aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  amz_theme()+
  theme(axis.title.x = ggtext::element_markdown(size=10))

cairo_pdf(file=here::here("figures/minf_scan.pdf"), width = 6, height = 2)
print(
p1
)
dev.off()

## Final runs (minf scan)
df_zeta <- tibble(zeta = seq(0.1, 0.3, by = 0.05)) %>%
  mutate(expt_dir = sprintf("AmzMortality_AMB_rs0.5e-3m_inf0.0020_K_leaf0.5e-16_zeta%gnspecies2", zeta)) |>
  mutate(data = purrr::map(expt_dir, ~as_tibble_row(read_scan(.x, ymin = -500, ymax = 0)))) %>%
  tidyr::unnest(cols = data)

df_zeta %>%
  write.csv(here::here("summarized_outputs/scan_zeta_summary.csv"), row.names = F)

## ------------

df_zeta = read.csv(here::here("summarized_outputs/scan_zeta_summary.csv"))

p2 <- df_zeta %>%
  select(ZETA=zeta, AGB, MORT, WD) %>%
  mutate(MORT = MORT*365.2425) %>%  # Convert kg m-2 day-1 --> kg m-2 yr-1
  pivot_longer(-ZETA) %>% 
  mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)])) %>%
  ggplot(aes(x=ZETA, y=value)) +
  geom_point(size=1)+
  geom_smooth(method = "lm", se=F, linewidth=0.5)+
  facet_wrap(~name, scales="free_y", 
             strip.position = "left") +
  labs(x = labels2["ZETA"], y="")+
  geom_label(data = . %>% count(name) %>% mutate(label = letters[1:n()+3]),
              aes(x=-Inf, y=-Inf, label=label), inherit.aes = F, hjust=0, vjust=0, label.size = 0, size = 4.5) +
  amz_theme()+
  theme(axis.title.x = ggtext::element_markdown(size=10))

cairo_pdf(file=here::here("figures/zeta_scan.pdf"), width = 6, height = 2)
print(
p2
)
dev.off()

library(patchwork)
cairo_pdf(file=here::here("figures/scan_m_inf_zeta.pdf"), width = 7, height = 4.5)
print(
p1/p2
)
dev.off()



## Sapwood respiration rate visualization
ggplot(data=data.frame(x=c(exp(-3.75), exp(-1.25), exp(1.25))*1e-9*12*86400*365), aes(y=x,x=1))+geom_violin()+geom_point(x=1,y=0.4)

