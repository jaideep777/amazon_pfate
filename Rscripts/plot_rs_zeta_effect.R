rm(list=ls())
library(tidyverse)

input_dir = here::here("input_data/")

source("~/codes/Plant-FATE/R/process_outputs.R")
source(here::here("Rscripts/definitions.R"))


read_supp = function(output_dir, expt_dir, ymin, ymax){

  l1 = pf_read_outputs(output_dir, expt_dir)
  
  l1$dat_d %>% 
    filter(YEAR > ymin & YEAR <= ymax) %>% 
    select(-IYEAR, -MON, -DAY) %>% 
    pivot_longer(-YEAR) %>% 
    group_by(name) %>% 
    summarize(value = mean(value)) %>% 
    bind_rows(
      l1$dat3 %>% 
        mutate(AGB = CL+CW) %>% 
        filter(YEAR > ymin & YEAR <= ymax) %>% 
        pivot_longer(-YEAR) %>% 
        group_by(name) %>% 
        summarize(value = mean(value))    
    ) %>% 
    bind_rows(
      l1$traits %>% 
        filter(!grepl("probe", SPP)) %>% 
        select(-SPP) %>% 
        filter(YEAR > ymin & YEAR <= ymax) %>% 
        pivot_longer(-YEAR) %>% 
        group_by(name) %>% 
        summarize(value = mean(value))    
    ) %>% 
    bind_rows(
      l1$Zp %>% select(YEAR, V2:V6) %>%
        setNames(c("YEAR", paste0("z", 1:5))) %>% 
        filter(YEAR > ymin & YEAR <= ymax) %>% 
        pivot_longer(-YEAR) %>% 
        group_by(name) %>% 
        summarize(value=mean(value, na.rm=T))
    )
}

list = tibble(
output_dir = c(
  here::here("pfate_output_co2scan/"),
  here::here("pfate_output_supplementary/"),
  here::here("pfate_output_supplementary/")
),
expt_dir = c(
  "AmzMIP_HIST_614.2_evol_20ky_4",
  "AmzMIP_ELE_614.2_evol_20ky_4_plusZeta",
  "AmzMIP_ELE_614.2_evol_20ky_4_plusRs"
  )
) %>% 
  pmap(.f=read_supp, ymin=19900, ymax=20000)

list[[length(list)+1]] = read_supp( here::here("pfate_output_co2scan/"),
                    "AmzMIP_HIST_614.2_evol_20ky_4",
                    1900, 2000)

df = list %>% plyr::join_all(by="name") %>% 
  setNames(c("name", "eco2", "eco2.zeta", "eco2.rs", "aco2")) %>%
  filter(name %in% c("GPP", "NPP", "RAU", "LAI", "AGB", "CFR", "BA", "WD", "P50X", "HMAT", "z2", "z3")) %>% 
  mutate(diff_co2 = (eco2-aco2)/aco2*100,
         diff_rs = (eco2.rs-aco2)/aco2*100,
         diff_zeta = (eco2.zeta-aco2)/aco2*100)

df %>% write_csv(here::here("summarized_outputs/rs_zeta_comparison.csv"))

cairo_pdf(here::here("figures/eco2_rs_zeta.pdf"), width = 6.5, height = 4.5)
print(
df %>%
  mutate(name = factor(name, levels=unique(name), labels=paste0("**",letters[1:length(unique(name))], "**. ", labels_nounit[unique(name)]))) %>%
  select(name, diff_co2, diff_rs, diff_zeta) %>%
  pivot_longer(-name, names_to = "scenario") %>% 
  ggplot()+
  geom_col(aes(y=value, x=scenario, fill=scenario))+
  facet_wrap(~name, scales="free_y", axes = "margins")+
  labs(y="Percent change", x="Scenario")+
  scale_x_discrete(labels = c("+CO<sub>2</sub>",
                              "+CO<sub>2</sub>,<br>+Rs", 
                              "+CO<sub>2</sub>,<br>+&zeta;"))+
  scale_fill_manual(values = c(
    diff_co2  = col_ele, 
    diff_rs   = "plum", 
    diff_zeta = "aquamarine4"
  ))+
  guides(fill="none")+
  amz_theme()+
  theme(axis.text.x = ggtext::element_markdown(lineheight=1.0))
)
dev.off()



