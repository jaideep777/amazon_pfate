rm(list=ls())
library(tidyverse)

input_dir = here::here("input_data/")

output_dir = here::here("pfate_output_evol_vs_comm/")
expt_dir = "AMB_wd_comm_20ky"

source("~/codes/Plant-FATE/R/process_outputs.R")
source(here::here("Rscripts/definitions.R"))

l = pf_read_outputs(output_dir, expt_dir)

n_species = l$dat2 %>% filter(!grepl("probe", .$PID)) %>% pull(PID) %>% unique() %>% length()
n_year = length(unique(l$dat2$YEAR))
col_species = rainbow(n = n_species, start = 0, end = 0.85, alpha = min(10/n_species, 1))
filter_span = span=30/n_year

cairo_pdf(here::here("figures/succesional_dynamics.pdf"), width = 6, height=2.25)
print(
l$dat2 %>% 
  filter(!grepl("probe", .$PID)) %>% 
  filter(YEAR > -20000 & YEAR < -19500) %>% 
  left_join(l$traits %>% rename(PID=SPP)) %>% 
  select(YEAR, PID, SEEDS, WD, BA) %>% 
  mutate(BA = BA*1e4,  # Convert m-2 to ha-1
         SEEDS = SEEDS*365.2425  # convert day-1 to yr-1
         ) %>% 
  pivot_longer(c(SEEDS,BA)) %>% 
  mutate(name = factor(name, levels=unique(name), labels=labels2[unique(name)])) %>% 
  ggplot(aes(x=YEAR, y=value, group=PID, color=WD))+
  stat_smooth(geom = "line", span=0.1, se=F, alpha=0.7, linewidth=0.4)+
  # geom_line(alpha=0.1)+
  facet_wrap(~name, scales="free_y", strip.position="left")+
  scale_x_continuous(n.breaks=3)+
  labs(y="", x="Year", col=labels3["WD"])+
  amz_theme()+
  theme(legend.title = ggtext::element_markdown())
)
dev.off()





