rm(list=ls())
library(tidyverse)

input_dir = here::here("input_data/")
output_dir = here::here("pfate_output_evol_vs_comm/")


source("~/codes/Plant-FATE/R/process_outputs.R")
source(here::here("Rscripts/definitions.R"))

expt_dir = "AMB_ELE_wd_evol_20ky"

l1 = pf_read_outputs(output_dir, expt_dir) %>% 
  pf_subsample_outputs(interval = 20)

expt_dir = "AMB_wd_comm_20ky"

l_amb = pf_read_outputs(output_dir, expt_dir) %>% 
  pf_subsample_outputs(interval = 20)

expt_dir = "ELE_wd_comm_20ky"

l_ele = pf_read_outputs(output_dir, expt_dir) %>% 
  pf_subsample_outputs(interval = 20)

df_amb = l_amb$dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(l_amb$traits %>% 
              select(YEAR, SPP, LMA:P50X), 
            by=c("PID"="SPP", "YEAR"="YEAR")) 

df_evol = l1$dat2 %>% select(YEAR, PID, BA) %>% 
  left_join(l1$traits %>% 
              select(YEAR, SPP, LMA:P50X), 
            by=c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(YEAR < 2020)


p1 = df_amb %>% mutate(type="CWM") %>% 
  bind_rows(
    df_evol %>% mutate(type="ESS")
  ) %>% 
  group_by(type, YEAR) %>% 
  summarize(WD = sum(WD*BA)/sum(BA)) %>% 
  ggplot(aes(x=YEAR, y=WD, col=type, group=type))+
  geom_line()+
  scale_x_continuous(n.breaks = 3)+
  amz_theme()+
  labs(y=paste("Mean", labels2["WD"]) %>% str_to_sentence(),
       x = "Year", col="Type")

p2 = data.frame(
  co2 = c(414, 828),
  cwm = c(814, 645),
  evol = c(803, 613)
) %>% 
  ggplot(aes(x=cwm, y=evol, col=factor(co2)))+
  geom_point()+
  geom_abline(slope = 1, intercept = 0, col=col_obs, alpha=0.5)+
  labs(y=paste("ESS", labels2["WD"] %>% str_to_lower()),
       x=paste("CWM", labels2["WD"] %>% str_to_lower()),
       col="CO<sub>2</sub>")+
  scale_colour_manual(values = c(`414`=col_amb, `828`=col_ele))+
  theme_bw()+
  amz_theme()+
  scale_x_continuous(expand=expansion(mult=0.5))+
  scale_y_continuous(expand=expansion(mult=0.5))+
  theme(legend.title = ggtext::element_markdown(lineheight=1.2))

library(patchwork)

cairo_pdf(here::here("figures/evol_vs_comm.pdf"), width=6, height=2.6)
print(
  p1+p2 + plot_layout(guides="collect")
)
dev.off()



