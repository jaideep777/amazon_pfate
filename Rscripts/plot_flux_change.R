library(tidyverse)

input_dir = here::here("input_data/")

output_dir = here::here("pfate_output_co2scan/")
expt_dir = "AmzMIP_HIST_614.2_evol_20ky_4"

source("~/codes/Plant-FATE/R/process_outputs.R")
source(here::here("Rscripts/definitions.R"))

l1 = pf_read_outputs(output_dir, expt_dir) %>% 
  pf_slice_time(1900, 2100)

dat_flux = l1$dat_d %>% 
  select(YEAR, GPP, NPP, GS, VCMAX, MORT) %>% 
  mutate(GPP = GPP*365.2425,
         NPP = NPP*365.2425,
         MORT = MORT*365.2425) %>% 
  pivot_longer(-YEAR) %>% 
  # Create facet labels
  mutate(name = factor(name, levels = unique(name), labels = labels[unique(name)]))

dat_z = l1$Zp %>% select(YEAR, V2:V6) %>%
  setNames(c("YEAR", paste0("z", 1:5))) %>% 
  pivot_longer(-YEAR, names_to="level") %>%
  mutate(name=factor("Z", levels="Z", labels=labels["Z"]))
         
cairo_pdf(here::here("figures/flux_change.pdf"), width=5.5, height=4.5)            
print(
dat_flux %>% 
  ggplot(aes(x=YEAR, y=value))+
  geom_line(alpha=0.3)+
  geom_smooth(method = "loess", span=0.15, se=F, col="black", linewidth=0.5)+
  geom_line(data = dat_z,
            aes(colour=level, group=level))+
  geom_rect(data = . %>% group_by(name) %>% slice(1), aes(xmin = 2020, xmax=2100, ymin=-Inf, ymax=Inf), fill="yellow", alpha=0.1)+
  geom_rect(data = dat_z %>% group_by(name) %>% slice(1), aes(xmin = 2020, xmax=2100, ymin=-Inf, ymax=Inf), fill="yellow", alpha=0.1)+
  facet_wrap(~name, scales="free_y", 
             strip.position = "left",
             labeller = label_parsed, nrow=3)+
  scale_colour_viridis_d(direction = -1, end=0.95)+
  labs(y="", x="Year", colour="Canopy\nlevel")+
  theme_bw()+
  amz_theme()+
  theme(strip.placement = "outside",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())
)
dev.off()

l1$dat_d %>% 
  select(YEAR, GPP, NPP, GS, VCMAX, MORT) %>% 
  mutate(GPP = GPP*365.2425,
         NPP = NPP*365.2425,
         MORT = MORT*365.2425) %>% 
  pivot_longer(-YEAR) %>% 
  bind_rows(
    l1$dat3 %>% 
      mutate(AGB = CL+CW,
             BGB = CCR+CFR,
             BA = BA*1e4) %>% 
      select(YEAR, BA, AGB, BGB, LAI) %>% 
      pivot_longer(-YEAR)
  ) %>% 
  mutate(hist = cut(YEAR, breaks = c(1900,1970,2000,2070,2100), labels = c("prehist","hist", "mid", "ele"))) %>%
  group_by(hist, name) %>% 
  summarize(value = mean(value)) %>% 
  pivot_wider(names_from = hist) %>% 
  mutate(pc_change = (ele-hist)/hist*100) %>% 
  write_csv(here::here("summarized_outputs/flux_change.csv"))



  