rm(list=ls())
library(tidyverse)

input_dir = here::here("input_data/")

output_dir = here::here("pfate_output_co2scan/")
expt_dir = "AmzMIP_HIST_614.2_evol_20ky_4"

source("~/codes/Plant-FATE/R/process_outputs.R")
source(here::here("Rscripts/definitions.R"))

l = pf_read_outputs(output_dir, expt_dir)

l1 = l %>% 
  pf_slice_time(1900, 2100)

dat_flux = l1$dat_d %>% 
  select(YEAR, GPP, NPP, GS, VCMAX, RAU, MORT) %>% 
  mutate(GPP = GPP*365.2425,
         NPP = NPP*365.2425,
         RAU = RAU*365.2425,
         MORT = MORT*365.2425) %>% 
  pivot_longer(-YEAR) %>%
  bind_rows(
    l1$dat3 %>% 
      mutate(AGB=CL+CW) %>% 
      select(YEAR, LAI, AGB) %>% 
      pivot_longer(-YEAR)
  ) %>% 
  # Create facet labels
  mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)]))

dat_z = l1$Zp %>% select(YEAR, V2:V6) %>%
  setNames(c("YEAR", paste0("z", 1:5))) %>% 
  pivot_longer(-YEAR, names_to="level") %>%
  mutate(name=factor("Z", levels="Z", labels=labels2["Z"]))

dat_traits = l$traits %>% 
  filter(!grepl("probe", SPP)) %>% 
  select(YEAR, SPP, WD, HMAT, P50X) %>% 
  pivot_longer(-c(YEAR, SPP)) %>% 
  mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)]))

p_flux = dat_flux %>% 
    ggplot(aes(x=YEAR, y=value))+
    geom_line(alpha=0.3)+
    geom_smooth(method = "loess", span=0.15, se=F, col="black", linewidth=0.5)+
    geom_line(data = dat_z,
              aes(colour=level, group=level))+
    geom_rect(data = . %>% group_by(name) %>% slice(1), aes(xmin = 2020, xmax=2100, ymin=-Inf, ymax=Inf), fill=col_ele, alpha=0.1)+
    geom_rect(data = dat_z %>% group_by(name) %>% slice(1), aes(xmin = 2020, xmax=2100, ymin=-Inf, ymax=Inf), fill=col_ele, alpha=0.1)+
    geom_label(data = . %>% count(name) %>% mutate(label = letters[1:n()]),
               aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
    geom_label(data = dat_z %>% count(name) %>% mutate(label = letters[1:n() + (dat_flux %>% pull(name) %>% unique() %>% length())] ),
               aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
    facet_wrap(~name, scales="free_y", 
               strip.position = "left",
               ncol=3)+
    scale_colour_viridis_d(direction = -1, end=0.95)+
    scale_x_continuous(n.breaks=3)+
    labs(y="", x="", colour="Canopy\nlevel")+
    amz_theme()


p_traits = dat_traits %>% 
  mutate(YEAR=as.integer(YEAR)) %>% 
  filter(YEAR %in% seq(-20000,20000, by=10)) %>% 
  ggplot(aes(x=YEAR, y=value))+
  geom_line(alpha=0.3)+
  geom_smooth(method = "loess", span=0.15, se=F, col="black", linewidth=0.5)+
  geom_rect(data = . %>% group_by(name) %>% slice(1), aes(xmin = 2020, xmax=Inf, ymin=-Inf, ymax=Inf), fill=col_ele, alpha=0.1)+
  geom_label(data = . %>% count(name) %>% mutate(label = letters[row_number()+9]),
             aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  facet_wrap(~name, scales="free_y", 
             strip.position = "left",
             ncol=3)+
  scale_x_continuous(n.breaks=3)+
  labs(y="", x="Year")+
  amz_theme()

cairo_pdf(here::here("figures/flux_change.pdf"), width=7, height=6.4)            
print(
  p_flux/p_traits + plot_layout(heights=c(3,1))
)
dev.off()


l$dat_d %>% 
  select(YEAR, GPP, NPP, GS, VCMAX, RAU, MORT) %>% 
  mutate(GPP = GPP*365.2425,
         NPP = NPP*365.2425,
         MORT = MORT*365.2425) %>% 
  pivot_longer(-YEAR) %>% 
  bind_rows(
    l$dat3 %>% 
      mutate(AGB = CL+CW,
             BGB = CCR+CFR,
             BA = BA*1e4) %>% 
      select(YEAR, BA, AGB, BGB, LAI) %>% 
      pivot_longer(-YEAR)
  ) %>% 
  mutate(hist = cut(YEAR, 
                    breaks = c(-Inf, 1970, 2000,   2070,    2100,     19900,     Inf), 
                    labels = c("prehist",  "hist", "mid", "ele_st", "midele", "ele_lt")
                    )
         ) %>%
  group_by(hist, name) %>% 
  summarize(value = mean(value)) %>% 
  pivot_wider(names_from = hist) %>% 
  mutate(pc_change_st = (ele_st-hist)/hist*100,
         pc_change_lt = (ele_lt-hist)/hist*100) %>% 
  write_csv(here::here("summarized_outputs/flux_change.csv"))



  