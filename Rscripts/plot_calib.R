rm(list=ls())
library(tidyverse)

input_dir = here::here("input_data/")
output_dir = here::here("pfate_output_co2scan/")

## Final runs (MIP)
expt_dir = "AmzMIP_HIST_614.2_evol_20ky_4"

source("~/codes/Plant-FATE/R/process_outputs.R")
source(here::here("Rscripts/definitions.R"))

traits_obs = read.csv(file = paste0(input_dir, "/Amz_trait_orig.csv"))

calib_fluxes = 
  read.csv(text = gsub(pattern = " ", replacement="",
                       "name, mean, min, max, type 
      GPP,  NA, 3,  3.5, fluxes     
      NPP, 1.31, NA, NA, fluxes    
      AGB, NA, 16.9, 20.7, structure  
      GS, 0.16, NA, NA, fluxes     
      LAI, NA, 5.3, 6.2, structure    
      CFR, NA, 0.48, 0.66, structure  
      VCMAX, NA, 39.3, 45.7, fluxes"),
           header=T, sep=",", as.is = F) %>% 
  as_tibble()

calib_ba = 
  traits_obs %>% 
  select(Species, Total.BasalArea_2017.cm2., Height_Max.m., meanWoodDensity..g.cm3., Leaf.LMA..g.m2., P50..Mpa.) %>% 
  mutate(BA.m2.ha = Total.BasalArea_2017.cm2.*1e-4*1e4/(8*pi*(30/2)^2)) %>% 
  summarize(
    BA = sum(BA.m2.ha, na.rm=T),
  ) %>% 
  pivot_longer(everything(), values_to = "mean") %>% 
  mutate(type="structure")


calib = calib_fluxes %>% bind_rows(calib_ba)


l1 = pf_read_outputs(output_dir, expt_dir) %>% 
       pf_slice_time(1900, 2100)

pred_ts = l1$dat_d %>% 
  select(YEAR, GPP, NPP, GS, VCMAX) %>% 
  filter(YEAR >= 1999 & YEAR <= 2020) %>%
  mutate(
    GPP = GPP*365.2425,
    NPP = NPP*365.2425,
  ) %>% 
  pivot_longer(-YEAR) %>% 
  bind_rows(
    l1$dat3 %>% 
      select(YEAR, CL, CW, CCR, CFR, BA, LAI, TB) %>% 
      filter(YEAR >= 1999 & YEAR <= 2020) %>%
      mutate(
        AGB = CL+CW, 
        BGB = CCR+CFR,
        BA = BA*1e4) %>% 
      pivot_longer(-YEAR)
  ) %>%  
  bind_rows(
    l1$traits %>% 
      filter(YEAR >= 1999 & YEAR <= 2020) %>%
      filter(!grepl("probe", SPP)) %>% 
      select(YEAR, LMA, WD, HMAT, P50X) %>% 
      pivot_longer(-YEAR)
  ) %>% 
  bind_rows(data.frame(YEAR=NA, name="BA", value=20)) %>% 
  bind_rows(data.frame(YEAR=NA, name="AGB", value=10)) %>% 
  inner_join(calib) %>% 
  mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)]))
  
pred = pred_ts %>% 
  filter(YEAR >= 2000 & YEAR <= 2015) %>% 
  filter(!is.na(YEAR)) %>% 
  group_by(name, type) %>% 
  summarize(across(everything(), ~mean(.))) %>% 
  ungroup()

calib_yr = 2007



cairo_pdf(file=here::here("figures/calib_fluxes_structure.pdf"), width = 7, height = 3.5)
print(
pred %>% 
  ggplot() + 
  # geom_col(aes(y=value, x=calib_yr), width = 1.5, fill="grey")+
  scale_y_continuous(expand = expansion(mult=0.4))+
  geom_line(data = pred_ts,
            aes(x=YEAR, y=value), col=col_amb, linewidth=0.3) + 
  geom_point(aes(y=mean, x=calib_yr), pch=1, col=col_obs, size=3, stroke=1.2)+
  geom_errorbar(aes(x=calib_yr, y=mean, ymin=min, ymax=max), col=col_obs, width = 2.5, linewidth=0.8)+
  geom_point(aes(y=value, x=calib_yr), col="grey20")+
  facet_wrap(~name, 
             scales="free_y", 
             strip.position = "left",
             ncol=4, axes = "margins")+
  geom_label(data = . %>% count(name) %>% mutate(label = letters[row_number()]),
             aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  labs(x = "Year", y = "")+
  scale_x_continuous(n.breaks = 3)+
  amz_theme()+
  theme(plot.margin = margin(t = 5, b=5, r = 10, unit = "pt"))
)
dev.off()


gauss_mix = function(x, means, sds, wts){
  y = x*0
  if (length(sds)==1) sds = rep(sds, length(means))
  for (i in 1:length(means)){
    y = y + wts[i]*dnorm(x, mean = means[i], sd=sds[i])
  }
  y
}

l = pf_read_outputs(output_dir, expt_dir)

dist_amb1 = l$dist %>% 
  mutate(period = "MID") %>% 
  mutate(period = ifelse(YEAR > 2000 & YEAR < 2020, yes="AMB", no=period)) %>% 
  mutate(period = ifelse(YEAR > 19980 & YEAR < 20000, yes="ELE", no=period)) %>% 
  filter(period %in% c("AMB", "ELE")) %>% 
  select(-X3) %>% 
  pivot_longer(cols=-c(YEAR,SPP,period), names_to="size_class") %>% 
  # sum over species
  group_by(YEAR,size_class, period) %>% 
  summarize(de = sum(value, na.rm=T)) %>% 
  # mean over years
  group_by(size_class, period) %>% 
  summarize(density=mean(de)) %>% 
  mutate(density = density*1e-2*1e4) %>% 
  mutate(size = l$x[as.numeric(sub('.','',size_class))-3]) %>% 
  arrange(size)
  

# Data for Manaus from https://link.springer.com/article/10.1007/s00442-004-1598-z 
dist_obs = data.frame(
  xobs = c(15,25,35,45,55,65,75,85,95,105)/100,
  yobs=c(350.5221340921042,
         132.41927918860426,
         62.62503296462008,
         29.61724892214378,
         15.095996574802413,
         5.702923697662178,
         2.3219542502889836,
         1.5968055466971947,
         0.7006940913385968,
         0.5597156879584093)/10
)

p1 = dist_amb1 %>% 
  ggplot() +
  geom_line(aes(x=size, y=density, group=period, col=period), linewidth=0.8)+
  scale_y_log10(limits=c(1e-3, 1000))+
  xlim(c(0.01, 1.2))+
  geom_point(data=dist_obs, aes(x=xobs, y=yobs), shape = 21, col=col_obs, fill=alpha(col_obs, 0.2), size=2)+
  geom_label(data = tibble(label="a"),
             aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  xlab("Diameter<br>(m)")+
  ylab("Density<br>(stems cm<sup>&minus;1</sup> ha<sup>&minus;1</sup>)")+
  scale_colour_manual(values = c(AMB=col_amb, ELE=col_ele))+
  theme_bw()+
  amz_theme()

p1

year_sq = 2000

df_trait = l$dat2 %>% 
  mutate(period = "MID") %>% 
  mutate(period = ifelse(YEAR > 2000 & YEAR < 2020, yes="AMB", no=period)) %>% 
  mutate(period = ifelse(YEAR > 19980 & YEAR < 20000, yes="ELE", no=period)) %>% 
  filter(period %in% c("AMB", "ELE")) %>% 
  select(YEAR, PID, BA, period) %>% 
  left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(!grepl("probe", PID)) %>% 
  mutate(YEAR = as.integer(YEAR)) %>%
  group_by(PID, period) %>% 
  summarize(across(everything(), ~mean(.)))

# %>%
#   mutate(yeardiff = abs(YEAR-year_sq)) %>% 
#   filter(yeardiff == min(yeardiff))
# 
# p2 = tibble(x = seq(0,300, length.out=1000),
#        y_obs = traits_obs %>% 
#          select(Leaf.LMA..g.m2., Total.BasalArea_2017.cm2.) %>% 
#          drop_na %>% 
#          mutate(means =Leaf.LMA..g.m2., 
#                 wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.)) %>% 
#          with(gauss_mix(x=x, means, wts, sds=180*0.1)),
#        y_pred_amb = df_trait %>% 
#          filter(period == "AMB") %>% 
#          with(gauss_mix(x=x, means =LMA*1000, wts=BA/sum(BA), sds=180*0.14)),
#        y_pred_ele = df_trait %>% 
#          filter(period == "ELE") %>% 
#          with(gauss_mix(x=x, means =LMA*1000, wts=BA/sum(BA), sds=180*0.14))
#        ) %>% 
#   ggplot(aes(x=x))+
#   geom_line(aes(y=y_obs), col="grey50")+
#   geom_ribbon(aes(ymax=y_obs, ymin=0), fill="grey", alpha=0.5)+
#   geom_line(aes(y=y_pred_amb))+
#   geom_line(aes(y=y_pred_ele), col="orange2")+
#   geom_label(data = tibble(label="b"),
#              aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
#   theme_bw()+
#   amz_theme()+
#   labs(y="Density", x=labels2["LMA"])
#   

p3 = tibble(x = seq(200,1200, length.out=1000),
       y_obs = traits_obs %>% 
         select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% 
         drop_na %>% 
         mutate(means =meanWoodDensity..g.cm3.*1000, 
                wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.)) %>% 
         with(gauss_mix(x=x, means, wts, sds=800*0.1)),
       y_pred_amb = df_trait %>% 
         filter(period == "AMB") %>% 
         with(gauss_mix(x=x, means =WD, wts=BA/sum(BA), sds=800*0.14)),
       y_pred_ele = df_trait %>% 
         filter(period == "ELE") %>% 
         with(gauss_mix(x=x, means =WD, wts=BA/sum(BA), sds=800*0.14))
  ) %>% 
  ggplot(aes(x=x))+
  geom_line(aes(y=y_obs), col=col_obs)+
  geom_ribbon(aes(ymax=y_obs, ymin=0), fill=col_obs, alpha=0.2)+
  geom_line(aes(y=y_pred_amb), col=col_amb, linewidth=0.8)+
  geom_line(aes(y=y_pred_ele), col=col_ele, linewidth=0.8)+
  geom_label(data = tibble(label="b"),
             aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  amz_theme()+
  labs(y="Density", x=labels2["WD"])

p4 = tibble(x = seq(0,50, length.out=1000),
       y_obs = traits_obs %>% 
         select(Height_Max.m., Total.BasalArea_2017.cm2.) %>% 
         drop_na %>% 
         mutate(means =Height_Max.m., 
                wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.)) %>% 
         with(gauss_mix(x=x, means, wts, sds=25*0.1)),
       y_pred_amb = df_trait %>% 
         filter(period == "AMB") %>% 
         with(gauss_mix(x=x, means =HMAT, wts=BA/sum(BA), sds=25*0.14)),
       y_pred_ele = df_trait %>% 
         filter(period == "ELE") %>% 
         with(gauss_mix(x=x, means =HMAT, wts=BA/sum(BA), sds=25*0.14))
  ) %>% 
  ggplot(aes(x=x))+
  geom_line(aes(y=y_obs), col=col_obs)+
  geom_ribbon(aes(ymax=y_obs, ymin=0), fill=col_obs, alpha=0.2)+
  geom_line(aes(y=y_pred_amb), col=col_amb, linewidth=0.8)+
  geom_line(aes(y=y_pred_ele), col=col_ele, linewidth=0.8)+
  geom_label(data = tibble(label="c"),
             aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  theme_bw()+
  amz_theme()+
  labs(y="Density", x=labels2["HMAT"])

p5 = tibble(x = seq(-6,-0.1, length.out=1000),
       y_obs = traits_obs %>% 
         select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% 
         drop_na %>% 
         mutate(means = P50..Mpa., 
                wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.)) %>% 
         with(gauss_mix(x=x, means, wts, sds=2*0.31)),
       y_pred_amb = df_trait %>% 
         filter(period == "AMB") %>% 
         with(gauss_mix(x=x, means =P50X, wts=BA/sum(BA), sds=2*0.31)),
       y_pred_ele = df_trait %>% 
         filter(period == "ELE") %>% 
         with(gauss_mix(x=x, means =P50X, wts=BA/sum(BA), sds=2*0.31))
  ) %>% 
  ggplot(aes(x=x))+
  geom_line(aes(y=y_obs), col=col_obs)+
  geom_ribbon(aes(ymax=y_obs, ymin=0), fill=col_obs, alpha=0.2)+
  geom_line(aes(y=y_pred_amb), col=col_amb, linewidth=0.8)+
  geom_line(aes(y=y_pred_ele), col=col_ele, linewidth=0.8)+
  geom_label(data = tibble(label="d"),
             aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  theme_bw()+
  amz_theme()+
  labs(y="Density", x=labels2["P50X"])


library(patchwork)

cairo_pdf(file=here::here("figures/calib_traits_sizedist.pdf"), width = 6, height = 5)
print(
(p1 + p3) / (p4 + p5) + plot_layout(guides="collect") 
)
dev.off()
