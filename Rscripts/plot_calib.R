library(tidyverse)

input_dir = here::here("input_data/")
output_dir = here::here("pfate_output_co2scan/")

## Final runs (MIP)
expt_dir = "AmzMIP_HIST_414.2_evol_20ky_4"

source("~/codes/Plant-FATE/R/process_outputs.R")
source(here::here("Rscripts/definitions.R"))

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
  ) 

pred = pred_ts %>% 
  filter(YEAR >= 2000 & YEAR <= 2015) %>% 
  group_by(name) %>% 
  summarize(value = mean(value))

calib_yr = 2007

cairo_pdf(file=here::here("figures/calib_fluxes_structure.pdf"), width = 6, height = 6)
print(
pred %>% 
  inner_join(calib) %>% 
  mutate(name = factor(name, levels = unique(name), labels = labels[unique(name)])) %>% 
  ggplot() + 
  # geom_col(aes(y=value, x=calib_yr), width = 1.5, fill="grey")+
  scale_y_continuous(expand = expansion(mult=0.4))+
  theme_bw() + 
  theme(axis.line=element_line(colour="grey30"))+
  geom_line(data = pred_ts %>% 
              filter(name %in% calib$name) %>% 
              rbind(data.frame(YEAR=NA, name="BA", value=20)) %>% 
              rbind(data.frame(YEAR=NA, name="AGB", value=10)) %>% 
              left_join(calib) %>% 
              mutate(name = factor(name, levels = unique(name), labels = labels[unique(name)])),
            aes(x=YEAR, y=value), col="grey50", linewidth=0.3) + 
  geom_point(aes(y=mean, x=calib_yr), pch=1, col="dodgerblue1", size=3, stroke=1)+
  geom_errorbar(aes(x=calib_yr, y=mean, ymin=min, ymax=max), col="dodgerblue1", width = 2.5)+
  geom_point(aes(y=value, x=calib_yr))+
  facet_wrap(~factor(name, levels=labels[c("GPP", "NPP", "GS", "VCMAX", "AGB", "BA", "CFR", "LAI")]), 
             scales="free_y", 
             strip.position = "left",
             ncol=2, axes = "margins",
             labeller = label_parsed)+
  theme(strip.placement = "outside",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  labs(x = "Year", y = "")
)
dev.off()


gauss_mix = function(x, means, sds, wts, add=T, col, add_polygon=F, ...){
  y = x*0
  if (length(sds)==1) sds = rep(sds, length(means))
  for (i in 1:length(means)){
    y = y + wts[i]*dnorm(x, mean = means[i], sd=sds[i])
  }
  y
}

traits_obs = read.csv(file = paste0(input_dir, "/Amz_trait_orig.csv"))

l = l1

dist_amb1 = l$dist %>% 
  filter(YEAR > 2000 & YEAR < 2020) %>% 
  select(-X3) %>% 
  pivot_longer(cols=-(YEAR:SPP), names_to="size_class") %>% 
  # sum over species
  group_by(YEAR,size_class) %>% 
  summarize(de = sum(value, na.rm=T)) %>% 
  # mean over years
  group_by(size_class) %>% 
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
  geom_line(aes(x=size, y=density))+
  scale_y_log10(limits=c(1e-3, 1000))+
  xlim(c(0.01, 1.2))+
  geom_point(data=dist_obs, aes(x=xobs, y=yobs), col="grey30", alpha=0.4, size=2)+
  xlab("Diameter (m)")+
  ylab(exprlabel("Density", "(stems cm"^"-1"~"ha"^"-1"*")"))+
  theme_bw()+
  amz_theme()

df_trait = l$dat2 %>% 
  select(YEAR, PID, BA) %>% 
  left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>% 
  filter(!grepl("probe", PID)) %>% 
  mutate(YEAR = as.integer(YEAR)) %>%
  mutate(yeardiff = abs(YEAR-year_sq)) %>% 
  filter(yeardiff == min(yeardiff))

p2 = tibble(x = seq(0,300, length.out=1000),
       y_obs = traits_obs %>% 
         select(Leaf.LMA..g.m2., Total.BasalArea_2017.cm2.) %>% 
         drop_na %>% 
         mutate(means =Leaf.LMA..g.m2., 
                wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.)) %>% 
         with(gauss_mix(x=x, means, wts, sds=180*0.1, add=F, col="grey", add_polygon=T, type="l", lwd=1.5,  las=0, main="", xlab="LMA", ylab="density", ylim=c(0, 0.02))),
       y_pred = df_trait %>% 
         with(gauss_mix(x=x, means =LMA*1000, wts=BA/sum(BA), sds=180*0.14, add=T, col="black", type="l", lwd=1.5))
) %>% 
  ggplot(aes(x=x))+
  geom_line(aes(y=y_pred))+
  geom_line(aes(y=y_obs), col="grey50")+
  geom_ribbon(aes(ymax=y_obs, ymin=0), fill="grey", alpha=0.5)+
  theme_bw()+
  amz_theme()+
  labs(y="Density", x=labels["LMA"]$LMA)
  

p3 = tibble(x = seq(200,1200, length.out=1000),
       y_obs = traits_obs %>% 
         select(meanWoodDensity..g.cm3., Total.BasalArea_2017.cm2.) %>% 
         drop_na %>% 
         mutate(means =meanWoodDensity..g.cm3.*1000, 
                wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.)) %>% 
         with(gauss_mix(x=x, means, wts, sds=800*0.1, add=F, col="grey", add_polygon=T, type="l", lwd=1.5,  las=0, main="", xlab="LMA", ylab="density", ylim=c(0, 0.02))),
       y_pred = df_trait %>% 
         with(gauss_mix(x=x, means =WD, wts=BA/sum(BA), sds=800*0.14, add=T, col="black", type="l", lwd=1.5))
) %>% 
  ggplot(aes(x=x))+
  geom_line(aes(y=y_pred))+
  geom_line(aes(y=y_obs), col="grey50")+
  geom_ribbon(aes(ymax=y_obs, ymin=0), fill="grey", alpha=0.5)+
  theme_bw()+
  amz_theme()+
  labs(y="Density", x=labels["WD"]$WD)

p4 = tibble(x = seq(0,50, length.out=1000),
       y_obs = traits_obs %>% 
         select(Height_Max.m., Total.BasalArea_2017.cm2.) %>% 
         drop_na %>% 
         mutate(means =Height_Max.m., 
                wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.)) %>% 
         with(gauss_mix(x=x, means, wts, sds=25*0.1, add=F, col="grey", add_polygon=T, type="l", lwd=1.5,  las=0, main="", xlab="LMA", ylab="density", ylim=c(0, 0.02))),
       y_pred = df_trait %>% 
         with(gauss_mix(x=x, means =HMAT, wts=BA/sum(BA), sds=25*0.14, add=T, col="black", type="l", lwd=1.5))
) %>% 
  ggplot(aes(x=x))+
  geom_line(aes(y=y_pred))+
  geom_line(aes(y=y_obs), col="grey50")+
  geom_ribbon(aes(ymax=y_obs, ymin=0), fill="grey", alpha=0.5)+
  theme_bw()+
  amz_theme()+
  labs(y="Density", x=labels["HMAT"]$HMAT)

p5 = tibble(x = seq(-6,-0.1, length.out=1000),
       y_obs = traits_obs %>% 
         select(P50..Mpa., Total.BasalArea_2017.cm2.) %>% 
         drop_na %>% 
         mutate(means = P50..Mpa., 
                wts=Total.BasalArea_2017.cm2./sum(Total.BasalArea_2017.cm2.)) %>% 
         with(gauss_mix(x=x, means, wts, sds=2*0.31, add=F, col="grey", add_polygon=T, type="l", lwd=1.5,  las=0, main="", xlab="LMA", ylab="density", ylim=c(0, 0.02))),
       y_pred = df_trait %>% 
         with(gauss_mix(x=x, means =P50X, wts=BA/sum(BA), sds=2*0.31, add=T, col="black", type="l", lwd=1.5))
) %>% 
  ggplot(aes(x=x))+
  geom_line(aes(y=y_pred))+
  geom_line(aes(y=y_obs), col="grey50")+
  geom_ribbon(aes(ymax=y_obs, ymin=0), fill="grey", alpha=0.5)+
  theme_bw()+
  amz_theme()+
  labs(y="Density", x=labels["P50X"]$P50X)


library(patchwork)

cairo_pdf(file=here::here("figures/calib_traits_sizedist.pdf"), width = 6, height = 5)
print(
(p1 + p3) / (p4 + p5)
)
dev.off()
