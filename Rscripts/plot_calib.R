rm(list=ls())
library(tidyverse)

source("~/codes/Plant-FATE/R/process_outputs.R")
source(here::here("Rscripts/definitions.R"))

input_dir = here::here("input_data/")

## Final runs (manuscript)
output_dir = data_path
expt_dir = "AmzMortality_AMB_rs0.5e-3m_inf0.002_K_leaf0.5e-16_zeta0.2nspecies2timestep0.0416"

# ## Final runs (MIP)
# output_dir = here::here("pfate_output_mip")
# expt_dir = "AmzMIP_HIST_AMB_evol_20ky_c2_rs0.04"

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

l = pf_read_outputs(input_dir, output_dir, expt_dir)
l1 = l %>% pf_slice_time(1900, 2100)

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


p0 <- pred %>%
  ggplot() +
  # geom_col(aes(y=value, x=calib_yr), width = 1.5, fill="grey")+
  scale_y_continuous(expand = expansion(mult=0.4))+
  geom_rect(aes(ymin = min, ymax = max), xmin = -Inf, xmax = Inf, 
           fill = col_obs, alpha = 0.2)+
  # geom_errorbar(aes(x=calib_yr, y=mean, ymin=min, ymax=max), col=col_obs, width = 2.5, linewidth=0.8)+
  geom_hline(aes(yintercept = mean), col = col_obs, alpha=1)+
  geom_hline(aes(yintercept = min), col = col_obs, alpha=0.5, linewidth=0.3)+
  geom_hline(aes(yintercept = max), col = col_obs, alpha=0.5, linewidth=0.3)+
  # geom_point(aes(y=mean, x=calib_yr), shape=21, col=col_obs, fill = alpha(col_obs, 0.2), size=3, stroke=1.2)+
  # geom_point(aes(y=value, x=calib_yr, col="Predicted (mean)"))+
  geom_line(data = pred_ts,
            aes(x=YEAR, y=value, col="Predicted (timeseries)"), linewidth=0.3) +
  facet_wrap(~name,
             scales="free_y",
             strip.position = "left",
             ncol=4, axes = "margins")+
  geom_label(data = . %>% count(name) %>% mutate(label = letters[row_number()]),
             aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  labs(x = "Year", y = "")+
  scale_x_continuous(n.breaks = 3)+
  amz_theme()+
  scale_color_manual(values = c("Predicted (mean)"="grey20",
                                "Predicted (timeseries)" = col_amb))+
  theme(plot.margin = margin(t = 5, b=5, r = 10, unit = "pt"))

p0

# cairo_pdf(file=here::here("figures/calib_fluxes_structure.pdf"), width = 7, height = 3.5)
cairo_pdf(file=here::here("figures/calib_fluxes_structure_mip.pdf"), width = 7, height = 3.5)
print(p0)
dev.off()


gauss_mix = function(x, means, sds, wts){
  y = x*0
  if (length(sds)==1) sds = rep(sds, length(means))
  for (i in 1:length(means)){
    y = y + wts[i]*dnorm(x, mean = means[i], sd=sds[i])
  }
  y
}

# l = pf_read_outputs(input_dir, output_dir, expt_dir)

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


year_sq = 2000

df_trait = l$dat2 %>%
  mutate(period = "MID") %>%
  mutate(period = ifelse(YEAR > 2000 & YEAR < 2020, yes="AMB", no=period)) %>%
  mutate(period = ifelse(YEAR > 19980 & YEAR < 20000, yes="ELE", no=period)) %>%
  filter(period %in% c("AMB", "ELE")) %>%
  select(YEAR, PID, BA, period) %>%
  mutate(BA = BA*1e4) %>%  # convert m2 m-2 ---> m2 ha-1
  left_join(l$traits, by = c("PID"="SPP", "YEAR"="YEAR")) %>%
  filter(!grepl("probe", PID)) %>%
  mutate(YEAR = as.integer(YEAR)) %>%
  group_by(PID, period) %>%
  summarize(across(everything(), ~mean(.)))


df_td = dplyr::bind_rows(
  dist_amb1 %>% 
    ungroup() %>% 
    pivot_wider(names_from=period, values_from=density) %>% 
    rename(x=size, y_amb=AMB) %>% 
    # rename(x=size, y_amb=AMB, y_ele=ELE) %>% 
    select(-size_class) %>% 
    full_join(dist_obs %>% 
                 rename(x=xobs, y_obs=yobs)) %>% 
    mutate(name = "D") %>% 
    mutate(id = 1:n()),
  
  tibble(x = seq(200,1200, length.out=1000),
                 y_obs = traits_obs %>%
                   select(meanWoodDensity..g.cm3., BA) %>%
                   drop_na %>%
                   mutate(means =meanWoodDensity..g.cm3.*1000,
                          wts=BA/sum(BA)*sum(traits_obs$BA, na.rm=T)) %>%
                   with(gauss_mix(x=x, means, wts, sds=800*0.1)),
                 y_amb = df_trait %>%
                   filter(period == "AMB") %>%
                   with(gauss_mix(x=x, means =WD, wts=BA, sds=800*0.14))
                #  y_ele = df_trait %>%
                #    filter(period == "ELE") %>%
                #    with(gauss_mix(x=x, means =WD, wts=BA, sds=800*0.14))
  ) %>% 
    mutate(id=1:n()) %>% 
    mutate(name = "WD"),
  
  tibble(x = seq(0,50, length.out=1000),
         y_obs = traits_obs %>%
           select(Height_Max.m., BA) %>%
           drop_na %>%
           mutate(means =Height_Max.m.,
                  wts=BA/sum(BA)*sum(traits_obs$BA, na.rm=T)) %>%
           with(gauss_mix(x=x, means, wts, sds=25*0.1)),
         y_amb = df_trait %>%
           filter(period == "AMB") %>%
           with(gauss_mix(x=x, means =HMAT, wts=BA, sds=25*0.14))
        #  y_ele = df_trait %>%
        #    filter(period == "ELE") %>%
        #    with(gauss_mix(x=x, means =HMAT, wts=BA, sds=25*0.14))
  ) %>%
    mutate(id=1:n()) %>% 
    mutate(name = "HMAT"),
  
  tibble(x = seq(-6,-0.1, length.out=1000),
         y_obs = traits_obs %>%
           select(P50..Mpa., BA) %>%
           drop_na %>%
           mutate(means = P50..Mpa.,
                  wts=BA/sum(BA)*sum(traits_obs$BA, na.rm=T)) %>%
           with(gauss_mix(x=x, means, wts, sds=2*0.31)),
         y_amb = df_trait %>%
           filter(period == "AMB") %>%
           with(gauss_mix(x=x, means =P50X, wts=BA, sds=2*0.31))
        #  y_ele = df_trait %>%
        #    filter(period == "ELE") %>%
        #    with(gauss_mix(x=x, means =P50X, wts=BA, sds=2*0.31))
  ) %>%
    mutate(id=1:n()) %>% 
    mutate(name = "P50X"),
)   

p_td = list()
for (i in 1:3){
  label = letters[i+9]
  p_td[[i]] = 
    df_td %>% 
      filter(name == c("WD", "HMAT", "P50X")[i]) %>% 
      mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)])) %>% 
      ggplot(aes(x=x))+
      geom_line(aes(y=y_obs, col="Site observations"), linewidth=0.4)+
      geom_ribbon(aes(ymax=y_obs, ymin=0), fill=col_obs, linewidth = 0.4, alpha=0.2)+
      geom_line(aes(y=y_amb, col="Predicted (distribution)"), linewidth=0.8)+
      # geom_line(aes(y=y_ele), col=col_ele, linewidth=0.8)+
      amz_theme()+
      theme(axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)))+
      geom_label(data = . %>% count(name),
                 aes(x=-Inf, y=Inf), label=label, inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
      facet_wrap(~name, scales="free",
                 strip.position = "bottom",
                 ncol=3)+
      scale_x_continuous(n.breaks = 3)+
      scale_y_continuous(n.breaks = 4)+
      scale_color_manual(values=c(
        "Predicted (distribution)" = col_amb,
        "Site observations" = col_obs
      ))+
      labs(y="Basal-area density<br>(m<sup>2</sup> trait-unit<sup>&minus;1</sup> ha<sup>&minus;1</sup>)", x=NA)
}

p_sd = df_td %>% 
  filter(name == "D") %>% 
  mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)])) %>% 
  ggplot(aes(x=x))+
  geom_point(data=dist_obs, aes(x=xobs, y=yobs), shape = 21, col=col_obs, fill=alpha(col_obs, 0.2), size=2)+
  geom_line(aes(y=y_obs), col=col_obs)+
  geom_ribbon(aes(ymax=y_obs, ymin=0), fill=col_obs, alpha=0.2)+
  geom_line(data=. %>% filter(!is.na(y_amb)), aes(y=y_amb), col=col_amb, linewidth=0.8)+
  # geom_line(data=. %>% filter(!is.na(y_ele)), aes(y=y_ele), col=col_ele, linewidth=0.8)+
  amz_theme()+
  facet_wrap(~name, scales="free",
             strip.position = "bottom",
             ncol=3, axes = "margins")+
  scale_y_log10(limits=c(1e-3, 1000))+
  scale_x_continuous(n.breaks = 3, limits = c(0.01, 1.2))+
  geom_label(data = . %>% count(name) %>% mutate(label = letters[row_number()+8]),
             aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  labs(y="Density<br>(stems cm<sup>&minus;1</sup> ha<sup>&minus;1</sup>)", x=NA)



library(patchwork)

cairo_pdf(file=here::here("figures/calib_all_v3.1.pdf"), width = 8.5, height = 6.5)
print(
(p0 + theme(panel.spacing.y = unit(0.5, "cm"))) / (p_sd + p_td[[1]]+ plot_spacer() + p_td[[2]] + p_td[[3]] + plot_layout(nrow=1, widths=c(0.98,0.98,-0.03, 1,1))) + 
  plot_layout(guides="collect", heights =c(2,1)) & 
  theme(legend.position = 'top',
        axis.title.y = element_text(size=9)) &
  labs(color=NULL)
)
dev.off()


### SI Figures for trait and size-dist change

# cairo_pdf(file=here::here("figures/trait_shift_si.pdf"), width = 7, height = 3.5)
# print(
# df_td %>% 
#   filter(name != "D") %>% 
#   mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)])) %>% 
#   ggplot(aes(x=x))+
#   geom_line(aes(y=y_amb, col="AMB (CO2 = 414 ppm)"), linewidth=0.8)+
#   geom_line(aes(y=y_ele, col="ELE (CO2 = 614 ppm)") , linewidth=0.8)+
#   amz_theme()+
#   theme(axis.title.y = element_text(margin = margin(t = 0, r = 6, b = 0, l = 0)))+
#   geom_label(data = . %>% count(name) %>% mutate(label = letters[row_number()]),
#              aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
#   facet_wrap(~name, scales="free",
#              strip.position = "bottom",
#              ncol=3)+
#   scale_x_continuous(n.breaks = 3)+
#   scale_y_continuous(n.breaks = 4)+
#   scale_color_manual(values=c(
#     "AMB (CO2 = 414 ppm)" = col_amb,
#     "ELE (CO2 = 614 ppm)" = col_ele
#   ))+
#   theme(legend.position = "top")+
#   labs(y="Basal-area density<br>(m<sup>2</sup> trait-unit<sup>&minus;1</sup> ha<sup>&minus;1</sup>)", x=NA, color="")
# )
# dev.off()


# cairo_pdf(file=here::here("figures/size_shift_si.pdf"), width = 5, height = 3.5)
# print(
# df_td %>% 
#   filter(name == "D") %>% 
#   mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)])) %>% 
#   ggplot(aes(x=x))+
#   geom_line(data=. %>% filter(!is.na(y_amb)), aes(y=y_amb, col="AMB (CO2 = 414 ppm)"), linewidth=0.8)+
#   geom_line(data=. %>% filter(!is.na(y_ele)), aes(y=y_ele, col="ELE (CO2 = 614 ppm)"), linewidth=0.8)+
#   amz_theme()+
#   scale_y_log10(limits=c(1e-3, 1000))+
#   scale_x_continuous(n.breaks = 3, limits = c(0.01, 1.2))+
#   geom_label(data = . %>% count(name) %>% mutate(label = letters[row_number()]),
#              aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
#   scale_color_manual(values=c(
#     "AMB (CO2 = 414 ppm)" = col_amb,
#     "ELE (CO2 = 614 ppm)" = col_ele
#   ))+
#   labs(y="Density<br>(stems cm<sup>&minus;1</sup> ha<sup>&minus;1</sup>)", x=labels2["D"], color="")
# )
# dev.off()


##### Sample for PPTs #### 

q1 = l$dat_d %>% 
  mutate(
    GPP = GPP*365.2425,
    NPP = NPP*365.2425,
  ) %>% 
  filter(YEAR < -19880) %>% 
  mutate(YEAR = YEAR + 20000) %>% 
  ggplot(aes(x=YEAR))+
  geom_line(aes(y=GPP), linewidth=0.4)+
  geom_line(aes(y=NPP), linewidth=0.1)+
  annotate(geom="rect", xmin=-Inf, xmax=Inf, ymin=calib %>% filter(name=="GPP") %>% pull(min), ymax=calib %>% filter(name=="GPP") %>% pull(max), fill=alpha("seagreen2", 0.5))+
  geom_hline(yintercept = calib %>% filter(name=="NPP") %>% pull(mean), col=alpha("seagreen2", 0.5), linewidth=2)+
  amz_theme()+
  scale_x_continuous(n.breaks = 3)+
  # geom_label(data = tibble(label="a"),
  #            aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  labs(y="**Gross**/Net productivity<br>(kg m<sup>&minus;2</sup> yr<sup>&minus;1</sup>)", x="Years")+
  theme(plot.title = ggtext::element_markdown(lineheight=1.2, colour = "grey40", hjust=0.5))+
  ggtitle("A. CO<sub>2</sub> fluxes")

q2 = dist_amb1 %>% 
  filter(period == "AMB") %>% 
  ggplot() +
  geom_line(aes(x=size, y=density, group=period), linewidth=0.8)+
  scale_y_log10(limits=c(1e-3, 1000))+
  geom_point(data=dist_obs, aes(x=xobs, y=yobs), shape = 21, col="seagreen", fill=alpha("seagreen2", 0.5), size=3)+
  # geom_label(data = tibble(label="b"),
  #            aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  xlab("Diameter<br>(m)")+
  ylab("Density<br>(stems cm<sup>&minus;1</sup> ha<sup>&minus;1</sup>)")+
  scale_x_continuous(n.breaks = 3, limits = c(0.01, 1.2))+
  theme_bw()+
  amz_theme()+
  theme(plot.title = ggtext::element_markdown(lineheight=1.2, colour = "grey40", hjust=0.5))+
  ggtitle("B. Size distribution")

q3 = df_td %>% 
  filter(name == "WD") %>% 
  ggplot(aes(x=x))+
  geom_line(aes(y=y_obs, col="Observed"))+
  geom_line(aes(y=y_obs), col="seagreen")+
  geom_ribbon(aes(ymax=y_obs, ymin=0), fill=alpha("seagreen2", 0.5), alpha=0.2)+
  geom_line(aes(y=y_amb, col="Predicted"), linewidth=0.8)+
  # geom_label(data = tibble(label="c"),
  #            aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  amz_theme()+
  scale_x_continuous(n.breaks = 3)+
  labs(y="Basal-area density<br>(m<sup>2</sup> trait-unit<sup>&minus;1</sup> ha<sup>&minus;1</sup>)", x=labels2["WD"], col="")+
  theme(plot.title = ggtext::element_markdown(lineheight=1.2, colour = "grey40", hjust=0.5))+
  ggtitle("C. Trait strategy")+
  scale_color_manual(values = c("Observed"=alpha("seagreen2", 0.5), 
                                "Predicted" = "black"))


library(patchwork)

cairo_pdf(here::here("figures/sample_prediction.pdf"), height = 3.5, width=8)
q1+q2+q3 + plot_layout(guides="collect") & theme(legend.position="bottom")
dev.off()

