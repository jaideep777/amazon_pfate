library(tidyverse)
# rm(list=ls())

input_dir = here::here("input_data/")
output_dir = here::here("pfate_output_co2scan/")

read_scan = function(co2, ymin=14000, ymax=14500){
  expt_dir = paste0("AmzMIP_HIST_",co2,"_evol_20ky_4")
  wd_back = getwd()
  setwd(paste0(output_dir,"/",expt_dir))
  
  l = list(
    Zp = read.csv("z_star.csv", header=F, col.names = paste0("V", 1:50)),
    traits = read.csv("traits.csv"),
    dat_d = readr::read_csv("D_PFATE.csv"),
    dat3 = read.csv("Y_mean_PFATE.csv")
  )
  
  names(l$Zp)[1] = c("YEAR")

  l$dat = l$dat_d |> 
    dplyr::mutate(YEAR = as.integer(YEAR)) |> 
    dplyr::group_by(YEAR) |> 
    dplyr::summarize_all(mean)
  
  setwd(wd_back)  
  
  l$Zp %>% 
    mutate(YEAR=as.integer(YEAR)) %>% 
    filter(YEAR > ymin & YEAR < ymax) %>% select(YEAR:V10) %>% 
    left_join(l$dat %>%  mutate(YEAR=as.integer(YEAR)) %>% filter(YEAR > ymin & YEAR < ymax)) %>% 
    left_join(l$traits %>% mutate(YEAR=as.integer(YEAR)) %>% filter(SPP == "Spp_1") %>% filter(YEAR > ymin & YEAR < ymax) %>% select(-SPP)) %>% 
    left_join(l$dat3 %>% mutate(YEAR=as.integer(YEAR)) %>% filter(YEAR > ymin & YEAR < ymax)) %>% 
    as_tibble() %>% 
    colMeans(na.rm=T)
    
}
  

## Final runs (co2 scan)
co2 = seq(414.2, by=50, length.out=20)

read_scan(co2[1])

df = co2 %>% map_df(~read_scan(., ymin=19900, ymax=20000))

df2 = df %>% bind_rows(
  read_scan(414.2, ymin=1920, ymax=2020)
)

df2 %>% write.csv(here::here("summarized_outputs/co2_scan_100y_summary.csv"), row.names = F)

df2 = read.csv(here::here("summarized_outputs/co2_scan_100y_summary.csv"))

ssp_co2 = c(
  ssp1_1.9 = 393.5,
  ssp1_2.6 = 445.6,
  ssp2_4.5 = 602.8,
  ssp3_7.0 = 867.2,
  ssp5_8.5 = 1135.2
)

ssp_cols = rgb(
  red   = c(0, 23, 247, 231, 149)/255,
  green = c(173, 60, 148, 29, 27)/255,
  blue  = c(207, 102, 32, 37, 30)/255
)
names(ssp_cols) = names(ssp_co2)

library(ggnewscale)

source(here::here("Rscripts/definitions.R"))

cairo_pdf(file=here::here("figures/co2_scan.pdf"), width = 7, height = 6)
print(
df2 %>% 
  mutate(AGB=CL+CW) %>% 
  mutate(BA=BA*1e4,
         GPP=GPP*365.2425,
         NPP=NPP*365.2425,
         RAU=RAU*365.2425) %>% 
  select(CO2, GPP, NPP, LAI, WD, AGB, BA, HMAT, P50X, RAU) %>% 
  pivot_longer(-CO2) %>% 
  # Create facet labels
  mutate(name = factor(name, levels = unique(name), labels = labels[unique(name)])) %>% 
  # plot
  ggplot(aes(x=CO2, y=value)) +
  geom_point(col="black", size=0.5)+
  geom_smooth(se=F, col=scales::alpha("grey50", 0.5), linewidth=0.75)+
  geom_line(data= df2 %>% select(CO2, V2:V6) %>%
              setNames(c("CO2", paste0("z", 1:5))) %>% 
              pivot_longer(-CO2, names_to="level") %>%
              mutate(name=factor("Z", levels="Z", labels=labels["Z"])),
              aes(colour=level, group=level))+
  facet_wrap(~name, scales="free_y", 
             strip.position = "left",
             ncol=2, axes = "margins",
             labeller = label_parsed)+
  theme_bw()+
  scale_colour_viridis_d(direction = -1, end=0.95)+
  labs(colour="Canopy\nlayer", x=expression("CO"[2]~"(ppm)"))+
  new_scale_color() +
  labs(y="", colour="CO2 in year\n2100 in")+
  theme(strip.placement = "outside",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  geom_vline(data = data.frame(ssp = names(ssp_co2), co2=ssp_co2), 
             aes(xintercept=co2, col=ssp), linewidth=0.2)+
  scale_color_manual(values = ssp_cols)
)  
dev.off()


## Sapwood respiration rate visualization
ggplot(data=data.frame(x=c(exp(-3.75), exp(-1.25), exp(1.25))*1e-9*12*86400*365), aes(y=x,x=1))+geom_violin()+geom_point(x=1,y=0.4)

