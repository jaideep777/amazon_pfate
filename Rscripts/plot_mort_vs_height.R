library(tidyverse)

input_dir = here::here("input_data/")
output_dir = here::here("pfate_output_co2scan/")

## Final runs (MIP)
expt_dir = "AmzMIP_HIST_1314.2_evol_20ky_4"

source("~/codes/Plant-FATE/R/process_outputs.R")
source(here::here("Rscripts/definitions.R"))

summarize_outputs = function(input_dir, expt_dir, output_dir, ymin=14000, ymax=14500){
  wd_back = getwd()
  setwd(paste0(output_dir,"/",expt_dir))
  
  l = list(
    Zp = read.csv("z_star.csv", header=F, col.names = paste0("V", 1:50)),
    co = read.csv("canopy_openness.csv", header=F, col.names = paste0("V", 1:50)),
    traits = read.csv("traits.csv"),
    pfile = paste0(getwd(), "/p.ini")
  )
  
  names(l$Zp)[1] = c("YEAR")
  names(l$co)[1] = c("YEAR")
  
  setwd(wd_back)  
  
  zp = l$Zp %>% 
    mutate(YEAR=as.integer(YEAR)) %>% 
    filter(YEAR > ymin & YEAR < ymax) %>% 
    pivot_longer(-YEAR) %>% 
    mutate(name = sub('.', '', name) %>% as.numeric()) %>% 
    group_by(name) %>% 
    summarize(z = mean(value, na.rm=T)) %>% 
    drop_na()
  
  co = l$co %>% 
    mutate(YEAR=as.integer(YEAR)) %>% 
    filter(YEAR > ymin & YEAR < ymax) %>% 
    pivot_longer(-YEAR) %>% 
    mutate(name = sub('.', '', name) %>% as.numeric()) %>% 
    group_by(name) %>% 
    summarize(co = mean(value, na.rm=T)) %>% 
    drop_na()
  
  zp %>% left_join(co)
  
  traits = l$traits %>% 
    mutate(YEAR=as.integer(YEAR)) %>% 
    filter(!grepl("probe", SPP)) %>% 
    filter(YEAR > ymin & YEAR < ymax) %>% 
    pivot_longer(-c(YEAR, SPP)) %>% 
    group_by(name, SPP) %>% 
    summarize(value = mean(value, na.rm=T))
  
  list(
    env = zp %>% left_join(co),
    traits = traits %>% pivot_wider(names_from = name, values_from=value),
    pfile = l$pfile
  )
}

zp_amb = summarize_outputs(input_dir, expt_dir, output_dir, ymin = 1995, ymax=2000)
zp_ele = summarize_outputs(input_dir, expt_dir, output_dir, ymin = 19995, ymax=20000)

df_cohorts = readr::read_csv(paste0(output_dir,"/",expt_dir,"/cohort_props.csv"))

cairo_pdf(here::here("figures/mortality_with_height.pdf"), height = 2.5, width = 4.5)
print(
  df_cohorts %>%
    mutate(YEAR = as.integer(YEAR)) %>%
    filter(YEAR %in% c(2000, 20000)) %>%
    mutate(CO2 = factor(ifelse(YEAR == 2000, yes=368.9, no=1314.2))) %>% 
    ggplot(aes(x=height, y=mort*365.2425, colour = CO2, group=CO2))+
    geom_line()+
    geom_vline(xintercept = zp_amb$env$z, col=col_amb, alpha=0.2)+
    geom_vline(xintercept = zp_ele$env$z, col=col_ele, alpha=0.2)+
    scale_color_manual(values=c(`368.9`=col_amb, `1314.2`=col_ele))+
    labs(y="Mortality rate<br>(Yr<sup>&minus;1</sup>)", x="Height (m)", col="CO2")+
    amz_theme()    
)
dev.off()
