library(PlantFATE)
library(tidyverse)

print(getwd())

input_dir = here::here("input_data/")
output_dir = here::here("pfate_output_co2scan/")
expt_dir = "AmzMIP_HIST_1314.2_evol_20ky_4"

read_scan = function(co2, ymin=14000, ymax=14500){
  expt_dir = paste0()
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


# # Ambient Env
# zp_amb = c(19.22669, 15.71661, 3.15710, 0.00000)
# co_amb = c(1.0000000, 0.4712540, 0.2218492, 0.0711068)
# 
# # eCO2 env
# zp_ele = c(20.695389, 18.106550, 14.087510, 2.206985, 0.000000)
# co_ele = c(1.000000, 4.712543e-01, 2.220795e-01, 1.032055e-01, 2.763073e-02)
# Above vectors are from very old runs

# Numbers from CO2 = 368.9 run
zp_amb = c(19.237911, 15.482511, 2.995447, 0.000)
co_amb = c(1.000000e+00, 4.712536e-01, 2.175038e-01, 6.701001e-02)

# CO2 = 1014
zp_ele = c(21.51922, 17.97074, 8.56094, 0.00000)
co_ele = c(1.000000e+00, 4.712548e-01, 2.199178e-01, 7.424568e-02)

# CO2 = 1314
zp_ele = c(21.82550, 18.42982, 11.61380,  0.00000)
co_ele = c(1.00000000, 0.47125256, 0.22044767, 0.08341148)

params_file = here::here("config_files/p_amz_mip_final.ini")

fitness_wd = function(wd, zp, co, co2, pfile=params_file,
                      lma=0.15, hmat, p50){
  lho = new(LifeHistoryOptimizer, pfile)

  lho$set_i_metFile(here::here("input_data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE_new.csv"))
  lho$set_a_metFile(here::here("input_data/MetData_AmzFACE_Monthly_2000_2015_PlantFATE_new.csv"))
  lho$init_co2(co2)

  lho$traits0$lma = lma
  lho$traits0$hmat = hmat
  lho$traits0$p50_xylem = p50
  lho$traits0$wood_density = wd
  # lho$traits0$zeta = zeta

  lho$env$z_star = zp
  lho$env$canopy_openness = co

  lho$init()
  
  # lho$printPlant()
  lho$calcFitness()
}

wd = seq(350, 900, length.out=20)

traits = data.frame(i=1:20,
                    lma = runif(n = 20, 0.08, 0.25),
                    hmat = runif(n = 20, 6, 35),
                    p50 = runif(n = 20, -4, -0.5))

traits_wd = list(i=1:20, wd=wd) %>% 
  cross_df() %>% 
  left_join(traits) %>% 
  select(-i) 

## Trial with one species

dat_amb = wd %>% purrr::map_dbl(.f = fitness_wd, zp=zp_amb, co=co_amb, co2=368.9, hmat=21.4, p50 = -2.16, lma=0.15)
dat_ele_base = wd %>% purrr::map_dbl(fitness_wd, zp=zp_amb, co=co_amb, co2=1014, hmat=21.18, p50 = -2.29, lma=0.15)
dat_ele = wd %>% purrr::map_dbl(fitness_wd, zp=zp_ele, co=co_ele, co2=1014, hmat=21.18, p50 = -2.29, lma=0.15)
dat_ele_t = wd %>% purrr::map_dbl(fitness_wd, zp=zp_ele, co=co_ele, co2=1014, hmat=23.6, p50 = -1.7, lma=0.15)

wd_opt = function(wd, fitness){
  fitness = fitness/max(fitness)
  fitness = fitness^20
  sum(wd*fitness)/sum(fitness)
}

par(mfrow=c(1,1))
matplot(y=cbind(dat_amb/max(dat_amb), 
                dat_ele_base/max(dat_ele_base), 
                dat_ele/max(dat_ele),
                dat_ele_t/max(dat_ele_t)), 
        x = wd, type="l", lty=1, col=c("black", "grey", "yellow3", "brown"), ylab="Fitness", xlab="Wood density", cex.lab=1.3, lwd=2)
abline(v = c(wd_opt(wd, dat_amb), 
           wd_opt(wd, dat_ele_base),
           wd_opt(wd, dat_ele),
           wd_opt(wd, dat_ele_t)),
           col=c("black", "grey", "yellow3", "brown"))



#### Final plot with 20 species ####

fitness_aCaE = traits_wd %>% purrr::pmap(fitness_wd, zp=zp_amb, co=co_amb, co2=368.9) %>% unlist()
fitness_eCaE = traits_wd %>% purrr::pmap(fitness_wd, zp=zp_amb, co=co_amb, co2=1314.2) %>% unlist()
fitness_eCeE = traits_wd %>% purrr::pmap(fitness_wd, zp=zp_ele, co=co_ele, co2=1314.2) %>% unlist()

source(here::here("Rscripts/definitions.R"))

traits_wd  %>% 
  mutate(
    fitness_aCaE = fitness_aCaE,
    fitness_eCaE = fitness_eCaE,
    fitness_eCeE = fitness_eCeE
  ) %>% write.csv(here::here("summarized_outputs/fitness_landscape.csv"), row.names = F)

traits_fitness = read.csv(here::here("summarized_outputs/fitness_landscape.csv"))

cairo_pdf(here::here("figures/optimal_wd.pdf"), width = 3.7, height = 2.5)
print(
traits_fitness %>% 
  group_by(lma, hmat, p50) %>% 
  summarize(
    wd_opt_aCaE = wd_opt(wd, fitness_aCaE),
    wd_opt_eCaE = wd_opt(wd, fitness_eCaE),
    wd_opt_eCeE = wd_opt(wd, fitness_eCeE)
  ) %>% 
  ungroup() %>% 
  select(starts_with("wd_opt")) %>% 
  pivot_longer(everything()) %>% 
  ggplot()+
    geom_boxplot(aes(y=value, x=name, fill=name))+
  scale_fill_manual(values = c(
    wd_opt_aCaE = "lightgreen",
    wd_opt_eCaE = "greenyellow",
    wd_opt_eCeE = "darkgreen"
  )) + 
  scale_x_discrete(labels = c(
    wd_opt_aCaE = "aC+aI",
    wd_opt_eCaE = "eC+aI",
    wd_opt_eCeE = "eC+eI"
  ))+
  guides(fill="none")+
  ylab(exprlabel("Optimal", "wood density", "(kg m"^"-3"*")"))+
  xlab("Scenario")+
  theme_bw()
)
dev.off()



