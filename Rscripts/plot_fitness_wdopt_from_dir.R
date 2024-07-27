library(PlantFATE)
library(tidyverse)

print(getwd())

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


# input_dir = here::here("input_data/")
# output_dir = here::here("pfate_output_supplementary/")
# 
# expt_dir = "AmzMIP_ELE_614.2_evol_20ky_4_plusRs"
# summarize_outputs(input_dir, expt_dir, output_dir, ymin=6000, ymax=7000)
# 
# expt_dir = "AmzMIP_ELE_614.2_evol_20ky_4_plusZeta"
# summarize_outputs(input_dir, expt_dir, output_dir, ymin=6000, ymax=7000)


# # Evol vs comm Comparison run
# input_dir = here::here("input_data/")
# output_dir = here::here("pfate_output_evol_vs_comm/")
# expt_dir = "AMB_ELE_wd_evol_20ky"
# 
# dat_amb = summarize_outputs(input_dir, expt_dir, output_dir, ymin=1900, ymax=2000)
# dat_ele = summarize_outputs(input_dir, expt_dir, output_dir, ymin=19900, ymax=20000)
# 
# dat_amb$co2 = 414
# dat_ele$co2 = 828

# Actual run
input_dir = here::here("input_data/")
output_dir = here::here("pfate_output_co2scan/")
expt_dir = "AmzMIP_HIST_1314.2_evol_20ky_4"

dat_amb = summarize_outputs(input_dir, expt_dir, output_dir, ymin=1900, ymax=2000)
dat_ele = summarize_outputs(input_dir, expt_dir, output_dir, ymin=19900, ymax=20000)

dat_amb$co2 = 368.9
dat_ele$co2 = 1314.2

fitness_wd = function(wd, zp, co, co2, pfile,
                      lma, hmat, p50){
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

#### 1 species plot from real simulation ####

fitness_amb = wd %>% purrr::map_dbl(.f = fitness_wd, 
                                pfile = dat_amb$pfile,
                                zp=dat_amb$env$z, co=dat_amb$env$co, 
                                co2=dat_amb$co2, 
                                hmat=dat_amb$traits$HMAT, p50 = dat_amb$traits$P50X, lma=dat_amb$traits$LMA)

fitness_eC = wd %>% purrr::map_dbl(.f = fitness_wd, 
                                pfile = dat_amb$pfile,
                                zp=dat_amb$env$z, co=dat_amb$env$co, 
                                co2=dat_ele$co2, 
                                hmat=dat_amb$traits$HMAT, p50 = dat_amb$traits$P50X, lma=dat_amb$traits$LMA)

fitness_eCeE = wd %>% purrr::map_dbl(.f = fitness_wd, 
                                   pfile = dat_amb$pfile,
                                   zp=dat_ele$env$z, co=dat_ele$env$co, 
                                   co2=dat_ele$co2, 
                                   hmat=dat_amb$traits$HMAT, p50 = dat_amb$traits$P50X, lma=dat_amb$traits$LMA)

fitness_eCeEeT = wd %>% purrr::map_dbl(.f = fitness_wd, 
                                     pfile = dat_amb$pfile,
                                     zp=dat_ele$env$z, co=dat_ele$env$co, 
                                     co2=dat_ele$co2, 
                                     hmat=dat_ele$traits$HMAT, p50 = dat_ele$traits$P50X, lma=dat_ele$traits$LMA)


df_real = cbind(wd,
                fitness_amb, 
                fitness_eC, 
                fitness_eCeE,
                fitness_eCeEeT) %>% as.data.frame()

df_real %>% write_csv(here::here("summarized_outputs/fitness_landscape_hist_ele_co2_1314.2.csv"))

wd_opt = function(wd, fitness){
  fitness = fitness/max(fitness)
  fitness = fitness^20
  sum(wd*fitness)/sum(fitness)
}

p_real = df_real %>% 
  pivot_longer(-wd) %>% 
  # pull(name) %>% unique()
  mutate(name = factor(name, levels = unique(name), labels = c("aC+aE+aT", "eC+aE+aT", "eC+eE+aT", "eC+eE+eT") )) %>% 
  group_by(name) %>% 
  mutate(value = value/max(value),
         mean = wd_opt(wd, value)) %>% 
  ggplot(aes(x=wd, y=value))+
  geom_line(aes(group=name, color=name))+
  labs(x=labels["WD"]$WD, y="Fitness", col="Scenario")+
  scale_color_manual(values = c(`aC+aE+aT`="lightgreen",
                                `eC+aE+aT`="greenyellow",
                                `eC+eE+aT`="darkgreen",
                                `eC+eE+eT`="darkcyan"))+
  geom_vline(aes(xintercept = mean, col=name), linewidth=0.3)+
  theme_bw()+
  amz_theme()

#### Plot with 20 species ####

traits = data.frame(i=1:20,
                    lma = runif(n = 20, 0.08, 0.25),
                    hmat = runif(n = 20, 6, 35),
                    p50 = runif(n = 20, -4, -0.5))

traits_wd = list(i=1:20, wd=wd) %>% 
  cross_df() %>% 
  left_join(traits) %>% 
  select(-i) %>% 
  mutate(pfile = here::here("config_files/p_amz_mip_final.ini"))


fitness_aCaE = traits_wd %>% purrr::pmap(fitness_wd, zp=zp_amb$env$z, co=zp_amb$env$co, co2=368.9) %>% unlist()
fitness_eCaE = traits_wd %>% purrr::pmap(fitness_wd, zp=zp_amb$env$z, co=zp_amb$env$co, co2=1314.2) %>% unlist()
fitness_eCeE = traits_wd %>% purrr::pmap(fitness_wd, zp=zp_ele$env$z, co=zp_ele$env$co, co2=1314.2) %>% unlist()

traits_wd  %>% 
  mutate(
    fitness_aCaE = fitness_aCaE,
    fitness_eCaE = fitness_eCaE,
    fitness_eCeE = fitness_eCeE
  ) %>% write.csv(here::here("summarized_outputs/fitness_landscape_20spp.csv"), row.names = F)

traits_fitness = read.csv(here::here("summarized_outputs/fitness_landscape_20spp.csv"))

p_20spp = traits_fitness %>% 
  group_by(lma, hmat, p50) %>% 
  summarize(
    wd_opt_aCaE = wd_opt(wd, fitness_aCaE),
    wd_opt_eCaE = wd_opt(wd, fitness_eCaE),
    wd_opt_eCeE = wd_opt(wd, fitness_eCeE)
  ) %>% 
  ungroup() %>% 
  select(starts_with("wd_opt")) %>% 
  pivot_longer(everything()) %>% 
  mutate(name = factor(name, levels = unique(name), labels = c("aC+aE+aT", "eC+aE+aT", "eC+eE+aT") )) %>% 
  ggplot()+
  geom_boxplot(aes(x=value, y=name, fill=name))+
  scale_fill_manual(values = c(`aC+aE+aT`="lightgreen",
                               `eC+aE+aT`="greenyellow",
                               `eC+eE+aT`="darkgreen",
                               `eC+eE+eT`="darkcyan"))+
  scale_y_discrete(limits=rev)+
  xlab(exprlabel("Optimal", "wood density", "(kg m"^"-3"*")"))+
  labs(y="Scenario", fill="Scenario")+
  theme_bw()+
  amz_theme()

cairo_pdf(here::here("figures/optimal_wd.pdf"), width = 4.6, height = 4.6)
print(
p_real/p_20spp + plot_layout(guides="collect")
)
dev.off()



# par(mfrow=c(1,1))
# matplot(y=cbind(fitness_amb/max(fitness_amb), 
#                 fitness_eC/max(fitness_eC), 
#                 fitness_eCeE/max(fitness_eCeE),
#                 fitness_eCeEeT/max(fitness_eCeEeT)), 
#         x = wd, type="l", lty=1, col=c("black", "grey", "yellow3", "brown"), ylab="Fitness", xlab="Wood density", cex.lab=1.3, lwd=2)
# abline(v = c(wd_opt(wd, fitness_amb), 
#              wd_opt(wd, fitness_eC),
#              wd_opt(wd, fitness_eCeE),
#              wd_opt(wd, fitness_eCeEeT)),
#            col=c("black", "grey", "yellow3", "brown"))
# abline(v = c(dat_amb$traits$WD,
#              dat_ele$traits$WD),
#        col=c("black", "brown"), lty=2)
       




