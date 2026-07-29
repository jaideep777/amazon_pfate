library(rphydro)
library(tidyverse)
library(ggplot2)

## Set up Phydrp parameters and environmental variables identical to those used in PlantFATE simulations

kphio = 0.0593;        # quantum yield efficiency
ppfd = 410.4675365;          # umol/m2/s
ppfd_max = 1798.319161;
vpd  = 652.7145833;         # Pa
co2  = 368;          # ppm
elv  = 0;            # m.a.s.l.
fapar = 1-exp(-0.5*5.1);         # fraction
rdark = 0.015;
tc = 25.75380729;
vwind = 3;
netrad = ppfd/2
psi_soil = -0.05
pa = calc_patm(elv)

par_cost = list(alpha=0.1, gamma=0.5);
par_plant = list(conductivity=5e-17, psi50=-2, b=1);
options = list(gs_method = "GS_IGF", 
               et_method = "ET_DIFFUSION",
               ftemp_vj_method = "FV_kumarathunge19",
               ftemp_rd_method = "FR_heskel16",
               ftemp_br_method = "FB_atkin15",
               scale_alpha = F)

aco2 = 368.9
eco2 = 591.6

## Calculate Acclimated Vcmax/Jmax at baseline and elevated CO2
dat_acc_amz = tibble(co2 = c(aco2, eco2)) |> 
    mutate(dat = purrr::map(
        .x=co2, 
        .f = ~rphydro_analytical(tc, tc, ppfd_max, netrad, vpd, .x, pa, fapar, kphio, psi_soil, rdark, vwind, par_plant, par_cost, options)
        )) |> 
    unnest_wider(dat) |>
    mutate(type = "acc")
  
## Calculated instantaneous quantities with daytime (12 hr) mean light. This is what results in daily GPP
dat_inst_amz = tibble(co2 =c(aco2, eco2)) |> 
    left_join(dat_acc_amz |> select(co2, vcmax25, jmax25)) |>
    mutate(dat = purrr::pmap(
        .l=list(co2, vcmax25, jmax25), 
        .f = ~rphydro_instantaneous_analytical(..2, ..3, tc, tc, ppfd_max/2, netrad, vpd, ..1, pa, fapar, kphio, psi_soil, rdark, vwind, par_plant, par_cost, options)
        )) |>
    select(-vcmax25, -jmax25) |>
    unnest_wider(dat) |> 
    mutate(a=a/2) |> # divide by two to account for nighttime in daily average
    mutate(type = "inst_daily")

## Calculated instantaneous photosynthesis with saturated light. This will give beta for Asat
dat_inst_amz_ls = tibble(co2 =c(aco2, eco2)) |> 
    left_join(dat_acc_amz |> select(co2, vcmax25, jmax25)) |>
    mutate(dat = purrr::pmap(
        .l=list(co2, vcmax25, jmax25), 
        .f = ~rphydro_instantaneous_analytical(..2, ..3, tc, tc, 1500, netrad, vpd, ..1, pa, fapar, kphio, psi_soil, rdark, vwind, par_plant, par_cost, options)
        )) |>
    select(-vcmax25, -jmax25) |>
    unnest_wider(dat) |> 
    mutate(type = "inst_ls")

df <- dplyr::bind_rows(
    dat_acc_amz,
    dat_inst_amz,
    dat_inst_amz_ls
)

beta_a_df <- df |> 
    select(co2, a, type) |> 
    pivot_longer(a) |>
    pivot_wider(names_from = co2, values_from=value) |> 
    mutate(beta = log(`591.6`/`368.9`)/log(591.6/368.9)) |>
    mutate(pc_change = (`591.6`/`368.9`-1)*100)


