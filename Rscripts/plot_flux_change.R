rm(list=ls())
library(tidyverse)
library(patchwork)

source("~/codes/Plant-FATE/R/process_outputs.R")
source(here::here("Rscripts/definitions.R"))

input_dir = here::here("input_data/")

output_dir = data_path
expt_dir = "AmzMortality_ELE_rs0.5e-3m_inf0.002_K_leaf0.5e-16_zeta0.2nspecies2timestep0.0416"

co2_amb = read.csv(paste0(input_dir, "/", "CO2_AMB_AmzFACE2000_2100.csv"))
co2_ele = read.csv(paste0(input_dir, "/", "CO2_ELE_AmzFACE2000_2100.csv"))

co2_amb |> rename(CO2_AMB = CO2) |> 
  left_join(co2_ele |> rename(CO2_ELE = CO2)) |>
  pivot_longer(-Year) |>
  ggplot(aes(x=Year, y=name)) +
  geom_line(aes(y=value, colour=name, group=name)) +
  amz_theme()

co2_change_year = 2000

l = pf_read_outputs(input_dir, output_dir, expt_dir)

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

dat_traits = l$traits |> 
    rename(PID = SPP) |>
    dplyr::mutate(YEAR = as.integer(YEAR)) |>
    filter(!grepl("probe", PID)) |>
    select(YEAR, PID, WD, P50X, HMAT) |> 
    left_join(
      l$dat2 |> 
      dplyr::mutate(YEAR = as.integer(YEAR)) |> 
      filter(!grepl("probe", PID)) |>
      select(YEAR, PID, BA)
    ) |> 
    group_by(YEAR) |> 
    summarize(
      WD = weighted.mean(WD, BA, na.rm=TRUE),
      P50X = weighted.mean(P50X, BA, na.rm=TRUE),
      HMAT = weighted.mean(HMAT, BA, na.rm=TRUE)
    )

p_flux = dat_flux %>% 
    ggplot(aes(x=YEAR, y=value))+
    geom_line(alpha=0.3)+
    geom_smooth(method = "loess", span=0.15, se=F, col="black", linewidth=0.5)+
    geom_line(data = dat_z,
              aes(colour=level, group=level))+
    geom_rect(data = . %>% group_by(name) %>% slice(1), aes(xmin = co2_change_year, xmax=2100, ymin=-Inf, ymax=Inf), fill=col_ele, alpha=0.1)+
    geom_rect(data = dat_z %>% group_by(name) %>% slice(1), aes(xmin = co2_change_year, xmax=2100, ymin=-Inf, ymax=Inf), fill=col_ele, alpha=0.1)+
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
  pivot_longer(-YEAR) %>%
  mutate(YEAR=as.integer(YEAR)) %>% 
  mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)])) |>
  filter(YEAR %in% seq(-20000,20000, by=10)) %>% 
  ggplot(aes(x=YEAR, y=value))+
  geom_line(alpha=0.3)+
  geom_smooth(method = "loess", span=0.15, se=F, col="black", linewidth=0.5)+
  geom_rect(data = . %>% group_by(name) %>% slice(1), aes(xmin = co2_change_year, xmax=Inf, ymin=-Inf, ymax=Inf), fill=col_ele, alpha=0.1)+
  geom_label(data = . %>% count(name) %>% mutate(label = letters[row_number()+9]),
             aes(x=-Inf, y=Inf, label=label), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  facet_wrap(~name, scales="free_y", 
             strip.position = "left",
             ncol=3)+
  scale_x_continuous(n.breaks=3)+
  labs(y="", x="Year")+
  amz_theme()

cairo_pdf(here::here("figures/flux_change_ELE.pdf"), width=7, height=6.4)            
print(
  p_flux/p_traits + plot_layout(heights=c(3,1))
)
dev.off()


beta_dat = l$dat_d %>% 
  select(YEAR, GPP, NPP, GS, VCMAX, RAU, MORT) %>% 
  mutate(GPP = GPP*365.2425,
         NPP = NPP*365.2425,
         MORT = MORT*365.2425) %>% 
  mutate(IWUE = NPP/GS) %>% 
  pivot_longer(-YEAR) %>% 
  bind_rows(
    l$dat3 %>% 
      mutate(AGB = CL+CW,
             BGB = CCR+CFR,
             BA = BA*1e4) %>% 
      select(YEAR, BA, AGB, BGB, LAI) %>% 
      pivot_longer(-YEAR)
  ) %>% 
  bind_rows(
    l$Zp %>% 
      select(YEAR, z3=V4) %>% 
      pivot_longer(-YEAR)
  ) %>% 
  bind_rows(
    l$traits %>% 
      select(YEAR, WD, HMAT, P50X) %>% 
      pivot_longer(-YEAR)
  ) %>% 
  mutate(hist = cut(YEAR, 
                    breaks = c(-Inf, 1970, 2000,   2070,    2100,     19900,     Inf), 
                    labels = c("prehist",  "hist", "mid", "ele_st", "midele", "ele_lt")
                    )
         ) %>%
  group_by(hist, name) %>% 
  summarize(value = mean(value)) %>% 
  ungroup() %>% 
  bind_rows(
    data.frame(hist=c("hist", "ele_st"),
               name = c("A","A"),
               value = c(9.193998, 13.432466))
  ) %>% 
  pivot_wider(names_from = hist) %>% 
  mutate(pc_change_st = (ele_st-hist)/hist*100) %>% 
        #  pc_change_lt = (ele_lt-hist)/hist*100) %>% 
  mutate(beta_st = log(ele_st/hist)/log(614.4/368.9))
        #  beta_lt = log(ele_lt/hist)/log(614.4/368.9)) 

beta_dat %>% 
  write_csv(here::here("summarized_outputs/flux_change.csv"))

# Betas reported in https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.16866
# VCMAX, -0.38, -0.48, -0.28
beta_obs =
  read.csv(text = gsub(pattern = " ", replacement="",
     "name, mean, min, max, type, source
      A,  0.75, 0.5,  1.0, fluxes, stocker 
      GPP,  0.49, 0.17,  0.82, fluxes, stocker 
      VCMAX,  -0.13, -0.27,  0.0, fluxes, stocker 
      AGB, 0.40, 0.34, 0.47, structure, stocker
      IWUE, 1.1, 0.65, 1.1, fluxes, walker"
     ),
     header=T, sep=",", as.is = F) %>%
  as_tibble()


p1 = beta_dat %>% 
  ungroup() %>% 
  select(name, beta_st) %>% 
  filter(name %in% c("GPP", "VCMAX", "NPP", "IWUE", "AGB", "MORT", "A")) %>% 
  arrange(match(name, c("AGB", "MORT", "IWUE", "NPP", "VCMAX", "GPP", "A"))) %>% 
  left_join(beta_obs) %>% 
  # Create facet labels
  mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)])) %>% 
  ggplot()+
  geom_col(aes(x=beta_st, y=name, fill="Predicted"), alpha=0.5)+
  geom_errorbar(aes(y=name, x=mean, xmin=min, xmax=max, col="Observed"), width = 0.2, linewidth=0.8)+
  amz_theme()+
  theme(axis.text = ggtext::element_markdown(lineheight=1.2))+
  labs(y="", x="Response ratio")+
  scale_x_continuous(limits = c(-0.6,1.2), breaks=c(-0.5, 0, 0.5, 1))+
  geom_label(data = . %>% slice(1),
             aes(x=-Inf, y=Inf, label="a"), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  scale_fill_manual(values = c("Predicted"=col_amb))+
  scale_color_manual(values = c("Observed"=col_obs))+
  labs(color="", fill="")

p1

cairo_pdf(here::here("figures/flux_change_beta.pdf"), width=5, height=4)
p1 + plot_layout(guides="collect")&theme(legend.position = "top")
dev.off()



