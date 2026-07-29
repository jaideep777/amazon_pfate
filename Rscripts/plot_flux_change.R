rm(list=ls())
library(tidyverse)
library(patchwork)

source("~/codes/Plant-FATE/R/process_outputs.R")
source(here::here("Rscripts/definitions.R"))
source(here::here("Rscripts/calc_asat_beta.R"))

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
    labs(y="", x="", colour="Canopy/nlevel")+
    amz_theme()


p_traits = dat_traits %>% 
  pivot_longer(-YEAR) %>%
  mutate(YEAR=as.integer(YEAR)) %>% 
  mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)])) |>
  filter(YEAR %in% seq(1900,2100, by=10)) %>% 
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

cairo_pdf(here::here("figures/flux_change_ELE.pdf"), width=7, height=6.6)            
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
                    breaks = c(-Inf,  1970, 2000,   2020,    2030,   2070,      2100,     19900,     Inf), 
                    labels = c("prehist",  "hist", "y2000s", "y2020s" , "mid", "ele_eoc", "midele", "ele_lt")
                    )
         ) %>%
  group_by(hist, name) %>% 
  summarize(value = mean(value)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = hist) %>% 
  mutate(pc_change_eoc = (ele_eoc-hist)/hist*100) %>% 
  mutate(pc_change_2020s = (y2020s-hist)/hist*100) %>% 
  mutate(pc_change_2000s = (y2000s-hist)/hist*100) %>% 
        #  pc_change_lt = (ele_lt-hist)/hist*100) %>% 
  mutate(beta_eoc = log(ele_eoc/hist)/log(614.4/368.9)) %>% 
  mutate(beta_2020s = log(y2020s/hist)/log(614.4/368.9)) |>
  mutate(beta_2000s = log(y2000s/hist)/log((614.4+568.9)/2/368.9)) |>
  dplyr::bind_rows(
    data.frame(
      name = "A", 
      beta_eoc = beta_a_df |> filter(type=="inst_ls") |> pull(beta),
      beta_2000s = beta_a_df |> filter(type=="inst_ls") |> pull(beta),
      beta_2020s = beta_a_df |> filter(type=="inst_ls") |> pull(beta)
    )
  )
        #  beta_lt = log(ele_lt/hist)/log(614.4/368.9)) 

beta_dat %>% 
  write_csv(here::here("summarized_outputs/flux_change.csv"))

# Digitized from Stocker et al: https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.20178
# iWUE from walker et al: https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.16866
beta_obs = read.csv(here::here("MESI_betas.csv")) |> 
      pivot_wider(names_from=metric, values_from=beta)

p1 = beta_dat %>% 
  ungroup() %>% 
  select(name, beta_eoc) %>% 
  pivot_wider(names_from=name, values_from=beta_eoc) |> 
  mutate(ANPP=NPP, BNPP=NPP, obs=1) |> 
  pivot_longer(-obs, values_to = "beta_eoc") |> 
  select(-obs) |> 
  filter(name %in% c("GPP", "VCMAX", "ANPP", "BNPP", "RAU", "LAI", "IWUE", "AGB", "MORT", "A")) %>% 
  arrange(match(name, c("AGB", "LAI", "MORT", "ANPP", "BNPP", "IWUE", "RAU", "VCMAX", "GPP", "A"))) %>% 
  left_join(beta_obs) %>% 
  # Create facet labels
  mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)])) %>% 
  ggplot()+
  geom_col(aes(x=beta_eoc, y=name, fill="Predicted"), alpha=0.5)+
  geom_errorbar(aes(y=name, x=mean, xmin=min, xmax=max, col="Observed"), width = 0.2, linewidth=0.8)+
  amz_theme()+
  theme(axis.text = ggtext::element_markdown(lineheight=1.2))+
  labs(y="", x="Response ratio")+
  scale_x_continuous(limits = c(-0.6,1.6), breaks=c(-0.5, 0, 0.5, 1))+
  geom_label(data = . %>% slice(1),
             aes(x=-Inf, y=Inf, label="a"), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  scale_fill_manual(values = c("Predicted"=col_amb))+
  scale_color_manual(values = c("Observed"=col_obs))+
  labs(color="", fill="")

p1

cairo_pdf(here::here("figures/flux_change_beta_eoc.pdf"), width=5, height=5)
p1 + plot_layout(guides="collect")&theme(legend.position = "top")
dev.off()


p2 = beta_dat %>% 
  ungroup() %>% 
  select(name, beta_2000s) %>% 
  pivot_wider(names_from=name, values_from=beta_2000s) |> 
  mutate(ANPP=NPP, BNPP=NPP, obs=1) |> 
  pivot_longer(-obs, values_to = "beta_2000s") |> 
  select(-obs) |> 
  filter(name %in% c("GPP", "VCMAX", "ANPP", "BNPP", "RAU", "LAI", "IWUE", "AGB", "MORT", "A")) %>% 
  arrange(match(name, c("AGB", "MORT", "ANPP", "BNPP", "RAU", "IWUE", "VCMAX", "LAI", "GPP", "A"))) %>% 
  left_join(beta_obs) %>% 
  # Create facet labels
  mutate(name = factor(name, levels = unique(name), labels = labels_nounit_oneline[unique(name)])) %>% 
  ggplot()+
  # geom_col(aes(x=beta_2000s, y=name, fill="Predicted"), alpha=0.5)+
  geom_errorbar(aes(y=name, x=mean, xmin=min, xmax=max, col="Observed"), width = 0.2, linewidth=0.8, alpha=0.7)+
  geom_point(aes(x=beta_2000s, y=name, col="Predicted"), alpha=1, size=3)+
  geom_vline(xintercept = 0, col="grey60")+
  amz_theme()+
  theme(axis.text = ggtext::element_markdown(lineheight=1.2))+
  labs(y="", x="Log-response ratio")+
  scale_x_continuous(limits = c(-0.6,1.6), breaks=c(-0.5, 0, 0.5, 1))+
  # geom_label(data = . %>% slice(1),
  #            aes(x=-Inf, y=Inf, label="a"), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  scale_fill_manual(values = c("Predicted"="orange2"))+
  scale_color_manual(values = c("Observed"=col_obs, "Predicted"="orange2"))+
  labs(color="", fill="")

p2

cairo_pdf(here::here("figures/flux_change_beta_2000s_orange_point.pdf"), width=6, height=4)
print(
p2 + plot_layout(guides="collect")&theme(legend.position = "top")
)
dev.off()


# df_coh <- readr::read_csv("c:/Users/Jaideep/OneDrive - IIASA/RESIST - Documents/Plant-FATE output_newmort/AmzMortality_AMB_rs0.5e-3m_inf0.01_K_leaf0.5e-16_zeta0.2nspecies2/cohort_props.csv")


p3 = beta_dat %>% 
  ungroup() %>% 
  select(name, beta_2020s) %>% 
  pivot_wider(names_from=name, values_from=beta_2020s) |> 
  mutate(ANPP=NPP, BNPP=NPP, obs=1) |> 
  pivot_longer(-obs, values_to = "beta_2020s") |> 
  select(-obs) |> 
  filter(name %in% c("GPP", "VCMAX", "ANPP", "BNPP", "RAU", "LAI", "IWUE", "AGB", "MORT", "A")) %>% 
  arrange(match(name, c("AGB", "LAI", "MORT", "ANPP", "BNPP", "IWUE", "RAU", "VCMAX", "GPP", "A"))) %>% 
  left_join(beta_obs) %>% 
  # Create facet labels
  mutate(name = factor(name, levels = unique(name), labels = labels2[unique(name)])) %>% 
  ggplot()+
  geom_col(aes(x=beta_2020s, y=name, fill="Predicted"), alpha=0.5)+
  geom_errorbar(aes(y=name, x=mean, xmin=min, xmax=max, col="Observed"), width = 0.2, linewidth=0.8)+
  amz_theme()+
  theme(axis.text = ggtext::element_markdown(lineheight=1.2))+
  labs(y="", x="Response ratio")+
  scale_x_continuous(limits = c(-0.6,1.6), breaks=c(-0.5, 0, 0.5, 1))+
  geom_label(data = . %>% slice(1),
             aes(x=-Inf, y=Inf, label="a"), inherit.aes = F, hjust=0, vjust=1, label.size = 0, size = 4.5) +
  scale_fill_manual(values = c("Predicted"=col_amb))+
  scale_color_manual(values = c("Observed"=col_obs))+
  labs(color="", fill="")

p3

cairo_pdf(here::here("figures/flux_change_beta_2020s.pdf"), width=5, height=5)
p3 + plot_layout(guides="collect")&theme(legend.position = "top")
dev.off()
