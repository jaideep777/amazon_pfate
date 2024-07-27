##  generate_co2_scenarios

dat = read.csv(here::here("input_data/CO2_AMB_AmzFACE2000_2100.csv"))

for (eco2 in seq(414.2, by=50, length.out=20)){
  print(eco2)
  dat %>% 
    mutate(CO2=ifelse(Year>2020, yes=eco2, no=CO2)) %>% 
    write.csv(file=here::here(paste0("input_data/CO2_",eco2,"_AmzFACE2000_2100.csv")), row.names = F)
}
