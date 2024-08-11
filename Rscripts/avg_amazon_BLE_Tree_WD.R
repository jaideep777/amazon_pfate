library(tidyverse)
library(BIOMASS)

try_cat = readr::read_csv("/home/jjoshi/Downloads/Try20248111693TRY_Categorical_Traits_Lookup_Table_2012_03_17_TestRelease/TRY_Categorical_Traits_Lookup_Table_2012_03_17_TestRelease.csv")

data("wdData")

wdData %>%
  filter(regionId == "SouthAmericaTrop") %>% dim()

wd_withmeta = wdData %>%
  # filter(regionId == "SouthAmericaTrop") %>%
  left_join(try_cat,
            by=c("family"="Family",
                 "genus"="Genus",
                 "species"="SpeciesEpithet"))

wd_withmeta %>%
  # filter(regionId == "SouthAmericaTrop") %>%
  filter(LeafType == "broadleaved") %>%
  filter(PlantGrowthForm == "tree") %>%
  filter(LeafPhenology == "evergreen") %>%
  group_by(family, genus, species) %>%
  summarize(species_wd = mean(wd)) %>%
  pull(species_wd) %>%
  mean()
