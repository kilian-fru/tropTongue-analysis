library(tidyverse)

rm(list=ls())

traits_melaug <- read.csv2("01_raw_data/traits_mel_aug.csv") %>%
  filter(Group=="Apidae: Apinae: Meliponini" | 
           Group=="Halictidae: Halictinae: Augochlorini") %>%
  select(Individual, Group, Genus, Species, Trait, Value) %>%
  filter(Trait=="intertegularDistance" | Trait=="proboscisLength") %>%
  mutate(subgenus=NA) %>%
  rename_with(tolower)

traits_eug <- read.csv2("01_raw_data/traits_eug.csv") %>%
  mutate(Genus=word(species, 1)) %>%
  rename_with(tolower)

traits <- rbind(traits_melaug, traits_eug) %>%
  mutate(family=str_extract(group, "[^:]+"), 
         tribe=word(group, -1)) %>%
  select(individual, family, tribe, genus, subgenus, species, trait, value)

rm(traits_eug, traits_melaug)

filtered_traits <- traits[traits$individual %in% 
                            traits$individual[duplicated(traits$individual)],]

write.csv2(filtered_traits, "02_processed_data/traits.csv")

mean_traits <- filtered_traits %>%
  pivot_wider(id_cols=individual:species, names_from=trait, 
              values_from=value)

mean_traits <- mean_traits %>%
  group_by(family, tribe, genus, species) %>%
  summarise(num=n(),
            mean_IT=round(mean(intertegularDistance), 2),
            sd_IT=round(sd(intertegularDistance), 2),
            mean_PB=round(mean(proboscisLength), 2),
            sd_PB=round(sd(proboscisLength), 2))

write.csv2(mean_traits, "03_final_data/mean_traits.csv")
