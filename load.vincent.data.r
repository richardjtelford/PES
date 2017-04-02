#load libraries
library("vegan")
library("tidyverse")
library("entropy")
library("lubridate")

#functions
div <- function(x) {
  apply(x, 1, entropy.ChaoShen)
}

### access data from database
con <- src_sqlite("data/PES_DB_8Feb2011.sqlite", create = FALSE)

##stations
stations <- tbl(con, "PES_DB_stations") %>%
  collect()

##macrofauna
macro <- tbl(con, "PES_DB_macrofauna_at_forams_stations") %>%
  collect() %>%
  mutate(DAT = ymd(DAT))

#replicate level
macro3r <- macro %>% 
  filter(DAT < "2005-01-01") %>% 
  select(-DAT, -SpeciesMacrofauna) %>%
  spread(key = CodeMacrofauna, value = `Number*1`, fill = 0)

RC  <- c("RC5", "RC9")# sites with multiple samples.
macro8r <- macro %>% 
  filter(DAT > "2005-01-01") %>%
  filter((Station_code %in% RC & month(DAT) == 9) | !Station_code %in% RC) %>% #only September samples from RC5/9
  select(-DAT, -SpeciesMacrofauna) %>% 
  spread(key = CodeMacrofauna, value = `Number*1`, fill = 0)

#site level
macro3g <- macro %>% 
  filter(DAT < "2005-01-01") %>% 
  group_by(Station_code, CodeMacrofauna) %>%
  summarise(count = sum(`Number*1`)) %>%
  spread(key = CodeMacrofauna, value = count, fill = 0)

macro8g <- macro %>% 
  filter(DAT > "2005-01-01") %>% 
  filter((Station_code %in% RC & month(DAT) == 9) | !Station_code %in% RC) %>% #only September samples from RC5/9
  group_by(Station_code, CodeMacrofauna) %>%
  summarise(count = sum(`Number*1`)) %>%
  spread(key = CodeMacrofauna, value = count, fill = 0)


##live forams
forams <- tbl(con, sql("select * from PES_DB_foraminifera_species_data where slice_numeric>0 and not Size  = '>500' and slice_numeric<3")) %>%
  collect() %>%
  mutate(DAT = ymd(DAT))

#replicate level
foram3r <- forams %>% 
  filter(DAT < "2005-01-01") %>% 
  group_by(Station_code, Replicate, SpeciesForam) %>% 
  summarise(count = sum(`Number*1`)) %>%
  spread(key = SpeciesForam, value = count, fill = 0)

foram8r <- forams %>% 
  filter(DAT > "2005-01-01") %>% 
  filter((Station_code %in% RC & month(DAT) == 9) | !Station_code %in% RC) %>% #only September samples from RC5/9
  group_by(Station_code, Replicate, SpeciesForam) %>% 
  summarise(count = sum(`Number*1`)) %>% # lump multiple depths
  spread(key = SpeciesForam, value = count, fill = 0)

#site level
foram3g <- forams %>% 
  filter(DAT < "2005-01-01") %>%
  group_by(Station_code, SpeciesForam) %>%
  summarise(count = sum(`Number*1`)) %>%
  spread(key = SpeciesForam, value = count, fill = 0)

foram8g <- forams %>% 
  filter(DAT > "2005-01-01") %>% 
  filter((Station_code %in% RC & month(DAT) == 9) | !Station_code %in% RC) %>% #only September samples from RC5/9
  group_by(Station_code, SpeciesForam) %>%
  summarise(count = sum(`Number*1`)) %>%
  spread(key = SpeciesForam, value = count, fill = 0)


##dead forams
dead <- tbl(con, "PES_DB_foraminifera_species_dead_data") %>%
 rename(Number = `Number*1`) %>%
 collect()

##foram species
fspp <- tbl(con,"PES_DB_foraminifera_species_list") %>%
  rename(TestStructure = `Test structure`, Sensitivity = `Sensitivity (for ISI)`) %>%
  collect()

##chemistry
#SQL was "SELECT STAS, Chemical_species, Avg(IIf([Value]<-98,0,[value])) AS Val FROM PES_db_chemistry_data_flat WHERE (((DATE)>20050000) AND ((Slice_numeric)<1 Or (Slice_numeric) Is Null)) GROUP BY STAS, Chemical_species"

chem0 <- tbl(con, "PES_db_chemistry_data_flat") %>% 
  filter(DATE > 20050000, Slice_numeric < 1 | is.na(Slice_numeric)) %>%
  filter(!Chemical_species %in% c("Zn","Cu","Pb", "Cd","H2S", "Chl c3")) %>% #unwanted variables
  group_by(STAS, Chemical_species) %>%
  mutate(Value = ifelse(Value < -98, 0, Value)) %>% #recode -99s as zero
  summarise(Val = mean(Value), min = min(Value)) %>% # mean chemistry per site
  mutate(Val = ifelse(Chemical_species == "O2", min, Val)) %>%
  select(-min) %>%
  collect()

pigments <- c("allo-xanthin","aphanizophyll","beta-carotene","cantha-xanthin","chl a","chl a total (a+allom)","Chl b","Chl c2","diadino-xanthin","diato-xanthin","fuco-xanthin","lutein","peridinin","pheo-phytin a","Pheophytin b","Pyropheophytin b", "violaxanthin","zea-xanthin")

chem0 <- chem0 %>%
  ungroup() %>%
  mutate(STAS = trimws(STAS)) %>%# zap trailing spaces
  filter(STAS %in% c(macro$Station_code, forams$Station_code)) %>%# only sites with macro/foram data
  filter(!Chemical_species %in% c("Pheophorbide", "Chl a allomer")) %>% #too many missing values
  left_join(
    y = filter(., Chemical_species == "TOC") %>% select(-Chemical_species) %>% rename(TOC = Val),
    by = c("STAS" = "STAS")
  ) %>%
  mutate(Val = ifelse(Chemical_species %in% pigments, Val/TOC, Val)) %>%#standardise pigments by TOC
  select(-TOC)

O2 <- chem0 %>% filter(Chemical_species %in% c("O2", "MinO2_2_years")) %>% 
  spread(key = Chemical_species, value = Val) %>%
  mutate(Val = ifelse(is.na(MinO2_2_years), O2, MinO2_2_years)) %>% 
  mutate(Chemical_species = "mO2") %>%
  select(-O2, -MinO2_2_years)

chem0 <- bind_rows(chem0, O2)

chem <- spread(chem0, key = Chemical_species, value = Val)
