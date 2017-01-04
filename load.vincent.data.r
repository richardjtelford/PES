#load libraries
library("vegan")
library("dplyr")
library("entropy")
library("tidyr")

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
  collect()

#replicate level
macro3r <- macro %>% 
  filter(DAT < 20050000) %>% 
  select(-DAT, -SpeciesMacrofauna) %>%
  spread(key = CodeMacrofauna, value = `Number*1`, fill = 0)
###FAILS because of RC5/9 repeated sampling
stop()
macro8r <- macro %>% 
  filter(DAT > 20050000) %>% group_by(Station_code, Replicate, CodeMacrofauna) %>% count() %>% arrange(desc(n))
  select(-DAT, -SpeciesMacrofauna) %>% 
  spread(key = CodeMacrofauna, value = `Number*1`, fill = 0)

#site level
macro3g <- macro %>% 
  filter(DAT < 20050000) %>% 
  group_by(Station_code, CodeMacrofauna) %>%
  summarise(count = sum(`Number*1`)) %>%
  spread(key = CodeMacrofauna, value = count, fill = 0)
macro8g <- macro %>% 
  filter(DAT > 20050000) %>% 
  group_by(Station_code, CodeMacrofauna) %>%
  summarise(count = sum(`Number*1`)) %>%
  spread(key = CodeMacrofauna, value = count, fill = 0)


##live forams
forams <- tbl(con, sql("select * from PES_DB_foraminifera_species_data where slice_numeric>0 and not Size  = '>500' and dat>20050000 and slice_numeric<3")) %>%
  collect()

#replicate level
foram3r <- forams %>% 
  filter(DAT < 20050000) %>% #obviously fails because of about sql
  group_by(Station_code, Replicate, SpeciesForam) %>% 
  summarise(count = sum(`Number*1`)) %>%
  spread(key = SpeciesForam, value = count, fill = 0)

foram8r <- forams %>% 
  filter(DAT > 20050000) %>% 
  group_by(Station_code, Replicate, SpeciesForam) %>% 
  summarise(count = sum(`Number*1`)) %>% # lump multiple depths
  spread(key = SpeciesForam, value = count, fill = 0)

#site level
foram3g <- forams %>% 
  filter(DAT < 20050000) %>%#fails because sql 
  group_by(Station_code, SpeciesForam) %>%
  summarise(count = sum(`Number*1`)) %>%
  spread(key = SpeciesForam, value = count, fill = 0)
foram8g <- forams %>% 
  filter(DAT > 20050000) %>% 
  group_by(Station_code, SpeciesForam) %>%
  summarise(count = sum(`Number*1`)) %>%
  spread(key = SpeciesForam, value = count, fill = 0)


##dead forams
dead <- tbl(con, "PES_DB_foraminifera_species_dead_data") %>%
 rename(Number = `Number*1`) %>%
 collect()

##foram species
fspp <- tbl(con,"PES_DB_foraminifera_species_list") %>%
  rename(fspp, TestStructure = `Test structure`, Sensitivity = `Sensitivity (for ISI)`) %>%
  collect()

##chemistry
#SQL was "SELECT STAS, Chemical_species, Avg(IIf([Value]<-98,0,[value])) AS Val FROM PES_db_chemistry_data_flat WHERE (((DATE)>20050000) AND ((Slice_numeric)<1 Or (Slice_numeric) Is Null)) GROUP BY STAS, Chemical_species"

chem0 <- tbl(con, "PES_db_chemistry_data_flat") %>% 
  filter(DATE > 20050000, Slice_numeric < 1 | is.na(Slice_numeric)) %>%
  filter(!Chemical_species %in% c("Zn","Cu","Pb", "Cd","H2S", "Chl c3")) %>% #unwanted variables
  group_by(STAS, Chemical_species) %>%
  mutate(Value = ifelse(Value < -98, 0, Value)) %>% #recode -99s as zero
  summarise(Val = mean(Value)) %>% # mean chemistry per site
  collect()

chem0 <- chem0 %>%
  ungroup() %>%
  mutate(STAS = trimws(STAS)) %>%# zap trailing spaces
  filter(STAS %in% c(macro$Station_code, forams$Station_code)) %>%# only sites with macro/foram data
  filter(!Chemical_species %in% c("Pheophorbide", "Chl a allomer")) #too many missing values

chem <- spread(chem0, key = Chemical_species, value = Val)

#standardise pigments by TOC
names(chem)         
pig <- names(chem) %in% c("allo-xanthin","aphanizophyll","beta-carotene","cantha-xanthin","chl a","chl a total (a+allom)","Chl b","Chl c2","diadino-xanthin","diato-xanthin","fuco-xanthin","lutein","peridinin","pheo-phytin a","Pheophytin b","Pyropheophytin b", "violaxanthin","zea-xanthin")

chem[, pig] <- chem[, pig] / chem$TOC

plot(chem$MinO2_2_years, chem$O2)

chem$mO2 <- chem$MinO2_2_years
chem$mO2[is.na(chem$mO2)] <- chem$O2[is.na(chem$mO2)]

