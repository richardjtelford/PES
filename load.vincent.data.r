#load libraries
library("vegan")
library("tidyverse")
library("entropy")
library("lubridate")
library("assertr")
library("assertthat")

#functions
div <- function(x) {
  apply(x, 1, entropy.ChaoShen)
}

cleanStationCodes <- function(codes){
  codes <- trimws(codes) # remove whitespace
  codes <- gsub("_", "", codes) # remove underscores
  codes <- gsub("^0", "", codes) # remove leading zeros
  codes <- recode(codes, R060 = "R60", R080 = "R80", GRO = "G69", GRO50 = "G50", GRO40 = "G40", GRO60 = "G60") # replace with alternative code
  codes
}

#species processing functions
sppProcess <- . %>%
  summarise(count = sum(N)) %>%
  spread(key = species, value = count, fill = 0)
group_replicate <- . %>% 
  group_by(Station_code, species, Replicate) %>% 
  sppProcess()
group_station <- . %>% 
  group_by(Station_code, species) %>% 
  sppProcess()
group_station_minCount <- . %>% 
  group_by(Station_code, species) %>% 
  filter(n > minCount) %>% 
  sppProcess()

### access data from database
con <- src_sqlite("data/PES_DB_8Feb2011.sqlite", create = FALSE)

##stations
stations <- tbl(con, "PES_DB_stations") %>%
  rename(Station_code = STATION) %>% 
  collect() %>% 
  mutate(Station_code = cleanStationCodes(Station_code))


##macrofauna
macro <- tbl(con, "PES_DB_macrofauna_at_forams_stations") %>%
  rename(species = CodeMacrofauna, N = `Number*1`) %>% 
  collect() %>%
  mutate(DAT = ymd(DAT)) %>% 
  mutate(Station_code = cleanStationCodes(Station_code))

#replicate level
macro3r <- macro %>% 
  filter(DAT < "2005-01-01") %>% 
  group_replicate() 

RC  <- c("RC5", "RC9")# sites with multiple samples.
macro8r <- macro %>% 
  filter(DAT > "2005-01-01") %>%
  filter((Station_code %in% RC & month(DAT) == 9) | !Station_code %in% RC) %>% #only September samples from RC5/9
  group_replicate()

#site level
macro3g <- macro %>% 
  filter(DAT < "2005-01-01")

macro8g <- macro %>% 
  filter(DAT > "2005-01-01") %>% 
  filter((Station_code %in% RC & month(DAT) == 9) | !Station_code %in% RC) #only September samples from RC5/9

##live forams
forams <- tbl(con, sql("select * from PES_DB_foraminifera_species_data where slice_numeric>0 and not Size  = '>500' and slice_numeric<3")) %>%
  rename(species = SpeciesForam, N = `Number*1`) %>% 
  collect() %>%
  mutate(DAT = ymd(DAT)) %>% 
  mutate(Station_code = cleanStationCodes(Station_code))

#replicate level
foram3r <- forams %>% 
  filter(DAT < "2005-01-01") %>% 
  group_replicate()

foram8r <- forams %>% 
  filter(DAT > "2005-01-01") %>% 
  filter((Station_code %in% RC & month(DAT) == 9) | !Station_code %in% RC) %>% #only September samples from RC5/9
  group_replicate()

#site level
foram3g <- forams %>% 
  filter(DAT < "2005-01-01") 

foram8g <- forams %>% 
  filter(DAT > "2005-01-01") %>% 
  filter((Station_code %in% RC & month(DAT) == 9) | !Station_code %in% RC) #only September samples from RC5/9

##dead forams
dead <- tbl(con, "PES_DB_foraminifera_species_dead_data") %>%
  rename(N = `Number*1`) %>%
  collect() %>% 
  mutate(Station_code = cleanStationCodes(Station_code))

##foram species
fspp <- tbl(con,"PES_DB_foraminifera_species_list") %>%
  rename(TestStructure = `Test structure`, Sensitivity = `Sensitivity (for ISI)`) %>%
  collect()

##chemistry
#SQL was "SELECT STAS, Chemical_species, Avg(IIf([Value]<-98,0,[value])) AS Val FROM PES_db_chemistry_data_flat WHERE (((DATE)>20050000) AND ((Slice_numeric)<1 Or (Slice_numeric) Is Null)) GROUP BY STAS, Chemical_species"

chem0 <- tbl(con, "PES_db_chemistry_data_flat") %>% 
  filter(DATE > 20050000, Slice_numeric < 1 | is.na(Slice_numeric)) %>%
  filter(!Chemical_species %in% c("Zn","Cu","Pb", "Cd","H2S", "Chl c3")) %>% #unwanted variables
  rename(Station_code = STAS) %>%
  group_by(Station_code, Chemical_species) %>%
  mutate(Value = ifelse(Value < -98, 0, Value)) %>% #recode -99s as zero
  summarise(Val = mean(Value), min = min(Value)) %>% # mean chemistry per site
  mutate(Val = ifelse(Chemical_species == "O2", min, Val)) %>%
  select(-min) %>%
  collect() %>% 
  ungroup() %>% 
  mutate(Station_code = cleanStationCodes(Station_code)) %>%
  group_by(Station_code, Chemical_species) %>% 
  filter(!Station_code %in% c("HV16", "KV01", "KRG")) # remove unwanted chemistry stations



pigments <- c("allo-xanthin","aphanizophyll","beta-carotene","cantha-xanthin","chl a","chl a total (a+allom)","Chl b","Chl c2","diadino-xanthin","diato-xanthin","fuco-xanthin","lutein","peridinin","pheo-phytin a","Pheophytin b","Pyropheophytin b", "violaxanthin","zea-xanthin")

min0 <- function(x) {
  min(x[x > 0], na.rm = TRUE)
}

#unified site list
allSites <- list(
  chem = unique(chem0$Station_code),
  macro = unique(macro8g$Station_code),
  foram = unique(foram8g$Station_code)
)
unified_site_list <- Reduce(intersect, allSites) # sites with all variables # would miss abiotic sites

chem0 <- chem0 %>%
  ungroup() %>%
  filter(Station_code %in% unified_site_list) %>%# only sites with 2008 macro AND foram data
  filter(!Chemical_species %in% c("Pheophorbide", "Chl a allomer")) %>% #too many missing values
  left_join(
    y = filter(., Chemical_species == "TOC") %>% select(-Chemical_species) %>% rename(TOC = Val),
    by = c("Station_code" = "Station_code")
  ) %>% #left join to TOC data
  mutate(Val = ifelse(Chemical_species %in% pigments, Val/TOC, Val)) %>%#standardise pigments by TOC
  group_by(Chemical_species) %>%
  select(-TOC) 
  

O2 <- chem0 %>% filter(Chemical_species %in% c("O2", "MinO2_2_years")) %>% 
  spread(key = Chemical_species, value = Val) %>%
  mutate(Val = ifelse(is.na(MinO2_2_years), O2, MinO2_2_years)) %>% 
  mutate(Chemical_species = "mO2") %>%
  select(-O2, -MinO2_2_years)

chem0 <- bind_rows(chem0, O2)

#remove unwanted chemistry sites
#HV16 and KV01 have no O2 data. Not sure why KRG was deleted
chem0 <- chem0 %>%
  filter(!Chemical_species %in% c("MinO2_2_years", "mO2"))

#spread

chem <- spread(chem0, key = Chemical_species, value = Val) %>% 
  mutate(ppna = `pheo-phytin a` + `chl a total (a+allom)`) %>%
  select(-`pheo-phytin a`, -`chl a total (a+allom)`) %>% 
  select(Station_code, `%<63`, `allo-xanthin`,`beta-carotene`,`diato-xanthin`,lutein, O2, TN, TOC,`zea-xanthin`, ppna) %>% #"Cd","Cu""Pb","Zn"
  assert(not_na, -Station_code) #check no NAs

chem <- chem %>% 
  left_join(
    stations %>% 
      filter(Station_code %in% unified_site_list) %>% 
      select(Station_code, DEPTH_BELOW_THRESHOLD) %>%
      rename(`Depth Below Threshold` = DEPTH_BELOW_THRESHOLD)
  ) %>% 
  assert(not_na, -Station_code)

pigments <- c(pigments, "ppna")
chem[names(chem) %in% pigments] <- sapply(chem[names(chem) %in% pigments], function(x){
    log(x + min0(x)/2) #log(x + k) where k is half minimum value
  })

#tidyup
rm(O2)

#Generate data sets with harmonised site lists f/m & f/m/c

m83 <- macro8r %>% filter(Station_code %in% macro3r$Station_code)
m38 <- macro3r %>% filter(Station_code %in% macro8r$Station_code)

foram83<-foram8r %>% filter(Station_code %in% foram3r$Station_code) 
foram38<-foram3r %>% filter(Station_code %in% foram8r$Station_code) 

macro8f<-macro8g %>% filter(Station_code %in% unified_site_list) %>% group_station()
foram8m<-foram8g %>% filter(Station_code %in% unified_site_list) %>% group_station()

assert_that(identical(rownames(macro8f), rownames(foram8m)))

#subset with > 30 individuals
minCount <- 30
macro8f30<-macro8g %>% filter(Station_code %in% unified_site_list) %>% group_station_minCount()
foram8m30<-foram8g %>% filter(Station_code %in% unified_site_list) %>% group_station_minCount()



rowSums(macro8f30[, -1])
assert_that(min(rowSums(foram8m[, -1])) > 50)


#tidy up
rm(group_station_minCount, group_station, group_replicate, cleanStationCodes)