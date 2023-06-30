### Processing Data from Direct Litter and Lateral litter inputs data
### Danielle Harre
### 11/1/2021
### Datasource: G:/My Drive/CREWS_official/300_wholestream/310_data_wholestream/321_organicmatterbudget

###
####
#Read RAW Data
file_direct <- "data/direct_OCinputs_CoweetaNC.csv"
file_lateral <- "data/lateral_OCinputs_CoweetaNC.csv"
POM_direct <- read.csv(file_direct)
POM_direct$sample.date <- as.Date(POM_direct$sample.date, "%Y-%m-%d")
POM_lateral <- read.csv(file_lateral)
POM_lateral$sample.date <- as.Date(POM_lateral$sample.date, "%Y-%m-%d")

### Process Data
#Summarize direct data into avg, but keep two basins separate
PD_day <- POM_direct%>%
  filter(organic.matter.fraction == "LEAF") %>%
  dplyr::select("sample.date", "sample.mo", "trap.no","stream", "areal.AFDM.input.g.m2.d")%>%
  mutate(Jdate = yday(sample.date)) %>%
  group_by(Jdate) %>%
  summarise(
    direct_gm2d_avg = round(mean(areal.AFDM.input.g.m2.d),2),
    direct_gm2d_sd = round(sd(areal.AFDM.input.g.m2.d),2)
  )%>%
  drop_na(direct_gm2d_avg)


saveRDS(PD_day, "./data/POM_direct_day.RDS")

#Summarize lateral data into avg, but keep basins separate
PL_day <- POM_lateral%>%
  filter(organic.matter.fraction == "LEAF") %>%
  dplyr::select("sample.date", "sample.mo", "trap.no","stream", "areal.AFDM.input.g.m.d")%>%
  mutate(Jdate = yday(sample.date)) %>%
  group_by(Jdate) %>%
  summarise(
    lateral_gmd_avg = round(mean(areal.AFDM.input.g.m.d),2),na.rm=TRUE,
    lateral_gmd_sd = round(sd(areal.AFDM.input.g.m.d),2),na.rm=TRUE
  )%>%
  drop_na(lateral_gmd_avg)

saveRDS(PL_day, "./data/POM_lateral_day.RDS")


