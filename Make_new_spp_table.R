## ---------------------------------------------------------------------##
## Working script to modify the pipeline to generate indices by regions ##
## ---------------------------------------------------------------------##

#setwd("~/NINA ONGOING projects/2023_Pipeline_for_data_access/birdIndicators/data")

require(tidyverse)
library(read.table)

#sppRegions <- read.csv2('~/NINA ONGOING projects/2023_Pipeline_for_data_access/birdIndicators/data/Hekkefugl_NI2020_SummaryUTF8.csv', dec = '.', sep = ',', header =  T)

#sppTable_ALL <- readRDS("C:\\Users\\diego.pavon-jordan\\OneDrive - NINA\\Documents\\NINA ONGOING projects\\2023_Pipeline_for_data_access\\birdIndicators\\data\\species_table_NI.rds")

sppRegions <- read.csv2("data/Hekkefugl_NI2020_SummaryUTF8.csv", dec = '.', sep = ',', header =  T) %>%
  dplyr::mutate(dataUse_expert_N = ifelse(dataUse_direct_N == "Expert", TRUE, FALSE),
                dataUse_direct_SN = ifelse(dataUse_direct_SN == "x", TRUE, FALSE),
                dataUse_direct_NN = ifelse(dataUse_direct_NN == "x", TRUE, FALSE),
                dataUse_direct_N = ifelse(dataUse_direct_N == "x", TRUE, FALSE))

sppTable_ALL <- readRDS("data/species_table.rds")

sppTable_ALL # from 'makeSpeciesLists.R'

# Add missing NI species (Cinclus cinclus and expert-judged species)
expertJudge <- c("Acrocephalus schoenobaenus",
                 "Anthus petrosus",
                 "Aythya marila",
                 "Calidris alpina",
                 "Charadrius morinellus",
                 "Clangula hyemalis",
                 "Eremophila alpestris",
                 "Gallinago media",
                 "Melanitta fusca",
                 "Melanitta nigra",
                 "Phalaropus lobatus",
                 "Plectrophenax nivalis")

expertJudge <- expertJudge[which(!(expertJudge %in% sppTable_ALL$Species))]

sppTable_additions <- tibble(
  EURINGCode = c(10500, # Cinclus cinclus
                 12430, # Acrocephalus schoenobaenus
                 10142, # Anthus petrosus
                 2040, # Aythya marila
                 5120, # Calidris alpina
                 4820, # Charadrius morinellus
                 2120, # Clangula hyemalis
                 9780, # Eremophila alpestris
                 5200, # Gallinago media
                 2150, # Melanitta fusca
                 2130, # Melanitta nigra
                 5640), # Phalaropus lobatus
  Species = c("Cinclus cinclus", expertJudge),
  WEBSITE = c(TRUE, rep(FALSE, length(expertJudge))),
  MSI = FALSE,
  PECBMS = FALSE,
  NI = TRUE
)

sppTable_ALL <- sppTable_ALL %>%
  dplyr::bind_rows(sppTable_additions)

sppTable_ALL %>% filter(Species == 'Pyrrhula pyrrhula')
#sppRegions %>% filter(Species == 'Pyrrhula pyrrhula') #CRN: There is no species column in sppRegions

# CRN: I guess I have no choice but retrieve the latin names and merge them in
NameList <- readxl::read_xlsx("data/NNKF142-20241111komplett.xlsx") %>%
  dplyr::select(IOC14.2, innhentet141) %>%
  dplyr::rename(Species = IOC14.2,
                indicatorName = innhentet141) %>%
  dplyr::mutate(indicatorName = stringr::str_to_sentence(indicatorName))

# Change norwegian names for species names with different spelling in NI
sppRegions <- sppRegions %>%
  dplyr::mutate(indicatorName = dplyr::case_when(indicatorName == "Grå fluesnapper" ~ "Gråfluesnapper",
                                                 indicatorName == "Fiskemåke ferskvann" ~ "Fiskemåke",
                                                 TRUE ~ indicatorName))
# Match latin names
sppTable_ALL <- sppTable_ALL %>%
  dplyr::left_join(NameList, by = "Species")

# Check for NAs and fix manually
sppTable_ALL %>%
  dplyr::filter(is.na(indicatorName))

sppTable_ALL <- sppTable_ALL %>%
  dplyr::mutate(indicatorName = dplyr::case_when(!is.na(indicatorName) ~ indicatorName,
                                                 Species == "Tetrao tetrix" ~ "Orrfugl",
                                                 Species == "Sylvia communis" ~ "Tornsanger",
                                                 Species == "Sylvia curruca" ~ "Møller",
                                                 Species == "Carduelis chloris" ~ "Grønnfink",
                                                 Species == "Carduelis spinus" ~ "Grønnsisik",
                                                 Species == "Carduelis flammea" ~ "Gråsisik",
                                                 Species == "Charadrius morinellus" ~ "Boltit"))



# Add information on select the regions for the species included in the Naturindeks
sppTable_ALL_NI <- sppTable_ALL %>%
  left_join(sppRegions)


# Check unmatched entries from sppRegions
sppRegions$indicatorName[which(!(sppRegions$indicatorName) %in% sppTable_ALL_NI$indicatorName)]
sppRegions[which(!(sppRegions$indicatorName) %in% sppTable_ALL_NI$indicatorName),]
# --> All species have a match

sppTable_ALL_NI %>% filter(Species == 'Pyrrhula pyrrhula')

sppTable_ALL_NI %>% filter(NI == TRUE) %>% print(n=100)

sppTable_ALL_NI %>% filter(dataUse_direct_N) %>% print(n=60) #CRN: This column dataUse_HeleNorge is something I do not have, but assuming it's the same as dataUse_direct_N



#saveRDS(sppTable_ALL_NI, 'C:\\Users\\diego.pavon-jordan\\OneDrive - NINA\\Documents\\NINA ONGOING projects\\2023_Pipeline_for_data_access\\birdIndicators\\data\\species_table_NI.rds')
saveRDS(sppTable_ALL_NI, file = "data/species_table_NI.rds")
readr::write_excel_csv(sppTable_ALL_NI, file = "data/species_table_NI.csv")

# Now, the table with the species information is updated and contains the
# four grouse species and the info on which regions to use.


## Now we need to modify the pipeline to do the different calculations by region 
## for the NI



