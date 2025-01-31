## ---------------------------------------------------------------------##
## Working script to modify the pipeline to generate indices by regions ##
## ---------------------------------------------------------------------##

#setwd("~/NINA ONGOING projects/2023_Pipeline_for_data_access/birdIndicators/data")

require(tidyverse)
library(read.table)

sppRegions <- read.csv2('~/NINA ONGOING projects/2023_Pipeline_for_data_access/birdIndicators/data/Hekkefugl_NI2020_SummaryUTF8.csv', dec = '.', sep = ',', header =  T)

sppTable_ALL <- readRDS("C:\\Users\\diego.pavon-jordan\\OneDrive - NINA\\Documents\\NINA ONGOING projects\\2023_Pipeline_for_data_access\\birdIndicators\\data\\species_table_NI.rds")


sppTable_ALL # from 'makeSpeciesLists.R'
sppTable_ALL %>% filter(Species == 'Pyrrhula pyrrhula')
sppRegions %>% filter(Species == 'Pyrrhula pyrrhula')


# select the regions for the species included in the Naturindeks
sppTable_ALL_NI <- sppTable_ALL %>%
  mutate(NI = case_when(Species == 'Tetrao tetrix' ~ TRUE,
                        Species == 'Tetrao urogallus' ~ TRUE,
                        Species == 'Lagopus lagopus' ~ TRUE,
                        Species == 'Lagopus muta' ~ TRUE,
                        TRUE ~ NI)) %>%
  mutate(StartDataHFT = case_when(Species == 'Tetrao tetrix' ~ '2008',
                        Species == 'Tetrao urogallus' ~ '2008',
                        Species == 'Lagopus lagopus' ~ '2008',
                        Species == 'Lagopus muta' ~ '2008',
                        TRUE ~ StartDataHFT)) %>%
  left_join(sppRegions)


sppTable_ALL_NI %>% filter(Species == 'Pyrrhula pyrrhula')

sppTable_ALL_NI %>% filter(NI == TRUE) %>% print(n=60)

sppTable_ALL_NI %>% filter(dataUse_HeleNorge == 1) %>% print(n=60)

saveRDS(sppTable_ALL_NI, 'C:\\Users\\diego.pavon-jordan\\OneDrive - NINA\\Documents\\NINA ONGOING projects\\2023_Pipeline_for_data_access\\birdIndicators\\data\\species_table_NI.rds')


# Now, the table with the species information is updated and contains the
# four grous species and the info on which regions to use.


## Now we need to modify the pipeline to do the different calculations by region 
## for the NI



