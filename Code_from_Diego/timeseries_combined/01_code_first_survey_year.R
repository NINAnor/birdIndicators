## This script requires the package dplyr or tidyverse, and the readxl to load the table of exemptions
library(dplyr)
library(readxl)

## We load the data from the backup folder, to prevent that some prior analyses affects the outcome
species_country<-read.csv2(paste0(wd, '/general+working folder/general_folder/Species_Countries.csv'),header = T, sep = ";",check.names = FALSE)


# Now we will add the exemption included in the table for those species which normally need to be truncated.
# None in our case, but has to be done...
exemptions<-read_excel(paste0(wd, '/general+working folder/general_folder/New_notorics_2022.xlsx'))


# This loop currently uses the combination of species code and country code,
# which are in the same row of the species_country file. This make the loop
# more robust since is not double or dependent of the position of the data. 

for (i in 1:nrow(species_country)){
  
  try({
    x<-read.csv(paste0(paste0(wd, '/general+working folder/working_folder/'),species_country$Species_nr[i],"_1_",species_country$Country_code[i],"_indices_TT.csv"),header = T, sep = ";")
    species_country$Year_first[i]<-min(x$Year)
    if(max(x$Year)>2022){ # This has to be change every year with the latest year. It should be easy by using the info of 'year_last'
      species_country$Year_last[i]=2022 # This has to be change every year with the latest year.
    }else{
      species_country$Year_last[i]<-max(x$Year)
    }
  },silent = T)
}


# This double loop looks through both files: the table with exemptions and the species_country file, 
# searching for a double match. 
# First it matches the species code.
# Second, it matches the country code, and if both matches are true, it replaces 
# the first Year_first previously filled in the species_country file with the year first.
# Additionally it prints for which species and which countries has the year changed.
# (Mostly to ensure it is working.)

for(i in 1:nrow(exemptions)){
  for (j in 1:nrow(species_country)) {
    try({ 
      if (exemptions$Species_name[i]== species_country$Species_nr[j]){
        if (exemptions$Country_code[i]== species_country$Country_code[j]){
          # print(paste0("the species ",exemptions$SciName[i]," in ",exemptions$Country[i], " has changed from ",species_country$Year_first[j] ," to ", exemptions$Year_First[i]))
          print(paste0("the species ",exemptions$Sci_name[i]," in ",exemptions$Country_name[i], " has changed from ",species_country$Year_last[j] ," to ", exemptions$Year_last[i]))
          
          species_country$Year_first[j]<-exemptions$Year_first[i]
          species_country$Year_last[j]<-exemptions$Year_last[i]
          
        }
      }
    },silent = T)
  }
}


str(species_country)
species_country$`Population size (geomean)`=as.numeric(species_country$`Population size (geomean)`)
species_country$`Population size_min`=as.numeric(species_country$`Population size_min`)
species_country$`Population size_max`=as.numeric(species_country$`Population size_max`)
species_country$Weight
# Now we are checking whether our loops worked properly or not. 
# Each loop is written to give you an error message locating the indices TT folder or 
# the species and the problematic years
# This loop compares the value of Year_first for each species and Year located in first position [1] in its particular indices_TT files
#The only printed errors should be the pairs "species&countries" within the exemption files.

#for (i in 1:nrow(species_country)) {
#  try({x<-read.csv(paste0("D:\\WORK\\RSWAN\\RSWAN_tests\\02_2020_PECBMS_data_till2017\\03_EUR\\V5\\working_folder\",species_country$Species_nr[i],"_1_",species_country$Country_code[i],"_indices_TT.csv"),header = T, sep = ";")
#  if (species_country$Year_first[i]!= x$Year[1]){
#print(paste0("Error in", species_country$Species_nr[i],"_1_",species_country$Country_code[i],"_indices_TT.csv" ))
#    print(paste0("The correct year for", species_country$SciName[i], " in ", species_country$Country[i], " is ", x$Year[1], " not ", species_country$Year_first[i]))

#     }  
#  },silent = T)
#}

#This loop compare the value of Year_first for each species and Year located in first position [nrow(x)] in its particular indices_TT files
#The only printed errors should be the pairs "species&countries" within the exemption files.
#for (i in 1:nrow(species_country)) {
#  try({x<-read.csv(paste0("D:\\WORK\\RSWAN\\RSWAN_tests\\02_2020_PECBMS_data_till2017\\03_EUR\\V5\\working_folder\",species_country$Species_nr[i],"_1_",species_country$Country_code[i],"_indices_TT.csv"),header = T, sep = ";")
#    if (species_country$Year_last[i]!= x$Year[nrow(x)]){
#      #print(paste0("Error in", species_country$Species_nr[i],"_1_",species_country$Country_code[i],"_indices_TT.csv" ))
#      print(paste0("The correct year for", species_country$SciName[i], " in ", species_country$Country[i], " is ", x$Year[nrow(x)], " not ", species_country$Year_last[i]))
#    }
#  },silent = T)
#}



#This exports the Species_Country file, directly to the general folder 
write.csv2(species_country,paste0(wd, '/general_folder/Species_Countries.csv'),row.names = F)

