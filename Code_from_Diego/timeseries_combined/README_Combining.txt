Combining time series
created: Agust 2nd 2023

A) information needed:

A.1) We need the trim indices from the 'old' and the 'new' schemes. The files from the old scheme are stored somwhere in NINA's database, perhaps John Atle has them (Diego also has them). 
These old indices are not changing as the monitoring schem was discontinued and thus we can use always the same files. The files from the 'new' scheme (hekkefuglovevåking, formerly TOV-E) are calculated every year using rtrim-shell (R package) as part of the breeding bird monitoring scheme and reportinf to PECBMS.
The trim files needed are: "_arg_output.csv", "_indices_TT.csv", "_ocv.csv". Files from the old scheme have a code == 62; whereas the new scheme has a code == 63.

A.2) Country population sizes. The estimates of population size has to overlap with atleast on of the time series to be merged. 
The estmiates used here are extracted from "BirdLife International 2021: European Red List of the Birds" (NOTE: this reference was given by PECBMS).


B) We need two folders:

B.1) geenral_folder: contains steering files and four R scripts for combination. 
Steering files have been provided by PECBMS and should not be modified to make sure that all countries use the same 'base information'. KEEP ALL FILES!
 - Countries.csv: list of countries and their country codes (schemes and scheme codes)
 - Pools.csv: list of supranational codes; Norway (combined) == 14
 - Schedules.csv: list of final (combined) schedule name, we use “COMB“ = old +new combination
 - Species.csv: list of species (Euring code, Latin name, English name)
 - Species_Countries.csv: list of country population sizes for species to be combined
 - Swan_function_list.csv: overview of functions to set for calculation
 - Swan_schedules_comb.csv: template for calculation schedule, must be prepared
 - Swan_schedules.csv: calculation schedule for all species automatically prepared from a template
 - New_notorics_2022.xlsx = list of species which index should be modified (for PECBMS needs but keep it)


And the functions (RSWAN) needed for combination are:
 - 00_Code_for_B_sel_schedules = multiplies the computation schedules for all species from a template, takes SWAN_schedule_comb.csv -> creates SWAN_schedule.csv

 - 01_CODE_adding_FirstSurvey_YEAR_OK = fills in the first and last year of each species index into Species_Countries.csv (add new columns Year_first and Year_last at the end)

 - 02_run_RSWAN_v7 = runs the RSWAN calculation (combination), operates the 4th script 03_Swan_Euromonitoring_v7

 - 03_Swan_Euromonitoring_v7 = calculation procedure itself


Before using the functions, you need to prepare the data, which is done with the code in "Data_preparation_combining_timeseries.R"


 
B.2) working_folder: contains all data files (outputs of RTRIM-shell) for both schemes differentiated by scheme code (old scheme == 62; new scheme == 63)
Example of the six files we need for merging for Anas platyrhynchos (species EURING code == 1860), Norway old == 62 and Norway new = 63: 
1860_1_62_arg_output.csv
1860_1_62_indices_TT.csv
1860_1_62_ocv.csv
1860_1_63_arg_output.csv
1860_1_63_indices_TT.csv
1860_1_63_ocv.csv

NOTE: Forget about the "_1_" in the name file, which is used internally by PECBMS for some checks.