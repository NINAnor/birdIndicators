#' Write schedule table for RSWAN Analyses
#'
#' @param working_folder character. Path to the working folder which contains
#' the PECBMS Trim results files for data from Hekkefuglovervåkingen and the 
#' predecessor monitoring. 
#' @param general_folder character. Path to the folder containing .csv files 
#' provided by PECBMS for RSWAN analyses. This function requires the files
#' "Swan_schedules.csv" and "Swan_schedules_comb.csv". 
#' @param MSI_speciesList character vector containing EURING codes of species to
#' include in RSWAN analyses. For Norway, this includes 55 species. 
#' @param loadSchedule logical. Whether to load a previously written schedule
#' table (TRUE) or instead reassemble a schedule table based on available data
#' files (FALSE).
#'
#' @return Writes a schedule table into the general_folder directory under the
#' names "Swan_schedules.csv" and "Swan_schedules_combine.csv".
#' @export
#'
#' @examples

writeSchedule_SWAN <- function(working_folder, general_folder, MSI_speciesList, loadSchedule){
  
  ## List all present data files
  patron <- (".*\\_indices_TT.csv")
  vector_files <- list.files(working_folder, pattern = patron)
  
  
  ## Generate a dataset for storing species-country combination 
  data_sp_country <- as.data.frame(matrix(NA, ncol = 2, nrow = 1))
  names(data_sp_country) = c("Sp_code", "Country_code")
  
  for(i in 1:length(vector_files)){
    x = strsplit(vector_files[i], "_")
    
    if(x[[1]][2] == 1){
      data_sp_country[i, 1] <- x[[1]][1]
      data_sp_country[i, 2] <- x[[1]][3]
    }
  }
  
  
  ## Drop species that are not part of the RSWAN analysis
  dropIdx <- which(!(data_sp_country$Sp_code %in% MSI_speciesList))
  data_sp_country <- data_sp_country[-dropIdx, ]
  
  
  ## Set up schedule table
  if(loadSchedule){
    
    # Load pre-written schedule table
    schedule.data <- read.table(paste0(general_folder,"/Swan_schedules.csv"), header = T, sep = ";", dec = ",")
    
  }else{
    
    # Load schedule primer
    schedule.primer <- read.table(paste0(general_folder,"/Swan_schedules_comb.csv"),header = T,sep = ";",dec=",")
    
    # Fill country and species codes into primer
    pp <- unique(data_sp_country$Sp_code)
    xx <- as.character(rep(pp, each = nrow(schedule.primer)))
    
    schedule.data <- as.data.frame(matrix(ncol = ncol(schedule.primer), 
                                              nrow = length(xx), NA))
    names(schedule.data) <- names(schedule.primer)
    
    for(i in 1:ncol(schedule.primer)){
      schedule.data[, i] <- rep(schedule.primer[, i], length(pp))
    }
    schedule.data$Species_nr <- xx

    # Check country codes and update column B_sel
    for (i in 1:nrow(schedule.data)) {
      sp.subset <- subset(data_sp_country, data_sp_country$Sp_code == schedule.data$Species_nr[i])
      
      for (j in 1:nrow(sp.subset)) {
        if(schedule.data$Level1[i] == sp.subset$Country_code[j]){
          schedule.data$B_sel[i] = 1
        }
      }
    }
  }
  
  
  ## Write schedule to CSV
  write.table(schedule.data, paste0(general_folder, "/Swan_schedules.csv"), row.names = F, col.names = T, sep = ";", dec = ".")
  write.table(schedule.data, paste0(general_folder, "/Swan_schedules_combine.csv"), row.names = F, col.names = T, sep = ";", dec = ".")
}

