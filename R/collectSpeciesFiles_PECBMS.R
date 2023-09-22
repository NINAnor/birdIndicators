#' Collects species and summary files in designated folder
#'
#' @param folder character. Path to directory which contains the RTRIM argument
#' and output files, and into which to save outputs.
#' @param subFolderName character. Desired name for sub-directory in which to 
#' collect renamed files. 
#'
#' @return This function returns no object but makes a new folder and copies
#' species and summary files into it. 
#' @export
#'
#' @examples

collectSpeciesFiles_PECBMS <- function(folder, subFolderName){

  ## Create new subfolder
  suppressWarnings(dir.create(file.path(paste0(folder, "/", subFolderName))))
  
  ## List names of files to move
  species_files <- dir(folder, pattern = ".csv")
  summary_tables <- dir(folder, pattern = ".txt")
  
  ## File copy function
  copyFiles <- function(x){
    file.copy(from = file.path(folder, x) ,
              to = file.path(paste0(folder, "/", subFolderName), x) )
  }
  
  ## Copy all files
  lapply(species_files, copyFiles)
  lapply(summary_tables, copyFiles)
  
  
  ## Remove "BMP_" from file names
  Old_names <- list.files(paste0(folder, "/", subFolderName))
  
  New_names <- Old_names %>%
    stringr::str_replace('BMP_', '')
  
  file.rename(paste0(folder, "/", subFolderName, "/", Old_names), 
              paste0(folder, "/", subFolderName, "/", New_names))
  
  ## Summarise actions
  message(paste0("Species and summary files have been saved in: ", folder, "/", subFolderName, "."))
}