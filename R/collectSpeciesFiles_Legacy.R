#' Collects legacy species and summary files in designated folder
#'
#' @param origin_folder character. Path to directory which contains the legacy
#' RTRIM output files (results of PECBMS Trim analyses of pre-2006 data).
#' @param target_folder character. Name of the directory into which to copy the
#' legacy files. Has to be the same directory that also holds the most recent
#' RTRIM output files from analysis of up-to-date data from Hekkefuglovervåkingen. 
#'
#' @return This function returns no object but copies filed from one folder to
#' another. 
#' @export
#'
#' @examples

collectSpeciesFiles_Legacy <- function(origin_folder, target_folder){
  
  ## List names of relevant files
  species_files <- dir(origin_folder, pattern = ".csv")

  ## Check whether files need to be moved
  species_files_target <- dir(target_folder, pattern = ".csv")
  
  if(all(species_files %in% species_files_target)){
    
    message("Legacy files are already present in target folder. Copying is not neccessary.")
  
  }else{
    
    ## File copy function
    copyFiles <- function(x){
      file.copy(from = file.path(origin_folder, x) ,
                to = file.path(target_folder, x))
    }
    
    ## Copy all files
    lapply(species_files, copyFiles)

    ## Summarise actions
    message(paste0("TRIM Legacy files have been copied to: ", target_folder, "."))
    
  }
  
}