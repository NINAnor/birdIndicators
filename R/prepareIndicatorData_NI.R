#' Prepare updated NI indicator data based on (area-specific) Trim analysis outputs
#'
#' @param sppNI a table containing information on species included in NI.
#' @param OldIndicator_data a list containing the NI indicator data currently 
#' stored in the NI database for each species in `sppNI`. 
#' @param updatePlans a list containing one table per species detailing which
#' Trim analysis output (area) should be used to update NI indicator values for
#' each area in the NI database. Output of makeUpdatePlan_NI(). 
#' @param inputFile_folders a data frame storing information on the filepaths
#' for the different subsets of Trim analysis outputs.
#' @param use_combTS logical. If TRUE, uses combined Trim time series (RSWAN runs)
#' if available for the selected species. This is so far only implemented for 
#' species for which Trim analyses are run at the national level (the sub-national
#' level Trim analyses do not include the RSWAN steps as per now).
#' @param NI_years a vector of years for which to calculate NI indicator values.
#' @param refAnchorYear integer. Which year in the current NI indicator data to 
#' use as the base for re-calibrating the references state. 
#' @param refAnchorYear_alt integer. Alternative year to use as a base for re-
#' calibrating reference state. Gets used if `refAnchorYear` is not present in 
#' Trim data. 
#' @param nsim integer. Number of simulations to run for calculating indicator
#' values as 3-year averages of Trim indices.  
#'
#' @returns a list containing updated NI indicator data for each species in `sppNI`.
#' Ready formatted for upload to NI database. 
#' @export
#'
#' @examples

prepareIndicatorData_NI <- function(sppNI, OldIndicator_data, updatePlans,
                                    inputFile_folders, use_combTS,
                                    NI_years, refAnchorYear, refAnchorYear_alt,
                                    nsim){
  
  ## Set up list for storing data
  UpdatedIndicator_data <- list()
  
  
  for(i in 1:nrow(sppNI)){
    
    message("")
    message(crayon::bold(paste0("Preparing data for indicator: ", sppNI$Species[i], " (", sppNI$indicatorName[i], ")")))
    
    ## Check areas used in NI database
    NI_areas <- updatePlans[[i]]$NI_areas
    message("NI database contains the following ", length(NI_areas), " area(s):")
    print(NI_areas)
    
    ## Check sub-areas of TRIM data to be used
    Trim_subAreas <- unique(updatePlans[[i]]$Trim_areaMatch)
    message(paste0("Trim analyses for the following (sub-)area(s) will be used for updating:"))
    print(Trim_subAreas)

    ## Check for availability of (combined) TRIM index data (only for indicators calculated for all of Norway)
    if(use_combTS & sppNI$StartDataHFT[i] == 1996){
      source_folder_COMB <- subset(inputFile_folders, ID == "PECBMS")$path
      spp_files <- list.files(source_folder_COMB)[which(grepl(paste0(sppNI$EURINGCode[i], "_"), list.files(source_folder_COMB)))]
      combTS <- ifelse(any(grepl(paste0("COMB_", sppNI$EURINGCode[i], "_"), spp_files)), TRUE, FALSE)
    }else{
      combTS <- FALSE
    }
    
    if(combTS){
      message("Using combined Trim time series (start 1996).")
    }else{
      message("Using TOV-E Trim time series (start 2008).")
    }
    
    
    ## Set up dataframe for temporarily storing area-specific indicator data
    Indicator_upload <- data.frame()
    
    for(x in 1:nrow(updatePlans[[i]])){
      
      ## Extract active area names
      NI_area_x <- as.character(updatePlans[[i]][x, "NI_areas"])
      Trim_area_x <- as.character(updatePlans[[i]][x, "Trim_areaMatch"])
      
      
      ## List files from relevant folder
      source_folder <- dplyr::case_when(
        combTS ~ subset(inputFile_folders, ID == "PECBMS")$path,
        Trim_area_x == "Norge" ~ subset(inputFile_folders, ID == "NI_Norge")$path,
        Trim_area_x == "Sør-Norge" & nrow(updatePlans[[i]]) == 1 ~ subset(inputFile_folders, ID == "NI_onlySNorge")$path,
        Trim_area_x == "Sør-Norge" & nrow(updatePlans[[i]]) > 1 ~ subset(inputFile_folders, ID == "NI_SNorge")$path,
        Trim_area_x == "Nord-Norge" ~ subset(inputFile_folders, ID == "NI_NNorge")$path
      )
      
      if(is.na(source_folder)){
        message("Source file folder could not be determined. Please check.")
        message("Jumping to next species-area")
        next()
      }
      
      spp_files <- list.files(source_folder)[which(grepl(paste0(sppNI$EURINGCode[i], "_"), list.files(source_folder)))]
      
      
      ## Load relevant TRIM index data 
      message(paste0("Loading TRIM index data for area: ", Trim_area_x))
      if(combTS){
        index_TS_file <- spp_files[which(grepl(paste0("COMB_", sppNI$EURINGCode[i], "_"), spp_files) & grepl("indices_TT", spp_files))]
      }else{
        index_TS_file <- spp_files[which(startsWith(spp_files, paste0(sppNI$EURINGCode[i], "_")) & grepl("indices_TT", spp_files))]
      }
      
      if(length(index_TS_file) > 1){
        stop("Identification of correct TRIM index data file failed. Check presence of file and detection criteria.")
      }
      
      TrimIndex_data <- read.csv(paste0(source_folder, "/", index_TS_file), sep = ";")
      
      ## Calculate 3-year averages for TRIM index data
      message("Calculating 3-year averages for TRIM data...")
      TrimIndex_averages <- data.frame()
      
      for(t in 3:nrow(TrimIndex_data)){
        idx <- matrix(NA, nrow = nsim, ncol = 3)
        
        for(s in 1:3){
          idx[, s] <- truncnorm::rtruncnorm(n = nsim, 
                                            mean = TrimIndex_data$Index_model[t-(s-1)],
                                            sd = TrimIndex_data$Index_model_SE[t-(s-1)],
                                            a = 0)
        }
        
        avg_3yr <- rowMeans(idx)
        
        avg_data <- data.frame(Year = TrimIndex_data$Year[t],
                               mean = mean(avg_3yr),
                               sd = sd(avg_3yr),
                               lowerQuart = unname(quantile(avg_3yr, probs = 0.25)),
                               upperQuart = unname(quantile(avg_3yr, probs = 0.75)))
        
        TrimIndex_averages <- rbind(TrimIndex_averages, avg_data)
      }
      
      ## Double-check that reference anchor year is in data (and use alternative if not)
      if(refAnchorYear %in% TrimIndex_averages$Year){
        refAnchorYear_use <- refAnchorYear
      }else{
        refAnchorYear_use <- refAnchorYear_alt
        message(crayon::cyan("NOTE: Switched to alternative reference anchor year (selected anchor year not in data)."))
      }
      
      ## Extract reference proportion for the reference anchor year in the relevant area from old NI data
      refProp_orig <- subset(OldIndicator_data[[i]]$indicatorValues, 
                             yearName == refAnchorYear_use & areaName == NI_area_x)$verdi/100
      
      ## Calculate and add reference value
      message(paste0("Recalibrate reference value for NI database area ", NI_area_x, "..."))
      
      ref <- subset(TrimIndex_averages, Year == refAnchorYear_use) %>%
        dplyr::mutate(mean = mean/refProp_orig,
                      sd = sd/refProp_orig,
                      lowerQuart = lowerQuart/refProp_orig,
                      upperQuart = upperQuart/refProp_orig)
      ref[,"Year"] <- "Referanseverdi"
      
      TrimIndex_averages <- rbind(TrimIndex_averages, ref)
      
      ## Do pre-scaling and drop values not relevant for upload to NI database
      message(paste0("Scale and format updated indicator data for NI database area ", NI_area_x, "..."))
      Indicator_prescaled <- TrimIndex_averages %>%
        dplyr::mutate(mean = 100*mean/ref$mean,
                      sd = 100*sd/ref$mean,
                      lowerQuart = 100*lowerQuart/ref$mean,
                      upperQuart = 100*upperQuart/ref$mean) %>%
        dplyr::filter(Year %in% c("Referanseverdi", NI_years)) %>%
        dplyr::rename(yearName = Year)

      row.names(Indicator_prescaled) <- NULL
      
      # Write updated indicator data in upload format
      Indicator_upload_new <- OldIndicator_data[[i]]$indicatorValues %>%
        dplyr::filter(areaName == NI_area_x) %>%
        dplyr::left_join(Indicator_prescaled, by = "yearName") %>%
        dplyr::mutate(update = ifelse(is.na(mean), FALSE, TRUE)) %>%
        dplyr::mutate(
          
          # (Over)write indicator values
          verdi = ifelse(update, mean, verdi),
          nedre_Kvartil = ifelse(update, lowerQuart, nedre_Kvartil),
          ovre_Kvartil = ifelse(update, upperQuart, ovre_Kvartil),
          
          # Set data type for updated values
          datatypeId = ifelse(update & yearName != "Referanseverdi", 3, datatypeId),
          datatypeName = ifelse(update & yearName != "Referanseverdi", "Beregnet fra modeller", datatypeName),
          
          # Remove no longer valid information on distributions
          distributionName = ifelse(update, NA, distributionName),
          distributionId = ifelse(update, NA, distributionId),
          distParam1 = ifelse(update, NA, distParam1),
          distParam2 = ifelse(update, NA, distParam2)) %>%
        
        # Delete auxiliary columns
        dplyr::select(-mean, -sd, -lowerQuart, -upperQuart, -update)
      
      ## Append
      Indicator_upload <- rbind(Indicator_upload, Indicator_upload_new)
    }
    
    UpdatedIndicator_data[[i]] <- OldIndicator_data[[i]]
    UpdatedIndicator_data[[i]]$indicatorValues <- Indicator_upload
    
    # Check for any legacy custom distribution information
    if(length(UpdatedIndicator_data[[i]]$customDistributions) > 0){
      warning(paste0("Legacy custom distribution information is present for species ", sppNI$Species_Norwegian_NI[i], ". Please review and delete information."))
    }
  }
  
  # Return data
  return(UpdatedIndicator_data)
}
