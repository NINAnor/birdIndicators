prepareIndicatorData_NI <- function(sppNI, inputFile_folders, use_combTS){
  
  ## Set up lists for storing data
  TrimIndex_data <- list()
  OldIndicator_data <- list()
  UpdatedIndicator_data <- list()
  
  
  for(i in 1:nrow(sppNI)){
    
    message("")
    message(paste0("Preparing data for indicator: ", sppNI$Species[i], " (", sppNI$indicatorName[i], ")"))
    
    ## Download old indicator data from NI database
    message("Retrieving data from NI database...")
    OldIndicator_data[[i]] <- NIcalc::getIndicatorValues(indicatorID = sppNI$indicatorId[i])
    
    ## Check areas used in NI database
    areas <- unique(OldIndicator_data[[i]]$indicatorValues$areaName)
    message("NI database contains the following ", length(areas), " area(s):")
    print(areas)
    
    ## Check whether TRIM data have been divided into sub-areas
    subAreas <- ifelse(sppNI$dataUse_direct_N[i] | sppNI$dataUse_expert_N[i], FALSE, TRUE)
    
    #if(length(areas) > 1 & subAreas){
    #  message("Distinct indicator data has to be reported for mulitple areas. Functionality for this is not implemented yet, but will be added soon. For now, this indicator will be skipped.")
    #  next()
    #}
    
    ## Check for availability of (combined) TRIM index data (only for indicators calculated for all of Norway)
    if(use_combTS & sppNI$StartDataHFT[i] == 1996){
      source_folder_COMB <- subset(inputFile_folders, ID == "PECBMS")$path
      spp_files <- list.files(source_folder_COMB)[which(grepl(paste0(sppNI$EURINGCode[i], "_"), list.files(source_folder_COMB)))]
      combTS <- ifelse(any(grepl(paste0("COMB_", sppNI$EURINGCode[i], "_"), spp_files)), TRUE, FALSE)
    }else{
      combTS <- FALSE
    }
    
    ## Count number of areas with separate TRIM data
    if(sppNI$dataUse_direct_SN[i] & sppNI$dataUse_direct_NN[i] & length(areas) > 1){
      areas_Trim <- c("Sør-Norge", "Nord-Norge")
    }else{
      areas_Trim <- ifelse(sppNI$dataUse_direct_SN[i] & !sppNI$dataUse_direct_NN[i], "Sør-Norge", "Norge")
    }
    
    nAreas_Trim <- length(areas_Trim)
    
    ## Set up dataframe for temporarily storing area-specific indicator data
    Indicator_upload <- data.frame()
    
    
    for(x in 1:nAreas_Trim){
      
      ## List files from relevant folder
      source_folder <- dplyr::case_when(
        combTS ~ subset(inputFile_folders, ID == "PECBMS")$path,
        sppNI$dataUse_expert_N[i] ~ subset(inputFile_folders, ID == "NI_Norge")$path,
        sppNI$dataUse_direct_N[i] ~ subset(inputFile_folders, ID == "NI_Norge")$path,
        sppNI$dataUse_direct_SN[i] & !sppNI$dataUse_direct_NN[i] ~ subset(inputFile_folders, ID == "NI_onlySNorge")$path,
        areas_Trim[x] == "Sør-Norge" ~ subset(inputFile_folders, ID == "NI_SNorge")$path,
        areas_Trim[x] == "Nord-Norge" ~ subset(inputFile_folders, ID == "NI_NNorge")$path
      )
      
      if(is.na(source_folder)){
        message("Source file folder could not be determined. Please check.")
        message("Jumping to next species-area")
        next()
      }
      
      spp_files <- list.files(source_folder)[which(grepl(paste0(sppNI$EURINGCode[i], "_"), list.files(source_folder)))]
      
      
      ## Load relevant TRIM index data 
      message(paste0("Loading TRIM index data for area: ", areas_Trim[x]))
      if(combTS){
        index_TS_file <- spp_files[which(grepl(paste0("COMB_", sppNI$EURINGCode[i], "_"), spp_files) & grepl("indices_TT", spp_files))]
      }else{
        index_TS_file <- spp_files[which(startsWith(spp_files, paste0(sppNI$EURINGCode[i], "_")) & grepl("indices_TT", spp_files))]
      }
      
      if(length(index_TS_file) > 1){
        stop("Identification of correct TRIM index data file failed. Check presence of file and detection criteria.")
      }
      
      TrimIndex_data[[i]] <- read.csv(paste0(source_folder, "/", index_TS_file), sep = ";")
      
      ## Calculate 3-year averages for TRIM index data
      message("Calculating 3-year averages for TRIM data...")
      TrimIndex_averages <- data.frame()
      
      for(t in 3:nrow(TrimIndex_data[[i]])){
        idx <- matrix(NA, nrow = nsim, ncol = 3)
        
        for(s in 1:3){
          idx[, s] <- truncnorm::rtruncnorm(n = nsim, 
                                            mean = TrimIndex_data[[i]]$Index_model[t-(s-1)],
                                            sd = TrimIndex_data[[i]]$Index_model_SE[t-(s-1)],
                                            a = 0)
        }
        
        avg_3yr <- rowMeans(idx)
        
        avg_data <- data.frame(Year = TrimIndex_data[[i]]$Year[t],
                               mean = mean(avg_3yr),
                               sd = sd(avg_3yr),
                               lowerQuart = unname(quantile(avg_3yr, probs = 0.25)),
                               upperQuart = unname(quantile(avg_3yr, probs = 0.75)))
        
        TrimIndex_averages <- rbind(TrimIndex_averages, avg_data)
      }
      
      Indicator_prescaled_areas <- list()
      
      ## Determine whether we need to use the same Trim data for setting values for multiple areas
      areas_count <- ifelse(nAreas_Trim == 1 & length(areas) > 1, length(areas), 1)
      
      for(a in 1:areas_count){
        
        # Add area information to TRIM index data
        TrimIndex_averages_add <- TrimIndex_averages %>%
          dplyr::mutate(areaName = areas[a])
        
        # Extract reference proportion for the reference anchor year in the relevant area from old NI data
        refProp_orig <- subset(OldIndicator_data[[i]]$indicatorValues, 
                               yearName == refAnchorYear & areaName == areas[a])$verdi/100
        
        # Calculate and add reference value
        message(paste0("Recalibrate reference value for NI database area ", areas[a], "..."))
        ref <- subset(TrimIndex_averages_add, Year == refAnchorYear) %>%
          dplyr::mutate(mean = mean/refProp_orig,
                        sd = sd/refProp_orig,
                        lowerQuart = lowerQuart/refProp_orig,
                        upperQuart = upperQuart/refProp_orig, 
                        areaName = areas[a])
        ref[,"Year"] <- "Referanseverdi"
        
        TrimIndex_averages_add <- rbind(TrimIndex_averages_add, ref)
        
        # Do pre-scaling and drop values not relevant for upload to NI database
        message(paste0("Scale and format updated indicator data for NI database area ", areas[a], "..."))
        Indicator_prescaled_areas[[a]] <- TrimIndex_averages_add %>%
          dplyr::mutate(mean = 100*mean/ref$mean,
                        sd = 100*sd/ref$mean,
                        lowerQuart = 100*lowerQuart/ref$mean,
                        upperQuart = 100*upperQuart/ref$mean) %>%
          dplyr::filter(Year %in% c("Referanseverdi", NI_years)) %>%
          dplyr::rename(yearName = Year)
      }
      
      # Write list as a dataframe
      Indicator_prescaled <- dplyr::bind_rows(Indicator_prescaled_areas)
      row.names(Indicator_prescaled) <- NULL
      
      # Write updated indicator data in upload format
      Indicator_upload_new <- OldIndicator_data[[i]]$indicatorValues %>%
        dplyr::left_join(Indicator_prescaled, by =c("yearName", "areaName")) %>%
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
  
  # Collate and return data
  results <- list(TrimIndex_data = TrimIndex_data,
                  OldIndicator_data = OldIndicator_data,
                  UpdatedIndicator_data = UpdatedIndicator_data)
  
  return(results)
}
