#' Run PECBMS' Rtrim shell script
#'
#' @param folder character. Path to directory which contains the RTRIM argument
#' and output files, and into which to save outputs.
#'
#' @return
#' @export
#'
#' @examples

runRtrimShell_PECBMS <- function(folder){
  
  ## Check whether required helper functions are available
  rtrim_shell_funNames <- c("fill_All_Indices_All_Trends",
                            "fill_arg_output_file",
                            "fill_Indices_TT_file",
                            "make_All_Indices_All_Trends",
                            "make_arg_output_file",
                            "make_Indices_TT_file",
                            "makeOverview")
  
  rtrim_shell_funExist <- rep(NA, length(rtrim_shell_funNames))
  
  for(i in 1:length(rtrim_shell_funNames)){
    rtrim_shell_funExist[i] <- exists(rtrim_shell_funNames[i])
  }
  
  if(any(!rtrim_shell_funExist)){
    missingFun <- rtrim_shell_funNames[which(!rtrim_shell_funExist)]
    stop(paste0("The following Rtrim shell helper functions are missing: ", 
                paste(missingFun, collapse = ", "),
                ". \nThese functions are provided in the script x3_functions_Shell_Rtrim.R"))
  }
  
  # ----------------------------------------------------------------------------
  # START OF CODE FROM PECBMS
  #-----------------------------------------------------------------------------
  
  #####################################################################################################################
  # RTRIM-shell. Version RTRIM-shell_v1.7
  # Tool to run rtrim (= TRIM in R) for multiple species and subsets of sites.
  # Marnix de Zeeuw rtrim@cbs.nl Statistics Netherlands 2019. 
  #####################################################################################################################
  
  #####################################################################################################################
  # EXPLANATION OF THE FILES USED AS INPUT. 
  #####################################################################################################################
  
  # The code in this shell reads .csv tables: _counts.csv & _arg_input_stratum.csv. These tables contain counts for a particular combination of species and stratum (stratum is e.g. a region).
  ### _counts.csv ### 
  # These counts are analysed using rtrim. 
  # The most simple format of _counts.csv (file containing the counts) is:
  #
  # site year count
  #    3 1994    NA
  #    3 1995     1
  #    3 1996     2
  #
  # Use "NA" for missing values instead of -1 as in the old trim version TRIM/BirdStats (in delphi).
  #
  # If weights are available, the format is:
  #
  # site year count weights
  #    3 1994    NA       1
  #    3 1995     1       1
  #    3 1996     2       1
  #
  # When no weights are used, then don't include that column.
  #
  # File name for proper data processing must be a combination of letters (e.g. BMP), 
  # followed by an underscore and a number that represents the species number,  
  # followed by an underscore and a number that represents the stratum type (1=standard; 2=strata combined, processing the latter however requires a different tool),  
  # followed by an underscore and a number that represents the stratum number, 
  # followed by "_counts.csv".
  # e.g. BMP_11870_1_4_counts.csv
  # In the example, BMP stands for Bird Monitoring Program, 11870 is the species code of blackbird, 1 refers to stratumtype and 4 is a country code
  #
  ### _arg_input_stratum.csv###
  # In addition to a file with count data, rtrim requires arguments to inform how the analysis should be done. 
  # For each file with count data, an associated argument file is needed, with a similar name, except for the last part: "_arg_input_stratum.csv" instead of "_counts.csv" .
  # Thus, datafile BMP_11870_1_4_counts.csv needs to be accompanied by arguments file BMP_11870_1_4_arg_input_stratum.csv.
  
  # Format of _arg_input_stratum.csv is as follows:
  #
  # File        Base_year_first_year Base_year_last_year Changepoints Serial_correlation Overdispersion Presence_weights Presence_monthfactors Year_from Save_fitted_values  
  # BMP_11870_1_4        1990                 1990          ALL          TRUE               TRUE               TRUE             FALSE	             2007         TRUE               
  #
  # File: short name of file to be analysed, i.e., datafilename without _"counts.csv"
  
  # Rather than choosing one particular base year (as in the old TRIM), it is also possible to use a number of years as base period.
  # This is NOT recommended for the PECBMS data, use first year as Base_year_first_year and Base_year_last_year, please. # 28-02-2022 Eva Silarova
  # Base_year_first_year: Calendar year used as base year for indices. 
  # Base_year_last_year: Calendar year used as base year for indices.                                                                      
  # If Base_year_first_year equals Base_year_last_year a single year is used as base year
  # If Base_year_first_year < Base_year_last_year, a period is used as base time
  # In the latter case, Base_year_first_year is the first year of the period.
  
  # Changepoints: Changepoints used to calculate the indices are usually set to "all". !!!BE AWARE!!! If any of changepoints shall be excluded e.g. due to covid, 
  #    changepoints shall be listed in _arg_input_stratum, do not include last year of time series!!!
  # Serial correlation: Set to "TRUE" when serial correlation needs to be taken into account. Otherwise set to "FALSE".
  # Overdispersion: Set to "TRUE" when overdisperion needs to be taken into account. Otherwise set to "FALSE".
  # Presence_weights: TRUE when weights are available, FALSE when not.
  # Presence_monthfactors: TRUE when available, FALSE when not. NOT recommended for PECBMS.
  # Year_from: first calendar year of the subperiod over which also a slope needs to be computed.
  #    Last year of the subperiod equals the last year of the entire period.
  # Save_fitted_values: Set to "TRUE" when the fitted values need to be saved. Otherwise set to "FALSE".
  
  # rtrim starts with trying the most elaborate model.
  # When the analysis of that model fails, a less elaborate model is tried.
  # Output for each species is stored in a RData-file, in the working directory with the count files.
  # That output needs to be processed with a separate script ("processingOutputSpeciesStratumCombinations.r").
  
  # rtrim enables to run models with covariates, but this rtrim shell does NOT support that option.  
  # For more information we refer to the trim-vignettes in R.
  
  # !!! IMPORTANT!!! CHANGES YOU have to do in the script:
  # 1. Line 94: Adjust the path to the file called functions_Shell_Rtrim.R according to your computer (try different slashes, if double backslash does not work)
  # 2. Line 101: Adjust the path to your working folder according to your computer (try different slashes, if double backslash does not work)
  
  ############################################################################
  # Changes done in the script for the version 6 (Alena Jechumtal Skalova & Martin Stjernman 30/05/2022):
  # Lines 151-152 were modified to prevent mismatches between the species-codes and their first-last year whenever species does not have similar first_year.
  #       Moreover, this modifications makes the program to analyze the first year of detection not only the first year of scheme.
  
  #####################################################################################################################
  # PREPARATION.
  #####################################################################################################################
  # 1/06/2020 Eva Silarova unified upper and lower cases in variable names - v1.1

  # Select all files in the directory which contain arguments.
  # The files can be recognized by the pattern "arg_input_arguments".
  listSpeciesStratumCombinations <- dir(folder, pattern = "arg_input_stratum.csv")
  
  # The number of files in the directory.                                                
  numberSpeciesStratumCombinations <- length(listSpeciesStratumCombinations)
  
  # Make a table for an overview. 
  # The table lists which analyses were successful or not for all runs.
  
  overview <- makeOverview(listSpeciesStratumCombinations)
  # The code of this function can be found in functions_Shell_Rtrim.R.
  # The source command in the main script makes it possible to use this function in this code.
  
  #####################################################################################################################
  # READING ARGUMENT AND COUNT FILES. 
  #####################################################################################################################
  
  for (j in 1:numberSpeciesStratumCombinations) {
    
    cat("Processing", listSpeciesStratumCombinations[j], "\n") # 19/07/2022 John Kennedy inserted instrumentation to find out what species/attempt each warning was for - v1.6
    

    # The file with arguments contains the information to analyse the counts for a particular combination of species and stratum. 
    # The arguments are used when calling the function 'rtrim'
    arguments <- read.csv2(paste0(folder, "/", listSpeciesStratumCombinations[j]), header = TRUE, stringsAsFactors = FALSE)
    
    arguments<-as.list(arguments)                                                  # Enables the script to read the change points properly, 
    if(!arguments$Changepoints%in%c("all","auto")) {                               # when they are specified in the ?_arg_input_stratum.csv? file as comma-separated numbers.
      Changepoints<-as.integer(unlist(strsplit(arguments$Changepoints,"-")))       # Function was created 20/03/2022 by Dario Massimino -> modified 11/04/2022 by Meelis Leivits -> modified 20/07/2022 by John Kennedy - v1.6
      arguments$Changepoints<-Changepoints                                         # Function was adjusted to use of "-" in the changepoints # 9/02/2023 Javier Rivas Salvador - v1.6
    } 
    
    # The file with counts contains the counts for a particular combination of species and stratum.
    # Weights may also be present in this file.
    #13/01/2023 Javier Rivas: I have modified the following lines to ensure that independently of the decimal symbol the counts$count variable is loaded as numeric
    counts <- read.table(paste0(folder, "/", arguments$File, "_counts.csv"), stringsAsFactors = FALSE, header = TRUE, dec = ".",sep=";")   
    
    if(is.character(counts$count)==T){
      counts <- read.table(paste0(folder, "/", arguments$File, "_counts.csv"), stringsAsFactors = FALSE, header = TRUE, dec = ",",sep=";")   
    }
    
    # 13/01/2023 end of modifications 
    # 30/05/22 Alena Jechumtal Skalova & Martin Stjernman corrected mismatches in the species-codes and their first last year.-v1.6 
    #   Moreover, the code evaluates the first year of detection and not the first year of the scheme.-v1.6
    overview$first_year[overview$ss_combinations == arguments$File] <- min(counts$year[counts$count > 0], na.rm = TRUE) 
    overview$last_year[overview$ss_combinations == arguments$File]  <- max(counts$year[counts$count > 0], na.rm = TRUE)
    # Original code
    # overview$first_year[j] <- min(counts$year, na.rm = TRUE)
    # overview$last_year[j]  <- max(counts$year, na.rm = TRUE)
    # End of the editation 30/05/2022 
    
    # First and last year for each combination of stratum and species is stored for later use.
    
    #####################################################################################################################
    # RUNNING RTRIM.
    #####################################################################################################################
    
    # Start with the most elaborate model and switch automatically to a more simple model when needed.
    
    result <- tryCatch(
      {
        # Attempt 1
        cat(" Attempt 1\n") #	19/07/2022 John Kennedy inserted instrumentation to find out what species/attempt each warning was for. - v1.6
        if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == TRUE) {
          trim(count ~ site + (year + month), data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
        } else {
          if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
            trim(count ~ site + year, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
          } else{
            if(arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == TRUE) {
              trim(count ~ site + (year + month), data = counts,            model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
            } else{
              trim(count ~ site + year, data = counts,                      model = 2, changepoints = arguments$Changepoints, serialcor = arguments$Serial_correlation, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
            }
          }
        }
      }
      , error = warning)     
    
    if (class(result) == "trim") {
      
      save(x = result,  file = paste0(folder, "/", arguments$File, ".RData"))
      overview$attempt_1[overview$ss_combinations == arguments$File] <- "success"
      overview$success[overview$ss_combinations == arguments$File] <- "yes"
      
    } else {
      
      # First attempt failed, try a less elaborate model by setting serial correlation off. 
      # Also, when month factors are available, no changepoints are estimated in the next model.
      
      overview$attempt_1[overview$ss_combinations == arguments$File] <- "error"
      overview$error_1[overview$ss_combinations == arguments$File] <- result
      
      result <- tryCatch(
        {
          # attempt 2 
          cat(" Attempt 2\n")          # 19/07/2022 John Kennedy inserted instrumentation to find out what species/attempt each warning was for. - v1.6
          if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == TRUE) {
            trim(count ~ site + (year + month), data = counts, weights = "weights", model = 2,                               serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
          } else {
            if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
              trim(count ~ site + year, data = counts, weights = "weights", model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
            } else{
              if(arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == TRUE) {
                trim(count ~ site + (year + month), data = counts,                      model = 2,                               serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
              } else{
                trim(count ~ site + year, data = counts,                      model = 2, changepoints = arguments$Changepoints, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
              }
            }
          }
        }
        , error = warning)
      
      if (class(result) == "trim") {
        
        save(x = result,  file = paste0(folder, "/", arguments$File, ".RData"))
        overview$attempt_2[overview$ss_combinations == arguments$File] <- "success"
        overview$success[overview$ss_combinations == arguments$File] <- "yes"
        
      } else {
        cat(" Attempt 3\n")          # 19/07/2022 John Kennedy inserted instrumentation to find out what species/attempt each warning was for. - v1.6
        
        # Second attempt also failed. Now try an even more simple model: no month factors, no changepoints, but serial correlation switched on. 
        # Note that no further options are available to include month factors in the model. 
        
        overview$attempt_2[overview$ss_combinations == arguments$File] <- "error"
        overview$error_2[overview$ss_combinations == arguments$File] <- result
        
        if (arguments$Presence_monthfactors == TRUE) {
          
          cat("Analysis failed for this combination of species and stratum:", arguments$File, "\n")
          
        }
        
        result <- tryCatch( 
          {
            # attempt 3
            
            if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
              
              trim(count ~ site + year, data = counts, weights = "weights", model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
              
            } else{
              
              if (arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == FALSE) {
                trim(count ~ site + year, data = counts,                      model = 2, serialcor = TRUE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
              }
            }
          }
          , error = warning)
        
        if (class(result) == "trim") {
          
          save(x = result,  file = paste0(folder, "/", arguments$File, ".RData"))
          overview$attempt_3[overview$ss_combinations == arguments$File] <- "success"
          overview$success[overview$ss_combinations == arguments$File] <- "yes"
          
        } else {
          if (arguments$Presence_monthfactors == FALSE) {
            overview$attempt_3[overview$ss_combinations == arguments$File] <- "error"
            overview$error_3[overview$ss_combinations == arguments$File] <- result
          }
          
          # Third attempt also failed. Now try the most simple model: no changepoints at all and no serial correlation. 
          
          result <- tryCatch(
            {
              # Final attempt
              cat(" Attempt 4\n")          # 19/07/2022 John Kennedy inserted instrumentation to find out what species/attempt each warning was for. - v1.6
              
              if (arguments$Presence_weights == TRUE & arguments$Presence_monthfactors == FALSE) {
                trim(count ~ site + year, data = counts, weights = "weights", model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)
              } else{
                if (arguments$Presence_weights == FALSE & arguments$Presence_monthfactors == FALSE) {
                  trim(count ~ site + year, data = counts,                      model = 2, serialcor = FALSE, overdisp = arguments$Overdispersion, max_iter = 200, conv_crit = 1e-5)    
                } 
              }
              
            }
            
            , error = warning)
          
          if (class(result) == "trim") {
            
            save(x = result,  file = paste0(folder, "/", arguments$File, ".RData"))
            overview$attempt_4[overview$ss_combinations == arguments$File] <- "success"
            overview$success[overview$ss_combinations == arguments$File] <- "yes"
            
          } else {
            if (arguments$Presence_monthfactors == FALSE) {
              
              # When this analysis also fails, send a error message to screen.
              
              overview$attempt_4[overview$ss_combinations == arguments$File] <- "error"
              overview$error_4[overview$ss_combinations == arguments$File] <- result
              
              cat("Analysis failed for this combination of species and stratum:", arguments$File, "\n")
            }        
          }
        }  
      } 
    }  
  }
  
  #####################################################################################################################
  # WRITING OVERVIEW OF RTRIM SUCCESSES and FAILURES.  
  #####################################################################################################################
  overview <- overview[order(overview$species_number, overview$stratum_number), ]
  write.table(overview, paste0(folder, "/overview.csv"), row.names = FALSE, sep = ";", dec = ".")   #13/01/2023 Javier Rivas: Changes to produce outputs with dots
  
  
}
