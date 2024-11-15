#####################################################################################################################
# RTRIM-shell functions. Version RTRIM-shell_2.0
# Functions to support the  rtrim shell scripts.
# Marnix de Zeeuw rtrim@cbs.nl Statistics Netherlands 2019. 

# Last: modifications J.Rivas 13/03/2024
#####################################################################################################################
# 1/06/2020 Eva Silarova unified upper and lower cases in variable names - v1.1
make_Indices_TT_file <- function(result) {
  # Aim: this function creates an empty dataframe for indices and time totals for a species-stratum combination.
  
  indices_TT_file <- data.frame(Year = integer(result$nyear), 
                                # Years for which indices and time totals have been calculated.
                                
                                TT_model = numeric(result$nyear), 
                                # Time totals based on the estimated model.
                                
                                TT_model_SE = numeric (result$nyear), 
                                # Standard errors of model time totals.
                                
                                TT_imputed = numeric(result$nyear), 
                                # Imputed time totals.
                                
                                TT_imputed_SE = numeric(result$nyear), 
                                # Standard errors of imputed time totals.
                                
                                Index_model = numeric(result$nyear), 
                                # Indices based on the estimated model.
                                
                                Index_model_SE = numeric(result$nyear), 
                                # Standard errors of model indices.
                                
                                Index_imputed = numeric(result$nyear), 
                                # Indices based on the imputed time totals.
                                
                                Index_imputed_SE = numeric(result$nyear))
  # Standard errors of imputed indices.
  
  return(indices_TT_file)
}

#####################################################################################################################
make_arg_output_file <- function(arguments) {
  # Aim: this function creates an empty dataframe for slopes spanning the entire study period and a subperiod for a species-stratum combination.
  
  arg_output_file <- data.frame(N_sites = integer(1), 
                                # Number of unique sites.
                                
                                N_time_values = integer(1), 
                                # Number of unique years.
                                
                                N_observed_zero_counts = integer(1), 
                                # Number of zero counts.
                                
                                N_observed_positive_counts = integer(1), 
                                # Number of positive counts.
                                
                                N_missing_counts = integer(1), 
                                # Number of missing counts.
                                
                                N_counts_total = integer(1), 
                                # Total number of counts.
                                
                                Base_year_first_year = numeric(1),
                                Base_year_last_year = numeric(1),
                                
                                # Calendar year used as base year for indices. 
                                # If Base_year_first_year equals Base_year_last_year a single year is used as base year.
                                # If Base_year_first_year < Base_year_last_year, a period is used as base time.
                                # In the latter case, Base_year_first_year is the first year of the period.
                                
                                Changepoints = numeric(1), 
                                # Changepoints used.
                                
                                Overdispersion = numeric(1),
                                # Estimated overdispersion.
                                
                                Serial_correlation = numeric(1),
                                # Estimated serial correlation.
                                
                                Slope_imputed_mul = numeric(1),
                                # Multiplicative imputed slope for the entire period.
                                
                                Slope_imputed_mul_SE = numeric(1),
                                # Standard error of multiplicative imputed slope for the entire period.
                                
                                Slope_imputed_classification = character(1),
                                # Trend classification of multiplicative imputed slope for the entire period.
                                
                                Year_from = numeric(1),
                                # First year of the subperiod from which a slope has been calculated. 
                                # Last year of the subperiod is equal to the last year of the entire period.
                                # Note that it is assumed that the first year of the subperiod is closer to the present than the first year of the entire period.
                                
                                Slope_from_imputed_mul = numeric(1),
                                # Multiplicative imputed slope for the subperiod.
                                
                                Slope_from_imputed_mul_SE = numeric(1),
                                # Standard error of multiplicative imputed slope for the subperiod.
                                
                                Slope_from_imputed_classification = character(1)
                                # Trend classification of multiplicative imputed slope for the subperiod.
                                
  )
  
  arg_output_file$Year_from <- arguments$Year_from  
  #   Here, a value is already assigned to Year_from, to make the function that fills the arg_output file more simple. 
  #   The variabele 'arguments' is not needed due to this, when calling the function.
  
  arg_output_file$Base_year_first_year <- arguments$Base_year_first_year
  arg_output_file$Base_year_last_year  <- arguments$Base_year_last_year
  
  return(arg_output_file)
}

#####################################################################################################################
fill_Indices_TT_file <- function(indices_TT_file, result, arguments) {
  # Aim: this function fills the dataframe with indices and time totals for a species-stratum combination.    
  
  indices <- index(result, which = "both", covars = FALSE, base = arguments$Base_year_first_year:arguments$Base_year_last_year)
  time_totals <- totals(result, which = "both")
  
  indices_TT_file$Year <- indices$time
  # The years over which the indices and time totals have been calculated.
  
  indices_TT_file$TT_model <- time_totals$fitted 
  # Time totals based on the estimated model. 
  
  indices_TT_file$TT_model_SE <- time_totals$se_fit
  # Standard errors of model time totals.
  
  indices_TT_file$TT_imputed <- time_totals$imputed
  # Imputed time totals.
  
  indices_TT_file$TT_imputed_SE <- time_totals$se_imp
  # Standard errors of imputed time totals.
  
  indices_TT_file$Index_model <- indices$fitted
  # Indices based on the estimated model.
  
  indices_TT_file$Index_model_SE <- indices$se_fit
  # Standard errors of model indices.
  
  indices_TT_file$Index_imputed <- indices$imputed
  # Indices based on the imputed time totals.
  
  indices_TT_file$Index_imputed_SE <- indices$se_imp
  # Standard errors of imputed indices.
  
  return(indices_TT_file)
}

#####################################################################################################################
fill_arg_output_file <- function(arg_output_file, result, counts) {
  # Aim: this function fills the dataframe with additional information.  
  
  overviewCounts <- count_summary(counts, count_col = "count", site_col = "site", year_col = "time") 
  # count_summary is a RTRIM function.
  
  slopes_imputed <- overall(result, which = "imputed")
  slopes_fitted  <- overall(result, which = "fitted")
  # The slopes over the entire period. Both variables are lists, a specific class of R-objects.
  
  slopes_subperiod_imputed <- overall(result, which = "imputed", changepoints = c(arg_output_file$Year_from)) 
  # Year_from is already available in the arg_output_file.
  slopes_subperiod_fitted <- overall(result, which = "fitted", changepoints = c(arg_output_file$Year_from))
  # The slopes over a subperiod. Both variables are lists.
  
  arg_output_file$N_sites <- result$nsite
  # Number of unique sites in the counts file used. 
  
  arg_output_file$N_time_values <- result$ntime
  # Number of unique years in the counts file used.
  
  arg_output_file$N_observed_zero_counts <- overviewCounts$zero_counts
  # Number of zero counts in the entire period, in the counts file used.
  
  arg_output_file$N_observed_positive_counts <- overviewCounts$positive_counts
  # Number of positive counts in the entire period, in the counts file used.
  
  arg_output_file$N_missing_counts <- overviewCounts$missing_counts
  # Number of missing counts in the entire period, in the counts file used.
  
  arg_output_file$N_counts_total <- overviewCounts$total_counts
  # Total number of counts in the entire period, in the counts file used.
  
  arg_output_file$Changepoints <- paste(result$changepoints, collapse = ", ")
  # Changepoints used to calculate the indices.
  
  arg_output_file$Overdispersion <- overdispersion(result)
  # Estimated overdispersion.
  
  arg_output_file$Serial_correlation <- serial_correlation(result)
  # Estimated serial correlation.
  
  arg_output_file$Slope_imputed_mul <- slopes_imputed$slope$mul
  # Multiplicative imputed slope over the entire period.
  
  arg_output_file$Slope_imputed_mul_SE <- slopes_imputed$slope$se_mul
  # Standard error of multiplicative imputed slope over the entire period.
  
  arg_output_file$Slope_imputed_classification <- slopes_imputed$slope$meaning
  # Trend classification of multiplicative imputed slope over the entire period.
  
  arg_output_file$Slope_from_imputed_mul <- slopes_subperiod_imputed$slope$mul[slopes_subperiod_fitted$slope$from == arg_output_file$Year_from]
  # Multiplicative imputed slope over the subperiod.
  
  arg_output_file$Slope_from_imputed_mul_SE <- slopes_subperiod_imputed$slope$se_mul[slopes_subperiod_fitted$slope$from == arg_output_file$Year_from]
  # Standard error of multiplicative imputed slope over the subperiod.
  
  arg_output_file$Slope_from_imputed_classification <- slopes_subperiod_imputed$slope$meaning[slopes_subperiod_imputed$slope$from == arg_output_file$Year_from]
  # Trend classification of multiplicative imputed slope over the subperiod.
  
  return(arg_output_file)
}

#####################################################################################################################
makeOverview <- function(listSpeciesStratumCombinations) {
  # Aim: this function makes an empty dataframe to list failure or success of rtrim runs.
  
  numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)
  
  overview <- data.frame(ss_combinations = character(numberSpeciesStratumCombinations), 
                         # All files of combination of species and stratum present in the working directory.
                         
                         species_group = character(numberSpeciesStratumCombinations),
                         # Short for specific species group.
                         
                         species_number = integer(numberSpeciesStratumCombinations),
                         # Unique species number.
                         
                         stratum_type = character(numberSpeciesStratumCombinations),
                         # Stratumtype. 1 = standard stratum (use this), 2 = combination of strata (not relevant here).
                         
                         stratum_number = integer(numberSpeciesStratumCombinations),
                         # Unique number for a stratum.
                         
                         first_year = integer(numberSpeciesStratumCombinations),
                         # The first year in which counts are available for this specific combination of species and stratum.
                         
                         last_year = integer(numberSpeciesStratumCombinations),
                         # The last year in which counts are available for this specific combination of species and stratum.
                         
                         success = character(numberSpeciesStratumCombinations), 
                         # yes if the analysis were successful, no if not.
                         
                         attempt_1 = character(numberSpeciesStratumCombinations), 
                         # Outcome of the first attempt. The attempt may be successful or not.
                         
                         attempt_2 = character (numberSpeciesStratumCombinations), 
                         # Outcome of the second attempt.
                         # Outcome is only available if the first attempt was not successful, else the value will be n.a.
                         
                         attempt_3 = character(numberSpeciesStratumCombinations), 
                         # Outcome of the third attempt.
                         # Outcome is only avaialble if earlier attempts were not succesful, else the value will be n.a.
                         
                         attempt_4 = character(numberSpeciesStratumCombinations),
                         # Outcome of the fourth attempt.
                         # Outcome is only avaialble if earlier attempts were not succesful, else the value will be n.a.
                         
                         error_1 = character(numberSpeciesStratumCombinations),
                         # Error message of the first attempt.
                         
                         error_2 = character(numberSpeciesStratumCombinations),
                         # Error message of the second attempt.
                         
                         error_3 = character(numberSpeciesStratumCombinations),
                         # Error message of the third attempt.
                         
                         error_4 = character(numberSpeciesStratumCombinations))
  # Error message of the last attempt.
  
  overview$ss_combinations <- gsub(listSpeciesStratumCombinations, pattern = "_arg_input_stratum.csv", replacement = "")
  
  overview$success <- "no"
  overview$attempt_1 <- "n.a."                                                                                       
  overview$attempt_2 <- "n.a."
  overview$attempt_3 <- "n.a."
  overview$attempt_4 <- "n.a."
  overview$error_1 <- "n.a."
  overview$error_2 <- "n.a."
  overview$error_3 <- "n.a."
  overview$error_4 <- "n.a."
  
  # retrieve species number, stratumtype and statumnumber from the filename 
  
  overview$species_group <- gsub(listSpeciesStratumCombinations, pattern = "_[0-9]+_[0-9]+_[0-9]+_arg_input_stratum.csv", replacement = "")
  
  without_Species_Group <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_", replacement = "")
  
  overview$species_number <- as.integer(gsub(without_Species_Group, pattern = "_[0-9]+_[0-9]+_arg_input_stratum.csv", replacement = ""))
  
  without_Species_Group_without_Species_number <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_[0-9]+_", replacement = "")
  
  overview$stratum_type <- gsub(without_Species_Group_without_Species_number, pattern = "_[0-9]_arg_input_stratum.csv", replacement = "")
  
  without_Species_Group_without_Species_number_without_Stratum_type <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_[0-9]+_[0-9]+_", replacement = "")
  
  overview$stratum_number <- gsub(without_Species_Group_without_Species_number_without_Stratum_type, pattern = "_arg_input_stratum.csv", replacement = "")
  
  overview <- overview[order(overview$species_number, overview$stratum_number), ]
  
  return(overview)
}

#####################################################################################################################
make_All_Indices_All_Trends <- function(overview, listSpeciesStratumCombinations){
  # Aim: this function makes an empty dataframe to list indices and time totals of all species-stratum combinations.    
  
  numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)
  
  first_year_over_all_species <- min(overview$first_year)
  last_year_over_all_species <- max(overview$last_year)
  
  complete_period <- first_year_over_all_species:last_year_over_all_species
  columnnames <- as.character(complete_period)
  number_of_columns <- length(complete_period)
  number_of_rows <- numberSpeciesStratumCombinations * 4
  # Four rows for each combination of species and stratum (one row for each record type)
  
  temporary_matrix <- matrix(nrow = number_of_rows, ncol = number_of_columns)
  temporary_dataframe <- as.data.frame(temporary_matrix)
  colnames(temporary_dataframe) <- columnnames
  
  all_Indices_All_Trends <- data.frame(Year_of_analysis = integer(numberSpeciesStratumCombinations * 4),
                                       # Dataframe, named all_Indices_All_Trends, is created.
                                       # First column (Year_of_analysis): Most recent year for which counts are available.
                                       
                                       Species_number = integer(numberSpeciesStratumCombinations),
                                       # Unique number for a species.
                                       
                                       Stratum_number = integer(numberSpeciesStratumCombinations),
                                       # Unique number for a stratum.
                                       
                                       Recordtype_number = integer(numberSpeciesStratumCombinations),
                                       # Unique number for record types. There are four record types: 
                                       # 1: indices. 
                                       # 2: standard errors of indices.
                                       # 3: time totals.
                                       # 4: standard errors of time totals.
                                       
                                       Recordtype_name = character(numberSpeciesStratumCombinations),	
                                       # Name of the record type.
                                       
                                       N_sites = integer(numberSpeciesStratumCombinations),	
                                       # Number of unique sites in the counts file used.
                                       
                                       Slope_imputed_mul =	numeric(numberSpeciesStratumCombinations),
                                       # Multiplicative imputed slope over the entire period.
                                       
                                       Slope_imputed_mul_SE = numeric(numberSpeciesStratumCombinations),	
                                       # Standard error of the multiplicative imputed slope over the entire period.
                                       
                                       Slope_imputed_classification = character(numberSpeciesStratumCombinations),	
                                       # Trend classification of the imputed slope over the entire period.
                                       
                                       Slope_from_imputed_mul = numeric(numberSpeciesStratumCombinations),
                                       # Multiplicative imputed slope over the subperiod.
                                       
                                       Slope_from_imputed_mul_SE = numeric(numberSpeciesStratumCombinations),	
                                       # Standard error of the multiplicative imputed slope over the subperiod.
                                       
                                       Slope_from_imputed_classification = character(numberSpeciesStratumCombinations),	
                                       # Trend classification of the imputed slope over the subperiod.
                                       
                                       Year_from = integer(numberSpeciesStratumCombinations),	
                                       # First year of the subperiod over which a slope has been calculated. 
                                       
                                       Date_analysis = character(numberSpeciesStratumCombinations))
  # Date at which the analysis was executed.
  
  # retrieve species code, stratumtype and statumnumber from the filename 
  
  without_Species_Group <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_", replacement = "")
  
  species_Code  <- as.integer(gsub(without_Species_Group, pattern = "_[0-9]_[0-9]+_arg_input_[a-zA-Z]+.csv", replacement = ""))
  
  without_Species_Group_Without_Species_Number <- gsub(listSpeciesStratumCombinations, pattern = "[A-Z]{3,4}_[0-9]+_", replacement = "")
  
  stratumtype_and_Stratumnumber <- gsub(without_Species_Group_Without_Species_Number, pattern = "_arg_input_[a-zA-Z]+.csv", replacement = "")
  
  stratumtype <- gsub(stratumtype_and_Stratumnumber, pattern = "_[0-9]+", replacement = "")
  
  stratum_Number <- gsub(stratumtype_and_Stratumnumber, pattern = "[0-9]+_", replacement = "")
  
  all_Indices_All_Trends$Recordtype_number <- rep(1:4, each = numberSpeciesStratumCombinations)
  # Unique number for record types. There are four record types: 
  # 1: indices. 
  # 2: standard errors of indices.
  # 3: time totals.
  # 4: standard errors of time totals.
  number_of_record_types <- length(unique(all_Indices_All_Trends$Recordtype_number))
  
  all_Indices_All_Trends$Species_number <- rep(species_Code, number_of_record_types)
  all_Indices_All_Trends$Stratum_number <- rep(stratum_Number, number_of_record_types)
  
  names_record_types <- c("indices", "se_indices", "time_totals", "se_time_totals")
  all_Indices_All_Trends$Recordtype_name <- rep(names_record_types, each = numberSpeciesStratumCombinations)	
  # Name of the record type.
  
  all_Indices_All_Trends$Slope_imputed_classification <- ""
  all_Indices_All_Trends$Slope_from_imputed_classification <- ""
  all_Indices_All_Trends$Date_analysis <- format(Sys.Date(), "%d-%m-%Y")
  # Date at which the analysis was done.
  
  all_Indices_All_Trends <- cbind(all_Indices_All_Trends, temporary_dataframe)
  
  return(all_Indices_All_Trends)  
}

#####################################################################################################################
fill_All_Indices_All_Trends <- function(result, arguments, j, listSpeciesStratumCombinations, all_Indices_All_Trends){
  # Aim: this function fills the dataframe with indices and time totals of all species-strata combinations.      
  
  numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)
  
  indices <- index(result, which = "both", covars = FALSE, base = arguments$Base_year_first_year:arguments$Base_year_last_year)
  time_totals <- totals(result, which = "both")
  slopes_imputed <- overall(result, which = "imputed")
  slopes_subperiod_imputed <- overall(result, which = "imputed", changepoints = c(arguments$Year_from))
  
  names_columns <- as.character(indices$time)
  position_colums <- colnames(all_Indices_All_Trends) %in% names_columns
  
  all_Indices_All_Trends$Year_of_analysis <- max(indices$time)
  # Most recent year in which counts are available.
  
  all_Indices_All_Trends$N_sites[j] <- result$nsite
  all_Indices_All_Trends$N_sites[j + numberSpeciesStratumCombinations] <- result$nsite
  all_Indices_All_Trends$N_sites[j + numberSpeciesStratumCombinations * 2] <- result$nsite
  all_Indices_All_Trends$N_sites[j + numberSpeciesStratumCombinations * 3] <- result$nsite
  # Number of unique sites in the counts file used to call rtrim.
  
  all_Indices_All_Trends$Slope_imputed_mul[j] <- slopes_imputed$slope$mul
  # Multiplicative imputed slope over the entire period.
  
  all_Indices_All_Trends$Slope_imputed_mul_SE[j] <-	slopes_imputed$slope$se_mul
  # Standard error of the multiplicative imputed slope over the entire period.
  
  all_Indices_All_Trends$Slope_imputed_classification[j] <-	slopes_imputed$slope$meaning
  # Trend classification of the imputed slope over the entire period.
  
  all_Indices_All_Trends$Slope_from_imputed_mul[j] <- slopes_subperiod_imputed$slope$mul[slopes_subperiod_imputed$slope$from == arguments$Year_from]
  # Multiplicative imputed slope over the subperiod.
  
  all_Indices_All_Trends$Slope_from_imputed_mul_SE[j] <- slopes_subperiod_imputed$slope$se_mul[slopes_subperiod_imputed$slope$from == arguments$Year_from]	
  # Standard error of the multiplicative imputed slope over the subperiod.
  
  all_Indices_All_Trends$Slope_from_imputed_classification[j] <- slopes_subperiod_imputed$slope$meaning[slopes_subperiod_imputed$slope$from == arguments$Year_from]	
  # Trend classification of the imputed slope over the subperiod.
  
  all_Indices_All_Trends$Year_from[j] <-	arguments$Year_from
  # First year of the subperiod over which a slope has been calculated. 
  
  all_Indices_All_Trends[j, position_colums] <- round(100 * indices$imputed, 1)
  all_Indices_All_Trends[j + numberSpeciesStratumCombinations, position_colums] <- round(100 * indices$se_imp, 2)    #12/6/2020 Eva Silarova added 100 into <- round(indices$se_imp, 2) to have SE corresponding with index in % - v1.2
  all_Indices_All_Trends[j + numberSpeciesStratumCombinations * 2, position_colums] <- round(time_totals$imputed, 2)
  all_Indices_All_Trends[j + numberSpeciesStratumCombinations * 3, position_colums] <- round(time_totals$se_imp, 2)
  
  return(all_Indices_All_Trends)  
}

#####################################################################################################################

agreggation<-function(datos,weight){
  if(weight=="NONE"){
    output <- datos                             
    names(output)<-c("species", "euring", "site", "year","count","weights")
    includeWeights <<- "TRUE"                             
    Presence_weights <<- "TRUE"
    
  }
  #MAX: selects the maximun count for each location/year combination.  
  if(weight=="MAX"){
    output <- aggregate(x=datos$count, by = datos[,c("species", "euring", "site", "year")], FUN = max)                             
    names(output)<-c("species", "euring", "site", "year","count")
    output$weights <- rep(1,length(output$count))                 
    includeWeights <<- "FALSE"                             
    Presence_weights <<- "FALSE"
  }
  #AVE_DIRECT: provides mean number of observed counts among all the visits for each location/year combination.
  if(weight=="AVE_DIRECT"){
    output <- aggregate(x=datos$count, by = datos[,c("species", "euring", "site", "year")], FUN = mean)      
    names(output)<-c("species", "euring", "site", "year","count")
    output$weights <- rep(1,length(output$count))                 
    includeWeights <<- "FALSE"                             
    Presence_weights <<- "FALSE"
  }
  #AVE_WEIGHTS: calculates the sum of count for each location/year combination, it also calculates weight in the form:  
  #             1/number of visits per year
  if(weight=="AVE_WEIGHTS"){
    output <- aggregate(x=datos$count, by = datos[,c("species", "euring", "site", "year")], FUN = sum)                 
    names(output)<-c("species", "euring", "site", "year","count")
    output_ncensus <- aggregate(x=datos$year, by = datos[,c("species", "euring", "site","year")], FUN = length) 
    output$ncensus<-as.numeric(output_ncensus$x)
    output$weights <- 1/output$ncensus                      
    output<-output[,-6]
    includeWeights <<- "TRUE"                              
    Presence_weights <<- "TRUE"
    
  }
  #FIRST: selects the earliest visit for each location/year combination
  
  if(weight=="FIRST"){
    datos1<-datos[,c("species", "euring", "site", "year","count","date")]
    output <- aggregate(x=datos1$date, by = datos1[,c("species", "euring", "site", "year")],FUN=min) 
    for(j in 1:nrow(output)){
      
      output$count[j]<-datos1[datos1$site ==output$site[j] & datos1$date==output$x[j],"count"]
    }   
    output<-output[,c(1:4,6)]
    names(output)<-c("species", "euring", "site", "year","count")
    output$weights <- rep(1,length(output$count))
    includeWeights <<- "FALSE"                             
    Presence_weights <<- "FALSE"
  }
  #LAST: selects the latest visit for each location/year combination
  if(weight=="LAST"){
    datos1<-datos[,c("species", "euring", "site", "year","count","date")]
    output <- aggregate(x=datos1$date, by = datos1[,c("species", "euring", "site", "year")],FUN=max) 
    for(j in 1:nrow(output)){
      
      output$count[j]<-datos1[datos1$site ==output$site[j] & datos1$date==output$x[j],"count"]
    }   
    output<-output[,c(1:4,6)]
    names(output)<-c("species", "euring", "site", "year","count")
    output$weights <- rep(1,length(output$count))
    includeWeights <<- "FALSE"                             
    Presence_weights <<- "FALSE"
  }
  return(output)
}

#####################################################################################################################
visited_sites<-function(datos,scheme){
  site <-unique(datos[,c("site")])
  cat("The total number of locations for the scheme",scheme ,"is: ",length(site),"\n")
  datos <- datos[order(datos$site, datos$year),] # ordered by location and year
  location_sp <- unique(datos[,c("species", "euring", "site")]) # all species for each location
  visit_location <- unique(datos[,c("site", "year")]) # all years for each location
  minyear <- min(datos$year)
  maxyear <- max(datos$year)
  d <-as.data.frame(seq(from = minyear, to = maxyear, by = 1)) # this line is to add missing years if you have no data at all
  colnames(d) <- "year"
  visit_location <- merge(visit_location, d, by.x = "year", by.y = "year", all.y = TRUE)
  dd <- merge(location_sp, visit_location, by = "site", all = TRUE) # adds species for missing location and missing year
  ddd <- merge(dd, datos, by = c("species", "euring", "year", "site"), all = TRUE) # adds number "n" for each species. If there is no observation, then n = NA
  ddd <- ddd[,c("species","euring", "year", "site", "count", "weights")]
  ddd[ddd$count %in% NA, "count"] <- 0 # we change NA into 0 : if a species is not seen/heard during the sampling we consider that it is not present
  year_species <- unique(datos[,c("species", "year"),]) # list of each year for each species
  year <- as.vector(unique(d$year)) # this vector shows all the years
  yy <- length(year) # this is the number of years for the project
  loca <- as.data.frame(unique(datos$site)) # dataframe of the locations
  loca_all_year <- as.data.frame(loca[rep(1:nrow(loca),each = c(yy)),]) # creates a table with a repetition of all the locations for each year
  years <- as.data.frame(unique(visit_location$year))
  loca_unique <- as.vector(unique(datos$site))
  ll <- length(loca_unique)
  year_all_loca <- as.data.frame(years[rep(seq_len(nrow(years)), times = c(ll)),]) # replication of the years as long as the location list
  loca_year <- cbind(loca_all_year, year_all_loca) #each location, each year
  colnames(loca_year) <- c("site", "year")
  sp_year_loca <- merge(loca_year, location_sp, by = "site", all = TRUE) # each species for each location and each year
  
  #Creation of the worktable
  worktable <- merge(ddd, sp_year_loca, by = c("site", "year", "species", "euring"), all.y = TRUE) # adds species and location for the missing years, n is NA for the missing years
  # 20240220 John Kennedy Populate weights for species with a count of 0 in a location/year, or in alocation not counted in a year
  if(weight=="NONE"){
    # extract the weight for each location for each year. Assumes weight does not vary by species in a location/year.
    location_year_weights <- unique(datos[,c("site", "year", "weights")]) # all locations, years and weights
    
    # assign relevant weight to every 0 or non-0 count in every location/year visited
    # drop weights column from worktable before merging in weights
    worktable <- merge(worktable[,-6], location_year_weights, by = c("site", "year"), all.x = TRUE)    
    
  }
  #13/03/2024  J.Rivas we extracted the lane john kennedy out of if (weight==none){} to deal with NA weights
  # if a location was not visited at all in a year, it will have a count of NA and a weight of NA.
  # RTRIM cannot handle weight of NA so assign a weight of 1
  worktable$weights <- replace(worktable$weights, is.na(worktable$weights), 1) 
  colnames(worktable) <- c("site", "year", "species", "euring", "count", "weights") # renames the columns
  worktable$count <- as.numeric(worktable$count)
  worktable$year <- as.numeric(worktable$year)
  worktable$species <- as.factor(worktable$species)
  worktable$euring <- as.numeric(worktable$euring)
  worktable$weights <- as.numeric(worktable$weights)
  return(worktable)
}

#####################################################################################################################
#This function read the control files , the raw data and also starts the logfile.
read.files<-function(){
  #Obtaining the version from the RTRIM folder
  pattern <- "(?:\\d+(?:\\.\\d+)?(?:\\.\\d+)?(?:\\.\\d+)?(?:\\.\\d+)?(?:\\.\\d+)?(?:\\.\\d+)?(?:\\.\\d+)?)" # this regular expression shall select up to 3 numbers separated by dots in the folder path
  version<- regmatches(rtrim_folder, regexpr(pattern, rtrim_folder))
  
  #Creation of logfile
  sink("Logfile RTRIM-shell.txt")
  cat("##################################","\n")
  cat(date(),"\n")  
  cat("\n")
  cat("Start TRIM-shell v",version," in R", "\n")  
  cat("\n")
  #Reading Analyses file
  analyses<-read.table("analyses.csv",sep=";",dec=".",header = F)
  
  if (ncol(analyses)==1){
    analyses <- read.table("analyses.csv",sep=",",dec=".",header = F)
  }
  Country <<- analyses[analyses$V1=="country",2] 
  country_code <<- analyses[analyses$V1=="country_code",2] 
  changepoints <<- analyses[analyses$V1=="changepoints",2] 
  data_file_name <<- analyses[analyses$V1=="data_file_name",2] 
  if(data_file_name== "pipeline"){
    execute<<-0
  }else{
    execute<<-1
  }
  decimal<<-analyses[analyses$V1=="decimal",2] 
  separator<<-analyses[analyses$V1=="separator",2] 
  prefix<<-analyses[analyses$V1=="monitoring_prefix",2]
  working_directory<<-getwd()
  cat("Analyses start for ", Country,"\n")
  cat("\n")
  if(execute==1){
    #Reading scheme file
    scheme<<-read.table("scheme.csv",sep=";",dec=".",header = T)
    if (ncol(scheme)==1){
      scheme<<-read.table("scheme.csv",sep=",",dec=".",header = T)
    }
    cat("RTRIM will run for:")
    cat("\n")
    for(i in 1:nrow(scheme)){
      cat(scheme$Species_nr[i]," species with ", scheme$Count_aggregation[i], "aggregation","\n" )
    }
    cat("\n")
    #generating the the exemption object (additional runs for RTRIM)
    exemptions<-as.data.frame(matrix(NA,ncol = 8,nrow = 1))
    names(exemptions)<-c("Species_nr","Monitoring_scheme","Base_year_first_year","Base_year_last_year","year_from","Count_aggregation","Serial_correlation","Overdispersion")
    exemptions<<-exemptions
    n_schemes<<-unique(scheme$Monitoring_scheme)
    main_scheme<<-scheme[scheme$Species_nr=="all","Monitoring_scheme"]
    if(nrow(scheme) >1){
      Base_year_first_year <<- scheme[scheme$Species_nr =="all","Base_year_first_year"]
      Base_year_last_year <<- scheme[scheme$Species_nr =="all","Base_year_last_year"]
      year_from <<- scheme[scheme$Species_nr =="all","year_from"]
      weight <<- scheme[scheme$Species_nr =="all","Count_aggregation"]
      Serial_correlation_status<<-scheme[scheme$Species_nr =="all","Serial_correlation"]
      Overdispersion_status<<-scheme[scheme$Species_nr =="all","Overdispersion"]
      extra_schemes<-scheme[scheme$Species_nr!="all",]
      extra_spp<-extra_schemes$Species_nr
      
      
      z=1
      for(i in 1:length(extra_spp)){
        extra_sp<-unlist(strsplit(extra_spp[i], split = "-"))
        for(j in 1:length(extra_sp)){
          exemptions[z,"Species_nr"]<-extra_sp[j]
          exemptions[z,"Monitoring_scheme"]<-extra_schemes[i,"Monitoring_scheme"]
          exemptions[z,"Base_year_first_year"]<-extra_schemes[i,"Base_year_first_year"]
          exemptions[z,"Base_year_last_year"]<-extra_schemes[i,"Base_year_last_year"]
          exemptions[z,"year_from"]<-extra_schemes[i,"year_from"]
          exemptions[z,"Count_aggregation"]<-extra_schemes[i,"Count_aggregation"]
          exemptions[z,"Serial_correlation"]<-extra_schemes[i,"Serial_correlation"]
          exemptions[z,"Overdispersion"]<-extra_schemes[i,"Overdispersion"]
          z=z+1
        }
      }
      exemptions<<-exemptions
    } else{
      Base_year_first_year <<- scheme[scheme$Species_nr =="all","Base_year_first_year"]
      Base_year_last_year <<- scheme[scheme$Species_nr =="all","Base_year_last_year"]
      year_from <<- scheme[scheme$Species_nr =="all","year_from"]
      weight <<- scheme[scheme$Species_nr =="all","Count_aggregation"]
      Serial_correlation_status<<-scheme[scheme$Species_nr =="all","Serial_correlation"]
      Overdispersion_status<<-scheme[scheme$Species_nr =="all","Overdispersion"]
    }
    #Reading the raw data
    datos <<- read.table(paste0(working_directory,"/00_Raw_data/", data_file_name), stringsAsFactors = FALSE, header = TRUE, dec = decimal,sep=separator)   
    names_columns<<-c("species","euring","site","count","year")
    if(is.character(datos$n)==T){
      datos <<- read.table(paste0(rtrim_folder,"00_Raw_data/", data_file_name), stringsAsFactors = FALSE, header = TRUE, dec = decimal,sep=separator)  
    }
    
    if(is.element("FIRST",scheme$Count_aggregation)){
      names_columns<<-c("species","euring","site","count","year","date")
      if(!is.element("date",names(datos))){
        cat ("To use First visit your raw data must have a date column")
      }
      datos$date<-as.Date(datos$date,format = "%d.%m.%Y")
    }
    if(is.element("LAST",scheme$Count_aggregation)){
      names_columns<<-c("species","euring","site","count","year","date")
      if(!is.element("date",names(datos))){
        cat ("To use Last visit your raw data must have a date column")
      }
      datos$date<-as.Date(datos$date,format = "%d.%m.%Y")
    }
  }
  #END of read.files function
}


#This function transfor raw data into RTRIM-shell inputs
data.preparation <- function(datos,weight,execute){
  if (execute==1){
    cat("##############################################","\n")
    cat("#Data preparation function:","\n")
    cat("##############################################","\n")
    datos$species <- as.factor(datos$species)
    datos$site <- as.factor(datos$site)
    datos<- datos[!is.na(datos$count),]
    
    if (length(n_schemes)>1){
      data_list<-vector("list",length= length(n_schemes))
      for(i in 1:length(data_list)){
        if(i==1){
          select_spp<-exemptions[exemptions$Monitoring_scheme!=main_scheme,"Species_nr"]
          data_list[[i]]<-datos[!datos$euring %in% select_spp,]
          names(data_list)[i]<-main_scheme
        }else{  
          select_spp<-exemptions[exemptions$Monitoring_scheme==n_schemes[i],"Species_nr"]
          data_list[[i]]<-datos[datos$euring %in% select_spp,]
          names(data_list)[i]<-n_schemes[i]
        }
      }
    }else{
      data_list<-vector("list",length= length(n_schemes))
      names(data_list)[1]<-main_scheme
      data_list[[1]]<-datos
    }
    
    
    main_exemptions<-exemptions[exemptions$Monitoring_scheme==main_scheme,]
    other_exemptions<-exemptions[exemptions$Monitoring_scheme!=main_scheme,]
    data.test<-as.data.frame(matrix(NA))
    worktable<-as.data.frame(matrix(NA))
    
    for(i in 1:length(data_list)){
      if(i==1){
        datos<-as.data.frame(data_list[1])
        names(datos)<-names_columns
        data.test<-agreggation(datos,weight)
        if(!is.na(main_exemptions$Species_nr[1])==T){
          for(j in 1:nrow(main_exemptions)){
            data.test<-data.test[data.test$euring!=main_exemptions$Species_nr[j],]
            datos1<- datos[datos$euring==main_exemptions$Species_nr[j],]
            data.test1<-agreggation(datos1,main_exemptions$Count_aggregation[j])
            data.test<-rbind(data.test,data.test1)
          }
        }  
        #aqui agregamos los sites con 0 y NA
        worktable<-visited_sites(data.test,names(data_list)[i])
      }else{
        if(!is.na(other_exemptions$Species_nr[1])==T){
          datos1<-as.data.frame(data_list[i])
          names(datos1)<-names_columns
          exemptions_related<-other_exemptions[other_exemptions$Monitoring_scheme==names(data_list)[i],]
          data.test2<-as.data.frame(matrix(NA))
          for (j in 1:nrow(exemptions_related)) {
            
            datos_exemption<-datos1[datos1$euring==exemptions_related$Species_nr[j],]
            if(j==1){
              data.test2<-agreggation(datos_exemption,exemptions_related$Count_aggregation[j])
            }else{
              data.test3<-agreggation(datos_exemption,exemptions_related$Count_aggregation[j])
              data.test2<-rbind(data.test2,data.test3)
            }
          }
          
          #Aqui agrefgamos los sites con 0 y NA
          worktable2<-visited_sites(data.test2,names(data_list)[i])
          worktable<-rbind(worktable,worktable2)
        }
      }
    }  
    
    
    
    setwd(paste0(working_directory,"/02_Inputs/"))		### AJS, 16.5.2022 - enables saving the outputs into different folder
    
    write.table(worktable, file = "worktable.csv",quote=FALSE, row.names=FALSE, col.names=TRUE, sep = ";" )
    cat("Worktable was successfully created","\n")
    
    # creates a list to store the species results
    ListeSp <- as.list(unique(worktable$euring)) 
    dataPozit <- worktable[data.test$count > 0,]	### AJS, 16.5.2022 - subset of positive n's in data
    
    # a) Create the _counts.csv files for a species called "i" in the script. This is a loop that picks all the data for each species and creates the _counts.csv files.
    for(i in 1:length(ListeSp)){
      sprecordsubset<-subset(dataPozit,dataPozit$euring==ListeSp[i])
      firstRec<-min(unique(sprecordsubset$year))
      #firstRec <- min(unique(dataPozit[dataPozit$euring == ListeSp[i],]$year))	### AJS, 16.5.2022 - finds the 1st year when species was recorded
      subsetSP <- worktable[worktable$euring == ListeSp[i],] # subset of the data for the i species
      subsetSPY <- subsetSP[subsetSP$year >= firstRec,]		### AJS, 16.5.2022 - subset of the data from the firstRec-year further (drops out years at the beginning of the series with no record of the species i)
      
      ### AJS - changes -- include weights
      if (includeWeights) {
        subsetSp <- subsetSPY[,c("site", "year", "count", "weights")]         ### AJS, 16.5.2022 - changes in table-names (headers) - V5weights
        colnames(subsetSp) <- c("site", "year", "count", "weights")
      } else {
        subsetSp <- subsetSPY[,c("site", "year", "count")]                    ### AJS, 16.5.2022 - changes in table-names (headers) - V5weights
        colnames(subsetSp) <- c("site", "year", "count")
      }   ### AJS - end of change due to weights
      
      spp<-unlist(ListeSp[i])
      write.table (subsetSp, file = paste0 (prefix,"_",ListeSp[i],"_1_",country_code,"_counts.csv"), row.names=FALSE, col.names=TRUE, sep = ";", dec = ".")	### AJS, 17.5.2022, "dec" added.
      cat(prefix,"_",spp,"_1_",country_code,"_counts.csv was successfully created and contains:","\n")
      missing_data<-sum(is.na(subsetSp$count))  
      if(is.na(missing_data)==T){missing_data<-0}
      zeroes<-nrow(subset(subsetSp,subsetSp$count==0))
      if(is.na(zeroes)==T){zeroes<-0}
      positive<-nrow(subsetSp)-(missing_data+zeroes) 
      
      cat("NAs: ",missing_data ,"\n")
      cat("Zero values: ",zeroes,"\n")
      cat("Counts different from zero: ",positive,"\n")
      cat("\n")
      
    }
    # b) Create input stratum files with a loop 
    
    for(i in 1:length(ListeSp)){
      sprecordsubset<-subset(dataPozit,dataPozit$euring==ListeSp[i])
      firstRec<-min(unique(sprecordsubset$year))
      #firstRec <- min(unique(dataPozit[dataPozit$euring == ListeSp[i],]$year))                            ### AJS, 16.5.2022 - finds the 1st year when species was recorded
      
      aa <- c("File","Base_year_first_year","Base_year_last_year","Changepoints","Serial_correlation","Overdispersion","Presence_weights","Presence_monthfactors","Year_from","Save_fitted_values") 
      
      #This fragment modify the base year(first or last) and the year from according to the first register of the species
      spp<-unlist(ListeSp[i])
      if(is.element(ListeSp[[i]],exemptions$Species_nr)){
        exemption_data<-exemptions[exemptions$Species_nr==ListeSp[[i]],]
        # 20/02/2024 John Kennedy Presence_weights needs to be saved in the _arg_input_stratum.csv file
        #bb <- c(paste(prefix,"_",ListeSp[i],"_1_", country_code),exemption_data$Base_year_first_year,exemption_data$Base_year_last_year,changepoints,Serial_correlation_status,Overdispersion_status, "FALSE", "FALSE", exemption_data$year_from,"TRUE") # set columns for arg_input_stratum
        #13/03/2024  J.Rivas fixed the issue of why exemptions weren't loading the right overdispersion and serial correlation
        bb <- c(paste(prefix,"_",ListeSp[i],"_1_", country_code),exemption_data$Base_year_first_year,exemption_data$Base_year_last_year,changepoints,exemption_data$Serial_correlation,exemption_data$Overdispersion, Presence_weights, "FALSE", exemption_data$year_from,"TRUE") # set columns for arg_input_stratum
        if (firstRec > exemption_data$Base_year_first_year) {                    
          cat("The first base year of species ",spp, "has been modified to ",firstRec,"\n")
          bb[2] <- firstRec
        }
        if (firstRec > exemption_data$Base_year_last_year) {                     
          cat("The last base year of species ",spp, "has been modified to ",firstRec,"\n")
          bb[3] <- firstRec
        }
        if (firstRec > exemption_data$year_from) {                                                      
          cat("The year from of species ",spp, "has been modified to ",firstRec,"\n")
          bb[9] <- firstRec
        }  
      }else{
        # 20/02/2024 John Kennedy Presence_weights needs to be saved in the _arg_input_stratum.csv file
        #bb <- c(paste(prefix,"_",ListeSp[i],"_1_", country_code),Base_year_first_year,Base_year_last_year,changepoints,Serial_correlation_status,Overdispersion_status, "FALSE", "FALSE", year_from,"TRUE") # set columns for arg_input_stratum
        bb <- c(paste(prefix,"_",ListeSp[i],"_1_", country_code),Base_year_first_year,Base_year_last_year,changepoints,Serial_correlation_status,Overdispersion_status, Presence_weights, "FALSE", year_from,"TRUE") # set columns for arg_input_stratum
        if (firstRec > Base_year_first_year) {                    
          cat("The first base year of species ",spp, "has been modified to ",firstRec ,"\n")
          bb[2] <- firstRec
        }
        if (firstRec > Base_year_last_year) {                     
          cat("The last base year of species ",spp, "has been modified to ",firstRec,"\n")
          bb[3] <- firstRec
        }
        if (firstRec > year_from) {                                                      
          cat("The year from of species ",spp, "has been modified to ",firstRec,"\n")
          bb[9] <- firstRec
        }  
      }
      #Creation of argç_input file
      argInput <- rbind(aa,bb)
      argInput <- as.data.frame(argInput)
      colnames(argInput) <- c("File","Base_year_first_year","Base_year_last_year","Changepoints","Serial_correlation","Overdispersion","Presence_weights","Presence_monthfactors","Year_from","Save_fitted_values")
      argInput <- argInput[c(2),]
      argInput$File <- gsub(pattern = " ", replacement = "", x = argInput$File)
      write.table (argInput, file = paste0 (prefix,"_",ListeSp[i],"_1_",country_code, "_arg_input_stratum.csv"), row.names=FALSE, col.names=TRUE, sep = ";")
      cat(prefix,"_",spp,"_1_",country_code, "_arg_input_stratum.csv was successfully created","\n")
    }
    cat("##############################################","\n")
    cat(" End of Data preparation","\n")
    cat("##############################################","\n")
    
  }else{
    setwd(paste0(working_directory,"/02_Inputs/"))
    cat("##############################################","\n")
    cat(" Data preparation was skipped","\n")
    cat("##############################################","\n")
  }
}  

#This function run the RTRIM models for each species and save the result as .R file
rtrim.strata<- function(){
  cat("##############################################","\n")
  cat("Starting Rtrim Strata","\n")
  cat("##############################################","\n")
  for (j in 1:numberSpeciesStratumCombinations) {
    
    cat("Processing", listSpeciesStratumCombinations[j], "\n") # 19/07/2022 John Kennedy inserted instrumentation to find out what species/attempt each warning was for - v1.6
    
    #j<- 1
    
    # The file with arguments contains the information to analyse the counts for a particular combination of species and stratum. 
    # The arguments are used when calling the function 'rtrim'
    arguments <- read.csv2(listSpeciesStratumCombinations[j], header = TRUE, stringsAsFactors = FALSE)
    
    arguments<-as.list(arguments)                                                  # Enables the script to read the change points properly, 
    if(!arguments$Changepoints%in%c("all","auto")) {                               # when they are specified in the ?_arg_input_stratum.csv? file as comma-separated numbers.
      Changepoints<-as.integer(unlist(strsplit(arguments$Changepoints,"-")))       # Function was created 20/03/2022 by Dario Massimino -> modified 11/04/2022 by Meelis Leivits -> modified 20/07/2022 by John Kennedy - v1.6
      arguments$Changepoints<-Changepoints                                         # Function was adjusted to use of "-" in the changepoints # 9/02/2023 Javier Rivas Salvador - v1.6
    } 
    
    # The file with counts contains the counts for a particular combination of species and stratum.
    # Weights may also be present in this file.
    #13/01/2023 Javier Rivas: I have modified the following lines to ensure that independently of the decimal symbol the counts$count variable is loaded as numeric
    counts <- read.table(paste(arguments$File, "_counts.csv", sep = ""), stringsAsFactors = FALSE, header = TRUE, dec = ".",sep=";")   
    
    if(is.character(counts$count)==T | is.character(counts$weight)==T){
      counts <- read.table(paste(arguments$File, "_counts.csv", sep = ""), stringsAsFactors = FALSE, header = TRUE, dec = ",",sep=";")   
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
      
      save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
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
        
        save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
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
          
          save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
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
            
            save(x = result,  file = paste(arguments$File, ".RData", sep = ""))
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
  write.table(overview, "overview.csv", row.names = FALSE,sep=";",dec=".")   #13/01/2023 Javier Rivas: Changes to produce outputs with dots
  cat("Overview table generated successfully","\n")
  cat("##############################################","\n")
  cat("End of Rtrim Strata","\n")
  cat("##############################################","\n")
}

#This function creates the outputs of RTRIM
processing.output<-function(){
  
  cat("##############################################","\n")
  cat("Start of Processing Outputs","\n")
  cat("##############################################","\n")
  outputs<-paste0(rtrim_folder,"03_Outputs\\")
  #####################################################################################################################
  # SELECTING SUCCESSFUL RUNS.
  #####################################################################################################################
  # Determines which datafiles (combinations of species and stratum) have been analysed successfully. 
  # This information is found in "overview.csv", which has been produced by the script that called rtrim.
  
  overview <- read.table("overview.csv",sep=";",dec=".",header = T)   #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
  listsuccessfulAnalyses <- overview$ss_combinations[overview$success == "yes"]
  listSpeciesStratumCombinations <- paste(listsuccessfulAnalyses, "_arg_input_stratum.csv", sep = "")
  
  # Determine how many combinations of species and stratum have been analysed successfully.
  numberSpeciesStratumCombinations <- length (listSpeciesStratumCombinations)
  
  # File with trends and indices for each combination of species and stratum
  all_Indices_All_Trends <- make_All_Indices_All_Trends(overview = overview, listSpeciesStratumCombinations = listSpeciesStratumCombinations)
  
  #####################################################################################################################
  # PROCESSING OUTPUT OF SUCCESFUL RUNS.
  #####################################################################################################################
  
  for (j in 1:numberSpeciesStratumCombinations) {
    cat("Generating outputs for sp: ",listSpeciesStratumCombinations[j],"\n" )
    # The file with arguments contains the arguments to run the analysis for a particular combination of species and stratum (stratum is e.g. a region). 
    # These arguments are used when calling the function 'rtrim'.
    arguments <- read.table(listSpeciesStratumCombinations[j], header = TRUE, stringsAsFactors = FALSE,sep=";",dec=".")
    counts <- read.table(paste(arguments$File, "_counts.csv", sep = ""), header = TRUE,dec=".",sep=";") 
    if(is.character(counts$count)==T){
      counts <- read.table(paste(arguments$File, "_counts.csv", sep = ""), stringsAsFactors = FALSE, header = TRUE, dec = ",",sep=";")   
    } 
    load(paste(arguments$File, ".RData", sep = "")) # produces object with name "result"
    
    #####################################################################################################################
    # CREATING FILES FOR LATER USE. 
    #####################################################################################################################
    
    # Several output files are created (dataframes). Filling the files with output is done at a later stage.
    # To create the file containing indices and time totals (indices_TT_file) requires the file with results.
    # The file arg_output contains the slopes and arguments used to run the rtrim function.
    
    indices_TT_file <- make_Indices_TT_file(result = result)
    
    arg_output_file <- make_arg_output_file(arguments = arguments)
    
    #####################################################################################################################
    # FILLING OUTPUT.
    # Several functions are used to fill the output files.
    # These functions can be found in "functions_Shell_Rtrim.r".
    # This script goes to that file to include them here, which makes it easier to see what is done in this script. 
    #####################################################################################################################
    
    indices_TT_file <- fill_Indices_TT_file(indices_TT_file = indices_TT_file, result = result, arguments = arguments)
    
    arg_output_file <- fill_arg_output_file(arg_output_file = arg_output_file, result = result, counts = counts)
    arg_output_file$Changepoints<-gsub(",\\ ", "-", arg_output_file$Changepoints) # changes the separator from function "results" from comma to dash # Javier Rivas Salvador 9/2/2023 - v1.6 & sub changed to gsub by Javier Rivas 23/3/2023 according to Meelis Leivits
    all_Indices_All_Trends <- fill_All_Indices_All_Trends(all_Indices_All_Trends = all_Indices_All_Trends, result = result, arguments = arguments, j = j, listSpeciesStratumCombinations = listSpeciesStratumCombinations)
    
    covariant_matrix <- vcov(result)
    
    if (arguments$Save_fitted_values){
      
      FI <- results(result)
      
    }
    
    
    
    #####################################################################################################################
    # WRITING OUTPUT FILES. 
    #####################################################################################################################
    # Indices and time totals.
    name_Indices_TT_file <- paste(outputs,arguments$File, "_indices_TT.csv", sep = "")
    write.table(indices_TT_file, name_Indices_TT_file, row.names = FALSE,sep=";",dec=".")   #13/01/2023 Javier Rivas: Changes to produce outputs with dots  - v1.6
    cat("Indices_TT generated.","\n")
    ############################################################ 
    # Slopes entire period and arguments.
    name_arg_output_file <- paste(outputs,arguments$File, "_arg_output.csv", sep = "") 
    write.table(arg_output_file, name_arg_output_file, row.names = FALSE,sep=";",dec=".")   #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
    cat("Arg_output generated","\n")
    ############################################################
    # Covariant matrix.
    name_covariant_matrix <- paste(outputs,arguments$File, "_ocv.csv", sep = "")
    write.table(covariant_matrix, name_covariant_matrix, row.names = FALSE,sep=";",dec=".")   #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
    cat("ocv generated","\n")
    ############################################################
    # File with fitted values.
    if (arguments$Save_fitted_values){
      name_Fitted_Values_File <- paste(outputs,arguments$File, "_fitted_values.csv", sep = "")
      write.table(FI, name_Fitted_Values_File, row.names = FALSE,sep=";",dec=".")  #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
      cat("fitted values generated","\n")
    }
    ############################################################
    # 25/01/2022 Javier Rivas: Generating Tables 1-3 for the PECBMS coordinators - v1.4 IMPORTANT!!: version is extracted out of the path indicated in folder, if coordinators change the name it wont work.
    name_rawdata=paste(outputs,arguments$File, "_summarizing_tables.txt", sep = "")
    pattern <- "(?:\\d+(?:\\.\\d+)?(?:\\.\\d+)?(?:\\.\\d+)?(?:\\.\\d+)?(?:\\.\\d+)?(?:\\.\\d+)?(?:\\.\\d+)?)" # this regular expression shall select up to 3 numbers separated by dots in the folder path
    version<- regmatches(rtrim_folder, regexpr(pattern, rtrim_folder))
    write(paste0("\n", arguments$File,"- Rtrim-shell version: ",version, "\n"), file = name_rawdata)
    
    #18/08/2023 Javier Rivas: Code to include the version of RTRIM-shell into the summarazing tables
    
    # write(paste0("\n Rtrim-shell version: ",version,"\n"), file = name_rawdata, append = T)
    
    # END 18/08/2023  Javier Rivas
    
    
    # We are filling the  table on numbers of sites and counts
    write(paste0("\n", arguments$File, "\n"), file = name_rawdata,append = T)
    
    Table1=data.frame(Parameter=c("Sites","Times","Number of observerd zero counts","Number of observed positive counts",
                                  "Total number of observed counts", "Number of missing counts", "Total number of counts"),
                      Value= c(arg_output_file$N_sites,arg_output_file$N_time_values,arg_output_file$N_observed_zero_counts,
                               arg_output_file$N_observed_positive_counts,(arg_output_file$N_observed_positive_counts+arg_output_file$N_observed_zero_counts),
                               arg_output_file$N_missing_counts,arg_output_file$N_counts_total))
    
    
    
    # Table showing those sites containing more than the 10% of observations
    
    total_count=sum(counts$count[!is.na(counts$count)])
    tenpercent=(total_count*10)/100
    
    Table2=as.data.frame(matrix(NA,ncol = 3))
    names(Table2)=c("Site Number","Observed total","%")
    xx=subset(counts,!is.na(counts$count)==T)
    y=unique(xx$site)
    for(i in 1:length(unique(xx$site))){
      x=subset(xx,xx$site==y[i])
      s=sum(x$count)
      if (s>tenpercent){
        if(is.na(Table2[1,1])==T){
          Table2$`Site Number`=x$site[1]
          Table2$`Observed total`=s
          percentage=round((s*100/total_count),digits = 2)  # 18/03/22 Javier Rivas rounding up to 2 decimal - v1.5
          Table2$`%`=percentage
        }else{
          percentage=round((s*100/total_count),digits = 2)  # 18/03/22 Javier Rivas rounding up to 2 decimal - v1.5
          vector=c(x$site[1],s,percentage)
          Table2=rbind(Table2,vector)
        }
      }
    }  
    
    # Table 3 is Timepoint average table
    
    Table3=as.data.frame(matrix(NA,ncol=5,nrow=length(unique(counts$year[!is.na(counts$count)]))))
    
    names(Table3)=c("Timepoint","Observations","Average","Index","Real_Number")
    years=unique(counts$year[!is.na(counts$count)])# 15/03/2022 Javier Rivas changed the years, to exclude the years with full NA to prevent mismatches between the table and the loop - v1.5
    
    sites= unique(xx$site) #18/03/2022 Javier Rivas - To imitate the clever old trim we have to add this so we exclude full zero sites prior building the table - v1.5
    for(i in 1:length(sites)){
      xy=subset(xx,xx$site==sites[i])
      xz=sum(xy$count,na.rm=T)
      
      if(xz==0){ 
        xx=subset(xx,xx$site!=sites[i])
      }
    }#End of 18/03/22 edition
    
    for(i in 1:length(years)){
      Table3$Timepoint[i]=years[i]
      visitedsites=subset(xx,xx$year==years[i])
      Table3$Observations[i]=nrow(visitedsites)
      Table3$Real_Number[i]=sum(visitedsites$count)
      Table3$Average[i]=round(Table3$Real_Number[i]/Table3$Observations[i],digits = 2)  #18/03/22 Javier Rivas rounding up to 2 decimal - v1.5
      if(i==1){
        Table3$Index[i]=round(1.00,digits=2)  #18/03/22 Javier Rivas rounding up to 2 decimal - v1.5
      }else{
        Table3$Index[i]=round(Table3$Average[i]/Table3$Average[1],digits=2) #18/03/22 Javier Rivas rounding up to 2 decimal - v1.5
      }
      
    } 
    #This is the last calculation for Table 1 that requires data from table 3
    x=0
    for(i in 1:nrow(Table3)){
      s=Table3$Observations[i]*Table3$Average[i]
      x=x+s
    }
    Table3=Table3[order(Table3$Timepoint,decreasing = F),] #30/05/23 Javier Rivas command sorts the _indices_TT.csv for correct order of the years (the years were disordered and did not start by first year of time series)
    
    Total_count=c("Total count",x)
    Total_count
    Table1=rbind(Table1,Total_count)
    #Filling the txt!
    write("\n 1.	Summarizing table on numbers of sites and counts \n", file = name_rawdata, append = T)
    capture.output( print.data.frame(Table1, row.names=F, print.gap=3, quote=F, right=F),  file=name_rawdata,append = T  )
    write("\n 2.	Sites containing more than 10% of the total count \n", file = name_rawdata, append = T)
    capture.output( print.data.frame(Table2, row.names=F, print.gap=3, quote=F, right=F),  file=name_rawdata,append = T  )
    write("\n 3.	Time Point Averages table \n", file = name_rawdata, append = T)
    capture.output( print.data.frame(Table3, row.names=F, print.gap=3, quote=F, right=F),  file=name_rawdata,append = T  )
  }
  cat("Summarizing tables successfully created","\n")
  all_Indices_All_Trends <- all_Indices_All_Trends[order(all_Indices_All_Trends$Species_number, all_Indices_All_Trends$Stratum_number, all_Indices_All_Trends$Recordtype_number), ]
  write.table(all_Indices_All_Trends, paste0(outputs,"All_Indices_All_Trends.csv"), row.names = FALSE,sep=";",dec=".")           # 09/04/2021 Eva Silarova changed "all" to capital "All": "all_Indices_All_Trends.csv" was changed to "All_Indices_All_Trends.csv" - v1.3
  cat("All indices all trends successfully created","\n")
  #18/09/2023 Javi: Creation of the Schedule table
  #Monitoring_scheme table
  if(execute==1){
    cat("creating monitoring scheme table","\n")
    species_list_pcbms=read.table(paste0(rtrim_folder,"01_Scripts\\DO_NOT_MODIFY\\species_list.csv"),sep=";",dec=".",header = T)
    names(species_list_pcbms)<-c("Euring","Latin.Name","English.Name","Is.PECBMS.Sp.")
    sp<-unique(all_Indices_All_Trends$Species_number)
    monitoring_scheme<-as.data.frame(matrix(NA, ncol = ncol(scheme)-3,nrow=length(sp)))
    nombres<-c("Species_nr","Sci_name","Country_code","Monitoring_scheme","Start_year","Number_visits_spring","Census_method","Number_of_points_per_square_or_transect",
               "Transect_length","Mapping_area","Count_unit","Count_aggregation","Recorded_species","All_localities_used_every_year")
    names(monitoring_scheme)<-nombres
    monitoring_scheme$Species_nr<-sp
    monitoring_scheme$Country_code<-country_code
    for(i in 1:nrow(monitoring_scheme)){
      if(is.element(monitoring_scheme$Species_nr[i],species_list_pcbms$Euring)){
        monitoring_scheme$Sci_name[i] <- species_list_pcbms[species_list_pcbms$Euring==monitoring_scheme$Species_nr[i],2] 
      }else{
        cat(paste0("Species: ",monitoring_scheme$Species_nr[i]," is not included in the species list" ))
        stop()
      }
      
    }
    monitoring_scheme[4:14]<- scheme[scheme$Species_nr=="all",c(2:9,13:15)]
    names(scheme)
    iterations<-nrow(scheme)
    for (i in 1:iterations){
      
      if (scheme$Species_nr[i]!="all"){
        species_vector<-scheme$Species_nr[i]
        species_list <-unlist(strsplit(species_vector,"-"))
        for (j in 1:length(species_list)){
          monitoring_scheme[monitoring_scheme$Species_nr==species_list[j],4:14] <- scheme[i,c(2:9,13:15)] 
        }
      }
    }
    
    
    write.table(monitoring_scheme,paste0(outputs, "Monitoring_scheme_",country_code,".csv"),sep = ";", dec = ".",row.names = F)
    cat("Monitoring scheme table successfully created","\n")
  }else{
    cat("Monitoring scheme shall be supplied manually","\n")
  }
  
  cat("##############################################","\n")
  cat("End of Processing output","\n")
  cat("##############################################","\n")
}

#This functions prepare the .zip file to be uploaded to ONLT
preparation.onlt<-function(){
  cat("##############################################","\n")
  cat("Start of preparation ONLT","\n")
  cat("##############################################","\n") 
  ONLT<-paste0(rtrim_folder,"04_ONLT\\")
  country<-Country
  country_code<-country_code
  prefix_to_remove=paste0(prefix,"_")
  datafolder <- paste0(rtrim_folder,"03_Outputs\\")
  datafolder2<- paste0(rtrim_folder,"02_Inputs\\")
  # Create directory to gather and prepare files to upload
  dir.create(paste0(ONLT,country,"_",country_code))
  outputfolder <- paste0(ONLT,country,"_",country_code) 
  
  # Empty the prep directory
  files_to_remove <- list.files(outputfolder, full.names = TRUE)
  if (length(files_to_remove) > 0) {
    file.remove(files_to_remove)
  }
  
  # copy files to prep directory
  
  files_to_copy <- list.files(datafolder,
                              str_c("^(All_Indices_All_Trends.csv)|",
                                    "(overview.csv)|",
                                    "(Monitoring_scheme_.+\\.csv)|",
                                    "(.*arg_input_stratum.csv)|",
                                    "(.*arg_output.csv)|",
                                    "(.*indices_TT.csv)|",
                                    "(.*ocv.csv)|",
                                    "(.*summarizing_tables.txt)$"),
                              full.names = TRUE)
  files_to_copy2 <- list.files(datafolder2,
                               str_c("^(overview.csv)|",
                                     "(.*arg_input_stratum.csv)$"),
                               full.names = TRUE)
  file.copy(files_to_copy, outputfolder)
  file.copy(files_to_copy2, outputfolder)
  ?file.copy
  # remove prefix
  files_to_rename <- list.files(outputfolder, str_c("^", prefix_to_remove, ".*$"))
  files_new_name <- str_replace(files_to_rename, prefix_to_remove, "")
  file.rename(from = str_c(outputfolder, "/", files_to_rename), 
              to = str_c(outputfolder, "/", files_new_name))
  overview=read.table(paste0(rtrim_folder,"02_Inputs\\overview.csv"),sep=";",dec=".",header = T)
  
  sp_to_delete<-overview[overview$success=="no",]
  sp_to_delete<-sp_to_delete[1]
  
  for ( i in 1:nrow(sp_to_delete)){
    sp_to_delete[i,] <- str_replace(sp_to_delete[i,],"BMP_", "")
    unlink(paste0(outputfolder,"\\",sp_to_delete[i,],"_arg_input_stratum.csv"))
  }
  # !!! LAST CHANCE TO ADD Monitoring_scheme_{country_code}.csv to the directory "Country_year" before you zip the folder!!! 
  # If you forget, add this table manually to the dataset just uploaded to the ONLT.
  
  # zip the contents of the prep directory
  files_to_zip <- list.files(outputfolder, full.names = TRUE)
  zip(str_c(outputfolder, ".zip"),
      files_to_zip,
      mode = "cherry-pick")
  cat("Zip file successfully generated","\n")
  cat("##############################################","\n")
  cat("End of preparation ONLT","\n")
  cat("##############################################","\n")
}

#This function executes in proper order all the RTRIM-shell functions
run_rtrim<-function(){
  read.files()
  data.preparation(datos,weight,execute)
  # Select all files in the directory which contain arguments.
  # The files can be recognized by the pattern "arg_input_arguments".
  listSpeciesStratumCombinations <<- dir(paste0(working_directory,"//02_Inputs"), pattern = "arg_input_stratum.csv")
  # The number of files in the directory.                                                
  numberSpeciesStratumCombinations <<- length (listSpeciesStratumCombinations)
  # Make a table for an overview. 
  # The table lists which analyses were successful or not for all runs.
  overview <<- makeOverview(listSpeciesStratumCombinations)
  rtrim.strata()
  processing.output()
  preparation.onlt()
  
}

createandclean<-function(){
  folders<-c("02_Inputs","03_Outputs","04_ONLT")
  for(i in 1:3){  
    if (file.exists(folders[i]) == TRUE){ 
      unlink(folders[i], recursive=TRUE) 
      dir.create(folders[i])
    } else{
      dir.create(folders[i])
    }
  }
  cat("Folders '02_Inputs', '03_Outputs' and '04_ONLT' have been successfully created")
}
