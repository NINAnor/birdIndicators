#Libraries and functions
#####################################################################################################################
library(rtrim)

#####################################################################################################################
#Paths to the data
#####################################################################################################################
combine_path<-"C:\\Users\\diego.pavon-jordan\\OneDrive - NINA\\Documents\\NINA ONGOING projects\\2023_Pipeline_for_data_access\\birdIndicators\\data"
base_year<-1990
year_from<-2000

data_path<-"C:\\Users\\diego.pavon-jordan\\OneDrive - NINA\\Documents\\NINA ONGOING projects\\2023_Pipeline_for_data_access\\birdIndicators\\PECBMS_Files\\Species_files\\"


dictionary<-read.table(paste0(combine_path,"\\combine_dictionary.csv"),sep=";", dec=".",header = T)
sp_countries<-read.table(paste0(combine_path,"\\Species_Countries_combine.csv"),sep=",", dec=".",header = T)
countries_file<-read.table(paste0(combine_path,"\\Countries.csv"),sep=";", dec=".",header = T)
species_file<-read.table(paste0(combine_path,"\\Species.csv"),sep=";", dec=".",header = T,quote = "\"")
level1<-dictionary$Level1
level2<-unique(dictionary$Level2)



output<-paste0(data_path,"\\output")
dir.create(output)

#####################################################################################################################
#Reading the files
#####################################################################################################################

for (i in 1:length(level2)){
  print(paste0("Country: ",level2[i]))
  schemes<-subset(dictionary,dictionary$Level2==level2[i])
  scheme<-schemes$Level1
  Country      <- level2[i]
  #We identify which species are in the scheme
  spp_scheme<-NA
  for(j in 1:length(scheme)){
    
    indices_vector<-list.files(data_path, pattern= paste0(scheme[j],"_indices_"))
    if(j==1){
      spp_scheme<-indices_vector
    }else{
      spp_scheme<-c(spp_scheme,indices_vector) 
    }
  }
  #Now we create the species vector
  ss<- gsub("_.*", "", spp_scheme)
  spp_vector<<-unique(ss)
  
  for(j in 1:length(spp_vector)){
 
    print(paste0("Species: ",spp_vector[j]))
    data1<-read.table(paste0(data_path,spp_vector[j],"_1_",min(scheme),"_indices_TT.csv"),sep=";",dec=".",header=T)
    data2<-read.table(paste0(data_path,spp_vector[j],"_1_",max(scheme),"_indices_TT.csv"),sep=";",dec=".",header=T)
    #####################################################################################################################
    #Stablishing old and new scheme
    ##################################################################################################################### 
  
    classificator<-schemes[schemes$Level1==min(scheme),]
    
    if (classificator$class=="new"){
      newscheme<-data1
      newscheme_stratum<-min(scheme)
      oldscheme<-data2
      oldscheme_stratum<-max(scheme)
    }else{
      newscheme<-data2
      newscheme_stratum<-max(scheme)
      oldscheme<-data1
      oldscheme_stratum<-min(scheme)
    }

    not_overlaping_years_old<-oldscheme[!is.element(oldscheme$Year,newscheme$Year),"Year"]
    not_overlaping_years_new<-newscheme[!is.element(newscheme$Year,oldscheme$Year),"Year"]
    if(!is.na(not_overlaping_years_new[1])==T){
      not_overlaping_years<-c(not_overlaping_years_old,not_overlaping_years_new)
    }else{
      not_overlaping_years<-not_overlaping_years_old
    }
    overlaping_years<-oldscheme[is.element(oldscheme$Year,newscheme$Year),"Year"]
    if(level2[i]==17){
      old_scheme<-oldscheme[,c("Year","TT_imputed")]
      names(old_scheme)<-c("time","count")
      old_scheme$site<-oldscheme_stratum
      
      new_scheme<-newscheme[,c("Year","TT_imputed")]
      names(new_scheme)<-c("time","count")
      new_scheme$site<-newscheme_stratum
      
      
      sp_country_subset<-as.data.frame(subset(sp_countries,sp_countries$Species_nr==spp_vector[j] & sp_countries$Country_code== min(scheme)))
      names(sp_country_subset)<-c("Species_nr","Country_code","Country","SciName","Population_size_mean","Population_size_min","Population_size_max","Weight","Weight_Year_First","Weight_Year_Last","Year_first","Year_last")
      weight_year<-sp_country_subset[,c("Country_code","Population_size_mean","Weight_Year_First","Weight_Year_Last")]
      mean_pop<-weight_year$Population_size_mean
      TT_subset<-subset(old_scheme,old_scheme$time >= weight_year$Weight_Year_First & old_scheme$time <= weight_year$Weight_Year_Last)
      mean_TT<-mean(TT_subset$count)
      weight<-((mean_pop/2)/mean_TT)
      old_scheme$weight<-weight
    
      sp_country_subset<-as.data.frame(subset(sp_countries,sp_countries$Species_nr==spp_vector[j] & sp_countries$Country_code== min(scheme)))
      names(sp_country_subset)<-c("Species_nr","Country_code","Country","SciName","Population_size_mean","Population_size_min","Population_size_max","Weight","Weight_Year_First","Weight_Year_Last","Year_first","Year_last")
      weight_year<-sp_country_subset[,c("Country_code","Population_size_mean","Weight_Year_First","Weight_Year_Last")]
      mean_pop<-weight_year$Population_size_mean
      TT_subset<-subset(new_scheme,new_scheme$time >= weight_year$Weight_Year_First & new_scheme$time <= weight_year$Weight_Year_Last)
      mean_TT<-mean(TT_subset$count)
      weight<-((mean_pop/2)/mean_TT)
      new_scheme$weight<-weight
      
      TimetotalsAllStrata<<-rbind(old_scheme,new_scheme)
      #################################################################################################################
      #Gathering ocv's
      #################################################################################################################
      icv = list() # defines file in list format
      
      
      
      for (st in 1:length(scheme)) {
        filenameCovmatrix <- paste(data_path, spp_vector[j], "_1_", scheme[st], "_ocv.csv", sep="")   # pool file input, with schedule in name 
        icv_1site = unname(as.matrix(read.table(filenameCovmatrix,header=TRUE,sep=";",dec="."))) # unname is required to delete colnames and rownames
        
        icv[[st]] = icv_1site # combine ocv's in a list needed for RTRIM
      }
      ###################
      
      weight <- TimetotalsAllStrata$weight
      
      result <- tryCatch(trim(count ~ site + time , data=TimetotalsAllStrata, weights = "weight", model = 2, changepoints = "all", autodelete = TRUE,covin=icv ),error=warning)
      slopes_subperiod_imputed <- overall(result, which = "imputed", changepoints = year_from)  
      slopes_subperiod_fitted <- overall(result, which = "fitted", changepoints = year_from)
    }else{
      ######################################################################################################################
      #We need to define which schedule  is within the weighting period since it will be used to estimate
      #the species population for the one out of the weighting period
      sp_country_subset_min<-as.data.frame(subset(sp_countries,sp_countries$Species_nr==spp_vector[j] & sp_countries$Country_code== min(scheme)))
      
      if(is.element(sp_country_subset_min$Weight_Year_First,newscheme$Year)==T & is.element(sp_country_subset_min$Weight_Year_Last ,newscheme$Year)==T){
        weighting_scheme<-"new"
        non_weighting<-"old"
      }else{
        weighting_scheme<-"old"
        non_weighting<-"new"
      }
      
      
      ######################################################################################################################
      #First we will run RSWAN like trim model
      ##################################################################################################################### 
      #Let's format the data
      ##################################################################################################################### 
      names(oldscheme)
      old_scheme<-oldscheme[,c("Year","TT_imputed")]
      names(old_scheme)<-c("time","count")
      old_scheme$site<-oldscheme_stratum
      
      new_scheme<-newscheme[,c("Year","TT_imputed")]
      names(new_scheme)<-c("time","count")
      new_scheme$site<-newscheme_stratum
      
      
      ##################################################################################################################### 
      #We need to calculate the weight for the weighting scheme
      ##################################################################################################################### 
      if (weighting_scheme == "old"){
        sp_country_subset<-as.data.frame(subset(sp_countries,sp_countries$Species_nr==spp_vector[j] & sp_countries$Country_code== min(scheme)))
        names(sp_country_subset)<-c("Species_nr","Country_code","Country","SciName","Population_size_mean","Population_size_min","Population_size_max","Weight","Weight_Year_First","Weight_Year_Last","Year_first","Year_last")
        weight_year<-sp_country_subset[,c("Country_code","Population_size_mean","Weight_Year_First","Weight_Year_Last")]
        mean_pop<-weight_year$Population_size_mean
        TT_subset<-subset(old_scheme,old_scheme$time >= weight_year$Weight_Year_First & old_scheme$time <= weight_year$Weight_Year_Last)
        mean_TT<-mean(TT_subset$count)
        
        weight<-(mean_pop/mean_TT)
        old_scheme_first<-old_scheme
        old_scheme_first$weight<-weight
        first_scheme<-old_scheme_first
        ##################################################################################################################### 
        #Loading .ocv
        #####################################################################################################################      
        icv = list() # defines file in list format
        
        filenameCovmatrix <- paste(data_path, spp_vector[j], "_1_", oldscheme_stratum, "_ocv.csv", sep="")   # pool file input, with schedule in name 
        icv_1site = unname(as.matrix(read.table(filenameCovmatrix,header=TRUE,sep=";",dec="."))) # unname is required to delete colnames and rownames
        icv[[1]] = icv_1site # combine ocv's in a list needed for RTRIM
        
      }
      ##################################################################################################################### 
      #We need to calculate the weight for the weighting scheme
      ##################################################################################################################### 
      if(weighting_scheme=="new"){
        sp_country_subset<-as.data.frame(subset(sp_countries,sp_countries$Species_nr==spp_vector[j] & sp_countries$Country_code== min(scheme)))
        names(sp_country_subset)<-c("Species_nr","Country_code","Country","SciName","Population_size_mean","Population_size_min","Population_size_max","Weight","Weight_Year_First","Weight_Year_Last","Year_first","Year_last")
        weight_year<-sp_country_subset[,c("Country_code","Population_size_mean","Weight_Year_First","Weight_Year_Last")]
        mean_pop<-weight_year$Population_size_mean
        
        TT_subset<-subset(new_scheme,new_scheme$time >= weight_year$Weight_Year_First & new_scheme$time <= weight_year$Weight_Year_Last)
        mean_TT<-mean(TT_subset$count)
        weight<-(mean_pop/mean_TT)
        new_scheme_first<-new_scheme
        new_scheme_first$weight<-weight
        first_scheme<-new_scheme_first
        ##################################################################################################################### 
        #Loading .ocv
        #####################################################################################################################   
        icv = list() # defines file in list format
        
        filenameCovmatrix <- paste(data_path, spp_vector[j], "_1_", newscheme_stratum, "_ocv.csv", sep="")   # pool file input, with schedule in name 
        icv_1site = unname(as.matrix(read.table(filenameCovmatrix,header=TRUE,sep=";",dec="."))) # unname is required to delete colnames and rownames
        icv[[1]] = icv_1site # combine ocv's in a list needed for RTRIM
      }
      ##################################################################################################################### 
      #We run the first trim model (RSWAN like), to estimate population on the non weighting scheme
      #####################################################################################################################  
      TimetotalsAllStrata<-first_scheme
      weight <- TimetotalsAllStrata$weight
      
      result1 <- tryCatch(trim(count ~ site + time , data=TimetotalsAllStrata, weights = "weight", model = 2, changepoints = "all", autodelete = TRUE,covin=icv ),error=warning)
      
      ##################################################################################################################### 
      #Now is time to select which years are overlapping and will be used to estimate the populations
      #####################################################################################################################    
      tt_first <- totals(result1, which = "both")
      if(length(overlaping_years)==1){
        xx<-tt_first[tt_first$time==overlaping_years,]
        min_ov_year<-overlaping_years[1]
        max_ov_year<-overlaping_years[1]
        mean_tt<-mean(xx$fitted)
      }
      if(length(overlaping_years)>1 & length(overlaping_years)<=5){
        xx<-tt_first[tt_first$time >= min(overlaping_years) & tt_first$time <= max(overlaping_years) ,]
        min_ov_year<-min(overlaping_years)
        max_ov_year<-max(overlaping_years)
        mean_tt<-mean(xx$fitted)
      }
      if (length(overlaping_years)==6){
        xx<-tt_first[tt_first$time >= (max(overlaping_years)-5) & tt_first$time <= (max(overlaping_years)-1) ,]
        min_ov_year<-(min(overlaping_years)+1)
        max_ov_year<-(min(overlaping_years)+5)
        mean_tt<-mean(xx$fitted)
      }
    
    if (length(overlaping_years)==7){
      xx<-tt_first[tt_first$time >= (max(overlaping_years)-6) & tt_first$time <= (max(overlaping_years)-2) ,]
      min_ov_year<-(min(overlaping_years)+1)
      max_ov_year<-(min(overlaping_years)+5)
      mean_tt<-mean(xx$fitted)
    }
    if (length(overlaping_years)>=8){
      xx<-tt_first[tt_first$time >= (max(overlaping_years)-7) & tt_first$time <= (max(overlaping_years)-3) ,]
      min_ov_year<-(min(overlaping_years)+1)
      max_ov_year<-(min(overlaping_years)+5)
      mean_tt<-mean(xx$fitted)
    }
    
    if (weighting_scheme == "old"){
      sp_countries[sp_countries$Species_nr==spp_vector[j] & sp_countries$Country_code==newscheme_stratum, "Weight_Year_First"]<-min_ov_year
      sp_countries[sp_countries$Species_nr==spp_vector[j] & sp_countries$Country_code==newscheme_stratum, "Weight_Year_Last"]<-max_ov_year
      sp_countries[sp_countries$Species_nr==spp_vector[j] & sp_countries$Country_code==newscheme_stratum, "Population.size..geomean."]<-mean_tt
    }
    if (weighting_scheme == "new"){
      sp_countries[sp_countries$Species_nr==spp_vector[j] & sp_countries$Country_code==oldscheme_stratum, "Weight_Year_First"]<-min_ov_year
      sp_countries[sp_countries$Species_nr==spp_vector[j] & sp_countries$Country_code==oldscheme_stratum, "Weight_Year_Last"]<-max_ov_year
      sp_countries[sp_countries$Species_nr==spp_vector[j] & sp_countries$Country_code==oldscheme_stratum, "Population.size..geomean."]<-mean_tt
    }
    
    ##################################################################################################################### 
    #We recalculate the wheights for both schedules and format the data for the second RSWAN like trim model
    #####################################################################################################################
    #################################################################################################################
    #Old scheme
    #################################################################################################################
    sp_country_subset_old<-as.data.frame(subset(sp_countries,sp_countries$Species_nr==spp_vector[j] & sp_countries$Country_code== oldscheme_stratum))
    names(sp_country_subset_old)<-c("Species_nr","Country_code","Country","SciName","Population_size_mean","Population_size_min","Population_size_max","Weight","Weight_Year_First","Weight_Year_Last","Year_first","Year_last")
    weight_year_old<-sp_country_subset_old[,c("Country_code","Population_size_mean","Weight_Year_First","Weight_Year_Last")]
    mean_pop_old<-weight_year_old$Population_size_mean
    TT_subset_old<-subset(old_scheme,old_scheme$time >= weight_year_old$Weight_Year_First & old_scheme$time <= weight_year_old$Weight_Year_Last)
    mean_TT_old<-mean(TT_subset_old$count)
    weight_old<-(mean_pop_old/2)/mean_TT_old
    
    old_scheme$weight<-weight_old
    #################################################################################################################
    #New scheme
    #################################################################################################################
    sp_country_subset_new<-as.data.frame(subset(sp_countries,sp_countries$Species_nr==spp_vector[j] & sp_countries$Country_code== newscheme_stratum))
    names(sp_country_subset_new)<-c("Species_nr","Country_code","Country","SciName","Population_size_mean","Population_size_min","Population_size_max","Weight","Weight_Year_First","Weight_Year_Last","Year_first","Year_last")
    weight_year_new<-sp_country_subset_new[,c("Country_code","Population_size_mean","Weight_Year_First","Weight_Year_Last")]
    mean_pop_new<-weight_year_new$Population_size_mean
    TT_subset_new<-subset(new_scheme,new_scheme$time >= weight_year_new$Weight_Year_First & new_scheme$time <= weight_year_new$Weight_Year_Last)
    mean_TT_new<-mean(TT_subset_new$count)
    weight_new<-(mean_pop_new/2)/mean_TT_new
    
    new_scheme$weight<-weight_new
    
    TimetotalsAllStrata<<-rbind(old_scheme,new_scheme)
    #################################################################################################################
    #Gathering ocv's
    #################################################################################################################
    icv = list() # defines file in list format
    
    
    
    for (st in 1:length(scheme)) {
      filenameCovmatrix <- paste(data_path, spp_vector[j], "_1_", scheme[st], "_ocv.csv", sep="")   # pool file input, with schedule in name 
      icv_1site = unname(as.matrix(read.table(filenameCovmatrix,header=TRUE,sep=";",dec="."))) # unname is required to delete colnames and rownames
      
      icv[[st]] = icv_1site # combine ocv's in a list needed for RTRIM
    }
    ###################
    
    weight <- TimetotalsAllStrata$weight
    
    result <- tryCatch(trim(count ~ site + time , data=TimetotalsAllStrata, weights = "weight", model = 2, changepoints = "all", autodelete = TRUE,covin=icv ),error=warning)
    }
    
    slopes_subperiod_imputed <- overall(result, which = "imputed", changepoints = year_from)
    slopes_subperiod_fitted <- overall(result, which = "fitted", changepoints = year_from)
    
    # Generating index_TT file:
    tt <- totals(result, which = "both")
    if(min(TimetotalsAllStrata$time)< base_year){        #Javi 21.09.2021
      indices <- index(result, which = "both",base= base_year)
    }else{
      indices <- index(result, which = "both")
    }  
    indices$time <- NULL                     
    tt_indices <- cbind(tt, indices)         
    
    # RTRIM requires terms site and time, and colnames should be renamed again as Year en TT_imputed
    colnames(tt_indices) <- c("Year", "TT_model", "TT_model_SE", "TT_imputed", "TT_imputed_SE", "Index_model", "Index_model_SE", "Index_imputed", "Index_imputed_SE")
    
    # Generating covariance matrix
    covmatrix <- vcov(result)
    
    
    #Generating  empty arg_output
    arg_output_file <<- data.frame(N_sites = integer(1),        # Number of unique sites.
                                   N_time_values = integer(1),                      # Number of unique years.
                                   N_observed_zero_counts = integer(1),             # Number of zero counts.
                                   N_observed_positive_counts = integer(1),         # Number of positive counts.
                                   N_missing_counts = integer(1),                   # Number of missing counts.
                                   N_counts_total = integer(1),                     # Total number of counts.
                                   Base_year_first_year = numeric(1),
                                   Base_year_last_year = numeric(1),
                                   # Calendar year used as base year for indices. 
                                   # If Base_year_first_year equals Base_year_last_year a single year is used as base year.
                                   # If Base_year_first_year < Base_year_last_year, a period is used as base time.
                                   # In the latter case, Base_year_first_year is the first year of the period.
                                   
                                   Changepoints = numeric(1),                        # Changepoints used.
                                   Overdispersion = numeric(1),                      # Estimated overdispersion.
                                   Serial_correlation = numeric(1),                  # Estimated serial correlation.
                                   Slope_imputed_mul = numeric(1),                   # Multiplicative imputed slope for the entire period.
                                   Slope_imputed_mul_SE = numeric(1),                # Standard error of multiplicative imputed slope for the entire period.
                                   Slope_imputed_classification = character(1),      # Trend classification of multiplicative imputed slope for the entire period.
                                   Year_from = numeric(1),                           # First year of the subperiod from which a slope has been calculated. 
                                   # Last year of the subperiod is equal to the last year of the entire period.
                                   # Note that it is assumed that the first year of the subperiod is closer to the present than the first year of the entire period.
                                   
                                   Slope_from_imputed_mul = numeric(1),              # Multiplicative imputed slope for the subperiod.
                                   Slope_from_imputed_mul_SE = numeric(1),           # Standard error of multiplicative imputed slope for the subperiod.
                                   Slope_from_imputed_classification = character(1)  # Trend classification of multiplicative imputed slope for the subperiod.
                                   
    )
    #Filling Arg_output
    
    arg_output_file$N_sites                      <- result$nsite               # Number of unique sites in the counts file used. 
    arg_output_file$N_time_values                <- result$ntime               # Number of unique years in the counts file used.
    counts                                       <- tt_indices[, "TT_imputed"]
    arg_output_file$N_observed_zero_counts       <- sum(counts == 0)           # Number of zero counts in the entire period, in the counts file used.
    arg_output_file$N_observed_positive_counts   <- sum(counts > 0)            # Number of positive counts in the entire period, in the counts file used.
    arg_output_file$N_missing_counts             <- sum(is.na(counts))         # Number of missing counts in the entire period, in the counts file used.
    arg_output_file$N_counts_total               <- sum(counts > 0)              # Total number of counts in the entire period, in the counts file used.
    arg_output_file$Base_year_first_year         <- base_year  # for case of accidental more occurences of 100% Arco 16-2-2021
    arg_output_file$Base_year_last_year          <- base_year # for case of accidental more occurences of 100% Arco 16-2-2021
    arg_output_file$Changepoints                 <- paste(result$changepoints, collapse = "-")    # Changepoints used to calculate the indices. #29/03/2023 Javi Change points shall have - instead , as separator, function modified accordingly.
    
    slopes_imputed <- overall(result, which = "imputed")
    arg_output_file$Slope_imputed_mul            <- slopes_imputed$slope$mul                       # Multiplicative imputed slope over the entire period.
    arg_output_file$Slope_imputed_mul_SE         <- slopes_imputed$slope$se_mul                    # Standard error of multiplicative imputed slope over the entire period.
    arg_output_file$Slope_imputed_classification <- slopes_imputed$slope$meaning                   # Trend classification of multiplicative imputed slope over the entire period.
    arg_output_file$Slope_from_imputed_mul <- slopes_subperiod_imputed$slope$mul[slopes_subperiod_fitted$slope$from == year_from]
    arg_output_file$Slope_from_imputed_mul_SE <- slopes_subperiod_imputed$slope$se_mul[slopes_subperiod_fitted$slope$from == year_from]
    arg_output_file$Slope_from_imputed_classification <- slopes_subperiod_imputed$slope$meaning[slopes_subperiod_imputed$slope$from == year_from]
    #Output file names
    filenameArg_out <- paste(output,"\\", spp_vector[j], "_",  1, "_",  level2[i], "_arg_output.csv", sep="")   
    filenameTimeTotals_out <- paste(output,"\\",spp_vector[j], "_",  1, "_",  level2[i], "_indices_TT.csv", sep="")
    filenameOcv_out <- paste(output,"\\", spp_vector[j], "_",  1, "_",  level2[i], "_ocv.csv", sep="")   
    
    #Saving output files
    write.table(tt_indices, file=filenameTimeTotals_out, dec=".",sep=";",row.names=FALSE)
    write.table(arg_output_file, filenameArg_out, dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
    write.table(covmatrix,file=filenameOcv_out, dec=".",sep=";",row.names=FALSE)
  }
  

  
  ###########################################################################
  # assess number of records in results file                                # 
  # = combinations of species * countries in active schedule                #
  ###########################################################################
  
   
  total_number_records    <- length(spp_vector)
  lowest_year_of_schedule <-min(TimetotalsAllStrata$time)
  highest_year_of_schedule<-max(TimetotalsAllStrata$time)
  uyear <- lowest_year_of_schedule:highest_year_of_schedule
  nyear <- length(uyear)
  
  Swan_ind <- array(NA,dim=c(total_number_records, (13+nyear)))  
  colnames(Swan_ind) <- c( "Year_of_analysis", "Species_number", "Recordtype_number", "Recordtype_name", "N_sites", "Slope_imputed_mul", "Slope_imputed_mul_SE", "Slope_imputed_classification", "Slope_from_imputed_mul", "Slope_from_imputed_mul_SE", "Slope_from_imputed_classification", "Year_from", "Date_analysis", uyear) 
  Swan_ind = data.frame(Swan_ind, check.names = FALSE) # check.names = FALSE prevents setting X before calendar year as colnames
  Swan_ind_se    <- Swan_ind
  Swan_tt        <- Swan_ind
  Swan_tt_se     <- Swan_ind
  
  ###########################################################################
  # loops species * countries                                               #
  ###########################################################################
  # Note: loops of species and countries here instead of loop of records from Swan_schedule; else the same species & countrynames need to be handeld multiple times
  
  species_stratum_rec <- 1  # record counter
  

  
  for (s in 1:length(spp_vector)) { 
    
    Species      <- spp_vector[s]
    Species_name <- as.character(species_file[species_file$Species_nr == Species, "Sci_name"] )

    

    Country_name <- as.character(countries_file[countries_file$Country_code == Country, "Country"] ) 
      
    RTRIM_timetotals   <- paste(output,"\\",Species, "_1_", Country, "_indices_TT.csv",  sep="") 
    TimetotalsStratum  <- read.table(RTRIM_timetotals, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol 
    RTRIM_arguments    <- paste(output,"\\",Species, "_1_", Country, "_arg_output.csv",  sep="") 
    TimeargStratum     <- read.table(RTRIM_arguments, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")#13/02/2023 Javi: change "," for "." as decimal symbol 
    nyear_stratum      <- TimeargStratum[, "N_time_values"]
    nplot_stratum      <- TimeargStratum[, "N_sites"]
    slope_imputed_mul            <- TimeargStratum[, "Slope_imputed_mul"]
    slope_imputed_mul_se         <- TimeargStratum[, "Slope_imputed_mul_SE"]
    slope_imputed_classification <- TimeargStratum[, "Slope_imputed_classification"]
      
    #     Slope_from ignored as Year_from varies in the files which countries have provided
      
    ###########################################################################
    #     fill in tables                                                      #
    ###########################################################################
      
    #     fill in table with indices

    Swan_ind[species_stratum_rec, "Year_of_analysis"]           <- max(TimetotalsStratum[, "Year"])  # Year_of_analysis = most recent year with data
    Swan_ind[species_stratum_rec, "Species_number"]             <- Species
    Swan_ind[species_stratum_rec, "Recordtype_number"]          <- 1  
    Swan_ind[species_stratum_rec, "Recordtype_name"]            <- "indices"
    Swan_ind[species_stratum_rec, "N_sites"]                    <- nplot_stratum
      
    #     fill in table with indices se's

    Swan_ind_se[species_stratum_rec, "Year_of_analysis"]        <- max(TimetotalsStratum[, "Year"])
    Swan_ind_se[species_stratum_rec, "Species_number"]          <- Species
    Swan_ind_se[species_stratum_rec, "Recordtype_number"]       <- 2   
    Swan_ind_se[species_stratum_rec, "Recordtype_name"]         <- "se_indices"
    Swan_ind_se[species_stratum_rec, "N_sites"]                 <- nplot_stratum
      
    #     fill in table with time totals

    Swan_tt[species_stratum_rec, "Year_of_analysis"]            <- max(TimetotalsStratum[, "Year"])
    Swan_tt[species_stratum_rec, "Species_number"]              <- Species
    Swan_tt[species_stratum_rec, "Recordtype_number"]           <- 3
    Swan_tt[species_stratum_rec, "Recordtype_name"]             <- "time_totals"
    Swan_tt[species_stratum_rec, "N_sites"]                     <- nplot_stratum
      
    #     fill in table with time totals se's

    Swan_tt_se[species_stratum_rec, "Year_of_analysis"]         <- max(TimetotalsStratum[, "Year"])
    Swan_tt_se[species_stratum_rec, "Species_number"]           <- Species
    Swan_tt_se[species_stratum_rec, "Recordtype_number"]        <- 4
    Swan_tt_se[species_stratum_rec, "Recordtype_name"]          <- "se_time_totals"
    Swan_tt_se[species_stratum_rec, "N_sites"]                  <- nplot_stratum
      
    #     select the proper column of calendar years taking into account that countries have different years in their data
    for (j in 1:nyear_stratum) { 
      target_year <- TimetotalsStratum[j, "Year"]  
      col_year    <- which(colnames(Swan_ind) == target_year) # column is similar in all Swan output arrays
      Swan_ind[species_stratum_rec, col_year]    <- round(TimetotalsStratum[j, "Index_imputed"] * 100)
      Swan_ind_se[species_stratum_rec, col_year] <- round(TimetotalsStratum[j, "Index_imputed_SE"] * 100)
      Swan_tt[species_stratum_rec, col_year]     <- round(TimetotalsStratum[j, "TT_imputed"])
      Swan_tt_se[species_stratum_rec, col_year]  <- round(TimetotalsStratum[j, "TT_imputed_SE"])
    }  
      
    #     overall slope entire period; possibly only multiplicative slope available in country output 
    Swan_ind[species_stratum_rec, "Slope_imputed_mul"]            <- slope_imputed_mul
    Swan_ind[species_stratum_rec, "Slope_imputed_mul_SE"]         <- slope_imputed_mul_se
    Swan_ind[species_stratum_rec, "Slope_imputed_classification"] <- as.character(slope_imputed_classification)
    
    #     Slope_from ignored as Year_from varies between countries
    Swan_ind[species_stratum_rec, "Slope_from_imputed_mul"]             <- TimeargStratum[, "Slope_from_imputed_mul"]
    Swan_ind[species_stratum_rec, "Slope_from_imputed_mul_SE"]          <- TimeargStratum[, "Slope_from_imputed_mul_SE"]
    Swan_ind[species_stratum_rec, "Slope_from_imputed_classification"]  <- TimeargStratum[, "Slope_from_imputed_classification"]
    Swan_ind[species_stratum_rec, "Year_from"]                          <- year_from
      
    species_stratum_rec <- species_stratum_rec + 1 
    
  } # end loop species
  
  cat("All species*country files processed", "\n")
  
  # write output in horizontal format in working_folder, not in output folder
  Swan_total                    <- rbind(Swan_ind, Swan_ind_se, Swan_tt, Swan_tt_se) 
  Swan_total[, "Date_analysis"] <- Sys.Date()
  Swan_total                    <- Swan_total[order(Swan_total[, "Species_number"], decreasing=FALSE),] 
  
  Swan_countries_name           <- paste(output,"\\All_Indices_All_Trends.csv", sep="")  
  if(length(level2)!=1){
    Swan_countries_name           <- paste(output,"\\",level2[i],"_All_Indices_All_Trends.csv", sep="")  
  }
  
  write.table(Swan_total, file = Swan_countries_name, dec=".",sep=";",row.names=FALSE) #13/02/2023 Decimal symbol
  
  trimind_countries <- NULL
  trimtot_countries <- NULL  
  for(s in 1:length(spp_vector)){
    Species<-spp_vector[s]
    Species_name <- as.character(species_file[species_file$Species_nr == Species, "Sci_name"] )
    Country_name <- as.character(countries_file[countries_file$Country_code == Country, "Country"] ) 

    RTRIM_timetotals   <- paste(output,"\\",Species, "_1_", Country, "_indices_TT.csv",  sep="") 
    TimetotalsStratum  <- read.table(RTRIM_timetotals, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")
    RTRIM_arguments    <- paste(output,"\\",Species, "_1_", Country, "_arg_output.csv",  sep="") 
    TimeargStratum     <- read.table(RTRIM_arguments, header=T,check.names = FALSE,dec=".",sep=";", comment.char = "",quote = "\"")
    
    TimetotalsStratum$Code          <- Species
    TimetotalsStratum$Species       <- Species_name
    TimetotalsStratum$StratumType   <- 1
    TimetotalsStratum$Stratum       <- Country
    TimetotalsStratum$CountryGroup  <- Country_name
    
    TimetotalsStratum$Index         <- TimetotalsStratum$Index_imputed
    TimetotalsStratum$Index_SE      <- TimetotalsStratum$Index_imputed_SE
    trimind_part                    <- subset(TimetotalsStratum, select = c(Code, Species, StratumType, Stratum, CountryGroup, Year, Index, Index_SE))
    trimind_countries               <- rbind(trimind_countries, trimind_part)
    
    TimetotalsStratum$TT_Model      <- TimetotalsStratum$TT_model
    TimetotalsStratum$TT_Model_SE   <- TimetotalsStratum$TT_model_SE
    TimetotalsStratum$TT_Imputed    <- TimetotalsStratum$TT_imputed
    TimetotalsStratum$TT_Imputed_SE <- TimetotalsStratum$TT_imputed_SE
    
    trimtot_part                    <- subset(TimetotalsStratum, select = c(Code, Species, StratumType, Stratum, CountryGroup, Year, TT_Model, TT_Model_SE, TT_Imputed, TT_Imputed_SE))
    trimtot_countries               <- rbind(trimtot_countries, trimtot_part) 
  }
  trimind_countries             <- subset(trimind_countries, !duplicated(trimind_countries))   
  trimind_countries             <- subset(trimind_countries, !duplicated(trimind_countries))
  
  trimind_countries_name           <- paste(output,"\\Trimind.csv", sep="")
  trimtot_countries_name           <- paste(output,"\\Trimtot.csv", sep="") 
  if(length(level2)!=1){
      
    trimind_countries_name           <- paste(output,"\\",level2[i],"Trimind.csv", sep="")
    trimtot_countries_name           <- paste(output,"\\",level2[i],"Trimtot.csv", sep="") 
  }
  write.table(trimind_countries, file = trimind_countries_name, dec=".",sep=";",row.names=FALSE)
  write.table(trimtot_countries, file = trimtot_countries_name, dec=".",sep=";",row.names=FALSE)
}


