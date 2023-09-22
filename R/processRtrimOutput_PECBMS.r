#' Process results from PECBMS Trim analysis
#'
#' @param folder character. Path to directory which contains the RTRIM argument
#' and output files, and into which to save outputs.
#'
#' @return a dataframe containing collated results from PECBMS Trim analyses
#' for all species and years.
#' @export
#'
#' @examples

processRtrimOutput_PECBMS <- function(folder){
  
  # ----------------------------------------------------------------------------
  # START OF CODE FROM PECBMS
  #-----------------------------------------------------------------------------
  
  #####################################################################################################################
  # RTRIM-shell second script. Version RTRIM-shell_v1.7
  # Tool to process rtrim results (= TRIM in R) for multiple species and subsets of sites.
  # Marnix de Zeeuw rtrim@cbs.nl Statistics Netherlands 2019.
  #####################################################################################################################
  
  # !!! IMPORTANT!!! CHANGES YOU have to do in the script:
  # 1. Line 17: Adjust the path to the file called functions_Shell_Rtrim.R according to your computer (try different slashes, if double backslash does not work)
  # 2. Line 24: Adjust the path to your working folder according to your computer (try different slashes, if double backslash does not work)
  
  #####################################################################################################################
  # PREPARATION. 
  #####################################################################################################################

  
  #####################################################################################################################
  # SELECTING SUCCESSFUL RUNS.
  #####################################################################################################################
  # Determines which datafiles (combinations of species and stratum) have been analysed successfully. 
  # This information is found in "overview.csv", which has been produced by the script that called rtrim.
  
  overview <- read.table(paste0(folder, "/overview.csv"),sep=";",dec=".",header = T)   #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
  listsuccessfulAnalyses <- overview$ss_combinations[overview$success == "yes"]
  listSpeciesStratumCombinations <- paste(listsuccessfulAnalyses, "_arg_input_stratum.csv", sep = "")
  
  # Determine how many combinations of species and stratum have been analysed successfully.
  numberSpeciesStratumCombinations <- length(listSpeciesStratumCombinations)
  
  # File with trends and indices for each combination of species and stratum
  all_Indices_All_Trends <- make_All_Indices_All_Trends(overview = overview, listSpeciesStratumCombinations = listSpeciesStratumCombinations)
  
  #####################################################################################################################
  # PROCESSING OUTPUT OF SUCCESFUL RUNS.
  #####################################################################################################################
  
  for (j in 1:numberSpeciesStratumCombinations) {
    
    # The file with arguments contains the arguments to run the analysis for a particular combination of species and stratum (stratum is e.g. a region). 
    # These arguments are used when calling the function 'rtrim'.
    arguments <- read.table(paste0(folder, "/", listSpeciesStratumCombinations[j]), header = TRUE, stringsAsFactors = FALSE,sep=";",dec=".")
    counts <- read.table(paste0(folder, "/", arguments$File, "_counts.csv"), header = TRUE,dec=".",sep=";")   #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
    #13/01/2023 : changes to accept either comma or dot as decimal symbol  - v1.6
    if(is.character(counts$count)==T){
      counts <- read.table(paste0(folder, "/", arguments$File, "_counts.csv"), stringsAsFactors = FALSE, header = TRUE, dec = ",",sep=";")   
    } 
    #13/01/2023 End of  changes
    # Also the result of the rtrim-analysis is required.
    load(paste0(folder, "/", arguments$File, ".RData")) # produces object with name "result"
    
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
    
    if(arguments$Save_fitted_values){
      
      FI <- results(result)
      
    }
    
    #####################################################################################################################
    # WRITING OUTPUT FILES. 
    #####################################################################################################################
    # Indices and time totals.
    name_Indices_TT_file <- paste0(folder, "/", arguments$File, "_indices_TT.csv")
    write.table(indices_TT_file, name_Indices_TT_file, row.names = FALSE,sep=";",dec=".")   #13/01/2023 Javier Rivas: Changes to produce outputs with dots  - v1.6
    ############################################################ 
    # Slopes entire period and arguments.
    name_arg_output_file <- paste0(folder, "/", arguments$File, "_arg_output.csv") 
    write.table(arg_output_file, name_arg_output_file, row.names = FALSE,sep=";",dec=".")   #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
    ############################################################
    # Covariant matrix.
    name_covariant_matrix <- paste0(folder, "/", arguments$File, "_ocv.csv")
    write.table(covariant_matrix, name_covariant_matrix, row.names = FALSE,sep=";",dec=".")   #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
    ############################################################
    # File with fitted values.
    if(arguments$Save_fitted_values){
      name_Fitted_Values_File <- paste0(folder, "/", arguments$File, "_fitted_values.csv")
      write.table(FI, name_Fitted_Values_File, row.names = FALSE,sep=";",dec=".")  #13/01/2023 Javier Rivas: Changes to produce outputs with dots - v1.6
    }
    ############################################################
    # 25/01/2022 Javier Rivas: Generating Tables 1-3 for the PECBMS coordinators - v1.4
    name_rawdata=paste0(folder, "/", arguments$File, "_summarizing_tables.txt")
    # We are filling the  table on numbers of sites and counts
    write(paste0("\n", arguments$File, "\n"), file = name_rawdata)
    
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
  
  all_Indices_All_Trends <- all_Indices_All_Trends[order(all_Indices_All_Trends$Species_number, all_Indices_All_Trends$Stratum_number, all_Indices_All_Trends$Recordtype_number), ]
  write.table(all_Indices_All_Trends, paste0(folder, "/All_Indices_All_Trends.csv"), row.names = FALSE,sep=";",dec=".")           # 09/04/2021 Eva Silarova changed "all" to capital "All": "all_Indices_All_Trends.csv" was changed to "All_Indices_All_Trends.csv" - v1.3
  ############################################################                                                  # 13/01/2023 Javier Rivas : Changes to produce outputs with dots
  
  # ----------------------------------------------------------------------------
  # END OF CODE FROM PECBMS
  #-----------------------------------------------------------------------------
  
  ## Return index trend data
  return(all_Indices_All_Trends)
}

