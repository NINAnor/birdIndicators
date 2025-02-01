#' Write newly calculated NI indicator data to CSV for manual expert review
#'
#' @param dir character. Name of directory into which to save CSV files.
#' @param UpdatedIndicator_data a list containing updated indicator data ready
#' formatted for upload to NI database. 
#' @param sppNI table containing information on species included in NI.
#'
#' @return one .csv file per species within the specified directory containing
#' updated indicator data ready for review. 
#' @export
#'
#' @examples
writeIndicatorData_forReview <- function(dir, UpdatedIndicator_data, sppNI){
  
  ## Create directory if needed
  if(!dir.exists(file.path(dir))){
    dir.create(file.path(dir))
  }
  
  ## Save overview table into directory
  readr::write_excel_csv(sppNI, file = paste0(dir, "/NI_indicators_overview.csv"))
  
  
  for(i in 1:length(UpdatedIndicator_data)){
    
    if(!is.null(UpdatedIndicator_data[[i]])){
      
      ## Extract indicator id
      indID <- UpdatedIndicator_data[[i]]$indicatorValue$indicatorId[1]
      
      ## Write updated indicator data as csv
      readr::write_excel_csv(UpdatedIndicator_data[[i]]$indicatorValues,
                             file = paste0(dir, "/IndID_", indID, ".csv"))
      
      ## Create and save simple visualization of indicator data
      plot_data <- UpdatedIndicator_data[[i]]$indicatorValues %>%
        dplyr::filter(!is.na(verdi) & verdi!=-1 & yearName != "Referanseverdi") %>%
        dplyr::mutate(Year = as.numeric(yearName)) %>%
        dplyr::filter(Year >= 1990)
      
      plot_1 <- ggplot(plot_data, aes(x = Year, y = verdi, group = areaName)) + 
        geom_point(aes(color = areaName), size = 2) + 
        geom_errorbar(aes(ymin = nedre_Kvartil, ymax = ovre_Kvartil, color = areaName)) +
        geom_hline(yintercept = 100, color = "grey30", linetype = "dashed") +
        ylab("Value") +
        scale_color_brewer(palette = "Dark2") +
        theme_classic()
      
      plot_2 <- ggplot(plot_data, aes(x = Year, y = verdi, group = areaName)) + 
        geom_line(aes(color = areaName)) + 
        geom_ribbon(aes(ymin = nedre_Kvartil, ymax = ovre_Kvartil, fill = areaName), alpha = 0.5) +
        geom_hline(yintercept = 100, color = "grey30", linetype = "dashed") +
        ylab("Value") +
        scale_color_brewer(palette = "Dark2") +
        scale_fill_brewer(palette = "Dark2") + 
        theme_classic()
      
      pdf(paste0(dir, "/IndID_", indID, "_graphs.pdf"), width = 6, height = 5)
      gridExtra::grid.arrange(plot_1, plot_2, ncol = 1, top = UpdatedIndicator_data[[i]]$indicatorValue$indicatorName[1])
      dev.off()
    }
  }
}