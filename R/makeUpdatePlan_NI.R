#' Match up correct Trim analysis outputs for each species-area in the NI database
#'
#' @param sppNI a table containing information on species included in NI.
#' @param OldIndicator_data a list containing the NI indicator data currently 
#' stored in the NI database for each species in `sppNI`. 
#'
#' @returns a list containing one table per species detailing which
#' Trim analysis output (area) should be used to update NI indicator values for
#' each area in the NI database. 
#' @export
#'
#' @examples

makeUpdatePlan_NI <- function(sppNI, OldIndicator_data){
  
  updatePlan <- list()
  
  for(i in 1:nrow(sppNI)){
    
    ## Check areas used in NI database
    NI_areas <- unique(OldIndicator_data[[i]]$indicatorValues$areaName)
    n_NI_areas <- length(NI_areas)

    ## Check whether TRIM data have been divided into sub-areas
    Trim_areaDivide <- ifelse(sppNI$dataUse_direct_N[i] | sppNI$dataUse_expert_N[i], FALSE, TRUE)
    
    ## List relevant Trim sub-areas
    if(Trim_areaDivide){
      Trim_subAreas <- ifelse(n_NI_areas > 1, c("Sør-Norge", "Nord-Norge"), "Sør-Norge")
    }else{
      Trim_subAreas <- ifelse(sppNI$dataUse_direct_SN[i] & !sppNI$dataUse_direct_NN[i], "Sør-Norge", "Norge")
    }
    
    ## Determine correct matches for NI and Trim areas
    Ind_update <- tibble::tibble(NI_areas = NI_areas) %>%
      dplyr::mutate(Trim_areaMatch = dplyr::case_when(
        !Trim_areaDivide ~ Trim_subAreas,
        NI_areas %in% c("Nord-Norge", "Nordnorge") ~ "Nord-Norge",
        TRUE ~ "Sør-Norge"
      ))

    updatePlan[[i]] <- Ind_update
  }
  
  names(updatePlan) <- sppNI$indicatorId
  return(updatePlan)
}  
  

  