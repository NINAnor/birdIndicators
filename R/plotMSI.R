
plotTimeSeries_MSI <- function(MSI, results_folder, plot_name, displayPlots = FALSE, savePDF = TRUE){
  
  ## Define custom colors
  plotCols <- colorspace::divergingx_hcl(palette = "Temps", n = length(unique(MSI$IndexName)))

  ## Make plot (values only)
  p1 <- ggplot(MSI, aes(x = year, y = MSI)) + 
    geom_line(aes(colour = IndexName)) +
    geom_ribbon(aes(ymin = MSI - sd_MSI, ymax = MSI + sd_MSI, fill = IndexName), alpha = 0.2) + 
    scale_x_continuous(name = "Year", breaks = min(MSI$year):max(MSI$year)) + 
    scale_color_manual(name = "Index", values = plotCols) +
    scale_fill_manual(name = "Index", values = plotCols) + 
    ylab("MSI Value (% of baseline)") + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  p2 <- ggplot(MSI, aes(x = year, y = MSI)) + 
    geom_point(aes(colour = IndexName), size = 3) + 
    geom_pointrange(aes(ymin = MSI - sd_MSI, ymax = MSI + sd_MSI, colour = IndexName)) + 
    geom_line(aes(y = Trend, colour = IndexName), linetype = "dashed") +
    geom_ribbon(aes(ymin = lower_CL_trend, ymax = upper_CL_trend, fill = IndexName), alpha = 0.2) + 
    scale_x_continuous(name = "Year", breaks = min(MSI$year):max(MSI$year)) + 
    scale_color_manual(name = "Index", values = plotCols) +
    scale_fill_manual(name = "Index", values = plotCols) + 
    ylab("MSI Value (% of baseline)") + 
    theme_bw() + 
    theme(panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  ## Optional: plot to PDF
  if(savePDF){
    pdf(paste0(results_folder, "/", plot_name, ".pdf"), width = 10, height = 5)
    print(p1)
    dev.off()
    
    pdf(paste0(results_folder, "/", plot_name, "_smooth.pdf"), width = 10, height = 5)
    print(p2)
    dev.off()
  }

  ## Optional: display
  if(displayPlots){
    print(p1)
    print(p2)
  }
  
}
