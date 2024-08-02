saturation_plot <- function(input_patristic_df, distance, patristic, distance_label = "", patristic_label = "", maxfillval = NA){
  plot_saturation =ggplot(input_patristic_df, aes_string(x=distance, y=patristic))+
          geom_point(alpha=0)+
          geom_hex(bins=50)+
          xlab(paste(distance_label, '', sep = " "))+
          ylab(paste(patristic_label, '', sep = " "))+
          theme_light()+
          scale_fill_gradientn(colors = c("#d4d4d4", rgb(0, 0.05, 0.3)), trans = "log", limits = c(NA, maxfillval),labels = function(x) format(round(x, 0), scientific = FALSE)) +
          theme(axis.text = element_text(color = "black", size = 10),
                axis.text.x = element_text(size = 10),
                axis.title.y = element_text(size = 15),
                axis.title.x = element_text(size = 15),
                axis.ticks.y = element_line(color = "black"))+
          expand_limits(y=0)+
          expand_limits(x=0)

  #ggsave(paste(plots,al,'_',tr,'_',x_label,'_', y_label,'_saturation.png', sep=''), plot=plot_saturation)
  return(plot_saturation)
}

saturation_plot_smooth <- function(input_patristic_df, distance, patristic, distance_label = "", patristic_label = ""){
  plot_saturation =ggplot(input_patristic_df, aes_string(x=distance, y=patristic))+
          geom_point(alpha=0.3)+
          geom_smooth(method = "lm")+
          xlab(paste(distance_label, '', sep = " "))+
          ylab(paste(patristic_label, '', sep = " "))+
          theme_light()+
          theme(axis.text = element_text(color = "black", size = 10),
                axis.text.x = element_text(size = 10),
                axis.title.y = element_text(size = 15),
                axis.title.x = element_text(size = 15),
                axis.ticks.y = element_line(color = "black"))+
          expand_limits(y=0)+
          expand_limits(x=0)
  #ggsave(paste(plots,al,'_',tr,'_',distance_label,'_', patristic_label,'_saturation_smooth.png', sep=''), plot=plot_saturation)
  return(plot_saturation)
}


get_before_after <- function(df, y, threshold, threshold_filtering =0){
  # IMD vs ML patristic distance
  df_imd_ml = df[,c(y, "pML", "family", "seq1", "seq2")]
  df_imd_ml_before = df_imd_ml[df_imd_ml$pML < threshold,]
  df_imd_ml_after = df_imd_ml[df_imd_ml$pML >= threshold,]
  # Check if there are enough points to fit a line
  if(nrow(df_imd_ml_before) < 10 | nrow(df_imd_ml_after) < 10){
    return( list(df = df_imd_ml, line_before = NA, line_after = NA, slope_before = NA, slope_after = NA, r2_before = NA, r2_after = NA))
  }

  # # only keep going if the number of points on the left and right are at least 49% of the total
  if(nrow(df_imd_ml_before) < nrow(df_imd_ml)*threshold_filtering | nrow(df_imd_ml_after) < nrow(df_imd_ml)*threshold_filtering){
    return( list(df = df_imd_ml, line_before = NA, line_after = NA, slope_before = NA, slope_after = NA, r2_before = NA, r2_after = NA))
  }


  line_imd_ml_before = lm(as.formula(paste(y, "~", "pML")), data = df_imd_ml_before)
  slope_imd_ml_before <- coef(line_imd_ml_before)[2]
  line_imd_ml_after = lm(as.formula(paste(y, "~", "pML")), data = df_imd_ml_after)
  slope_imd_ml_after <- coef(line_imd_ml_after)[2]
  # round
  slope_imd_ml_before = round(slope_imd_ml_before, 2)
  slope_imd_ml_after = round(slope_imd_ml_after, 2)
  # calculate R²
  r_imd_ml_before = cor(df_imd_ml_before[[y]], df_imd_ml_before$pML, method = "pearson")
  r2_imd_ml_before = round(r_imd_ml_before^2,2)
  r_imd_ml_after = cor(df_imd_ml_after[[y]], df_imd_ml_after$pML, method = "pearson")
  r2_imd_ml_after = round(r_imd_ml_after^2,2)

  return( list(df = df_imd_ml, line_before = line_imd_ml_before, line_after = line_imd_ml_after, slope_before = slope_imd_ml_before, slope_after = slope_imd_ml_after, r2_before = r2_imd_ml_before, r2_after = r2_imd_ml_after))

}


plot_distribution_across_families <- function(slopes_per_family, y){
  # plot the distribution of R² values before
  slopes_per_family$y = slopes_per_family[[y]]
  medians = aggregate(y ~ method, data = slopes_per_family, FUN = median)
  p <- ggplot(slopes_per_family, aes(x = y, fill = method, color = method)) +
              geom_density(alpha = 0.5) +
              theme_minimal() +
              xlab(y) + 
              ylab("Density") +
              theme(legend.title=element_blank())+ 
              geom_vline(data = medians, aes(xintercept = y, color = method), linetype = "dashed", linewidth = 0.9)+
              scale_fill_manual(values = palette)+
              scale_color_manual(values = palette)+
              theme(legend.title = element_blank())+
              theme(axis.text = element_text(color = "black", size = 10))
  return(p)
}


plot_boxplot_across_families <- function(slopes_per_family, y){
  correlations_melted = melt(slopes_per_family[,c(y, "method")])
  # plot the distribution of R² values before
  correlations_melted$method = factor(correlations_melted$method, levels = c("pdist", "TM ", "IMD "))
  slopes_per_family$y = slopes_per_family[[y]]
  p <- ggplot(correlations_melted, aes(x = method, y = value, fill = method, color = method)) +
              geom_boxplot(alpha = 0.5) +
              theme_light() +
              xlab("") + 
              ylab(paste(y)) +
              theme(legend.title=element_blank())+ theme(legend.position = "right")+
              scale_y_continuous(labels = function(x) format(x, scientific = FALSE))+
              scale_fill_manual(values = palette)+
              scale_color_manual(values = palette)+
              scale_y_continuous(labels = function(x) format(x, scientific = FALSE))
  return(p)
}


plot_saturation_with_fitted_lines <- function(df_imd_ml,line_imd_ml_before, line_imd_ml_after, y, ylabel,  threshold, maxfillval = NA){
  threshold_x <- threshold
  IMD_intercept <- coef(line_imd_ml_before)[1]
  IMD_slope <- coef(line_imd_ml_before)[2]
  IMD_threshold_y <- IMD_intercept + IMD_slope * threshold_x
  IMD_max_y = max(df_imd_ml[[y]])
  IMD_x_for_y_100 = (IMD_max_y - IMD_intercept)/IMD_slope

  # if IMD x_for_y_100 is greater than max x
  if(IMD_x_for_y_100 > max(df_imd_ml$pML)){
    IMD_x_for_y_100 = max(df_imd_ml$pML)
    IMD_max_y = IMD_intercept + IMD_slope * IMD_x_for_y_100
  }
  
  col_left = "#5de3f5"
  col_right = "#22bb74"
  line_size = 0.7

  IMD_intercept_after <- coef(line_imd_ml_after)[1]
  IMD_slope_after <- coef(line_imd_ml_after)[2]
  IMD_threshold_y_after <- IMD_intercept_after + IMD_slope_after * threshold_x

  pml_dimd <- saturation_plot(df_imd_ml, "pML", y, "ML patristic distance", ylabel, maxfillval) +
    geom_segment(aes(x = min(df_imd_ml$pML), y = IMD_intercept + IMD_slope * min(df_imd_ml$pML), 
                    xend = threshold_x, yend = IMD_threshold_y), color = col_left, size =line_size) +
    geom_segment(aes(x = threshold_x, y = IMD_threshold_y, 
                    xend = IMD_x_for_y_100, yend = IMD_max_y), 
                color = col_left, size =line_size)+
    geom_segment(aes(x = min(df_imd_ml$pML), y = IMD_intercept_after + IMD_slope_after * min(df_imd_ml$pML), 
                    xend = threshold_x, yend = IMD_threshold_y_after), color = col_right, size =line_size) +
    geom_segment(aes(x = threshold_x, y = IMD_threshold_y_after, 
                    xend = max(df_imd_ml$pML), yend = IMD_intercept_after + IMD_slope_after * max(df_imd_ml$pML)), 
                color = col_right, size =line_size)

  r = cor(df_imd_ml[[y]], df_imd_ml$pML, method = "pearson")
  r2 = r^2
  pml_dimd = pml_dimd + ggtitle(paste("R² = ", round(r2, 2)))
  pml_dimd = pml_dimd + theme(plot.title = element_text(hjust = 0.5))+theme(legend.position = "none")
  pml_dimd = pml_dimd + geom_vline(xintercept = threshold, color = "red", linetype = "dashed", linewidth = 0.7)
  return(pml_dimd)

}



