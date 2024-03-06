#принимает на вход df, c колонками
# cb.x номер баркода  из df отсортированного для knee_plot для образца 1
# umi_per_CB.x количество UMI для этго баркода в образце 1
# dist_s1  - отношение количества UMI для этого баркода в образце 1 к количеству всех UMI для этого баркода

knee_plot_UMI_sample_dist <- function(df){
  colors <- brewer.pal(n = 11, name = "PRGn")
  knee_plot_UMI_sample_dist <- ggplot(df, aes(x = df$cb.x, y = df$umi_per_CB.x, colour = df$dist_s1)) +
    geom_point(size = 1) +
    labs(y = 'UMI count, log10', x = 'cell barcodes, log10') +
    theme_linedraw() +
    theme(axis.text=element_text(size=14),
          axis.title=element_text(size=16,face="bold")) +
    scale_x_log10() +
    scale_y_log10() +
    scale_colour_gradientn(colours = colors,
                           guide = guide_colorbar(title = "UMI in Sample1 presense"))+
    ggtitle('Solo1')
  return (knee_plot_UMI_sample_dist)
}
  
#DF_for_knee_s1_without1 <-subset(DF_for_knee_s1, dist_s1 < 1) 
#plot <- knee_plot_UMI_sample_dist(DF_for_knee_s1_without1)
#plot

