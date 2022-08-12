#' Boxplot of expression data
#' Plotting expression data of Coprinopsis cinerea
#'
#' @param gene A character of a gene name: CopciAB_ and a number. E.g., CopciAB_411205
#' @param log A logical indicating whether expression values should be log 2 transformed (Default is TRUE)
#' @import ggplot2
#' @export


copci_exp_plot <- function(gene, log = TRUE, expdata = exp_data, coldata = df_combined, order = "order",
                           facetting = "plot_facets", color = "fill"){
  df <- data.frame(expression = expdata[rownames(expdata) %in% gene,],
                   samples = coldata[,colnames(coldata) %in% order],
                   facets = coldata[,colnames(coldata) %in% facetting],
                   fill = coldata[,colnames(coldata) %in% color])
  if(log){
    ylab_text <- "Log 2 Normalized Counts"
  }else{
    ylab_text <- "Normalized Counts"
    df$expression <- 2^df$expression
  }
  ggplot(df, aes(x = samples, y = expression, fill = fill)) +
    geom_boxplot() +
    ylab(ylab_text) +
    xlab("Samples") +
    labs(fill = "Tissue types") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          axis.title = element_text(size=14, face="bold"),
          panel.grid.major.y = element_line(color = "black",
                                            size = 0.5,
                                            linetype = 1),
          legend.position = "top",
          strip.text = element_text(size=12, face = "bold"),
          strip.background = element_blank()) +
    facet_wrap(~ facets, scales = "free_x", nrow = 2, ncol = 2, drop = TRUE) +
    guides(fill=guide_legend(nrow=1,byrow=TRUE)) +
    scale_fill_brewer(palette = "Set3")
}
