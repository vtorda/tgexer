#' Heatmap plot
#'
#' Heatmap plot of the expression data through the developmental stage of Coprinopsis cinerea
#'
#' @param gene A character of a gene name: CopciAB_ and a number. E.g., CopciAB_411205
#' @param image_data A list which is the output of the copci_image_data function
#' @param norm_val NULL or a character vector length of 2. If NULL the expression values will be normalized and centered between 0 and 1 based on the values of the selected gene. If it is a vector then the minimum and the maximum values should be added which will be used to normalize and center the data.
#' @param color_n Integer. Number of colors to pick from the color palette.
#' @param colfun Function. A color palette function.
#' @import magick
#' @import RColorBrewer
#' @export
#'
copci_dev_heatmap <- function(gene, image_data,  norm_val = NULL, color_n = 75,
                              colfun = colorRampPalette(brewer.pal(8, "YlOrRd"))){

  all_val <- image_data$mean_exp_data[rownames(mean_exp_data) %in% gene,]
  all_val_stages <- sapply(str_split(names(all_val), "_"), function(x) x[1])
  # calculate color palette
  if(is.null(norm_val)){
    all_val <- (all_val - min(all_val)) / (max(all_val) - min(all_val))
  }else{
    all_val <- (all_val - norm_val[1]) / (norm_val[2] - norm_val[1])
  }
  colors_v <- colfun(color_n)
  bins <- max(all_val) / color_n
  color_df <- data.frame(color = colors_v, value_borders = 0:(color_n - 1) * bins, stringsAsFactors = FALSE)
  plot_list <- list()
  stages <- unique(image_data$image_df$images)
  for(k in 1:length(stages)){
    stage <- stages[k]
    tissue_values <- all_val[all_val_stages %in% stage]
    tissue_colors <- tissue_values
    tissue_colors <- color_df$color[sapply(tissue_colors, function(x) max(which(x >= color_df$value_borders)))]
    names(tissue_colors) <- names(tissue_values)
    image_temp <- image_data$image_list[[stage]]
    for(j in 1:length(tissue_colors)){
      tissue <- names(tissue_colors)[j]
      col <- tissue_colors[j]
      coord <- image_data$image_df$coords[str_detect(image_data$image_df$areas, tissue)]
      if(length(coord) == 1){
        image_temp <- image_fill(image_temp, col, point = coord, fuzz = 20)
      }else{
        for(i in 1:length(coord)){
          image_temp <- image_fill(image_temp, col, point = coord[i], fuzz = 20)
        }
      }
    }
    plot_list[[k]] <- image_temp
  }
  img <- c(image_scale(image_rotate(image_trim(plot_list[[1]]), 180), "150"),
           image_scale(image_rotate(image_trim(plot_list[[2]]), 180), "200"),
           image_scale(image_rotate(image_trim(plot_list[[3]]), 180), "350"),
           image_scale(image_rotate(image_trim(plot_list[[4]]), 180), "x450"),
           image_scale(image_rotate(image_trim(plot_list[[5]]), 180), "x750"),
           image_scale(image_rotate(image_trim(plot_list[[6]]), 180), "x900"),
           image_scale(image_rotate(image_trim(plot_list[[7]]), 180), "x1200"))
  merged <- image_append(img)
  image_scale(image_flop(image_rotate(merged, 180)), "1200")
}
