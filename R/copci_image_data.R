#' Copci data image processing
#'
#' Creates a data based on Coprinopsis cinerea tissue-specific expression data that can be used for image heat map plotting.
#'
#' @param HyphalKnot Two kind of Hyphal knot samples exist. You can choose of which expressiond data should be used: H1_VM or H1_cap
#' @return A list of three elements. First element contains the image_df with three columns. The second element contains the image pointer created by image_read. The third element contains the expression data.
#' @import magick stringr
#' @export

copci_image_data <- function(HyphalKnot = "H1_VM"){
  #edit the image_df
  H1_samples <- image_df$areas[image_df$images %in% "H1"]
  to_remove <- H1_samples[!H1_samples %in% HyphalKnot]
  image_df <- image_df[!image_df$areas %in% to_remove,]
  images <- unique(image_df$images)
  path <- paste0(system.file("extdata", package = "tgexer"), "/")
  # order svg_files
  svg_files <- list.files(path)
  svg_files_name <- svg_files
  svg_files_name <- str_remove(svg_files_name, "\\.svg")
  names(svg_files) <- str_remove(svg_files_name, "_.*")
  svg_files <- svg_files[match(images, names(svg_files))]
  image_list <- list()
  for(i in 1:length(svg_files)){
    image_temp <- image_read(paste0(path, svg_files[i]))
    image_list[[i]] <- image_trim(image_temp)
  }
  names(image_list) <- images
  # edit the mean_exp_data
  mean_exp_data <- mean_exp_data[,!colnames(mean_exp_data) %in% to_remove]
  image_data <- list(image_df = image_df, image_list = image_list, mean_exp_data = mean_exp_data)
  return(image_data)
}
