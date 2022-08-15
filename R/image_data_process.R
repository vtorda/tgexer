#' Image data processing
#'
#' Creates a data that can be used for image heat map plotting. It uses interactive image import area locator.
#'
#' @param image_df A dataframe of two columns. Column 1 should be named images. Column 2 should be named areas. Both columns should contain characters defining the images that will be imported, and the areas on the images that should be coloured.
#' @param coldata A dataframe containing infromation for the columns of the expression_data
#' @return A list of two elements. First element contains the image_df with a third column defining the coordinates of the areas. The second element contains the image pointer created by image_read
#' @import magick
#' @export

image_data_process <- function(image_df){
  image_df$coords <- NA
  images <- unique(image_df$images)
  image_list <- list()
  count <- 1
  for(j in 1:length(images)){
    message(paste0("Choose the file of the image:\n ", images[j], "\n"))
    Sys.sleep(2)
    path_temp <- file.choose()
    image_temp <- image_read(path_temp)
    image_temp <- image_trim(image_temp)
    areas <- image_df$areas[image_df$images %in% images[j]]
    plot(as.raster(image_temp))
    ## the height coordinates will be handled inversely
    ## so in the function I am going to subtract the retrieved y data from the total height
    imageinfo <- image_info(image_temp)
    tot_height <- imageinfo$height
    for(i in 1:length(areas)){
      message(paste0("Point to the area:\n ", areas[i], "\n"))
      coords <- locator(n = 1)
      image_df[count, "coords"] <- paste0("+", coords$x, "+", tot_height - coords$y)
      count <- count + 1
    }
    image_list[[j]] <- image_temp
  }
  names(image_list) <- images
  image_data <- list(image_df = image_df, image_list = image_list)
  return(image_data)
}
