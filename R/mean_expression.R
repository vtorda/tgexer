#' Average expression
#'
#' Calculate the mean of expession data replicates
#'
#' @param expression_data A matrix containing gene expressin data. Column names should be sample id and row names should be gene ids
#' @param coldata A dataframe containing infromation for the columns of the expression_data
#' @param sample_col Character. The name of the column in the coldata defining samples that has multiple replicates
#' @param sample_col Character. The name of the column in the coldata defining the sample id that can be found in the column nams of the expression_data
#' @return A dataframe containing mean expression data with rownames of gene id and colnams of samples
#' @export


# df_combined <- df_combined[df_combined$groups %in% "TissueSpecFBDevelopment", ]
# expression_data <- exp_data[,colnames(exp_data) %in% df_combined$Library_ID_final]
# coldata <- df_combined
# sample_col <- "sample"
# id_col <- "Library_ID_final"
# i <- 1
mean_expression <- function(expression_data, coldata, sample_col, id_col){
  samples <- unique(coldata[,sample_col])
  mean_expression_df <- as.data.frame(matrix(nrow = nrow(expression_data)))
  for(i in 1:length(samples)){
    idx <- coldata[,sample_col] %in% samples[i]
    ids <- coldata[idx,id_col]
    temp <- expression_data[,colnames(expression_data) %in% ids]
    mean_expression_df <- cbind(mean_expression_df, apply(temp, 1, mean))
  }
  mean_expression_df <- mean_expression_df[,-1]
  colnames(mean_expression_df) <- samples
  return(mean_expression_df)
}
