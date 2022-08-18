#' TPM calculator
#'
#' This function uses feature count tables to calculate TPM values for each gene and sample.
#'
#' @param count_table A CSV file containing feature counts for each sample.
#' @param complete_ann A GFF3 annotation file or SAF dataframe
#' @param feature_type A string indicating desired feature type(s) from annotation.
#' @param is_gff A boolean indicating whether annotation is gff file, default=T
#' @param excl_rna A boolean indicating if misc RNA features (rRNA, tRNA) are counted in addition to CDS and predicted RNA elements. (Defaults=T)
#' @param output_file A string indicating the name of the output file.
#'
#' @return A dataframe with TPM values for each gene and sample; the same is written into the output file.
#'
#'
#' @importFrom utils read.delim write.table
#' @importFrom assertthat assert_that
#'
#' @export
tpm_normalisation <- function(count_table, complete_ann, feature_type = c("putative_sRNA", "putative_UTR"), is_gff = T, output_file = NA, excl_rna = T) {
  ## check output directory exists
  if (is.na(output_file)==FALSE) {
    out_dir <- dirname(output_file)
    assert_that(dir.exists(out_dir), msg="Output directory doesn't exist.")
  }
  ## make saf from gff (uses make_saf function)
  if (is_gff == T){
    nsaf_df <- make_saf(ann_file=complete_ann, exclude = excl_rna)
  }else{
    nsaf_df <- complete_ann
  }
  feature_names <- nsaf_df$GeneID

  ## Load in the count table and filter for feature types
  count_df <- read.delim(count_table)

  ## Calculate the length of the features, if they are in the right order.
  assert_that(all(rownames(count_df) == feature_names), msg="Wrong feature order")
  feature_lengths <- c()
  feature_lengths <- (nsaf_df$End-nsaf_df$Start+1)/1000

  ## Calculate RPK by dividing the feature count of each gene (per feature) by its length in kilobases.
  rpk_df <- data.frame( count_df[,1] / feature_lengths )
  if (ncol(count_df) > 1){
    for (i in 2:ncol(count_df)) {
      sample_rpk <- count_df[ , i] / feature_lengths
      rpk_df     <- data.frame(rpk_df, sample_rpk)
    }
  }
  ## Calculate the scaling factor for each sample by summing up all RPKs per sample and dividing by a million.
  sample_rpk_sum <- apply(rpk_df, 2, sum)
  scaling_fact <- sample_rpk_sum/1000000

  ## Divide all RPK by the corresponding sampling factor.
  tpm_df <- data.frame(rpk_df[,1]/scaling_fact[1])
  if (ncol(rpk_df) > 1){
    for (n in 2:ncol(rpk_df)) {
      sample_tpm <- rpk_df[,n]/scaling_fact[n]
      tpm_df <- data.frame(tpm_df, sample_tpm)
    }
  }

  colnames(tpm_df) <- colnames(count_df)
  rownames(tpm_df) <- feature_names

  ## If the output file is set by the user, write the TPM dataframe into it.
  if (is.na(output_file)==FALSE) {
    write.table(tpm_df, output_file, sep = "\t")
  }
  return(tpm_df)
}


#' Flagging features depending on TPM value profile
#'
#' A helper function to analyse each row of the TPM table. Each feature gets allocated a flag depending on the expression profile.
#'
#' @param tpm_data A CSV file containing TPM values for each feature in each sample.
#' @param complete_annotation A GFF3 annotation file.
#' @param output_file A string indicating the name of the output file.
#'
#' @return A Gff3 file, where a target feature type has an expression flag added to its attribute column.
#'
#'
#' @importFrom utils read.delim write.table
#' @export
tpm_flagging <- function(tpm_data, complete_annotation, output_file) {


  tpm_analyser <- function(num_vec) {
    if (all(num_vec<0.5)) {
      return("expression_below_cutoff")
    } else  if (any(num_vec > 1000)) {
      return("high_expression_hit")
    } else if (any((num_vec > 10) & (num_vec <= 1000))) {
      return("medium_expression_hit")
    } else if (any((num_vec >= 0.5) & (num_vec <= 10))) {
      return("low_expression_hit")
    }
  }
  ## Read the TPM table and create a flag vector in the corresponding gene order.
  norm_data <- read.delim(tpm_data, header = TRUE)
  # return if the data frame contains non-numeric data (comparisons will fail otherwise)
  if (!is.numeric(as.matrix(norm_data))) {
    stop("Input file contains unexpected non-numeric data\n")
  }
  flags <- apply(norm_data, 1, function(x) tpm_analyser(as.vector(x)))
  flag_names <- names(flags)
  ann_file <- readLines(complete_annotation)
  ## Add the flag to the corresponding feature's attribute column.
  new_annot <- c()
  for (i in ann_file) {
    feature_name <- sub(".*?ID=(.*?:.*?);.*", "\\1", i)

    if (feature_name %in% flag_names) {
      new_line <- paste(i, ";expression_flag=", flags[feature_name], sep = "")
      new_annot <- c(new_annot, new_line)
    } else {
      new_annot <- c(new_annot, i)
    }
  }

  write.table(new_annot, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

}


#' Filtering features by a target flag.
#'
#' A function to filter the marked features selected by the user by the flag of choice.
#'
#' @param complete_annotation_file A flagged GFF3 annotation file.
#' @param target_features A string indicating feature type to filter (optional).
#' @param target_flag A string indicating a flag for filtering.
#' @param output_file A string indicating the name of the output file.
#'
#' @return A GFF3 file where a target feature type is filtered by the expression flag of interest.
#'
#'
#' @importFrom utils read.delim write.table
#' @export
tpm_flag_filtering <- function(complete_annotation_file, target_features = c("putative_sRNA", "putative_UTR"), target_flag, output_file) {

  ##load in annotation data.
  annot_data <- read.delim(complete_annotation_file, header = FALSE, comment.char = "#")

  ## An internal function to go examine a table row: all target features are checked and filtered by the desired flag; all the other features are kept.
  selection <- function(table_row, target_features, target_flag) {
    if ((as.character(table_row[[3]]) %in% target_features)==TRUE) {
      if (grepl(target_flag, as.character(table_row[[9]]), ignore.case = TRUE)==TRUE){
        return(table_row)
      }
    } else if ((as.character(table_row[[3]]) %in% target_features)==FALSE) {
      return(table_row)
    }
  }
  ## Select the features by the flag and create a new dataframe.
  filtered_selection <- apply(annot_data, 1, function(x) selection(x, target_features, target_flag))
  selection_not_null <- filtered_selection[-which(sapply(filtered_selection, is.null))]
  #IN df <- data.frame(matrix(unlist(selection_not_null), nrow=length(d), byrow=T))
  df <- data.frame(matrix(unlist(selection_not_null), nrow=length(selection_not_null), byrow=T))

  ## Restore the original header.
  f <- readLines(complete_annotation_file)
  header <- c()
  i <- 1
  while (grepl("#",f[i])==TRUE) {
    header <- c(header,f[i])
    i <- i+1
  }
  ## Write a new GFF3 file.
  write.table(header, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(df, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)

}
