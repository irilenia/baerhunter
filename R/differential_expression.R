#' Differential gene expression analysis
#'
#' A wrapper function to perform analysis of differential gene expression using DESeq2 method.
#'
#' @param feature_count_file A CSV file containing a table of counts for a particular feature type.
#' @param metadata_file A file outlining the experimental conditions for each sample.
#' @param cutoff_value Retain only the genes that have at least half the samples with the counts above the cut-off value.
#' @param multiple_variables A boolean indicating if the experimental design contained more than one condition.
#' @param main_condition A string indicating which of the conditions is the most important; should coincide with the name in the metadata table.
#' @param output_file_name A string containing the name of the output file.
#'
#' @return DESeq result table and the result summary are written into separate files. Summary file name is derived from the provided output file name with an addition of "summary.txt".
#'
#'
#' @import DESeq2
#' @importFrom stats as.formula
#' @importFrom utils capture.output read.csv read.delim write.table
#'
#' @export
differential_expression <- function(feature_count_file, metadata_file, cutoff_value, multiple_variables, main_condition, output_file_name){

  ##Load in the count and condition tables
  count_data <- as.matrix(read.csv(feature_count_file,sep="\t"))
  condition_metadata <- read.csv(metadata_file, row.names=1)
  coldata <- data.frame()
  ## Check if the samples in the condition table are the same and in the same order as the count tables.
  if (all(rownames(condition_metadata) == colnames(count_data))) {
    coldata <- condition_metadata
  }
  samples <- colnames(count_data)
  coldata <- as.matrix(condition_metadata[samples,])
  rownames(coldata) <- samples
  colnames(coldata) <- colnames(condition_metadata)

  coldata <- as.data.frame(coldata)
  conditions <- colnames(coldata)

  ## Perform the condition model formula.
  DESeq_data <- c()
  ## If there is only one experimental condition.
  if (length(conditions)==1 & multiple_variables==FALSE){
    if (main_condition==colnames(coldata)){
      dds_formula <- as.formula(paste('design=', main_condition, sep = " ~ "))
      DESeq_data <- DESeqDataSetFromMatrix(countData = round(count_data), colData = coldata,design= dds_formula)
    } else {
      return("Argument 'main_condition' does not match to the variable in the metadata table.")
    }
  } else if (length(conditions)>1 & multiple_variables==TRUE){
    ## If there are multiple experimental conditions, incorporate them in the model.
    extra_conditions <- conditions[! conditions %in% main_condition]
    dds_formula <- as.formula(paste('design=', paste(paste(extra_conditions, sep = ' +'), main_condition, sep = " + "), sep = " ~ "))
    DESeq_data <- DESeqDataSetFromMatrix(countData = count_data, colData = coldata,design= dds_formula)
  } else {
    return("Mismatch in variable input.")
  }

  ## Retain only the genes that have at least half the samples with the counts above the cut-off value.
  keep <- rowSums(counts(DESeq_data) > cutoff_value) >=(length(samples)/2)
  DESeq_data <- DESeq_data[keep,]
  ## Perform DESeq2 analysis.
  DESeq_data <- DESeq(DESeq_data)
  DESeq_results <- results(DESeq_data)
  ## Write the analysis results into a file.
  write.table(DESeq_results, output_file_name, sep = "\t")
  ## Write the analysis summary statistics into another file.
  summary_file <- paste(c(sub("^([^.]*).*", "\\1", output_file_name), "_summary.txt"), collapse = "")
  capture.output(summary(DESeq_results), file = summary_file)
  return(DESeq_results)
}
