#' Assigning reads to genomic features
#' 
#' A wrapper function allowing to count reads to different typed of genomic features in one go. Designed to incorporate newly-discovered sRNAs and UTRs into the analysis.
#' 
#' @param bam_dir The directory where bam files can be found.
#' @param annotation_dir The directory where the annotation file can be found.
#' @param output_dir The directory where the resulting CSV files will be deposited.
#' @param annotation_file A GFF3 or GTF genome annotation file.
#' @param chromosome_alias_file A comma-delimited TXT file containing a character string with the chromosome names. This file has to have two columns: first with the chromosome name in the annotation file, second with the chromosome name in the BAM file.
#' @param target_features A vector of strings outlining all feature types to be examined.
#' @param strandedness A string outlining the type of the sequencing library: unstranded, stranded, or reversely stranded.
#' @param is_paired_end A boolean indicating if the reads are paired-end.
#' @param ... Optional parameters passed on to featureCounts()  
#' @return Count tables for each feature are written into separate files, as well as the result summary.
#' 
#'
#' @import Rsubread
#' @importFrom tools file_path_sans_ext
#' @importFrom utils write.table
#' 
#' @export
count_features <- function(bam_dir=".", annotation_dir=".", output_dir=".", annotation_file, chromosome_alias_file, target_features, strandedness, is_paired_end, ...){
  ## Compile a list of BAM files present in the Data folder.
  bam_files <- list.files(path = bam_dir, pattern = ".BAM$", full.names = TRUE, ignore.case = TRUE)
  #IN bam_files <- list.files(pattern = ".BAM$", full.names = TRUE, ignore.case = TRUE)
  ## Feature types that will be used for counting.
  feature_types <- target_features
    
  ## Attribute type is set to a default of GTF file format, but if the annotation provided is in GFF3 format, the type is adjusted accordingly.
  attribute_type <- "gene_id"
  gff3_files <- list.files(path=annotation_dir, pattern = ".gff3")
  #IN gff3_files <- list.files(pattern = ".gff3")
  if(is.element(annotation_file, gff3_files)){
    attribute_type <- "ID"
  }
  
  ## The strandedness set by the user is translated into the integer for strandSpecific argument of the featureCounts function.
  strand_specific <- integer()
  if(strandedness=="unstranded"){
    strand_specific <- 0
  } else if(strandedness=="stranded"){
    strand_specific <- 1
  } else if(strandedness=="reversely_stranded"){
    strand_specific <- 2
  } else {
    strand_specific <- 0
  }
  
  ## Paired-end is FALSE by default unless specified by the user otherwise.
  paired_end <- FALSE
  if(is_paired_end==TRUE){
    paired_end <- TRUE
  }
  
  ## Extract BAM file names without extension.
  sample_names <- c(file_path_sans_ext(basename(bam_files)))
  
  ## Execute featureCounts function for each of the feature types. Feature counts and summary stats are written into separate TAB delimited files.
  for (i in 1:length(feature_types)) {
    fc <- featureCounts(bam_files, annot.ext = paste(annotation_dir, annotation_file, sep = "/"),  isGTFAnnotationFile = TRUE, GTF.featureType = feature_types[i], GTF.attrType = attribute_type,chrAliases = chromosome_alias_file, strandSpecific = strand_specific, isPairedEnd = paired_end, ...)
    colnames(fc$counts) <- sample_names
    count_file_name <- paste(output_dir, feature_types[i], "_Counts.csv", sep = "")
    write.table(fc$counts, count_file_name, sep = "\t")
    colnames(fc$stat) <- c('Status',sample_names)
    #IN added output_dir to name
    summary_file_name <- paste(output_dir, feature_types[i], "Count_summary.csv", sep = "")
    write.table(fc$stat, summary_file_name, sep = "\t")
  }
  return("Done!")
}

