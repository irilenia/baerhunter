#' read_annotation_file function
#'
#' This function pastes together the path and filename of annotation file (if not in current directory) and tests file for existence.
#'
#' @param annot_dir Directory where genome annotation file is located
#' @param annot_file  GFF3 or GTF genome annotation file
#'
#' @return annot_file_loc Complete path of existing annotation file
#'
#' @importFrom assertthat assert_that
#' @export
read_annotation_file <- function(annot_dir, annot_file){
  #read in annotation file
  annot_file_loc <- c()
  if (annot_dir==".") {
    annot_file_loc <- annot_file
    assert_that(file.access(annot_file_loc, mode=0) == 0, msg="Annotation file not found")
  } else {
    annot_file_loc <- paste(annot_dir, annot_file, sep = "/")
    assert_that(file.access(annot_file_loc, mode=0) == 0, msg="Annotation file not found")
  }
  return(annot_file_loc)
}


#' find_strandedness function
#'
#' This function translates the user-inputted strandedness parameter to the required integer input for strandSpecific arg of featureCounts.
#'
#' @param strand_param user input 'stranded' or 'reversely-stranded'
#'
#' @return strand_sp integer value
#'
#' @importFrom assertthat assert_that
#' @export
find_strandedness <- function(strand_param){
  strand_sp <- integer()
  strand_strings <- c("unstranded", "reversely_stranded", "stranded")
  assert_that(strand_param %in% strand_strings,
              msg="Invalid strandedness parameter: must be either 'unstranded', 'stranded' or 'reversely_stranded'")
  if(strand_param=="unstranded"){
    strand_sp <- 0
  } else if(strand_param=="stranded"){
    strand_sp <- 1
  } else if(strand_param=="reversely_stranded"){
    strand_sp <- 2
  } #else {
  #  strand_sp <- 0
  #}
  return(strand_sp)
}


#' SAF converter
#'
#' This function converts gff file into simplified annotation format with 5 columns, appropriate for use in featureCounts/tpm_norm_flagging.
#'
#' @param ann_file A GFF3 annotation file.
#' @param exclude A boolean to indicate whether or not to include rRNA/tRNA features
#'
#' @return A dataframe in Simplified Annotation Format
#'
#' @importFrom stringr str_match
#' @importFrom assertthat assert_that
#' @export
make_saf <- function(ann_file, exclude=F){
  #function to create SAF file from gff from feature file editor
  gff <- read.delim(ann_file, header = FALSE, comment.char = "#")
  ## check that file is correct format (9 cols)
  assert_that(ncol(gff)==9, msg="annotation file format is invalid")
  ## Select only the major genomic features: remove all child features (like CDS, mRNA etc.) and extra features
  if (exclude==F){
    major_f <- gff[grepl("Parent", gff[,9], ignore.case = TRUE)==FALSE &
                     gff[,3]!='chromosome' & gff[,3]!='biological_region' &
                     gff[,3]!='region' & gff[,3]!='sequence_feature',]
  }else{ # if you want to exclude rRNA and tRNA features)
    major_f <- gff[grepl("Parent", gff[,9], ignore.case = TRUE)==FALSE &
                     gff[,3]!='chromosome' & gff[,3]!='biological_region' &
                     gff[,3]!='region' &
                     gff[,3]!='sequence_feature' &
                     gff[,3]!='tRNA_gene' &
                     gff[,3]!='rRNA_gene',]
  }
  saf_df <- as.data.frame(matrix(0, ncol=5, nrow=nrow(major_f)))
  colnames(saf_df) <- c("GeneID", "Chr", "Start", "End", "Strand")
  saf_df$GeneID   <- str_match(major_f$V9, "ID=(.*?);")[,2]
  saf_df$Chr      <- major_f$V1
  saf_df$Start    <- major_f$V4
  saf_df$End      <- major_f$V5
  saf_df$Strand   <- major_f$V7
  return(saf_df)
}

#' count_features function
#'
#' This is a function to employ Rsubread featureCounts to quantify expression of annotated and predicted elements
#'
#' @param bam_dir The directory where bam files located
#' @param annotation_dir The directory where annotation file is located
#' @param annotation_file The complete annotation file: GFF3 or GTF genome annotation file
#' @param output_dir The full directory path for CSV output files to be written
#' @param output_filename The name for the output files--for example dataset name
#' @param chromosome_alias_file A comma-delimited TXT file containing a character string with the chromosome names. This file has to have two columns: first with the chromosome name in the annotation file, second with the chromosome name in the BAM file.
#' @param strandedness A string outlining the type of the sequencing library: unstranded, stranded, or reversely stranded.
#' @param is_paired_end A boolean indicating if the reads are paired-end.
#' @param excl_rna A boolean indicating if misc RNA features (rRNA, tRNA) are excluded from quantification. (Defaults=T)
#' @param ... Optional parameters passed on to featureCounts()  Default: allowMultiOverlap = T, fraction = T
#'
#' @return Count tables for each feature are written into separate files, as well as the result summary.
#'
#' @import Rsubread
#' @importFrom tools file_path_sans_ext
#' @importFrom utils write.table
#' @importFrom assertthat assert_that
#' @export
count_features <- function(bam_dir=".",
                           annotation_dir=".",
                           annotation_file,
                           output_dir=".",
                           output_filename="dataset",
                           chromosome_alias_file,
                           strandedness,
                           is_paired_end,
                           excl_rna = T,
                           ...){
  ## function to call rsubread featureCounts

  ##read in annotation file (in gff3 form)
  annot_file_loc <- read_annotation_file(annotation_dir, annotation_file)

  ## Compile a list of BAM files present in the bam directory
  bam_files <- list.files(path = bam_dir, pattern = ".BAM$", full.names = TRUE, ignore.case = TRUE)
  assert_that(length(bam_files) > 0 , msg="Empty bam directory")

  ## Convert gff to SAF dataframe
  nsaf_df <- make_saf(ann_file=annot_file_loc, exclude=excl_rna)

  ## if output directory exists, create filenames and path to output file
  assert_that(dir.exists(output_dir), msg="Output directory doesn't exist")
  output_file  <- paste(output_dir, output_filename, sep = "/")
  count_file_name <- paste(output_file, "_Counts.csv", sep = "")
  summary_file_name <- paste(output_file, "_Count_summary.csv", sep="")

  ## The strandedness set by the user is translated into the integer for strandSpecific argument of the featureCounts function.
  strand_specific <- find_strandedness(strandedness)

  ## Paired-end is FALSE by default unless specified by the user otherwise.
  paired_end <- FALSE
  if(is_paired_end==TRUE){
    paired_end <- TRUE
  }

  ## Extract BAM file names without extension.
  sample_names <- c(file_path_sans_ext(basename(bam_files)))

  ## Execute featureCounts function. Feature counts and summary stats are written into separate TAB delimited files.
  fc <- featureCounts(bam_files,
                      annot.ext = nsaf_df,
                      chrAliases = chromosome_alias_file,
                      strandSpecific = strand_specific,
                      isPairedEnd = paired_end,
                      allowMultiOverlap = T,
                      fraction = T,
                      ...)
  colnames(fc$counts) <- sample_names
  write.table(fc$counts, count_file_name, sep = "\t")
  colnames(fc$stat) <- c('Status',sample_names)
  write.table(fc$stat, summary_file_name, sep = "\t")

  return("Done!")
}
