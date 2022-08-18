#' Peak union calculation
#' 
#' 
#' The function goes over each BAM file in the directory and finds the expression peaks that satisfy the coverage boundary and length criteria in each file. Then it unifies the peak information to obtain a single set of peak genomic coordinates.
#' 
#' @param bam_location The directory containing BAM files.
#' @param bam_txt_list Optional newline separated text file of filenames of bam files. File must be located in bam_location directory. 
#' @param target_strand A character string indicating the strand. Supports two values; '+' and '-'.
#' @param low_coverage_cutoff An integer indicating the low coverage threshold value.
#' @param high_coverage_cutoff An integer indicating the high coverage threshold value.
#' @param peak_width An integer indicating the minimum peak width.
#' @param paired_end_data A boolean indicating if the reads are paired-end.
#' @param strandedness A string outlining the type of the sequencing library: stranded, or reversely stranded.
#' 
#' @return An object of IRanges class, containing the genomic coordinates of selected and unified peaks.
#' 
#' @import IRanges
#' @import GenomicAlignments
#' @import Rsamtools
#' @importFrom utils capture.output read.csv read.delim write.table
#'
#' @export
peak_union_calc <- function(bam_location = ".", bam_txt_list = "", target_strand, low_coverage_cutoff, high_coverage_cutoff,  peak_width, paired_end_data = FALSE, strandedness  = "unstranded") {
  ## Find all BAM files in the directory.
  if (bam_txt_list != ""){
    bam_files <- readLines(bam_txt_list)
    bam_files <- lapply(bam_files, function(x) paste(bam_location, x, sep="/"))
  }else{
    bam_files <- list.files(path = bam_location, pattern = ".BAM$", full.names = TRUE, ignore.case = TRUE)
  }
  peak_union <- IRanges()
  ## Go over each BAM file to extract coverage peaks for a target strand and gradually build a union of all peak sets.
  for (f in bam_files) {
    ## Read a BAM file in accordance with its type and select only the reads aligning to a target strand.
    strand_alignment <- c()
    if (paired_end_data == FALSE & strandedness  == "stranded") {
      file_alignment <- readGAlignments(f)
      strand_alignment <- file_alignment[strand(file_alignment)==target_strand,]
    } else if (paired_end_data == TRUE & strandedness  == "stranded") {
      file_alignment <- readGAlignmentPairs(f, strandMode = 1)
      strand_alignment_unmerged <- file_alignment[strand(file_alignment)==target_strand,]
      #coerce to get single ranges
      strand_alignment <- granges(strand_alignment_unmerged)
    } else if (paired_end_data == FALSE & strandedness  == "reversely_stranded") {
      file_alignment <- readGAlignments(f)
      relevant_strand <- c()
      if (target_strand=="+") {
        relevant_strand <- "-"
      } else {
        relevant_strand <- "+"
      }
      strand_alignment <- file_alignment[strand(file_alignment)==relevant_strand,]
    } else if (paired_end_data == TRUE & strandedness  == "reversely_stranded") {
      file_alignment <- readGAlignmentPairs(f, strandMode = 2)
      strand_alignment_unmerged <- file_alignment[strand(file_alignment)==target_strand,]
      #coerce to get single ranges
      strand_alignment <- granges(strand_alignment_unmerged)
    }
    ## Create a strand coverage vector and extract it
    strand_cvg <- coverage(strand_alignment)
    list_components <- names(strand_cvg)
    target <- c()
    if (length(list_components)==1) {
      target <- list_components
    } else {
      return(paste("Invalid BAM file:",f, sep = " "))
    }
    ## Cut the coverage vector to obtain the expression peaks with the coverage above the low cut-off values.
    peaks <- IRanges::slice(strand_cvg[[target]], lower = low_coverage_cutoff, includeLower=TRUE)
    ## Examine the peaks for the stretches of coverage above the high cut-off. The stretches have to be a defined width.
    test <- viewApply(peaks, function(x) peak_analysis(x,high_coverage_cutoff,peak_width))
    ## Select only the peaks that satisfy the high cut-off condition.
    selected_peaks <- peaks[lapply(test, function(x) !is.null(x))==TRUE]
    ## Convert peak coordinates into IRanges.
    peaks_IRange <- IRanges(start = start(selected_peaks), end = end(selected_peaks))
    ## Calculate the peak union in with the previous peak sets.
    peak_union <- union(peak_union,peaks_IRange)
  }
  return(peak_union)
}

#' Peak checking for the second coverage threshold and width.
#' 
#' This is a helper function that is used to examine if the peak had a continuous stretch of a given width that has coverage above the high cut-off value.
#' 
#' @param View_line A line from a RleViews object.
#' @param high_cutoff An integer indicating the high coverage threshold value.
#' @param min_sRNA_length An integer indicating the minimum sRNA length (peak width).
#' 
#' @return Returns a RleViews line if it satisfies conditions.
#' 
#' @export
#' 
peak_analysis <- function(View_line, high_cutoff, min_sRNA_length) {
  ## This is a helper function that is used to examine if the peak had a continuous stretch of a given width that has coverage above the high cut-off value.
  cvg_string <- as.vector(View_line)
  target_peak <- which(cvg_string>high_cutoff)
  if (length(target_peak)>min_sRNA_length) {
    return(View_line)
  }
}


#' Extract major features from the annotation file
#' 
#' The function extracts parent features only; it also excludes all non-coding RNAs that are already annotated in the file.
#' 
#' @param annotation_file  GFF3 genome annotation file.
#' @param annot_file_directory The directory path for the annotation file (default is '.')
#' @param target_strand A character string indicating the strand. Supports two valies; '+' and '-'.
#' @param original_sRNA_annotation A string indicating how the biotype of pre-annotated ncRNA, which can be found in the attribute column.In case if the user does not know how the sRNA is annotated, it can be set as "unknown". In this case, all RNAs apart from tRNAs and rRNAs will be removed from the selection.
#' 
#' @return A dataframe with the major features from a set strand.
#' 
#' @importFrom utils read.delim
#' @export
major_features <- function(annotation_file, annot_file_directory = ".", target_strand, original_sRNA_annotation) {
  annot_file_loc <- c()
  if (annot_file_directory==".") {
    annot_file_loc <- annotation_file
  } else {
    annot_file_loc <- paste(annot_file_directory, annotation_file, sep = "/")
  }
  
  gff <- read.delim(annot_file_loc, header = FALSE, comment.char = "#")
  ## Pre-annotated sRNAs always have a defined biotype, which can be found in the attribute column.
  ## The following code creates a regex that will recognise pre-annotated sRNAs
  ori_sRNA_biotype <- c()
  ## In case if the user does not know how the sRNA is annotated, it can be set as "unknown". In this case, all RNAs apart from tRNAs and rRNAs will be removed from the selection.
  if (original_sRNA_annotation=="unknown") {
    ori_sRNA_biotype <- "biotype=.*?[^tr]RNA;"
  } else {
    ori_sRNA_biotype <- paste("biotype=", original_sRNA_annotation, sep = "")
  }
  
  ## Select only the major genomic features: remove all child features (like CDS, mRNA etc.), previously annotated sRNAs and extra features
  major_f <- gff[grepl("Parent", gff[,9], ignore.case = TRUE)==FALSE & gff[,3]!='chromosome' & gff[,3]!='biological_region' & grepl(ori_sRNA_biotype, gff[,9], ignore.case = TRUE)==FALSE & gff[,3]!='region' & gff[,3]!='sequence_feature',]
  ## Select only major features for the target strand.
  m_strand_features <- data.frame()
  if (target_strand=="+") {
    m_strand_features <- major_f[major_f[,7]=="+",]
  } else if (target_strand=="-") {
    m_strand_features <- major_f[major_f[,7]=="-",]
  } else {
    return(major_f)
  }
  
  return(m_strand_features)
  
}



#' sRNA prediction
#' 
#' The function for prediction and annotation of sRNA.
#' 
#' @param major_strand_features A dataframe containing the major features for a particular strand.
#' @param target_strand A character string indicating the strand. Supports two values; '+' and '-'.
#' @param union_peak_ranges An IRanges object containing genomic coordinated for all peaks detected on the target strand.
#' 
#' @return An IRanges object containing coordinates and names of the predicted sRNA.
#' 
#' @export
sRNA_calc <- function(major_strand_features, target_strand, union_peak_ranges) {
  ## This function predicts sRNAs.
  ## define function to make sure match is IRanges::match and not base
  match <- function(x, table) IRanges::match(x, table, nomatch = 0)
  ## Convert strand feature coordinates into IRanges.
  strand_IRange <- IRanges(start = major_strand_features[,4], end = major_strand_features[,5])
  ## Select only the ranges that do not overlap the annotated features
  ## Also, disregard the ranges that finish/start 1 position before the genomic feature, because they should be considered as UTRs.
  IGR_sRNAs <- union_peak_ranges[match(union_peak_ranges, subsetByOverlaps(union_peak_ranges, strand_IRange, maxgap = 1L)) == 0,]
  ## Construct the IDs for the new sRNAs to be added into the attribute column of the annotation.
  if (target_strand=="+") {
    names(IGR_sRNAs) <- apply(as.data.frame(IGR_sRNAs),1, function(x) paste("ID=putative_sRNA:p", x[1], "_", x[2], ";", sep = ''))
  } else if (target_strand== "-") {
    names(IGR_sRNAs) <- apply(as.data.frame(IGR_sRNAs),1, function(x) paste("ID=putative_sRNA:m", x[1], "_", x[2], ";", sep = ''))
  } else {
    return("Select strand")
  }
  return(IGR_sRNAs)
}



#' UTR prediction
#' 
#' Function for prediction and annotation of UTRs.
#' 
#' @param major_strand_features A dataframe containing the major features for a particular strand.
#' @param target_strand A character string indicating the strand. Supports two valies; '+' and '-'.
#' @param union_peak_ranges An IRanges object containing genomic coordinated for all peaks detected on the target strand.
#' @param min_UTR_length An integer indicating the minimum UTR length.
#' 
#' @return An IRanges object containing coordinates and names of the predicted UTRs.
#' 
#' @export
UTR_calc <- function(major_strand_features, target_strand, union_peak_ranges, min_UTR_length) {
  ## This function predicts UTRs.
  ## Define function to make sure match is IRanges::match and not base
  match <- function(x, table) IRanges::match(x, table, nomatch = 0)
  ## Convert strand feature coordinates into IRanges.
  strand_IRange <- IRanges(start = major_strand_features[,4], end = major_strand_features[,5])
  ## Find the peak union ranges that overlap with genomic features. Also, include the ranges that do not overlap the features but start/finish 1 position away from it.
  overapping_features <- subsetByOverlaps(union_peak_ranges, strand_IRange, maxgap = 1L)
  ## Join the overlapping features with the genomic ones for further cutting.
  overapping_features <- c(overapping_features, strand_IRange)
  ## Cut the overlapping features on teh border where they overlap with the genomic features.
  split_features <- disjoin(overapping_features)
  ## Now select only the UTR "overhangs" that are created by cutting overlapping features on the border.
  UTRs <- split_features[match(split_features, subsetByOverlaps(split_features, strand_IRange)) == 0]
  ## Select only UTRs that satisfy the minimum length condition.
  UTRs <- UTRs[width(UTRs)>=min_UTR_length,]
  ## Construct the IDs for the new UTRs to be added into the attribute column of the annotation.
  if (target_strand=="+") {
    names(UTRs) <- apply(as.data.frame(UTRs),1, function(x) paste("ID=putative_UTR:p", x[1], "_", x[2],";", sep = ''))
  } else if (target_strand== "-") {
    names(UTRs) <- apply(as.data.frame(UTRs),1, function(x) paste("ID=putative_UTR:m", x[1], "_", x[2],";", sep = ''))
  } else {
    return("Select strand")
  }
  return(UTRs)
  
}



#' Strand annotation
#' 
#' This function constructs the full annotation for the strand.
#' 
#' @param target_strand A character string indicating the strand. Supports two values; '+' and '-'.
#' @param sRNA_IRanges An IRanges object containing coordinates and names of the predicted sRNAs.
#' @param UTR_IRanges An IRanges object containing coordinates and names of the predicted UTRs.
#' @param major_strand_features A dataframe containing the major features for a particular strand.
#' 
#' @return A dataframe containing strand annotation populated with the prediction features, build in accordance with GFF3 file format.
#' 
#' @export
strand_feature_editor <- function(target_strand, sRNA_IRanges, UTR_IRanges, major_strand_features) {
  
  ## Join the sRNA and UTR ranges together.
  sRNA_UTR <- c(sRNA_IRanges, UTR_IRanges)
  
  ## Create information to go into corresponding columns.
  seqid <- rep(major_strand_features[1,1],nrow(as.data.frame(sRNA_UTR)))
  empty_col <- rep(".",nrow(as.data.frame(sRNA_UTR)))
  ## Create a dataframe for sRNAs and UTRs with all GFF3 file columns.
  cmp_strand <- data.frame()
  if (target_strand=="+") {
    cmp_strand <- data.frame(seqid, empty_col, empty_col, as.integer(start(sRNA_UTR)), as.integer(end(sRNA_UTR)), empty_col, rep("+", nrow(as.data.frame(sRNA_UTR))), empty_col, names(sRNA_UTR), stringsAsFactors = FALSE)
  } else if (target_strand== "-") {
    cmp_strand <- data.frame(seqid, empty_col, empty_col, as.integer(start(sRNA_UTR)), as.integer(end(sRNA_UTR)), empty_col, rep("-", nrow(as.data.frame(sRNA_UTR))), empty_col, names(sRNA_UTR), stringsAsFactors = FALSE)
  } else {
    return("Select strand")
  }
  
  names(cmp_strand) <-names(major_strand_features)
  ## Join the sRNA/UTR dataframe with the dataframe for the major features.
  cmp_strand <- rbind(cmp_strand, major_strand_features)
  ## Order the dataframe by the feature start position.
  cmp_strand <- cmp_strand[order(cmp_strand[,4]),]
  
  ## Set the previous feature name to be the ID of the last feature in the chromosome, accounting for the fact that bacterial genomes are circular.
  previous_feature_name <- sub("ID=.*?:(.*?);.*", "\\1", cmp_strand[nrow(cmp_strand),9])
  
  for (i in 1:nrow(cmp_strand)) {
    ## Determine feature type from the attribute column if the third column is empty.
    feature_name <- sub("ID=.*?:(.*?);.*", "\\1", cmp_strand[i,9])
    if (cmp_strand[i,3]==".") {
      feature_type <- sub("ID=(.*?):.*?;.*", "\\1", cmp_strand[i,9])
      cmp_strand[i,3] <- feature_type
      ## Find the name of the next feature in the annotation.
      next_feature_name <- c()
      if (i+1 <= nrow(cmp_strand)) {
        next_feature_name <- sub("ID=.*?:(.*?);.*", "\\1", cmp_strand[i+1,9])
      } else {
        next_feature_name <- sub("ID=.*?:(.*?);.*", "\\1", cmp_strand[1,9])
      }
      
      ## Build feature attribute column information for sRNAs and UTRs, including upstream and downstream features (with reagards to the strand).
      feature_attribute <- c()
      
      if (target_strand=="+") {
        feature_attribute <- paste(cmp_strand[i,9],"upstream_feature=", previous_feature_name, ";downstream_feature=", next_feature_name, sep = "")
      } else if (target_strand== "-") {
        feature_attribute <- paste(cmp_strand[i,9],"upstream_feature=", next_feature_name, ";downstream_feature=", previous_feature_name, sep = "")
      } else {
        return("Select strand")
      }
      cmp_strand[i,9] <- feature_attribute
      
    }
    previous_feature_name <- feature_name
  }
  
  return(cmp_strand)
  
}





#' Prediction and annotation of sRNAs and UTRs from RNA-seq data
#' 
#' A wrapper function that executes all prediction steps for each strand and builds the final GFF3 annotation.
#' 
#' @param bam_directory The directory containing BAM files.
#' @param bam_list Optional newline separated text file of filenames of bam files. File must be located in bam_location directory. 
#' @param original_annotation_file GFF3 genome annotation file.
#' @param annot_file_dir The directory containing the GFF3 annotation file.
#' @param output_file A string containing the name of an output file.
#' @param original_sRNA_annotation A string indicating how the biotype of pre-annotated ncRNA, which can be found in the attribute column.In case if the user does not know how the sRNA is annotated, it can be set as "unknown". In this case, all RNAs apart from tRNAs and rRNAs will be removed from the selection.
#' @param low_coverage_cutoff An integer indicating the low coverage threshold value.
#' @param high_coverage_cutoff An integer indicating the high coverage threshold value.
#' @param min_sRNA_length An integer indicating the minimum peak width/sRNA length.
#' @param min_UTR_length An integer indicating the minimum UTR length.
#' @param paired_end_data A boolean indicating if the reads are paired-end.
#' @param strandedness A string outlining the type of the sequencing library: stranded, or reversely stranded.
#' 
#' @return Outputs a new GFF3 file populated with predicted sRNAs and UTRs.
#'
#' 
#' @export
feature_file_editor <- function(bam_directory = ".", bam_list = "", original_annotation_file, annot_file_dir = ".", output_file, original_sRNA_annotation, low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, min_UTR_length, paired_end_data = FALSE, strandedness  = "stranded") {
  test <- list.files(path = bam_directory, pattern = ".BAM$", full.names = TRUE, ignore.case = TRUE)
  if (length(test) > 0){
    ## Plus strand
    plus_strand_peaks <- peak_union_calc(bam_location = bam_directory, bam_txt_list = bam_list, "+", low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, paired_end_data, strandedness)
    print("Extracted plus strand data from BAM files")
    maj_plus_features <- major_features(original_annotation_file, annot_file_directory = annot_file_dir, "+", original_sRNA_annotation)
    plus_sRNA <- sRNA_calc(maj_plus_features, "+", plus_strand_peaks)
    plus_UTR <- UTR_calc(maj_plus_features, "+", plus_strand_peaks, min_sRNA_length)
    plus_annot_dataframe <- strand_feature_editor("+", plus_sRNA, plus_UTR, maj_plus_features)
    print("Built plus strand annotation dataframe")
    ## Minus strand
    minus_strand_peaks <- peak_union_calc(bam_location = bam_directory, bam_txt_list = bam_list, "-", low_coverage_cutoff, high_coverage_cutoff, min_sRNA_length, paired_end_data, strandedness)
    print("Extracted minus strand data from BAM files")
    maj_minus_features <- major_features(original_annotation_file, annot_file_directory = annot_file_dir, "-", original_sRNA_annotation)
    minus_sRNA <- sRNA_calc(maj_minus_features, "-", minus_strand_peaks)
    minus_UTR <- UTR_calc(maj_minus_features, "-", minus_strand_peaks, min_UTR_length)
    minus_annot_dataframe <- strand_feature_editor("-", minus_sRNA, minus_UTR, maj_minus_features)
    print("Built minus strand annotation dataframe")
  
    ## Creating the final annotation dataframe by combining both strand dataframe and adding missing information like child features from the original GFF3 file.
    annot_file_loc <- c()
    if (annot_file_dir==".") {
      annot_file_loc <- original_annotation_file
    } else {
      annot_file_loc <- paste(annot_file_dir, original_annotation_file, sep = "/")
    }
  
    gff <- read.delim(annot_file_loc, header = FALSE, comment.char = "#")
    annotation_dataframe <- rbind(gff, plus_annot_dataframe, minus_annot_dataframe)
    ## Remove all teh repeating information.
    annotation_dataframe <- unique(annotation_dataframe)
    ## Order the dataframe by feature start coordinates.
    annotation_dataframe <- annotation_dataframe[order(annotation_dataframe[,4]),]
    print("Prepared complete annotation dataframe")
  
    ## Restore the original header.
    f <- readLines(annot_file_loc)
    header <- c()
    i <- 1
    while (grepl("#",f[i])==TRUE) {
      f_line <- f[i]
      header <- c(header,f_line)
      i <- i+1
    }
    # add a line to indicate the origin of the file (single # commas should be ignored by programs)
    header <- c(header, "# produced by baerhunter")
  
    print("Building output file now")
  
    ## Create the final GFF3 file.
    write.table(header, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    write.table(annotation_dataframe, output_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
  
    return("Done!")
  }else{
    return("No BAMs in bam directory!")
  }
}