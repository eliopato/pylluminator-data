# load all versions of the genome and save relevant information for pylluminator as csv files, in the _generated_data/genome_info folder

# package installation
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("sesame")

library('sesame')
library('sesameData')

sesameDataCache() # first time only

# if using r studio
current_directory <- dirname(rstudioapi::getSourceEditorContext()$path)
# otherwise, manually specify the script location ..
# current_directory <- "C:/path/to/your/script"

print(paste('source directory is', current_directory))

for(genome_version in c("hg38", "hg19", "mm10", "mm39")) {
  print(paste('getting ', genome_version))
  
  # define and create output directory
  directory <- paste0(current_directory, "/_generated_data/genome_info/", genome_version)
  dir.create(directory, showWarnings = FALSE, recursive=TRUE)

  # get genome data
  genomeInfo <- sesameDataGet(paste0("genomeInfo.", genome_version))

  # reformat and save chromosome regions
  chromosome_regions <- genomeInfo$cytoBand
  names(chromosome_regions) <- c('chromosome', 'start', 'end', 'name', 'giemsa_staining')
  chromosome_regions$chromosome <- sub("^chr", "", chromosome_regions$chromosome) # Remove 'chr' prefix
  write.csv(chromosome_regions, paste0(directory, '/chromosome_regions.csv'), row.names=FALSE, quote=FALSE)

  # reformat and save chromosome sequence length
  seq_length <- as.data.frame(genomeInfo$seqLength)
  seq_length$Chromosome <- row.names(seq_length)
  names(seq_length) <- c('seq_length', 'chromosome')
  seq_length$chromosome <- sub("^chr", "", seq_length$chromosome) # Remove 'chr' prefix
  write.csv(seq_length[c('chromosome', 'seq_length')], paste0(directory, '/seq_length.csv'), row.names=FALSE, quote=FALSE)

  # reformat and save chromosome gap information
  gap_info <- as.data.frame(genomeInfo$gapInfo)
  names(gap_info) <- c('chromosome','start','end','width','strand','type','bridge')
  gap_info$chromosome <- sub("^chr", "", gap_info$chromosome) # Remove 'chr' prefix
  write.csv(gap_info, paste0(directory, '/gap_info.csv'), row.names=FALSE, quote=FALSE)

  # reformat and save transcript information
  transcript_list <- as.data.frame(genomeInfo$txns)[c('group_name', 'start', 'end', 'width', 'exon_number')]
  write.csv(transcript_list, paste0(directory, '/transcripts_list.csv'), row.names=FALSE, quote=FALSE)
  transcript_exons <- as.data.frame(GenomicRanges::mcols(genomeInfo$txns))
  names(transcript_exons) <- c("chromosome", "transcript_start", "transcript_end",  "transcript_strand", "transcript_id", "transcript_type",
                               "transcript_name", "gene_name", "gene_id", "gene_type", "source", "level", "cds_start", "cds_end")
  transcript_exons$chromosome <- sub("^chr", "", transcript_exons$chromosome) # Remove 'chr' prefix
  write.csv(transcript_exons, paste0(directory, '/transcripts_exons.csv'), row.names=FALSE, quote=FALSE)

  # not useful ?
  # write.csv(genomeInfo$cen_info, paste0(directory, '/cen_info.csv'), row.names=FALSE)

}
