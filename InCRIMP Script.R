#=============================================================================
#
#    Code chunk 0: Pre-requisites
#
#=============================================================================

# devtools::install_github("r-dbi/DBI")
# devtools::install_github("r-dbi/RMariaDB")
# devtools::install_github("truecluster/bit64")
# devtools::install_github("tidyverse/blob")
# devtools::install_github("tidyverse/lubridate")
# devtools::install_github("tidyverse/hms")
# devtools::install_github("RcppCore/Rcpp")
library(DBI)
library(RMariaDB)
library(bigrquery)
library(dplyr)
library(tidyverse)
library(influential)
library(httr)
library(jsonlite)
library(VariantAnnotation) ## A very good package for preprocessing and preparing the VCF files for downstream analyses
library(parallel)

#=============================================================================
#
#    Code chunk 1: Prepare the settings fro retrieving data from Database
#
#=============================================================================

# Connecting to the database

## Connect to the database
db <- dbConnect(
  drv = RMariaDB::MariaDB(),
  default.file=".my.cnf",
  groups="azure_lkcgp_curated" # must be the same as the name between [] in .my.cnf file
)

#=============================================================================
#
#    Code chunk 2: Get the sample data from Database
#
#=============================================================================

## Get Patient and Sample data
lkcgp_sample <- tbl(db, "lkcgp_sample") %>% as.data.frame()

### Add diagnosis instead of those cancer_types that are NA
lkcgp_sample$cancer_type[which(is.na(lkcgp_sample$cancer_type))] <- lkcgp_sample$diagnosis[which(is.na(lkcgp_sample$cancer_type))]

## Retain only PRISM study data
lkcgp_sample <- subset(lkcgp_sample, study == "PRISM")

###########################

## Prepare the sample tables for publication
somatic_samples_tbl <- somatic_snv %>% dplyr::select(sample_id, cancer_type) %>% dplyr::distinct(sample_id, .keep_all = T)
germline_samples_tbl <- germline_snv %>% dplyr::select(sample_id, cancer_type) %>% dplyr::distinct(sample_id, .keep_all = T)

somatic_samples_tbl <- table(somatic_samples_tbl$cancer_type) %>% data.frame()
colnames(somatic_samples_tbl) <- c("Cancer_type", "Number_of_samples")
somatic_samples_tbl$Context <- 'Somatic'

germline_samples_tbl <- table(germline_samples_tbl$cancer_type) %>% data.frame()
colnames(germline_samples_tbl) <- c("Cancer_type", "Number_of_samples")
germline_samples_tbl$Context <- 'Germline'

combined_samples_tbl <- rbind(somatic_samples_tbl, germline_samples_tbl)

## Unify names
combined_samples_tbl <- combined_samples_tbl[-which(combined_samples_tbl$Cancer_type == "Glioma Other"),]
combined_samples_tbl[which(combined_samples_tbl$Cancer_type == "Glioma other" & combined_samples_tbl$Context == 'Somatic'),2] <- 15

write_csv(x = combined_samples_tbl, file = "Results/combined_samples_tbl.csv")


#=============================================================================
#
#    Code chunk 3: Get and prepare the germline SNV data
#
#=============================================================================

## Get the SNV data

### Germline variants
germline_snv <- tbl(db, "lkcgp_curated_germline_snv") %>% as.data.frame()
sample_germline_snv <- tbl(db, "lkcgp_curated_sample_germline_snv") %>% as.data.frame()

### Complement the dataset
germline_snv <- germline_snv[c(germline_snv$variant_id %in% sample_germline_snv$variant_id),]
germline_snv$exon <- NULL

sample_germline_snv <- sample_germline_snv[,c(1:12)]
sample_germline_snv <- cbind(germline_snv[match(sample_germline_snv$variant_id, germline_snv$variant_id),],
                             sample_germline_snv[,-3])

germline_snv <- sample_germline_snv
rm(sample_germline_snv)

germline_snv <- germline_snv[c(germline_snv$sample_id %in% lkcgp_sample$matched_normal_id),]
germline_snv$study <- lkcgp_sample$study[match(germline_snv$sample_id, lkcgp_sample$matched_normal_id)]

## Retain only PRISM study data
germline_snv <- subset(germline_snv, study == "PRISM")

## Split the variants that have two alternative alterations (denoted by a comma) into two
germline_snv <- tidyr::separate_rows(germline_snv, alt)

## Calculate nalt for all variants
germline_snv_nalt <- lapply(unique(germline_snv$variant_id), function(i) {
  tmp_data <- subset(germline_snv, variant_id == i)
  ref <- tmp_data$ref[1]
  genotypes <- tmp_data$genotype %>% 
    stringr::str_split(pattern = "/", n = 2) %>% unlist()
  count <- length(genotypes) - grep(pattern = paste0("^", ref, "$"), x = genotypes) %>% length()
  tmp_tbl <- data.frame(variant_id = i, nalt = count)
})

germline_snv_nalt <- do.call(rbind, germline_snv_nalt)

germline_snv$nalt <- germline_snv_nalt$nalt[match(germline_snv$variant_id, germline_snv_nalt$variant_id)]

#### Add patient ID, cancer type, Study, and lastUpdated date
germline_snv$patient_id <- lkcgp_sample$patient_id[match(germline_snv$sample_id, lkcgp_sample$matched_normal_id)]
germline_snv$cancer_type <- lkcgp_sample$cancer_type[match(germline_snv$sample_id, lkcgp_sample$matched_normal_id)]
germline_snv$lastUpdated <- lkcgp_sample$lastUpdated[match(germline_snv$sample_id, lkcgp_sample$matched_normal_id)]

## Remove variants with 0 or NA depth
germline_snv <- germline_snv[-c(which(germline_snv$depth == 0 | is.na(germline_snv$depth))),]

## Remove variants with <= 3 altAD or NA altAD
germline_snv <- germline_snv[-c(which(germline_snv$altad <= 3 | is.na(germline_snv$altad))),]

## Remove the stars (*) from alt column as these denote spanning deletions and the the effects of the their other corresponding variants will be captured by their own data record
germline_snv$alt <- gsub(pattern = ",[*]|[*],", replacement = "", x = germline_snv$alt)

# Process downloaded gnomAD pop freq. data

### First, get min and max genomic positions of all chromosomes and create batches of 100,000 bps (genomic positions) for each chromosome
germline_snv_chr <- 
  unique(germline_snv$chr)

germline_snv_chr_tbl <-
parallel::mclapply(1:length(germline_snv_chr), FUN = function(i) {
  
  ## Calculate the min and max of chromosome i
  germline_snv_chr_min <- min(germline_snv$pos[which(germline_snv$chr %in% germline_snv_chr[i])]) - 10
  germline_snv_chr_max <- max(germline_snv$pos[which(germline_snv$chr %in% germline_snv_chr[i])]) + 10
  
  germline_snv_chr_tbl_tmp <- data.frame(chr = germline_snv_chr[i], 
                                         min = germline_snv_chr_min, 
                                         max = germline_snv_chr_max)
  
}, mc.cores = parallel::detectCores() - 3)

germline_snv_chr_tbl <-  do.call(rbind, germline_snv_chr_tbl)

############################

## Get the header of VCF file
genomes.r2.1.1.exome_calling_header <- VariantAnnotation::scanVcfHeader(file = "Datasets/gnomAD/V2/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz")

## Filter out chromosomes of germline_snv_chr_tbl that are not available in the genomes.r2.1.1.exome_calling_header
genomes.r2.1.1.exome_calling_chrs <- GenomeInfoDb::seqlevels(genomes.r2.1.1.exome_calling_header)
germline_snv_chr_gnomAD_tbl <- germline_snv_chr_tbl[which(germline_snv_chr_tbl$chr %in% genomes.r2.1.1.exome_calling_chrs),]

## Remove Y chromosome as it is not present in the Tabix index file
germline_snv_chr_gnomAD_tbl <- germline_snv_chr_gnomAD_tbl[-which(germline_snv_chr_gnomAD_tbl$chr == "Y"),]

## Get the info column names
infoCols <- rownames(VariantAnnotation::info(genomes.r2.1.1.exome_calling_header))

## Separate the AF column names
AFcols <- infoCols[grep("^AF", infoCols)]

## Filter out the first and third column names (overall AF and AF_raw) and AF_popmax as well as male/female sub-populations as we are not looking for sex-specific variant
AFcols <- AFcols[-c(1, grep("male|raw|popmax", AFcols))]

germline_snv_popFreq_list <- parallel::mclapply(1:nrow(germline_snv_chr_gnomAD_tbl), function(i) {
  
  ## Define a genomic ranges
  chr = germline_snv_chr_gnomAD_tbl[i,1]
  start.loc <- germline_snv_chr_gnomAD_tbl[i,2]
  end.loc <- germline_snv_chr_gnomAD_tbl[i,3]
  gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(start.loc, end.loc))
  
  ## restrict VCF INFO columns to AC and AN values
  vcfPar <- VariantAnnotation::ScanVcfParam(geno = NA,
                                            fixed = c("ALT", "FILTER"),
                                            info = AFcols,
                                            which = gr)
  
  # Load the desired chunk of VCF file
  vcf <- VariantAnnotation::readVcf(
    Rsamtools::TabixFile("Datasets/gnomAD/V2/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz"),
    param=vcfPar)
  
  # Discard variants not passing all FILTERS
  mask <- VariantAnnotation::filt(vcf) == "PASS"
  vcf <- vcf[mask, ]
  
  # prepare the variants freq data
  ref <- VariantAnnotation::ref(vcf) %>% as.character()
  alt <- VariantAnnotation::alt(vcf) %>% unlist() %>% as.character()
  pos <- MatrixGenerics::rowRanges(vcf) %>% BiocGenerics::start()
  af <- VariantAnnotation::info(vcf) %>% as.data.frame()
  
  # Generate the gnomAD pop. freq. table
  gnomAD_af_tbl <- cbind(data.frame(chr = chr, pos = pos, ref = ref, alt = alt), af)
  
  gnomAD_af_tbl

   }, mc.cores = parallel::detectCores() - 3)

germline_snv_popFreq_list <- do.call(rbind, germline_snv_popFreq_list)

############################3

## Filter out high gnomAD pop freq vars

### identify gnomAD vars that are present in germline_snv
germline_snv_gnomAD_pop_freq_index <- which(paste(germline_snv_popFreq_list$chr, 
                                                  germline_snv_popFreq_list$pos, 
                                                  germline_snv_popFreq_list$ref,
                                                  "_",
                                                  germline_snv_popFreq_list$alt, sep = "") %in% 
                                              paste(germline_snv$chr, 
                                                    germline_snv$pos, 
                                                    germline_snv$ref,
                                                    "_",
                                                    germline_snv$alt, sep = ""))

### identify germline_snv_gnomAD_pop_freq_index which have high frequency (FC > 0.01)
germline_snv_gnomAD_pop_freq_high_freq_index <- (((germline_snv_popFreq_list[germline_snv_gnomAD_pop_freq_index,c(5:ncol(germline_snv_popFreq_list))] > 0.01) %>% 
                                                    rowSums(na.rm = TRUE)) > 0)

### Filter out rare vars of germline_snv_gnomAD_pop_freq_index
germline_snv_gnomAD_pop_freq_index <- germline_snv_gnomAD_pop_freq_index[germline_snv_gnomAD_pop_freq_high_freq_index]

### identify index of high freq vars in germline_snv
germline_snv_gnomAD_index <- match(paste(germline_snv_popFreq_list$chr, 
                                         germline_snv_popFreq_list$pos, 
                                         germline_snv_popFreq_list$ref,
                                         germline_snv_popFreq_list$alt, sep = "_")[germline_snv_gnomAD_pop_freq_index],
                                   paste(germline_snv$chr, 
                                         germline_snv$pos, 
                                         germline_snv$ref,
                                         germline_snv$alt, sep = "_"))

# filter the germline_snv
germline_snv_flt <- germline_snv[-germline_snv_gnomAD_index,]

# remove the germline_snv_popFreq_list to free up space
rm(germline_snv_popFreq_list)

################################################

# Shell script for compressing the CADD data (adapted from https://gist.github.com/arq5x/9378020)
# We will perform the same pipeline on the InDel data as well to compresse it.

## 1. Download the raw CADD TSV and Tabix index (no annotations, just scores)

### https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz
### https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh37/whole_genome_SNVs.tsv.gz.tbi

##########

## 2. Add the required tools to the path (downloaded from https://vcf.iobio.io/help.html; https://vcf.iobio.io/tabix/TabixBinary.zip)
### export PATH=$PATH:/Users/asalavaty/OneDrive\ -\ Childrens\ Cancer\ Institute\ Australia/Apps\ and\ Libraries/TabixBinary/

##########

## 3. Have a look at the data and its size (remove `<` from after zcat when working within linux)

### ls -ltrh whole_genome_SNVs.tsv.gz
### zcat < whole_genome_SNVs.tsv.gz | head

##########

## 4. Do the following for compressing the file:
### Remove the 5th column (raw CADD scores)

# zcat < whole_genome_SNVs.tsv.gz \
# | cut -f 1,2,3,4,6 \
# | bgzip \
# > whole_genome_SNVs.compressed.tsv.gz

##########

## 5. As the original file is changed, we need to re-index the file

# tabix -p vcf -f whole_genome_SNVs.compressed.tsv.gz

##########

## 6. Move the Compressed directory to your desired location (OneDrive)

##########

######################

# Define a function for retrieving the CADD scores and reformat them to a proper and readable style for downstream analyses

scanCADD2df <-   function(variants_table, # Table of variants to check their CADD scores (Chromosome number/name on the first column, and POS, REFs and ALTs on the second third, and forth columns, respectively).
                          tabix_snv_file = NULL, # The tabix indexed file of SNVs (located in the same directory that .tbi file is located)
                          tabix_indel_file = NULL, # The tabix indexed file of InDels (located in the same directory that .tbi file is located)
                          gr_table, # A dataframe of genomic regions to scan; Chromosome number/name on the first column, and start and end positions on the second and third columns, respectively.
                          verbose = TRUE,
                          ...){
  
  # Check if a tabix file is provided
  if(is.null(tabix_snv_file) & is.null(tabix_indel_file)) {
    stop("Either of tabix_snv_file or tabix_indel_file, or both should be provided!")
  }
  
  library(magrittr)
  library(dplyr)

  if(verbose) {
    print(paste("There are ", nrow(gr_table), " chunks to be processed!", sep = ""))
  }
  
  # Make sure variants_table and are of data.frame class
  variants_table <- as.data.frame(variants_table)
  gr_table <- as.data.frame(gr_table)
  
  # Initialize a vector of variants
  variants <- paste(variants_table[,1], as.integer(variants_table[,2]), 
                    variants_table[,3], variants_table[,4], sep = "_")
  
  # Define regions for the seqminer package
  gr_table$region <- paste(gr_table[,1], ":", as.integer(gr_table[,2]), "-", as.integer(gr_table[,3]), sep = "")
  
  # Retrieve the tabix data based on gr_table
  cadd_df_list <-
    lapply(1:nrow(gr_table), function(i) {
      
      # Retrieve tabix SNV data
      tabix_snv_data <- if(is.null(tabix_snv_file)) {
        NULL
      } else {
        seqminer::tabix.read(tabixFile = tabix_snv_file, 
                             tabixRange = gr_table$region[i])
      }
      
      # Retrieve tabix InDel data
      tabix_indel_data <- if(is.null(tabix_indel_file)) {
        NULL
      } else {
        seqminer::tabix.read(tabixFile = tabix_indel_file, 
                             tabixRange = gr_table$region[i])
      }
      
      # Combine SNV and InDel results
      tabix_data <- c(tabix_snv_data, tabix_indel_data)
      
      # Prepare the tabix_data and re-format it
      tabix_data <- tabix_data %>% 
        stringr::str_split_fixed(pattern = "\t", n = 5) %>% 
        as.data.frame()
      
      # Set column names
      colnames(tabix_data) <- c("CHROM", "POS", "REF", "ALT", "CADD_PHRED")
      
      # Correct the class of POS and CADD_PHRED
      tabix_data$POS <- as.integer(tabix_data$POS)
      tabix_data$CADD_PHRED <- as.numeric(tabix_data$CADD_PHRED)
      
      # Match the variants in the tabix_data with variants vector generated from the variants_table
      tabix_data_variants <- paste(tabix_data[,1], as.integer(tabix_data[,2]),
                                   tabix_data[,3], tabix_data[,4], sep = "_")
      
      variants_index <- which(variants %in% tabix_data_variants)
      
      cadd_phred_scores <- tabix_data$CADD_PHRED[match(variants[variants_index], tabix_data_variants)]
      
      rm(tabix_data, tabix_data_variants)
      
      if(verbose) {
        print(paste("The chunk number ", i, " is done!", sep = ""))
      }
      
      return(data.frame(index = variants_index, score = cadd_phred_scores))
      
    })
  
  # Merge the lists of results
  cadd_df_list <- base::do.call(rbind, cadd_df_list)
  
  # Remove duplicates in case there are multiple duplicates of the same variant (e.g. corresponding to different tissues/cancers)
  cadd_df_list <- cadd_df_list %>% dplyr::distinct(index, .keep_all = TRUE)
  
  cadd_df_list
  
}

######################

# Define a function for creating intervals

interval_gen <- function(min, max, size) {
  
  # Generate sequence of intervals
  interval_seq <- seq(min, max, by = size)
  
  # Add the last item if the modulus of min-max by size does not equal zero
  if(!identical(max(interval_seq), max)) {
    interval_seq <- c(interval_seq, max)
  }
  
  # Generate the intervals table
  interval_tbl <- data.frame(min = vector(mode = "integer"), max = vector(mode = "integer"))
  
  lapply(0:(length(interval_seq) - 2), FUN = function(i) {
    
    if(i == 0) {
      interval_tbl <<- rbind(interval_tbl, c(interval_seq[i + 1], interval_seq[i + 2]))
    } else {
      interval_tbl <<- rbind(interval_tbl, c((interval_seq[i + 1] + 1), interval_seq[i + 2]))
    } 
    
  })
  
  colnames(interval_tbl) <- c("min", "max")
  
  return(interval_tbl)
}

######################

## Prepare germline_snv_chr_cadd_tbl table
genomic_intervals <- 50000

germline_snv_chr_cadd_tbl <- 
  lapply(1:length(germline_snv_chr), FUN = function(i) {
    
    ## Calculate the min and max of chromosome i
    germline_snv_chr_min <- min(germline_snv_flt$pos[which(germline_snv_flt$chr %in% germline_snv_chr[i])])
    germline_snv_chr_max <- max(germline_snv_flt$pos[which(germline_snv_flt$chr %in% germline_snv_chr[i])])
    
    ## Define the intervals for chromosome i
    if((germline_snv_chr_max - germline_snv_chr_min) <= genomic_intervals) {
      germline_snv_chr_cadd_tbl_tmp <- data.frame(min = germline_snv_chr_min, max = germline_snv_chr_max)
    } else {
      germline_snv_chr_cadd_tbl_tmp <- 
        interval_gen(min = germline_snv_chr_min, 
                     max = germline_snv_chr_max, 
                     size = genomic_intervals)
    }
    
    ## Add chromosome and refs columns to the table
    germline_snv_chr_cadd_tbl_tmp <- cbind(chr = germline_snv_chr[i],
                                           germline_snv_chr_cadd_tbl_tmp,
                                           refs = vector(mode = "logical", length = nrow(germline_snv_chr_cadd_tbl_tmp)))
    
    ## Define the refs for each interval
    sapply(as.list(1:nrow(germline_snv_chr_cadd_tbl_tmp)), function(j) {
      germline_snv_chr_cadd_tbl_tmp[j, 4] <<- ifelse(any(germline_snv_flt$chr %in% germline_snv_chr[i] &
                                                           germline_snv_flt$pos >= germline_snv_chr_cadd_tbl_tmp[j, 2] &
                                                           germline_snv_flt$pos < germline_snv_chr_cadd_tbl_tmp[j, 3]),
                                                     TRUE, FALSE)
      
      
    })
    
    germline_snv_chr_cadd_tbl_tmp$refs <- as.logical(germline_snv_chr_cadd_tbl_tmp$refs)
    
    germline_snv_chr_cadd_tbl_tmp
    
  })

## Merge the list of tables
germline_snv_chr_cadd_tbl <- do.call(rbind, germline_snv_chr_cadd_tbl)

## Filter out genomic regions that have no variants in them
germline_snv_chr_cadd_tbl <- germline_snv_chr_cadd_tbl[germline_snv_chr_cadd_tbl$refs,]

## Remove MT rows corresponding to mitochondrial chromosome
germline_snv_chr_cadd_tbl <- germline_snv_chr_cadd_tbl[-grep("MT", germline_snv_chr_cadd_tbl$chr),]

germline_snv_chr_cadd_tbl$refs <- NULL

## Retrieve CADD scores for all chromosomes of germline_snv (each interval of 50,000 last ~ 0.7-0.8 sec to be retrieved)

germline_snv_cadd_phred_scores <- scanCADD2df(variants_table = germline_snv_flt[,-1], 
                                      tabix_snv_file = "Datasets/CADD/v1.6-hg19 (compressed)/whole_genome_SNVs.compressed.tsv.gz",
                                      tabix_indel_file = "Datasets/CADD/v1.6-hg19 (compressed)/InDels.compressed.tsv.gz", 
                                      gr_table = germline_snv_chr_cadd_tbl)

germline_snv_cadd_phred_scores <- do.call(rbind, germline_snv_cadd_phred_scores)

## Add CADD scores to the germline_snv_flt
germline_snv_flt$cadd_phred[germline_snv_cadd_phred_scores$index] <- germline_snv_cadd_phred_scores$score

## Remove variants with Ref == Alt (These are wrong)
germline_snv_flt <- germline_snv_flt[which(germline_snv_flt$ref != germline_snv_flt$alt), ]

### Filter out intergenic variants
germline_snv_flt <- germline_snv_flt[-which(germline_snv_flt$consequence == "intergenic_variant"),]
dim(germline_snv_flt)
length(unique(germline_snv_flt$variant_id))
germline_snv_flt$sample_id %>% unique() %>% length()

germline_snv_flt_withCADD <- germline_snv_flt

# Retrieve missing CADD scores
germline_snv_cadd_phred_scores_missing_indices <- which(is.na(germline_snv_flt_withCADD$cadd_phred))

germline_snv_flt_withCADD_missed <- germline_snv_flt_withCADD[germline_snv_cadd_phred_scores_missing_indices,]

germline_snv_chr_cadd_tbl_missed <- germline_snv_flt_withCADD_missed[,c(2,3)]
colnames(germline_snv_chr_cadd_tbl_missed) <- c("chr", "min")
germline_snv_chr_cadd_tbl_missed$max <- germline_snv_chr_cadd_tbl_missed$min

germline_snv_cadd_phred_scores_missed <- scanCADD2df(variants_table = germline_snv_flt_withCADD_missed[,-1],
                                                     tabix_snv_file = "Datasets/CADD/v1.6-hg19 (compressed)/whole_genome_SNVs.compressed.tsv.gz",
                                                     tabix_indel_file = "Datasets/CADD/v1.6-hg19 (compressed)/InDels.compressed.tsv.gz",
                                                     gr_table = germline_snv_chr_cadd_tbl_missed)

germline_snv_cadd_phred_scores_missed <- do.call(rbind, germline_snv_cadd_phred_scores_missed)

germline_snv_flt_withCADD$cadd_phred[germline_snv_cadd_phred_scores_missed$index] <- germline_snv_cadd_phred_scores_missed$score

### Filter out SNV variants (length of ref and alt == 1) with NA CADD scores (these SNVs have wrong Refs in our data)
germline_snv_flt_withCADD <- germline_snv_flt_withCADD[-which(is.na(germline_snv_flt_withCADD$cadd_phred) & 
                                                                nchar(germline_snv_flt_withCADD$ref) == 1 & 
                                                                nchar(germline_snv_flt_withCADD$alt) == 1),]

#################################################

### Visualize CADD score distribution for getting an idea for imputation of the remaining NA CADD scores

vep_impact_tbl <- read.delim("Datasets/Ensembl_VEP_IMPACT_tbl.txt", sep = "\t")
vep_impact_tbl$IMPACT_score <- 1
vep_impact_tbl$IMPACT_score[vep_impact_tbl$IMPACT == "LOW"] <- 2
vep_impact_tbl$IMPACT_score[vep_impact_tbl$IMPACT == "MODERATE"] <- 3
vep_impact_tbl$IMPACT_score[vep_impact_tbl$IMPACT == "HIGH"] <- 4
 
germline_snv_flt_withCADD_noNA <- germline_snv_flt_withCADD[-which(is.na(germline_snv_flt_withCADD$cadd_phred)),]

germline_snv_flt_withCADD_noNA$consequence <-
sapply(germline_snv_flt_withCADD_noNA$consequence, function(i) {
  tmp_consequence <-
  if(length(grep("&|,", i)) > 0) {
    split_cons <- strsplit(x = i, split = "&|,") %>% unlist()
    split_cons_impact <- vep_impact_tbl$IMPACT_score[match(split_cons, vep_impact_tbl$SO.term)]
    
    split_cons_final <-
      split_cons[which(split_cons_impact == max(split_cons_impact))[1]]
    
    split_cons_final

  } else {
    i
  }
  tmp_consequence
}) %>% unlist()

#### Add Variant type to the table (SNV vs InDel)
germline_snv_flt_withCADD_noNA$Type <- "InDel"
germline_snv_flt_withCADD_noNA$Type[which(nchar(germline_snv_flt_withCADD_noNA$ref) == 1 & nchar(germline_snv_flt_withCADD_noNA$alt) == 1)] <- "SNV"

germline_snv_flt_withCADD_noNA_cons_plot <-
ggplot(germline_snv_flt_withCADD_noNA, aes(x = consequence, y = cadd_phred, fill = consequence)) +
  geom_boxplot() +
  facet_wrap(~Type) +
  scale_fill_discrete(name = "Consequence Type", type = fishualize::fish(n = 22)) +
  labs(x = NULL, y= "CADD Phred Score") +
  coord_flip() +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(filename = "Results/Germline SNV/germline_snv_flt_withCADD_noNA_cons_plot.pdf", 
       plot = germline_snv_flt_withCADD_noNA_cons_plot, device = "pdf", 
       width = 16, height = 16, units = "in")

#################################################

### Assign Ensembl VEP IMPACT type to our variants for imputing the CADD scores of variants with no CADD score

#### Add Variant type to the table (SNV vs InDel)
germline_snv_flt_withCADD$Type <- "InDel"
germline_snv_flt_withCADD$Type[which(nchar(germline_snv_flt_withCADD$ref) == 1 & nchar(germline_snv_flt_withCADD$alt) == 1)] <- "SNV"

#### Convert multiple consequence variants to single consequence by keeping only the consequence with the highest VEP impact
germline_snv_flt_withCADD$consequence <-
  sapply(germline_snv_flt_withCADD$consequence, function(i) {
    tmp_consequence <-
      if(length(grep("&|,", i)) > 0) {
        split_cons <- strsplit(x = i, split = "&|,") %>% unlist()
        split_cons_impact <- vep_impact_tbl$IMPACT_score[match(split_cons, vep_impact_tbl$SO.term)]
        
        split_cons_final <-
          split_cons[which(split_cons_impact == max(split_cons_impact))[1]]
        
        split_cons_final
        
      } else {
        i
      }
    tmp_consequence
  }) %>% unlist()

#### Remove SNV Type variants that belong to frameshift or synonymous (with high CADD phred scores) consequence class:
##### frameshift:their Alt is wrong and there are only 18 of them
##### synonymous: these variants strangely have high CADD phred scores (above 15) and there are only three of them
germline_snv_flt_withCADD <- germline_snv_flt_withCADD[-which(germline_snv_flt_withCADD$Type == "SNV" & 
                                                                germline_snv_flt_withCADD$consequence == "frameshift_variant" |
                                                                germline_snv_flt_withCADD$consequence == "synonymous_variant" &
                                                                germline_snv_flt_withCADD$cadd_phred > 15),]


germline_snv_flt_withCADD$VEP_IMPACT <- ""

for(i in c("MODIFIER", "LOW", "MODERATE", "HIGH")) {
  so_terms_index <- which(vep_impact_tbl$IMPACT == i)
  so_terms <- c(vep_impact_tbl$SO.term[so_terms_index], vep_impact_tbl$Display.term[so_terms_index])
  vep_impact_index <- grep(paste0(so_terms, collapse = "|"), germline_snv_flt_withCADD$consequence, ignore.case = TRUE)
  germline_snv_flt_withCADD$VEP_IMPACT[vep_impact_index] <- i
}

#### Impute the missing CADD scores

##### Calculate the median of CADD scores of each class of VEP Impacts
germline_VEP_HIGH_median <- median(germline_snv_flt_withCADD$cadd_phred[germline_snv_flt_withCADD$VEP_IMPACT == "HIGH"], na.rm = TRUE)
germline_VEP_MODERATE_median <- median(germline_snv_flt_withCADD$cadd_phred[germline_snv_flt_withCADD$VEP_IMPACT == "MODERATE"], na.rm = TRUE)
germline_VEP_LOW_median <- median(germline_snv_flt_withCADD$cadd_phred[germline_snv_flt_withCADD$VEP_IMPACT == "LOW"], na.rm = TRUE)
germline_VEP_MODIFIER_median <- median(germline_snv_flt_withCADD$cadd_phred[germline_snv_flt_withCADD$VEP_IMPACT == "MODIFIER"], na.rm = TRUE)

germline_snv_flt_withCADD$cadd_phred <- 
sapply(1:nrow(germline_snv_flt_withCADD), function(i) {
  imputed_score <- 
    if(!is.na(germline_snv_flt_withCADD$cadd_phred[i])) {
      germline_snv_flt_withCADD$cadd_phred[i]
    } else if(germline_snv_flt_withCADD$VEP_IMPACT[i] == "HIGH") {
      germline_VEP_HIGH_median
    } else if(germline_snv_flt_withCADD$VEP_IMPACT[i] == "MODERATE") {
      germline_VEP_MODERATE_median
    } else if(germline_snv_flt_withCADD$VEP_IMPACT[i] == "LOW") {
      germline_VEP_LOW_median
    } else if(germline_snv_flt_withCADD$VEP_IMPACT[i] == "MODIFIER") {
      germline_VEP_MODIFIER_median
    }
})

#################################################

## Calibrating the CADD scores

### First define the binary labels (pathogenic (1) vs benign (0)) based on clinVar data for training a calibration model

#### Download the clinVar data file (variant_summary.txt.gz)
download.file(url = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz",
              destfile = "~/Downloads/variant_summary.txt.gz")

#### Read in the file
clinVar_summary_con <- gzfile("~/Downloads/variant_summary.txt.gz")
clinVar_summary <- readr::read_delim(clinVar_summary_con, delim = "\t")
rm(clinVar_summary_con)

#### Remove redundant columns
#### and retain only rows with Benign and Pathogenic ClinicalSignificance
#### and remove rows with unknown Position/Ref/Alt ("na" ReferenceAlleleVCF)
clinVar_summary <- clinVar_summary %>% select(Type, Name, GeneSymbol, 
                                              ClinicalSignificance, ReviewStatus, Assembly, 
                                              Chromosome, PositionVCF, 
                                              ReferenceAlleleVCF, AlternateAlleleVCF) %>% 
  dplyr::filter(grepl("practice guideline|reviewed by expert panel", ReviewStatus)) %>% 
  dplyr::filter(grepl('^Benign$|^Pathogenic$', ClinicalSignificance)) %>% 
  subset(ReferenceAlleleVCF != "na")

#### Add coding/non-coding class to the table
clinVar_summary$coding_type <- "noncoding"
clinVar_summary$coding_type[grep(" \\(p\\.", clinVar_summary$Name)] <- "coding"

#### Save the dataset for later use
clinVar_summary_con <-  gzfile("Datasets/clinVar_summary.txt.gz")
readr::write_delim(clinVar_summary, clinVar_summary_con, delim = "\t")
rm(clinVar_summary_con)

#################

clinVar_summary_con <-  gzfile("Datasets/clinVar_summary.txt.gz")
clinVar_summary <- readr::read_delim(clinVar_summary_con, delim = "\t")
rm(clinVar_summary_con)

### Retain only hg19 for our current project
clinVar_summary <- subset(clinVar_summary, Assembly == "GRCh37")

### Remove germline and somatic variants (testing/calibrating set) from clinVar_summary to make it ready for training
clinVar_summary_omit_index <- which(c(paste(clinVar_summary$Chromosome, as.integer(clinVar_summary$PositionVCF),
                                            clinVar_summary$ReferenceAlleleVCF, clinVar_summary$AlternateAlleleVCF, sep = "_")) %in% 
                                      c(paste(germline_snv_flt_withCADD$chr, as.integer(germline_snv_flt_withCADD$pos), 
                                            germline_snv_flt_withCADD$ref, germline_snv_flt_withCADD$alt, sep = "_"),
                                      paste(somatic_snv_flt_withCADD$chr, as.integer(somatic_snv_flt_withCADD$pos), 
                                            somatic_snv_flt_withCADD$ref, somatic_snv_flt_withCADD$alt, sep = "_")))


clinVar_summary <- clinVar_summary[-clinVar_summary_omit_index,]

### Randomly select 1500 pathogenic (some of them will not have CADD scores) and 1500 benign variants
set.seed(1234)
clinVar_summary_pathogenic <- clinVar_summary[sample(which(clinVar_summary$ClinicalSignificance == "Pathogenic"), size = 1500),]

set.seed(1234)
clinVar_summary_benign <- clinVar_summary[sample(which(clinVar_summary$ClinicalSignificance == "Benign"), size = 1500),]

clinVar_summary_for_training <- rbind(clinVar_summary_pathogenic, clinVar_summary_benign)
rm(clinVar_summary_pathogenic, clinVar_summary_benign)

### Get the CADD scores of the training set
clinVar_summary_for_training$ClinicalSigLabel <- 0
clinVar_summary_for_training$ClinicalSigLabel[clinVar_summary_for_training$ClinicalSignificance == "Pathogenic"] <- 1

#### Shuffle the rows
clinVar_summary_for_training <- clinVar_summary_for_training[sample(nrow(clinVar_summary_for_training)),]

#### Prepare clinVar_summary_chr_cadd_tbl table
clinVar_summary_chr_cadd_tbl <- clinVar_summary_for_training[,c(7,8)]
colnames(clinVar_summary_chr_cadd_tbl) <- c("chr", "min")
clinVar_summary_chr_cadd_tbl$max <- clinVar_summary_chr_cadd_tbl$min

## Retrieve CADD scores for all chromosomes of clinVar_summary

clinVar_summary_cadd_phred_scores <- scanCADD2df(variants_table = clinVar_summary_for_training[,c(7:10)], 
                                                 tabix_snv_file = "Datasets/CADD/v1.6-hg19 (compressed)/whole_genome_SNVs.compressed.tsv.gz",
                                                 tabix_indel_file = "Datasets/CADD/v1.6-hg19 (compressed)/InDels.compressed.tsv.gz", 
                                                 gr_table = clinVar_summary_chr_cadd_tbl)

## Add cadd data to the table
clinVar_summary_for_training$cadd_phred[clinVar_summary_cadd_phred_scores$index] <- clinVar_summary_cadd_phred_scores$score

# Retrieve missing CADD scores
clinVar_summary_cadd_phred_scores_missing_indices <- which(is.na(clinVar_summary_for_training$cadd_phred))

clinVar_summary_for_training_missed <- clinVar_summary_for_training[clinVar_summary_cadd_phred_scores_missing_indices,]

clinVar_summary_chr_cadd_tbl_missed <- clinVar_summary_for_training_missed[,c(7,8)]
colnames(clinVar_summary_chr_cadd_tbl_missed) <- c("chr", "min")
clinVar_summary_chr_cadd_tbl_missed$max <- clinVar_summary_chr_cadd_tbl_missed$min

clinVar_summary_cadd_phred_scores_missed <- scanCADD2df(variants_table = clinVar_summary_for_training_missed[,c(7:10)],
                                                        tabix_snv_file = "Datasets/CADD/v1.6-hg19 (compressed)/whole_genome_SNVs.compressed.tsv.gz",
                                                        tabix_indel_file = "Datasets/CADD/v1.6-hg19 (compressed)/InDels.compressed.tsv.gz",
                                                        gr_table = clinVar_summary_chr_cadd_tbl_missed)


## Remove the variants with NA cadd scores
clinVar_summary_for_training <- clinVar_summary_for_training[-which(is.na(clinVar_summary_for_training$cadd_phred)),]

## Get the number of remaining Benign and Pathogenic variants
table(clinVar_summary_for_training$ClinicalSignificance) 
#Benign Pathogenic 
#1485       1443 

## Equalize the number of Benign and Pathogenic
diff_ClinicalSignificance <- table(clinVar_summary_for_training$ClinicalSignificance) %>% sort() %>% diff()
set.seed(1234)
clinVar_summary_for_training <- clinVar_summary_for_training[-sample(which(clinVar_summary_for_training$ClinicalSignificance == "Benign"), size = diff_ClinicalSignificance),]

#####################

## Calibrate using GAM

### Visualize coding vs non-coding
library(mgcv)
clinVar_summary_for_training_gam <- gam(data = clinVar_summary_for_training,
                                        formula = ClinicalSigLabel ~ s(cadd_phred, by=as.factor(coding_type)), 
                                        family=binomial())

plot(clinVar_summary_for_training_gam,pages=1,residuals=TRUE)  ## show partial residuals

### Train for coding and noncoding together
clinVar_summary_for_training_gam <- gam(data = clinVar_summary_for_training,
                                        formula = ClinicalSigLabel ~ s(cadd_phred), 
                                        family=binomial())

### We can either use the type="link" to get (-Inf,Inf) predictions then use logistic regression to transform it to (0,1) probabilities or use type="response" to directly get (0,1) probabilities.
#### with type="response" predictions on the scale of the response are returned. So, type="link" would be redundant.
# germline_snv_flt_withCADD_gam_prob <- 1/(1+exp(-predict(clinVar_summary_for_training_gam, 
#                                                            newdata = germline_snv_flt_withCADD[,"cadd_phred"], 
#                                                            type="link")))

### get the gam probabilities using the logistic transformation
germline_snv_flt_withCADD_gam_prob <- predict(clinVar_summary_for_training_gam, 
                                              newdata = germline_snv_flt_withCADD[,"cadd_phred"], 
                                              type="response")

### Visualize calibration probabilities vs original cadd
germline_snv_flt_withCADD$calibCADD_GAM_prob <- as.numeric(germline_snv_flt_withCADD_gam_prob)

germline_snv_flt_GAM_probvsPhred_plot <-
  ggplot(germline_snv_flt_withCADD, aes(x = cadd_phred, y = calibCADD_GAM_prob, color = VEP_IMPACT)) +
  geom_point() +
  scale_color_discrete(name = "Ensembl\nVEP IMPACT", type = fishualize::fish(n = 4)) +
  labs(x = "CADD Phred Score", y= "GAM-based Calibrated CADD Score (Pathogenicity Probability)") +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

germline_snv_flt_GAM_probvsPhred_plot

ggsave(filename = "Results/Germline SNV/ML-based calibration/GAM/germline_snv_flt_GAM_probvsPhred_plot.pdf", 
       plot = germline_snv_flt_GAM_probvsPhred_plot, 
       device = "pdf", width = 8, height = 8, units = "in")

# Log transform GAM probabilities
germline_snv_flt_withCADD$calibCADD_GAM_prob_log <- (-10)*log10(1-germline_snv_flt_withCADD$calibCADD_GAM_prob)

###################

# Define Min-Max (range) normalization
range_normalize <- function(data, min = 0, max = 1, na.rm = TRUE, ...) {
  
  min+(((data-min(data, na.rm = na.rm, ...))*(max-min))/
         (max(data, na.rm = na.rm, ...)-min(data, na.rm = na.rm, ...)))
  
}

#################################################

## Retain only PRISM study data
germline_snv_flt_withCADD <- subset(germline_snv_flt_withCADD, study == "PRISM")

## Calculate nalt for all variants
germline_snv_flt_withCADD_nalt <- lapply(unique(germline_snv_flt_withCADD$variant_id), function(i) {
  tmp_data <- subset(germline_snv_flt_withCADD, variant_id == i)
  ref <- tmp_data$ref[1]
  genotypes <- tmp_data$genotype %>% 
    stringr::str_split(pattern = "/", n = 2) %>% unlist()
  count <- length(genotypes) - grep(pattern = paste0("^", ref, "$"), x = genotypes) %>% length()
  tmp_tbl <- data.frame(variant_id = i, nalt = count)
})

germline_snv_flt_withCADD_nalt <- do.call(rbind, germline_snv_flt_withCADD_nalt)

germline_snv_flt_withCADD$nalt <- germline_snv_flt_withCADD_nalt$nalt[match(germline_snv_flt_withCADD$variant_id, germline_snv_flt_withCADD_nalt$variant_id)]

######################################################################################
######################################################################################

# Evaluating whether the variants that have recurrences fall within the bad genomic regions or not

library(rtracklayer)
library(GenomicRanges)
library(parallel)

## Import bad (difficult to sequence regions)
GRCh37_difficult_regions <- rtracklayer::import("Datasets/GRCh37_alldifficultregions.bed.gz")

germline_query_granges <- GenomicRanges::GRanges(seqnames = germline_snv_flt_withCADD$chr, 
                                                 ranges = IRanges(germline_snv_flt_withCADD$pos, width = 1))

germline_bad_region <- GenomicRanges::findOverlaps(germline_query_granges, GRCh37_difficult_regions) %>% 
  as.data.frame() %>% select(queryHits) %>% unlist() %>% unname()

### Add to the dataset whether each variant is within a bad region or not

germline_snv_flt_withCADD$bad_region <- FALSE
germline_snv_flt_withCADD$bad_region[germline_bad_region] <- TRUE

### Add to the dataset whether each variant is recurrent or not

germline_recurrrent_vars <- # name of the vector will be the name of recurrent variants and the value of each element will be its number of recurrence
  germline_snv_flt_withCADD %>% 
  janitor::get_dupes(variant_id) %>% 
  select(variant_id) %>% 
  unlist() %>% 
  unname() %>% 
  table()

germline_snv_flt_withCADD <- 
  germline_snv_flt_withCADD %>% 
  mutate(recurrent = ifelse(variant_id %in% names(germline_recurrrent_vars), TRUE, FALSE))

#### Add the number of recurrence
germline_snv_flt_withCADD$n_recurrence <- 1

lapply(1:length(germline_recurrrent_vars), function(i) {
  germline_snv_flt_withCADD$n_recurrence[germline_snv_flt_withCADD$variant_id == names(germline_recurrrent_vars)[i]] <<- 
    germline_recurrrent_vars[i]
  return(NULL)
})

germline_snv_flt_withCADD$recurrent[germline_snv_flt_withCADD$recurrent == TRUE] <- "Recurrent"
germline_snv_flt_withCADD$recurrent[germline_snv_flt_withCADD$recurrent == FALSE] <- "Non-recurrent"

germline_snv_flt_withCADD$bad_region[germline_snv_flt_withCADD$bad_region == TRUE] <- "Difficult Region"
germline_snv_flt_withCADD$bad_region[germline_snv_flt_withCADD$bad_region == FALSE] <- "Non-difficult Region"

### Visualizing the results

germline_snv_bad_region_vis1 <- ggplot(data = germline_snv_flt_withCADD, aes(x = recurrent, fill = bad_region)) +
  geom_bar(position = "dodge") +
  labs(fill = "Region Quality Type") +
  ylab("Number of variants") +
  xlab(NULL) +
  fishualize::scale_fill_fish_d() +
  theme_classic()

germline_snv_bad_region_vis2 <- ggplot(data = germline_snv_flt_withCADD, aes(x = recurrent, fill = bad_region)) +
  geom_bar(position = "fill") + 
  ylab("Proportion") +
  labs(fill = "Region Quality Type") +
  xlab(NULL) +
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_fill(vjust=0.5), colour="white") +
  fishualize::scale_fill_fish_d() +
  theme_classic()

######################################################################################
######################################################################################

# filtering out variants seen in > 5% of the cohort (filtering out platform errors)
## Although it is conceivable that real biological signal might be present at 5% rate within a cohort (specific cancer type), it's less likely that such signal would be present at 5% across all cancers. But platform error will be present at the same rate in all.

### Save the data that include the highly recurrent variants as well
saveRDS(object = germline_snv_flt_withCADD, file = "Germline data with recurrences/germline_snv_flt_withCADD.rds")

### Also keep as an object
germline_snv_withRecurr_withCADD <- germline_snv_flt_withCADD

germline_snv_flt_withCADD <- subset(germline_snv_flt_withCADD, n_recurrence <= round(0.05*length(unique(germline_snv_flt_withCADD$patient_id))))

## Plot that highly frequent variants are removed.
table(janitor::get_dupes(germline_snv_flt_withCADD, variant_id)$variant_id) %>% table() %>% plot()

### Visualize the results

germline_snv_recurrflt_bad_region_vis1 <- ggplot(data = germline_snv_flt_withCADD, aes(x = recurrent, fill = bad_region)) +
  geom_bar(position = "dodge") +
  labs(fill = "Region Quality Type") +
  ylab("Number of variants") +
  xlab(NULL) +
  fishualize::scale_fill_fish_d() +
  theme_classic()

germline_snv_recurrflt_bad_region_vis2 <- ggplot(data = germline_snv_flt_withCADD, aes(x = recurrent, fill = bad_region)) +
  geom_bar(position = "fill") + 
  ylab("Proportion") +
  labs(fill = "Region Quality Type") +
  xlab(NULL) +
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_fill(vjust=0.5), colour="white") +
  fishualize::scale_fill_fish_d() +
  theme_classic()

### We see that Generally, the recurrent variants are mostly present in the "Difficult regions". 
### After filtration, we still have recurrent variants in the difficult regions but the proportion of "Difficult region"- vs "Non-difficult regions"-recurrent variants has changed after the filtration.

######################################################################################
######################################################################################

# checking the region quality type of ribosomal protein genes

ribo_p_pattern <- "^RPS\\d+$|^RPL\\d+$"

germline_snv_ribosomal_genes <- germline_snv_flt_withCADD$gene[stringr::str_detect(germline_snv_flt_withCADD$gene, ribo_p_pattern)]

ggplot(data = germline_snv_flt_withCADD %>% filter(gene %in% germline_snv_ribosomal_genes), aes(x = recurrent, fill = bad_region)) +
  geom_bar(position = "fill") + 
  ylab("Proportion") +
  labs(fill = "Region Quality Type") +
  xlab(NULL) +
  stat_count(geom = "text", 
             aes(label = stat(count)),
             position=position_fill(vjust=0.5), colour="white") +
  fishualize::scale_fill_fish_d() +
  theme_classic()

germline_snv_flt_withCADD %>% filter(gene %in% germline_snv_ribosomal_genes) %>% filter(bad_region == "Non-difficult Region") %>% select(variant_id, chr, pos, ref, alt, gene, transcript, cadd_phred, bad_region, n_recurrence) %>% View()

################################

## check if the ribosomal protein gene variants fall within repeatitive regions
library(rtracklayer)
library(GenomicRanges)
library(UCSCRepeatMasker)
library(AnnotationHub)

germline_snv_flt_withCADD_ribo <- germline_snv_flt_withCADD %>% filter(gene %in% germline_snv_ribosomal_genes)

### Load AnnotationHub
ah <- AnnotationHub()

### find the reference version you wanna use
query(ah, c("RepeatMasker", "Homo sapiens"))

### retrieve RepeatMaskers
rmskhg19 <- ah[["AH99002"]]

### create the GRanges of germline_snv_flt_withCADD_ribo
germline_snv_ribo_granges <- GenomicRanges::GRanges(seqnames = paste("chr", germline_snv_flt_withCADD_ribo$chr, sep = ""), 
                                                 ranges = IRanges(germline_snv_flt_withCADD_ribo$pos, width = 1))

### find the overlap between RepeatMaskers and our variants
germline_snv_ribo_repeats <- GenomicRanges::findOverlaps(germline_snv_ribo_granges, rmskhg19) %>% 
  as.data.frame() %>% select(queryHits) %>% unlist() %>% unname()

### Add to the dataset whether each variant is within a bad region or not

germline_snv_flt_withCADD_ribo$repeat_region <- FALSE
germline_snv_flt_withCADD_ribo$repeat_region[germline_snv_ribo_repeats] <- TRUE

germline_snv_flt_withCADD_ribo$repeat_region[germline_snv_flt_withCADD_ribo$repeat_region == TRUE] <- "Repeat Region"
germline_snv_flt_withCADD_ribo$repeat_region[germline_snv_flt_withCADD_ribo$repeat_region == FALSE] <- "Non-repeat Region"

### Visualizing the results

germline_snv_flt_withCADD_ribo_vis1 <- ggplot(data = germline_snv_flt_withCADD_ribo, aes(x = recurrent, fill = repeat_region)) +
  geom_bar(position = "dodge") +
  labs(fill = "Region Type") +
  ylab("Number of variants") +
  xlab(NULL) +
  fishualize::scale_fill_fish_d(begin = 0.5) +
  theme_classic()

######################################################################################
######################################################################################

# # removing ribosomal protein, HLA, and olfactory receptor genes from our data

ribo_p_pattern <- "^RPS\\d+$|^RPL\\d+$"
or_pattern <- "^OR\\w+"
hla_pattern <- "^HLA-|^DRB"

germline_snv_flt_withCADD <- germline_snv_flt_withCADD[-which(stringr::str_detect(string = germline_snv_flt_withCADD$gene, 
                                                                           pattern = paste(ribo_p_pattern, or_pattern, hla_pattern, sep = "|"))),]

######################################################################################
######################################################################################

## Summarize the data by summing the calibrated CADD scores of each gene within each sample

### Summarize
germline_snv_flt_withCADD_summarized <- germline_snv_flt_withCADD %>% 
  group_by(chr, gene, cancer_type, study, sample_id, patient_id) %>% 
  dplyr::summarise(cadd_phred = sum(calibCADD_GAM_prob))

######################

### Re-transform the Phred/log scaled pathogenicities to linear scale probabilities (Not required as we used calibCADD_GAM_prob not calibCADD_GAM_prob_log)
# germline_snv_flt_withCADD_summarized$cadd_phred <- 1-10^(germline_snv_flt_withCADD_summarized$cadd_phred/-10)

######################

## We should create a list of germline SNV data with different cancer types and perform subsequent analyses separately for each cancer type
germline_snv_flt_withCADD_summarized <- lapply(unique(germline_snv_flt_withCADD_summarized$cancer_type), function(i) {
  germline_snv_flt_withCADD_summarized[which(germline_snv_flt_withCADD_summarized$cancer_type == i),]
})

### Set the names of list
germline_snv_flt_cancer_types <- sapply(germline_snv_flt_withCADD_summarized, function(i) unique(i$cancer_type))
names(germline_snv_flt_withCADD_summarized) <- germline_snv_flt_cancer_types

## Prepare tables in the form of genes on columns and samples on rows for corr analysis
germline_snv_flt_withCADD_summarized <- lapply(germline_snv_flt_withCADD_summarized, function(i) {
  tmp.data <- i[, c(2,5,7)]
  tmp.data <- tidyr::pivot_wider(data = tmp.data, 
                                 id_cols = "sample_id", 
                                 names_from = "gene", 
                                 values_from = "cadd_phred",
                                 values_fn = mean) %>% as.data.frame()
  
  rownames(tmp.data) <- tmp.data$sample_id
  tmp.data <- tmp.data[,-1]
  
  # Replace NA with 0
  tmp.data <- apply(tmp.data, 1, function(i) {i[is.na(i)]  <- 0; i}) %>% t() %>% as.data.frame()
  
  # Report the samples that their number of corresponding genes with data is less than < 1 percent of the total number of genes in the table of its corresponding cancer type
  row_index <- apply(tmp.data, 1, function(i) {sum(i > 0)})/ncol(tmp.data) < 0.01
  if(sum(row_index) > 0) {
    tmp_text <- paste("\nThe following samples have mutations in very few number of genes. You may further check their quality.\n",
                      paste0(rownames(tmp.data)[row_index], collapse = "\n"), sep = "")
    cat(tmp_text, sep = "")
  }

  # Remove noise genes
  tmp.data <- tmp.data[,colSums(tmp.data > 0.01) != 0]
  
  tmp.data
})

# Remove tables with less than 3 samples
germline_snv_flt_withCADD_summarized[which(sapply(germline_snv_flt_withCADD_summarized, nrow) < 3)] <- NULL

################################################################
################################################################

# Inspect the statistics of Deleteriousness data

## Mean Del across all genes in each cancer
sapply(1:length(germline_snv_flt_withCADD_summarized), function(i) {
  cat(paste0(names(germline_snv_flt_withCADD_summarized)[i], ": ", mean(unlist(germline_snv_flt_withCADD_summarized[[i]]))))
  cat("\n********************\n")
})

findoutlier <- function(x) {
  return(x < quantile(x, .25) - 1.5*IQR(x) | x > quantile(x, .75) + 1.5*IQR(x))
}

## Mean Del for each gene (trends)
sapply(1:length(germline_snv_flt_withCADD_summarized), function(i) {
  
  tmp_data = apply(germline_snv_flt_withCADD_summarized[[i]], 2, mean)
  
  png(filename = paste0("Results/Germline Deleteriousness stats/0.95/", names(germline_snv_flt_withCADD_summarized)[i], ".png"))
  
  plot(y = tmp_data, 
       x = 1:ncol(germline_snv_flt_withCADD_summarized[[i]]),
       xlab = "Genes", ylab = "Mean Del Score", 
       main = names(germline_snv_flt_withCADD_summarized)[i])
  
  qt_tmp = quantile(tmp_data, probs=c(.25, .95), na.rm = FALSE)
  iqr_tmp = IQR(tmp_data)
  lower_tmp = qt_tmp[1] - 1.5*iqr_tmp
  upper_tmp = qt_tmp[2] + 1.5*iqr_tmp
  
    abline(h = c(lower_tmp, upper_tmp), col = "red", lwd = 2)
    
    dev.off()
    
})

## Number of outlier genes based on IQR of 0.25 and 0.75
sapply(1:length(germline_snv_flt_withCADD_summarized), function(i) {
  
  tmp_data = apply(germline_snv_flt_withCADD_summarized[[i]], 2, mean)
  
  quartiles = quantile(tmp_data, probs=c(.25, .75), na.rm = FALSE)
  
  IQR = IQR(tmp_data)
  
  Lower = quartiles[1] - 1.5*IQR
  Upper = quartiles[2] + 1.5*IQR
  
  cat(paste0("Number of genes in ", names(germline_snv_flt_withCADD_summarized)[i], ": ", ncol(germline_snv_flt_withCADD_summarized[[i]])))
  cat("\n")
  cat(paste0("Number of outliers: ", length(which(tmp_data < Lower | tmp_data > Upper))))
  cat("\n********************\n\n")
})

############

# Visualize all cancers together

germline_snv_Del_data_4stats <-
  lapply(1:length(germline_snv_flt_withCADD_summarized), function(i) {
    
    tmp_data = data.frame(Del = apply(germline_snv_flt_withCADD_summarized[[i]], 2, mean))
    tmp_data = cbind(Gene = rownames(tmp_data), tmp_data)
    rownames(tmp_data) = NULL
    tmp_data = tmp_data %>% mutate(outlier = ifelse(findoutlier(Del), Gene, NA), 
                                   Cancer_Type = names(germline_snv_flt_withCADD_summarized)[i])
    
    return(tmp_data)

})

germline_snv_Del_data_4stats <- do.call(rbind, germline_snv_Del_data_4stats)

germline_snv_Del_data_4stats_plot <- 
ggplot(germline_snv_Del_data_4stats, aes(x = Cancer_Type, y = Del, fill = Cancer_Type)) + 
  geom_boxplot(show.legend = F, outlier.size = 0.1) +
  theme_classic() + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), legend.position = "none") +
  xlab("Cancer Type") +
  ylab("Mean Del Score") +
  facet_wrap(facets = vars(Cancer_Type), nrow = 3) +
  # scale_fill_manual(values = )
  geom_text_repel(aes(label=outlier), na.rm=TRUE, hjust=-.5, max.overlaps = 15,
                   size = 2, segment.alpha = 0.3)

ggsave(filename = "Results/Germline Deleteriousness stats/germline_snv_Del_data_4stats_plot.pdf", 
       plot = germline_snv_Del_data_4stats_plot, device = "pdf", width = 17, height = 22, units = "cm")

################################################################

## Perform correlation analysis

germline_snv_flt_cor <- lapply(germline_snv_flt_withCADD_summarized, function(i){
  tmp_cor <- fcor(i)
  tmp_cor <- subset(tmp_cor, mr < 20)
  tmp_cor
})

##################

## Create a sample sub-graph of co-del data for visualization

library(influential)
library(igraph)
library(visNetwork)

hgg_germline_subnet <- influential::graph_from_data_frame(germline_snv_flt_cor$HGG, directed = F)
hgg_germline_subnet <- igraph::subgraph(graph = hgg_germline_subnet, vids = V(hgg_germline_subnet)[1:2000])
hgg_germline_subnets_genes <- igraph::cliques(graph = hgg_germline_subnet, min = 15, max = 20)
hgg_germline_subnet_final <- subgraph(hgg_germline_subnet, vids = V(hgg_germline_subnet)[hgg_germline_subnets_genes[[6]]])

hgg_germline_subnet_final_vis <-
visIgraph(igraph = hgg_germline_subnet_final) %>% 
  visOptions(
    # nodesIdSelection = TRUE, 
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '30px arial black',
    size = 25,
    # color = "orange",
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "circle" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
    ) %>% 
  visEdges(color = "lightgrey")

visExport(graph = hgg_germline_subnet_final_vis, 
          name = "hgg_germline_subnet_final_vis.png", 
          label = "Subset of HGG co-deleteriousness network")



######################## For pan cancer analysis ######################## 

## Prepare table in the form of genes on columns and samples on rows
pan_germline_snv_tbl <- germline_snv_flt_withCADD %>% 
  group_by(chr, gene, cancer_type, study, sample_id, patient_id) %>% 
  dplyr::summarise(cadd_phred = sum(calibCADD_GAM_prob))

pan_germline_snv_n_cancer <- pan_germline_snv_tbl$cancer_type %>% unique() %>% length()

pan_germline_snv_tbl <- pan_germline_snv_tbl[,c(2,5,7)]
pan_germline_snv_tbl <- tidyr::pivot_wider(data = pan_germline_snv_tbl, 
                                           id_cols = "sample_id", 
                                           names_from = "gene", 
                                           values_from = "cadd_phred",
                                           values_fn = mean) %>% as.data.frame()

rownames(pan_germline_snv_tbl) <- pan_germline_snv_tbl$sample_id
pan_germline_snv_tbl <- pan_germline_snv_tbl[,-1]

# Replace NA with 0
pan_germline_snv_tbl <- apply(pan_germline_snv_tbl, 1, function(i) {i[is.na(i)]  <- 0; i}) %>% t() %>% as.data.frame()

# Report the samples that their number of corresponding genes with data is less than < 1 percent of the total number of genes in the table
pan_germline_snv_tbl_row_index <- apply(pan_germline_snv_tbl, 1, function(i) {sum(i > 0)})/ncol(pan_germline_snv_tbl) < 0.01

if(sum(pan_germline_snv_tbl_row_index) > 0) {
  tmp_text <- paste("\nThe following samples have mutations in very few number of genes. You may further check their quality.\n",
                    paste0(rownames(pan_germline_snv_tbl)[pan_germline_snv_tbl_row_index], collapse = "\n"), sep = "")
  cat(tmp_text, sep = "")
}

# Remove noise genes
pan_germline_snv_tbl <- pan_germline_snv_tbl[,colSums(pan_germline_snv_tbl > 0.01) != 0]

# Perform corr analysis
pan_germline_snv_cor <-  fcor(data = pan_germline_snv_tbl, 
                              method = "spearman", 
                              mutualRank = TRUE, pvalue = FALSE, flat = TRUE)

# Filter top 20 (MR < 20) correlations
pan_germline_snv_cor <- subset(pan_germline_snv_cor, mr < 20)

#######################################################################

#=============================================================================
#
#    Code chunk 4: Get and prepare the somatic SNV data
#
#=============================================================================

### Somatic Variant

somatic_snv <- tbl(db, "lkcgp_curated_somatic_snv") %>% as.data.frame()

sample_somatic_snv <- tbl(db, "lkcgp_curated_sample_somatic_snv") %>% as.data.frame()

#### Remove exon (redundant)
somatic_snv$exon <- NULL

### Complement the dataset
somatic_snv <- somatic_snv[c(somatic_snv$variant_id %in% sample_somatic_snv$variant_id),]

sample_somatic_snv <- sample_somatic_snv[,c(1:12)]
sample_somatic_snv <- cbind(somatic_snv[match(sample_somatic_snv$variant_id, somatic_snv$variant_id),],
                            sample_somatic_snv[,-3])

somatic_snv <- sample_somatic_snv
rm(sample_somatic_snv)

#### Filter variants that are NOT present in lkcgp_sample$sample_id (this table is used to extract the patient ID and cancer type)
somatic_snv <- somatic_snv[c(somatic_snv$sample_id %in% lkcgp_sample$sample_id),]

#### Add patient ID, cancer type, Study, and lastUpdated date
somatic_snv$patient_id <- lkcgp_sample$patient_id[match(somatic_snv$sample_id, lkcgp_sample$sample_id)]
somatic_snv$cancer_type <- lkcgp_sample$cancer_type[match(somatic_snv$sample_id, lkcgp_sample$sample_id)]
somatic_snv$study <- lkcgp_sample$study[match(somatic_snv$sample_id, lkcgp_sample$sample_id)]
somatic_snv$lastUpdated <- lkcgp_sample$lastUpdated[match(somatic_snv$sample_id, lkcgp_sample$sample_id)]

## Remove variants with 0 or NA depth
somatic_snv <- somatic_snv[-c(which(somatic_snv$depth == 0 | is.na(somatic_snv$depth))),]

## Remove variants with <= 3 altAD or NA altAD
somatic_snv <- somatic_snv[-c(which(somatic_snv$altad <= 3 | is.na(somatic_snv$altad))),]

#### Add Other as the cancer_type of variants that correspond to NA cancer_type
somatic_snv[c(which(is.na(somatic_snv$cancer_type))), "cancer_type"] <- "Other"

#### Remove redundant columns
somatic_snv$copynumber <- NULL
somatic_snv$sjc_medal <- NULL

################################################
################################################

## Get the population frequencies of somatic SNVs

# Process downloaded gnomAD pop freq. data

### First, get min and max genomic positions of all chromosomes and create batches of 100,000 bps (genomic positions) for each chromosome
somatic_snv_chr <- 
  unique(somatic_snv$chr)

somatic_snv_chr_tbl <-
parallel::mclapply(1:length(somatic_snv_chr), FUN = function(i) {
  
  ## Calculate the min and max of chromosome i
  somatic_snv_chr_min <- min(somatic_snv$pos[which(somatic_snv$chr %in% somatic_snv_chr[i])]) - 10
  somatic_snv_chr_max <- max(somatic_snv$pos[which(somatic_snv$chr %in% somatic_snv_chr[i])]) + 10
  
  somatic_snv_chr_tbl_tmp <- data.frame(chr = somatic_snv_chr[i], 
                                        min = somatic_snv_chr_min, 
                                        max = somatic_snv_chr_max)
  
}, mc.cores = parallel::detectCores() - 3)

somatic_snv_chr_tbl <- do.call(rbind, somatic_snv_chr_tbl)

############################

## Filter out chromosomes of somatic_snv_chr_tbl that are not available in the genomes.r2.1.1.exome_calling_header
somatic_snv_chr_gnomAD_tbl <- somatic_snv_chr_tbl[which(somatic_snv_chr_tbl$chr %in% genomes.r2.1.1.exome_calling_chrs),]

## Remove Y chromosome as it is not present in the Tabix index file
somatic_snv_chr_gnomAD_tbl <- somatic_snv_chr_gnomAD_tbl[-which(somatic_snv_chr_gnomAD_tbl$chr == "Y"),]

somatic_snv_popFreq_list <- parallel::mclapply(1:nrow(somatic_snv_chr_gnomAD_tbl), function(i) {
  
  ## Define a genomic ranges
  chr = somatic_snv_chr_gnomAD_tbl[i,1]
  start.loc <- somatic_snv_chr_gnomAD_tbl[i,2]
  end.loc <- somatic_snv_chr_gnomAD_tbl[i,3]
  gr <- GRanges(chr, IRanges(start.loc, end.loc))
  
  ## restrict VCF INFO columns to AC and AN values
  vcfPar <- VariantAnnotation::ScanVcfParam(geno = NA,
                                            fixed = c("ALT", "FILTER"),
                                            info = AFcols,
                                            which = gr)
  
  # Load the desired chunk of VCF file
  vcf <- VariantAnnotation::readVcf(
    Rsamtools::TabixFile("Datasets/gnomAD/V2/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz"),
    param=vcfPar)
  
  # Discard variants not passing all FILTERS
  mask <- VariantAnnotation::filt(vcf) == "PASS"
  vcf <- vcf[mask, ]
  
  # prepare the variants freq data
  ref <- VariantAnnotation::ref(vcf) %>% as.character()
  alt <- VariantAnnotation::alt(vcf) %>% unlist() %>% as.character()
  pos <- MatrixGenerics::rowRanges(vcf) %>% BiocGenerics::start()
  af <- VariantAnnotation::info(vcf) %>% as.data.frame()
  
  # Generate the gnomAD pop. freq. table
  gnomAD_af_tbl <- cbind(data.frame(chr = chr, pos = pos, ref = ref, alt = alt), af)
  
  gnomAD_af_tbl
  
}, mc.cores = parallel::detectCores() - 3)

somatic_snv_popFreq_list <- do.call(rbind, somatic_snv_popFreq_list)

############################3

## Filter out high gnomAD pop freq vars

### identify gnomAD vars that are present in somatic_snv
somatic_snv_gnomAD_pop_freq_index <- which(paste(somatic_snv_popFreq_list$chr, 
                                                 somatic_snv_popFreq_list$pos, 
                                                 somatic_snv_popFreq_list$ref,
                                                 "_",
                                                 somatic_snv_popFreq_list$alt, sep = "") %in% 
                                             paste(somatic_snv$chr, 
                                                   somatic_snv$pos, 
                                                   somatic_snv$ref,
                                                   "_",
                                                   somatic_snv$alt, sep = ""))

### identify somatic_snv_gnomAD_pop_freq_index which have high frequency (FC > 0.01)
somatic_snv_gnomAD_pop_freq_high_freq_index <- (((somatic_snv_popFreq_list[somatic_snv_gnomAD_pop_freq_index,c(5:ncol(somatic_snv_popFreq_list))] > 0.01) %>% 
                                                   rowSums(na.rm = TRUE)) > 0)

### Filter out rare vars of somatic_snv_gnomAD_pop_freq_index
somatic_snv_gnomAD_pop_freq_index <- somatic_snv_gnomAD_pop_freq_index[somatic_snv_gnomAD_pop_freq_high_freq_index]

### identify index of high freq vars in somatic_snv
somatic_snv_gnomAD_index <- match(paste(somatic_snv_popFreq_list$chr, 
                                        somatic_snv_popFreq_list$pos, 
                                        somatic_snv_popFreq_list$ref,
                                        "_",
                                        somatic_snv_popFreq_list$alt, sep = "")[somatic_snv_gnomAD_pop_freq_index],
                                  paste(somatic_snv$chr, 
                                        somatic_snv$pos, 
                                        somatic_snv$ref,
                                        "_",
                                        somatic_snv$alt, sep = ""))

# filter the somatic_snv
somatic_snv_flt <- somatic_snv[-somatic_snv_gnomAD_index,]

# remove the somatic_snv_popFreq_list to free up space
rm(somatic_snv_popFreq_list)

####################################

## Get CADD scores for somatic SNVs

### First, get min and max genomic positions of all chromosomes and create batches of 50,000 bps (genomic positions) for each chromosome
somatic_snv_chr <- 
  unique(somatic_snv$chr)

## Prepare somatic_snv_chr_cadd_tbl table
somatic_snv_chr_cadd_tbl <- somatic_snv_flt[,c(2,3)]
colnames(somatic_snv_chr_cadd_tbl) <- c("chr", "min")
somatic_snv_chr_cadd_tbl$max <- somatic_snv_chr_cadd_tbl$min

somatic_snv_cadd_phred_scores <- scanCADD2df(variants_table = somatic_snv_flt[,-1], 
                                             tabix_snv_file = "Datasets/CADD/v1.6-hg19 (compressed)/whole_genome_SNVs.compressed.tsv.gz",
                                             tabix_indel_file = "Datasets/CADD/v1.6-hg19 (compressed)/InDels.compressed.tsv.gz", 
                                             gr_table = somatic_snv_chr_cadd_tbl)

somatic_snv_cadd_phred_scores <- do.call(rbind, somatic_snv_cadd_phred_scores)

## Add CADD scores to the germline_snv_flt
somatic_snv_flt$cadd_phred[somatic_snv_cadd_phred_scores$index] <- somatic_snv_cadd_phred_scores$score

somatic_snv_flt_withCADD <- somatic_snv_flt
rm(somatic_snv_flt)

# Retrieve missing CADD scores
somatic_snv_cadd_phred_scores_missing_indices <- which(is.na(somatic_snv_flt_withCADD$cadd_phred))

somatic_snv_flt_withCADD_missed <- somatic_snv_flt_withCADD[somatic_snv_cadd_phred_scores_missing_indices,]

somatic_snv_chr_cadd_tbl_missed <- somatic_snv_flt_withCADD_missed[,c(2,3)]
colnames(somatic_snv_chr_cadd_tbl_missed) <- c("chr", "min")
somatic_snv_chr_cadd_tbl_missed$max <- somatic_snv_chr_cadd_tbl_missed$min

somatic_snv_cadd_phred_scores_missed <- scanCADD2df(variants_table = somatic_snv_flt_withCADD_missed[,-1],
                                                    tabix_snv_file = "Datasets/CADD/v1.6-hg19 (compressed)/whole_genome_SNVs.compressed.tsv.gz",
                                                    tabix_indel_file = "Datasets/CADD/v1.6-hg19 (compressed)/InDels.compressed.tsv.gz",
                                                    gr_table = somatic_snv_chr_cadd_tbl_missed)

### Filter out all Mitochondrial variants as there is no CADD scores for such variants
somatic_snv_flt_withCADD <- subset(somatic_snv_flt_withCADD, chr != "MT")

### Assign Ensembl VEP IMPACT type to our variants for imputing the CADD scores of variants with no CADD score

#### Add Variant type to the table (SNV vs InDel)
somatic_snv_flt_withCADD$Type <- "InDel"
somatic_snv_flt_withCADD$Type[which(nchar(somatic_snv_flt_withCADD$ref) == 1 & nchar(somatic_snv_flt_withCADD$alt) == 1)] <- "SNV"

#### Convert multiple consequence variants to single consequence by keeping only the consequence with the highest VEP impact
somatic_snv_flt_withCADD$consequence <-
  sapply(somatic_snv_flt_withCADD$consequence, function(i) {
    tmp_consequence <-
      if(length(grep("&|,", i)) > 0) {
        split_cons <- strsplit(x = i, split = "&|,") %>% unlist()
        split_cons_impact <- vep_impact_tbl$IMPACT_score[match(split_cons, vep_impact_tbl$SO.term)]
        
        split_cons_final <-
          split_cons[which(split_cons_impact == max(split_cons_impact))[1]]
        
        split_cons_final
        
      } else {
        i
      }
    tmp_consequence
  }) %>% unlist()

#### Remove SNV Type variants that belong to frameshift or synonymous (with high CADD phred scores) consequence class:
##### frameshift:their Alt is wrong (there is none in somatic_snv_flt_withCADD)
##### synonymous: these variants strangely have high CADD phred scores (above 15); there is none in somatic_snv_flt_withCADD

somatic_snv_flt_withCADD$VEP_IMPACT <- ""

for(i in c("MODIFIER", "LOW", "MODERATE", "HIGH")) {
  so_terms_index <- which(vep_impact_tbl$IMPACT == i)
  so_terms <- c(vep_impact_tbl$SO.term[so_terms_index], vep_impact_tbl$Display.term[so_terms_index])
  vep_impact_index <- grep(pattern = paste0(so_terms, collapse = "|"), 
                           x = somatic_snv_flt_withCADD$consequence, 
                           ignore.case = TRUE)
  somatic_snv_flt_withCADD$VEP_IMPACT[vep_impact_index] <- i
}

#### Impute the missing CADD scores

##### Calculate the median of CADD scores of each class of VEP Impacts
somatic_VEP_HIGH_median <- median(somatic_snv_flt_withCADD$cadd_phred[somatic_snv_flt_withCADD$VEP_IMPACT == "HIGH"], na.rm = TRUE)
somatic_VEP_MODERATE_median <- median(somatic_snv_flt_withCADD$cadd_phred[somatic_snv_flt_withCADD$VEP_IMPACT == "MODERATE"], na.rm = TRUE)
somatic_VEP_LOW_median <- median(somatic_snv_flt_withCADD$cadd_phred[somatic_snv_flt_withCADD$VEP_IMPACT == "LOW"], na.rm = TRUE)
somatic_VEP_MODIFIER_median <- median(somatic_snv_flt_withCADD$cadd_phred[somatic_snv_flt_withCADD$VEP_IMPACT == "MODIFIER"], na.rm = TRUE)

somatic_snv_flt_withCADD$cadd_phred <- 
  sapply(1:nrow(somatic_snv_flt_withCADD), function(i) {
    imputed_score <- 
      if(!is.na(somatic_snv_flt_withCADD$cadd_phred[i])) {
        somatic_snv_flt_withCADD$cadd_phred[i]
      } else if(somatic_snv_flt_withCADD$VEP_IMPACT[i] == "HIGH") {
        somatic_VEP_HIGH_median
      } else if(somatic_snv_flt_withCADD$VEP_IMPACT[i] == "MODERATE") {
        somatic_VEP_MODERATE_median
      } else if(somatic_snv_flt_withCADD$VEP_IMPACT[i] == "LOW") {
        somatic_VEP_LOW_median
      } else if(somatic_snv_flt_withCADD$VEP_IMPACT[i] == "MODIFIER") {
        somatic_VEP_MODIFIER_median
      }
  })

####################################

## Calibrating the CADD scores (based on the GAM model)

### get the gam probabilities using the logistic transformation
somatic_snv_flt_withCADD_gam_prob <- predict(clinVar_summary_for_training_gam, 
                                             newdata = somatic_snv_flt_withCADD["cadd_phred"], 
                                             type="response")

### Visualize calibration probabilities vs original cadd
somatic_snv_flt_withCADD$calibCADD_GAM_prob <- as.numeric(somatic_snv_flt_withCADD_gam_prob)

somatic_snv_flt_GAM_probvsPhred_plot <-
  ggplot(somatic_snv_flt_withCADD, aes(x = cadd_phred, y = calibCADD_GAM_prob, color = VEP_IMPACT)) +
  geom_point() +
  scale_color_discrete(name = "Ensembl\nVEP IMPACT", type = fishualize::fish(n = 4)) +
  labs(x = "CADD Phred Score", y= "GAM-based Calibrated CADD Score (Pathogenicity Probability)") +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

ggsave(filename = "Results/Somatic SNV/ML-based calibration/somatic_snv_flt_GAM_probvsPhred_plot.pdf", 
       plot = somatic_snv_flt_GAM_probvsPhred_plot, 
       device = "pdf", width = 8, height = 8, units = "in")

# Log transform GAM probabilities
somatic_snv_flt_withCADD$calibCADD_GAM_prob_log <- (-10)*log10(1-somatic_snv_flt_withCADD$calibCADD_GAM_prob)

#################################################

## Retain only PRISM study data
somatic_snv_flt_withCADD <- subset(somatic_snv_flt_withCADD, study == "PRISM")

#### Discard purity NA and amber_qc FAIL, WARN, and NA samples
somatic_snv_pass_index <- which(somatic_snv_flt_withCADD$sample_id %in% subset(lkcgp_sample, !is.na(purity) & amber_qc == "PASS")$sample_id)
somatic_snv_flt_withCADD <- somatic_snv_flt_withCADD[somatic_snv_pass_index,]

#### Retain only samples with purity > 0.1
somatic_snv_high_purity_index <- which(somatic_snv_flt_withCADD$sample_id %in% subset(lkcgp_sample, purity > 0.1)$sample_id)
somatic_snv_flt_withCADD <- somatic_snv_flt_withCADD[somatic_snv_high_purity_index,]

dim(somatic_snv_flt_withCADD)
length(unique(somatic_snv_flt_withCADD$variant_id))
somatic_snv_flt_withCADD$sample_id %>% unique() %>% length()

# # removing ribosomal protein, HLA, and olfactory receptor genes from our data
somatic_snv_flt_withCADD <- somatic_snv_flt_withCADD[-which(stringr::str_detect(string = somatic_snv_flt_withCADD$gene, 
                                                                            pattern = paste(ribo_p_pattern, or_pattern, hla_pattern, sep = "|"))),]

## Summarize the data by summing the calibrated CADD scores of each gene within each sample

#### To evaluate somatic SNV data we need to combine them with germline data as their germline variants have been filtered out
somatic_snv_flt_withCADD_summarized <- rbind(somatic_snv_flt_withCADD[,-c(10,11,17:21)], 
                                             germline_snv_flt_withCADD[,c(colnames(somatic_snv_flt_withCADD[,-c(10,11,17:21)]))])

### Summarize
somatic_snv_flt_withCADD_summarized <- somatic_snv_flt_withCADD_summarized %>% 
  group_by(chr, gene, cancer_type, patient_id) %>% 
  dplyr::summarise(cadd_phred = sum(calibCADD_GAM_prob))

######################

### Re-transform the Phred/log scaled pathogenicities to linear scale probabilities
# somatic_snv_flt_withCADD_summarized$cadd_phred <- 1-10^(somatic_snv_flt_withCADD_summarized$cadd_phred/-10)

######################

# Define the selectGenes function

selectGenes <- function(tbl,                  # table to be inspected (genes on rows and samples on columns)
                        gene_axis = "row",    # are the genes on rows or columns?
                        min.abs=0.1,          # min absolute value to be selected as true variation
                        min.sample=0.05){     # min percentage of samples required for a gene to be considered as true altered gene

  # Check if there is any NA in the table
  if(any(is.na(tbl))) {
    message("There are NA values in the input table. Remove them and re-run the selectGenes function!")
  } else {

  #detect genes that have min.abs in at least samples.length
  keep <- vector(mode = "integer")
  samples.length <- round(min.sample * ifelse(gene_axis == "row", ncol(tbl), nrow(tbl)))
  temp.keep <- apply(abs(tbl), ifelse(gene_axis == "row", 1, 2),
                     function(x, n = samples.length){
                       t = sum(x >= min.abs) >= n
                       t
                     }
  )
  temp.keep <- which(temp.keep)
  keep <- append(keep, temp.keep)

  keep <- unique(keep)
  }
  return(keep)
}

######################

## We should create a list of somatic SNV data with different cancer types and perform subsequent analyses separately for each cancer type
somatic_snv_flt_withCADD_summarized <- lapply(unique(somatic_snv_flt_withCADD_summarized$cancer_type), function(i) {
  somatic_snv_flt_withCADD_summarized[which(somatic_snv_flt_withCADD_summarized$cancer_type == i),]
})

### Set the names of list
somatic_snv_flt_cancer_types <- sapply(somatic_snv_flt_withCADD_summarized, function(i) unique(i$cancer_type))
names(somatic_snv_flt_withCADD_summarized) <- somatic_snv_flt_cancer_types

## Prepare tables in the form of genes on columns and patients on rows for corr analysis (we prepare the table based on patients so that germline and somatic SNV data are combined)
somatic_snv_flt_withCADD_summarized <- lapply(somatic_snv_flt_withCADD_summarized, function(i) {
  tmp.data <- i[, c(2,4,5)]
  tmp.data <- tidyr::pivot_wider(data = tmp.data, 
                                 id_cols = "patient_id", 
                                 names_from = "gene", 
                                 values_from = "cadd_phred",
                                 values_fn = mean) %>% as.data.frame()
  
  rownames(tmp.data) <- tmp.data$patient_id
  tmp.data <- tmp.data[,-1]
  
  # Replace NA with 0
  tmp.data <- apply(tmp.data, 1, function(i) {i[is.na(i)]  <- 0; i}) %>% t() %>% as.data.frame()

  # Report the samples that their number of corresponding genes with data is less than < 1 percent of the total number of genes in the table
  row_index <- apply(tmp.data, 1, function(i) {sum(i > 0)})/ncol(tmp.data) < 0.01
  if(sum(row_index) > 0) {
    tmp_text <- paste("\nThe following samples have mutations in very few number of genes. You may further check their quality.\n",
                      paste0(rownames(tmp.data)[row_index], collapse = "\n"), sep = "")
    cat(tmp_text, sep = "")
  }
  
  # Remove noise genes
  tmp.data <- tmp.data[,colSums(tmp.data > 0.01) != 0]

  tmp.data
})

# Remove tables with less than 3 samples
somatic_snv_flt_withCADD_summarized[which(sapply(somatic_snv_flt_withCADD_summarized, nrow) < 3)] <- NULL

## Get the number of genes in each sample in each dataset
somatic_snv_flt_withCADD_geneCount <- lapply(1:length(somatic_snv_flt_withCADD_summarized), function(i) {
  
  tmp_tbl <- somatic_snv_flt_withCADD_summarized[[i]]
  
  tmp_geneCount <- apply(tmp_tbl, 1, function(j) {sum(j > 0)}) %>% as.data.frame()
  tmp_geneCount <- cbind(cancer_type = somatic_snv_flt_cancer_types[i], 
                         sample_id = rownames(tmp_geneCount), 
                         gene_count = tmp_geneCount[,1]) %>% as.data.frame()
  rownames(tmp_geneCount) <- NULL
  
  tmp_geneCount
})

somatic_snv_flt_withCADD_geneCount <- do.call(rbind, somatic_snv_flt_withCADD_geneCount)

## Perform correlation analysis
somatic_snv_flt_cor <- lapply(somatic_snv_flt_withCADD_summarized, function(i){
  tmp_cor <- fcor(i)
  tmp_cor <- subset(tmp_cor, mr < 20)
  tmp_cor
})

#################

## Create a sample sub-graph of co-del data for visualization

library(influential)
library(igraph)
library(visNetwork)

hgg_somatic_subnet <- influential::graph_from_data_frame(somatic_snv_flt_cor$HGG, directed = F)
hgg_somatic_subnet <- igraph::subgraph(graph = hgg_somatic_subnet, vids = V(hgg_somatic_subnet)[1:2000])
hgg_somatic_subnets_genes <- igraph::cliques(graph = hgg_somatic_subnet, min = 15, max = 20)
hgg_somatic_subnet_final <- subgraph(hgg_somatic_subnet, vids = V(hgg_somatic_subnet)[hgg_somatic_subnets_genes[[63]]])

hgg_somatic_subnet_final_vis <-
  visIgraph(igraph = hgg_somatic_subnet_final) %>% 
  visOptions(
    # nodesIdSelection = TRUE, 
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '30px arial black',
    size = 25,
    # color = "orange",
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "circle" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = "lightgrey")

visExport(graph = hgg_somatic_subnet_final_vis, 
          name = "hgg_somatic_subnet_final_vis.png", 
          label = "Subset of HGG co-deleteriousness network")

######################## For pan cancer analysis ######################## 

## Prepare the table in the form of genes on columns and samples on rows
pan_somatic_snv_tbl <- rbind(somatic_snv_flt_withCADD[,-c(10,11,17:21)], 
                             germline_snv_flt_withCADD[,c(colnames(somatic_snv_flt_withCADD[,-c(10,11,17:21)]))])

pan_somatic_snv_n_cancer <- pan_somatic_snv_tbl$cancer_type %>% unique() %>% length()

pan_somatic_snv_tbl <- pan_somatic_snv_tbl %>% 
  group_by(chr, gene, patient_id) %>% 
  dplyr::summarise(cadd_phred = sum(calibCADD_GAM_prob))

pan_somatic_snv_tbl <- tidyr::pivot_wider(data = pan_somatic_snv_tbl, 
                                          id_cols = "patient_id", 
                                          names_from = "gene", 
                                          values_from = "cadd_phred",
                                          values_fn = mean) %>% as.data.frame()

rownames(pan_somatic_snv_tbl) <- pan_somatic_snv_tbl$patient_id
pan_somatic_snv_tbl <- pan_somatic_snv_tbl[,-1]

# Replace NA with 0
pan_somatic_snv_tbl <- apply(pan_somatic_snv_tbl, 1, function(i) {i[is.na(i)]  <- 0; i}) %>% t() %>% as.data.frame()

# Report the samples that their number of corresponding genes with data is less than < 0.5 percent of the total number of genes in the table
pan_somatic_snv_tbl_row_index <- apply(pan_somatic_snv_tbl, 1, function(i) {sum(i > 0)})/ncol(pan_somatic_snv_tbl) < 0.01

if(sum(pan_somatic_snv_tbl_row_index) > 0) {
  tmp_text <- paste("\nThe following samples have mutations in very few number of genes. You may further check their quality.\n",
                    paste0(rownames(pan_somatic_snv_tbl)[pan_somatic_snv_tbl_row_index], collapse = "\n"), sep = "")
  cat(tmp_text, sep = "")
}

# Remove noise genes
pan_somatic_snv_tbl <- pan_somatic_snv_tbl[,colSums(pan_somatic_snv_tbl > 0.01) != 0]

# Perform corr analysis
pan_somatic_snv_cor <-  fcor(data = pan_somatic_snv_tbl, 
                             method = "spearman", 
                             mutualRank = TRUE, pvalue = FALSE, flat = TRUE)

# Filter top 20 (MR < 20) correlations
pan_somatic_snv_cor <- subset(pan_somatic_snv_cor, mr < 20)

#=============================================================================
#
#    Code chunk 5: Get and prepare the transcriptome data
#
#=============================================================================

###### Redundant ######
somatic_rnaexp <- tbl(db, "lkcgp_curated_sample_somatic_rnaexp") %>% as.data.frame()
###### Redundant ######

rnaseq <- tbl(db, "rnaseq") %>% as.data.frame()
rnaseq <- rnaseq[,c(1:3, 7, 8, 11, 12, 13, 14)]

## Keep only PRISM samples (present in lkcgp_sample)
rnaseq <- rnaseq[which(rnaseq$rnaseq_id %in% lkcgp_sample$rnaseq_id),]

## Add cancer type to the table
rnaseq$cancer_type <- lkcgp_sample$cancer_type[match(rnaseq$rnaseq_id, lkcgp_sample$rnaseq_id)]

#### Discard purity NA and amber_qc FAIL, WARN, and NA
rnaseq_pass_index <- which(rnaseq$rnaseq_id %in% subset(lkcgp_sample, !is.na(purity) & amber_qc == "PASS")$rnaseq_id)
rnaseq <- rnaseq[rnaseq_pass_index,]

#### Retain only rnaseq data with purity > 0.1
rnaseq_high_purity_index <- which(rnaseq$rnaseq_id %in% subset(lkcgp_sample, purity > 0.1)$rnaseq_id)
rnaseq <- rnaseq[rnaseq_high_purity_index,]

#### # removing ribosomal protein, HLA, and olfactory receptor genes from our data
rnaseq <- rnaseq[-which(stringr::str_detect(string = rnaseq$gene, pattern = paste(ribo_p_pattern, or_pattern, hla_pattern, sep = "|"))),]

#### We should create a list of transcriptome data with different cancer types and perform subsequent analyses separately for each cancer type

rnaseq_expr_list <- rnaseq
rnaseq_expr_list <- lapply(unique(rnaseq_expr_list$cancer_type), function(i) {
  rnaseq_expr_list[which(rnaseq_expr_list$cancer_type == i),]
})

### Set the names of list
rnaseq_expr_cancer_types <- sapply(rnaseq_expr_list, function(i) unique(i$cancer_type))
names(rnaseq_expr_list) <- rnaseq_expr_cancer_types

#######################################################

# Check the accuracy of co-expressions based on external rna-seq data
nbl_target_2018 <- read.delim("~/Downloads/data_mrna_seq_rpkm.txt", sep = "\t")

## Remove duplicate genes
nbl_target_2018 <- nbl_target_2018[match(unique(nbl_target_2018$Hugo_Symbol), nbl_target_2018$Hugo_Symbol),]
rownames(nbl_target_2018) <- nbl_target_2018$Hugo_Symbol 
nbl_target_2018 <- nbl_target_2018[,-c(1,2)]

nbl_target_2018 <- as.data.frame(t(nbl_target_2018))

# Replace NA with 0
nbl_target_2018 <- apply(nbl_target_2018, 1, function(i) {i[is.na(i)]  <- 0; i}) %>% t() %>% as.data.frame()

nbl_target_2018_keep_index <- selectGenes(tbl = nbl_target_2018, gene_axis = "column",
                                          min.abs = 1, # At least RPKM of 1
                                          min.sample = 0.1) # In at least 10 percent of the samples
nbl_target_2018 <- nbl_target_2018[,nbl_target_2018_keep_index]

nbl_target_2018_cor <- fcor(data = nbl_target_2018, 
                            method = "spearman", 
                            mutualRank = TRUE, pvalue = FALSE, flat = TRUE)

nbl_target_2018_cor <- subset(nbl_target_2018_cor, mr < 10) # Get highly reliable correlations

# Split our data based on purity
rnaseq_expr_nbl <- rnaseq_expr_list$NBL
rnaseq_expr_nbl <- list(
  rnaseq_expr_nbl_0.1_0.3 = rnaseq_expr_nbl[which(rnaseq_expr_nbl$rnaseq_id %in% subset(lkcgp_sample, purity > 0.1 & purity <= 0.3)$rnaseq_id),],
  rnaseq_expr_nbl_0.3_0.5 = rnaseq_expr_nbl[which(rnaseq_expr_nbl$rnaseq_id %in% subset(lkcgp_sample, purity > 0.3 & purity <= 0.5)$rnaseq_id),],
  rnaseq_expr_nbl_0.5_0.8 = rnaseq_expr_nbl[which(rnaseq_expr_nbl$rnaseq_id %in% subset(lkcgp_sample, purity > 0.5 & purity <= 0.8)$rnaseq_id),],
  rnaseq_expr_nbl_0.8_1.0 = rnaseq_expr_nbl[which(rnaseq_expr_nbl$rnaseq_id %in% subset(lkcgp_sample, purity > 0.8 & purity <= 1)$rnaseq_id),]
  
)

## Prepare tables in the form of genes on columns and samples on rows for corr analysis
rnaseq_expr_nbl <- lapply(rnaseq_expr_nbl, function(i) {
  tmp.data <- i[,c(1,3,4)]
  tmp.data <- tidyr::pivot_wider(data = tmp.data, 
                                 id_cols = "rnaseq_id", 
                                 names_from = "gene", 
                                 values_from = "tpm") %>% as.data.frame()
  
  rownames(tmp.data) <- tmp.data$rnaseq_id
  tmp.data <- tmp.data[,-1]
  
  # Replace NA with 0
  tmp.data <- apply(tmp.data, 1, function(i) {i[is.na(i)]  <- 0; i}) %>% t() %>% as.data.frame()
  
  # Report the samples that their number of corresponding genes with data is less than < 1 percent of the total number of genes in the table
  row_index <- apply(tmp.data, 1, function(i) {sum(i > 0)})/ncol(tmp.data) < 0.01
  if(sum(row_index) > 0) {
    tmp_text <- paste("\nThe following samples have mutations in very few number of genes. You may further check their quality.\n",
                      paste0(rownames(tmp.data)[row_index], collapse = "\n"), sep = "")
    cat(tmp_text, sep = "")
  }
  
  # Remove the genes that might be noise (their expression is very low in very low number of samples)
  rnaseq_expr_tbl_keep_index <- selectGenes(tbl = tmp.data, gene_axis = "column",
                                            min.abs = 1, # At least TPM of 1
                                            min.sample = 0.1) # In at least 10 percent of the samples
  tmp.data <- tmp.data[,rnaseq_expr_tbl_keep_index]
  
  tmp.data
})

## Perform correlation analysis
rnaseq_expr_nbl_cor <- lapply(rnaseq_expr_nbl, function(i){
  
  tmp_cor <- fcor(data = i, 
                  method = "spearman", 
                  mutualRank = TRUE, pvalue = FALSE, flat = TRUE)
  
  tmp_cor <- subset(tmp_cor, mr < 20)
  tmp_cor
})

# Check intersection
nbl_target_intersect <-
  which(paste(nbl_target_2018_cor$row, nbl_target_2018_cor$column, sep = "_") %in% 
          c(paste(rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.8_1.0$row, rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.8_1.0$column, sep = "_"), 
            paste(rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.8_1.0$column, rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.8_1.0$row, sep = "_"))) %>% length()

(nbl_target_intersect/nrow(nbl_target_2018_cor)) *100 
# 0.1_0.3 ~ 0 % overlap with cBioPortal
# 0.3_0.5 ~ 0.36 % overlap with cBioPortal
# 0.5_0.8 ~ 3 % overlap with cBioPortal
# 0.8_1.0 ~ 3.5 % overlap with cBioPortal

## Check the overlap of correlations of each section of our data with others
nbl_cor_vec_list <-  list(
  cor_0.1_0.3 = c(paste(rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.1_0.3$row, rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.1_0.3$column, sep = "_"), 
                  paste(rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.1_0.3$column, rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.1_0.3$row, sep = "_")),
  cor_0.3_0.5 = c(paste(rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.3_0.5$row, rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.3_0.5$column, sep = "_"), 
                  paste(rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.3_0.5$column, rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.3_0.5$row, sep = "_")),
  cor_0.5_0.8 = c(paste(rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.5_0.8$row, rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.5_0.8$column, sep = "_"), 
                  paste(rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.5_0.8$column, rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.5_0.8$row, sep = "_")),
  cor_0.8_1.0 = c(paste(rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.8_1.0$row, rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.8_1.0$column, sep = "_"), 
                  paste(rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.8_1.0$column, rnaseq_expr_nbl_cor$rnaseq_expr_nbl_0.8_1.0$row, sep = "_"))
    )

nbl_cor_vec_list_combinations <- combn(x = c(1:length(nbl_cor_vec_list)), m = 2) %>% as.data.frame() %>% as.list()

nbl_cor_intersect <- lapply(nbl_cor_vec_list_combinations, function(i) {
  intersect(x = unlist(nbl_cor_vec_list[i[[1]]]), y = unlist(nbl_cor_vec_list[i[[2]]])) %>% length()
})

#######################################################

## Prepare tables in the form of genes on columns and samples on rows for corr analysis
rnaseq_expr_list <- lapply(rnaseq_expr_list, function(i) {
  tmp.data <- i[,c(1,3,4)]
  tmp.data <- tidyr::pivot_wider(data = tmp.data, 
                                 id_cols = "rnaseq_id", 
                                 names_from = "gene", 
                                 values_from = "tpm") %>% as.data.frame()
  
  rownames(tmp.data) <- tmp.data$rnaseq_id
  tmp.data <- tmp.data[,-1]
  
  # Replace NA with 0
  tmp.data <- apply(tmp.data, 1, function(i) {i[is.na(i)]  <- 0; i}) %>% t() %>% as.data.frame()
  
  # Report the samples that their number of corresponding genes with data is less than < 1 percent of the total number of genes in the table
  row_index <- apply(tmp.data, 1, function(i) {sum(i > 0)})/ncol(tmp.data) < 0.01
  if(sum(row_index) > 0) {
    tmp_text <- paste("\nThe following samples have mutations in very few number of genes. You may further check their quality.\n",
                      paste0(rownames(tmp.data)[row_index], collapse = "\n"), sep = "")
    cat(tmp_text, sep = "")
  }
  
  # Remove the genes that might be noise (their expression is very low in very low number of samples)
  rnaseq_expr_tbl_keep_index <- selectGenes(tbl = tmp.data, gene_axis = "column",
                                            min.abs = 1, # At least TPM of 1
                                            min.sample = 0.1) # In at least 10 percent of the samples
  tmp.data <- tmp.data[,rnaseq_expr_tbl_keep_index]
  
  tmp.data
})

# Remove tables with less than 3 samples
rnaseq_expr_list[which(sapply(rnaseq_expr_list, nrow) < 3)] <- NULL

## Perform correlation analysis
rnaseq_expr_cor <- lapply(rnaseq_expr_list, function(i){
  
  tmp_cor <- fcor(data = i, 
                  method = "spearman", 
                  mutualRank = TRUE, pvalue = FALSE, flat = TRUE)
  
  tmp_cor <- subset(tmp_cor, mr < 20)
  tmp_cor
})

#################

## Create a sample sub-graph of co-exp data for visualization

library(influential)
library(igraph)
library(visNetwork)

hgg_somatic_tx_subnet <- influential::graph_from_data_frame(rnaseq_expr_cor$HGG, directed = F)
hgg_somatic_tx_subnet <- igraph::subgraph(graph = hgg_somatic_tx_subnet, vids = V(hgg_somatic_tx_subnet)[5000:15000])
hgg_somatic_tx_subnets_genes <- igraph::cliques(graph = hgg_somatic_tx_subnet, min = 10, max = 20)
hgg_somatic_tx_subnet_final <- subgraph(hgg_somatic_tx_subnet, vids = V(hgg_somatic_tx_subnet)[hgg_somatic_tx_subnets_genes[[1]]])

hgg_somatic_tx_subnet_final_vis <-
  visIgraph(igraph = hgg_somatic_tx_subnet_final) %>% 
  visOptions(
    # nodesIdSelection = TRUE, 
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '25px arial black',
    # size = 25,
    color = "orange",
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "box" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = "lightgrey")

visExport(graph = hgg_somatic_tx_subnet_final_vis, 
          name = "hgg_somatic_tx_subnet_final_vis.png", 
          label = "Subset of HGG co-deleteriousness network")

######################## For pan cancer analysis ######################## 

## Prepare RNA-seq table in the form of genes on columns and samples on rows
pan_somatic_rnaseq_tbl <- rnaseq

pan_somatic_rnaseq_n_cancer <- pan_somatic_rnaseq_tbl$cancer_type %>% unique() %>% length()

pan_somatic_rnaseq_tbl <- pan_somatic_rnaseq_tbl[,c(1,3,4)]

pan_somatic_rnaseq_tbl <- tidyr::pivot_wider(data = pan_somatic_rnaseq_tbl, 
                                      id_cols = "rnaseq_id", 
                                      names_from = "gene", 
                                      values_from = "tpm") %>% as.data.frame()

rownames(pan_somatic_rnaseq_tbl) <- pan_somatic_rnaseq_tbl$rnaseq_id
pan_somatic_rnaseq_tbl <- pan_somatic_rnaseq_tbl[,-1]

# Replace NA with 0
pan_somatic_rnaseq_tbl <- apply(pan_somatic_rnaseq_tbl, 1, function(i) {i[is.na(i)]  <- 0; i}) %>% t() %>% as.data.frame()

# Remove the genes that might be noise (their rnaseq is very low in very high number of samples)
pan_somatic_rnaseq_tbl_keep_index <- selectGenes(tbl = pan_somatic_rnaseq_tbl, gene_axis = "column",
                                                 min.abs = 1, min.sample = 0.1)
pan_somatic_rnaseq_tbl <- pan_somatic_rnaseq_tbl[,pan_somatic_rnaseq_tbl_keep_index]

# Perform corr analysis
pan_somatic_rnaseq_cor <-  fcor(data = pan_somatic_rnaseq_tbl, 
                                     method = "spearman", 
                                     mutualRank = TRUE, pvalue = FALSE, flat = TRUE)

# Filter top 20 (MR < 20) correlations
pan_somatic_rnaseq_cor <- subset(pan_somatic_rnaseq_cor, mr < 20)

rm(pan_somatic_rnaseq_tbl)

#######################################################################

#=============================================================================
#
#    Code chunk 6: Import and pre-process normal RNA-seq data (GTEx data)
#
#=============================================================================

# Read in the expression TPM file
gtex_data_lines <- vroom::vroom_lines("~/Downloads/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz")
buf_size <- sum(nchar(gtex_data_lines[1:109]))
buf_size
Sys.setenv(VROOM_CONNECTION_SIZE = as.integer(buf_size))

gtex_data <- vroom::vroom("~/Downloads/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz", delim = "\t", skip = 2)

# Set the gene Ensembl accessions as rownames and remove it
rownames(gtex_data) <- gtex_data$Name
gtex_data$Name <- NULL

# Summarize duplicate gene names by getting their average
gtex_data_dup_index <- duplicated(gtex_data$Description) %>% which()
gtex_data_dup_index <- which(gtex_data$Description %in% gtex_data$Description[gtex_data_dup_index])

gtex_data_dup <- gtex_data[gtex_data_dup_index,]
gtex_data <- gtex_data[-gtex_data_dup_index,]

gtex_data_dup_geneNames <- unique(gtex_data_dup$Description)
gtex_data_dup <- lapply(unique(gtex_data_dup$Description), FUN = function(i) {
  tmp_data <- gtex_data_dup[which(gtex_data_dup$Description ==  i),]
  tmp_data <- colMeans(tmp_data[,-1])
  tmp_data
})

names(gtex_data_dup) <- gtex_data_dup_geneNames
gtex_data_dup <- do.call(rbind, gtex_data_dup)
gtex_data_dup <- cbind(Description = gtex_data_dup_geneNames, gtex_data_dup)

gtex_data <- rbind(gtex_data, gtex_data_dup)
gtex_data <- as.data.frame(gtex_data)

rownames(gtex_data) <- gtex_data$Description
gtex_data_geneNames <- gtex_data$Description
gtex_data$Description <- NULL

# Import the annotation file
gtex_data_annot <- vroom::vroom("~/Downloads/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt", delim = "\t")

# Add tissue type to the expression table
gtex_data_tissue_index <- match(colnames(gtex_data), gtex_data_annot$SAMPID)
gtex_data_tissues <- gtex_data_annot$SMTS[gtex_data_tissue_index]

# create a list of germline gtex data with different tissue types
gtex_data_list <- lapply(unique(gtex_data_tissues), function(i) {
  gtex_data[,which(gtex_data_tissues == i)]
})

# Remove original data to free up memory
rm(gtex_data)

names(gtex_data_list) <- unique(gtex_data_tissues)

# Add another table (combination of Adipose Tissue, Muscle, and Blood Vessel for Sarcoma)
gtex_data_list <- append(gtex_data_list, list(do.call(cbind, gtex_data_list[c("Adipose Tissue", "Muscle", "Blood Vessel")])))
names(gtex_data_list)[31] <- "Adipose Tissue; Muscle; Blood Vessel"

# Convert the data to numeric
gtex_data_list <- lapply(gtex_data_list, function(i) {
  geneNames_tmp <- rownames(i)
  table_tmp <- apply(i, 2, as.numeric) %>% as.data.frame()
  rownames(table_tmp) <- geneNames_tmp
  table_tmp
})

# Remove noise genes
gtex_data_list <- lapply(gtex_data_list, function(i) {
  keep_tmp <- selectGenes(tbl = i[-nrow(i),], gene_axis = "row",
                          min.abs = 1, # At least TPM of 1
                          min.sample = 0.1) # In at least 10 percent of the samples
  table_tmp <- i[keep_tmp,]
  table_tmp
})

gtex_data_list_names <- names(gtex_data_list)

#### # removing ribosomal protein, HLA, and olfactory receptor genes from our data
gtex_data_list <- lapply(1:length(gtex_data_list), function(i) {
  gtex_data_list[[i]] <<- gtex_data_list[[i]][-which(stringr::str_detect(string = rownames(gtex_data_list[[i]]), 
                                                                         pattern = paste(ribo_p_pattern, or_pattern, hla_pattern, sep = "|"))),]
})

names(gtex_data_list) <- gtex_data_list_names
  
## Save the dataset
gtex_data_list_con <-  gzfile("Datasets/GTEx/gtex_V8_preprocessed.rds.gz")
saveRDS(gtex_data_list, gtex_data_list_con)
close(gtex_data_list_con)
rm(gtex_data_list_con)

# Perform correlation analysis
gtex_data_cor_list <- lapply(gtex_data_list, function(i) {
  tmp_data <- fcor(data = t(i))
  tmp_data <- subset(tmp_data, mr < 20)
  tmp_data
})

## Save the correlations
gtex_data_cor_list_con <-  gzfile("Datasets/GTEx/gtex_V8_correlations.rds.gz")
saveRDS(gtex_data_cor_list, gtex_data_cor_list_con)
close(gtex_data_cor_list_con)
rm(gtex_data_cor_list_con)

rm(gtex_data_list)

######################## For pan cancer analysis ######################## 

## Prepare RNA-seq table in the form of genes on columns and samples on rows

pan_normal_rnaseq_tbl <- read_rds("Datasets/GTEx/gtex_V8_preprocessed.rds.gz")

## Remove the last manually created table (Adipose Tissue;Muscle; Blood Vessel)
pan_normal_rnaseq_tbl$`Adipose Tissue;Muscle; Blood Vessel` <- NULL

pan_normal_rnaseq_n_tissue <- pan_normal_rnaseq_tbl %>% length()

pan_normal_rnaseq_tbl <- lapply(pan_normal_rnaseq_tbl, function(i) as.data.frame(t(i)))

pan_normal_rnaseq_tbl <- data.table::rbindlist(pan_normal_rnaseq_tbl, fill = TRUE)

# Replace NA with 0
pan_normal_rnaseq_tbl <- apply(pan_normal_rnaseq_tbl, 1, function(i) {i[is.na(i)]  <- 0; i}) %>% t() %>% as.data.frame()

# Remove the genes that might be noise (their rnaseq is very low in very high number of samples)
pan_normal_rnaseq_tbl_keep_index <- selectGenes(tbl = pan_normal_rnaseq_tbl, gene_axis = "column",
                                                min.abs = 1, min.sample = 0.1)
pan_normal_rnaseq_tbl <- pan_normal_rnaseq_tbl[,pan_normal_rnaseq_tbl_keep_index]

# Save the table
pan_normal_rnaseq_tbl_con <-  gzfile("Datasets/GTEx/pan_tissue_gtex_V8_preprocessed.rds.gz")
saveRDS(pan_normal_rnaseq_tbl, pan_normal_rnaseq_tbl_con)
close(pan_normal_rnaseq_tbl_con)
rm(pan_normal_rnaseq_tbl_con)

# Perform corr analysis
pan_normal_rnaseq_cor <-  influential::fcor(data = pan_normal_rnaseq_tbl, 
                               method = "spearman", 
                               mutualRank = TRUE, pvalue = FALSE, flat = TRUE)

# Filter top 20 (MR < 20) correlations
pan_normal_rnaseq_cor <- subset(pan_normal_rnaseq_cor, mr < 20)

# Save the table
pan_normal_rnaseq_cor_con <- gzfile("Datasets/GTEx/pan_tissue_gtex_V8_correlations.rds.gz")
saveRDS(pan_normal_rnaseq_cor, pan_normal_rnaseq_cor_con)
close(pan_normal_rnaseq_cor_con)
rm(pan_normal_rnaseq_cor_con)

rm(pan_normal_rnaseq_tbl)

#=============================================================================
#
#    Code chunk 7: Reconstruction and analysis of the Risk Net
#
#=============================================================================

# Create the multi-layer Risk network

## First remove "Other" from the lists of correlations as the cancer type is not known and we will perform a pan-cancer analysis anyway
germline_snv_flt_cor$Other <- NULL

risk_net <- germline_snv_flt_cor

# Prepare the tables for merging with other data
risk_net <- lapply(risk_net, function(i) {
  tmp_tbl <- i[,c(1,2)]
  colnames(tmp_tbl) <- c("from", "to")
  tmp_tbl$type <- "co-deletriousness"
  tmp_tbl
})

# Add the co-expression data to the tables
## We complement the 
### CNS other, HGG (high-grade glioma), NBL (Neuroblastoma), DMG (Diffuse Midline Glioma), CNS embryonal, Glioma other, MB (medulloblastoma), and EPD (ependymoma) with Brain coexpression data
### Rhabdoid with Kidney coexpression data
### Leukaemia and Lymphoma with Blood coexpression data
### Sarcoma with Adipose Tissue; Muscle; Blood Vessel coexpression data

for(i in 1:length(risk_net)) {
  
  cancer_type = names(risk_net[i]) 
  
  genes <- risk_net[[i]][,c(1,2)] %>% 
    unlist() %>% 
    unname() %>% 
    unique()
  
  if(cancer_type %in% c("CNS other", "HGG", "NBL", "DMG", "CNS embryonal", "Glioma other", "MB", "EPD")) {
    
    first_level_coex_index <- c(which(gtex_data_cor_list$Brain$row %in% genes), 
                                which(gtex_data_cor_list$Brain$column %in% genes)) %>% unique()
    
    first_level_coex_genes <- c(genes, 
                                gtex_data_cor_list$Brain$row[first_level_coex_index],
                                gtex_data_cor_list$Brain$column[first_level_coex_index]) %>% unique()
    
    first_two_levels_coex_index <- c(which(gtex_data_cor_list$Brain$row %in% first_level_coex_genes), 
                                     which(gtex_data_cor_list$Brain$column %in% first_level_coex_genes)) %>% unique()
    
    tmp_tbl <- gtex_data_cor_list$Brain[first_two_levels_coex_index,c(1,2)]
    colnames(tmp_tbl) <- c("from", "to")
    tmp_tbl$type <- "co-expression"
    
    risk_net[[i]] <- rbind(risk_net[[i]], tmp_tbl)
    
  } else if(cancer_type %in% c("Leukaemia", "Lymphoma")) {
    
    first_level_coex_index <- c(which(gtex_data_cor_list$Blood$row %in% genes), 
                                which(gtex_data_cor_list$Blood$column %in% genes)) %>% unique()
    
    first_level_coex_genes <- c(genes, 
                                gtex_data_cor_list$Blood$row[first_level_coex_index],
                                gtex_data_cor_list$Blood$column[first_level_coex_index]) %>% unique()
    
    first_two_levels_coex_index <- c(which(gtex_data_cor_list$Blood$row %in% first_level_coex_genes), 
                                     which(gtex_data_cor_list$Blood$column %in% first_level_coex_genes)) %>% unique()
    
    tmp_tbl <- gtex_data_cor_list$Blood[first_two_levels_coex_index,c(1,2)]
    colnames(tmp_tbl) <- c("from", "to")
    tmp_tbl$type <- "co-expression"
    
    risk_net[[i]] <- rbind(risk_net[[i]], tmp_tbl)
    
  } else if(cancer_type %in% c("Rhabdoid")) {
    
    first_level_coex_index <- c(which(gtex_data_cor_list$Kidney$row %in% genes), 
                                which(gtex_data_cor_list$Kidney$column %in% genes)) %>% unique()
    
    first_level_coex_genes <- c(genes, 
                                gtex_data_cor_list$Kidney$row[first_level_coex_index],
                                gtex_data_cor_list$Kidney$column[first_level_coex_index]) %>% unique()
    
    first_two_levels_coex_index <- c(which(gtex_data_cor_list$Kidney$row %in% first_level_coex_genes), 
                                     which(gtex_data_cor_list$Kidney$column %in% first_level_coex_genes)) %>% unique()
    
    tmp_tbl <- gtex_data_cor_list$Kidney[first_two_levels_coex_index,c(1,2)]
    colnames(tmp_tbl) <- c("from", "to")
    tmp_tbl$type <- "co-expression"
    
    risk_net[[i]] <- rbind(risk_net[[i]], tmp_tbl)
    
  } else if(cancer_type %in% c("Sarcoma")) {
    
    first_level_coex_index <- c(which(gtex_data_cor_list$`Adipose Tissue; Muscle; Blood Vessel`$row %in% genes), 
                                which(gtex_data_cor_list$`Adipose Tissue; Muscle; Blood Vessel`$column %in% genes)) %>% unique()
    
    first_level_coex_genes <- c(genes, 
                                gtex_data_cor_list$`Adipose Tissue; Muscle; Blood Vessel`$row[first_level_coex_index],
                                gtex_data_cor_list$`Adipose Tissue; Muscle; Blood Vessel`$column[first_level_coex_index]) %>% unique()
    
    first_two_levels_coex_index <- c(which(gtex_data_cor_list$`Adipose Tissue; Muscle; Blood Vessel`$row %in% first_level_coex_genes), 
                                     which(gtex_data_cor_list$`Adipose Tissue; Muscle; Blood Vessel`$column %in% first_level_coex_genes)) %>% unique()
    
    tmp_tbl <- gtex_data_cor_list$`Adipose Tissue; Muscle; Blood Vessel`[first_two_levels_coex_index,c(1,2)]
    colnames(tmp_tbl) <- c("from", "to")
    tmp_tbl$type <- "co-expression"
    
    risk_net[[i]] <- rbind(risk_net[[i]], tmp_tbl)
    
  }
}

# Add PPI data to the Risk net
library(STRINGdb)

# Map gene names to stringDB ids for PPI analysis

## get the gene names in each risk net
risk_net_genes <- lapply(risk_net, function(i) {
  data.frame(gene_symbol = i[,c(1,2)] %>% 
               unlist() %>% 
               unname() %>% 
               unique())
})

## getSTRINGdb for human
string_db <- STRINGdb$new(species=9606, version = "11.5", 
                          score_threshold=900, # very high confidence (>0.9)
                          input_directory=""
                          )

risk_net_genes <- lapply(risk_net_genes, function(i) {
  
  string_db$map(my_data_frame = i, 
                my_data_frame_id_col_names = "gene_symbol", 
                takeFirst = TRUE, removeUnmappedRows = TRUE)

})

# Get the PPIs
risk_net_ppi <- lapply(risk_net_genes, function(i) {
  
  string_db$get_interactions(string_ids = i$STRING_id)
  
})

# Convert stringdb ids to gene symbols
risk_net_ppi <- lapply(1:length(risk_net_ppi), function(i) {
  tmp_row_index <- match(risk_net_ppi[[i]]$from, risk_net_genes[[i]]$STRING_id)
  risk_net_ppi[[i]]$from <- risk_net_genes[[i]]$gene_symbol[tmp_row_index]
  
  tmp_column_index <- match(risk_net_ppi[[i]]$to, risk_net_genes[[i]]$STRING_id)
  risk_net_ppi[[i]]$to <- risk_net_genes[[i]]$gene_symbol[tmp_column_index]
  
  tmp_tbl <- risk_net_ppi[[i]][,c(1,2)]
  tmp_tbl$type <- "PPI"
  
  ## Remove duplicates
  tmp_tbl <- tmp_tbl[-which(duplicated(paste(tmp_tbl$from,
                                             tmp_tbl$to, 
                                             sep = "_"))),]
  
  return(tmp_tbl)
})

names(risk_net_ppi) <- names(risk_net_genes)

#################

## Create a sample sub-graph of ppi data for visualization

library(influential)
library(igraph)
library(visNetwork)

hgg_ppi_subnet <- influential::graph_from_data_frame(risk_net_ppi$HGG, directed = F)
hgg_ppi_subnet <- igraph::subgraph(graph = hgg_ppi_subnet, vids = V(hgg_ppi_subnet)[1:500])
hgg_ppi_subnets_genes <- igraph::cliques(graph = hgg_ppi_subnet, min = 15, max = 20)
hgg_ppi_subnet_final <- subgraph(hgg_ppi_subnet, vids = V(hgg_ppi_subnet)[hgg_ppi_subnets_genes[[12]]])

hgg_ppi_subnet_final_vis <-
  visIgraph(igraph = hgg_ppi_subnet_final) %>% 
  visOptions(
    # nodesIdSelection = TRUE, 
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '25px arial black',
    # size = 25,
    color = "lightblue",
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "ellipse" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = "lightgrey")

visExport(graph = hgg_ppi_subnet_final_vis, 
          name = "hgg_ppi_subnet_final_vis.png", 
          label = "Subset of HGG co-deleteriousness network")

####################

# Combine PPI with the original risk net
risk_net <- lapply(1:length(risk_net), function(i) {
  rbind(risk_net[[i]], risk_net_ppi[[i]])
})

names(risk_net) <- names(risk_net_genes)

# Reconstruct the networks
risk_net <- lapply(risk_net, function(i) {
  igraph::graph_from_data_frame(i, directed=FALSE)
})

## Check the network layers
table(igraph::E(risk_net$`CNS other`)$type)

###############################

# Calculate the IVI values of genes within each net
library(influential)
risk_net_ivi <- lapply(risk_net, ivi, verbose = T)

###############################

# Calculate the primitive gene risk scores
risk_net_gene_tbl <- lapply(risk_net_ivi, function(i) {
  data.frame(gene = names(i), ivi = i)
})

## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer

### First remove other
germline_snv_flt_withCADD_summarized$Other <- NULL

germline_snv_mean_deletriousness <- lapply(germline_snv_flt_withCADD_summarized, function(i) {
  colMeans(i)
})

## Generate the gene table for each net and calculate the primitive scores
risk_net_gene_tbl <- lapply(1:length(risk_net_gene_tbl), function(i) {
  tmp_tbl <- risk_net_gene_tbl[[i]]
  tmp_tbl$mean_deletriousness <- 0
  match_index <- which(names(germline_snv_mean_deletriousness[[i]]) %in% rownames(tmp_tbl))
  tmp_tbl[names(germline_snv_mean_deletriousness[[i]])[match_index], "mean_deletriousness"] <- 
    germline_snv_mean_deletriousness[[i]][match_index]
  
  tmp_tbl$primitive_risk_score <- (tmp_tbl$ivi)*(tmp_tbl$mean_deletriousness) # node weights
  
  return(tmp_tbl)
  
})

names(risk_net_gene_tbl) <- names(germline_snv_mean_deletriousness)

# add gene scores
risk_net <- lapply(1:length(risk_net), function(i) {
  set_vertex_attr(graph = risk_net[[i]], name = "score", value = risk_net_gene_tbl[[i]]$primitive_risk_score)
})

names(risk_net) <- names(germline_snv_mean_deletriousness)

# Convert the un-weighted to a node-weighted network
risk_net <- lapply(1:length(risk_net), function(i) {
  set_vertex_attr(graph = risk_net[[i]], name = "weight", value = (risk_net_gene_tbl[[i]]$ivi)^(risk_net_gene_tbl[[i]]$mean_deletriousness))
})

names(risk_net) <- names(germline_snv_mean_deletriousness)

## Check the network node weights
summary(igraph::V(risk_net$Sarcoma)$weight)

# Calculate the mediator scores by multiplying the scores to the mean neighborhood weight of nodes
risk_net_gene_tbl <- lapply(1:length(risk_net), function(i) {
  tmp_tbl <- risk_net_gene_tbl[[i]]
  tmp_tbl$mean_neighborhood_weight <- sapply(1:nrow(tmp_tbl), function(j) {
    igraph::neighbors(graph = risk_net[[i]], 
                      v = V(risk_net[[i]])[j])$score %>% mean()
  })
  
  tmp_tbl$risk_mediator_score <- tmp_tbl$primitive_risk_score * tmp_tbl$mean_neighborhood_weight
  
  return(tmp_tbl)
})

names(risk_net_gene_tbl) <- names(germline_snv_mean_deletriousness)

########################

## List of HGG drivers suggested by https://pecan.stjude.cloud/
hgg_driver_genes <- c("TP53", "H3F3A", "TTN", "ATRX", "PDGFRA", "C1ORF112", "PIK3CA", "ACVR1", 
                      "ZBTB3", "LRRC7", "HIST1H3B", "NF1", "RNF213", "BCOR", "NEB", "LRCH3", 
                      "MUC16", "BRD1", "OBSCN", "CREBBP")

  ## Top 20 Risk genes of HGG
risk_net_gene_tbl$HGG %>% 
  slice_max(order_by = primitive_risk_score, n = 20) %>% 
  mutate(Rank = rank(-(primitive_risk_score))) %>% 
  select(Rank, gene) %>% 
  rename(Gene = gene) %>% 
  View()

## Top 20 mediators of HGG
risk_net_gene_tbl$HGG %>% 
  slice_max(order_by = risk_mediator_score, n = 20) %>% 
  mutate(Rank = rank(-(risk_mediator_score))) %>% 
  select(Rank, gene) %>% 
  rename(Gene = gene) %>% 
  View()

########################

## Create a dataframe of top 20 risk genes of all cancers
risk_net_top20_genes <- data.frame(Rank = c(1:20))

for(i in 1:length(risk_net_gene_tbl)) {
  tmp_data = risk_net_gene_tbl[[i]] %>% 
    slice_max(order_by = primitive_risk_score, n = 20) %>% 
    select(gene) %>% 
    rename(cancer_type = gene)
  
  risk_net_top20_genes <- cbind(risk_net_top20_genes, tmp_data)
}

colnames(risk_net_top20_genes)[2:ncol(risk_net_top20_genes)] <- names(risk_net_gene_tbl)

########################

## Sarcoma drivers suggested by https://pecan.stjude.cloud/
sarcoma_driver_genes <- unique(c("TP53", "VNN3", "RB1", "ATRX", "DMD", "TTN", "CSMD1", "MUC16", "FER", "DGKB", "LAMA5", "AUTS2", "DLG2", "NAALADL2", "SYT16", 
                                 "MIRLET7BHG", "EYS", "LRRC7", "OTOGL", "C1ORF112", "EWSR1", "FLI1", "STAG2", "ERG", "AP1B1", "TTC28", "TP53", "EEF2", "ZNRF3", 
                                 "RUNX1", "EMID1", "ERF", "CSMD1", "SLC35F4", "RYR3", "MUC16", "DEPDC5", "ETS1", "LRRC7", "CACNG2", "ASCC1", "FOXO1", "SYT16", 
                                 "TP53", "LRP1B", "PAX3", "PSMG3-AS1", "FGFR4", "GRM1", "PIK3CA", "SYNE1", "BCOR", "SCRIB", "MUC6", "NOTCH1", "LRRC7", "TUG1", 
                                 "LOC100499484-C9ORF174", "SEC24A", "ACSM2A", "PERM1", "KRTAP4-8", "UBXN11", "ITGAM", "ZNF891", "NF1", "ADCY8", "PTPN11", "C17ORF97", 
                                 "BCOR", "TTBK1", "PCCA", "MN1", "OTOP1", "DMD", "KDM5C", "GPR161", "POLR1A", "ADGRV1", "TENM2"))

## Top 20 Risk genes of Sarcoma
risk_net_gene_tbl$Sarcoma %>% 
  slice_max(order_by = primitive_risk_score, n = 20) %>% 
  mutate(Rank = rank(-(primitive_risk_score))) %>% 
  select(Rank, gene) %>% 
  rename(Gene = gene) %>% 
  View()

## Top 20 mediators of Sarcoma
risk_net_gene_tbl$Sarcoma %>% 
  slice_max(order_by = risk_mediator_score, n = 20) %>% 
  mutate(Rank = rank(-(risk_mediator_score))) %>% 
  select(Rank, gene) %>% 
  rename(Gene = gene) %>% 
  View()

############################################# 

# Detect communities/modules/clusters
risk_net_modules <- lapply(risk_net, function(i) {
  set.seed(3847)
  igraph::cluster_leiden(
    graph = i,
    objective_function = "CPM",
    weights = NULL,
    resolution_parameter = 0.5,
    beta = 0.05,
    initial_membership = NULL,
    n_iterations = 100,
    vertex_weights = igraph::V(i)$weight
  )
})

sizes(risk_net_modules$`CNS other`) %>% as.vector() %>% summary()
length(risk_net_modules$`CNS other`)

## Visualize the entire network and color by modules (NOT Recommended; It will take a while and will be too crowded)
plot(
  x = risk_net_modules$`CNS other`, 
  y = risk_net$`CNS other`,
  col = membership(risk_net_modules$`CNS other`),
  vertex.label=NA,
  mark.groups = communities(risk_net_modules$`CNS other`),
  edge.color = c("black", "red")[crossing(risk_net_modules$`CNS other`, risk_net$`CNS other`) + 1]
)

###########################################

## Inspect and filter modules
risk_net_modules_flt <- lapply(1:length(risk_net_modules), function(i) {
  
  ### Create a table of modules and their genes
  tmp_tbl <- data.frame(
    module = membership(risk_net_modules[[i]]) %>% as.integer(),
    gene = risk_net_modules[[i]]$names
  )
  
  ### Add if the gene is mutated or not (could be a risk gene or not)
  tmp_tbl$mutated <- 0
  tmp_tbl$mutated[which(risk_net_gene_tbl[[i]]$mean_deletriousness > 0)] <- 1
  
  ### Separate modules
  tmp_tbl <- lapply(unique(tmp_tbl$module), function(j){
    subset(tmp_tbl, module == j)
  })
  
  names(tmp_tbl) <- sapply(tmp_tbl, function(m) {
    unique(m$module)
  })
  
  ### Remove modules that have less than 2 Risk genes or their size is less than 4
  tmp_tbl <- lapply(tmp_tbl, function(k) {
    if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
      k
    } else {
      NULL
    }
  })
  
  ### Remove NULL modules
  tmp_tbl[which(sapply(tmp_tbl, is.null))] <- NULL
  
  ### Calculate the module scores
  tmp_tbl <- lapply(tmp_tbl, function(l) {
    tmp_tbl4score <- cbind(l, mean_score = 0)
    tmp_tbl4score$mean_score <- mean(risk_net_gene_tbl[[i]]$primitive_risk_score[which(risk_net_gene_tbl[[i]]$gene %in% l$gene)])
    tmp_tbl4score
  })
  
  tmp_tbl
})

names(risk_net_modules_flt) <- names(risk_net_modules)

## Sort the mean scores of all modules of each cancer type
risk_net_modules_flt_scores <-
lapply(risk_net_modules_flt, function(i) {
  tmp_tbl <-
    sapply(i, function(j) {unique(j$mean_score)})
  tmp_tbl <- rev(sort(tmp_tbl))
  tmp_tbl
})

## Visualize the modules
plot.igraph(subgraph(graph = risk_net$`CNS other`, 
                     vids = which(as_ids(V(risk_net$`CNS other`)) %in%  risk_net_modules_flt$`CNS other`$`37`$gene)), 
            vertex.size =  3)

#####################

# Create a dataframe of 1st-ranked modules of all cancer types

risk_net_1st_modules <- data.frame(Number = vector("integer"))

for(i in 1:length(risk_net_modules_flt)) {
  tmp_module_genes = risk_net_modules_flt[[i]][names(risk_net_modules_flt_scores[[i]])[1]][[1]]$gene
  risk_net_1st_modules <- merge(risk_net_1st_modules, 
                                data.frame(Number = c(1:length(tmp_module_genes)),
                                           "cancer_type" = tmp_module_genes), 
                                by = "Number", all = TRUE)
}

risk_net_1st_modules$Number <- NULL

colnames(risk_net_1st_modules) <- names(risk_net_modules_flt)

#####################

## Check Rhabdoid top genes
rhabdoid_risk_genes <- c("ARID1A", "ARID1B",
                         "SMARCA2", "SMARCA4",
                         "SMARCC1", "SMARCC2",
                         "SMARCD1", "SMARCD2", "SMARCD3",
                         "ACTL6A", "ACTL6B",
                         "SMARCB1")

risk_net_gene_tbl$Rhabdoid$Rank <- rank(-1*(risk_net_gene_tbl$Rhabdoid$primitive_risk_score))
risk_net_gene_tbl$Rhabdoid %>% filter(gene %in% rhabdoid_risk_genes) %>% view()

############

# Check the module rhabdoid_risk_genes are involved in

Rhabdoid_risk_including_module_names <- vector(mode = "character")

for(i in 1:length(risk_net_modules_flt$Rhabdoid)) {
  module_tmp <-
    grep(paste0(paste(paste("^", rhabdoid_risk_genes, sep = ""), "$", sep = ""), collapse = "|"), risk_net_modules_flt$Rhabdoid[[i]]$gene)
  if(length(module_tmp >= 1)) { 
    print(names(risk_net_modules_flt$Rhabdoid[i]))
    print(risk_net_modules_flt$Rhabdoid[[i]])
    cat("\n###########################################\n")
    Rhabdoid_risk_including_module_names <<- append(Rhabdoid_risk_including_module_names,
                                                   names(risk_net_modules_flt$Rhabdoid[i]))
  }
}


Rhabdoid_risk_including_module_ranks <- which(names(risk_net_modules_flt$Rhabdoid) %in% Rhabdoid_risk_including_module_names)
rhabdoid_risk_genes[which(rhabdoid_risk_genes %in% risk_net_modules_flt$Rhabdoid$`130`$gene)] %>% cat()


#####################

## Create module 1 of HGG for visualization

library(influential)
library(igraph)
library(visNetwork)

risk_hgg_module1 <- subgraph(graph = risk_net$HGG, 
                             vids = which(as_ids(V(risk_net$HGG)) %in%  risk_net_modules_flt$HGG$`590`$gene))

risk_hgg_module1_vis <-
  visIgraph(igraph = risk_hgg_module1) %>% 
  visOptions(
    # nodesIdSelection = TRUE, 
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '25px arial black',
    size = 25,
    # color = c("lightgrey", "orange"),
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "box" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = c("lightgrey"))

visExport(graph = risk_hgg_module1_vis, 
          name = "risk_hgg_module1_vis", 
          label = "")

## It's also possible to create a subset of two modules and color the nodes based on mark.groups arg
### https://stackoverflow.com/questions/26913419/plot-communities-with-igraph

#####################

## Create module 1 of Sarcoma for visualization

risk_net_modules_flt_scores$Sarcoma %>% head()

library(influential)
library(igraph)
library(visNetwork)

risk_Sarcoma_module1 <- subgraph(graph = risk_net$Sarcoma, 
                                 vids = which(as_ids(V(risk_net$Sarcoma)) %in%  risk_net_modules_flt$Sarcoma$`2375`$gene))

risk_Sarcoma_module1_vis <-
  visIgraph(igraph = risk_Sarcoma_module1) %>% 
  visOptions(
    # nodesIdSelection = TRUE, 
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '25px arial black',
    size = 25,
    # color = c("lightgrey", "orange"),
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "box" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = c("lightgrey"))

visExport(graph = risk_Sarcoma_module1_vis, 
          name = "risk_Sarcoma_module1_vis", 
          label = "")

#####################

## Create Risk module 1 of Lymphoma for visualization

risk_net_modules_flt_scores$Lymphoma %>% head()

library(influential)
library(igraph)
library(visNetwork)

risk_Lymphoma_module1 <- subgraph(graph = risk_net$Lymphoma, 
                                  vids = which(as_ids(V(risk_net$Lymphoma)) %in%  risk_net_modules_flt$Lymphoma$`974`$gene))

risk_Lymphoma_module1_node_colors <- rep("red", length(V(risk_Lymphoma_module1)))
risk_Lymphoma_module1_node_colors[c(1,4,5)] <- "lightgrey"
risk_Lymphoma_module1_node_colors[2] <- "lightgreen"
risk_Lymphoma_module1_node_colors[6] <- "yellow"
risk_Lymphoma_module1 <- set_vertex_attr(graph = risk_Lymphoma_module1, name = "color", value = risk_Lymphoma_module1_node_colors)

risk_Lymphoma_module1_vis <-
  visIgraph(igraph = risk_Lymphoma_module1) %>% 
  visOptions(
    # nodesIdSelection = TRUE, 
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '25px arial black',
    size = 25,
    # color = c("lightgrey", "orange"),
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "box" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = c("lightgrey"))

visExport(graph = risk_Lymphoma_module1_vis, 
          name = "risk_Lymphoma_module1_vis", 
          label = "")

#####################

## Create driver module 1 of Lymphoma for visualization

driver_marker_net_modules_flt_scores$Lymphoma %>% head()

library(influential)
library(igraph)
library(visNetwork)

driver_marker_Lymphoma_module1 <- subgraph(graph = driver_marker_net$Lymphoma, 
                                           vids = which(as_ids(V(driver_marker_net$Lymphoma)) %in%  driver_marker_net_modules_flt$Lymphoma$`1035`$gene))

driver_marker_Lymphoma_module1_node_colors <- rep("yellow", length(V(driver_marker_Lymphoma_module1)))
driver_marker_Lymphoma_module1_node_colors[1] <- "violet"
driver_marker_Lymphoma_module1_node_colors[5] <- "red"
driver_marker_Lymphoma_module1 <- set_vertex_attr(graph = driver_marker_Lymphoma_module1, name = "color", value = driver_marker_Lymphoma_module1_node_colors)

driver_marker_Lymphoma_module1_vis <-
  visIgraph(igraph = driver_marker_Lymphoma_module1) %>% 
  visOptions(
    # nodesIdSelection = TRUE, 
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '25px arial black',
    size = 25,
    # color = c("lightgrey", "orange"),
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "box" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = c("lightgrey"))

visExport(graph = driver_marker_Lymphoma_module1_vis, 
          name = "driver_marker_Lymphoma_module1_vis", 
          label = "")

#############################################

# Create a function for inspecting the edge type between specific nodes

edge_type <- function(g, v) {
  library(igraph)
  node_indices <- grep(paste0(paste(paste("^", v, sep = ""), "$", sep = ""), collapse = "|"), as_ids(V(g)))
  if(length(v) == length(node_indices)) {
    possible_combn <- combn(x = node_indices, m = 2)
    for(i in 1:ncol(possible_combn)) {
      edge_type_tmp <- E(g)[possible_combn[1,i] %--% possible_combn[2,i]]$type
      cat(as_ids(V(g)[possible_combn[1,i]]), 
          "-",
          as_ids(V(g)[possible_combn[2,i]]),
          ":", edge_type_tmp, "\n")
    }
  } else {
    message("Some input node names are not present in the network")
  }
}

#############################################

edge_type(g = risk_net$HGG, v = c("VENTX", "CLEC5A","SFRP4", "CCND1"))

############

# Check the module a specific gene is involved in
for (i in 1:length(risk_net_modules_flt$HGG)) {
  module_tmp <-
    grep("^TP53$", risk_net_modules_flt$HGG[[i]]$gene)
  if(length(module_tmp >= 1)) print(i)
}

##########################

# Check the module TP53 is involved in
for (i in 1:length(risk_net_modules_flt$HGG)) {
  module_tmp <-
    grep("^TP53$", risk_net_modules_flt$HGG[[i]]$gene)
  if(length(module_tmp >= 1)) { 
    print(names(risk_net_modules_flt$HGG[i]))
    print(risk_net_modules_flt$HGG[[i]])
    cat("\n###########################################\n")
  }
}

##########################

#### Separate the 5 first ranked modules for the visualization
hgg_first_ranked_risk_module_net <- subgraph(graph = risk_net$HGG, 
                                             vids = which(as_ids(V(risk_net$HGG)) %in% risk_net_modules_flt$HGG$`617`$gene))
plot(hgg_first_ranked_risk_module_net)

hgg_first_ranked_risk_module_netTbl <- igraph::as_data_frame(hgg_first_ranked_risk_module_net)
hgg_first_ranked_risk_module_netTbl$mutated <- ""
hgg_first_ranked_risk_module_netTbl <- rbind(hgg_first_ranked_risk_module_netTbl, 
                                             data.frame(from = risk_net_modules_flt$HGG$`617`$gene, 
                                                        to = "", type = "", 
                                                        mutated = risk_net_modules_flt$HGG$`617`$mutated))

write.csv(x = hgg_first_ranked_risk_module_netTbl, 
          file = "Results/Networks/HGG top 5 risk modules/hgg_first_ranked_risk_module.csv", 
          quote = F, row.names = F)

####################################################
####################################################

# Create a heatmap of HLA and CTS Risk genes

## Get the desired genes from the 1st rank risk modules of all cancer types
hla_cts_genes <- vector(mode = "character")

for(i in 1:length(risk_net_modules_flt_scores)) {
  tmp_genes = risk_net_modules_flt[[i]][names(risk_net_modules_flt_scores[[i]][1])][[1]] %>% 
    select(gene) %>% 
    unlist() %>% 
    grep(pattern = "HLA|CTS", value = TRUE)
  
  hla_cts_genes <- append(hla_cts_genes, tmp_genes)
}

hla_cts_genes <- hla_cts_genes %>% unique() %>% unname()

############

## Prepare the Del data

### Somatic  + Germline
somatic_plus_germline_snv_CADD_for_eval <- rbind(somatic_snv_flt_withCADD[,-c(10,11,17:21)], 
                                                 germline_snv_flt_withCADD[,c(colnames(somatic_snv_flt_withCADD[,-c(10,11,17:21)]))])

#### Summarize
somatic_plus_germline_snv_CADD_for_eval <- somatic_plus_germline_snv_CADD_for_eval %>% 
  group_by(gene, cancer_type, patient_id) %>% 
  dplyr::summarise(cadd_phred = sum(calibCADD_GAM_prob))

somatic_plus_germline_snv_CADD_for_eval <- somatic_plus_germline_snv_CADD_for_eval %>% 
  filter(gene %in% hla_cts_genes)

#### Order the table by cancer type (for visualization)
somatic_plus_germline_snv_CADD_for_eval <- somatic_plus_germline_snv_CADD_for_eval[order(somatic_plus_germline_snv_CADD_for_eval$cancer_type),]

#### Transform the table for visualization

somatic_plus_germline_snv_CADD_for_eval <- 
somatic_plus_germline_snv_CADD_for_eval %>% 
pivot_wider(
  names_from = gene,
  values_from = cadd_phred,
  values_fill = NA
) %>% 
  data.frame()

somatic_plus_germline_snv_CADD_for_eval_cancers <- somatic_plus_germline_snv_CADD_for_eval$cancer_type
rownames(somatic_plus_germline_snv_CADD_for_eval) <- somatic_plus_germline_snv_CADD_for_eval$patient_id
somatic_plus_germline_snv_CADD_for_eval <- somatic_plus_germline_snv_CADD_for_eval[,-c(1,2)]
somatic_plus_germline_snv_CADD_for_eval <- t(somatic_plus_germline_snv_CADD_for_eval)

#### Correct the order of genes
somatic_plus_germline_snv_CADD_for_eval <- somatic_plus_germline_snv_CADD_for_eval[order(rownames(somatic_plus_germline_snv_CADD_for_eval)),]

########

### Somatic only
somatic_only_snv_CADD_for_eval <- somatic_snv_flt_withCADD[,-c(10,11,17:21)]

#### Summarize
somatic_only_snv_CADD_for_eval <- somatic_only_snv_CADD_for_eval %>% 
  group_by(gene, cancer_type, patient_id) %>% 
  dplyr::summarise(cadd_phred = sum(calibCADD_GAM_prob))

somatic_only_snv_CADD_for_eval <- somatic_only_snv_CADD_for_eval %>% 
  filter(gene %in% hla_cts_genes)

#### Order the table by cancer type (for visualization)
somatic_only_snv_CADD_for_eval <- somatic_only_snv_CADD_for_eval[order(somatic_only_snv_CADD_for_eval$cancer_type),]

#### Transform the table for visualization

somatic_only_snv_CADD_for_eval <- 
  somatic_only_snv_CADD_for_eval %>% 
  pivot_wider(
    names_from = gene,
    values_from = cadd_phred,
    values_fill = NA
  ) %>% 
  data.frame()

somatic_only_snv_CADD_for_eval_cancers <- somatic_only_snv_CADD_for_eval$cancer_type
rownames(somatic_only_snv_CADD_for_eval) <- somatic_only_snv_CADD_for_eval$patient_id
somatic_only_snv_CADD_for_eval <- somatic_only_snv_CADD_for_eval[,-c(1,2)]
somatic_only_snv_CADD_for_eval <- t(somatic_only_snv_CADD_for_eval)

#### Correct the order of genes
somatic_only_snv_CADD_for_eval <- somatic_only_snv_CADD_for_eval[order(rownames(somatic_only_snv_CADD_for_eval)),]

########

### Germline only
germline_only_snv_CADD_for_eval <- germline_snv_flt_withCADD[,c(colnames(somatic_snv_flt_withCADD[,-c(10,11,17:21)]))]

#### Summarize
germline_only_snv_CADD_for_eval <- germline_only_snv_CADD_for_eval %>% 
  group_by(gene, cancer_type, patient_id) %>% 
  dplyr::summarise(cadd_phred = sum(calibCADD_GAM_prob))

germline_only_snv_CADD_for_eval <- germline_only_snv_CADD_for_eval %>% 
  filter(gene %in% hla_cts_genes)

#### Order the table by cancer type (for visualization)
germline_only_snv_CADD_for_eval <- germline_only_snv_CADD_for_eval[order(germline_only_snv_CADD_for_eval$cancer_type),]

#### Tranforme the table for visualization

germline_only_snv_CADD_for_eval <- 
  germline_only_snv_CADD_for_eval %>% 
  pivot_wider(
    names_from = gene,
    values_from = cadd_phred,
    values_fill = NA
  ) %>% 
  data.frame()

germline_only_snv_CADD_for_eval_cancers <- germline_only_snv_CADD_for_eval$cancer_type
rownames(germline_only_snv_CADD_for_eval) <- germline_only_snv_CADD_for_eval$patient_id
germline_only_snv_CADD_for_eval <- germline_only_snv_CADD_for_eval[,-c(1,2)]
germline_only_snv_CADD_for_eval <- t(germline_only_snv_CADD_for_eval)

#### Correct the order of genes
germline_only_snv_CADD_for_eval <- germline_only_snv_CADD_for_eval[order(rownames(germline_only_snv_CADD_for_eval)),]

################################

## Visualize

library(pheatmap)

### Somatic  + Germline

#### Prepare the annotation data and colors

somatic_plus_germline_snv_CADD_for_eval_metaData <- data.frame(Cancer_Type = somatic_plus_germline_snv_CADD_for_eval_cancers)
rownames(somatic_plus_germline_snv_CADD_for_eval_metaData) <- colnames(somatic_plus_germline_snv_CADD_for_eval)

#### Draw the heatmap

somatic_plus_germline_snv_CADD_for_eval_heatmap <-
pheatmap(somatic_plus_germline_snv_CADD_for_eval, 
         cluster_cols = F, cluster_rows = F, show_colnames = F, 
         main = "The Mutation Burden of HLA and CTS genes in Combined Somatic and Germline Data", 
         legend_labels = c("2", "4", "6", "8", "10", "12", "Deleteriousness\n\n"),
         annotation_col = somatic_plus_germline_snv_CADD_for_eval_metaData)

somatic_plus_germline_snv_CADD_for_eval_heatmap

ggsave(somatic_plus_germline_snv_CADD_for_eval_heatmap, device = "pdf", width = 11, height = 6, 
       filename = "Results/Mut burden of HLA and CTS/somatic_plus_germline.pdf")

################

### Somatic only

#### Prepare the annotation data and colors

somatic_only_snv_CADD_for_eval_metaData <- data.frame(Cancer_Type = somatic_only_snv_CADD_for_eval_cancers)
rownames(somatic_only_snv_CADD_for_eval_metaData) <- colnames(somatic_only_snv_CADD_for_eval)

#### Draw the heatmap

somatic_only_snv_CADD_for_eval_heatmap <-
  pheatmap(somatic_only_snv_CADD_for_eval, 
           cluster_cols = F, cluster_rows = F, show_colnames = F, 
           main = "The Mutation Burden of HLA and CTS genes in Somatic Data", 
           legend_labels = c("2", "4", "6", "8", "10", "12", "Deleteriousness\n\n"),
           annotation_col = somatic_only_snv_CADD_for_eval_metaData)

somatic_only_snv_CADD_for_eval_heatmap

ggsave(somatic_only_snv_CADD_for_eval_heatmap, device = "pdf", width = 11, height = 6, 
       filename = "Results/Mut burden of HLA and CTS/somatic_only.pdf")

################

### Germline only

#### Prepare the annotation data and colors

germline_only_snv_CADD_for_eval_metaData <- data.frame(Cancer_Type = germline_only_snv_CADD_for_eval_cancers)
rownames(germline_only_snv_CADD_for_eval_metaData) <- colnames(germline_only_snv_CADD_for_eval)

#### Draw the heatmap

germline_only_snv_CADD_for_eval_heatmap <-
  pheatmap(germline_only_snv_CADD_for_eval, 
           cluster_cols = F, cluster_rows = F, show_colnames = F, 
           main = "The Mutation Burden of HLA and CTS genes in Combined Somatic and Germline Data", 
           legend_labels = c("2", "4", "6", "8", "10", "12", "Deleteriousness\n\n"),
           annotation_col = germline_only_snv_CADD_for_eval_metaData)

germline_only_snv_CADD_for_eval_heatmap

ggsave(germline_only_snv_CADD_for_eval_heatmap, device = "pdf", width = 11, height = 6, 
       filename = "Results/Mut burden of HLA and CTS/germline_only.pdf")

##########################################################################################
##########################################################################################

# Evaluation of the first ranked Risk module of Sarcoma

risk_net_modules_flt_scores$Sarcoma %>% head()

## First ranked Risk Module corresponding mutated samples
sarcoma_frist_risk_module_mutated <- risk_net_modules_flt$Sarcoma$`2375`
sarcoma_frist_risk_module_mutated_genes = sarcoma_frist_risk_module_mutated$gene[sarcoma_frist_risk_module_mutated$mutated == 1]

## subset the SNV table based on mutated genes
frist_risk_module_snv <- germline_snv_flt_withCADD_summarized$Sarcoma[,sarcoma_frist_risk_module_mutated_genes]

### Add cancer subtype
frist_risk_module_snv$cancer_subtype <- lkcgp_sample$diagnosis[match(rownames(frist_risk_module_snv), lkcgp_sample$matched_normal_id)]
frist_risk_module_snv$patient_id <- lkcgp_sample$patient_id[match(rownames(frist_risk_module_snv), lkcgp_sample$matched_normal_id)]

## samples including mutated genes
sarcoma_frist_risk_module_mutated_samples <- (apply(frist_risk_module_snv[,c(1:4)], 
                                                    1, 
                                                    function(i) {
                                                      return(sum(i > 0))
                                                    }
) > 0) %>% 
  which() %>% 
  names()

sarcoma_frist_risk_module_mutated_samples <- frist_risk_module_snv[sarcoma_frist_risk_module_mutated_samples,]
sarcoma_frist_risk_module_mutated_samples_tbl <- sarcoma_frist_risk_module_mutated_samples
sarcoma_frist_risk_module_mutated_samples_tbl <- sarcoma_frist_risk_module_mutated_samples_tbl[order(sarcoma_frist_risk_module_mutated_samples_tbl$cancer_subtype),]
sarcoma_frist_risk_module_mutated_samples_tbl <- sarcoma_frist_risk_module_mutated_samples_tbl[,c(1,3,2,4:6)]

rownames(sarcoma_frist_risk_module_mutated_samples_tbl) <- sarcoma_frist_risk_module_mutated_samples_tbl$patient_id

################

# Oncoplot for mutations of First ranked Risk Module genes in their corresponding samples

sarcoma_frist_risk_module_mutated_samples_cancer_subtype <- data.frame(cancer_subtype = sarcoma_frist_risk_module_mutated_samples_tbl$cancer_subtype)
rownames(sarcoma_frist_risk_module_mutated_samples_cancer_subtype) <- rownames(sarcoma_frist_risk_module_mutated_samples_tbl)

sarcoma_frist_risk_module_mutated_samples_cancer_subtype_colors <-
  list(cancer_subtype = c(EWS = fish(n = 6)[1], MPNST = fish(n = 6)[2], OST = fish(n = 6)[3], 
                          `RMS FN` = fish(n = 6)[4], `RMS FP` = fish(n = 6)[5], `Sarcoma other` = fish(n = 6)[6]))

sarcoma_frist_risk_module_mutated_samples_oncoplot <- pheatmap(t(sarcoma_frist_risk_module_mutated_samples_tbl[,-c(5,6)]),
                                                               scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
                                                               color = colorRampPalette(colors = c("lightgrey", "darkgreen"))(100),
                                                               annotation_col = sarcoma_frist_risk_module_mutated_samples_cancer_subtype, 
                                                               annotation_colors = sarcoma_frist_risk_module_mutated_samples_cancer_subtype_colors,
                                                               fontsize_row = 6, fontsize_col = 6,
                                                               main = "Oncoplot of Sarcoma First ranked Risk Module Genes")

ggsave(filename = "Results/Evaluations/Sarcoma first ranked Risk module/sarcoma_frist_risk_module_mutated_samples_oncoplot.pdf", 
       plot = sarcoma_frist_risk_module_mutated_samples_oncoplot,
       device = "pdf", 
       width = 18, height = 3.8, units = "cm")

#################

# Visualization of mutated genes and their genomic information
library(ggpubr)

sarcoma_frist_risk_module_mutated_samples_genomics_tbl <-
germline_snv_flt_withCADD %>% filter(gene %in% sarcoma_frist_risk_module_mutated_genes & sample_id %in% rownames(sarcoma_frist_risk_module_mutated_samples))

sarcoma_frist_risk_module_mutated_samples_genomics_tbl <- sarcoma_frist_risk_module_mutated_samples_genomics_tbl[order(sarcoma_frist_risk_module_mutated_samples_genomics_tbl$gene),]

## add the genome strand of genes
sarcoma_frist_risk_module_mutated_samples_genomics_tbl$strand <- 1
sarcoma_frist_risk_module_mutated_samples_genomics_tbl$strand[grep("ESR2", sarcoma_frist_risk_module_mutated_samples_genomics_tbl$gene)] <- -1
sarcoma_frist_risk_module_mutated_samples_genomics_tbl$strand[grep("NCOR1", sarcoma_frist_risk_module_mutated_samples_genomics_tbl$gene)] <- -1

sarcoma_frist_risk_module_mutated_samples_genomics_tbl_plot <- ggtexttable(sarcoma_frist_risk_module_mutated_samples_genomics_tbl[,c("patient_id", "gene", "chr", "pos", "ref", "alt", "consequence")],
                                                                           rows = NULL, theme = ttheme("lBlack"),
                                                                           cols = c("Patient ID", "Gene", "Chr", "Pos", "Ref", "Alt", "Consequence")) %>% 
  table_cell_bg(row = grep(sarcoma_frist_risk_module_mutated_genes[1], 
                           sarcoma_frist_risk_module_mutated_samples_genomics_tbl$gene) + 1, column = 2, fill = "darkolivegreen1") %>% 
  table_cell_bg(row = grep(sarcoma_frist_risk_module_mutated_genes[2], 
                           sarcoma_frist_risk_module_mutated_samples_genomics_tbl$gene) + 1, column = 2, fill = "firebrick1") %>% 
  table_cell_bg(row = grep(sarcoma_frist_risk_module_mutated_genes[3], 
                           sarcoma_frist_risk_module_mutated_samples_genomics_tbl$gene) + 1, column = 2, fill = "cyan") %>% 
  table_cell_bg(row = grep(sarcoma_frist_risk_module_mutated_genes[4], 
                           sarcoma_frist_risk_module_mutated_samples_genomics_tbl$gene) + 1, column = 2, fill = "deeppink")
  # table_cell_font(row = 2, column = 1, face = "bold")

sarcoma_frist_risk_module_mutated_samples_genomics_tbl_plot <- sarcoma_frist_risk_module_mutated_samples_genomics_tbl_plot %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(nrow(sarcoma_frist_risk_module_mutated_samples_genomics_tbl_plot) + 1), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(sarcoma_frist_risk_module_mutated_samples_genomics_tbl_plot), column.side = "left", from.row = 2, linetype = 2) %>%
  tab_add_title(text = "Sarcoma first ranked module genomic information", face = "bold", padding = unit(1.5, "line"))

sarcoma_frist_risk_module_mutated_samples_genomics_tbl_plot

# Save the visualization
ggsave(filename = "Results/Evaluations/Sarcoma first ranked Risk module/sarcoma_frist_risk_module_mutated_samples_genomics_tbl_plot.pdf", 
       plot = sarcoma_frist_risk_module_mutated_samples_genomics_tbl_plot,
       device = "pdf", 
       width = 22, height = 26, units = "cm")

##########################################

# Evaluation of the effect of mutations

## for SIFT (https://sift.bii.a-star.edu.sg/www/Extended_SIFT_chr_coords_submit.html)

paste(sarcoma_frist_risk_module_mutated_samples_genomics_tbl$chr, 
      sarcoma_frist_risk_module_mutated_samples_genomics_tbl$pos,
      sarcoma_frist_risk_module_mutated_samples_genomics_tbl$strand,
      paste(sarcoma_frist_risk_module_mutated_samples_genomics_tbl$ref,
            paste(sarcoma_frist_risk_module_mutated_samples_genomics_tbl$alt, "\n", sep = ""), sep = "/"), 
      sep = ",") %>% cat()

#################################################
#################################################

# Oncoplot for the entire dataset

entire_cohort_sarcoma_frist_risk_module_snv <- lapply(germline_snv_flt_withCADD_summarized, function(i) {
  tmp_genes = sarcoma_frist_risk_module_mutated_genes[which(sarcoma_frist_risk_module_mutated_genes %in% colnames(i))]
  if(length(tmp_genes) > 0) {
    i[,tmp_genes, drop = FALSE]
  } else {
    NULL
  }
})

### Add cancer subtype
entire_cohort_sarcoma_frist_risk_module_snv <- lapply(1:length(entire_cohort_sarcoma_frist_risk_module_snv), function(i) {
  entire_cohort_sarcoma_frist_risk_module_snv[[i]]$cancer_type <- names(entire_cohort_sarcoma_frist_risk_module_snv)[i]
  entire_cohort_sarcoma_frist_risk_module_snv[[i]]
})

entire_cohort_sarcoma_frist_risk_module_snv <- do.call(plyr::rbind.fill, entire_cohort_sarcoma_frist_risk_module_snv)
entire_cohort_sarcoma_frist_risk_module_snv <- entire_cohort_sarcoma_frist_risk_module_snv[,c(sarcoma_frist_risk_module_mutated_genes, "cancer_type")]

entire_cohort_sarcoma_frist_risk_module_snv[is.na(entire_cohort_sarcoma_frist_risk_module_snv)] <- 0

## samples including mutated genes
# entire_cohort_sarcoma_frist_risk_module_snv_mutated_samples <- (apply(entire_cohort_sarcoma_frist_risk_module_snv[,c(1:4)], 
#                                                     1, 
#                                                     function(i) {
#                                                       return(sum(i > 0))
#                                                     }
# ) > 0) %>% 
#   which()
# 
# entire_cohort_sarcoma_frist_risk_module_snv <- entire_cohort_sarcoma_frist_risk_module_snv[entire_cohort_sarcoma_frist_risk_module_snv_mutated_samples,]

##############

## for the purpose of using pheatmpa the input dataframe should have rownames
rownames(entire_cohort_sarcoma_frist_risk_module_snv) <- paste("row_", c(1:nrow(entire_cohort_sarcoma_frist_risk_module_snv)), sep = "")

sarcoma_frist_risk_module_mutated_samples_cancer_type <- data.frame(cancer_type = entire_cohort_sarcoma_frist_risk_module_snv$cancer_type)
rownames(sarcoma_frist_risk_module_mutated_samples_cancer_type) <- rownames(entire_cohort_sarcoma_frist_risk_module_snv)

sarcoma_frist_risk_module_mutated_samples_cancer_type_colors <-
  list(cancer_type = set_names(fish(n = 12, option = "Acanthurus_sohal"), nm = names(germline_snv_flt_withCADD_summarized)))

entire_cohort_sarcoma_frist_risk_module_mutated_samples_oncoplot <- pheatmap(t(entire_cohort_sarcoma_frist_risk_module_snv[,-c(5)]),
                                                               scale = "none", cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = F,
                                                               border_color = NA,
                                                               color = colorRampPalette(colors = c("lightgrey", "darkgreen"))(100),
                                                               annotation_col = sarcoma_frist_risk_module_mutated_samples_cancer_type,
                                                               annotation_colors = sarcoma_frist_risk_module_mutated_samples_cancer_type_colors,
                                                               fontsize_row = 6, fontsize_col = 6,
                                                               main = "Oncoplot of Sarcoma First ranked Risk Module Genes in the entire cohort")

ggsave(filename = "Results/Evaluations/Sarcoma first ranked Risk module/entire_cohort_sarcoma_frist_risk_module_mutated_samples_oncoplot1.pdf", 
       plot = entire_cohort_sarcoma_frist_risk_module_mutated_samples_oncoplot,
       device = "pdf", 
       width = 45, height = 12, units = "cm")

###########################################

## Fischer' exact test

################

### Method
# Set a threshold on mean_del_score (say 0.5).
# For each cancer type, count:
# The number of samples with any gene exceeding the threshold
# The number of samples with no genes exceeding the threshold
# Put all the numbers above into an n x 2 table (n being the number of cancer types).
# Perform Fisher's test on the table.

################

# Check if any of the genes exceed the threshold (0.5 del level) for each sample
entire_cohort_sarcoma_frist_risk_module_snv$any_sig_mutation <- 
sapply(1:nrow(entire_cohort_sarcoma_frist_risk_module_snv), function(i) {
  (entire_cohort_sarcoma_frist_risk_module_snv[i,c(1:4)] > 0.5) %>% sum() %>% as.logical()
})

entire_cohort_frist_risk_module_cancer_types <- unique(entire_cohort_sarcoma_frist_risk_module_snv$cancer_type)

### Defining a dataframe including
#### The number of samples with any gene exceeding the threshold
#### The number of samples with no genes exceeding the threshold
sarcoma_frist_risk_module_fischer <- data.frame(Num_mutated = vector(mode = "numeric", length = length(entire_cohort_frist_risk_module_cancer_types)), 
                                                Num_non_mutated = vector(mode = "numeric", length = length(entire_cohort_frist_risk_module_cancer_types)))

rownames(sarcoma_frist_risk_module_fischer) <- entire_cohort_frist_risk_module_cancer_types

lapply(1:length(entire_cohort_frist_risk_module_cancer_types), function(i) {
  tmp_table = entire_cohort_sarcoma_frist_risk_module_snv %>% filter(cancer_type == entire_cohort_frist_risk_module_cancer_types[i])
  sarcoma_frist_risk_module_fischer[entire_cohort_frist_risk_module_cancer_types[i],"Num_mutated"] <<- sum(tmp_table$any_sig_mutation)
  sarcoma_frist_risk_module_fischer[entire_cohort_frist_risk_module_cancer_types[i],"Num_non_mutated"] <<- nrow(tmp_table) - sum(tmp_table$any_sig_mutation)
})

## Fisher’s exact test is used when there is at least one cell in the contingency table of the expected frequencies below 5
### To retrieve the expected frequencies, use the chisq.test() function together with $expected
chisq.test(sarcoma_frist_risk_module_fischer)$expected

sarcoma_frist_risk_module_fischer_test <- fisher.test(sarcoma_frist_risk_module_fischer, workspace=2e8)

sarcoma_frist_risk_module_fischer_test

##########################################################################################
##########################################################################################
##########################################################################################

# Evaluation of the Risk results.

## Save the gene data and modules as supplementary files

tmp_risk_net_gene_tbl <- risk_net_gene_tbl[-c(1,9)]
names(tmp_risk_net_gene_tbl)[names(tmp_risk_net_gene_tbl) == "Rhabdoid"] <- "Rhabdoid carcinoma"

for(i in 1:10) {
  colnames(tmp_risk_net_gene_tbl[[i]]) <- c("Gene", "IVI", "Mean_deletriousness_score", "Risk_score", "Mean_neighborhood_weight", "Risk_mediator_score")
}

library(openxlsx)
write.xlsx(x = tmp_risk_net_gene_tbl, file = "Results/Evaluations/risk_net_gene_tbl.xlsx")

rm(tmp_risk_net_gene_tbl) 

############################

### save the modules

tmp_risk_net_modules_flt <- risk_net_modules_flt[-c(1,9)]
names(tmp_risk_net_modules_flt)[names(tmp_risk_net_modules_flt) == "Rhabdoid"] <- "Rhabdoid carcinoma"

tmp_risk_net_modules_flt_scores <- risk_net_modules_flt_scores[-c(1,9)]

tmp_risk_net_modules_flt <- lapply(tmp_risk_net_modules_flt, function(i) {
  
  tmp_data <- i
  tmp_data <- tmp_data[match(names(tmp_risk_net_modules_flt_scores[[which(sapply(tmp_risk_net_modules_flt_scores, length) == length(tmp_data))]]),
                             names(i))]
  
  tmp_data
    
})

for(i in 1:10) {
  tmp_risk_net_modules_flt[[i]] <-
  lapply(1:length(tmp_risk_net_modules_flt[[i]]), function(j) {
    tmp_data <- tmp_risk_net_modules_flt[[i]][[j]]
    colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
    tmp_data[tmp_data == 1] <- "Yes"
    tmp_data[tmp_data == 0] <- "No"
    tmp_data$Module <- paste0("Module number ", j)
    tmp_data <- data.frame(Module = unique(tmp_data$Module),
                           Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                           Mean_score = unique(tmp_data$Mean_score))
    tmp_data
  })
  tmp_risk_net_modules_flt[[i]] <- do.call(rbind, tmp_risk_net_modules_flt[[i]])
}

library(openxlsx)
write.xlsx(x = tmp_risk_net_modules_flt, file = "Results/Evaluations/risk_net_modules_flt.xlsx")

rm(tmp_risk_net_modules_flt, tmp_risk_net_modules_flt_scores)

##########################################################################################

## PeCan-based childhood cancer drivers

##########################################################################################

# convert Pedican hgg entrez gene names to symbols

pedican_hgg_driver_genes <-
gprofiler2::gconvert(query = c("596", "598", "581", "7157", "7157", "11186", "4609", "1956", "5159", "5156", "3265", "3845", 
                               "673", "2048", "1441", "5395", "2956", "4292", "4436", "142", "3417", "3418", "3020", 
                               "4763", "4603", "2735", "5728"), 
                     organism = "hsapiens", numeric_ns = "ENTREZGENE_ACC",
                     target = "HGNC", mthreshold = 1, filter_na = F)$target

hgg_driver_genes <- unique(c("TP53", "H3F3A", "TTN", "ATRX", "PDGFRA", "C1ORF112", "PIK3CA", "ACVR1", 
                      "ZBTB3", "LRRC7", "HIST1H3B", "NF1", "RNF213", "BCOR", "NEB", "LRCH3", 
                      "MUC16", "BRD1", "OBSCN", "CREBBP", pedican_hgg_driver_genes))

####################

pedican_sarcoma_driver_genes <-
  gprofiler2::gconvert(query = c("2308", "7852", "5077", "5081", "7021", "1958", "5156", "1029", "7157", "1653", "4613", "407975", "6839", "1956", "5584", 
                                 "472", "2264", "1001", "7434", "7054", "4613", "1630", "4267", "2130", "1029", "1026", "5925", "2113", "4137", "5502", "4751", 
                                 "595", "2064", "466", "2065", "3084", "2254", "2260", "3678", "1017", "843", "7157", "23095", "3439", "3458", "7124", "8743", 
                                 "10370", "2719", "5395", "11200", "885", "4790", "190", "3298", "6657", "3479", "3480", "23435", "1959", "332", "2737", "4193", 
                                 "3265", "4893", "3845", "406937", "2308", "4613", "4763", "4193", "1019", "3481", "545", "5726", "7157", "1029", "6598", "5925", 
                                 "7054", "5081", "5727", "1630", "7248", "7249", "574", "7490", "11186", "4267", "23095", "1030", "8312", "4255", "7852", "5077", 
                                 "5081", "2308", "2309", "114625", "6007", "11200", "6495", "595", "4609", "7430", "5077", "4609", "3918", "5058", "3815", "2262", 
                                 "2247", "3082", "5243", "4363", "4035", "7021", "7157", "1029", "1026", "5155", "7157", "5925", "2064", "843", "5925", "1029", 
                                 "1026", "581", "3815", "11186", "7078", "4255", "1612", "3479", "355", "356", "5396", "7158", "9402", "9401", "6696", "860", "3381", 
                                 "80005", "4045", "596", "4170", "1869", "1871", "2099", "4193", "2736", "407021", "2735", "7157", "7158", "6495", "595", "4609", 
                                 "7430", "7157", "7030", "4613"), 
                       organism = "hsapiens", numeric_ns = "ENTREZGENE_ACC",
                       target = "HGNC", mthreshold = 1, filter_na = F)$target

sarcoma_driver_genes <- unique(c("TP53", "VNN3", "RB1", "ATRX", "DMD", "TTN", "CSMD1", "MUC16", "FER", "DGKB", "LAMA5", "AUTS2", "DLG2", "NAALADL2", "SYT16", 
                                 "MIRLET7BHG", "EYS", "LRRC7", "OTOGL", "C1ORF112", "EWSR1", "FLI1", "STAG2", "ERG", "AP1B1", "TTC28", "TP53", "EEF2", "ZNRF3", 
                                 "RUNX1", "EMID1", "ERF", "CSMD1", "SLC35F4", "RYR3", "MUC16", "DEPDC5", "ETS1", "LRRC7", "CACNG2", "ASCC1", "FOXO1", "SYT16", 
                                 "TP53", "LRP1B", "PAX3", "PSMG3-AS1", "FGFR4", "GRM1", "PIK3CA", "SYNE1", "BCOR", "SCRIB", "MUC6", "NOTCH1", "LRRC7", "TUG1", 
                                 "LOC100499484-C9ORF174", "SEC24A", "ACSM2A", "PERM1", "KRTAP4-8", "UBXN11", "ITGAM", "ZNF891", "NF1", "ADCY8", "PTPN11", "C17ORF97", 
                                 "BCOR", "TTBK1", "PCCA", "MN1", "OTOP1", "DMD", "KDM5C", "GPR161", "POLR1A", "ADGRV1", "TENM2", pedican_sarcoma_driver_genes))

####################

pedican_leukaemia_driver_genes <-
  gprofiler2::gconvert(query = c("4297", "4299", "8842", "1029", "595", "4869", "2324", "3270", "2120", "1027", "2842", "7490", "7490", "3492", "9", "10", 
                                 "2944", "1543", "4297", "7157", "2322", "6929", "4613", "5087", "10397", "7124", "4049", "7490", "1728", "4353", "1571", 
                                 "9760", "23532", "4193", "4790", "3077", "2950", "373156", "9429", "2353", "5243", "5578", "843", "23532", "54111", "3953", 
                                 "3195", "3251", "2263", "4524", "4609", "5111", "2067", "3726", "596", "7187", "5310", "581", "836", "1030", "472", "6545", 
                                 "861", "10962", "3845", "4893", "2120", "4088", "2120", "861", "1565", "2952", "100507249", "3439", "7157", "7520", "5925", 
                                 "4683", "23368", "7163", "7163", "6929", "1029", "1030", "1026", "1027", "5243", "10628", "7704", "5594", "2957", "6821", 
                                 "4297", "4524", "7298", "4397", "4552", "55294", "4851", "1747", "3479", "4363", "4035", "978", "6598", "892", "3315", 
                                 "6275", "1026", "5079", "27086", "2149", "51561", "7172", "1442", "2068", "7515", "30012", "3197", "166824", "644943", 
                                 "9429", "207", "4050", "6382", "862", "3717", "3305", "3303", "3304", "7535", "54363", "2037", "6936", "4124", "84295", 
                                 "841", "355", "3480", "7159", "84159", "4193", "8841", "51564", "9734", "7491", "3110", "23286", "2200", "3490", "2623", 
                                 "5925", "79718", "1387", "4585", "8623", "105", "84296", "578", "3574", "7012", "92912", "22877", "7490", "356", "9759", 
                                 "51564", "5836", "5142", "1493", "912", "3575", "4040", "79370", "80824", "1389", "4602", "5243", "2272", "1788", "10320", 
                                 "440", "1053", "64109", "3717", "7153", "6886", "1031", "3251", "4297", "4763", "5243", "3458", "1438", "3562", "947", 
                                 "7157", "1464", "1029", "3845", "5925", "406938", "406912"), 
                       organism = "hsapiens", numeric_ns = "ENTREZGENE_ACC",
                       target = "HGNC", mthreshold = 1, filter_na = F)$target

leukaemia_driver_genes <- unique(c("ETV6", "PAX5", "RUNX1", "IKZF1", "NRAS", "KRAS", "TCF3", "CDKN2A", "PBX1", "ERG", "CRLF2", "TP53", 
                                   "ABL1", "FLT3", "CREBBP", "DGKB", "TBL1XR1", "KMT2A", "JAK2", "BCR", "NOTCH1", "CDKN2A", "CDKN2B", "TAL1", 
                                   "FBXW7", "PHF6", "PTEN", "WT1", "MYB", "MLLT3", "BCL11B", "DNM2", "STIL", "LEF1", "FLT3", "USP7", "ETV6", "MYC", 
                                   "KDM6A", "NRAS", "RUNX1", "NRAS", "RUNX1T1", "KMT2A", "MYH11", "CBFB", "KIT", "CSNKA2IP", "FLT3", "WT1", "MLLT3", 
                                   "CBL", "KRAS", "MLLT10", "ASXL2", "NUP98", "MYC", "TTN", "TET2", "NF1", "KMT2A", "AFF1", "BTD", "FLT3", "NRAS", "MLLT10", 
                                   "MN1", "FOXO1", "PRDM16", "NOTCH1", "TJP1", "CP", "WT1", "NECTIN2", "TET3", "RPSAP58", "NUDT6", "FGF2", "CAND2", "PARD3B", 
                                   pedican_leukaemia_driver_genes))

####################

pedican_nbl_driver_genes <-
  gprofiler2::gconvert(query = c("7161", "207", "4613", "7054", "8643", "5272", "841", "4763", "4192", "2051", "1948", "1949", "6752", "7481", "2535", 
                                 "324", "4904", "5204", "9168", "8797", "8326", "5979", "2668", "2674", "1910", "1908", "1996", "1995", "8568", "7490", 
                                 "843", "54111", "4613", "23095", "274", "1029", "8312", "7057", "4914", "1270", "4363", "5728", "1755", "429", "3280", 
                                 "8929", "23261", "4763", "29108", "894", "3205", "378708", "7157", "1026", "960", "150465", "650", "4601", "4153", "9063", 
                                 "3066", "4100", "4102", "1485", "11337", "1030", "7474", "5588", "581", "836", "11186", "51364", "1641", "3265", "11200", 
                                 "5243", "5290", "4609", "3845", "4893", "673", "1385", "3479", "84432", "3297", "1499", "7471", "5925", "7709", "5789", 
                                 "3815", "3091", "4193", "4915", "238", "401", "92140", "2247", "3399", "3398", "3397", "1017", "100302292", "3458", 
                                 "148753", "604", "407975", "1969", "4487", "22943", "27123", "27122", "6422", "4610", "4790", "51752", "64167", "580", 
                                 "2064", "79923", "5460", "4751", "3266", "2778", "8842", "4804", "4004", "54897", "6337", "112464", "3880", "10563", 
                                 "643", "1027", "596", "10257", "8714", "538", "2146", "1191", "864", "3135", "1601", "129807", "6790", "2066", "2309", 
                                 "472", "5629", "57448", "1789", "57531", "389421", "4000", "50614", "5727", "4853", "9787", "79075", "6608", "6632", 
                                 "6742", "11065", "2034", "3481", "558", "4854", "1441", "599", "842", "55665", "5465", "2597", "283120", "3558", 
                                 "1956", "5978", "3480", "4105", "4832", "2026", "6285", "5409", "27087", "7345", "4684", "9464", "407040"), 
                       organism = "hsapiens", numeric_ns = "ENTREZGENE_ACC",
                       target = "HGNC", mthreshold = 1, filter_na = F)$target

nbl_driver_genes <- unique(c("MYCN", "ALK", "TERT", "ATRX", "FAM49A", "MUC16", "PTPRD", "TTN", "DDX1", "SHANK2", "LRRC7", 
                             "ODF2", "ASCC1", "SLC12A7", "MYCNOS", "AUTS2", "MUC17", "SLIT3", "SLC6A18", "LRP1B", pedican_nbl_driver_genes))

####################

pedican_embryonal_cns_driver_genes <-
  gprofiler2::gconvert(query = c("8522", "7157", "4609", "1026", "1956"), 
                       organism = "hsapiens", numeric_ns = "ENTREZGENE_ACC",
                       target = "HGNC", mthreshold = 1, filter_na = F)$target

embryonal_cns_driver_genes <- unique(c("HLA-DQB1", "OR4A16", "ATP8B2", "HLA-B", "LILRA6", "SF3A2", "MESP2", "TMEM235", 
                                       "GBP4", "KRTAP4-5", "IGFN1", "TNRC6A", "PCDHA8", "HLA-DRB1", "KRTAP5-2", 
                                       "FBLN1", "MUC17", "FAM174B", "ADGRB1", "MMP17", pedican_embryonal_cns_driver_genes))

####################

pedican_epd_driver_genes <-
  gprofiler2::gconvert(query = c("7157", "4613", "4609", "1027", "5156", "1956", "2066", "7015", "4916", "7422", "4771", 
                                 "23136", "2035", "2037", "6598", "595", "1019", "664", "637", "596", "5341", "4791", "374491", 
                                 "324", "7471", "3371", "4851", "3280", "23493", "55294", "10481", "4747", "5156", 
                                 "5599", "5468", "144", "7249"), 
                       organism = "hsapiens", numeric_ns = "ENTREZGENE_ACC",
                       target = "HGNC", mthreshold = 1, filter_na = F)$target

epd_driver_genes <- unique(c("C11ORF95", "RELA", "IGF2R", "PISD", "DNAH8", "TUG1", "MICA", "CDC5L", "BRCA1", "NF1", "KRTAP4-5", 
                             "HLA-DQA1", "MID1", "LAMA5", "TSC2", "NOTCH2", "TTN", "LINC00935", "TUBGCP2", "SOX13", pedican_epd_driver_genes))

####################

pedican_lymphoma_driver_genes <-
  gprofiler2::gconvert(query = c("596", "3492", "7157", "604", "7054", "472", "4683", "472", "4613", "7157", 
                                 "238", "3251", "4609", "7157", "4193", "1029", "8793", "22974", "4524", "5395", 
                                 "541466", "4609", "596", "604", "2956", "4609", "4684", "3662", "604", "238", 
                                 "472", "5396", "596", "6693", "596", "9260", "1029", "1031", "4763"), 
                       organism = "hsapiens", numeric_ns = "ENTREZGENE_ACC",
                       target = "HGNC", mthreshold = 1, filter_na = F)$target

lymphoma_driver_genes <- unique(c("MYC", "ID3", "TP53", "FBXO11", "CCND3", "IGLL5", "SMARCA4", "RHOA", "HLA-DQB1", 
                                  "P2RY8", "PCBP1", "HLA-DRB1", "TFAP4", "BMS1P20", "RYR2", "RFX7", "CARD11", "FOXO1", "ZAN", "DDX3X", 
                                  pedican_lymphoma_driver_genes))

####################

pedican_mb_driver_genes <-
  gprofiler2::gconvert(query = c("3156", "6598", "3213", "3214", "3223", "4609", "5015", "7545", "7157", "1029", "5727", 
                                 "324", "1030", "7249", "1387", "2737", "6469", "843", "11186", "22931", "8314", 
                                 "59336", "1027", "1032", "1031", "6422", "8643", "474", "100133941", "841", "6790", 
                                 "5925", "4255", "5290", "4771", "3815", "5728", "7490", "3728", "1499", "3090", 
                                 "4916", "147040", "4613", "1936", "6156", "6224", "4601", "51684", "8945"), 
                       organism = "hsapiens", numeric_ns = "ENTREZGENE_ACC",
                       target = "HGNC", mthreshold = 1, filter_na = F)$target

mb_driver_genes <- unique(c("PTCH1", "DDX3X", "KMT2D", "CTNNB1", "KMT2C", "SMARCA4", "MYCN", "TTN", "MYC", "KBTBD4", 
                            "KDM6A", "TP53", "OTX2", "CREBBP", "MUC16", "GLI2", "ZMYM3", "SMO", "ZIC1", "TCF4", "PTEN", 
                            "GSE1", "FBXW7", "PIK3CA", "CTDNEP1", "GFI1B", "BCOR", "SYNE1", "TERT", "PRKAR1A", "NEB", 
                            "MUC17", "STAG2", "CDK6", "OBSCN", "RYR2", "RYR3", "SUFU", "ARID1A", "SCN4A", "KDM3B", "CHD7", 
                            "DNAH14", "FAT1", "FSIP2", "LRP1B", "PCDH19", "CSF2RA", "MED12", "LDB1", "FMR1", 
                            "FCGBP", "APOB", "IDH1", "MACF1", pedican_mb_driver_genes))

####################

pedican_rhabdoid_driver_genes <-
  gprofiler2::gconvert(query = c("6598", "7157", "5925", "1029", "5728", "4609", "595", "841", "6598", 
                                 "6599", "9074", "5274", "5460", "79923", "6657", "207", "3417", "6598", "6598"), 
                       organism = "hsapiens", numeric_ns = "ENTREZGENE_ACC",
                       target = "HGNC", mthreshold = 1, filter_na = F)$target

rhabdoid_driver_genes <- unique(c("MYCN", "SMARCB1", "ATP8B2", "FMN2", "PABPC1", "RPTN", "TAF5L", "STARD9", "RP1L1", "PCBP2", 
                                  "OR4C6", "MUC16", "KRTAP1-3", "ARMCX4", "KRT39", "POU3F3", "ZAK", "MSH3", "LGSN", "SCARA5",
                                  pedican_rhabdoid_driver_genes))


##########################################################################################

## Evaluation of the top 20 genes

### DMG
hgg_driver_genes[
which(
  hgg_driver_genes 
  %in% 
(risk_net_gene_tbl$DMG %>% 
  slice_max(order_by = primitive_risk_score, n = 20) %>% 
  select(gene) %>% 
  unname() %>% 
  unlist())
)
] %>% cat()

#####################

### HGG
hgg_driver_genes[
  which(
    hgg_driver_genes 
    %in% 
      (risk_net_gene_tbl$HGG %>% 
         slice_max(order_by = primitive_risk_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### Leukaemia
leukaemia_driver_genes[
  which(
    leukaemia_driver_genes 
    %in% 
      (risk_net_gene_tbl$Leukaemia %>% 
         slice_max(order_by = primitive_risk_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### NBL
nbl_driver_genes[
  which(
    nbl_driver_genes 
    %in% 
      (risk_net_gene_tbl$NBL %>% 
         slice_max(order_by = primitive_risk_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### Sarcoma
sarcoma_driver_genes[
  which(
    sarcoma_driver_genes 
    %in% 
      (risk_net_gene_tbl$Sarcoma %>% 
         slice_max(order_by = primitive_risk_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### CNS embryonal
embryonal_cns_driver_genes[
  which(
    embryonal_cns_driver_genes 
    %in% 
      (risk_net_gene_tbl$`CNS embryonal` %>% 
         slice_max(order_by = primitive_risk_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### EPD
epd_driver_genes[
  which(
    epd_driver_genes 
    %in% 
      (risk_net_gene_tbl$EPD %>% 
         slice_max(order_by = primitive_risk_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### Lymphoma
lymphoma_driver_genes[
  which(
    lymphoma_driver_genes 
    %in% 
      (risk_net_gene_tbl$Lymphoma %>% 
         slice_max(order_by = primitive_risk_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### MB (Medulloblastoma)
mb_driver_genes[
  which(
    mb_driver_genes 
    %in% 
      (risk_net_gene_tbl$MB %>% 
         slice_max(order_by = primitive_risk_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### Rhabdoid
rhabdoid_driver_genes[
  which(
    rhabdoid_driver_genes 
    %in% 
      (risk_net_gene_tbl$Rhabdoid %>% 
         slice_max(order_by = primitive_risk_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

##########################################################################################

## Evaluation of the 1st-ranked modules of all cancer types

### DMG
hgg_driver_genes[
  which(
    hgg_driver_genes 
    %in% 
      na.omit(risk_net_1st_modules$DMG)
  )
] %>% cat()

#####################

### HGG
hgg_driver_genes[
  which(
    hgg_driver_genes 
    %in% 
      na.omit(risk_net_1st_modules$HGG)
  )
] %>% cat()

#####################

### Leukaemia
leukaemia_driver_genes[
  which(
    leukaemia_driver_genes 
    %in% 
      na.omit(risk_net_1st_modules$Leukaemia)
  )
] %>% cat()

#####################

### NBL
nbl_driver_genes[
  which(
    nbl_driver_genes 
    %in% 
      na.omit(risk_net_1st_modules$NBL)
  )
] %>% cat()

#####################

### Sarcoma
sarcoma_driver_genes[
  which(
    sarcoma_driver_genes 
    %in% 
      na.omit(risk_net_1st_modules$Sarcoma)
  )
] %>% cat()

#####################

### CNS embryonal
embryonal_cns_driver_genes[
  which(
    embryonal_cns_driver_genes 
    %in% 
      na.omit(risk_net_1st_modules$`CNS embryonal`)
  )
] %>% cat()

#####################

### EPD
epd_driver_genes[
  which(
    epd_driver_genes 
    %in% 
      na.omit(risk_net_1st_modules$EPD)
  )
] %>% cat()

#####################

### Lymphoma
lymphoma_driver_genes[
  which(
    lymphoma_driver_genes 
    %in% 
      na.omit(risk_net_1st_modules$Lymphoma)
  )
] %>% cat()

#####################

### MB
mb_driver_genes[
  which(
    mb_driver_genes 
    %in% 
      na.omit(risk_net_1st_modules$MB)
  )
] %>% cat()

#####################

### Rhabdoid
rhabdoid_driver_genes[
  which(
    rhabdoid_driver_genes 
    %in% 
      na.omit(risk_net_1st_modules$Rhabdoid)
  )
] %>% cat()

##########################################################################################
##########################################################################################

## Literature mining of cancer risk results
library(pubmedR)

lit_risk_data <- tibble(`Cancer type` = vector("character"), 
                        `Gene group` = vector("character"),
                        Gene = vector("character"),
                        Num_PMIDs = vector("integer")
)

for(i in 1:length(risk_net_gene_tbl[-c(1,9)])) {
  
  if(names(risk_net_gene_tbl[-c(1,9)])[i] == "HGG") {
    main_search_text <- "((HGG[Title/Abstract]) OR (HGGs[Title/Abstract]) OR (High-Grade Glioma[Title/Abstract]) OR (High-Grade Gliomas[Title/Abstract]))"
    
  } else if(names(risk_net_gene_tbl[-c(1,9)])[i] == "DMG") {
    main_search_text <- "((DMG[Title/Abstract]) OR (DMGs[Title/Abstract]) OR (Diffuse Midline Glioma[Title/Abstract]) OR (Diffuse Midline Gliomas[Title/Abstract]))"
    
  } else if(names(risk_net_gene_tbl[-c(1,9)])[i] == "Leukaemia") {
    main_search_text <- "((leukemia[Title/Abstract]) OR (leukemias[Title/Abstract]) OR (Leukaemias[Title/Abstract]) OR (Leukaemia[Title/Abstract]))"
    
  } else if(names(risk_net_gene_tbl[-c(1,9)])[i] == "Lymphoma") {
    main_search_text <- "((Lymphoma[Title/Abstract]) OR (Lymphomas[Title/Abstract]))"
    
  } else if(names(risk_net_gene_tbl[-c(1,9)])[i] == "NBL") {
    main_search_text <- "((NBL[Title/Abstract]) OR (NBLs[Title/Abstract]) OR (Neuroblastoma[Title/Abstract]) OR (Neuroblastomas[Title/Abstract]))"
    
  } else if(names(risk_net_gene_tbl[-c(1,9)])[i] == "Sarcoma") {
    main_search_text <- "((Sarcoma[Title/Abstract]) OR (Sarcomas[Title/Abstract]))"
    
  } else if(names(risk_net_gene_tbl[-c(1,9)])[i] == "CNS embryonal") {
    main_search_text <- "((atypical teratoid/rhabdoid tumor[Title/Abstract]) OR (AT/RT[Title/Abstract]) OR (embryonal tumors with multilayered rosettes[Title/Abstract]) OR (ETMR[Title/Abstract]) OR (ependymoblastoma[Title/Abstract]))"
    
  } else if(names(risk_net_gene_tbl[-c(1,9)])[i] == "EPD") {
    main_search_text <- "((EPD[Title/Abstract]) OR (Ependymoma[Title/Abstract]) OR (Ependymomas[Title/Abstract]) OR (Subependymoma[Title/Abstract]))"
    
  } else if(names(risk_net_gene_tbl[-c(1,9)])[i] == "MB") {
    main_search_text <- "((MB[Title/Abstract]) OR (medulloblastoma[Title/Abstract]) OR (medulloblastomas[Title/Abstract]))"
    
  } else if(names(risk_net_gene_tbl[-c(1,9)])[i] == "Rhabdoid") {
    main_search_text <- "((Rhabdoid[Title/Abstract]) OR (Renal medullary carcinoma[Title/Abstract]) OR (RMC[Title/Abstract]))"
  }
  
  ### Add top 20 genes
  tmp_top20_tbl <- tibble(`Cancer type` = names(risk_net_gene_tbl[-c(1,9)])[i],
                          `Gene group` = "Top 20 genes",
                          Gene = (risk_net_gene_tbl[-c(1,9)][[i]] %>%
                                    slice_max(order_by = primitive_risk_score, n = 20) %>%
                                    select(gene) %>%
                                    unname() %>%
                                    unlist()))
  
  tmp_top20_pmid_list <- lapply(1:20, function(j) {
    query_text <- paste0("(", tmp_top20_tbl$Gene[j], "[Title/Abstract]) AND ", main_search_text)
    tmp_entrez <- pmQueryTotalCount(query = query_text, api_key = "1a1d2b704b3bad87d1ded4f2d95d611e1909")
    tmp_entrez_num <- tmp_entrez$total_count
    print(paste0("Top gene ", j, " is done!"))
    return(tmp_entrez_num)
  })
  
  tmp_top20_tbl$Num_PMIDs <- unlist(tmp_top20_pmid_list)
  
  lit_risk_data <- rbind(lit_risk_data, tmp_top20_tbl)
  
  ### Add the 1st ranked module
  tmp_1st_module_tbl <- tibble(`Cancer type` = names(risk_net_gene_tbl[-c(1,9)])[i],
                               `Gene group` = "1st-ranked module",
                               Gene = (na.omit(risk_net_1st_modules[,-c(1,9)][,i])))
  
  tmp_1st_module_pmid_list <- lapply(1:length(na.omit(risk_net_1st_modules[,-c(1,9)][,i])), function(m) {
    query_text <- paste0("(", tmp_1st_module_tbl$Gene[m], "[Title/Abstract]) AND ", main_search_text)
    tmp_entrez <- pmQueryTotalCount(query = query_text, api_key = "1a1d2b704b3bad87d1ded4f2d95d611e1909")
    tmp_entrez_num <- tmp_entrez$total_count
    print(paste0("Module 1 gene ", m, " is done!"))
    return(tmp_entrez_num)
  })
  
  tmp_1st_module_tbl$Num_PMIDs <- unlist(tmp_1st_module_pmid_list)
  
  lit_risk_data <- rbind(lit_risk_data, tmp_1st_module_tbl)
  
  print(paste0(names(risk_net_gene_tbl[-c(1,9)])[i], " is done!"))
  
}

#################################

## Visualization of the results

lit_risk_data[lit_risk_data == "Rhabdoid"] <- "Rhabdoid carcinoma"

risk_lit_mining_plot <- ggplot(data = lit_risk_data, aes(x = `Cancer type`, 
                                                         y = log2(Num_PMIDs+1), 
                                                         fill = `Gene group`)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.01, width=0.58, size = 0.2) +
  labs(y= 'Log2(number of papers + 1)',
       title = "Literature-based Association of Top\nRisk Genes with Their Respective Cancers") +
  coord_flip() +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_fill_fish_d(option = "Elagatis_bipinnulata") +
  # scale_fill_viridis_d(begin = 0.32, end = 0.83) +
  theme_bw() +
  theme(text = element_text(size = 7),
        legend.position = "right", 
        legend.box.background = element_rect(colour = "black"))

risk_lit_mining_plot

##########################################################################################
##########################################################################################

## ORA on the first ranked risk modules

library(enrichR)

### Using enrichR

#Get the list of all databases
enrichR_dbs <- enrichR::listEnrichrDbs()
enrichR_dbs %>% as.data.frame() %>% view()

#Select you desired databases
dbs.GO.BP <- "GO_Biological_Process_2023"
# dbs.KEGG.pathway <- "KEGG_2021_Human"

#ORA on GO-BP terms for 1st risk module of each cancer type

risks_1st_module_ORA <- lapply(1:10, function(i) {
  enrichR::enrichr(genes = as.character(na.omit(risk_net_1st_modules[,-c(1,9)][,i])),
                   databases = dbs.GO.BP)[[1]]
})

names(risks_1st_module_ORA) <- colnames(risk_net_1st_modules[,-c(1,9)])

###########################

### save the results
library(openxlsx)

write.xlsx(x = risks_1st_module_ORA, file = "Results/Evaluations/ORA/risks/risks_1st_module_ORA.xlsx")

###########################

### visualize the results
risks_1st_module_ORA_vis_list <- lapply(1:10, function(i) {
  
  #Filter the data based on their Combined Score
  tmp_ora <- risks_1st_module_ORA[[i]] %>% 
    slice_max(Combined.Score, n = 5) 
  
  #Correcting the levels
  tmp_ora$Term <- factor(tmp_ora$Term, levels = c(as.character(tmp_ora$Term)[order(tmp_ora$Combined.Score)]))
  
  #Visualization of terms
  tmp_ora_plot <- ggplot(tmp_ora, aes(x = Term, y = Combined.Score,
                                      color = Adjusted.P.value)) +  
    geom_point(size = 2) +
    ylab("Combined Score") + xlab(NULL) + 
    scale_color_viridis_c(name="Adjusted\nP-value", breaks= scales::pretty_breaks()) + 
    coord_flip() +
    labs(title = expression("ORA of the"~"1"^"st"~ "Ranked Risk Module of"),
         subtitle = colnames(risk_net_1st_modules[,-c(1,9)])[i]) +
    scale_y_continuous(breaks= scales::pretty_breaks(), labels = function(x) format(x, scientific = TRUE)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) +
    theme_classic() + theme(
      legend.position = "right",
      # legend.position = c(.95, 0.02),
      # legend.justification = c("right", "bottom"),
      legend.box.background = element_rect(color="black", size=1),
      text = element_text(size = 7))
  
  #Saving the plot
  ggsave(filename = paste0("Results/Evaluations/ORA/risks/", colnames(risk_net_1st_modules[,-c(1,9)])[i], "_risks_1st_module_ORA.pdf"), 
         plot = tmp_ora_plot,
         device = "pdf",width = 9.5, height = 4.6, units = "cm")
  
  tmp_ora_plot
  
})

risks_1st_module_ORA_vis_list[[3]]

#=============================================================================
#
#    Code chunk 8: Reconstruction and analysis of the pan-cancer Risk Net
#
#=============================================================================

pan_cancer_risk_net <- pan_germline_snv_cor

# Prepare the table for merging with other data
pan_cancer_risk_net <- pan_cancer_risk_net[,c(1,2)]
colnames(pan_cancer_risk_net) <- c("from", "to")
pan_cancer_risk_net$type <- "co-deletriousness"

# Add the co-expression data to the table
pan_cancer_risk_net_genes <- pan_cancer_risk_net[,c(1,2)] %>% 
  unlist() %>% 
  unname() %>% 
  unique()
  
pan_cancer_risk_net_first_level_coex_index <- c(which(pan_normal_rnaseq_cor$row %in% pan_cancer_risk_net_genes), 
                                                which(pan_normal_rnaseq_cor$column %in% pan_cancer_risk_net_genes)) %>% unique()
  
  
pan_cancer_risk_net_first_level_coex_genes <- c(pan_cancer_risk_net_genes, 
                                                pan_normal_rnaseq_cor$row[pan_cancer_risk_net_first_level_coex_index],
                                                pan_normal_rnaseq_cor$column[pan_cancer_risk_net_first_level_coex_index]) %>% unique()

pan_cancer_risk_net_first_two_level_coex_index <- c(which(pan_normal_rnaseq_cor$row %in% pan_cancer_risk_net_first_level_coex_genes), 
                                                    which(pan_normal_rnaseq_cor$column %in% pan_cancer_risk_net_first_level_coex_genes)) %>% unique()

pan_normal_rnaseq_cor_flt <- pan_normal_rnaseq_cor[pan_cancer_risk_net_first_two_level_coex_index,c(1,2)]
colnames(pan_normal_rnaseq_cor_flt) <- c("from", "to")
pan_normal_rnaseq_cor_flt$type <- "co-expression"

pan_cancer_risk_net <- rbind(pan_cancer_risk_net, pan_normal_rnaseq_cor_flt)

rm(pan_normal_rnaseq_cor_flt)

# Add PPI data to the Risk net
library(STRINGdb)

# Map gene names to stringDB ids for PPI analysis

## get the gene names in each risk net
pan_cancer_risk_net_genes_plus_coex <- 
  data.frame(gene_symbol = pan_cancer_risk_net[,c(1,2)] %>% 
               unlist() %>% 
               unname() %>% 
               unique())

## getSTRINGdb for human
string_db <- STRINGdb$new(species=9606, version = "11.5", 
                          score_threshold=900, # very high confidence (>0.9)
                          input_directory=""
)

pan_cancer_risk_net_genes_plus_coex <- string_db$map(my_data_frame = pan_cancer_risk_net_genes_plus_coex, 
                                                     my_data_frame_id_col_names = "gene_symbol", 
                                                     takeFirst = TRUE, removeUnmappedRows = TRUE)

# Get the PPIs
pan_cancer_risk_net_ppi <- string_db$get_interactions(string_ids = pan_cancer_risk_net_genes_plus_coex$STRING_id)

# Convert stringdb ids to gene symbols

pan_cancer_risk_net_ppi_row_index <- match(pan_cancer_risk_net_ppi$from, pan_cancer_risk_net_genes_plus_coex$STRING_id)
pan_cancer_risk_net_ppi$from <- pan_cancer_risk_net_genes_plus_coex$gene_symbol[pan_cancer_risk_net_ppi_row_index]

pan_cancer_risk_net_ppi_column_index <- match(pan_cancer_risk_net_ppi$to, pan_cancer_risk_net_genes_plus_coex$STRING_id)
pan_cancer_risk_net_ppi$to <- pan_cancer_risk_net_genes_plus_coex$gene_symbol[pan_cancer_risk_net_ppi_column_index]

pan_cancer_risk_net_ppi <- pan_cancer_risk_net_ppi[,c(1,2)]
pan_cancer_risk_net_ppi$type <- "PPI"

## Remove duplicates
pan_cancer_risk_net_ppi <- pan_cancer_risk_net_ppi[-which(duplicated(paste(pan_cancer_risk_net_ppi$from,
                                                                           pan_cancer_risk_net_ppi$to, 
                                                                           sep = "_"))),]

#####################################

# Combine PPI with the original risk net
pan_cancer_risk_net <- rbind(pan_cancer_risk_net, pan_cancer_risk_net_ppi)

# Reconstruct the networks
pan_cancer_risk_net <- igraph::graph_from_data_frame(pan_cancer_risk_net, directed=FALSE)

## Check the network layers
table(igraph::E(pan_cancer_risk_net)$type)

###############################

# Calculate the IVI values of genes within each net
library(influential)
pan_cancer_risk_net_ivi <- ivi(pan_cancer_risk_net)

###############################

# Calculate the primitive gene risk scores
pan_cancer_risk_net_gene_tbl <- data.frame(gene = names(pan_cancer_risk_net_ivi), ivi = pan_cancer_risk_net_ivi)

## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer

pan_cancer_germline_snv_mean_deletriousness <- colMeans(pan_germline_snv_tbl)

## Generate the gene table and calculate the primitive scores
pan_cancer_risk_net_gene_tbl$mean_deletriousness <- 0
pan_cancer_risk_net_gene_tbl_match_index <- which(names(pan_cancer_germline_snv_mean_deletriousness) %in% rownames(pan_cancer_risk_net_gene_tbl))

pan_cancer_risk_net_gene_tbl[names(pan_cancer_germline_snv_mean_deletriousness)[pan_cancer_risk_net_gene_tbl_match_index], "mean_deletriousness"] <- 
pan_cancer_germline_snv_mean_deletriousness[pan_cancer_risk_net_gene_tbl_match_index]
  
pan_cancer_risk_net_gene_tbl$primitive_risk_score <- (pan_cancer_risk_net_gene_tbl$ivi)*(pan_cancer_risk_net_gene_tbl$mean_deletriousness) # node weights

# add gene scores
pan_cancer_risk_net <- set_vertex_attr(graph = pan_cancer_risk_net, name = "score", value = pan_cancer_risk_net_gene_tbl$primitive_risk_score)

# Convert the un-weighted to a node-weighted network
pan_cancer_risk_net <- set_vertex_attr(graph = pan_cancer_risk_net, name = "weight", value = (pan_cancer_risk_net_gene_tbl$ivi)^(pan_cancer_risk_net_gene_tbl$mean_deletriousness))

## Check the network node weights
summary(igraph::V(pan_cancer_risk_net)$weight)

# Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
  pan_cancer_risk_net_gene_tbl$mean_neighborhood_weight <- sapply(1:nrow(pan_cancer_risk_net_gene_tbl), function(j) {
    igraph::neighbors(graph = pan_cancer_risk_net, 
                      v = V(pan_cancer_risk_net)[j])$score %>% mean()
  })
  
pan_cancer_risk_net_gene_tbl$risk_mediator_score <- pan_cancer_risk_net_gene_tbl$primitive_risk_score * pan_cancer_risk_net_gene_tbl$mean_neighborhood_weight

## Top 20 Risk genes
pan_cancer_risk_net_gene_tbl %>% 
  slice_max(order_by = primitive_risk_score, n = 20) %>% 
  mutate(Rank = rank(-(primitive_risk_score))) %>% 
  select(Rank, gene) %>% 
  rename(Gene = gene) %>% 
  View()


## Top 30 mediators
pan_cancer_risk_net_gene_tbl %>% 
  slice_max(order_by = risk_mediator_score, n = 30) %>% 
  mutate(Rank = rank(-(risk_mediator_score))) %>% 
  select(Rank, gene) %>% 
  rename(Gene = gene) %>% 
  View()

# Detect communities/modules/clusters
set.seed(3847)
pan_cancer_risk_net_modules <- 
  igraph::cluster_leiden(
    graph = pan_cancer_risk_net,
    objective_function = "CPM",
    weights = NULL,
    resolution_parameter = 0.5,
    beta = 0.05,
    initial_membership = NULL,
    n_iterations = 100,
    vertex_weights = igraph::V(pan_cancer_risk_net)$weight
  )

sizes(pan_cancer_risk_net_modules) %>% as.vector() %>% summary()
length(pan_cancer_risk_net_modules)

## Inspect and filter modules

### Create a table of modules and their genes
  pan_cancer_risk_net_modules_tbl <- data.frame(
    module = membership(pan_cancer_risk_net_modules) %>% as.integer(),
    gene = pan_cancer_risk_net_modules$names
  )

### Add if the gene is mutated or not (could be a risk gene or not)
pan_cancer_risk_net_modules_tbl$mutated <- 0
pan_cancer_risk_net_modules_tbl$mutated[which(pan_cancer_risk_net_gene_tbl$mean_deletriousness > 0)] <- 1

### Separate modules
pan_cancer_risk_net_modules_tbl <- lapply(unique(pan_cancer_risk_net_modules_tbl$module), function(j){
  subset(pan_cancer_risk_net_modules_tbl, module == j)
})

names(pan_cancer_risk_net_modules_tbl) <- sapply(pan_cancer_risk_net_modules_tbl, function(m) {
  unique(m$module)
})

### Remove modules that have less than 2 Risk genes or their size is less than 4
pan_cancer_risk_net_modules_tbl <- lapply(pan_cancer_risk_net_modules_tbl, function(k) {
  if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
    k
  } else {
    NULL
  }
})

### Remove NULL modules
pan_cancer_risk_net_modules_tbl[which(sapply(pan_cancer_risk_net_modules_tbl, is.null))] <- NULL

### Calculate the module scores
pan_cancer_risk_net_modules_tbl <- lapply(pan_cancer_risk_net_modules_tbl, function(l) {
  pan_cancer_risk_net_modules_tbl4score <- cbind(l, mean_score = 0)
  pan_cancer_risk_net_modules_tbl4score$mean_score <- mean(pan_cancer_risk_net_gene_tbl$primitive_risk_score[which(pan_cancer_risk_net_gene_tbl$gene %in% l$gene)])
  pan_cancer_risk_net_modules_tbl4score
})

## Sort the mean scores of all modules
pan_cancer_risk_net_modules_flt_scores <- sapply(pan_cancer_risk_net_modules_tbl, function(j) {unique(j$mean_score)})
pan_cancer_risk_net_modules_flt_scores <- rev(sort(pan_cancer_risk_net_modules_flt_scores))

## Create module 1 for visualization

library(influential)
library(igraph)
library(visNetwork)

pan_cancer_risk_module1 <- subgraph(graph = pan_cancer_risk_net, 
                             vids = which(as_ids(V(pan_cancer_risk_net)) %in%  pan_cancer_risk_net_modules_tbl$`3843`$gene))

pan_cancer_risk_module1_vis <-
  visIgraph(igraph = pan_cancer_risk_module1) %>% 
  visOptions(
    # nodesIdSelection = TRUE, 
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '25px arial black',
    size = 25,
    # color = c("lightgrey", "orange"),
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "box" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = c("lightgrey"))

visExport(graph = pan_cancer_risk_module1_vis, 
          name = "pan_cancer_risk_module1_vis", 
          label = "")

#=============================================================================
#
#    Code chunk 9: Reconstruction and analysis of the Driver Net
#
#=============================================================================

# Remove other as the cancer type is not known and we will do a pan-cancer analysis any way!

somatic_snv_flt_cor$Other <- NULL
# somatic_cnv_cor$Other <- NULL
rnaseq_expr_cor$Other <- NULL

## Check if the cancer types of SNV, CNV and RNA-seq are identical
# identical(sort(names(somatic_snv_flt_cor)), sort(names(somatic_cnv_cor))) #TRUE
identical(sort(names(somatic_snv_flt_cor)), sort(names(rnaseq_expr_cor))) #TRUE

# Make the order of cancer types identical
# somatic_cnv_cor <- somatic_cnv_cor[match(names(somatic_snv_flt_cor), names(somatic_cnv_cor))]
rnaseq_expr_cor <- rnaseq_expr_cor[match(names(somatic_snv_flt_cor), names(rnaseq_expr_cor))]

# Prepare the tables for merging with other data
somatic_snv_flt_cor <- lapply(somatic_snv_flt_cor, function(i) {
  tmp_tbl <- i[,c(1,2)]
  colnames(tmp_tbl) <- c("from", "to")
  tmp_tbl$type <- "co-deletriousness"
  tmp_tbl
})

# somatic_cnv_cor <- lapply(somatic_cnv_cor, function(i) {
#   tmp_tbl <- i[,c(1,2)]
#   colnames(tmp_tbl) <- c("from", "to")
#   tmp_tbl$type <- "co-CNV"
#   tmp_tbl
# })

rnaseq_expr_cor <- lapply(rnaseq_expr_cor, function(i) {
  tmp_tbl <- i[,c(1,2)]
  colnames(tmp_tbl) <- c("from", "to")
  tmp_tbl$type <- "co-expression"
  tmp_tbl
})

# Merge the similarity tables
driver_marker_net <- somatic_snv_flt_cor
# for(i in 1:length(somatic_cnv_cor)) {
#   driver_marker_net[[i]] <- rbind(somatic_snv_flt_cor[[i]], somatic_cnv_cor[[i]])
# }
# names(driver_marker_net) <- names(somatic_cnv_cor)

# Add the co-expression data to the tables

for(i in 1:length(driver_marker_net)) {
  
  genes <- driver_marker_net[[i]][,c(1,2)] %>% 
    unlist() %>% 
    unname() %>% 
    unique()
  
  first_level_coex_index <- c(which(rnaseq_expr_cor[[i]]$from %in% genes), 
                              which(rnaseq_expr_cor[[i]]$to %in% genes)) %>% unique()
  
  first_level_coex_genes <- c(genes, 
                              rnaseq_expr_cor[[i]]$from[first_level_coex_index],
                              rnaseq_expr_cor[[i]]$to[first_level_coex_index]) %>% unique()
  
  first_two_levels_coex_index <- c(which(rnaseq_expr_cor[[i]]$from %in% first_level_coex_genes), 
                                   which(rnaseq_expr_cor[[i]]$to %in% first_level_coex_genes)) %>% unique()
  
  driver_marker_net[[i]] <- rbind(driver_marker_net[[i]], rnaseq_expr_cor[[i]][first_two_levels_coex_index,])
  
}

# Add PPI data to the driver/marker net

# Map gene names to stringDB ids for PPI analysis

## get the gene names in each risk net
driver_marker_net_genes <- lapply(driver_marker_net, function(i) {
  data.frame(gene_symbol = i[,c(1,2)] %>% 
               unlist() %>% 
               unname() %>% 
               unique())
})

## getSTRINGdb for human
string_db <- STRINGdb$new(species=9606, version = "11.5", 
                          score_threshold=900, # very high confidence (>0.9)
                          input_directory=""
)

driver_marker_net_genes <- lapply(driver_marker_net_genes, function(i) {
  
  string_db$map(my_data_frame = i, 
                my_data_frame_id_col_names = "gene_symbol", 
                takeFirst = TRUE, removeUnmappedRows = TRUE)
  
})

# Get the PPIs
driver_marker_net_ppi <- lapply(driver_marker_net_genes, function(i) {
  
  string_db$get_interactions(string_ids = i$STRING_id)
  
})

# Convert stringdb ids to gene symbols
driver_marker_net_ppi <- lapply(1:length(driver_marker_net_ppi), function(i) {
  tmp_row_index <- match(driver_marker_net_ppi[[i]]$from, driver_marker_net_genes[[i]]$STRING_id)
  driver_marker_net_ppi[[i]]$from <- driver_marker_net_genes[[i]]$gene_symbol[tmp_row_index]
  
  tmp_column_index <- match(driver_marker_net_ppi[[i]]$to, driver_marker_net_genes[[i]]$STRING_id)
  driver_marker_net_ppi[[i]]$to <- driver_marker_net_genes[[i]]$gene_symbol[tmp_column_index]
  
  tmp_tbl <- driver_marker_net_ppi[[i]][,c(1,2)]
  tmp_tbl$type <- "PPI"
  
  ## Remove duplicates
  tmp_tbl <- tmp_tbl[-which(duplicated(paste(tmp_tbl$from,
                                             tmp_tbl$to, 
                                             sep = "_"))),]
  
  return(tmp_tbl)
})

names(driver_marker_net_ppi) <- names(driver_marker_net_genes)

# Combine PPI with the original driver/marker net
driver_marker_net <- lapply(1:length(driver_marker_net), function(i) {
  rbind(driver_marker_net[[i]], driver_marker_net_ppi[[i]])
})

names(driver_marker_net) <- names(driver_marker_net_genes)

# Reconstruct the networks
driver_marker_net <- lapply(driver_marker_net, function(i) {
  igraph::graph_from_data_frame(i, directed=FALSE)
})

## Check the network layers
table(igraph::E(driver_marker_net$`CNS other`)$type)

###############################

# Calculate the IVI values of genes within each net
library(influential)
driver_marker_net_ivi <- lapply(driver_marker_net, ivi)

###############################

# Calculate the primitive gene driver scores
driver_marker_net_gene_tbl <- lapply(driver_marker_net_ivi, function(i) {
  data.frame(gene = names(i), ivi = i)
})

## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer

### First remove other
somatic_snv_flt_withCADD_summarized$Other <- NULL



somatic_snv_mean_deletriousness <- lapply(somatic_snv_flt_withCADD_summarized, function(i) {
  colMeans(i)
})

## Generate the gene table for each net and calculate the primitive scores
driver_marker_net_gene_tbl <- lapply(1:length(driver_marker_net_gene_tbl), function(i) {
  tmp_tbl <- driver_marker_net_gene_tbl[[i]]
  tmp_tbl$mean_deletriousness <- 0
  match_index <- which(names(somatic_snv_mean_deletriousness[[i]]) %in% rownames(tmp_tbl))
  tmp_tbl[names(somatic_snv_mean_deletriousness[[i]])[match_index], "mean_deletriousness"] <- 
    somatic_snv_mean_deletriousness[[i]][match_index]
  
  tmp_tbl$primitive_driver_score <- (tmp_tbl$ivi)*(tmp_tbl$mean_deletriousness) # node weights
  
  return(tmp_tbl)
  
})

names(driver_marker_net_gene_tbl) <- names(somatic_snv_mean_deletriousness)

# add gene score
driver_marker_net <- lapply(1:length(driver_marker_net), function(i) {
  set_vertex_attr(graph = driver_marker_net[[i]], name = "score", value = driver_marker_net_gene_tbl[[i]]$primitive_driver_score)
})

names(driver_marker_net) <- names(somatic_snv_mean_deletriousness)

# Convert the un-weighted to a node-weighted network
driver_marker_net <- lapply(1:length(driver_marker_net), function(i) {
  set_vertex_attr(graph = driver_marker_net[[i]], name = "weight", value = (driver_marker_net_gene_tbl[[i]]$ivi)^(driver_marker_net_gene_tbl[[i]]$mean_deletriousness))
})

names(driver_marker_net) <- names(somatic_snv_mean_deletriousness)

## Check the network node weights
summary(igraph::V(driver_marker_net$`CNS other`)$weight)

# Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
driver_marker_net_gene_tbl <- lapply(1:length(driver_marker_net), function(i) {
  tmp_tbl <- driver_marker_net_gene_tbl[[i]]
  tmp_tbl$mean_neighborhood_weight <- sapply(1:nrow(tmp_tbl), function(j) {
    igraph::neighbors(graph = driver_marker_net[[i]], 
                      v = V(driver_marker_net[[i]])[j])$score %>% mean()
  })
  
  tmp_tbl$mediator_score <- tmp_tbl$primitive_driver_score * tmp_tbl$mean_neighborhood_weight
  
  return(tmp_tbl)
})

names(driver_marker_net_gene_tbl) <- names(somatic_snv_mean_deletriousness)

########################

## Create a dataframe of top 20 drivers of all cancers
driver_net_top20_genes <- data.frame(Rank = c(1:20))

for(i in 1:length(driver_marker_net_gene_tbl)) {
  tmp_data = driver_marker_net_gene_tbl[[i]] %>% 
    slice_max(order_by = primitive_driver_score, n = 20) %>% 
    select(gene) %>% 
    rename(cancer_type = gene)
  
  driver_net_top20_genes <- cbind(driver_net_top20_genes, tmp_data)
}

colnames(driver_net_top20_genes)[2:ncol(driver_net_top20_genes)] <- names(driver_marker_net_gene_tbl)

#######################

# Top 20 drivers of HGG

driver_marker_net_gene_tbl$HGG %>% 
  slice_max(order_by = primitive_driver_score, n = 20) %>% 
  mutate(Rank = rank(-(primitive_driver_score))) %>% 
  select(Rank, gene) %>% 
  rename(Gene = gene) %>% 
  View()

#######################

# Top 20 drivers of Sarcoma

driver_marker_net_gene_tbl$Sarcoma %>% 
  slice_max(order_by = primitive_driver_score, n = 20) %>% 
  mutate(Rank = rank(-(primitive_driver_score))) %>% 
  select(Rank, gene) %>% 
  rename(Gene = gene) %>% 
  View()

####################

# Lymphoma top for poster

## Top 10 Risk of Lymhoma
risk_net_gene_tbl_Lymphoma <- 
risk_net_gene_tbl$Lymphoma %>% 
  slice_max(order_by = primitive_risk_score, n = 10) %>% 
  mutate(Rank = rank(-(primitive_risk_score))) %>% 
  select(Rank, gene) %>% 
  rename(Gene = gene)


## Top 10 Driver of Lymhoma
driver_marker_net_gene_tbl_Lymphoma <- 
  driver_marker_net_gene_tbl$Lymphoma %>% 
  slice_max(order_by = primitive_driver_score, n = 10) %>% 
  mutate(Rank = rank(-(primitive_driver_score))) %>% 
  select(Rank, gene) %>% 
  rename(Gene = gene)

net_gene_tbl_Lymphoma <- cbind(risk_net_gene_tbl_Lymphoma, driver_marker_net_gene_tbl_Lymphoma)
net_gene_tbl_Lymphoma <- net_gene_tbl_Lymphoma[,-3]
colnames(net_gene_tbl_Lymphoma)[c(2,3)] <- c("Risk", "Driver")

net_gene_tbl_Lymphoma_full_tbl <- ggtexttable(net_gene_tbl_Lymphoma,
                                                         rows = NULL, theme = ttheme("blank")) %>% 
  table_cell_bg(row = c(2:11), column = 1, fill = "darkolivegreen1") %>%
  table_cell_font(row = c(2:11), column = 1, face = "bold")

net_gene_tbl_Lymphoma_full_tbl <- net_gene_tbl_Lymphoma_full_tbl %>%
  tab_add_hline(at.row = c(1, 2), row.side = "top", linewidth = 3, linetype = 1) %>%
  tab_add_hline(at.row = c(nrow(net_gene_tbl_Lymphoma_full_tbl) + 1), row.side = "bottom", linewidth = 3, linetype = 1) %>%
  tab_add_vline(at.column = 2:tab_ncol(net_gene_tbl_Lymphoma_full_tbl), column.side = "left", from.row = 2, linetype = 2)

net_gene_tbl_Lymphoma_full_tbl






####################

# Detect communities/modules/clusters
driver_marker_net_modules <- lapply(driver_marker_net, function(i) {
  set.seed(3847)
  igraph::cluster_leiden(
    graph = i,
    objective_function = "CPM",
    weights = NULL,
    resolution_parameter = 0.5,
    beta = 0.05,
    initial_membership = NULL,
    n_iterations = 100,
    vertex_weights = igraph::V(i)$weight
  )
})

sizes(driver_marker_net_modules$`CNS other`) %>% as.vector() %>% summary()
length(driver_marker_net_modules$`CNS other`)

## Visualize the entire network and color by modules (NOT Recommended; It will take a while and will be too crowded)
plot(
  x = driver_marker_net_modules$`CNS other`, 
  y = driver_marker_net$`CNS other`,
  col = membership(driver_marker_net_modules$`CNS other`),
  vertex.label=NA,
  mark.groups = communities(driver_marker_net_modules$`CNS other`),
  edge.color = c("black", "red")[crossing(driver_marker_net_modules$`CNS other`, driver_marker_net$`CNS other`) + 1]
)

## Inspect and filter modules
driver_marker_net_modules_flt <- lapply(1:length(driver_marker_net_modules), function(i) {
  
  ### Create a table of modules and their genes
  tmp_tbl <- data.frame(
    module = membership(driver_marker_net_modules[[i]]) %>% as.integer(),
    gene = driver_marker_net_modules[[i]]$names
  )
  
  ### Add if the gene is mutated or not (could be a risk gene or not)
  tmp_tbl$mutated <- 0
  tmp_tbl$mutated[which(driver_marker_net_gene_tbl[[i]]$mean_deletriousness > 0)] <- 1
  
  ### Separate modules
  tmp_tbl <- lapply(unique(tmp_tbl$module), function(j){
    subset(tmp_tbl, module == j)
  })
  
  names(tmp_tbl) <- sapply(tmp_tbl, function(m) {
    unique(m$module)
  })
  
  ### Remove modules that have less than 2 Risk genes or their size is less than 4
  tmp_tbl <- lapply(tmp_tbl, function(k) {
    if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
      k
    } else {
      NULL
    }
  })
  
  ### Remove NULL modules
  tmp_tbl[which(sapply(tmp_tbl, is.null))] <- NULL
  
  ### Calculate the module scores
  tmp_tbl <- lapply(tmp_tbl, function(l) {
    tmp_tbl4score <- cbind(l, mean_score = 0)
    tmp_tbl4score$mean_score <- mean(driver_marker_net_gene_tbl[[i]]$primitive_driver_score[which(driver_marker_net_gene_tbl[[i]]$gene %in% l$gene)])
    tmp_tbl4score
  })
  
  tmp_tbl
})

names(driver_marker_net_modules_flt) <- names(driver_marker_net_modules)

## Sort the mean scores of all modules of each cancer type
driver_marker_net_modules_flt_scores <-
  lapply(driver_marker_net_modules_flt, function(i) {
    tmp_tbl <-
      sapply(i, function(j) {unique(j$mean_score)})
    tmp_tbl <- rev(sort(tmp_tbl))
    tmp_tbl
  })

### View the number of first ranked modules
sapply(driver_marker_net_modules_flt_scores, function(i) {i[1]}) %>% as.data.frame() %>% View()

## Visualize the modules
plot.igraph(subgraph(graph = driver_marker_net$`CNS other`, 
                     vids = which(as_ids(V(driver_marker_net$`CNS other`)) %in%  driver_marker_net_modules_flt$`CNS other`$`2`$gene)), 
            vertex.size =  3)

#####################

# Inspect the HGG results

## Create module 1 of HGG for visualization

library(influential)
library(igraph)
library(visNetwork)

driver_hgg_module1 <- subgraph(graph = driver_marker_net$HGG, 
                               vids = which(as_ids(V(driver_marker_net$HGG)) %in%  driver_marker_net_modules_flt$HGG$`1204`$gene))

driver_hgg_module1_vis <-
  visIgraph(igraph = driver_hgg_module1) %>% 
  visOptions(
    # nodesIdSelection = TRUE, 
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '40px arial black',
    size = 25,
    color = "orange",
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "box" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = c("lightgrey"))

visExport(graph = driver_hgg_module1_vis, 
          name = "driver_hgg_module1_vis", 
          label = "")

## It's also possible to create a subset of two modules and color the nodes based on mark.groups arg
### https://stackoverflow.com/questions/26913419/plot-communities-with-igraph

#####################

# Inspect the Sarcoma results

## Create module 1 of Sarcoma for visualization

driver_marker_net_modules_flt_scores$Sarcoma %>% head()

library(influential)
library(igraph)
library(visNetwork)

driver_Sarcoma_module1 <- subgraph(graph = driver_marker_net$Sarcoma, 
                                   vids = which(as_ids(V(driver_marker_net$Sarcoma)) %in%  driver_marker_net_modules_flt$Sarcoma$`3590`$gene))


driver_Sarcoma_module1_vis <-
  visIgraph(igraph = driver_Sarcoma_module1) %>% 
  visOptions(
    # nodesIdSelection = TRUE, 
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '40px arial black',
    size = 25,
    color = "orange",
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "box" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = c("lightgrey"))

visExport(graph = driver_Sarcoma_module1_vis, 
          name = "driver_Sarcoma_module1_vis", 
          label = "")

############################################
############################################

#### Inspect HGG drivers suggested by https://pecan.stjude.cloud/
driver_marker_net_gene_tbl$HGG[hgg_driver_genes,] %>% View()

driver_marker_net_gene_tbl$HGG %>% 
  slice_max(order_by = primitive_driver_score, n = 20) %>% 
  select(gene) %>% 
  slice(which(gene %in% hgg_driver_genes)) %>% 
  unlist()

driver_marker_net_modules_flt$HGG$`685` %>% 
  select(gene) %>% 
  filter(gene %in% hgg_driver_genes)

risk_net_modules_flt_scores$HGG %>% head()
risk_net_modules_flt$HGG$`2830` %>% 
  select(gene) %>% 
  filter(gene %in% hgg_driver_genes)

# Check the module TP53 is involved in
for (i in 1:length(driver_marker_net_modules_flt$HGG)) {
  module_tmp <-
    grep("^TP53$", driver_marker_net_modules_flt$HGG[[i]]$gene)
  if(length(module_tmp >= 1)) { 
    print(names(driver_marker_net_modules_flt$HGG[i]))
    print(driver_marker_net_modules_flt$HGG[[i]])
    cat("\n###########################################\n")
  }
}

head(driver_marker_net_modules_flt_scores$HGG, n=10)
driver_marker_net_modules_flt$HGG$`685`

View(driver_marker_net_gene_tbl$HGG)

View(somatic_snv_flt_withCADD_summarized$HGG[,"TP53", drop = F])

##############

#### Check the rank of modules with driver genes in HGG

hgg_driver_including_module_names <- vector(mode = "character")
for (i in 1:length(driver_marker_net_modules_flt$HGG)) {
  module_tmp <-
    grep(paste0(paste(paste("^", hgg_driver_genes, sep = ""), "$", sep = ""), collapse = "|"), driver_marker_net_modules_flt$HGG[[i]]$gene)
  if(length(module_tmp >= 1)) { 
    print(names(driver_marker_net_modules_flt$HGG[i]))
    hgg_driver_including_module_names <- append(hgg_driver_including_module_names,
                                                names(driver_marker_net_modules_flt$HGG[i]))
  }
}

hgg_driver_including_module_ranks <- which(names(driver_marker_net_modules_flt_scores$HGG) %in% hgg_driver_including_module_names)
hgg_driver_including_module_names <- names(driver_marker_net_modules_flt_scores$HGG)[hgg_driver_including_module_ranks]
hgg_driver_genes[which(hgg_driver_genes %in% driver_marker_net_modules_flt$HGG$`2411`$gene)]

############################################

# Inspect the Sarcoma results based on Sarcoma drivers suggested by https://pecan.stjude.cloud/

sarcoma_driver_genes[which(sarcoma_driver_genes %in% (driver_marker_net_gene_tbl$Sarcoma %>% slice_max(order_by = primitive_driver_score, n = 20) %>% select(gene) %>% unlist() %>% unname()))]
sarcoma_driver_genes[which(sarcoma_driver_genes %in% (risk_net_gene_tbl$Sarcoma %>% slice_max(order_by = primitive_risk_score, n = 20) %>% select(gene) %>% unlist() %>% unname()))]

driver_marker_net_modules_flt_scores$Sarcoma %>% head()
sarcoma_driver_genes[which(sarcoma_driver_genes %in% (driver_marker_net_modules_flt$Sarcoma$`152` %>% select(gene) %>% unlist() %>% unname()))]

risk_net_modules_flt_scores$Sarcoma %>% head()
risk_net_modules_flt$Sarcoma$`44` %>% 
  select(gene) %>% 
  filter(gene %in% sarcoma_driver_genes)

##########################

# Checking the first ranked module that includes co-deleteriousness data

head(driver_marker_net_modules_flt_scores$HGG, n=11)

driver_marker_net_modules_flt$HGG$`135`$gene

edge_type(g = driver_marker_net$HGG, v = driver_marker_net_modules_flt$HGG$`135`$gene)

##########################################################################################
##########################################################################################
##########################################################################################

# Evaluation of the Driver results.

## Save the gene data and modules as supplementary files

tmp_driver_marker_net_gene_tbl <- driver_marker_net_gene_tbl[-c(1,10)]
names(tmp_driver_marker_net_gene_tbl)[names(tmp_driver_marker_net_gene_tbl) == "Rhabdoid"] <- "Rhabdoid carcinoma"

for(i in 1:10) {
  colnames(tmp_driver_marker_net_gene_tbl[[i]]) <- c("Gene", "IVI", "Mean_deletriousness_score", "Driver_score", "Mean_neighborhood_weight", "Driverness_mediator_score")
}

library(openxlsx)
write.xlsx(x = tmp_driver_marker_net_gene_tbl, file = "Results/Evaluations/driver_net_gene_tbl.xlsx")

rm(tmp_driver_marker_net_gene_tbl)

############################

### save the modules

tmp_driver_marker_net_modules_flt <- driver_marker_net_modules_flt[-c(1,10)]
names(tmp_driver_marker_net_modules_flt)[names(tmp_driver_marker_net_modules_flt) == "Rhabdoid"] <- "Rhabdoid carcinoma"

tmp_driver_marker_net_modules_flt_scores <- driver_marker_net_modules_flt_scores[-c(1,10)]

tmp_driver_marker_net_modules_flt <- lapply(tmp_driver_marker_net_modules_flt, function(i) {
  
  tmp_data <- i
  tmp_data <- tmp_data[match(names(tmp_driver_marker_net_modules_flt_scores[[which(sapply(tmp_driver_marker_net_modules_flt_scores, length) == length(tmp_data))]]),
                             names(i))]
  
  tmp_data
  
})

for(i in 1:10) {
  tmp_driver_marker_net_modules_flt[[i]] <-
    lapply(1:length(tmp_driver_marker_net_modules_flt[[i]]), function(j) {
      tmp_data <- tmp_driver_marker_net_modules_flt[[i]][[j]]
      colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
      tmp_data[tmp_data == 1] <- "Yes"
      tmp_data[tmp_data == 0] <- "No"
      tmp_data$Module <- paste0("Module number ", j)
      tmp_data <- data.frame(Module = unique(tmp_data$Module),
                             Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                             Mean_score = unique(tmp_data$Mean_score))
      tmp_data
    })
  tmp_driver_marker_net_modules_flt[[i]] <- do.call(rbind, tmp_driver_marker_net_modules_flt[[i]])
}

library(openxlsx)
write.xlsx(x = tmp_driver_marker_net_modules_flt, file = "Results/Evaluations/driver_net_modules_flt.xlsx")

rm(tmp_driver_marker_net_modules_flt, tmp_driver_marker_net_modules_flt_scores)

##########################################################################################

# Create a dataframe of 1st-ranked driver modules of all cancer types

driver_net_1st_modules <- data.frame(Number = vector("integer"))

for(i in 1:length(driver_marker_net_modules_flt)) {
  tmp_module_genes = driver_marker_net_modules_flt[[i]][names(driver_marker_net_modules_flt_scores[[i]])[1]][[1]]$gene
  driver_net_1st_modules <- merge(driver_net_1st_modules, 
                                data.frame(Number = c(1:length(tmp_module_genes)),
                                           "cancer_type" = tmp_module_genes), 
                                by = "Number", all = TRUE)
}

driver_net_1st_modules$Number <- NULL

colnames(driver_net_1st_modules) <- names(driver_marker_net_modules_flt)

#####################

## Evaluation of the top 20 genes

### DMG
hgg_driver_genes[
  which(
    hgg_driver_genes 
    %in% 
      (driver_marker_net_gene_tbl$DMG %>% 
         slice_max(order_by = primitive_driver_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### HGG
hgg_driver_genes[
  which(
    hgg_driver_genes 
    %in% 
      (driver_marker_net_gene_tbl$HGG %>% 
         slice_max(order_by = primitive_driver_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### Leukaemia
leukaemia_driver_genes[
  which(
    leukaemia_driver_genes 
    %in% 
      (driver_marker_net_gene_tbl$Leukaemia %>% 
         slice_max(order_by = primitive_driver_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### NBL
nbl_driver_genes[
  which(
    nbl_driver_genes 
    %in% 
      (driver_marker_net_gene_tbl$NBL %>% 
         slice_max(order_by = primitive_driver_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### Sarcoma
sarcoma_driver_genes[
  which(
    sarcoma_driver_genes 
    %in% 
      (driver_marker_net_gene_tbl$Sarcoma %>% 
         slice_max(order_by = primitive_driver_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### CNS embryonal
embryonal_cns_driver_genes[
  which(
    embryonal_cns_driver_genes 
    %in% 
      (driver_marker_net_gene_tbl$`CNS embryonal` %>% 
         slice_max(order_by = primitive_driver_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### EPD
epd_driver_genes[
  which(
    epd_driver_genes 
    %in% 
      (driver_marker_net_gene_tbl$EPD %>% 
         slice_max(order_by = primitive_driver_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### Lymphoma
lymphoma_driver_genes[
  which(
    lymphoma_driver_genes 
    %in% 
      (driver_marker_net_gene_tbl$Lymphoma %>% 
         slice_max(order_by = primitive_driver_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### MB (Medulloblastoma)
mb_driver_genes[
  which(
    mb_driver_genes 
    %in% 
      (driver_marker_net_gene_tbl$MB %>% 
         slice_max(order_by = primitive_driver_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

#####################

### Rhabdoid
rhabdoid_driver_genes[
  which(
    rhabdoid_driver_genes 
    %in% 
      (driver_marker_net_gene_tbl$Rhabdoid %>% 
         slice_max(order_by = primitive_driver_score, n = 20) %>% 
         select(gene) %>% 
         unname() %>% 
         unlist())
  )
] %>% cat()

##########################################################################################

## Evaluation of the 1st-ranked modules of all cancer types

### DMG
hgg_driver_genes[
  which(
    hgg_driver_genes 
    %in% 
      na.omit(driver_net_1st_modules$DMG)
  )
] %>% cat()

#####################

### HGG
hgg_driver_genes[
  which(
    hgg_driver_genes 
    %in% 
      na.omit(driver_net_1st_modules$HGG)
  )
] %>% cat()

#####################

### Leukaemia
leukaemia_driver_genes[
  which(
    leukaemia_driver_genes 
    %in% 
      na.omit(driver_net_1st_modules$Leukaemia)
  )
] %>% cat()

#####################

### NBL
nbl_driver_genes[
  which(
    nbl_driver_genes 
    %in% 
      na.omit(driver_net_1st_modules$NBL)
  )
] %>% cat()

#####################

### Sarcoma
sarcoma_driver_genes[
  which(
    sarcoma_driver_genes 
    %in% 
      na.omit(driver_net_1st_modules$Sarcoma)
  )
] %>% cat()

#####################

### CNS embryonal
embryonal_cns_driver_genes[
  which(
    embryonal_cns_driver_genes 
    %in% 
      na.omit(driver_net_1st_modules$`CNS embryonal`)
  )
] %>% cat()

#####################

### EPD
epd_driver_genes[
  which(
    epd_driver_genes 
    %in% 
      na.omit(driver_net_1st_modules$EPD)
  )
] %>% cat()

#####################

### Lymphoma
lymphoma_driver_genes[
  which(
    lymphoma_driver_genes 
    %in% 
      na.omit(driver_net_1st_modules$Lymphoma)
  )
] %>% cat()

#####################

### MB
mb_driver_genes[
  which(
    mb_driver_genes 
    %in% 
      na.omit(driver_net_1st_modules$MB)
  )
] %>% cat()

#####################

### Rhabdoid
rhabdoid_driver_genes[
  which(
    rhabdoid_driver_genes 
    %in% 
      na.omit(driver_net_1st_modules$Rhabdoid)
  )
] %>% cat()

##########################################################################################

## Visualization of the results
library(fishualize)
library(patchwork)

pecan_prism_common <- readr::read_delim(file = "Results/Top Childhood Evaluation.txt", delim = "\t")
pecan_prism_common[pecan_prism_common == "None"] <- ""

pecan_prism_common_top20_plot <- ggplot(data = pecan_prism_common, aes(x = `Cancer type`, 
                                                                       y = `Number of Genes in PeCan/PediCan (out of top 20)`, 
                                                                       fill = Context)) +
  geom_bar(position = "dodge2", stat = 'identity') +
  geom_text(aes(x = `Cancer type`, 
                y = ifelse(`Number of Genes in PeCan/PediCan (out of top 20)` > 0, `Number of Genes in PeCan/PediCan (out of top 20)` - `Number of Genes in PeCan/PediCan (out of top 20)`/2, 0.5), 
                label = `Genes in PeCan/PediCan (out of top 20)`), position=position_dodge(width=0.9), size = 1.6, show.legend = FALSE) +
  labs(y= '# Genes common with PeCan/PediCan drivers',
       title = "Top 20 Genes") +
  coord_flip() +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_fill_fish_d(begin = 0.27, end = 0.51) +
  # scale_fill_viridis_d(begin = 0.7) +
  theme_bw() +
  theme(text = element_text(size = 7),
        legend.position = "right", 
        legend.box.background = element_rect(colour = "black"))

pecan_prism_common_top20_plot


#########################

pecan_prism_common_1st_module_plot <- ggplot(data = pecan_prism_common, aes(x = `Cancer type`, 
                                                                       y = `Number of Genes in PeCan/PediCan (out of module 1 genes)`, 
                                                                       fill = Context)) +
  geom_bar(position = "dodge2", stat = 'identity') +
  geom_text(aes(x = `Cancer type`, 
                y = ifelse(`Number of Genes in PeCan/PediCan (out of module 1 genes)` > 0, `Number of Genes in PeCan/PediCan (out of module 1 genes)` - `Number of Genes in PeCan/PediCan (out of module 1 genes)`/2, 0.5), 
                label = `Genes in PeCan/PediCan (out of module 1 genes)`), position=position_dodge(width=0.9), size = 1.6, show.legend = FALSE) +
  labs(y= '# Genes common with PeCan/PediCan drivers', 
       title = "1st-Ranked Module Genes") +
  coord_flip() +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_fill_fish_d(begin = 0.27, end = 0.51) +
  # scale_fill_viridis_d(begin = 0.7) +
  theme_bw() +
  theme(text = element_text(size = 7),
        legend.position = "right", 
        legend.box.background = element_rect(colour = "black"))

pecan_prism_common_1st_module_plot

#########################

pecan_prism_common_plot <- (pecan_prism_common_top20_plot + theme(legend.position = "none")) + 
  (pecan_prism_common_1st_module_plot + xlab(NULL)) +
  plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = 'bold'))

pecan_prism_common_plot

ggsave(filename = "Results/Evaluations/pecan_prism_common_plot.pdf", 
       plot = pecan_prism_common_plot, 
       device = "pdf", width = 19, height = 7, units = "cm")

##########################################################################################
##########################################################################################

## Literature mining of cancer driver results
library(pubmedR)

lit_driver_data <- tibble(`Cancer type` = vector("character"), 
                              `Gene group` = vector("character"),
                              Gene = vector("character"),
                              Num_PMIDs = vector("integer")
                              )

for(i in 1:length(driver_marker_net_gene_tbl[-c(1,10)])) {
  
  if(names(driver_marker_net_gene_tbl[-c(1,10)])[i] == "HGG") {
    main_search_text <- "((HGG[Title/Abstract]) OR (HGGs[Title/Abstract]) OR (High-Grade Glioma[Title/Abstract]) OR (High-Grade Gliomas[Title/Abstract]))"
    
  } else if(names(driver_marker_net_gene_tbl[-c(1,10)])[i] == "DMG") {
    main_search_text <- "((DMG[Title/Abstract]) OR (DMGs[Title/Abstract]) OR (Diffuse Midline Glioma[Title/Abstract]) OR (Diffuse Midline Gliomas[Title/Abstract]))"
    
  } else if(names(driver_marker_net_gene_tbl[-c(1,10)])[i] == "Leukaemia") {
    main_search_text <- "((leukemia[Title/Abstract]) OR (leukemias[Title/Abstract]) OR (Leukaemias[Title/Abstract]) OR (Leukaemia[Title/Abstract]))"
    
  } else if(names(driver_marker_net_gene_tbl[-c(1,10)])[i] == "Lymphoma") {
    main_search_text <- "((Lymphoma[Title/Abstract]) OR (Lymphomas[Title/Abstract]))"
    
  } else if(names(driver_marker_net_gene_tbl[-c(1,10)])[i] == "NBL") {
    main_search_text <- "((NBL[Title/Abstract]) OR (NBLs[Title/Abstract]) OR (Neuroblastoma[Title/Abstract]) OR (Neuroblastomas[Title/Abstract]))"
    
  } else if(names(driver_marker_net_gene_tbl[-c(1,10)])[i] == "Sarcoma") {
    main_search_text <- "((Sarcoma[Title/Abstract]) OR (Sarcomas[Title/Abstract]))"
    
  } else if(names(driver_marker_net_gene_tbl[-c(1,10)])[i] == "CNS embryonal") {
    main_search_text <- "((atypical teratoid/rhabdoid tumor[Title/Abstract]) OR (AT/RT[Title/Abstract]) OR (embryonal tumors with multilayered rosettes[Title/Abstract]) OR (ETMR[Title/Abstract]) OR (ependymoblastoma[Title/Abstract]))"
    
  } else if(names(driver_marker_net_gene_tbl[-c(1,10)])[i] == "EPD") {
    main_search_text <- "((EPD[Title/Abstract]) OR (Ependymoma[Title/Abstract]) OR (Ependymomas[Title/Abstract]) OR (Subependymoma[Title/Abstract]))"
    
  } else if(names(driver_marker_net_gene_tbl[-c(1,10)])[i] == "MB") {
    main_search_text <- "((MB[Title/Abstract]) OR (medulloblastoma[Title/Abstract]) OR (medulloblastomas[Title/Abstract]))"
    
  } else if(names(driver_marker_net_gene_tbl[-c(1,10)])[i] == "Rhabdoid") {
    main_search_text <- "((Rhabdoid[Title/Abstract]) OR (Renal medullary carcinoma[Title/Abstract]) OR (RMC[Title/Abstract]))"
  }
  
  ### Add top 20 genes
  tmp_top20_tbl <- tibble(`Cancer type` = names(driver_marker_net_gene_tbl[-c(1,10)])[i],
                          `Gene group` = "Top 20 genes",
                          Gene = (driver_marker_net_gene_tbl[-c(1,10)][[i]] %>%
                                    slice_max(order_by = primitive_driver_score, n = 20) %>%
                                    select(gene) %>%
                                    unname() %>%
                                    unlist()))
  
  tmp_top20_pmid_list <- lapply(1:20, function(j) {
    query_text <- paste0("(", tmp_top20_tbl$Gene[j], "[Title/Abstract]) AND ", main_search_text)
    tmp_entrez <- pmQueryTotalCount(query = query_text, api_key = "1a1d2b704b3bad87d1ded4f2d95d611e1909")
    tmp_entrez_num <- tmp_entrez$total_count
    print(paste0("Top gene ", j, " is done!"))
    return(tmp_entrez_num)
  })
  
  tmp_top20_tbl$Num_PMIDs <- unlist(tmp_top20_pmid_list)
  
  lit_driver_data <- rbind(lit_driver_data, tmp_top20_tbl)
  
  ### Add the 1st ranked module
  tmp_1st_module_tbl <- tibble(`Cancer type` = names(driver_marker_net_gene_tbl[-c(1,10)])[i],
                          `Gene group` = "1st-ranked module",
                          Gene = (na.omit(driver_net_1st_modules[,-c(1,10)][,i])))
  
  tmp_1st_module_pmid_list <- lapply(1:length(na.omit(driver_net_1st_modules[,-c(1,10)][,i])), function(m) {
    query_text <- paste0("(", tmp_1st_module_tbl$Gene[m], "[Title/Abstract]) AND ", main_search_text)
    tmp_entrez <- pmQueryTotalCount(query = query_text, api_key = "1a1d2b704b3bad87d1ded4f2d95d611e1909")
    tmp_entrez_num <- tmp_entrez$total_count
    print(paste0("Module 1 gene ", m, " is done!"))
    return(tmp_entrez_num)
  })
  
  tmp_1st_module_tbl$Num_PMIDs <- unlist(tmp_1st_module_pmid_list)
  
  lit_driver_data <- rbind(lit_driver_data, tmp_1st_module_tbl)
  
  print(paste0(names(driver_marker_net_gene_tbl[-c(1,10)])[i], " is done!"))
  
}

#################################

## Visualization of the results

lit_driver_data[lit_driver_data == "Rhabdoid"] <- "Rhabdoid carcinoma"

driver_lit_mining_plot <- ggplot(data = lit_driver_data, aes(x = `Cancer type`, 
                                                                y = log2(Num_PMIDs+1), 
                                                                fill = `Gene group`)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.01, width=0.58, size = 0.2) +
  # geom_boxplot(position = position_dodge(0.8), width=0.4, 
  #              aes(color = `Gene group`),
  #              fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0, show.legend = F) +
  labs(y= 'Log2(number of papers + 1)',
       title = "Literature-based Association of Top\nDrivers with Their Respective Cancers") +
  coord_flip() +
  scale_y_continuous(breaks = scales::pretty_breaks()) +
  scale_fill_fish_d(option = "Elagatis_bipinnulata") +
  # scale_color_fish_d(option = "Elagatis_bipinnulata") +
  # scale_fill_viridis_d(begin = 0.32, end = 0.83) +
  theme_bw() +
  theme(text = element_text(size = 7),
        legend.position = "right", 
        legend.box.background = element_rect(colour = "black"))

driver_lit_mining_plot

#########################

lit_mining_plot <- (driver_lit_mining_plot + theme(legend.position = "none")) + 
  (risk_lit_mining_plot + xlab(NULL)) +
  plot_annotation(tag_levels = list(c("C", "D"))) &
  theme(plot.tag = element_text(face = 'bold'))

lit_mining_plot

ggsave(filename = "Results/Evaluations/lit_mining_plot.pdf", 
       plot = lit_mining_plot, 
       device = "pdf", width = 19, height = 7, units = "cm")

### Save the tables
lit_mining_data <- lit_driver_data
lit_mining_data$Context <- "Cancer Driver"

lit_mining_data <- rbind(lit_mining_data, cbind(lit_risk_data, Context = "Cancer Risk"))

write.table(x = lit_mining_data, file = "Results/Evaluations/lit_mining_data.csv", append = F, quote = F, sep = ",", row.names = F)

##########################################################################################
##########################################################################################

## ORA on the first ranked driver modules

library(enrichR)

### Using enrichR

#Get the list of all databases
enrichR_dbs <- enrichR::listEnrichrDbs()
enrichR_dbs %>% as.data.frame() %>% view()

#Select you desired databases
dbs.GO.BP <- "GO_Biological_Process_2023"
# dbs.KEGG.pathway <- "KEGG_2021_Human"

#ORA on GO-BP terms for 1st driver module of each cancer type

drivers_1st_module_ORA <- lapply(1:10, function(i) {
  enrichR::enrichr(genes = as.character(na.omit(driver_net_1st_modules[,-c(1,10)][,i])),
                   databases = dbs.GO.BP)[[1]]
})

names(drivers_1st_module_ORA) <- colnames(driver_net_1st_modules[,-c(1,10)])

###########################

### save the results
library(openxlsx)

write.xlsx(x = drivers_1st_module_ORA, file = "Results/Evaluations/ORA/Drivers/drivers_1st_module_ORA.xlsx")

###########################

### visualize the results
drivers_1st_module_ORA_vis_list <- lapply(1:10, function(i) {
  
  #Filter the data based on their Combined Score
  tmp_ora <- drivers_1st_module_ORA[[i]] %>% 
    slice_max(Combined.Score, n = 5) 
  
  #Correcting the levels
  tmp_ora$Term <- factor(tmp_ora$Term, levels = c(as.character(tmp_ora$Term)[order(tmp_ora$Combined.Score)]))
  
  #Visualization of terms
  tmp_ora_plot <- ggplot(tmp_ora, aes(x = Term, y = Combined.Score,
                                      color = Adjusted.P.value)) +  
    geom_point(size = 2) +
    ylab("Combined Score") + xlab(NULL) + 
    scale_color_viridis_c(name="Adjusted\nP-value", breaks= scales::pretty_breaks()) + 
    coord_flip() +
    labs(title = expression("ORA of the"~"1"^"st"~ "Ranked Driver Module of"),
         subtitle = colnames(driver_net_1st_modules[,-c(1,10)])[i]) +
    scale_y_continuous(breaks= scales::pretty_breaks(), labels = function(x) format(x, scientific = TRUE, digits = 1)) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 35)) +
    theme_classic() + theme(
      legend.position = "right",
      # legend.position = c(.95, 0.02),
      # legend.justification = c("right", "bottom"),
      legend.box.background = element_rect(color="black", size=1),
      text = element_text(size = 7))
  
  #Saving the plot
  ggsave(filename = paste0("Results/Evaluations/ORA/drivers/", colnames(driver_net_1st_modules[,-c(1,10)])[i], "_drivers_1st_module_ORA.pdf"), 
         plot = tmp_ora_plot,
         device = "pdf",width = 9.7, height = 4.6, units = "cm")
  
  tmp_ora_plot
  
})

drivers_1st_module_ORA_vis_list[[4]]

###########################

### combine the driver and risk results of lymphoma

#Filter the data based on their Combined Score
combined_1st_module_ORA_lymphoma <- rbind(
  cbind(drivers_1st_module_ORA[[4]] %>% slice_max(Combined.Score, n = 5), Context= "Cancer Driver"),
  cbind(risks_1st_module_ORA[[4]] %>% slice_max(Combined.Score, n = 5) , Context= "Cancer Risk")
) 

#Correcting the levels
combined_1st_module_ORA_lymphoma$Term <- factor(combined_1st_module_ORA_lymphoma$Term, 
                                               levels = c(as.character(combined_1st_module_ORA_lymphoma$Term)[order(combined_1st_module_ORA_lymphoma$Combined.Score)]))

#Visualization of terms
combined_1st_module_ORA_lymphoma_plot <- ggplot(combined_1st_module_ORA_lymphoma, 
                                               aes(x = Term, y = Combined.Score, color = Adjusted.P.value)) +  
  geom_point(size = 2) +
  facet_wrap(~ Context, scales = "free") +
  ylab("Combined Score") + xlab(NULL) + 
  scale_color_viridis_c(name="Adjusted\nP-value", breaks= scales::pretty_breaks()) + 
  coord_flip() +
  labs(title = expression("GO-BP-based ORA of the"~"1"^"st"~ "Ranked Modules of Lymphoma")) +
  scale_y_continuous(breaks= scales::pretty_breaks(), labels = function(x) format(x, scientific = TRUE, digits = 1)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 35)) +
  theme_classic() + theme(
    legend.position = "right",
    # legend.position = c(.95, 0.02),
    # legend.justification = c("right", "bottom"),
    legend.box.background = element_rect(color="black", size=1),
    text = element_text(size = 7))

combined_1st_module_ORA_lymphoma_plot

#Saving the plot
ggsave(filename = paste0("Results/Evaluations/ORA/lymphoma_1st_modules_ORA.pdf"), 
       plot = combined_1st_module_ORA_lymphoma_plot,
       device = "pdf",width = 20, height = 5.5, units = "cm")

##########################################################################################
##########################################################################################

## Evaluation of the first ranked mediators

### Evaluate the first ranked mediators and their first-order connected modules in driver networks

first_ranked_driver_mediator_tbl <- tibble("Cancer type" = vector(mode = "character", length = 10), 
                                         "1st Ranked Driverness mediator" = vector(mode = "character", length = 10), 
                                         "1st order connections" = vector(mode = "character", length = 10), 
                                         "Commonalities with 1st module" = vector(mode = "character", length = 10), 
                                         "Commonalities with 2nd module" = vector(mode = "character", length = 10), 
                                         "Commonalities with 3rd module" = vector(mode = "character", length = 10), 
                                         "Commonalities with 4th module" = vector(mode = "character", length = 10), 
                                         "Commonalities with 5th module" = vector(mode = "character", length = 10))

first_ranked_driver_mediator_list <- list()

for(j in 1:10) {
  first_ranked_driver_mediator_tbl[j,1] <- names(driver_marker_net_gene_tbl[-c(1,10)])[j]
  
  tmp_first_ranked_mediator <- driver_marker_net_gene_tbl[-c(1,10)][[j]] %>%
    slice_max(mediator_score, n = 1) %>%
    select(gene) %>%
    unlist() %>%
    unname()
  first_ranked_driver_mediator_tbl[j,2] <- tmp_first_ranked_mediator

  tmp_first_ranked_mediator_Neighbors <- igraph::neighbors(graph = driver_marker_net[-c(1,10)][[j]], v = tmp_first_ranked_mediator, mode = "all") %>% as_ids()
  first_ranked_driver_mediator_tbl[j,3] <- paste0(tmp_first_ranked_mediator_Neighbors, collapse = ", ")

  tmp_first_ranked_mediator_moduleNeighbors <- vector(mode = "character")

  for(i in 1:5) {
    tmp_module <- driver_marker_net_modules_flt[-c(1,10)][[j]][[names(driver_marker_net_modules_flt_scores[-c(1,10)][[j]])[i]]]
    first_orders_tmp = tmp_module$gene[which(tmp_module$gene %in% tmp_first_ranked_mediator_Neighbors)]
    first_ranked_driver_mediator_tbl[j,i+3] <- paste0(first_orders_tmp, collapse = ", ")

    if(length(first_orders_tmp) > 0) {
      tmp_first_ranked_mediator_moduleNeighbors <- append(tmp_first_ranked_mediator_moduleNeighbors, first_orders_tmp)
    }
  }
  first_ranked_driver_mediator_list <- append(first_ranked_driver_mediator_list, list(tmp_first_ranked_mediator_moduleNeighbors))
}

first_ranked_driver_mediator_tbl[first_ranked_driver_mediator_tbl == "Rhabdoid"] <- "Rhabdoid carcinoma"
readr::write_delim(x = first_ranked_driver_mediator_tbl,
                   file = "Results/Evaluations/Mediators/first_ranked_driver_mediator_tbl.csv", 
                   delim = ",")

names(first_ranked_driver_mediator_list) <- names(driver_marker_net_gene_tbl[-c(1,10)])

###########

#### Visualization of results of DMG
DMG_first_ranked_driver_mediator_neighbors_net <- subgraph(graph = driver_marker_net$DMG, 
                                                           vids = which(as_ids(V(driver_marker_net$DMG)) %in%  c(first_ranked_driver_mediator_tbl[first_ranked_driver_mediator_tbl[,1] == "DMG",2], 
                                                                                                                 first_ranked_driver_mediator_list$DMG)))

DMG_first_ranked_driver_mediator_neighbors_net_tbl <- igraph::as_data_frame(DMG_first_ranked_driver_mediator_neighbors_net)
readr::write_delim(x = DMG_first_ranked_driver_mediator_neighbors_net_tbl,
                   file = "Results/Evaluations/Mediators/DMG_driver_mediator_neighbors_net_tbl.txt", 
                   delim = "\t")

############################

#### Visualization of results of Sarcoma

library(influential)
library(igraph)
library(visNetwork)

sarcoma_first_ranked_driver_mediator_neighbors_net <- subgraph(graph = driver_marker_net$Sarcoma, 
                                                  vids = which(as_ids(V(driver_marker_net$Sarcoma)) %in%  c(first_ranked_driver_mediator_tbl[first_ranked_driver_mediator_tbl[,1] == "Sarcoma",2], 
                                                                                                            first_ranked_driver_mediator_list$Sarcoma)))

sarcoma_first_ranked_driver_mediator_neighbors_net_tbl <- igraph::as_data_frame(sarcoma_first_ranked_driver_mediator_neighbors_net)
readr::write_delim(x = sarcoma_first_ranked_driver_mediator_neighbors_net_tbl,
                   file = "Results/Evaluations/Mediators/sarcoma_driver_mediator_neighbors_net_tbl.txt", 
                   delim = "\t")

sarcoma_first_ranked_driver_mediator_neighbors_net_vis <-
  visIgraph(igraph = sarcoma_first_ranked_driver_mediator_neighbors_net) %>% 
  visOptions(
    nodesIdSelection = TRUE,
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '25px arial black',
    size = 25,
    # color = c("lightgrey", "orange"),
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "circle" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = c("lightgrey"))

visExport(graph = sarcoma_first_ranked_driver_mediator_neighbors_net_vis, 
          name = "sarcoma_first_ranked_driver_mediator_neighbors_net_vis", 
          label = "")

##########################################

## Evaluate the first ranked mediators and their first-order connected modules in risk networks

first_ranked_risk_mediator_tbl <- tibble("Cancer type" = vector(mode = "character", length = 10), 
                                         "1st Ranked Risk mediator" = vector(mode = "character", length = 10), 
                                         "1st order connections" = vector(mode = "character", length = 10), 
                                         "Commonalities with 1st module" = vector(mode = "character", length = 10), 
                                         "Commonalities with 2nd module" = vector(mode = "character", length = 10), 
                                         "Commonalities with 3rd module" = vector(mode = "character", length = 10), 
                                         "Commonalities with 4th module" = vector(mode = "character", length = 10), 
                                         "Commonalities with 5th module" = vector(mode = "character", length = 10))

first_ranked_risk_mediator_list <- list()

for(j in 1:10) {
  first_ranked_risk_mediator_tbl[j,1] <- names(risk_net_gene_tbl[-c(1,9)])[j]
  
  tmp_first_ranked_mediator <- risk_net_gene_tbl[-c(1,9)][[j]] %>%
    slice_max(risk_mediator_score, n = 1) %>%
    select(gene) %>%
    unlist() %>%
    unname()
  first_ranked_risk_mediator_tbl[j,2] <- tmp_first_ranked_mediator
  
  tmp_first_ranked_mediator_Neighbors <- igraph::neighbors(graph = risk_net[-c(1,9)][[j]], v = tmp_first_ranked_mediator, mode = "all") %>% as_ids()
  first_ranked_risk_mediator_tbl[j,3] <- paste0(tmp_first_ranked_mediator_Neighbors, collapse = ", ")
  
  tmp_first_ranked_mediator_moduleNeighbors <- vector(mode = "character")
  
  for(i in 1:5) {
    tmp_module <- risk_net_modules_flt[-c(1,9)][[j]][[names(risk_net_modules_flt_scores[-c(1,9)][[j]])[i]]]
    first_orders_tmp = tmp_module$gene[which(tmp_module$gene %in% tmp_first_ranked_mediator_Neighbors)]
    first_ranked_risk_mediator_tbl[j,i+3] <- paste0(first_orders_tmp, collapse = ", ")
    
    if(length(first_orders_tmp) > 0) {
      tmp_first_ranked_mediator_moduleNeighbors <- append(tmp_first_ranked_mediator_moduleNeighbors, first_orders_tmp)
    }
  }
  first_ranked_risk_mediator_list <- append(first_ranked_risk_mediator_list, list(tmp_first_ranked_mediator_moduleNeighbors))
}

first_ranked_risk_mediator_tbl[first_ranked_risk_mediator_tbl == "Rhabdoid"] <- "Rhabdoid carcinoma"
readr::write_delim(x = first_ranked_risk_mediator_tbl,
                   file = "Results/Evaluations/Mediators/first_ranked_risk_mediator_tbl.csv", 
                   delim = ",")

names(first_ranked_risk_mediator_list) <- names(risk_net_gene_tbl[-c(1,9)])

###########

#### Visualization of results of DMG

DMG_first_ranked_risk_mediator_neighbors_net <- subgraph(graph = risk_net$DMG, 
                                                         vids = which(as_ids(V(risk_net$DMG)) %in%  c(first_ranked_risk_mediator_tbl[first_ranked_risk_mediator_tbl[,1] == "DMG",2], 
                                                                                                      first_ranked_risk_mediator_list$DMG)))

DMG_first_ranked_risk_mediator_neighbors_net_tbl <- igraph::as_data_frame(DMG_first_ranked_risk_mediator_neighbors_net)
readr::write_delim(x = DMG_first_ranked_risk_mediator_neighbors_net_tbl,
                   file = "Results/Evaluations/Mediators/DMG_risk_mediator_neighbors_net_tbl.txt", 
                   delim = "\t")

##########################

#### Visualization of results of Sarcoma

sarcoma_first_ranked_risk_mediator_neighbors_net <- subgraph(graph = risk_net$Sarcoma, 
                                                             vids = which(as_ids(V(risk_net$Sarcoma)) %in%  c(first_ranked_risk_mediator_tbl[first_ranked_risk_mediator_tbl[,1] == "Sarcoma",2], 
                                                                                                              first_ranked_risk_mediator_list$Sarcoma)))

sarcoma_first_ranked_risk_mediator_neighbors_net_tbl <- igraph::as_data_frame(sarcoma_first_ranked_risk_mediator_neighbors_net)
readr::write_delim(x = sarcoma_first_ranked_risk_mediator_neighbors_net_tbl,
                   file = "Results/Evaluations/Mediators/sarcoma_risk_mediator_neighbors_net_tbl.txt", 
                   delim = "\t")

sarcoma_first_ranked_risk_mediator_neighbors_net_vis <-
  visIgraph(igraph = sarcoma_first_ranked_risk_mediator_neighbors_net) %>% 
  visOptions(
    nodesIdSelection = TRUE,
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '25px arial black',
    size = 25,
    # color = c("lightgrey", "orange"),
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "circle" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = c("lightgrey"))

visExport(graph = sarcoma_first_ranked_risk_mediator_neighbors_net_vis, 
          name = "sarcoma_first_ranked_risk_mediator_neighbors_net_vis", 
          label = "")

#=============================================================================
#
#    Code chunk 10: Reconstruction and analysis of the pan-cancer Driver Net
#
#=============================================================================

pan_cancer_driver_net <- pan_somatic_snv_cor

# Prepare the table for merging with other data
pan_cancer_driver_net <- pan_cancer_driver_net[,c(1,2)]
colnames(pan_cancer_driver_net) <- c("from", "to")
pan_cancer_driver_net$type <- "co-deletriousness"

# Add the co-expression data to the table
pan_cancer_driver_net_genes <- pan_cancer_driver_net[,c(1,2)] %>% 
  unlist() %>% 
  unname() %>% 
  unique()

pan_cancer_driver_net_first_level_coex_index <- c(which(pan_somatic_rnaseq_cor$row %in% pan_cancer_driver_net_genes), 
                                                  which(pan_somatic_rnaseq_cor$column %in% pan_cancer_driver_net_genes)) %>% unique()


pan_cancer_driver_net_first_level_coex_genes <- c(pan_cancer_driver_net_genes, 
                                                  pan_somatic_rnaseq_cor$row[pan_cancer_driver_net_first_level_coex_index],
                                                  pan_somatic_rnaseq_cor$column[pan_cancer_driver_net_first_level_coex_index]) %>% unique()

pan_cancer_driver_net_first_two_level_coex_index <- c(which(pan_somatic_rnaseq_cor$row %in% pan_cancer_driver_net_first_level_coex_genes), 
                                                      which(pan_somatic_rnaseq_cor$column %in% pan_cancer_driver_net_first_level_coex_genes)) %>% unique()

pan_somatic_rnaseq_cor_flt <- pan_somatic_rnaseq_cor[pan_cancer_driver_net_first_two_level_coex_index,c(1,2)]
colnames(pan_somatic_rnaseq_cor_flt) <- c("from", "to")
pan_somatic_rnaseq_cor_flt$type <- "co-expression"

pan_cancer_driver_net <- rbind(pan_cancer_driver_net, pan_somatic_rnaseq_cor_flt)

rm(pan_somatic_rnaseq_cor_flt)

# Add PPI data to the driver net
library(STRINGdb)

# Map gene names to stringDB ids for PPI analysis

## get the gene names in each driver net
pan_cancer_driver_net_genes_plus_coex <- 
  data.frame(gene_symbol = pan_cancer_driver_net[,c(1,2)] %>% 
               unlist() %>% 
               unname() %>% 
               unique())

## getSTRINGdb for human
string_db <- STRINGdb$new(species=9606, version = "11.5", 
                          score_threshold=900, # very high confidence (>0.9)
                          input_directory=""
)

pan_cancer_driver_net_genes_plus_coex <- string_db$map(my_data_frame = pan_cancer_driver_net_genes_plus_coex, 
                                                       my_data_frame_id_col_names = "gene_symbol", 
                                                       takeFirst = TRUE, removeUnmappedRows = TRUE)

# Get the PPIs
pan_cancer_driver_net_ppi <- string_db$get_interactions(string_ids = pan_cancer_driver_net_genes_plus_coex$STRING_id)

# Convert stringdb ids to gene symbols

pan_cancer_driver_net_ppi_row_index <- match(pan_cancer_driver_net_ppi$from, pan_cancer_driver_net_genes_plus_coex$STRING_id)
pan_cancer_driver_net_ppi$from <- pan_cancer_driver_net_genes_plus_coex$gene_symbol[pan_cancer_driver_net_ppi_row_index]

pan_cancer_driver_net_ppi_column_index <- match(pan_cancer_driver_net_ppi$to, pan_cancer_driver_net_genes_plus_coex$STRING_id)
pan_cancer_driver_net_ppi$to <- pan_cancer_driver_net_genes_plus_coex$gene_symbol[pan_cancer_driver_net_ppi_column_index]

pan_cancer_driver_net_ppi <- pan_cancer_driver_net_ppi[,c(1,2)]
pan_cancer_driver_net_ppi$type <- "PPI"

## Remove duplicates
pan_cancer_driver_net_ppi <- pan_cancer_driver_net_ppi[-which(duplicated(paste(pan_cancer_driver_net_ppi$from,
                                                                               pan_cancer_driver_net_ppi$to, 
                                                                               sep = "_"))),]

#####################################

# Combine PPI with the original driver net
pan_cancer_driver_net <- rbind(pan_cancer_driver_net, pan_cancer_driver_net_ppi)

# Reconstruct the networks
pan_cancer_driver_net <- igraph::graph_from_data_frame(pan_cancer_driver_net, directed=FALSE)

## Check the network layers
table(igraph::E(pan_cancer_driver_net)$type)

###############################

# Calculate the IVI values of genes within each net
library(influential)
pan_cancer_driver_net_ivi <- ivi(pan_cancer_driver_net)

###############################

# Calculate the primitive gene driver scores
pan_cancer_driver_net_gene_tbl <- data.frame(gene = names(pan_cancer_driver_net_ivi), ivi = pan_cancer_driver_net_ivi)

## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer

pan_cancer_somatic_snv_mean_deletriousness <- colMeans(pan_somatic_snv_tbl)

## Generate the gene table and calculate the primitive scores
pan_cancer_driver_net_gene_tbl$mean_deletriousness <- 0
pan_cancer_driver_net_gene_tbl_match_index <- which(names(pan_cancer_somatic_snv_mean_deletriousness) %in% rownames(pan_cancer_driver_net_gene_tbl))

pan_cancer_driver_net_gene_tbl[names(pan_cancer_somatic_snv_mean_deletriousness)[pan_cancer_driver_net_gene_tbl_match_index], "mean_deletriousness"] <- 
  pan_cancer_somatic_snv_mean_deletriousness[pan_cancer_driver_net_gene_tbl_match_index]

pan_cancer_driver_net_gene_tbl$primitive_driver_score <- (pan_cancer_driver_net_gene_tbl$ivi)*(pan_cancer_driver_net_gene_tbl$mean_deletriousness) # node weights

# add gene scores
pan_cancer_driver_net <- set_vertex_attr(graph = pan_cancer_driver_net, name = "score", value = pan_cancer_driver_net_gene_tbl$primitive_driver_score)

# Convert the un-weighted to a node-weighted network
pan_cancer_driver_net <- set_vertex_attr(graph = pan_cancer_driver_net, name = "weight", value = (pan_cancer_driver_net_gene_tbl$ivi)^(pan_cancer_driver_net_gene_tbl$mean_deletriousness))

## Check the network node weights
summary(igraph::V(pan_cancer_driver_net)$weight)

# Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
pan_cancer_driver_net_gene_tbl$mean_neighborhood_weight <- sapply(1:nrow(pan_cancer_driver_net_gene_tbl), function(j) {
  igraph::neighbors(graph = pan_cancer_driver_net, 
                    v = V(pan_cancer_driver_net)[j])$score %>% mean()
})

pan_cancer_driver_net_gene_tbl$driver_mediator_score <- pan_cancer_driver_net_gene_tbl$primitive_driver_score * pan_cancer_driver_net_gene_tbl$mean_neighborhood_weight

## Top 20 driver genes
pan_cancer_driver_net_gene_tbl %>% 
  slice_max(order_by = primitive_driver_score, n = 20) %>% 
  mutate(Rank = rank(-(primitive_driver_score))) %>% 
  select(Rank, gene) %>% 
  rename(Gene = gene) %>% 
  View()


## Top 30 mediators
pan_cancer_driver_net_gene_tbl %>% 
  slice_max(order_by = driver_mediator_score, n = 30) %>% 
  mutate(Rank = rank(-(driver_mediator_score))) %>% 
  select(Rank, gene) %>% 
  rename(Gene = gene) %>% 
  View()

# Detect communities/modules/clusters
set.seed(3847)
pan_cancer_driver_net_modules <- 
  igraph::cluster_leiden(
    graph = pan_cancer_driver_net,
    objective_function = "CPM",
    weights = NULL,
    resolution_parameter = 0.5,
    beta = 0.05,
    initial_membership = NULL,
    n_iterations = 100,
    vertex_weights = igraph::V(pan_cancer_driver_net)$weight
  )

sizes(pan_cancer_driver_net_modules) %>% as.vector() %>% summary()
length(pan_cancer_driver_net_modules)

## Inspect and filter modules

### Create a table of modules and their genes
pan_cancer_driver_net_modules_tbl <- data.frame(
  module = membership(pan_cancer_driver_net_modules) %>% as.integer(),
  gene = pan_cancer_driver_net_modules$names
)

### Add if the gene is mutated or not (could be a driver gene or not)
pan_cancer_driver_net_modules_tbl$mutated <- 0
pan_cancer_driver_net_modules_tbl$mutated[which(pan_cancer_driver_net_gene_tbl$mean_deletriousness > 0)] <- 1

### Separate modules
pan_cancer_driver_net_modules_tbl <- lapply(unique(pan_cancer_driver_net_modules_tbl$module), function(j){
  subset(pan_cancer_driver_net_modules_tbl, module == j)
})

names(pan_cancer_driver_net_modules_tbl) <- sapply(pan_cancer_driver_net_modules_tbl, function(m) {
  unique(m$module)
})

### Remove modules that have less than 2 driver genes or their size is less than 4
pan_cancer_driver_net_modules_tbl <- lapply(pan_cancer_driver_net_modules_tbl, function(k) {
  if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
    k
  } else {
    NULL
  }
})

### Remove NULL modules
pan_cancer_driver_net_modules_tbl[which(sapply(pan_cancer_driver_net_modules_tbl, is.null))] <- NULL

### Calculate the module scores
pan_cancer_driver_net_modules_tbl <- lapply(pan_cancer_driver_net_modules_tbl, function(l) {
  pan_cancer_driver_net_modules_tbl4score <- cbind(l, mean_score = 0)
  pan_cancer_driver_net_modules_tbl4score$mean_score <- mean(pan_cancer_driver_net_gene_tbl$primitive_driver_score[which(pan_cancer_driver_net_gene_tbl$gene %in% l$gene)])
  pan_cancer_driver_net_modules_tbl4score
})

## Sort the mean scores of all modules
pan_cancer_driver_net_modules_flt_scores <- sapply(pan_cancer_driver_net_modules_tbl, function(j) {unique(j$mean_score)})
pan_cancer_driver_net_modules_flt_scores <- rev(sort(pan_cancer_driver_net_modules_flt_scores))

## Create module 1 for visualization

library(influential)
library(igraph)
library(visNetwork)

pan_cancer_driver_module1 <- subgraph(graph = pan_cancer_driver_net, 
                                      vids = which(as_ids(V(pan_cancer_driver_net)) %in%  pan_cancer_driver_net_modules_tbl$`267`$gene))

pan_cancer_driver_module1_vis <-
  visIgraph(igraph = pan_cancer_driver_module1) %>% 
  visOptions(
    # nodesIdSelection = TRUE, 
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '25px arial black',
    size = 25,
    color = "orange",
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "box" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = c("lightgrey"))

visExport(graph = pan_cancer_driver_module1_vis, 
          name = "pan_cancer_driver_module1_vis", 
          label = "")

#########################################################

# Evaluation of pan-cancer results

## checking the number of common top 20 and the first ranked module genes in the pecan/pedican reported driver genes

### defining a set of pan-cancer pediatric driver genes
pan_cancer_ped_drivers <-
  unique(c(rhabdoid_driver_genes,
           mb_driver_genes,
           lymphoma_driver_genes,
           epd_driver_genes,
           embryonal_cns_driver_genes,
           sarcoma_driver_genes,
           nbl_driver_genes,
           leukaemia_driver_genes,
           hgg_driver_genes))

### create a table of pan_cancer_pecan_prism_common genes

pan_cancer_driver_top_20 <- 
pan_cancer_driver_net_gene_tbl %>% 
  slice_max(order_by = primitive_driver_score, n = 20) %>% 
  select(gene) %>% 
  unlist() %>% 
  unname()

pan_cancer_risk_top_20 <- 
  pan_cancer_risk_net_gene_tbl %>% 
  slice_max(order_by = primitive_risk_score, n = 20) %>% 
  select(gene) %>% 
  unlist() %>% 
  unname()

pan_cancer_pecan_prism_common_top20 <- list("Cancer Driver" = pan_cancer_driver_top_20[which(pan_cancer_driver_top_20 %in% pan_cancer_ped_drivers)],
                                            "Cancer Risk" = pan_cancer_risk_top_20[which(pan_cancer_risk_top_20 %in% pan_cancer_ped_drivers)])

pan_cancer_driver_first_module <- pan_cancer_driver_net_modules_tbl[[names(pan_cancer_driver_net_modules_flt_scores)[2]]]$gene
pan_cancer_risk_first_module <- pan_cancer_risk_net_modules_tbl[[names(pan_cancer_risk_net_modules_flt_scores)[2]]]$gene

pan_cancer_pecan_prism_common_1st_module <- list("Cancer Driver" = pan_cancer_driver_first_module[which(pan_cancer_driver_first_module %in% pan_cancer_ped_drivers)],
                                                 "Cancer Risk" = pan_cancer_risk_first_module[which(pan_cancer_risk_first_module %in% pan_cancer_ped_drivers)])

pan_cancer_pecan_prism_common <- tibble(`Number of Genes in PeCan/PediCan` = c(sapply(pan_cancer_pecan_prism_common_top20, length),
                                                                               sapply(pan_cancer_pecan_prism_common_1st_module, length)),
                                        `Genes in PeCan/PediCan` = c(sapply(pan_cancer_pecan_prism_common_top20, function(i) paste0(i, collapse = ", ")),
                                                                     sapply(pan_cancer_pecan_prism_common_1st_module, function(i) paste0(i, collapse = ", "))),
                                        `Gene group` = rep(c("Top 20 Genes", "First Ranked Module"), each = 2),
                                        `Context` = rep(c("Cancer Driver", "Cancer Risk"), 2))

pan_cancer_pecan_prism_common_plot <- ggplot(data = pan_cancer_pecan_prism_common, aes(x = `Gene group`, 
                                                                       y = `Number of Genes in PeCan/PediCan`, 
                                                                       fill = Context)) +
  geom_bar(position = "dodge2", stat = 'identity') +
  geom_text(aes(x = `Gene group`, 
                y = `Number of Genes in PeCan/PediCan` - `Number of Genes in PeCan/PediCan`/2, 
                label = `Genes in PeCan/PediCan`), position=position_dodge(width=0.9), size = 1.6, show.legend = FALSE) +
  labs(y= '# Genes common with pan-cancer\nPeCan/PediCan drivers') +
  coord_flip() +
  scale_fill_fish_d(begin = 0.27, end = 0.51) +
  # scale_fill_viridis_d(begin = 0.7) +
  theme_bw() +
  theme(text = element_text(size = 7),
        legend.position = "right", 
        legend.box.background = element_rect(colour = "black"))

pan_cancer_pecan_prism_common_plot

ggsave(filename = "Results/Evaluations/pan-cancer/pan_cancer_pecan_prism_common_plot.pdf", 
       plot = pan_cancer_pecan_prism_common_plot, 
       device = "pdf", width = 11.5, height = 4, units = "cm")

########################################

## Literature based evaluation

library(pubmedR)

pan_cancer_lit_data <- tibble(`Gene group` = c(rep(c("Top 20 Genes"), 40),
                                                      rep(c("First Ranked Module"), 
                                                          (length(pan_cancer_driver_first_module) + length(pan_cancer_risk_first_module)))),
                                     Gene = c(pan_cancer_driver_top_20, 
                                              pan_cancer_risk_top_20, 
                                              pan_cancer_driver_first_module, 
                                              pan_cancer_risk_first_module),
                                     `Number of PMIDs` = vector("integer", length = 70),
                                     Context = c(rep(c("Cancer Driver", "Cancer Risk"), each = 20),
                                                 rep(c("Cancer Driver"), length(pan_cancer_driver_first_module)),
                                                 rep(c("Cancer Risk"), length(pan_cancer_risk_first_module)))
)

### Define the pan-cancer main search text
pan_cancer_main_search_text <- "((HGG[Title/Abstract]) OR (HGGs[Title/Abstract]) OR (High-Grade Glioma[Title/Abstract]) OR (High-Grade Gliomas[Title/Abstract]) OR (DMG[Title/Abstract]) OR (DMGs[Title/Abstract]) OR (Diffuse Midline Glioma[Title/Abstract]) OR (Diffuse Midline Gliomas[Title/Abstract]) OR (leukemia[Title/Abstract]) OR (leukemias[Title/Abstract]) OR (Leukaemias[Title/Abstract]) OR (Leukaemia[Title/Abstract]) OR (Lymphoma[Title/Abstract]) OR (Lymphomas[Title/Abstract]) OR (NBL[Title/Abstract]) OR (NBLs[Title/Abstract]) OR (Neuroblastoma[Title/Abstract]) OR (Neuroblastomas[Title/Abstract]) OR (Sarcoma[Title/Abstract]) OR (Sarcomas[Title/Abstract]) OR (atypical teratoid/rhabdoid tumor[Title/Abstract]) OR (AT/RT[Title/Abstract]) OR (embryonal tumors with multilayered rosettes[Title/Abstract]) OR (ETMR[Title/Abstract]) OR (ependymoblastoma[Title/Abstract]) OR (EPD[Title/Abstract]) OR (Ependymoma[Title/Abstract]) OR (Ependymomas[Title/Abstract]) OR (Subependymoma[Title/Abstract]) OR (MB[Title/Abstract]) OR (medulloblastoma[Title/Abstract]) OR (medulloblastomas[Title/Abstract]) OR (Rhabdoid[Title/Abstract]) OR (Renal medullary carcinoma[Title/Abstract]) OR (RMC[Title/Abstract]) OR (Cancer[Title/Abstract]) OR (Carcinoma[Title/Abstract]) OR (Tumor[Title/Abstract]) OR (Tumour[Title/Abstract]))"

pan_cancer_pmid_list <- lapply(1:nrow(pan_cancer_lit_data), function(j) {
  query_text <- paste0("(", pan_cancer_lit_data$Gene[j], "[Title/Abstract]) AND ", pan_cancer_main_search_text)
  tmp_entrez <- pmQueryTotalCount(query = query_text, api_key = "1a1d2b704b3bad87d1ded4f2d95d611e1909")
  tmp_entrez_num <- tmp_entrez$total_count
  return(tmp_entrez_num)
})

pan_cancer_lit_data$`Number of PMIDs` <- unlist(pan_cancer_pmid_list)

### Visualization of the results

pan_cancer_lit_mining_plot <- ggplot(data = pan_cancer_lit_data, aes(x = `Context`, 
                                                             y = log2(`Number of PMIDs`+1), 
                                                             fill = `Gene group`)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.01, width=0.58, size = 0.2) +
  # geom_boxplot(position = position_dodge(0.8), width=0.4, 
  #              aes(color = `Gene group`),
  #              fatten = NULL, fill = NA, coef = 0, outlier.alpha = 0, show.legend = F) +
  labs(y= 'Log2(number of papers + 1)',
       title = "Literature-based Association of Top Pan-cancer\nGenes with Their Respective Cancers") +
  coord_flip() +
  scale_fill_fish_d(option = "Elagatis_bipinnulata") +
  # scale_color_fish_d(option = "Elagatis_bipinnulata") +
  # scale_fill_viridis_d(begin = 0.32, end = 0.83) +
  theme_bw() +
  theme(text = element_text(size = 7),
        legend.position = "right", 
        legend.box.background = element_rect(colour = "black"))

pan_cancer_lit_mining_plot

ggsave(filename = "Results/Evaluations/pan-cancer/pan_cancer_lit_mining_plot.pdf", 
       plot = pan_cancer_lit_mining_plot, 
       device = "pdf", width = 10, height = 4, units = "cm")

########################################

## ORA on the first ranked driver modules

library(enrichR)

### Using enrichR

#Get the list of all databases
enrichR_dbs <- enrichR::listEnrichrDbs()
enrichR_dbs %>% as.data.frame() %>% view()

#Select you desired databases
dbs.GO.BP <- "GO_Biological_Process_2023"
# dbs.KEGG.pathway <- "KEGG_2021_Human"

#ORA on GO-BP terms for 1st driver module of pan-cancer results

pan_cancer_drivers_1st_module_ORA <- 
  enrichR::enrichr(genes = as.character(pan_cancer_driver_first_module),
                   databases = dbs.GO.BP)[[1]]

library(openxlsx)
write.xlsx(x = pan_cancer_drivers_1st_module_ORA, 
           file = "Results/Evaluations/pan-cancer/pan_cancer_drivers_1st_module_ORA.xlsx")

#### summarizing GO terms

# First filter the GO terms based on adjusted p-value
pan_cancer_drivers_1st_module_ORA <- subset(pan_cancer_drivers_1st_module_ORA, Adjusted.P.value < 0.05)

#Separate the GO terms for summarization
pan_cancer_drivers_1st_module_ORA$GO_ID <- ""

for (i in 1:nrow(pan_cancer_drivers_1st_module_ORA)) {
  pan_cancer_drivers_1st_module_ORA$GO_ID[i] <- substr(x = pan_cancer_drivers_1st_module_ORA$Term[i],
                                  start = nchar(pan_cancer_drivers_1st_module_ORA$Term[i])-10,
                                  stop = nchar(pan_cancer_drivers_1st_module_ORA$Term[i])-1)
}

#First check the name of your desired orgDB at:
#http://bioconductor.org/packages/release/BiocViews.html#___OrgDb

#Then install your desired orgDB

#to install all GO terms associated with the desired organism
BiocManager::install("org.Hs.eg.db")
BiocManager::install("GOSemSim")
GOSemSim::godata("org.Hs.eg.db", ont="BP")

#calculating the similarity between GO terms
pan_cancer_drivers_1st_module_ORA_simMatrix <- rrvgo::calculateSimMatrix(pan_cancer_drivers_1st_module_ORA$GO_ID,
                                                orgdb="org.Hs.eg.db",
                                                ont="BP",
                                                method="Rel")

#Reducing the list of GO terms based on their similarities

#First create a named vector of scores associated with each term
pan_cancer_drivers_1st_module_ORA_GO.scores <- setNames(pan_cancer_drivers_1st_module_ORA$Combined.Score, 
                                                        pan_cancer_drivers_1st_module_ORA$GO_ID)

#Now perform the summarization
pan_cancer_drivers_1st_module_ORA_reducedTerms <- rrvgo::reduceSimMatrix(pan_cancer_drivers_1st_module_ORA_simMatrix,
                                                pan_cancer_drivers_1st_module_ORA_GO.scores,
                                                threshold=0.9,
                                                orgdb="org.Hs.eg.db")

############################

#ORA on GO-BP terms for 1st risk module of pan-cancer results
pan_cancer_risks_1st_module_ORA <- 
  enrichR::enrichr(genes = as.character(pan_cancer_risk_first_module),
                   databases = dbs.GO.BP)[[1]]

write.xlsx(x = pan_cancer_risks_1st_module_ORA, 
           file = "Results/Evaluations/pan-cancer/pan_cancer_risks_1st_module_ORA.xlsx")

#### summarizing GO terms

# First filter the GO terms based on adjusted p-value
pan_cancer_risks_1st_module_ORA <- subset(pan_cancer_risks_1st_module_ORA, Adjusted.P.value < 0.05)

#Separate the GO terms for summarization
pan_cancer_risks_1st_module_ORA$GO_ID <- ""

for (i in 1:nrow(pan_cancer_risks_1st_module_ORA)) {
  pan_cancer_risks_1st_module_ORA$GO_ID[i] <- substr(x = pan_cancer_risks_1st_module_ORA$Term[i],
                                                     start = nchar(pan_cancer_risks_1st_module_ORA$Term[i])-10,
                                                     stop = nchar(pan_cancer_risks_1st_module_ORA$Term[i])-1)
}

#calculating the similarity between GO terms
pan_cancer_risks_1st_module_ORA_simMatrix <- rrvgo::calculateSimMatrix(pan_cancer_risks_1st_module_ORA$GO_ID,
                                                                       orgdb="org.Hs.eg.db",
                                                                       ont="BP",
                                                                       method="Rel")

#Reducing the list of GO terms based on their similarities

#First create a named vector of scores associated with each term
pan_cancer_risks_1st_module_ORA_GO.scores <- setNames(pan_cancer_risks_1st_module_ORA$Combined.Score, 
                                                      pan_cancer_risks_1st_module_ORA$GO_ID)

#Now perform the summarization
pan_cancer_risks_1st_module_ORA_reducedTerms <- rrvgo::reduceSimMatrix(pan_cancer_risks_1st_module_ORA_simMatrix,
                                                                       pan_cancer_risks_1st_module_ORA_GO.scores,
                                                                       threshold=0.9,
                                                                       orgdb="org.Hs.eg.db")

########################

## Visualizing a combined alluvial plot
library(ggalluvial)

pan_cancer_risks_1st_module_ORA_reducedTerms <- pan_cancer_risks_1st_module_ORA_reducedTerms[,c("parent", "score", "parentTerm")]

pan_cancer_risks_1st_module_ORA_reducedTerms <- pan_cancer_risks_1st_module_ORA_reducedTerms %>% 
  group_by(parent, parentTerm) %>% 
  summarise(score = mean(score)) %>% 
  ungroup()
pan_cancer_risks_1st_module_ORA_reducedTerms$Context <- "Cancer Risk"
pan_cancer_risks_1st_module_ORA_reducedTerms$`Biological process category` <- c("Lipid metabolism regulation", "Protein function regulation", 
                                                                                "Protein function regulation", "Lipid metabolism regulation",
                                                                                "Gene expression regulation", "Gene expression regulation", 
                                                                                "Lipid metabolism regulation", "Cell-cell interaction", 
                                                                                "Lipid metabolism regulation", "Lipid metabolism regulation",
                                                                                "Lipid metabolism regulation", "Cell-cell interaction", 
                                                                                "Metabolism regulation")

###############

pan_cancer_drivers_1st_module_ORA_reducedTerms <- pan_cancer_drivers_1st_module_ORA_reducedTerms[,c("parent", "score", "parentTerm")]

pan_cancer_drivers_1st_module_ORA_reducedTerms <- pan_cancer_drivers_1st_module_ORA_reducedTerms %>% 
  group_by(parent, parentTerm) %>% 
  summarise(score = mean(score)) %>% 
  ungroup()

pan_cancer_drivers_1st_module_ORA_reducedTerms$Context <- "Cancer Driver"

pan_cancer_drivers_1st_module_ORA_reducedTerms$`Biological process category` <- c("Immune system modulation", "Gene expression regulation",
                                                                                  "Protein function regulation", "Protein function regulation",
                                                                                  "Lipid metabolism regulation", "Lipid metabolism regulation",
                                                                                  "Gene expression regulation", "Lipid metabolism regulation",
                                                                                  "Lipid metabolism regulation", "Gene expression regulation", 
                                                                                  "Cell-cell interaction", "Cell death modulation",
                                                                                  "Lipid metabolism regulation", "Lipid metabolism regulation", 
                                                                                  "Metabolism regulation")

###############

pan_cancer_1st_module_ORA_reducedTerms <- rbind(pan_cancer_risks_1st_module_ORA_reducedTerms, 
                                                pan_cancer_drivers_1st_module_ORA_reducedTerms)


pan_cancer_1st_module_ORA_reducedTerms$`Major function category` <- "Non-metabolic biological process"
pan_cancer_1st_module_ORA_reducedTerms$`Major function category`[grep("metab", pan_cancer_1st_module_ORA_reducedTerms$`Biological process category`, ignore.case = T)] <- "Metabolic biological process"

pan_cancer_1st_module_ORA_reducedTerms$parent <- factor(pan_cancer_1st_module_ORA_reducedTerms$parent, 
                                                        levels = unique(c(pan_cancer_1st_module_ORA_reducedTerms$parent[which(pan_cancer_1st_module_ORA_reducedTerms$`Major function category` == "Metabolic biological process")],
                                                                   pan_cancer_1st_module_ORA_reducedTerms$parent[which(pan_cancer_1st_module_ORA_reducedTerms$`Major function category` == "Non-metabolic biological process")])))

pan_cancer_1st_module_ORA_reducedTerms_plot <-
ggplot(data = pan_cancer_1st_module_ORA_reducedTerms,
       aes(axis1 = `Major function category`,   # First variable on the X-axis
           axis2 = Context, # Second variable on the X-axis
           axis3 = parent,   # Third variable on the X-axis
           y = score)) +
  geom_alluvium(aes(fill = `Biological process category`)) +
  geom_stratum(linewidth = 0.1) +
  geom_text(stat = "stratum", 
            # angle=c(rep(90, 4),rep(0, 24)),
            size = 2, 
            aes(label = after_stat(stratum))) +
  scale_x_discrete(expand = c(0.15, 0.05)) +
  scale_fill_fish_d() +
  theme_void() +
  theme(text = element_text(size = 8), legend.box.background = element_rect(colour = "black"))

pan_cancer_1st_module_ORA_reducedTerms_plot

#Saving the plot
ggsave(filename = paste0("Results/Evaluations/pan-cancer/pan_cancer_1st_module_ORA_reducedTerms_plot.pdf"), 
       plot = pan_cancer_1st_module_ORA_reducedTerms_plot,
       device = "pdf",width = 17, height = 8, units = "cm")

###########

#### Visualization of the connections of the first ranked risk and driver mediators in the five first ranked modules

pan_cancer_risk_net_gene_tbl %>% slice_max(risk_mediator_score) %>% select(gene) %>% unlist() %>% unname() # TP53

pan_cancer_risk_first_ranked_mediator_Neighbors <- igraph::neighbors(graph = pan_cancer_risk_net, 
                                                                     v = "TP53", mode = "all") %>% as_ids()

# separate only the ones present in the top modules
pan_cancer_risk_first_ranked_mediator_Neighbors <-
  pan_cancer_risk_first_ranked_mediator_Neighbors[which(pan_cancer_risk_first_ranked_mediator_Neighbors %in% (sapply(pan_cancer_risk_net_modules_tbl[names(pan_cancer_risk_net_modules_flt_scores)[1:5]], function(i) i$gene) %>% unlist() %>% unname()))]

pan_cancer_first_ranked_risk_mediator_neighbors_net <- subgraph(graph = pan_cancer_risk_net, 
                                                         vids = which(as_ids(V(pan_cancer_risk_net)) %in%  c("TP53", pan_cancer_risk_first_ranked_mediator_Neighbors)))

## check connections with different modules
pan_cancer_risk_first_ranked_mediator_Neighbors[which(pan_cancer_risk_first_ranked_mediator_Neighbors %in% pan_cancer_risk_net_modules_tbl[[names(pan_cancer_risk_net_modules_flt_scores)[5]]]$gene)]

pan_cancer_first_ranked_risk_mediator_neighbors_net_tbl <- igraph::as_data_frame(pan_cancer_first_ranked_risk_mediator_neighbors_net)
readr::write_delim(x = pan_cancer_first_ranked_risk_mediator_neighbors_net_tbl,
                   file = "Results/Evaluations/pan-cancer/first_ranked_risk_mediator_neighbors_tbl.txt", 
                   delim = "\t")

############################

pan_cancer_driver_net_gene_tbl %>% slice_max(driver_mediator_score) %>% select(gene) %>% unlist() %>% unname() # TP53

pan_cancer_driver_first_ranked_mediator_Neighbors <- igraph::neighbors(graph = pan_cancer_driver_net, 
                                                                       v = "TP53", mode = "all") %>% as_ids()

# separate only the ones present in the top modules
pan_cancer_driver_first_ranked_mediator_Neighbors <-
  pan_cancer_driver_first_ranked_mediator_Neighbors[which(pan_cancer_driver_first_ranked_mediator_Neighbors %in% (sapply(pan_cancer_driver_net_modules_tbl[names(pan_cancer_driver_net_modules_flt_scores)[1:5]], function(i) i$gene) %>% unlist() %>% unname()))]

pan_cancer_first_ranked_driver_mediator_neighbors_net <- subgraph(graph = pan_cancer_driver_net, 
                                                                  vids = which(as_ids(V(pan_cancer_driver_net)) %in%  c("TP53", pan_cancer_driver_first_ranked_mediator_Neighbors)))

## check connections with different modules
pan_cancer_driver_first_ranked_mediator_Neighbors[which(pan_cancer_driver_first_ranked_mediator_Neighbors %in% pan_cancer_driver_net_modules_tbl[[names(pan_cancer_driver_net_modules_flt_scores)[5]]]$gene)]

pan_cancer_first_ranked_driver_mediator_neighbors_net_tbl <- igraph::as_data_frame(pan_cancer_first_ranked_driver_mediator_neighbors_net)
readr::write_delim(x = pan_cancer_first_ranked_driver_mediator_neighbors_net_tbl,
                   file = "Results/Evaluations/pan-cancer/first_ranked_driver_mediator_neighbors_tbl.txt", 
                   delim = "\t")


########################################

# save pan cancer tables

tmp_pan_cancer_driver_net_modules_tbl <-
  lapply(1:length(pan_cancer_driver_net_modules_tbl), function(j) {
    tmp_data <- pan_cancer_driver_net_modules_tbl[[j]]
    colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
    tmp_data[tmp_data == 1] <- "Yes"
    tmp_data[tmp_data == 0] <- "No"
    tmp_data$Module <- paste0("Module number ", j)
    tmp_data <- data.frame(Module = unique(tmp_data$Module),
                           Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                           Mean_score = unique(tmp_data$Mean_score))
    tmp_data
  })
tmp_pan_cancer_driver_net_modules_tbl <- do.call(rbind, tmp_pan_cancer_driver_net_modules_tbl)

library(openxlsx)
write.xlsx(x = tmp_pan_cancer_driver_net_modules_tbl, file = "Results/Evaluations/pan-cancer/pan_cancer_driver_net_modules_flt.xlsx")

rm(tmp_pan_cancer_driver_net_modules_tbl)

########################

tmp_pan_cancer_driver_net_gene_tbl <- pan_cancer_driver_net_gene_tbl
colnames(tmp_pan_cancer_driver_net_gene_tbl) <- c("Gene", "IVI", "Mean_deletriousness_score", "Driver_score", "Mean_neighborhood_weight", "Driverness_mediator_score")
write.xlsx(x = tmp_pan_cancer_driver_net_gene_tbl, file = "Results/Evaluations/pan-cancer/pan_cancer_driver_net_gene_tbl.xlsx")
rm(tmp_pan_cancer_driver_net_gene_tbl)

########################

tmp_pan_cancer_risk_net_gene_tbl <- pan_cancer_risk_net_gene_tbl
colnames(tmp_pan_cancer_risk_net_gene_tbl) <- c("Gene", "IVI", "Mean_deletriousness_score", "Risk_score", "Mean_neighborhood_weight", "Risk_mediator_score")
write.xlsx(x = tmp_pan_cancer_risk_net_gene_tbl, file = "Results/Evaluations/pan-cancer/pan_cancer_risk_net_gene_tbl.xlsx")
rm(tmp_pan_cancer_risk_net_gene_tbl)

########################

tmp_pan_cancer_risk_net_modules_tbl <-
  lapply(1:length(pan_cancer_risk_net_modules_tbl), function(j) {
    tmp_data <- pan_cancer_risk_net_modules_tbl[[j]]
    colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
    tmp_data[tmp_data == 1] <- "Yes"
    tmp_data[tmp_data == 0] <- "No"
    tmp_data$Module <- paste0("Module number ", j)
    tmp_data <- data.frame(Module = unique(tmp_data$Module),
                           Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                           Mean_score = unique(tmp_data$Mean_score))
    tmp_data
  })
tmp_pan_cancer_risk_net_modules_tbl <- do.call(rbind, tmp_pan_cancer_risk_net_modules_tbl)

library(openxlsx)
write.xlsx(x = tmp_pan_cancer_risk_net_modules_tbl, file = "Results/Evaluations/pan-cancer/pan_cancer_risk_net_modules_flt.xlsx")

rm(tmp_pan_cancer_risk_net_modules_tbl)

#=============================================================================
#
#    Code chunk 11: Test the pipeline based on permutation
#
#=============================================================================

# Test the pipeline based on permutation (shuffling both deleteriousness and expression data)

## Prepare somatic snv data

permut_somatic_snv_flt_withCADD_summarized <- somatic_snv_flt_withCADD_summarized
## Remove some uncertain cancer types
permut_somatic_snv_flt_withCADD_summarized <- permut_somatic_snv_flt_withCADD_summarized[-c(1,10)]

## Shuffle the data
permut_somatic_snv_flt_withCADD_summarized <- lapply(permut_somatic_snv_flt_withCADD_summarized, function(i) {
  tmp_colnames <- colnames(i)
  
  ## shuffle
  tmp_data <- apply(i, 1, function(j) {
    set.seed(1234)
    j[sample(length(j))]
  })
  tmp_data <- as.data.frame(t(tmp_data))
  colnames(tmp_data) <- tmp_colnames
  
  ## re-shuffle
  tmp_data <- apply(i, 1, function(j) {
    set.seed(5678)
    j[sample(length(j))]
  })
  tmp_data <- as.data.frame(t(tmp_data))
  colnames(tmp_data) <- tmp_colnames
  
  tmp_data
})

## Perform correlation analysis

permut_somatic_snv_flt_cor <- lapply(permut_somatic_snv_flt_withCADD_summarized, function(i){
  tmp_cor <- fcor(i)
  tmp_cor <- subset(tmp_cor, mr < 20)
  tmp_cor
})

##################

## Prepare expression data
permut_rnaseq_expr_list <- rnaseq_expr_list
## Remove some uncertain cancer types
permut_rnaseq_expr_list <- permut_rnaseq_expr_list[names(permut_somatic_snv_flt_cor)]

## Shuffle the data
permut_rnaseq_expr_list <- lapply(permut_rnaseq_expr_list, function(i) {
  tmp_colnames <- colnames(i)
  
  ## shuffle
  tmp_data <- apply(i, 1, function(j) {
    set.seed(1234)
    j[sample(length(j))]
  })
  tmp_data <- as.data.frame(t(tmp_data))
  colnames(tmp_data) <- tmp_colnames
  
  ## re-shuffle
  tmp_data <- apply(i, 1, function(j) {
    set.seed(5678)
    j[sample(length(j))]
  })
  tmp_data <- as.data.frame(t(tmp_data))
  colnames(tmp_data) <- tmp_colnames
  
  tmp_data
})

## Perform correlation analysis
permut_rnaseq_expr_cor <- lapply(permut_rnaseq_expr_list, function(i){
  
  tmp_cor <- fcor(data = i, 
                  method = "spearman", 
                  mutualRank = TRUE, pvalue = FALSE, flat = TRUE)
  
  tmp_cor <- subset(tmp_cor, mr < 20)
  tmp_cor
})

###################

# Prepare the tables for merging with other data
permut_somatic_snv_flt_cor <- lapply(permut_somatic_snv_flt_cor, function(i) {
  tmp_tbl <- i[,c(1,2)]
  colnames(tmp_tbl) <- c("from", "to")
  tmp_tbl$type <- "co-deletriousness"
  tmp_tbl
})

permut_rnaseq_expr_cor <- lapply(permut_rnaseq_expr_cor, function(i) {
  tmp_tbl <- i[,c(1,2)]
  colnames(tmp_tbl) <- c("from", "to")
  tmp_tbl$type <- "co-expression"
  tmp_tbl
})

# Merge the similarity tables
permut_driver_marker_net <- permut_somatic_snv_flt_cor

# Add the co-expression data to the tables

for(i in 1:length(permut_driver_marker_net)) {
  
  genes <- permut_driver_marker_net[[i]][,c(1,2)] %>% 
    unlist() %>% 
    unname() %>% 
    unique()
  
  first_level_coex_index <- c(which(permut_rnaseq_expr_cor[[i]]$from %in% genes), 
                              which(permut_rnaseq_expr_cor[[i]]$to %in% genes)) %>% unique()
  
  first_level_coex_genes <- c(genes, 
                              permut_rnaseq_expr_cor[[i]]$from[first_level_coex_index],
                              permut_rnaseq_expr_cor[[i]]$to[first_level_coex_index]) %>% unique()
  
  first_two_levels_coex_index <- c(which(permut_rnaseq_expr_cor[[i]]$from %in% first_level_coex_genes), 
                                   which(permut_rnaseq_expr_cor[[i]]$to %in% first_level_coex_genes)) %>% unique()
  
  permut_driver_marker_net[[i]] <- rbind(permut_driver_marker_net[[i]], permut_rnaseq_expr_cor[[i]][first_two_levels_coex_index,])
  
}

# Add PPI data to the driver/marker net

# Map gene names to stringDB ids for PPI analysis (STRING  11.5)

## get the gene names in each risk net
permut_driver_marker_net_genes <- lapply(permut_driver_marker_net, function(i) {
  data.frame(gene_symbol = i[,c(1,2)] %>% 
               unlist() %>% 
               unname() %>% 
               unique())
})

## getSTRINGdb for human
### WARNING: You didn't specify a version of the STRING database to use. Hence we will use STRING  11.5
string_db <- STRINGdb$new(species=9606, version = "11.5", 
                          score_threshold=900, # very high confidence (>0.9)
                          input_directory=""
)

permut_driver_marker_net_genes <- lapply(permut_driver_marker_net_genes, function(i) {
  
  string_db$map(my_data_frame = i, 
                my_data_frame_id_col_names = "gene_symbol", 
                takeFirst = TRUE, removeUnmappedRows = TRUE)
  
})

# Get the PPIs
permut_driver_marker_net_ppi <- lapply(permut_driver_marker_net_genes, function(i) {
  
  string_db$get_interactions(string_ids = i$STRING_id)
  
})

# Convert stringdb ids to gene symbols
permut_driver_marker_net_ppi <- lapply(1:length(permut_driver_marker_net_ppi), function(i) {
  tmp_row_index <- match(permut_driver_marker_net_ppi[[i]]$from, permut_driver_marker_net_genes[[i]]$STRING_id)
  permut_driver_marker_net_ppi[[i]]$from <- permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_row_index]
  
  tmp_column_index <- match(permut_driver_marker_net_ppi[[i]]$to, permut_driver_marker_net_genes[[i]]$STRING_id)
  permut_driver_marker_net_ppi[[i]]$to <- permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_column_index]
  
  tmp_tbl <- permut_driver_marker_net_ppi[[i]][,c(1,2)]
  tmp_tbl$type <- "PPI"
  
  ## Remove duplicates
  tmp_tbl <- tmp_tbl[-which(duplicated(paste(tmp_tbl$from,
                                             tmp_tbl$to, 
                                             sep = "_"))),]
  
  return(tmp_tbl)
})

names(permut_driver_marker_net_ppi) <- names(permut_driver_marker_net_genes)

# Combine PPI with the original driver/marker net
permut_driver_marker_net <- lapply(1:length(permut_driver_marker_net), function(i) {
  rbind(permut_driver_marker_net[[i]], permut_driver_marker_net_ppi[[i]])
})

names(permut_driver_marker_net) <- names(permut_driver_marker_net_genes)

# Reconstruct the networks
permut_driver_marker_net <- lapply(permut_driver_marker_net, function(i) {
  igraph::graph_from_data_frame(i, directed=FALSE)
})

## Check the network layers
table(igraph::E(permut_driver_marker_net$DMG)$type)

###############################

# Calculate the IVI values of genes within each net
library(influential)
permut_driver_marker_net_ivi <- lapply(permut_driver_marker_net, function(i) ivi(i, verbose = T))

###############################

# Calculate the primitive gene driver scores
permut_driver_marker_net_gene_tbl <- lapply(permut_driver_marker_net_ivi, function(i) {
  data.frame(gene = names(i), ivi = i)
})

## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer

permut_somatic_snv_mean_deletriousness <- lapply(permut_somatic_snv_flt_withCADD_summarized, function(i) {
  colMeans(i)
})

## Generate the gene table for each net and calculate the primitive scores
permut_driver_marker_net_gene_tbl <- lapply(1:length(permut_driver_marker_net_gene_tbl), function(i) {
  tmp_tbl <- permut_driver_marker_net_gene_tbl[[i]]
  tmp_tbl$mean_deletriousness <- 0
  match_index <- which(names(permut_somatic_snv_mean_deletriousness[[i]]) %in% rownames(tmp_tbl))
  tmp_tbl[names(permut_somatic_snv_mean_deletriousness[[i]])[match_index], "mean_deletriousness"] <- 
    permut_somatic_snv_mean_deletriousness[[i]][match_index]
  
  tmp_tbl$primitive_driver_score <- (tmp_tbl$ivi)*(tmp_tbl$mean_deletriousness) # node weights
  
  return(tmp_tbl)
  
})

names(permut_driver_marker_net_gene_tbl) <- names(permut_somatic_snv_mean_deletriousness)

# add gene scores
permut_driver_marker_net <- lapply(1:length(permut_driver_marker_net), function(i) {
  set_vertex_attr(graph = permut_driver_marker_net[[i]], name = "score", value = permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score)
})

names(permut_driver_marker_net) <- names(permut_somatic_snv_mean_deletriousness)

# Convert the un-weighted to a node-weighted network
permut_driver_marker_net <- lapply(1:length(permut_driver_marker_net), function(i) {
  set_vertex_attr(graph = permut_driver_marker_net[[i]], name = "weight", value = (permut_driver_marker_net_gene_tbl[[i]]$ivi)^(permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness))
})

names(permut_driver_marker_net) <- names(permut_somatic_snv_mean_deletriousness)

## Check the network node weights
summary(igraph::V(permut_driver_marker_net$DMG)$weight)

# Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
permut_driver_marker_net_gene_tbl <- lapply(1:length(permut_driver_marker_net), function(i) {
  tmp_tbl <- permut_driver_marker_net_gene_tbl[[i]]
  tmp_tbl$mean_neighborhood_weight <- sapply(1:nrow(tmp_tbl), function(j) {
    igraph::neighbors(graph = permut_driver_marker_net[[i]], 
                      v = V(permut_driver_marker_net[[i]])[j])$score %>% mean()
  })
  
  tmp_tbl$mediator_score <- tmp_tbl$primitive_driver_score * tmp_tbl$mean_neighborhood_weight
  
  return(tmp_tbl)
})

names(permut_driver_marker_net_gene_tbl) <- names(permut_somatic_snv_mean_deletriousness)

# Detect communities/modules/clusters
permut_driver_marker_net_modules <- lapply(permut_driver_marker_net, function(i) {
  set.seed(3847)
  igraph::cluster_leiden(
    graph = i,
    objective_function = "CPM",
    weights = NULL,
    resolution_parameter = 0.5,
    beta = 0.05,
    initial_membership = NULL,
    n_iterations = 100,
    vertex_weights = igraph::V(i)$weight
  )
})

sizes(permut_driver_marker_net_modules$DMG) %>% as.vector() %>% summary()
length(permut_driver_marker_net_modules$DMG)

## Inspect and filter modules
permut_driver_marker_net_modules_flt <- lapply(1:length(permut_driver_marker_net_modules), function(i) {
  
  ### Create a table of modules and their genes
  tmp_tbl <- data.frame(
    module = membership(permut_driver_marker_net_modules[[i]]) %>% as.integer(),
    gene = permut_driver_marker_net_modules[[i]]$names
  )
  
  ### Add if the gene is mutated or not (could be a risk gene or not)
  tmp_tbl$mutated <- 0
  tmp_tbl$mutated[which(permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness > 0)] <- 1
  
  ### Separate modules
  tmp_tbl <- lapply(unique(tmp_tbl$module), function(j){
    subset(tmp_tbl, module == j)
  })
  
  names(tmp_tbl) <- sapply(tmp_tbl, function(m) {
    unique(m$module)
  })
  
  ### Remove modules that have less than 2 Risk genes or their size is less than 4
  tmp_tbl <- lapply(tmp_tbl, function(k) {
    if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
      k
    } else {
      NULL
    }
  })
  
  ### Remove NULL modules
  tmp_tbl[which(sapply(tmp_tbl, is.null))] <- NULL
  
  ### Calculate the module scores
  tmp_tbl <- lapply(tmp_tbl, function(l) {
    tmp_tbl4score <- cbind(l, mean_score = 0)
    tmp_tbl4score$mean_score <- mean(permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score[which(permut_driver_marker_net_gene_tbl[[i]]$gene %in% l$gene)])
    tmp_tbl4score
  })
  
  tmp_tbl
})

names(permut_driver_marker_net_modules_flt) <- names(permut_driver_marker_net_modules)

## Sort the mean scores of all modules of each cancer type
permut_driver_marker_net_modules_flt_scores <-
  lapply(permut_driver_marker_net_modules_flt, function(i) {
    tmp_tbl <-
      sapply(i, function(j) {unique(j$mean_score)})
    tmp_tbl <- rev(sort(tmp_tbl))
    tmp_tbl
  })

### View the number of first ranked modules
sapply(permut_driver_marker_net_modules_flt_scores, function(i) {i[1]}) %>% as.data.frame() %>% View()

permut_driver_marker_net_modules_flt$HGG$`1667`

########

# Top 20 gene 
permut_driver_marker_net_gene_tbl_top20 <-
lapply(permut_driver_marker_net_gene_tbl, function(i) {
  slice_max(i, n = 20, order_by = primitive_driver_score) %>% select(gene)
})

permut_driver_marker_net_gene_tbl_top20 <- do.call(cbind, permut_driver_marker_net_gene_tbl_top20)
colnames(permut_driver_marker_net_gene_tbl_top20) <- names(permut_driver_marker_net_gene_tbl)

###########################################
###########################################
###########################################

## Save the gene data and modules as supplementary files

tmp_permut_driver_marker_net_gene_tbl <- permut_driver_marker_net_gene_tbl
names(tmp_permut_driver_marker_net_gene_tbl)[names(tmp_permut_driver_marker_net_gene_tbl) == "Rhabdoid"] <- "Rhabdoid carcinoma"

for(i in 1:10) {
  colnames(tmp_permut_driver_marker_net_gene_tbl[[i]]) <- c("Gene", "IVI", "Mean_deletriousness_score", "Driver_score", "Mean_neighborhood_weight", "Driverness_mediator_score")
}

library(openxlsx)
write.xlsx(x = tmp_permut_driver_marker_net_gene_tbl, file = "Results/Evaluations/Permutation/permut_driver_marker_net_gene_tbl.xlsx")

rm(tmp_permut_driver_marker_net_gene_tbl)

############################

### save the modules

tmp_permut_driver_marker_net_modules_flt <- permut_driver_marker_net_modules_flt
names(tmp_permut_driver_marker_net_modules_flt)[names(tmp_permut_driver_marker_net_modules_flt) == "Rhabdoid"] <- "Rhabdoid carcinoma"

tmp_permut_driver_marker_net_modules_flt_scores <- permut_driver_marker_net_modules_flt_scores

tmp_permut_driver_marker_net_modules_flt <- lapply(tmp_permut_driver_marker_net_modules_flt, function(i) {
  
  tmp_data <- i
  tmp_data <- tmp_data[rev(order(sapply(tmp_data, function(j) unique(j$mean_score))))]
  
  tmp_data
  
})

for(i in 1:10) {
  tmp_permut_driver_marker_net_modules_flt[[i]] <-
    lapply(1:length(tmp_permut_driver_marker_net_modules_flt[[i]]), function(j) {
      tmp_data <- tmp_permut_driver_marker_net_modules_flt[[i]][[j]]
      colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
      tmp_data[tmp_data == 1] <- "Yes"
      tmp_data[tmp_data == 0] <- "No"
      tmp_data$Module <- paste0("Module number ", j)
      tmp_data <- data.frame(Module = unique(tmp_data$Module),
                             Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                             Mean_score = unique(tmp_data$Mean_score))
      tmp_data
    })
  tmp_permut_driver_marker_net_modules_flt[[i]] <- do.call(rbind, tmp_permut_driver_marker_net_modules_flt[[i]])
}

library(openxlsx)
write.xlsx(x = tmp_permut_driver_marker_net_modules_flt, file = "Results/Evaluations/Permutation/permut_driver_marker_net_modules_flt.xlsx")

rm(tmp_permut_driver_marker_net_modules_flt, tmp_permut_driver_marker_net_modules_flt_scores)

###############################################################################################################

# Test the pipeline based on permutation (shuffling everything) 

## Prepare cor data

### Setting snv cor to the ones derived from permutated data
plusPPI_permut_somatic_snv_flt_cor <- permut_somatic_snv_flt_cor

### Setting expr cor to the ones derived from permutated data
plusPPI_permut_rnaseq_expr_cor <- permut_rnaseq_expr_cor
plusPPI_permut_rnaseq_expr_cor <- plusPPI_permut_rnaseq_expr_cor[which(names(plusPPI_permut_rnaseq_expr_cor) %in% names(plusPPI_permut_somatic_snv_flt_cor))]

# Merge the similarity tables
plusPPI_permut_driver_marker_net <- plusPPI_permut_somatic_snv_flt_cor

# Add the co-expression data to the tables

for(i in 1:length(plusPPI_permut_driver_marker_net)) {
  
  genes <- plusPPI_permut_driver_marker_net[[i]][,c(1,2)] %>% 
    unlist() %>% 
    unname() %>% 
    unique()
  
  first_level_coex_index <- c(which(plusPPI_permut_rnaseq_expr_cor[[i]]$from %in% genes), 
                              which(plusPPI_permut_rnaseq_expr_cor[[i]]$to %in% genes)) %>% unique()
  
  first_level_coex_genes <- c(genes, 
                              plusPPI_permut_rnaseq_expr_cor[[i]]$from[first_level_coex_index],
                              plusPPI_permut_rnaseq_expr_cor[[i]]$to[first_level_coex_index]) %>% unique()
  
  first_two_levels_coex_index <- c(which(plusPPI_permut_rnaseq_expr_cor[[i]]$from %in% first_level_coex_genes), 
                                   which(plusPPI_permut_rnaseq_expr_cor[[i]]$to %in% first_level_coex_genes)) %>% unique()
  
  plusPPI_permut_driver_marker_net[[i]] <- rbind(plusPPI_permut_driver_marker_net[[i]], plusPPI_permut_rnaseq_expr_cor[[i]][first_two_levels_coex_index,])
  
}

# Permute PPI data
plusPPI_permut_driver_marker_net_ppi <- driver_marker_net_ppi[-c(1,10)]
plusPPI_permut_driver_marker_net_ppi <- lapply(plusPPI_permut_driver_marker_net_ppi, function(i) {
  tmp_data <- i
  
  ## Remove duplicates
  tmp_data <- tmp_data[-which(duplicated(paste(tmp_data$from,
                                               tmp_data$to, 
                                             sep = "_"))),]
  
  ## Shuffle the data
  set.seed(1234)
  tmp_data$from <- tmp_data$from[sample(nrow(tmp_data))]
  # set.seed(5678)
  # tmp_data$from <- tmp_data$from[sample(nrow(tmp_data))]
  tmp_data
})

# Add PPI data to the driver/marker net

# Combine PPI with the original driver/marker net
plusPPI_permut_driver_marker_net <- lapply(1:length(plusPPI_permut_driver_marker_net), function(i) {
  rbind(plusPPI_permut_driver_marker_net[[i]], plusPPI_permut_driver_marker_net_ppi[[i]])
})

names(plusPPI_permut_driver_marker_net) <- names(plusPPI_permut_driver_marker_net_ppi)

# Reconstruct the networks
plusPPI_permut_driver_marker_net <- lapply(plusPPI_permut_driver_marker_net, function(i) {
  igraph::graph_from_data_frame(i, directed=FALSE)
})

## Check the network layers
table(igraph::E(plusPPI_permut_driver_marker_net$DMG)$type)

###############################

# Calculate the IVI values of genes within each net
library(influential)
plusPPI_permut_driver_marker_net_ivi <- lapply(plusPPI_permut_driver_marker_net, function(i) ivi(i, verbose = T))

###############################

# Calculate the primitive gene driver scores
plusPPI_permut_driver_marker_net_gene_tbl <- lapply(plusPPI_permut_driver_marker_net_ivi, function(i) {
  data.frame(gene = names(i), ivi = i)
})

## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer

plusPPI_permut_somatic_snv_mean_deletriousness <- lapply(permut_somatic_snv_flt_withCADD_summarized, function(i) {
  colMeans(i)
})

## Generate the gene table for each net and calculate the primitive scores
plusPPI_permut_driver_marker_net_gene_tbl <- lapply(1:length(plusPPI_permut_driver_marker_net_gene_tbl), function(i) {
  tmp_tbl <- plusPPI_permut_driver_marker_net_gene_tbl[[i]]
  tmp_tbl$mean_deletriousness <- 0
  match_index <- which(names(plusPPI_permut_somatic_snv_mean_deletriousness[[i]]) %in% rownames(tmp_tbl))
  tmp_tbl[names(plusPPI_permut_somatic_snv_mean_deletriousness[[i]])[match_index], "mean_deletriousness"] <- 
    plusPPI_permut_somatic_snv_mean_deletriousness[[i]][match_index]
  
  tmp_tbl$primitive_driver_score <- (tmp_tbl$ivi)*(tmp_tbl$mean_deletriousness) # node weights
  
  return(tmp_tbl)
  
})

names(plusPPI_permut_driver_marker_net_gene_tbl) <- names(plusPPI_permut_somatic_snv_mean_deletriousness)

# add gene scores
plusPPI_permut_driver_marker_net <- lapply(1:length(plusPPI_permut_driver_marker_net), function(i) {
  set_vertex_attr(graph = plusPPI_permut_driver_marker_net[[i]], name = "score", value = plusPPI_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score)
})

names(plusPPI_permut_driver_marker_net) <- names(plusPPI_permut_somatic_snv_mean_deletriousness)

# Convert the un-weighted to a node-weighted network
plusPPI_permut_driver_marker_net <- lapply(1:length(plusPPI_permut_driver_marker_net), function(i) {
  set_vertex_attr(graph = plusPPI_permut_driver_marker_net[[i]], name = "weight", value = (plusPPI_permut_driver_marker_net_gene_tbl[[i]]$ivi)^plusPPI_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness)
})

names(plusPPI_permut_driver_marker_net) <- names(plusPPI_permut_somatic_snv_mean_deletriousness)

## Check the network node weights
summary(igraph::V(plusPPI_permut_driver_marker_net$DMG)$weight)

# Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
plusPPI_permut_driver_marker_net_gene_tbl <- lapply(1:length(plusPPI_permut_driver_marker_net), function(i) {
  tmp_tbl <- plusPPI_permut_driver_marker_net_gene_tbl[[i]]
  tmp_tbl$mean_neighborhood_weight <- sapply(1:nrow(tmp_tbl), function(j) {
    igraph::neighbors(graph = plusPPI_permut_driver_marker_net[[i]], 
                      v = V(plusPPI_permut_driver_marker_net[[i]])[j])$score %>% mean()
  })
  
  tmp_tbl$mediator_score <- tmp_tbl$primitive_driver_score * tmp_tbl$mean_neighborhood_weight
  
  return(tmp_tbl)
})

names(plusPPI_permut_driver_marker_net_gene_tbl) <- names(plusPPI_permut_somatic_snv_mean_deletriousness)

# Detect communities/modules/clusters
plusPPI_permut_driver_marker_net_modules <- lapply(plusPPI_permut_driver_marker_net, function(i) {
  set.seed(3847)
  igraph::cluster_leiden(
    graph = i,
    objective_function = "CPM",
    weights = NULL,
    resolution_parameter = 0.5,
    beta = 0.05,
    initial_membership = NULL,
    n_iterations = 100,
    vertex_weights = igraph::V(i)$weight
  )
})

sizes(plusPPI_permut_driver_marker_net_modules$DMG) %>% as.vector() %>% summary()
length(plusPPI_permut_driver_marker_net_modules$DMG)

## Inspect and filter modules
plusPPI_permut_driver_marker_net_modules_flt <- lapply(1:length(plusPPI_permut_driver_marker_net_modules), function(i) {
  
  ### Create a table of modules and their genes
  tmp_tbl <- data.frame(
    module = membership(plusPPI_permut_driver_marker_net_modules[[i]]) %>% as.integer(),
    gene = plusPPI_permut_driver_marker_net_modules[[i]]$names
  )
  
  ### Add if the gene is mutated or not (could be a risk gene or not)
  tmp_tbl$mutated <- 0
  tmp_tbl$mutated[which(plusPPI_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness > 0)] <- 1
  
  ### Separate modules
  tmp_tbl <- lapply(unique(tmp_tbl$module), function(j){
    subset(tmp_tbl, module == j)
  })
  
  names(tmp_tbl) <- sapply(tmp_tbl, function(m) {
    unique(m$module)
  })
  
  ### Remove modules that have less than 2 Risk genes or their size is less than 4
  tmp_tbl <- lapply(tmp_tbl, function(k) {
    if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
      k
    } else {
      NULL
    }
  })
  
  ### Remove NULL modules
  tmp_tbl[which(sapply(tmp_tbl, is.null))] <- NULL
  
  ### Calculate the module scores
  tmp_tbl <- lapply(tmp_tbl, function(l) {
    tmp_tbl4score <- cbind(l, mean_score = 0)
    tmp_tbl4score$mean_score <- mean(plusPPI_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score[which(plusPPI_permut_driver_marker_net_gene_tbl[[i]]$gene %in% l$gene)])
    tmp_tbl4score
  })
  
  tmp_tbl
})

names(plusPPI_permut_driver_marker_net_modules_flt) <- names(plusPPI_permut_driver_marker_net_modules)

sapply(plusPPI_permut_driver_marker_net_modules_flt, length)
sapply(driver_marker_net_modules_flt, length)

## Sort the mean scores of all modules of each cancer type
plusPPI_permut_driver_marker_net_modules_flt_scores <-
  lapply(plusPPI_permut_driver_marker_net_modules_flt, function(i) {
    tmp_tbl <-
      sapply(i, function(j) {unique(j$mean_score)})
    tmp_tbl <- rev(sort(tmp_tbl))
    tmp_tbl
  })

### View the number of first ranked modules
sapply(plusPPI_permut_driver_marker_net_modules_flt_scores, function(i) {i[1]}) %>% as.data.frame() %>% View()

plusPPI_permut_driver_marker_net_modules_flt$HGG$`1462` %>% View()

########

# Top 20 gene 
plusPPI_permut_driver_marker_net_gene_tbl_top20 <-
  lapply(plusPPI_permut_driver_marker_net_gene_tbl, function(i) {
    slice_max(i, n = 20, order_by = primitive_driver_score) %>% select(gene)
  })

plusPPI_permut_driver_marker_net_gene_tbl_top20 <- do.call(cbind, plusPPI_permut_driver_marker_net_gene_tbl_top20)
colnames(plusPPI_permut_driver_marker_net_gene_tbl_top20) <- names(plusPPI_permut_driver_marker_net_gene_tbl)

###########################################
###########################################
###########################################

## Save the gene data and modules as supplementary files

tmp_plusPPI_permut_driver_marker_net_gene_tbl <- plusPPI_permut_driver_marker_net_gene_tbl
names(tmp_plusPPI_permut_driver_marker_net_gene_tbl)[names(tmp_plusPPI_permut_driver_marker_net_gene_tbl) == "Rhabdoid"] <- "Rhabdoid carcinoma"

for(i in 1:10) {
  colnames(tmp_plusPPI_permut_driver_marker_net_gene_tbl[[i]]) <- c("Gene", "IVI", "Mean_deletriousness_score", "Driver_score", "Mean_neighborhood_weight", "Driverness_mediator_score")
}

library(openxlsx)
write.xlsx(x = tmp_plusPPI_permut_driver_marker_net_gene_tbl, file = "Results/Evaluations/Permutation/plusPPI_permut_driver_marker_net_gene_tbl.xlsx")

rm(tmp_plusPPI_permut_driver_marker_net_gene_tbl)

############################

### save the modules

tmp_plusPPI_permut_driver_marker_net_modules_flt <- plusPPI_permut_driver_marker_net_modules_flt
names(tmp_plusPPI_permut_driver_marker_net_modules_flt)[names(tmp_plusPPI_permut_driver_marker_net_modules_flt) == "Rhabdoid"] <- "Rhabdoid carcinoma"

tmp_plusPPI_permut_driver_marker_net_modules_flt_scores <- plusPPI_permut_driver_marker_net_modules_flt_scores

tmp_plusPPI_permut_driver_marker_net_modules_flt <- lapply(tmp_plusPPI_permut_driver_marker_net_modules_flt, function(i) {
  
  tmp_data <- i
  tmp_data <- tmp_data[rev(order(sapply(tmp_data, function(j) unique(j$mean_score))))]
  
  tmp_data
  
})

for(i in 1:10) {
  tmp_plusPPI_permut_driver_marker_net_modules_flt[[i]] <-
    lapply(1:length(tmp_plusPPI_permut_driver_marker_net_modules_flt[[i]]), function(j) {
      tmp_data <- tmp_plusPPI_permut_driver_marker_net_modules_flt[[i]][[j]]
      colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
      tmp_data[tmp_data == 1] <- "Yes"
      tmp_data[tmp_data == 0] <- "No"
      tmp_data$Module <- paste0("Module number ", j)
      tmp_data <- data.frame(Module = unique(tmp_data$Module),
                             Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                             Mean_score = unique(tmp_data$Mean_score))
      tmp_data
    })
  tmp_plusPPI_permut_driver_marker_net_modules_flt[[i]] <- do.call(rbind, tmp_plusPPI_permut_driver_marker_net_modules_flt[[i]])
}

library(openxlsx)
write.xlsx(x = tmp_plusPPI_permut_driver_marker_net_modules_flt, file = "Results/Evaluations/Permutation/plusPPI_permut_driver_marker_net_modules_flt.xlsx")

rm(tmp_plusPPI_permut_driver_marker_net_modules_flt, tmp_plusPPI_permut_driver_marker_net_modules_flt_scores)

###############################################################################################################
###############################################################################################################
###############################################################################################################

# Test the pipeline based on permutation (shuffling everything but deleteriousness data) 

## Prepare cor data

### Setting snv cor to the ones derived from original data
noSNV_permut_somatic_snv_flt_cor <- somatic_snv_flt_cor[-c(1,10)]

### Setting expr cor to the ones derived from permutated data
noSNV_permut_rnaseq_expr_cor <- permut_rnaseq_expr_cor
noSNV_permut_rnaseq_expr_cor <- noSNV_permut_rnaseq_expr_cor[which(names(noSNV_permut_rnaseq_expr_cor) %in% names(noSNV_permut_somatic_snv_flt_cor))]

# Merge the similarity tables
noSNV_permut_driver_marker_net <- noSNV_permut_somatic_snv_flt_cor

# Add the co-expression data to the tables

for(i in 1:length(noSNV_permut_driver_marker_net)) {
  
  genes <- noSNV_permut_driver_marker_net[[i]][,c(1,2)] %>% 
    unlist() %>% 
    unname() %>% 
    unique()
  
  first_level_coex_index <- c(which(noSNV_permut_rnaseq_expr_cor[[i]]$from %in% genes), 
                              which(noSNV_permut_rnaseq_expr_cor[[i]]$to %in% genes)) %>% unique()
  
  first_level_coex_genes <- c(genes, 
                              noSNV_permut_rnaseq_expr_cor[[i]]$from[first_level_coex_index],
                              noSNV_permut_rnaseq_expr_cor[[i]]$to[first_level_coex_index]) %>% unique()
  
  first_two_levels_coex_index <- c(which(noSNV_permut_rnaseq_expr_cor[[i]]$from %in% first_level_coex_genes), 
                                   which(noSNV_permut_rnaseq_expr_cor[[i]]$to %in% first_level_coex_genes)) %>% unique()
  
  noSNV_permut_driver_marker_net[[i]] <- rbind(noSNV_permut_driver_marker_net[[i]], noSNV_permut_rnaseq_expr_cor[[i]][first_two_levels_coex_index,])
  
}

### Setting PPI to the ones derived from permutated data
noSNV_permut_driver_marker_net_ppi <- plusPPI_permut_driver_marker_net_ppi

# Add PPI data to the driver/marker net

# Combine PPI with the original driver/marker net
noSNV_permut_driver_marker_net <- lapply(1:length(noSNV_permut_driver_marker_net), function(i) {
  rbind(noSNV_permut_driver_marker_net[[i]], noSNV_permut_driver_marker_net_ppi[[i]])
})

names(noSNV_permut_driver_marker_net) <- names(noSNV_permut_driver_marker_net_ppi)

# Reconstruct the networks
noSNV_permut_driver_marker_net <- lapply(noSNV_permut_driver_marker_net, function(i) {
  igraph::graph_from_data_frame(i, directed=FALSE)
})

## Check the network layers
table(igraph::E(noSNV_permut_driver_marker_net$DMG)$type)

###############################

# Calculate the IVI values of genes within each net
library(influential)
noSNV_permut_driver_marker_net_ivi <- lapply(noSNV_permut_driver_marker_net, function(i) ivi(i, verbose = T))

###############################

# Calculate the primitive gene driver scores
noSNV_permut_driver_marker_net_gene_tbl <- lapply(noSNV_permut_driver_marker_net_ivi, function(i) {
  data.frame(gene = names(i), ivi = i)
})

## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer

noSNV_permut_somatic_snv_mean_deletriousness <- lapply(somatic_snv_flt_withCADD_summarized[-c(1,10)], function(i) {
  colMeans(i)
})

## Generate the gene table for each net and calculate the primitive scores
noSNV_permut_driver_marker_net_gene_tbl <- lapply(1:length(noSNV_permut_driver_marker_net_gene_tbl), function(i) {
  tmp_tbl <- noSNV_permut_driver_marker_net_gene_tbl[[i]]
  tmp_tbl$mean_deletriousness <- 0
  match_index <- which(names(noSNV_permut_somatic_snv_mean_deletriousness[[i]]) %in% rownames(tmp_tbl))
  tmp_tbl[names(noSNV_permut_somatic_snv_mean_deletriousness[[i]])[match_index], "mean_deletriousness"] <- 
    noSNV_permut_somatic_snv_mean_deletriousness[[i]][match_index]
  
  tmp_tbl$primitive_driver_score <- (tmp_tbl$ivi)*(tmp_tbl$mean_deletriousness) # node weights
  
  return(tmp_tbl)
  
})

names(noSNV_permut_driver_marker_net_gene_tbl) <- names(noSNV_permut_somatic_snv_mean_deletriousness)

# add gene scores
noSNV_permut_driver_marker_net <- lapply(1:length(noSNV_permut_driver_marker_net), function(i) {
  set_vertex_attr(graph = noSNV_permut_driver_marker_net[[i]], name = "score", value = noSNV_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score)
})

names(noSNV_permut_driver_marker_net) <- names(noSNV_permut_somatic_snv_mean_deletriousness)

# Convert the un-weighted to a node-weighted network
noSNV_permut_driver_marker_net <- lapply(1:length(noSNV_permut_driver_marker_net), function(i) {
  set_vertex_attr(graph = noSNV_permut_driver_marker_net[[i]], name = "weight", value = (noSNV_permut_driver_marker_net_gene_tbl[[i]]$ivi)^(noSNV_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness))
})

names(noSNV_permut_driver_marker_net) <- names(noSNV_permut_somatic_snv_mean_deletriousness)

## Check the network node weights
summary(igraph::V(noSNV_permut_driver_marker_net$DMG)$weight)

# Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
noSNV_permut_driver_marker_net_gene_tbl <- lapply(1:length(noSNV_permut_driver_marker_net), function(i) {
  tmp_tbl <- noSNV_permut_driver_marker_net_gene_tbl[[i]]
  tmp_tbl$mean_neighborhood_weight <- sapply(1:nrow(tmp_tbl), function(j) {
    igraph::neighbors(graph = noSNV_permut_driver_marker_net[[i]], 
                      v = V(noSNV_permut_driver_marker_net[[i]])[j])$score %>% mean()
  })
  
  tmp_tbl$mediator_score <- tmp_tbl$primitive_driver_score * tmp_tbl$mean_neighborhood_weight
  
  return(tmp_tbl)
})

names(noSNV_permut_driver_marker_net_gene_tbl) <- names(noSNV_permut_somatic_snv_mean_deletriousness)

# Detect communities/modules/clusters
noSNV_permut_driver_marker_net_modules <- lapply(noSNV_permut_driver_marker_net, function(i) {
  set.seed(3847)
  igraph::cluster_leiden(
    graph = i,
    objective_function = "CPM",
    weights = NULL,
    resolution_parameter = 0.5,
    beta = 0.05,
    initial_membership = NULL,
    n_iterations = 100,
    vertex_weights = igraph::V(i)$weight
  )
})

sizes(noSNV_permut_driver_marker_net_modules$DMG) %>% as.vector() %>% summary()
length(noSNV_permut_driver_marker_net_modules$DMG)

## Inspect and filter modules
noSNV_permut_driver_marker_net_modules_flt <- lapply(1:length(noSNV_permut_driver_marker_net_modules), function(i) {
  
  ### Create a table of modules and their genes
  tmp_tbl <- data.frame(
    module = membership(noSNV_permut_driver_marker_net_modules[[i]]) %>% as.integer(),
    gene = noSNV_permut_driver_marker_net_modules[[i]]$names
  )
  
  ### Add if the gene is mutated or not (could be a risk gene or not)
  tmp_tbl$mutated <- 0
  tmp_tbl$mutated[which(noSNV_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness > 0)] <- 1
  
  ### Separate modules
  tmp_tbl <- lapply(unique(tmp_tbl$module), function(j){
    subset(tmp_tbl, module == j)
  })
  
  names(tmp_tbl) <- sapply(tmp_tbl, function(m) {
    unique(m$module)
  })
  
  ### Remove modules that have less than 2 Risk genes or their size is less than 4
  tmp_tbl <- lapply(tmp_tbl, function(k) {
    if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
      k
    } else {
      NULL
    }
  })
  
  ### Remove NULL modules
  tmp_tbl[which(sapply(tmp_tbl, is.null))] <- NULL
  
  ### Calculate the module scores
  tmp_tbl <- lapply(tmp_tbl, function(l) {
    tmp_tbl4score <- cbind(l, mean_score = 0)
    tmp_tbl4score$mean_score <- mean(noSNV_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score[which(noSNV_permut_driver_marker_net_gene_tbl[[i]]$gene %in% l$gene)])
    tmp_tbl4score
  })
  
  tmp_tbl
})

names(noSNV_permut_driver_marker_net_modules_flt) <- names(noSNV_permut_driver_marker_net_modules)

sapply(noSNV_permut_driver_marker_net_modules_flt, length)
sapply(driver_marker_net_modules_flt, length)

## Sort the mean scores of all modules of each cancer type
noSNV_permut_driver_marker_net_modules_flt_scores <-
  lapply(noSNV_permut_driver_marker_net_modules_flt, function(i) {
    tmp_tbl <-
      sapply(i, function(j) {unique(j$mean_score)})
    tmp_tbl <- rev(sort(tmp_tbl))
    tmp_tbl
  })

### View the number of first ranked modules
sapply(noSNV_permut_driver_marker_net_modules_flt_scores, function(i) {i[1]}) %>% as.data.frame() %>% View()

noSNV_permut_driver_marker_net_modules_flt$HGG$`1462` %>% View()

########

# Top 20 gene 
noSNV_permut_driver_marker_net_gene_tbl_top20 <-
  lapply(noSNV_permut_driver_marker_net_gene_tbl, function(i) {
    slice_max(i, n = 20, order_by = primitive_driver_score) %>% select(gene)
  })

noSNV_permut_driver_marker_net_gene_tbl_top20 <- do.call(cbind, noSNV_permut_driver_marker_net_gene_tbl_top20)
colnames(noSNV_permut_driver_marker_net_gene_tbl_top20) <- names(noSNV_permut_driver_marker_net_gene_tbl)

###########################################
###########################################
###########################################

## Save the gene data and modules as supplementary files

tmp_noSNV_permut_driver_marker_net_gene_tbl <- noSNV_permut_driver_marker_net_gene_tbl
names(tmp_noSNV_permut_driver_marker_net_gene_tbl)[names(tmp_noSNV_permut_driver_marker_net_gene_tbl) == "Rhabdoid"] <- "Rhabdoid carcinoma"

for(i in 1:10) {
  colnames(tmp_noSNV_permut_driver_marker_net_gene_tbl[[i]]) <- c("Gene", "IVI", "Mean_deletriousness_score", "Driver_score", "Mean_neighborhood_weight", "Driverness_mediator_score")
}

library(openxlsx)
write.xlsx(x = tmp_noSNV_permut_driver_marker_net_gene_tbl, file = "Results/Evaluations/Permutation/noSNV_permut_driver_marker_net_gene_tbl.xlsx")

rm(tmp_noSNV_permut_driver_marker_net_gene_tbl)

############################

### save the modules

tmp_noSNV_permut_driver_marker_net_modules_flt <- noSNV_permut_driver_marker_net_modules_flt
names(tmp_noSNV_permut_driver_marker_net_modules_flt)[names(tmp_noSNV_permut_driver_marker_net_modules_flt) == "Rhabdoid"] <- "Rhabdoid carcinoma"

tmp_noSNV_permut_driver_marker_net_modules_flt_scores <- noSNV_permut_driver_marker_net_modules_flt_scores

tmp_noSNV_permut_driver_marker_net_modules_flt <- lapply(tmp_noSNV_permut_driver_marker_net_modules_flt, function(i) {
  
  tmp_data <- i
  tmp_data <- tmp_data[rev(order(sapply(tmp_data, function(j) unique(j$mean_score))))]
  
  tmp_data
  
})

for(i in 1:10) {
  tmp_noSNV_permut_driver_marker_net_modules_flt[[i]] <-
    lapply(1:length(tmp_noSNV_permut_driver_marker_net_modules_flt[[i]]), function(j) {
      tmp_data <- tmp_noSNV_permut_driver_marker_net_modules_flt[[i]][[j]]
      colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
      tmp_data[tmp_data == 1] <- "Yes"
      tmp_data[tmp_data == 0] <- "No"
      tmp_data$Module <- paste0("Module number ", j)
      tmp_data <- data.frame(Module = unique(tmp_data$Module),
                             Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                             Mean_score = unique(tmp_data$Mean_score))
      tmp_data
    })
  tmp_noSNV_permut_driver_marker_net_modules_flt[[i]] <- do.call(rbind, tmp_noSNV_permut_driver_marker_net_modules_flt[[i]])
}

library(openxlsx)
write.xlsx(x = tmp_noSNV_permut_driver_marker_net_modules_flt, file = "Results/Evaluations/Permutation/noSNV_permut_driver_marker_net_modules_flt.xlsx")

rm(tmp_noSNV_permut_driver_marker_net_modules_flt, tmp_noSNV_permut_driver_marker_net_modules_flt_scores)

###############################################################################################################
###############################################################################################################
###############################################################################################################

# Test the pipeline based on permutation (shuffling only deleteriousness data) 

## Prepare cor data

### Setting snv cor to the ones derived from permutated data
onlySNV_permut_somatic_snv_flt_cor <- permut_somatic_snv_flt_cor

### Setting expr cor to the ones derived from NON-permutated data
onlySNV_permut_rnaseq_expr_cor <- rnaseq_expr_cor
onlySNV_permut_rnaseq_expr_cor <- onlySNV_permut_rnaseq_expr_cor[which(names(onlySNV_permut_rnaseq_expr_cor) %in% names(onlySNV_permut_somatic_snv_flt_cor))]

# Merge the similarity tables
onlySNV_permut_driver_marker_net <- onlySNV_permut_somatic_snv_flt_cor

# Add the co-expression data to the tables

for(i in 1:length(onlySNV_permut_driver_marker_net)) {
  
  genes <- onlySNV_permut_driver_marker_net[[i]][,c(1,2)] %>% 
    unlist() %>% 
    unname() %>% 
    unique()
  
  first_level_coex_index <- c(which(onlySNV_permut_rnaseq_expr_cor[[i]]$from %in% genes), 
                              which(onlySNV_permut_rnaseq_expr_cor[[i]]$to %in% genes)) %>% unique()
  
  first_level_coex_genes <- c(genes, 
                              onlySNV_permut_rnaseq_expr_cor[[i]]$from[first_level_coex_index],
                              onlySNV_permut_rnaseq_expr_cor[[i]]$to[first_level_coex_index]) %>% unique()
  
  first_two_levels_coex_index <- c(which(onlySNV_permut_rnaseq_expr_cor[[i]]$from %in% first_level_coex_genes), 
                                   which(onlySNV_permut_rnaseq_expr_cor[[i]]$to %in% first_level_coex_genes)) %>% unique()
  
  onlySNV_permut_driver_marker_net[[i]] <- rbind(onlySNV_permut_driver_marker_net[[i]], onlySNV_permut_rnaseq_expr_cor[[i]][first_two_levels_coex_index,])
  
}

# Add PPI data to the driver/marker net

# Map gene names to stringDB ids for PPI analysis (STRING  11.5)

## get the gene names in each risk net
onlySNV_permut_driver_marker_net_genes <- lapply(onlySNV_permut_driver_marker_net, function(i) {
  data.frame(gene_symbol = i[,c(1,2)] %>% 
               unlist() %>% 
               unname() %>% 
               unique())
})

## getSTRINGdb for human
string_db <- STRINGdb$new(species=9606, version = "11.5", 
                          score_threshold=900, # very high confidence (>0.9)
                          input_directory=""
)

onlySNV_permut_driver_marker_net_genes <- lapply(onlySNV_permut_driver_marker_net_genes, function(i) {
  
  string_db$map(my_data_frame = i, 
                my_data_frame_id_col_names = "gene_symbol", 
                takeFirst = TRUE, removeUnmappedRows = TRUE)
  
})

# Get the PPIs
onlySNV_permut_driver_marker_net_ppi <- lapply(onlySNV_permut_driver_marker_net_genes, function(i) {
  
  string_db$get_interactions(string_ids = i$STRING_id)
  
})

# Convert stringdb ids to gene symbols
onlySNV_permut_driver_marker_net_ppi <- lapply(1:length(onlySNV_permut_driver_marker_net_ppi), function(i) {
  tmp_row_index <- match(onlySNV_permut_driver_marker_net_ppi[[i]]$from, onlySNV_permut_driver_marker_net_genes[[i]]$STRING_id)
  onlySNV_permut_driver_marker_net_ppi[[i]]$from <- onlySNV_permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_row_index]
  
  tmp_column_index <- match(onlySNV_permut_driver_marker_net_ppi[[i]]$to, onlySNV_permut_driver_marker_net_genes[[i]]$STRING_id)
  onlySNV_permut_driver_marker_net_ppi[[i]]$to <- onlySNV_permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_column_index]
  
  tmp_tbl <- onlySNV_permut_driver_marker_net_ppi[[i]][,c(1,2)]
  tmp_tbl$type <- "PPI"
  
  ## Remove duplicates
  tmp_tbl <- tmp_tbl[-which(duplicated(paste(tmp_tbl$from,
                                             tmp_tbl$to, 
                                             sep = "_"))),]
  
  return(tmp_tbl)
})

names(onlySNV_permut_driver_marker_net_ppi) <- names(onlySNV_permut_driver_marker_net_genes)

# Combine PPI with the original driver/marker net
onlySNV_permut_driver_marker_net <- lapply(1:length(onlySNV_permut_driver_marker_net), function(i) {
  rbind(onlySNV_permut_driver_marker_net[[i]], onlySNV_permut_driver_marker_net_ppi[[i]])
})

names(onlySNV_permut_driver_marker_net) <- names(onlySNV_permut_driver_marker_net_genes)

# Reconstruct the networks
onlySNV_permut_driver_marker_net <- lapply(onlySNV_permut_driver_marker_net, function(i) {
  igraph::graph_from_data_frame(i, directed=FALSE)
})

## Check the network layers
table(igraph::E(onlySNV_permut_driver_marker_net$DMG)$type)

###############################

# Calculate the IVI values of genes within each net
library(influential)
onlySNV_permut_driver_marker_net_ivi <- lapply(onlySNV_permut_driver_marker_net, function(i) ivi(i, verbose = T))

###############################

# Calculate the primitive gene driver scores
onlySNV_permut_driver_marker_net_gene_tbl <- lapply(onlySNV_permut_driver_marker_net_ivi, function(i) {
  data.frame(gene = names(i), ivi = i)
})

## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer

onlySNV_permut_somatic_snv_mean_deletriousness <- lapply(permut_somatic_snv_flt_withCADD_summarized, function(i) {
  colMeans(i)
})

## Generate the gene table for each net and calculate the primitive scores
onlySNV_permut_driver_marker_net_gene_tbl <- lapply(1:length(onlySNV_permut_driver_marker_net_gene_tbl), function(i) {
  tmp_tbl <- onlySNV_permut_driver_marker_net_gene_tbl[[i]]
  tmp_tbl$mean_deletriousness <- 0
  match_index <- which(names(onlySNV_permut_somatic_snv_mean_deletriousness[[i]]) %in% rownames(tmp_tbl))
  tmp_tbl[names(onlySNV_permut_somatic_snv_mean_deletriousness[[i]])[match_index], "mean_deletriousness"] <- 
    onlySNV_permut_somatic_snv_mean_deletriousness[[i]][match_index]
  
  tmp_tbl$primitive_driver_score <- (tmp_tbl$ivi)*(tmp_tbl$mean_deletriousness) # node weights
  
  return(tmp_tbl)
  
})

names(onlySNV_permut_driver_marker_net_gene_tbl) <- names(onlySNV_permut_somatic_snv_mean_deletriousness)

# add gene scores
onlySNV_permut_driver_marker_net <- lapply(1:length(onlySNV_permut_driver_marker_net), function(i) {
  set_vertex_attr(graph = onlySNV_permut_driver_marker_net[[i]], name = "score", value = onlySNV_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score)
})

names(onlySNV_permut_driver_marker_net) <- names(onlySNV_permut_somatic_snv_mean_deletriousness)

# Convert the un-weighted to a node-weighted network
onlySNV_permut_driver_marker_net <- lapply(1:length(onlySNV_permut_driver_marker_net), function(i) {
  set_vertex_attr(graph = onlySNV_permut_driver_marker_net[[i]], name = "weight", value = (onlySNV_permut_driver_marker_net_gene_tbl[[i]]$ivi)^(onlySNV_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness))
})

names(onlySNV_permut_driver_marker_net) <- names(onlySNV_permut_somatic_snv_mean_deletriousness)

## Check the network node weights
summary(igraph::V(onlySNV_permut_driver_marker_net$DMG)$weight)

# Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
onlySNV_permut_driver_marker_net_gene_tbl <- lapply(1:length(onlySNV_permut_driver_marker_net), function(i) {
  tmp_tbl <- onlySNV_permut_driver_marker_net_gene_tbl[[i]]
  tmp_tbl$mean_neighborhood_weight <- sapply(1:nrow(tmp_tbl), function(j) {
    igraph::neighbors(graph = onlySNV_permut_driver_marker_net[[i]], 
                      v = V(onlySNV_permut_driver_marker_net[[i]])[j])$score %>% mean()
  })
  
  tmp_tbl$mediator_score <- tmp_tbl$primitive_driver_score * tmp_tbl$mean_neighborhood_weight
  
  return(tmp_tbl)
})

names(onlySNV_permut_driver_marker_net_gene_tbl) <- names(onlySNV_permut_somatic_snv_mean_deletriousness)

# Detect communities/modules/clusters
onlySNV_permut_driver_marker_net_modules <- lapply(onlySNV_permut_driver_marker_net, function(i) {
  set.seed(3847)
  igraph::cluster_leiden(
    graph = i,
    objective_function = "CPM",
    weights = NULL,
    resolution_parameter = 0.5,
    beta = 0.05,
    initial_membership = NULL,
    n_iterations = 100,
    vertex_weights = igraph::V(i)$weight
  )
})

sizes(onlySNV_permut_driver_marker_net_modules$DMG) %>% as.vector() %>% summary()
length(onlySNV_permut_driver_marker_net_modules$DMG)

## Inspect and filter modules
onlySNV_permut_driver_marker_net_modules_flt <- lapply(1:length(onlySNV_permut_driver_marker_net_modules), function(i) {
  
  ### Create a table of modules and their genes
  tmp_tbl <- data.frame(
    module = membership(onlySNV_permut_driver_marker_net_modules[[i]]) %>% as.integer(),
    gene = onlySNV_permut_driver_marker_net_modules[[i]]$names
  )
  
  ### Add if the gene is mutated or not (could be a risk gene or not)
  tmp_tbl$mutated <- 0
  tmp_tbl$mutated[which(onlySNV_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness > 0)] <- 1
  
  ### Separate modules
  tmp_tbl <- lapply(unique(tmp_tbl$module), function(j){
    subset(tmp_tbl, module == j)
  })
  
  names(tmp_tbl) <- sapply(tmp_tbl, function(m) {
    unique(m$module)
  })
  
  ### Remove modules that have less than 2 Risk genes or their size is less than 4
  tmp_tbl <- lapply(tmp_tbl, function(k) {
    if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
      k
    } else {
      NULL
    }
  })
  
  ### Remove NULL modules
  tmp_tbl[which(sapply(tmp_tbl, is.null))] <- NULL
  
  ### Calculate the module scores
  tmp_tbl <- lapply(tmp_tbl, function(l) {
    tmp_tbl4score <- cbind(l, mean_score = 0)
    tmp_tbl4score$mean_score <- mean(onlySNV_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score[which(onlySNV_permut_driver_marker_net_gene_tbl[[i]]$gene %in% l$gene)])
    tmp_tbl4score
  })
  
  tmp_tbl
})

names(onlySNV_permut_driver_marker_net_modules_flt) <- names(onlySNV_permut_driver_marker_net_modules)

sapply(onlySNV_permut_driver_marker_net_modules_flt, length)
sapply(driver_marker_net_modules_flt, length)

## Sort the mean scores of all modules of each cancer type
onlySNV_permut_driver_marker_net_modules_flt_scores <-
  lapply(onlySNV_permut_driver_marker_net_modules_flt, function(i) {
    tmp_tbl <-
      sapply(i, function(j) {unique(j$mean_score)})
    tmp_tbl <- rev(sort(tmp_tbl))
    tmp_tbl
  })

### View the number of first ranked modules
sapply(onlySNV_permut_driver_marker_net_modules_flt_scores, function(i) {i[1]}) %>% as.data.frame() %>% View()

onlySNV_permut_driver_marker_net_modules_flt$HGG$`1462` %>% View()

########

# Top 20 gene 
onlySNV_permut_driver_marker_net_gene_tbl_top20 <-
  lapply(onlySNV_permut_driver_marker_net_gene_tbl, function(i) {
    slice_max(i, n = 20, order_by = primitive_driver_score) %>% select(gene)
  })

onlySNV_permut_driver_marker_net_gene_tbl_top20 <- do.call(cbind, onlySNV_permut_driver_marker_net_gene_tbl_top20)
colnames(onlySNV_permut_driver_marker_net_gene_tbl_top20) <- names(onlySNV_permut_driver_marker_net_gene_tbl)

###########################################
###########################################
###########################################

## Save the gene data and modules as supplementary files

tmp_onlySNV_permut_driver_marker_net_gene_tbl <- onlySNV_permut_driver_marker_net_gene_tbl
names(tmp_onlySNV_permut_driver_marker_net_gene_tbl)[names(tmp_onlySNV_permut_driver_marker_net_gene_tbl) == "Rhabdoid"] <- "Rhabdoid carcinoma"

for(i in 1:10) {
  colnames(tmp_onlySNV_permut_driver_marker_net_gene_tbl[[i]]) <- c("Gene", "IVI", "Mean_deletriousness_score", "Driver_score", "Mean_neighborhood_weight", "Driverness_mediator_score")
}

library(openxlsx)
write.xlsx(x = tmp_onlySNV_permut_driver_marker_net_gene_tbl, file = "Results/Evaluations/Permutation/onlySNV_permut_driver_marker_net_gene_tbl.xlsx")

rm(tmp_onlySNV_permut_driver_marker_net_gene_tbl)

############################

### save the modules

tmp_onlySNV_permut_driver_marker_net_modules_flt <- onlySNV_permut_driver_marker_net_modules_flt
names(tmp_onlySNV_permut_driver_marker_net_modules_flt)[names(tmp_onlySNV_permut_driver_marker_net_modules_flt) == "Rhabdoid"] <- "Rhabdoid carcinoma"

tmp_onlySNV_permut_driver_marker_net_modules_flt_scores <- onlySNV_permut_driver_marker_net_modules_flt_scores

tmp_onlySNV_permut_driver_marker_net_modules_flt <- lapply(tmp_onlySNV_permut_driver_marker_net_modules_flt, function(i) {
  
  tmp_data <- i
  tmp_data <- tmp_data[rev(order(sapply(tmp_data, function(j) unique(j$mean_score))))]
  
  tmp_data
  
})

for(i in 1:10) {
  tmp_onlySNV_permut_driver_marker_net_modules_flt[[i]] <-
    lapply(1:length(tmp_onlySNV_permut_driver_marker_net_modules_flt[[i]]), function(j) {
      tmp_data <- tmp_onlySNV_permut_driver_marker_net_modules_flt[[i]][[j]]
      colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
      tmp_data[tmp_data == 1] <- "Yes"
      tmp_data[tmp_data == 0] <- "No"
      tmp_data$Module <- paste0("Module number ", j)
      tmp_data <- data.frame(Module = unique(tmp_data$Module),
                             Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                             Mean_score = unique(tmp_data$Mean_score))
      tmp_data
    })
  tmp_onlySNV_permut_driver_marker_net_modules_flt[[i]] <- do.call(rbind, tmp_onlySNV_permut_driver_marker_net_modules_flt[[i]])
}

library(openxlsx)
write.xlsx(x = tmp_onlySNV_permut_driver_marker_net_modules_flt, file = "Results/Evaluations/Permutation/onlySNV_permut_driver_marker_net_modules_flt.xlsx")

rm(tmp_onlySNV_permut_driver_marker_net_modules_flt, tmp_onlySNV_permut_driver_marker_net_modules_flt_scores)

#############################################################################################

# Visualize permutation results

########################################
# define a function for the evaluation of intersections

intersect_combn <- function(named_list) {
  
  comb_list = combn(seq(length(named_list)), 2)
  
  intersects =
  lapply(1:ncol(comb_list), function(i) {
    
    tmp_intersect =
    base::Reduce(intersect, named_list[comb_list[,i]])
    
    intersect_tbl = data.frame(Vector1 = names(named_list[comb_list[,i][1]]),
                               Vector2 = names(named_list[comb_list[,i][2]]),
                               Intersection_length = length(tmp_intersect),
                               Intersects = paste0(tmp_intersect, collapse = ", "))
    
    intersect_tbl
  })
  
  intersects <- do.call(rbind, intersects)
  intersects
}

########################################

## for the top 20 genes

onlySNV_permut_driver_marker_net_gene_tbl_top20
onlySNV_permut_driver_marker_net_modules_flt
driver_net_top20_genes

permut_eval_list <- lapply(1:10, function(i) {
  
  original = driver_net_top20_genes[,-c(1,2,11)][,i]
  onlySNV = onlySNV_permut_driver_marker_net_gene_tbl_top20[,i]
  SNV_RNA = permut_driver_marker_net_gene_tbl_top20[,i]
  SNV_RNA_PPI = plusPPI_permut_driver_marker_net_gene_tbl_top20[,i]
  RNA_PPI = noSNV_permut_driver_marker_net_gene_tbl_top20[,i]
  
  intersect_combn(named_list = list("Original data" = original,
                                    "Permuted mutations" = onlySNV,
                                    "Permuted mutations+RNA-Seq" = SNV_RNA,
                                    "Permuted mutations+RNA-Seq+PPIs" = SNV_RNA_PPI,
                                    "Permuted RNA-Seq+PPIs" = RNA_PPI))
})

names(permut_eval_list) <- c(names(driver_marker_net_modules_flt[-c(1,10,12)]), "Rhabdoid carcinoma")

permut_eval_tbl <- do.call(rbind, permut_eval_list)
permut_eval_tbl$Cancer_type <- rep(names(permut_eval_list), each = 10)

colnames(permut_eval_tbl)[c(1,2)] <- paste("Origin_of_gene_set_", c(1,2), sep = "")

#########################

## for the top 5 modules

permut_modules_eval_list <-
lapply(1:10, function(i) {
  
  ## prepare the top 5 modules
  original_modules <- driver_marker_net_modules_flt[-c(1,10)][[i]]
  original_module_scores <- driver_marker_net_modules_flt_scores[-c(1,10)][[i]]
  original_first_five_modules = sapply(1, function(j) {
    tmp_module =
      original_modules[[names(original_module_scores)[j]]]

    tmp_module$gene
  }) %>% unlist()
  ###########
  onlySNV_modules <- onlySNV_permut_driver_marker_net_modules_flt[[i]]
  onlySNV_module_scores <- onlySNV_permut_driver_marker_net_modules_flt_scores[[i]]
  onlySNV_first_five_modules = sapply(1, function(j) {
    tmp_module =
      onlySNV_modules[[names(onlySNV_module_scores)[j]]]
    
    tmp_module$gene
  }) %>% unlist()
  ###########
  SNV_RNA_modules <- permut_driver_marker_net_modules_flt[[i]]
  SNV_RNA_module_scores <- permut_driver_marker_net_modules_flt_scores[[i]]
  SNV_RNA_first_five_modules = sapply(1, function(j) {
    tmp_module =
      SNV_RNA_modules[[names(SNV_RNA_module_scores)[j]]]
    
    tmp_module$gene
  }) %>% unlist()
  ###########
  SNV_RNA_PPI_modules <- plusPPI_permut_driver_marker_net_modules_flt[[i]]
  SNV_RNA_PPI_module_scores <- plusPPI_permut_driver_marker_net_modules_flt_scores[[i]]
  SNV_RNA_PPI_first_five_modules = sapply(1, function(j) {
    tmp_module =
      SNV_RNA_PPI_modules[[names(SNV_RNA_PPI_module_scores)[j]]]
    
    tmp_module$gene
  }) %>% unlist()
  ###########
  RNA_PPI_modules <- noSNV_permut_driver_marker_net_modules_flt[[i]]
  RNA_PPI_module_scores <- noSNV_permut_driver_marker_net_modules_flt_scores[[i]]
  RNA_PPI_first_five_modules = sapply(1, function(j) {
    tmp_module =
      RNA_PPI_modules[[names(RNA_PPI_module_scores)[j]]]
    
    tmp_module$gene
  }) %>% unlist()
  ###########
  
  tmp_intersect_tbl =
  intersect_combn(named_list = list("Original data" = original_first_five_modules,
                                    "Permuted mutations" = onlySNV_first_five_modules,
                                    "Permuted mutations+RNA-Seq" = SNV_RNA_first_five_modules,
                                    "Permuted mutations+RNA-Seq+PPIs" = SNV_RNA_PPI_first_five_modules,
                                    "Permuted RNA-Seq+PPIs" = RNA_PPI_first_five_modules))
  
  tmp_intersect_tbl
})

names(permut_modules_eval_list) <- c(names(driver_marker_net_modules_flt[-c(1,10,12)]), "Rhabdoid carcinoma")

permut_modules_eval_tbl <- do.call(rbind, permut_modules_eval_list)
permut_modules_eval_tbl$Cancer_type <- rep(names(permut_modules_eval_list), each = 10)

colnames(permut_modules_eval_tbl)[c(1,2)] <- paste("Origin_of_gene_set_", c(1,2), sep = "")

## save tables of results
openxlsx::write.xlsx(x = list("Top 20 genes" = permut_eval_tbl, "First ranked module" = permut_modules_eval_tbl), 
                     file = "Results/Evaluations/Permutation/permut_eval_tbl.xlsx")

## Further evaluate the results
for(i in 1:4) {
  
  cat(paste0(permut_eval_list$DMG$Vector2[i], ": ",
  mean(sapply(permut_eval_list, function(j) {
    j[i,3]
  })), sep = "\n"
  ))
}

##############

for(i in 1:4) {
  
  cat(paste0(permut_modules_eval_list$DMG$Vector2[i], ": ",
             mean(sapply(permut_modules_eval_list, function(j) {
               j[i,3]
             })), sep = "\n"
  ))
}

###########################

## Chord diagram visualization

par(mfrow = c(2,10))

# Visualize top 20s
for(i in 1:10) {
  
  set.seed(1234)
  chordDiagram(x = permut_eval_list[[i]],annotationTrackHeight = mm_h(c(1, 1)),
               annotationTrack = c("grid"),
               grid.col = setNames(c("#8F1D1EFF", "#02734AFF", "#FFC737FF", "#3E52A1FF", "#000000FF"), 
                                   union(permut_eval_list[[i]][,1], permut_eval_list[[i]][,2])),
               link.border = rep("black", 4))
  
  # circos.xaxis(h = 0.5, sector.index = "Original data", 
  #              labels = rev(permut_eval_list[[i]]$Intersection_length[1:4]),
  #              major.tick = F,
  #              minor.ticks = 0, 
  #              labels.cex = 1)
  
  title(names(permut_eval_list)[i])
  
  circos.clear()
  
}

##################

# Visualize the first ranked module
## alaki <- permut_modules_eval_list
### alaki <- alaki[c(1:6, 10, 7:9)] # and use alaki instead of permut_modules_eval_list as the last one has intersection but not the three networks before that

for (i in 1:10) {
  
  set.seed(1234)
  chordDiagram(x = permut_modules_eval_list[[i]],annotationTrackHeight = mm_h(c(1, 1)),
               annotationTrack = c("grid"),
               grid.col = setNames(c("#8F1D1EFF", "#02734AFF", "#FFC737FF", "#3E52A1FF", "#000000FF"), 
                                   union(permut_modules_eval_list[[i]][,1], permut_modules_eval_list[[i]][,2])),
               link.border = rep("black", 4))
  
  # circos.xaxis(h = 0.5,
  #              labels = permut_modules_eval_list[[i]]$Intersection_length[permut_modules_eval_list[[i]]$Intersection_length > 0],
  #              major.tick = F,
  #              minor.ticks = 0, 
  #              labels.cex = 1)
  
  circos.clear()
}

par(mfrow = c(1,1))

## draw the legend
library(ComplexHeatmap)

Legend(at = union(permut_modules_eval_list[[i]][,1], permut_modules_eval_list[[i]][,2]), 
       type = "lines", 
       legend_gp = gpar(col = c("#8F1D1EFF", "#02734AFF", "#FFC737FF", "#3E52A1FF", "#000000FF"), 
                        lwd = 6), 
       title_position = "topleft", 
       title = "Permutation level") %>% draw()

dev.off()

#=============================================================================
#
#    Code chunk 12: Disconnect from the database
#
#=============================================================================

## Disconnect from the database
dbDisconnect(db)


#=============================================================================
#
#    Code chunk 13: Testing the pipeline on ICGC adult data
#
#=============================================================================


#### Read in the sample data
icgc_metadata <- readr::read_delim("~/Downloads/pcawg_specimen_histology_August2016_v9.txt", delim = "\t")
icgc_metadata <- icgc_metadata %>% filter(grepl("Liver-HCC|Skin-Melanoma|Stomach-AdenoCA|Breast-AdenoCA|Myeloid-AML", histology_abbreviation))

####################################################################################################
####################################################################################################

#### Read in the genomic file
icgc_snv_con <- gzfile("~/Downloads/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz")
icgc_snv_data <- readr::read_delim(icgc_snv_con, delim = "\t")
rm(icgc_snv_con)

## Keep only the most common cancers with known specific types for model evaluation and the required columns and remove the intergenic region (IGR) variants
icgc_snv_data <- icgc_snv_data %>% select("Chromosome", "Start_position", "Reference_Allele", "Tumor_Seq_Allele2",
                                  "Hugo_Symbol", "Variant_Classification", "End_position", "Tumor_Sample_Barcode",
                                  "Genome_Change", "i_VAF", "t_alt_count", "t_ref_count", "Project_Code", "Donor_ID") %>% 
  filter(grepl("Liver-HCC|Skin-Melanoma|Stomach-AdenoCA|Breast-AdenoCa|Myeloid-AML", Project_Code)) %>% 
  filter(!grepl("IGR", Variant_Classification))

### Correct the project code Breast-AdenoCa to Breast-AdenoCA so that it matches the metadata
icgc_snv_data$Project_Code[which(icgc_snv_data$Project_Code == "Breast-AdenoCa")] <- "Breast-AdenoCA"

############################

# Remove InDels and MNVs

## Remove InDels
icgc_snv_data <- icgc_snv_data[-c(grep("-", icgc_snv_data$Reference_Allele), grep("-", icgc_snv_data$Tumor_Seq_Allele2)),]

## Remove MNVs
icgc_snv_data <- icgc_snv_data[-which(nchar(icgc_snv_data$Reference_Allele) == nchar(icgc_snv_data$Tumor_Seq_Allele2) & 
                                        nchar(icgc_snv_data$Tumor_Seq_Allele2) > 1),]

############################

# # Split the consecutive SNVs (SNVs in a row, e.g. 2:189704688-AGG-TTT)
# 
# ## First separate the multiple SNV rows
# icgc_snv_data_multiple_snv_index <- which(nchar(icgc_snv_data$Reference_Allele) == nchar(icgc_snv_data$Tumor_Seq_Allele2) & 
#                                             nchar(icgc_snv_data$Tumor_Seq_Allele2) > 1)
# 
# icgc_snv_data_multiple_snv <- icgc_snv_data[icgc_snv_data_multiple_snv_index,]
# icgc_snv_data <- icgc_snv_data[-icgc_snv_data_multiple_snv_index,]
# 
# ## Now split the rows of icgc_snv_data_multiple_snv and clean it
# icgc_snv_data_multiple_snv <- icgc_snv_data_multiple_snv %>% tidyr::separate_rows(Reference_Allele, Tumor_Seq_Allele2, sep = "")
# icgc_snv_data_multiple_snv <- icgc_snv_data_multiple_snv[-which(icgc_snv_data_multiple_snv$Reference_Allele == "" & icgc_snv_data_multiple_snv$Tumor_Seq_Allele2 == ""),]
# 
# ## Now merge two tables
# icgc_snv_data <- rbind(icgc_snv_data, icgc_snv_data_multiple_snv)
# rm(icgc_snv_data_multiple_snv)

############################

# Process downloaded gnomAD pop freq. data

### First, get min and max genomic positions of all chromosomes and create batches of 100,000 bps (genomic positions) for each chromosome
icgc_snv_data_chr <- 
  unique(icgc_snv_data$Chromosome)

icgc_snv_data_chr_tbl <-
  lapply(1:length(icgc_snv_data_chr), FUN = function(i) {
    
    ## Calculate the min and max of chromosome i
    icgc_snv_data_chr_min <- min(icgc_snv_data$Start_position[which(icgc_snv_data$Chromosome %in% icgc_snv_data_chr[i])]) - 10
    icgc_snv_data_chr_max <- max(icgc_snv_data$Start_position[which(icgc_snv_data$Chromosome %in% icgc_snv_data_chr[i])]) + 10
    
    icgc_snv_data_chr_tbl_tmp <- data.frame(Chromosome = icgc_snv_data_chr[i], 
                                            min = icgc_snv_data_chr_min, 
                                            max = icgc_snv_data_chr_max)
    
  })

icgc_snv_data_chr_tbl <-  do.call(rbind, icgc_snv_data_chr_tbl)

############################

## Get the header of VCF file
genomes.r2.1.1.exome_calling_header %>% class() # already imported

## Filter out chromosomes of icgc_snv_data_chr_tbl that are not available in the genomes.r2.1.1.exome_calling_header
genomes.r2.1.1.exome_calling_chrs <- GenomeInfoDb::seqlevels(genomes.r2.1.1.exome_calling_header)
icgc_snv_data_chr_gnomAD_tbl <- icgc_snv_data_chr_tbl[which(icgc_snv_data_chr_tbl$Chromosome %in% genomes.r2.1.1.exome_calling_chrs),]

## Remove Y chromosome as it is not present in the Tabix index file
icgc_snv_data_chr_gnomAD_tbl <- icgc_snv_data_chr_gnomAD_tbl[-which(icgc_snv_data_chr_gnomAD_tbl$Chromosome == "Y"),]

icgc_snv_data_popFreq_list <- parallel::mclapply(1:nrow(icgc_snv_data_chr_gnomAD_tbl), function(i) {
  
  ## Define a genomic ranges
  chr = icgc_snv_data_chr_gnomAD_tbl[i,1]
  start.loc <- icgc_snv_data_chr_gnomAD_tbl[i,2]
  end.loc <- icgc_snv_data_chr_gnomAD_tbl[i,3]
  gr <- GenomicRanges::GRanges(chr, IRanges::IRanges(start.loc, end.loc))
  
  ## restrict VCF INFO columns to AC and AN values
  vcfPar <- VariantAnnotation::ScanVcfParam(geno = NA,
                                            fixed = c("ALT", "FILTER"),
                                            info = AFcols, # AFcols was previously prepared based on rownames(VariantAnnotation::info(genomes.r2.1.1.exome_calling_header))
                                            which = gr)
  
  # Load the desired chunk of VCF file
  vcf <- VariantAnnotation::readVcf(
    Rsamtools::TabixFile("Datasets/gnomAD/V2/gnomad.genomes.r2.1.1.exome_calling_intervals.sites.vcf.bgz"),
    param=vcfPar)
  
  # Discard variants not passing all FILTERS
  mask <- VariantAnnotation::filt(vcf) == "PASS"
  vcf <- vcf[mask, ]
  
  # prepare the variants freq data
  ref <- VariantAnnotation::ref(vcf) %>% as.character()
  alt <- VariantAnnotation::alt(vcf) %>% unlist() %>% as.character()
  pos <- MatrixGenerics::rowRanges(vcf) %>% BiocGenerics::start()
  af <- VariantAnnotation::info(vcf) %>% as.data.frame()
  
  # Generate the gnomAD pop. freq. table
  gnomAD_af_tbl <- cbind(data.frame(chr = chr, pos = pos, ref = ref, alt = alt), af)
  
  gnomAD_af_tbl
  
}, mc.cores = parallel::detectCores() - 4)

icgc_snv_data_popFreq_list <- do.call(rbind, icgc_snv_data_popFreq_list)

############################

# Correct the indels of gnomAD pop freq data as ICGC indel data are dashed

## Correct insertions
icgc_snv_data_popFreq_list_insertions_index <- which(nchar(icgc_snv_data_popFreq_list$ref) < nchar(icgc_snv_data_popFreq_list$alt))

if(length(icgc_snv_data_popFreq_list_insertions_index) > 0) {
  
  icgc_snv_data_popFreq_list$ref[icgc_snv_data_popFreq_list_insertions_index] <- "-"
  icgc_snv_data_popFreq_list$alt[icgc_snv_data_popFreq_list_insertions_index] <- substr(x = icgc_snv_data_popFreq_list$alt[icgc_snv_data_popFreq_list_insertions_index], 
                                                                                        start = 2, stop = nchar(icgc_snv_data_popFreq_list$alt[icgc_snv_data_popFreq_list_insertions_index]))
}

## Correct deletions
icgc_snv_data_popFreq_list_deletions_index <- which(nchar(icgc_snv_data_popFreq_list$ref) > nchar(icgc_snv_data_popFreq_list$alt))

if(length(icgc_snv_data_popFreq_list_deletions_index) > 0) {
  
  icgc_snv_data_popFreq_list$pos[icgc_snv_data_popFreq_list_deletions_index] <- icgc_snv_data_popFreq_list$pos[icgc_snv_data_popFreq_list_deletions_index] + 1
  icgc_snv_data_popFreq_list$ref[icgc_snv_data_popFreq_list_deletions_index] <- substr(x = icgc_snv_data_popFreq_list$ref[icgc_snv_data_popFreq_list_deletions_index], 
                                                                                       start = 2, stop = nchar(icgc_snv_data_popFreq_list$ref[icgc_snv_data_popFreq_list_deletions_index]))
  icgc_snv_data_popFreq_list$alt[icgc_snv_data_popFreq_list_deletions_index] <- "-"
}

############################

## Filter out high gnomAD pop freq vars
  
  ### identify gnomAD vars that are present in icgc_snv_data
  icgc_snv_data_gnomAD_pop_freq_index <- which(paste(icgc_snv_data_popFreq_list$chr, 
                                icgc_snv_data_popFreq_list$pos, 
                                icgc_snv_data_popFreq_list$ref,
                                "_",
                                icgc_snv_data_popFreq_list$alt, sep = "") %in% 
                            paste(icgc_snv_data$Chromosome, 
                                  icgc_snv_data$Start_position, 
                                  icgc_snv_data$Reference_Allele,
                                  "_",
                                  icgc_snv_data$Tumor_Seq_Allele2, sep = ""))
  
  ### identify icgc_snv_data_gnomAD_pop_freq_index which have high frequency (FC > 0.01)
  icgc_snv_data_gnomAD_pop_freq_high_freq_index <- (((icgc_snv_data_popFreq_list[icgc_snv_data_gnomAD_pop_freq_index,c(5:ncol(icgc_snv_data_popFreq_list))] > 0.01) %>% 
                                  rowSums(na.rm = TRUE)) > 0)
  
  ### Filter out rare vars of icgc_snv_data_gnomAD_pop_freq_index
  icgc_snv_data_gnomAD_pop_freq_index <- icgc_snv_data_gnomAD_pop_freq_index[icgc_snv_data_gnomAD_pop_freq_high_freq_index]
  
  ### identify index of high freq vars in icgc_snv_data
  icgc_snv_data_gnomAD_index <- match(paste(icgc_snv_data_popFreq_list$chr, 
                                            icgc_snv_data_popFreq_list$pos, 
                                            icgc_snv_data_popFreq_list$ref,
                                            "_",
                                            icgc_snv_data_popFreq_list$alt, sep = "")[icgc_snv_data_gnomAD_pop_freq_index],
                                      paste(icgc_snv_data$Chromosome, 
                                            icgc_snv_data$Start_position, 
                                            icgc_snv_data$Reference_Allele,
                                            "_",
                                            icgc_snv_data$Tumor_Seq_Allele2, sep = ""))


# filter the icgc_snv_data
icgc_snv_data_flt <- icgc_snv_data[-icgc_snv_data_gnomAD_index,]

# remove icgc_snv_data and icgc_snv_data_popFreq_list to free up space
rm(icgc_snv_data_popFreq_list, icgc_snv_data)

#####################################################

## Add CADD scores of SNVs

## Prepare icgc_snv_data_chr_cadd_tbl table
genomic_intervals <- 50000

icgc_snv_data_chr_cadd_tbl <- 
  lapply(1:length(icgc_snv_data_chr), FUN = function(i) {
    
    ## Calculate the min and max of chromosome i
    icgc_snv_data_chr_min <- min(icgc_snv_data_flt$Start_position[which(icgc_snv_data_flt$Chromosome %in% icgc_snv_data_chr[i])])
    icgc_snv_data_chr_max <- max(icgc_snv_data_flt$Start_position[which(icgc_snv_data_flt$Chromosome %in% icgc_snv_data_chr[i])])
    
    ## Define the intervals for chromosome i
    if((icgc_snv_data_chr_max - icgc_snv_data_chr_min) <= genomic_intervals) {
      icgc_snv_data_chr_cadd_tbl_tmp <- data.frame(min = icgc_snv_data_chr_min, max = icgc_snv_data_chr_max)
    } else {
      icgc_snv_data_chr_cadd_tbl_tmp <- 
        interval_gen(min = icgc_snv_data_chr_min, 
                     max = icgc_snv_data_chr_max, 
                     size = genomic_intervals)
    }
    
    ## Add chromosome and refs columns to the table
    icgc_snv_data_chr_cadd_tbl_tmp <- cbind(chr = icgc_snv_data_chr[i],
                                            icgc_snv_data_chr_cadd_tbl_tmp,
                                            refs = vector(mode = "logical", length = nrow(icgc_snv_data_chr_cadd_tbl_tmp)))
    
    ## Define the refs for each interval
    sapply(as.list(1:nrow(icgc_snv_data_chr_cadd_tbl_tmp)), function(j) {
      icgc_snv_data_chr_cadd_tbl_tmp[j, 4] <<- ifelse(any(icgc_snv_data_flt$Chromosome %in% icgc_snv_data_chr[i] &
                                                            icgc_snv_data_flt$Start_position >= icgc_snv_data_chr_cadd_tbl_tmp[j, 2] &
                                                            icgc_snv_data_flt$Start_position < icgc_snv_data_chr_cadd_tbl_tmp[j, 3]),
                                                      TRUE, FALSE)
      
      
    })
    
    icgc_snv_data_chr_cadd_tbl_tmp$refs <- as.logical(icgc_snv_data_chr_cadd_tbl_tmp$refs)
    
    icgc_snv_data_chr_cadd_tbl_tmp
    
  })

## Merge the list of tables
icgc_snv_data_chr_cadd_tbl <- do.call(rbind, icgc_snv_data_chr_cadd_tbl)

## Filter out genomic regions that have no variants in them
icgc_snv_data_chr_cadd_tbl <- icgc_snv_data_chr_cadd_tbl[icgc_snv_data_chr_cadd_tbl$refs,]

icgc_snv_data_chr_cadd_tbl$refs <- NULL

## Retrieve CADD scores for all chromosomes of icgc_snv_data (each interval of 50,000 last ~ 0.7-0.8 sec to be retrieved)

icgc_snv_data_cadd_phred_scores <- scanCADD2df(variants_table = icgc_snv_data_flt,
                                               tabix_snv_file = "Datasets/CADD/v1.6-hg19 (compressed)/whole_genome_SNVs.compressed.tsv.gz",
                                               tabix_indel_file = NULL, # The tabix indexed file of InDels (located in the same directory that .tbi file is located)
                                               gr_table = icgc_snv_data_chr_cadd_tbl)

icgc_snv_data_cadd_phred_scores <- do.call(rbind, icgc_snv_data_cadd_phred_scores)

# Checking the accuracy of the data (the results are correct, no idea why some CADD scores are missing!)
sample_icgc_snv_data_cadd_phred_scores <- icgc_snv_data_cadd_phred_scores[sample(seq(nrow(icgc_snv_data_cadd_phred_scores)), size = 30),]
View(icgc_snv_data_flt[sample_icgc_snv_data_cadd_phred_scores$index,])

# Retrieve missing CADD scores
icgc_snv_data_cadd_phred_scores_missing_indices <- seq(nrow(icgc_snv_data_flt))[-which(seq(nrow(icgc_snv_data_flt)) %in% icgc_snv_data_cadd_phred_scores$index)]

icgc_snv_data_flt_missed <- icgc_snv_data_flt[icgc_snv_data_cadd_phred_scores_missing_indices,]

icgc_snv_data_chr_cadd_tbl_missed <- icgc_snv_data_flt_missed[,c(1,2)]
colnames(icgc_snv_data_chr_cadd_tbl_missed) <- c("chr", "min")
icgc_snv_data_chr_cadd_tbl_missed$max <- icgc_snv_data_chr_cadd_tbl_missed$min

icgc_snv_data_cadd_phred_scores_missed <- scanCADD2df(variants_table = icgc_snv_data_flt_missed,
                                               tabix_snv_file = "Datasets/CADD/v1.6-hg19 (compressed)/whole_genome_SNVs.compressed.tsv.gz",
                                               tabix_indel_file = NULL, # The tabix indexed file of InDels (located in the same directory that .tbi file is located)
                                               gr_table = icgc_snv_data_chr_cadd_tbl_missed)

icgc_snv_data_cadd_phred_scores_missed <- do.call(rbind, icgc_snv_data_cadd_phred_scores_missed)

# Add CADD scores to the table
icgc_snv_data_flt$cadd_score[icgc_snv_data_cadd_phred_scores$index] <- icgc_snv_data_cadd_phred_scores$score
icgc_snv_data_flt$cadd_score[icgc_snv_data_cadd_phred_scores_missing_indices] <- icgc_snv_data_cadd_phred_scores_missed$score

icgc_snv_data_flt_withCADD <- icgc_snv_data_flt
rm(icgc_snv_data_flt)

## Change CADD score colname to match future ClinVar data for callibration
colnames(icgc_snv_data_flt_withCADD)[ncol(icgc_snv_data_flt_withCADD)] <- "cadd_phred"

#####################################

# Calibrating the CADD scores

library(mgcv)

### We can either use the type="link" to get (-Inf,Inf) predictions then use logistic regression to transform it to (0,1) probabilities or use type="response" to directly get (0,1) probabilities.
#### with type="response" predictions on the scale of the response are returned. So, type="link" would be redundant.
# icgc_snv_data_flt_withCADD_gam_prob <- 1/(1+exp(-predict(clinVar_summary_for_training_gam, 
#                                                            newdata = icgc_snv_data_flt_withCADD[,"cadd_phred"], 
#                                                            type="link")))

### get the gam probabilities using the logistic transformation
icgc_snv_data_flt_withCADD_gam_prob <- predict(clinVar_summary_for_training_gam, 
                                               newdata = icgc_snv_data_flt_withCADD[,"cadd_phred"], 
                                               type="response")

### Visualize calibration probabilities vs original cadd
icgc_snv_data_flt_withCADD$calibCADD_GAM_prob <- as.numeric(icgc_snv_data_flt_withCADD_gam_prob)

icgc_snv_data_flt_GAM_probvsPhred_plot <-
  ggplot(icgc_snv_data_flt_withCADD, aes(x = cadd_phred, y = calibCADD_GAM_prob)) +
  scale_color_gradientn(colours = c("blue", "red")) + 
  # guides(color = guide_legend(title = "Original CADD Phred")) +
  geom_point(aes(color = calibCADD_GAM_prob), show.legend = F) +
  labs(x = "CADD Phred Score", y= "CADD GAM Probability") +
  guides(fill=guide_legend(ncol=1)) +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

icgc_snv_data_flt_GAM_probvsPhred_plot

ggsave(filename = "Results/ICGC/GAM-based calibration/icgc_snv_data_flt_GAM_probvsPhred_plot.png", 
       plot = icgc_snv_data_flt_GAM_probvsPhred_plot, 
       device = "png", dpi = 300, width = 8, height = 8, units = "in")

# Log transform GAM probabilities
icgc_snv_data_flt_withCADD$calibCADD_GAM_prob_log <- (-10)*log10(1-icgc_snv_data_flt_withCADD$calibCADD_GAM_prob)

#######################

# removing ribosomal protein, HLA, and olfactory receptor genes from our data
icgc_snv_data_flt_withCADD <- icgc_snv_data_flt_withCADD[-which(stringr::str_detect(string = icgc_snv_data_flt_withCADD$Hugo_Symbol, 
                                                                                    pattern = paste(ribo_p_pattern, or_pattern, hla_pattern, sep = "|"))),]

#######################

## Summarize the data by summing the calibrated CADD scores of each gene within each sample

### Summarize
icgc_snv_data_flt_withCADD_summarized <- icgc_snv_data_flt_withCADD %>% 
  group_by(Chromosome, Hugo_Symbol, Project_Code, Tumor_Sample_Barcode) %>% 
  dplyr::summarise(cadd_phred = sum(calibCADD_GAM_prob))

######################

### Re-transform the Phred/log scaled pathogenicities to linear scale probabilities (Not required as we used calibCADD_GAM_prob instead of calibCADD_GAM_prob_log as the cadd_phred)
# icgc_snv_data_flt_withCADD_summarized$cadd_phred <- 1-10^(icgc_snv_data_flt_withCADD_summarized$cadd_phred/-10)

######################

## We should create a list of SNV data with different cancer types and perform subsequent analyses separately for each cancer type
icgc_snv_data_flt_withCADD_summarized <- lapply(unique(icgc_snv_data_flt_withCADD_summarized$Project_Code), function(i) {
  icgc_snv_data_flt_withCADD_summarized[which(icgc_snv_data_flt_withCADD_summarized$Project_Code == i),]
})

### Set the names of list
icgc_snv_data_flt_cancer_types <- sapply(icgc_snv_data_flt_withCADD_summarized, function(i) unique(i$Project_Code))
names(icgc_snv_data_flt_withCADD_summarized) <- icgc_snv_data_flt_cancer_types

## Prepare tables in the form of genes on columns and samples on rows for corr analysis
icgc_snv_data_flt_withCADD_summarized <- lapply(icgc_snv_data_flt_withCADD_summarized, function(i) {
  tmp.data <- i[, c(2,4,5)]
  tmp.data <- tidyr::pivot_wider(data = tmp.data, 
                                 id_cols = "Tumor_Sample_Barcode", 
                                 names_from = "Hugo_Symbol", 
                                 values_from = "cadd_phred",
                                 values_fn = mean) %>% as.data.frame()
  
  
  rownames(tmp.data) <- tmp.data$sample_id
  tmp.data <- tmp.data[,-1]
  
  # Replace NA with 0
  tmp.data <- apply(tmp.data, 1, function(i) {i[is.na(i)]  <- 0; i}) %>% t() %>% as.data.frame()
  
  # Report the samples that their number of corresponding genes with data is less than < 1 percent of the total number of genes in the table of its corresponding cancer type
  row_index <- apply(tmp.data, 1, function(i) {sum(i > 0)})/ncol(tmp.data) < 0.01
  if(sum(row_index) > 0) {
    tmp_text <- paste("\nThe following samples have mutations in very few number of genes. You may further check their quality.\n",
                      paste0(rownames(tmp.data)[row_index], collapse = "\n"), sep = "")
    cat(tmp_text, sep = "")
  }
  
  # Remove noise genes (very low probability of being pathogenic in all samples)
  print("I changed tmp.data > 0.05 to tmp.data > 0.01 to be consistent with the rest of the project")
  tmp.data <- tmp.data[,colSums(tmp.data > 0.01) != 0]
  
  tmp.data
})

## Get the number of genes in each sample in each dataset
icgc_snv_data_flt_withCADD_geneCount <- lapply(1:length(icgc_snv_data_flt_withCADD_summarized), function(i) {
  
  tmp_tbl <- icgc_snv_data_flt_withCADD_summarized[[i]]
  
  tmp_geneCount <- apply(tmp_tbl, 1, function(j) {sum(j > 0)}) %>% as.data.frame()
  tmp_geneCount <- cbind(cancer_type = icgc_snv_data_flt_cancer_types[i],
                         sample_id = rownames(tmp_geneCount),
                         gene_count = tmp_geneCount[,1]) %>% as.data.frame()
  rownames(tmp_geneCount) <- NULL
  
  tmp_geneCount
})

icgc_snv_data_flt_withCADD_geneCount <- do.call(rbind, icgc_snv_data_flt_withCADD_geneCount)

# icgc_snv_data_flt_withCADD_summarized <- lapply(icgc_snv_data_flt_withCADD_summarized, function(i) {
#   tmp_tbl <- i
#   tmp_tbl <- tmp_tbl[,colSums(tmp_tbl > 0.05) != 0]
#   tmp_tbl
#   }
#   )

######################

## Correlation analysis

icgc_snv_data_flt_cor <- lapply(icgc_snv_data_flt_withCADD_summarized, function(i){
  tmp_cor <- fcor(i)
  tmp_cor <- subset(tmp_cor, mr < 20)
  tmp_cor
})

####################################################################################################
####################################################################################################

#### Read in the transcriptomic data 
icgc_tx_con <- gzfile("~/Downloads/tophat_star_fpkm.v2_aliquot_gl.tsv.gz")
icgc_tx_data <- readr::read_delim(icgc_tx_con, delim = "\t")
rm(icgc_tx_con)

icgc_tx_data <- as.data.frame(icgc_tx_data)
rownames(icgc_tx_data) <- icgc_tx_data$feature

## Filter out the samples corresponding to other cancer types
icgc_tx_data <- icgc_tx_data[,which(colnames(icgc_tx_data) %in% icgc_metadata$tcga_sample_uuid)]

#######################

# Convert the Ensembl gene names to gene symbols

## Create a dictionary of gene name and symbols
icgc_tx_ens_genes <- rownames(icgc_tx_data)
icgc_tx_ens_genes_noVersion <- gsub(pattern = "\\..*", replacement = "", x = icgc_tx_ens_genes)

icgc_tx_gene_symbols <- gprofiler2::gconvert(query = icgc_tx_ens_genes_noVersion, organism = "hsapiens",
                                             target = "HGNC", mthreshold = 1, filter_na = F)

icgc_tx_gene_dict <- data.frame(ens = icgc_tx_ens_genes, hgnc = icgc_tx_gene_symbols$target)
which(is.na(icgc_tx_gene_dict$hgnc)) %>% length()

## Remove the genes without a name
icgc_tx_data <- icgc_tx_data[-which(is.na(icgc_tx_gene_dict$hgnc)),]
icgc_tx_gene_dict <- icgc_tx_gene_dict[-which(is.na(icgc_tx_gene_dict$hgnc)),]

## Check the duplicate gene names corresponding to different isoforms
length(which(duplicated(icgc_tx_gene_dict$hgnc))) # 8

### add the gene symbols to the main df
icgc_tx_data$symbol <- icgc_tx_gene_dict$hgnc

### remove the duplicates and average their expressions across all samples for each duplicated gene
icgc_tx_data <- icgc_tx_data %>% 
  dplyr::group_by(symbol) %>% 
  summarize(across(everything(), mean))

length(unique(icgc_tx_data$symbol)) == nrow(icgc_tx_data)
icgc_tx_data <- as.data.frame(icgc_tx_data)
rownames(icgc_tx_data) <- icgc_tx_data$symbol
icgc_tx_data$symbol <- NULL

#######################

# # removing ribosomal protein, HLA, and olfactory receptor genes from our data
icgc_tx_data <- icgc_tx_data[-which(stringr::str_detect(string = rownames(icgc_tx_data), 
                                                        pattern = paste(ribo_p_pattern, or_pattern, hla_pattern, sep = "|"))),]

#######################

## Add the cancer type of each sample
icgc_tx_data <- as.data.frame(t(icgc_tx_data))
icgc_tx_data$Project_Code <- icgc_metadata$histology_abbreviation[match(rownames(icgc_tx_data), icgc_metadata$tcga_sample_uuid)]

## Separate the tx of different cancer types
icgc_tx_data_Ca_names <- unique(icgc_tx_data$Project_Code)
icgc_tx_data <- lapply(icgc_tx_data_Ca_names, function(i) {
  tmp_data <- subset(icgc_tx_data, Project_Code == i)
  tmp_data <- tmp_data[,-ncol(tmp_data)]
  return(tmp_data)
})

names(icgc_tx_data) <- icgc_tx_data_Ca_names

## Filter out noise genes from each data set
# Remove noise genes
icgc_tx_data <- lapply(icgc_tx_data, function(i) {
  keep_tmp <- selectGenes(tbl = i[-nrow(i),], gene_axis = "column",
                          min.abs = 1, # At least FPKM of 1
                          min.sample = 0.1) # In at least 10 percent of the samples
  table_tmp <- i[,keep_tmp]
  table_tmp
})

# Perform correlation analysis
icgc_tx_cor_list <- lapply(icgc_tx_data, function(i) {
  tmp_data <- fcor(i)
  tmp_data <- subset(tmp_data, mr < 20)
  tmp_data
})

####################################################################################################
####################################################################################################
####################################################################################################

# Reconstruction and analysis of the Driver Net

## Unify the order of cancers in icgc_snv_data_flt_cor and icgc_tx_cor_list
icgc_tx_cor_list <- icgc_tx_cor_list[names(icgc_snv_data_flt_cor)]

# Prepare the tables for merging with other data
icgc_snv_data_flt_cor <- lapply(icgc_snv_data_flt_cor, function(i) {
  tmp_tbl <- i[,c(1,2)]
  colnames(tmp_tbl) <- c("from", "to")
  tmp_tbl$type <- "co-deletriousness"
  tmp_tbl
})

icgc_tx_cor_list <- lapply(icgc_tx_cor_list, function(i) {
  tmp_tbl <- i[,c(1,2)]
  colnames(tmp_tbl) <- c("from", "to")
  tmp_tbl$type <- "co-expression"
  tmp_tbl
})

# Merge the similarity tables
icgc_driver_net <- icgc_snv_data_flt_cor

# Add the co-expression data to the tables

for(i in 1:length(icgc_driver_net)) {
  
  genes <- icgc_driver_net[[i]][,c(1,2)] %>% 
    unlist() %>% 
    unname() %>% 
    unique()
  
  first_level_coex_index <- c(which(icgc_tx_cor_list[[i]]$from %in% genes), 
                              which(icgc_tx_cor_list[[i]]$to %in% genes)) %>% unique()
  
  first_level_coex_genes <- c(genes, 
                              icgc_tx_cor_list[[i]]$from[first_level_coex_index],
                              icgc_tx_cor_list[[i]]$to[first_level_coex_index]) %>% unique()
  
  first_two_levels_coex_index <- c(which(icgc_tx_cor_list[[i]]$from %in% first_level_coex_genes), 
                                   which(icgc_tx_cor_list[[i]]$to %in% first_level_coex_genes)) %>% unique()
  
  icgc_driver_net[[i]] <- rbind(icgc_driver_net[[i]], icgc_tx_cor_list[[i]][first_two_levels_coex_index,])
  
}

##################

# Add PPI data to the driver/marker net

# Map gene names to stringDB ids for PPI analysis

## get the gene names in each risk net
icgc_driver_net_genes <- lapply(icgc_driver_net, function(i) {
  data.frame(gene_symbol = i[,c(1,2)] %>% 
               unlist() %>% 
               unname() %>% 
               unique())
})

## getSTRINGdb for human
library(STRINGdb)
string_db <- STRINGdb$new(species=9606, version = "11.5", 
                          score_threshold=900, # very high confidence (>0.9)
                          input_directory=""
)

icgc_driver_net_genes <- lapply(icgc_driver_net_genes, function(i) {
  
  string_db$map(my_data_frame = i, 
                my_data_frame_id_col_names = "gene_symbol", 
                takeFirst = TRUE, removeUnmappedRows = TRUE)
  
})

# Get the PPIs
icgc_driver_net_ppi <- lapply(icgc_driver_net_genes, function(i) {
  
  string_db$get_interactions(string_ids = i$STRING_id)
  
})

# Convert stringdb ids to gene symbols
icgc_driver_net_ppi <- lapply(1:length(icgc_driver_net_ppi), function(i) {
  tmp_row_index <- match(icgc_driver_net_ppi[[i]]$from, icgc_driver_net_genes[[i]]$STRING_id)
  icgc_driver_net_ppi[[i]]$from <- icgc_driver_net_genes[[i]]$gene_symbol[tmp_row_index]
  
  tmp_column_index <- match(icgc_driver_net_ppi[[i]]$to, icgc_driver_net_genes[[i]]$STRING_id)
  icgc_driver_net_ppi[[i]]$to <- icgc_driver_net_genes[[i]]$gene_symbol[tmp_column_index]
  
  tmp_tbl <- icgc_driver_net_ppi[[i]][,c(1,2)]
  tmp_tbl$type <- "PPI"
  
  ## Remove duplicates
  tmp_tbl <- tmp_tbl[-which(duplicated(paste(tmp_tbl$from,
                                             tmp_tbl$to, 
                                             sep = "_"))),]
  
  return(tmp_tbl)
})

names(icgc_driver_net_ppi) <- names(icgc_driver_net_genes)

# Combine PPI with the original driver/marker net
icgc_driver_net <- lapply(1:length(icgc_driver_net), function(i) {
  rbind(icgc_driver_net[[i]], icgc_driver_net_ppi[[i]])
})

names(icgc_driver_net) <- names(icgc_driver_net_genes)

# Reconstruct the networks
icgc_driver_net <- lapply(icgc_driver_net, function(i) {
  igraph::graph_from_data_frame(i, directed=FALSE)
})

## Check the network layers
table(igraph::E(icgc_driver_net$`Breast-AdenoCA`)$type)

###############################

# Calculate the IVI values of genes within each net
library(influential)
icgc_driver_net_ivi <- lapply(icgc_driver_net, ivi)

###############################

# Calculate the primitive gene driver scores
icgc_driver_net_gene_tbl <- lapply(icgc_driver_net_ivi, function(i) {
  data.frame(gene = names(i), ivi = i)
})

## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer


icgc_snv_mean_deletriousness <- lapply(icgc_snv_data_flt_withCADD_summarized, function(i) {
  colMeans(i)
})

## Generate the gene table for each net and calculate the primitive scores
icgc_driver_net_gene_tbl <- lapply(1:length(icgc_driver_net_gene_tbl), function(i) {
  tmp_tbl <- icgc_driver_net_gene_tbl[[i]]
  tmp_tbl$mean_deletriousness <- 0
  match_index <- which(names(icgc_snv_mean_deletriousness[[i]]) %in% rownames(tmp_tbl))
  tmp_tbl[names(icgc_snv_mean_deletriousness[[i]])[match_index], "mean_deletriousness"] <- 
    icgc_snv_mean_deletriousness[[i]][match_index]
  
  tmp_tbl$primitive_driver_score <- (tmp_tbl$ivi)*((tmp_tbl$mean_deletriousness)) # gene scores
  
  return(tmp_tbl)
  
})

# Add node scores
icgc_driver_net <- lapply(1:length(icgc_driver_net), function(i) {
  set_vertex_attr(graph = icgc_driver_net[[i]], name = "score", 
                  value = (icgc_driver_net_gene_tbl[[i]]$ivi)*(icgc_driver_net_gene_tbl[[i]]$mean_deletriousness))
})

names(icgc_driver_net) <- names(icgc_snv_mean_deletriousness)

## Check the network node scores
summary(igraph::V(icgc_driver_net$`Breast-AdenoCA`)$score)

names(icgc_driver_net_gene_tbl) <- names(icgc_snv_mean_deletriousness)

# Convert the un-weighted to a node-weighted network
icgc_driver_net <- lapply(1:length(icgc_driver_net), function(i) {
  set_vertex_attr(graph = icgc_driver_net[[i]], name = "weight", 
                  value = (icgc_driver_net_gene_tbl[[i]]$ivi)^(icgc_driver_net_gene_tbl[[i]]$mean_deletriousness))
})

names(icgc_driver_net) <- names(icgc_snv_mean_deletriousness)

## Check the network node weights
summary(igraph::V(icgc_driver_net$`Breast-AdenoCA`)$weight)

# Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
icgc_driver_net_gene_tbl <- lapply(1:length(icgc_driver_net), function(i) {
  tmp_tbl <- icgc_driver_net_gene_tbl[[i]]
  tmp_tbl$mean_neighborhood_weight <- sapply(1:nrow(tmp_tbl), function(j) {
    igraph::neighbors(graph = icgc_driver_net[[i]], 
                      v = V(icgc_driver_net[[i]])[j])$score %>% mean()
  })
  
  tmp_tbl$mediator_score <- tmp_tbl$primitive_driver_score * tmp_tbl$mean_neighborhood_weight
  
  return(tmp_tbl)
})

names(icgc_driver_net_gene_tbl) <- names(icgc_snv_mean_deletriousness)

icgc_driver_net_gene_tbl$`Breast-AdenoCA` %>% 
  slice_max(order_by = primitive_driver_score, n = 20) %>% 
  mutate(Rank = rank(-(primitive_driver_score))) %>% 
  select(Rank, gene) %>% 
  rename(Gene = gene) %>% 
  View()

# Check the top 10 genes of breast cancer
icgc_driver_net_gene_tbl$`Breast-AdenoCA` %>% slice_max(order_by = primitive_driver_score, n = 10) %>% select(gene) %>% unlist() %>% unname()

###################################################

# Detect communities/modules/clusters
icgc_driver_net_modules <- lapply(icgc_driver_net, function(i) {
  set.seed(3847)
  igraph::cluster_leiden(
    graph = i,
    objective_function = "CPM",
    weights = NULL,
    resolution_parameter = 0.5,
    beta = 0.05,
    initial_membership = NULL,
    n_iterations = 100,
    vertex_weights = igraph::V(i)$weight
  )
})

sizes(icgc_driver_net_modules$`Breast-AdenoCA`) %>% as.vector() %>% summary()
length(icgc_driver_net_modules$`Breast-AdenoCA`)

## Visualize the entire network and color by modules (NOT Recommended; It will take a while and will be too crowded)
# plot(
#   x = icgc_driver_net_modules$`CNS other`, 
#   y = icgc_driver_net$`CNS other`,
#   col = membership(icgc_driver_net_modules$`CNS other`),
#   vertex.label=NA,
#   mark.groups = communities(icgc_driver_net_modules$`CNS other`),
#   edge.color = c("black", "red")[crossing(icgc_driver_net_modules$`CNS other`, icgc_driver_net$`CNS other`) + 1]
# )

## Inspect and filter modules
icgc_driver_net_modules_flt <- lapply(1:length(icgc_driver_net_modules), function(i) {
  
  ### Create a table of modules and their genes
  tmp_tbl <- data.frame(
    module = membership(icgc_driver_net_modules[[i]]) %>% as.integer(),
    gene = icgc_driver_net_modules[[i]]$names
  )
  
  ### Add if the gene is mutated or not (could be a driver gene or not)
  tmp_tbl$mutated <- 0
  tmp_tbl$mutated[which(icgc_driver_net_gene_tbl[[i]]$mean_deletriousness > 0)] <- 1
  
  ### Separate modules
  tmp_tbl <- lapply(unique(tmp_tbl$module), function(j){
    subset(tmp_tbl, module == j)
  })
  
  names(tmp_tbl) <- sapply(tmp_tbl, function(m) {
    unique(m$module)
  })
  
  ### Remove modules that have less than 2 driver genes or their size is less than 4
  tmp_tbl <- lapply(tmp_tbl, function(k) {
    if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
      k
    } else {
      NULL
    }
  })
  
  ### Remove NULL modules
  tmp_tbl[which(sapply(tmp_tbl, is.null))] <- NULL
  
  ### Calculate the module scores
  tmp_tbl <- lapply(tmp_tbl, function(l) {
    tmp_tbl4score <- cbind(l, mean_score = 0)
    # tmp_tbl4score$mean_score <- mean(V(icgc_driver_net[[i]])$score[which(as_ids(V(icgc_driver_net[[i]])) %in% l$gene)])
    tmp_tbl4score$mean_score <- mean(icgc_driver_net_gene_tbl[[i]]$primitive_driver_score[which(icgc_driver_net_gene_tbl[[i]]$gene %in% l$gene)])
    tmp_tbl4score
  })
  
  tmp_tbl
})

names(icgc_driver_net_modules_flt) <- names(icgc_driver_net_modules)
sapply(icgc_driver_net_modules_flt, length)

beep(8)

################

## Sort the mean scores of all modules of each cancer type
icgc_driver_net_modules_flt_scores <-
  lapply(icgc_driver_net_modules_flt, function(i) {
    tmp_tbl <-
      sapply(i, function(j) {unique(j$mean_score)})
    tmp_tbl <- rev(sort(tmp_tbl))
    tmp_tbl
  })

icgc_driver_net_modules_flt_scores$`Breast-AdenoCA` %>% head()

### View the number of first ranked modules
sapply(icgc_driver_net_modules_flt_scores, function(i) {i[1]}) %>% as.data.frame() %>% View()

####################################

### Evaluate the breast cancer results

### Top 10 modules

### View the genes of top 10 modules of breast cancer
icgc_top10_driver_modules <- data.frame(Module_Rank = c(1:10))
icgc_top10_driver_modules$Genes <- lapply(icgc_driver_net_modules_flt$`Breast-AdenoCA`[match(names(icgc_driver_net_modules_flt_scores$`Breast-AdenoCA`[1:10]),
                                                                                             names(icgc_driver_net_modules_flt$`Breast-AdenoCA`))], 
                                          function(i) i$gene)

### Import the drivers of BrCa
icgc_top10_driver_modules <- mutate(icgc_top10_driver_modules, Genes = sapply(Genes, toString))

####*********####


### Top 5 modules

### View the genes of top 5 modules of breast cancer
icgc_top5_driver_modules <- data.frame(Module_Rank = c(1:5))
icgc_top5_driver_modules$Genes <- lapply(icgc_driver_net_modules_flt$`Breast-AdenoCA`[match(names(icgc_driver_net_modules_flt_scores$`Breast-AdenoCA`[1:5]),
                                                                                            names(icgc_driver_net_modules_flt$`Breast-AdenoCA`))], 
                                         function(i) i$gene)

### Import the drivers of BrCa
icgc_top5_driver_modules <- mutate(icgc_top5_driver_modules, Genes = sapply(Genes, toString))

#############

## Breast cancer drivers
BreastCA_driver_genes <- readr::read_delim("Results/ICGC/IntOGen-DriverGenes_BRCA.tsv", delim = "\t") %>% select(Symbol) %>% unlist() %>% unname()

driverdb_brca <- read_csv("Results/ICGC/driverdb_brca.csv")$gene

BreastCA_driver_genes <- unique(c(BreastCA_driver_genes, driverdb_brca))

View(icgc_driver_net_gene_tbl$`Breast-AdenoCA`[BreastCA_driver_genes,])

##############

## Liver cancer drivers
LiverCA_driver_genes <- readr::read_delim("Results/ICGC/IntOGen-DriverGenes_HC.tsv", delim = "\t") %>% select(Symbol) %>% unlist() %>% unname()

driverdb_hc <- read_csv("Results/ICGC/driverdb_hc.csv")$gene

LiverCA_driver_genes <- unique(c(LiverCA_driver_genes, driverdb_hc))

##############

## Melanoma drivers
Melanoma_driver_genes <- readr::read_delim("Results/ICGC/IntOGen-DriverGenes_CM.tsv", delim = "\t") %>% select(Symbol) %>% unlist() %>% unname()

driverdb_melanoma <- read_csv("Results/ICGC/driverdb_melanoma.csv")$gene

Melanoma_driver_genes <- unique(c(Melanoma_driver_genes, driverdb_melanoma))

##############

## stomach cancer drivers
StomachCA_driver_genes <- readr::read_delim("Results/ICGC/IntOGen-DriverGenes_ST.tsv", delim = "\t") %>% select(Symbol) %>% unlist() %>% unname()

driverdb_stomach_ad <- read_csv("Results/ICGC/driverdb_stomach_ad.csv")$gene

StomachCA_driver_genes <- unique(c(StomachCA_driver_genes, driverdb_stomach_ad))

##############

## AML drivers
AML_driver_genes <- readr::read_delim("Results/ICGC/IntOGen-DriverGenes_AML.tsv", delim = "\t") %>% select(Symbol) %>% unlist() %>% unname()

##############

adult_drivers_list <- list(BreastCA_driver_genes, LiverCA_driver_genes, AML_driver_genes,
                           Melanoma_driver_genes, StomachCA_driver_genes)

names(adult_drivers_list) <- names(icgc_driver_net_modules_flt)

##############

#### Check the rank of modules with driver genes in BrCa

BreastCA_driver_including_module_names <- vector(mode = "character")
for (i in 1:length(icgc_driver_net_modules_flt$`Breast-AdenoCA`)) {
  module_tmp <-
    grep(paste0(paste(paste("^", BreastCA_driver_genes, sep = ""), "$", sep = ""), collapse = "|"), icgc_driver_net_modules_flt$`Breast-AdenoCA`[[i]]$gene)
  if(length(module_tmp >= 1)) { 
    print(names(icgc_driver_net_modules_flt$`Breast-AdenoCA`[i]))
    BreastCA_driver_including_module_names <- append(BreastCA_driver_including_module_names,
                                                     names(icgc_driver_net_modules_flt$`Breast-AdenoCA`[i]))
  }
}

BreastCA_driver_including_module_ranks <- which(names(icgc_driver_net_modules_flt_scores$`Breast-AdenoCA`) %in% BreastCA_driver_including_module_names)
BreastCA_driver_including_module_names <- names(icgc_driver_net_modules_flt_scores$`Breast-AdenoCA`)[BreastCA_driver_including_module_ranks]
BreastCA_driver_genes[which(BreastCA_driver_genes %in% icgc_driver_net_modules_flt$`Breast-AdenoCA`$`597`$gene)]

##################

## Create module 1 for visualization

library(influential)
library(igraph)
library(visNetwork)

icgc_brca_module1 <- subgraph(graph = icgc_driver_net$`Breast-AdenoCA`, 
                              vids = which(as_ids(V(icgc_driver_net$`Breast-AdenoCA`)) %in%  icgc_driver_net_modules_flt$`Breast-AdenoCA`$`1456`$gene))

icgc_brca_module1_node_colors <- rep("lightblue", length(V(icgc_brca_module1)))
icgc_brca_module1_node_colors[which(as_ids(V(icgc_brca_module1)) %in% BreastCA_driver_genes)] <- "orange"

icgc_brca_module1 <- set_vertex_attr(graph = icgc_brca_module1, name = "color", value = icgc_brca_module1_node_colors)

icgc_brca_module1_vis <-
  visIgraph(igraph = icgc_brca_module1) %>% 
  visOptions(
    # nodesIdSelection = TRUE, 
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '25px arial black',
    size = 25,
    # color = c("lightgrey", "orange"),
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "circle" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = c("lightgrey"))

visExport(graph = icgc_brca_module1_vis, 
          name = "icgc_brca_module1_vis", 
          label = "")

## It's also possible to create a subset of two modules and color the nodes based on mark.groups arg
### https://stackoverflow.com/questions/26913419/plot-communities-with-igraph

##########################

# Check the module TP53 is involved in
for (i in 1:length(icgc_driver_net_modules_flt$`Breast-AdenoCA`)) {
  module_tmp <-
    grep("^TP53$", icgc_driver_net_modules_flt$`Breast-AdenoCA`[[i]]$gene)
  if(length(module_tmp >= 1)) { 
    print(names(icgc_driver_net_modules_flt$`Breast-AdenoCA`[i]))
    print(icgc_driver_net_modules_flt$`Breast-AdenoCA`[[i]])
    cat("\n###########################################\n")
  }
}

##########################

## Visualize TP53 and its first-order connected modules

icgc_brca_driverNet_tp53Neighbors <- igraph::neighbors(graph = icgc_driver_net$`Breast-AdenoCA`, v = "TP53", mode = "all") %>% as_ids()
icgc_brca_driverNet_tp53_moduleNeighbors <- vector(mode = "character")

for(i in icgc_driver_net_modules_flt$`Breast-AdenoCA`[names(icgc_driver_net_modules_flt_scores$`Breast-AdenoCA`)[1:5]]) {
  first_orders_tmp = i$gene[which(i$gene %in% icgc_brca_driverNet_tp53Neighbors)]
  if(length(first_orders_tmp) > 0) {
    print(unique(i$module))
    icgc_brca_driverNet_tp53_moduleNeighbors <- append(icgc_brca_driverNet_tp53_moduleNeighbors, first_orders_tmp)
  }
}

library(influential)
library(igraph)
library(visNetwork)

icgc_brca_driverNet_tp53Neighbors_net <- subgraph(graph = icgc_driver_net$`Breast-AdenoCA`, 
                                                  vids = which(as_ids(V(icgc_driver_net$`Breast-AdenoCA`)) %in%  c("TP53", icgc_brca_driverNet_tp53_moduleNeighbors)))

icgc_brca_driverNet_tp53Neighbors_net_vis <-
  visIgraph(igraph = icgc_brca_driverNet_tp53Neighbors_net) %>% 
  visOptions(
    nodesIdSelection = TRUE,
    highlightNearest = TRUE) %>% 
  visNodes(
    font = '25px arial black',
    size = 25,
    # color = c("lightgrey", "orange"),
    shadow = TRUE, 
    labelHighlightBold = TRUE, 
    shape = "circle" # Shapes: "square", "triangle", "box", "circle", "dot", "star", "ellipse", "database", "text", "diamond"
  ) %>% 
  visEdges(color = c("lightgrey"))

visExport(graph = icgc_brca_driverNet_tp53Neighbors_net_vis, 
          name = "icgc_brca_driverNet_tp53Neighbors_net_vis", 
          label = "")

for(i in icgc_driver_net_modules_flt$`Breast-AdenoCA`[names(icgc_driver_net_modules_flt_scores$`Breast-AdenoCA`)[1:5]]) {
  tmp_d = which(i$gene %in% icgc_brca_driverNet_tp53_moduleNeighbors)
  if(length(tmp_d) > 0) {
    print(unique(i$module))
    print("*******")
    print(i$gene[tmp_d])
  }
}

##################################################################################################################################
##################################################################################################################################
##################################################################################################################################

# Evaluate the ICGC results

## Ground truth driver genes
adult_drivers_list

## Top 20 gene of ICGC cancers
icgc_driver_net_gene_tbl_top20 <-
  lapply(icgc_driver_net_gene_tbl, function(i) {
    slice_max(i, n = 20, order_by = primitive_driver_score) %>% select(gene)
  })

icgc_driver_net_gene_tbl_top20 <- do.call(cbind, icgc_driver_net_gene_tbl_top20)
colnames(icgc_driver_net_gene_tbl_top20) <- names(icgc_driver_net_gene_tbl)


########################################################

## Create a dataframe of 1st-ranked driver modules of all cancer types

icgc_net_1st_modules <- data.frame(Number = vector("integer"))

for(i in 1:length(icgc_driver_net_modules_flt)) {
  tmp_module_genes = icgc_driver_net_modules_flt[[i]][names(icgc_driver_net_modules_flt_scores[[i]])[1]][[1]]$gene
  icgc_net_1st_modules <- merge(icgc_net_1st_modules, 
                                data.frame(Number = c(1:length(tmp_module_genes)),
                                           "cancer_type" = tmp_module_genes), 
                                by = "Number", all = TRUE)
}

icgc_net_1st_modules$Number <- NULL

colnames(icgc_net_1st_modules) <- names(icgc_driver_net_modules_flt)

##########################

top_adulthood_evals <- tibble("Cancer type" = vector(mode = "character", length = 5),
                              "IntOGen/DriverDB Source Genes" = vector(mode = "character", length = 5),
                              "Number of Genes in IntOGen/DriverDB (out of top 20)" = vector(mode = "integer", length = 5),
                              "Genes in IntOGen/DriverDB (out of top 20)" = vector(mode = "character", length = 5),
                              "Number of Genes in IntOGen/DriverDB (out of module 1 genes)" = vector(mode = "integer", length = 5),
                              "Genes in IntOGen/DriverDB (out of module 1 genes)" = vector(mode = "character", length = 5),
                              "First Ranked Module Size" = vector(mode = "integer", length = 5))

for (i in 1:length(adult_drivers_list)) {
  top_adulthood_evals[i,1] <- names(adult_drivers_list)[i]
  top_adulthood_evals[i,2] <- paste0(adult_drivers_list[[i]], collapse = ", ")
  top_adulthood_evals[i,3] <- length(which(icgc_driver_net_gene_tbl_top20[[i]] %in% adult_drivers_list[[i]]))
  top_adulthood_evals[i,4] <- paste0(icgc_driver_net_gene_tbl_top20[[i]][which(icgc_driver_net_gene_tbl_top20[[i]] %in% adult_drivers_list[[i]])], collapse = ", ")
  top_adulthood_evals[i,5] <- length(which(na.omit(icgc_net_1st_modules[,i]) %in% adult_drivers_list[[i]]))
  top_adulthood_evals[i,6] <- paste0(na.omit(icgc_net_1st_modules[,i])[which(na.omit(icgc_net_1st_modules[,i]) %in% adult_drivers_list[[i]])], collapse = ", ")
  top_adulthood_evals[i,7] <- length(na.omit(icgc_net_1st_modules[,i]))
}

openxlsx::write.xlsx(top_adulthood_evals, file = "Results/Evaluations/ICGC/top_adulthood_evals.xlsx")

#########################

## visualization of the results

top_adulthood_evals_top20 <- top_adulthood_evals[,c(1,2,7)]
top_adulthood_evals_top20$`Gene group` <- "Top 20 genes"
top_adulthood_evals_top20$`Number of Genes in IntOGen/DriverDB` <- top_adulthood_evals$`Number of Genes in IntOGen/DriverDB (out of top 20)`
top_adulthood_evals_top20$`Genes in IntOGen/DriverDB` <- top_adulthood_evals$`Genes in IntOGen/DriverDB (out of top 20)`
top_adulthood_evals_top20$`First Ranked Module Size` <- ""

top_adulthood_evals_1st_module <- top_adulthood_evals[,c(1,2,7)]
top_adulthood_evals_1st_module$`Gene group` <- "1st ranked module"
top_adulthood_evals_1st_module$`Number of Genes in IntOGen/DriverDB` <- top_adulthood_evals$`Number of Genes in IntOGen/DriverDB (out of module 1 genes)`
top_adulthood_evals_1st_module$`Genes in IntOGen/DriverDB` <- top_adulthood_evals$`Genes in IntOGen/DriverDB (out of module 1 genes)`

top_adulthood_evals_forVis <- rbind(top_adulthood_evals_top20, top_adulthood_evals_1st_module)

###########

icgc_top_eval_plot <- ggplot(data = top_adulthood_evals_forVis, aes(x = `Cancer type`, 
                                                                    y = `Number of Genes in IntOGen/DriverDB`, 
                                                                    fill = `Gene group`)) +
  geom_bar(position = "dodge2", stat = 'identity') +
  geom_text(aes(x = `Cancer type`, 
                y = ifelse(`Number of Genes in IntOGen/DriverDB` > 0, `Number of Genes in IntOGen/DriverDB` - `Number of Genes in IntOGen/DriverDB`/2, 0.5), 
                label = `Genes in IntOGen/DriverDB`), position=position_dodge(width=0.9), size = 1.6, show.legend = FALSE) +
  geom_text(aes(label = paste0(c(rep("", 5), rep("Module size: ", 5)), `First Ranked Module Size`)), 
            hjust = -1, position=position_dodge(width=0.9), size = 1.6) +
  labs(y= '# Genes common with IntOGen/DriverDB drivers') +
  coord_flip() +
  scale_y_continuous(breaks= scales::pretty_breaks()) +
  scale_fill_fish_d(option = "Elagatis_bipinnulata") +
  # scale_color_fish_d(option = "Elagatis_bipinnulata") +
  # scale_fill_viridis_d(begin = 0.32, end = 0.83) +
  theme_bw() +
  theme(text = element_text(size = 7),
        legend.position = "right", 
        legend.box.background = element_rect(colour = "black"))

icgc_top_eval_plot

##############################################################################

## Literature mining of cancer driver results
library(pubmedR)

icgc_lit_driver_data <- tibble(`Cancer type` = vector("character"), 
                               `Gene group` = vector("character"),
                               Gene = vector("character"),
                               Num_PMIDs = vector("integer")
)

for(i in 1:length(icgc_driver_net_gene_tbl_top20)) {
  
  if(names(icgc_driver_net_gene_tbl_top20)[i] == "Breast-AdenoCA") {
    main_search_text <- "((breast cancer[Title/Abstract]) OR (breast carcinoma[Title/Abstract]) OR (breast adenocarcinoma[Title/Abstract]) OR (breast adeno-carcinoma[Title/Abstract]))"
    
  } else if(names(icgc_driver_net_gene_tbl_top20)[i] == "Liver-HCC") {
    main_search_text <- "((liver cancer[Title/Abstract]) OR (hepatic cancer[Title/Abstract]) OR (liver carcinoma[Title/Abstract]) OR (hepatocellular cancer[Title/Abstract]) OR (hepatocellular carcinoma[Title/Abstract]) OR (hepatic carcinoma[Title/Abstract]))"
    
  } else if(names(icgc_driver_net_gene_tbl_top20)[i] == "Myeloid-AML") {
    main_search_text <- "((AML[Title/Abstract]) OR (Acute Myeloid Leukaemia[Title/Abstract]) OR (Acute Myelogenous Leukemia[Title/Abstract]) OR (Acute Myelogenous Leukaemia[Title/Abstract]) OR (Acute Myeloid Leukemia[Title/Abstract]))"
    
  } else if(names(icgc_driver_net_gene_tbl_top20)[i] == "Skin-Melanoma") {
    main_search_text <- "((melanoma[Title/Abstract]) OR (skin carcinoma[Title/Abstract]) OR (skin cutaneous carcinoma[Title/Abstract]) OR (skin cancer[Title/Abstract]))"
    
  } else if(names(icgc_driver_net_gene_tbl_top20)[i] == "Stomach-AdenoCA") {
    main_search_text <- "((Stomach cancer[Title/Abstract]) OR (Stomach carcinoma[Title/Abstract]) OR (Stomach adenocarcinoma[Title/Abstract]) OR (Stomach adeno-carcinoma[Title/Abstract]) OR (Gastric cancer[Title/Abstract]) OR (Gastric carcinoma[Title/Abstract]) OR (Gastric adenocarcinoma[Title/Abstract]) OR (Gastric adeno-carcinoma[Title/Abstract]))"
    
  }
  
  ### Add top 20 genes
  tmp_top20_tbl <- tibble(`Cancer type` = names(icgc_driver_net_gene_tbl_top20)[i],
                          `Gene group` = "Top 20 genes",
                          Gene = (icgc_driver_net_gene_tbl_top20[,i]))
  
  tmp_top20_pmid_list <- lapply(1:20, function(j) {
    query_text <- paste0("(", tmp_top20_tbl$Gene[j], "[Title/Abstract]) AND ", main_search_text)
    tmp_entrez <- pmQueryTotalCount(query = query_text, api_key = "1a1d2b704b3bad87d1ded4f2d95d611e1909")
    tmp_entrez_num <- tmp_entrez$total_count
    print(paste0("Top gene ", j, " is done!"))
    return(tmp_entrez_num)
  })
  
  tmp_top20_tbl$Num_PMIDs <- unlist(tmp_top20_pmid_list)
  
  icgc_lit_driver_data <- rbind(icgc_lit_driver_data, tmp_top20_tbl)
  
  ### Add the 1st ranked module
  tmp_1st_module_tbl <- tibble(`Cancer type` = names(icgc_driver_net_gene_tbl_top20)[i],
                               `Gene group` = "1st-ranked module",
                               Gene = (na.omit(icgc_net_1st_modules[,i])))
  
  tmp_1st_module_pmid_list <- lapply(1:length(na.omit(icgc_net_1st_modules[,i])), function(m) {
    query_text <- paste0("(", tmp_1st_module_tbl$Gene[m], "[Title/Abstract]) AND ", main_search_text)
    tmp_entrez <- pmQueryTotalCount(query = query_text, api_key = "1a1d2b704b3bad87d1ded4f2d95d611e1909")
    tmp_entrez_num <- tmp_entrez$total_count
    print(paste0("Module 1 gene ", m, " is done!"))
    return(tmp_entrez_num)
  })
  
  tmp_1st_module_tbl$Num_PMIDs <- unlist(tmp_1st_module_pmid_list)
  
  icgc_lit_driver_data <- rbind(icgc_lit_driver_data, tmp_1st_module_tbl)
  
  print(paste0(names(icgc_driver_net_gene_tbl_top20)[i], " is done!"))
  
}

### Save the tables
write.table(x = icgc_lit_driver_data, file = "Results/Evaluations/ICGC/icgc_lit_driver_data.csv", append = F, quote = F, sep = ",", row.names = F)

######################

## Visualization of the results

icgc_driver_lit_mining_plot <- ggplot(data = icgc_lit_driver_data, aes(x = `Cancer type`, 
                                                                       y = log2(Num_PMIDs+1), 
                                                                       fill = `Gene group`)) +
  geom_boxplot(position = position_dodge(0.8), outlier.size = 0.01, width=0.58, size = 0.2) +
  labs(y= 'Log2(number of papers + 1)',
       title = "Literature-based Association of Top\nDrivers with Their Respective Cancers") +
  coord_flip() +
  scale_y_continuous(breaks= scales::pretty_breaks()) +
  scale_fill_fish_d(option = "Elagatis_bipinnulata") +
  # scale_color_fish_d(option = "Elagatis_bipinnulata") +
  # scale_fill_viridis_d(begin = 0.32, end = 0.83) +
  theme_bw() +
  theme(text = element_text(size = 7),
        legend.position = "right", 
        legend.box.background = element_rect(colour = "black"))

icgc_driver_lit_mining_plot

########

### Combine plots

icgc_all_eval_plot <- (icgc_top_eval_plot + theme(legend.position = "none")) + 
  (icgc_driver_lit_mining_plot + xlab(NULL)) + patchwork::plot_layout(widths = c(2,1)) +
  plot_annotation(tag_levels = list(c("A", "B"))) &
  theme(plot.tag = element_text(face = 'bold'))

icgc_all_eval_plot

ggsave(filename = "Results/Evaluations/ICGC/icgc_all_eval_plot.pdf", 
       plot = icgc_all_eval_plot, 
       device = "pdf", width = 22, height = 7, units = "cm")

##############################################################################
##############################################################################

## ORA on the first ranked risk modules

library(enrichR)

### Using enrichR

#Get the list of all databases
enrichR_dbs <- enrichR::listEnrichrDbs()
enrichR_dbs %>% as.data.frame() %>% view()

#Select you desired databases
dbs.GO.BP <- "GO_Biological_Process_2023"
dbs.KEGG.pathway <- "KEGG_2021_Human"

#ORA on GO-BP terms for 1st risk module of each cancer type

icgc_driver_1st_module_ORA <- lapply(1:5, function(i) {
  enrichR::enrichr(genes = as.character(na.omit(icgc_net_1st_modules[,i])),
                   databases = dbs.GO.BP)[[1]]
})

names(icgc_driver_1st_module_ORA) <- colnames(icgc_net_1st_modules)

###########################

### save the results
library(openxlsx)

write.xlsx(x = icgc_driver_1st_module_ORA, file = "Results/Evaluations/ICGC/icgc_driver_1st_module_ORA.xlsx")

###########################

### visualize the results

 # Create a table of the first ranked GO-BP of each cancer based on the combined scores
icgc_driver_1st_module_ORA_1st_ranked <- tibble("Term" = vector(mode = "character", length = 5),
                                                "Padj" = vector(mode = "numeric", length = 5),
                                                "Combined score" = vector(mode = "numeric", length = 5),
                                                "Genes" = vector(mode = "character", length = 5),
                                                "Cancer type" = vector(mode = "character", length = 5))

for (i in 1:length(icgc_driver_1st_module_ORA)) {
  icgc_driver_1st_module_ORA_1st_ranked[i,1] <- icgc_driver_1st_module_ORA[[i]] %>% slice_max(Combined.Score, with_ties = F) %>% select(Term) %>% unlist() %>% unname()
  icgc_driver_1st_module_ORA_1st_ranked[i,2] <- icgc_driver_1st_module_ORA[[i]] %>% slice_max(Combined.Score, with_ties = F) %>% select(Adjusted.P.value) %>% unlist() %>% unname()
  icgc_driver_1st_module_ORA_1st_ranked[i,3] <- icgc_driver_1st_module_ORA[[i]] %>% slice_max(Combined.Score, with_ties = F) %>% select(Combined.Score) %>% unlist() %>% unname()
  icgc_driver_1st_module_ORA_1st_ranked[i,4] <- icgc_driver_1st_module_ORA[[i]] %>% slice_max(Combined.Score, with_ties = F) %>% select(Genes) %>% unlist() %>% unname()
  icgc_driver_1st_module_ORA_1st_ranked[i,5] <- names(icgc_driver_1st_module_ORA)[i]
}

icgc_driver_1st_module_ORA_vis <- ggplot(icgc_driver_1st_module_ORA_1st_ranked, 
                                         aes(x = Term, y = -log10(Padj),
                                             color = `Cancer type`,
                                             size = log10(`Combined score`))) +  
    geom_point() +
  scale_size_identity(guide = "legend") +
    ylab("-log10(Padj)") + xlab(NULL) + 
    scale_color_fish_d(option = "Hexagrammos_lagocephalus", name="Cancer type") +
    coord_flip() +
    labs(title = expression("ORA of the"~"1"^"st"~ "Ranked ICGC Modules")) +
    scale_y_continuous(breaks= scales::pretty_breaks()) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 35)) +
    theme_classic() + theme(
      legend.position = "right",
      # legend.position = c(.95, 0.02),
      # legend.justification = c("right", "bottom"),
      legend.box.background = element_rect(color="black", size=1),
      text = element_text(size = 7))

icgc_driver_1st_module_ORA_vis
  
  #Saving the plot
  ggsave(filename = paste0("Results/Evaluations/ICGC/icgc_driver_1st_module_ORA_vis.pdf"), 
         plot = icgc_driver_1st_module_ORA_vis,
         device = "pdf",width = 14, height = 5, units = "cm")

  ##########################################################################################
  
  ## Evaluation of the first ranked mediators
  
  ### Evaluate the first ranked mediators and their first-order connected modules in driver networks
  
  icgc_driver_mediator_tbl <- tibble("Cancer type" = vector(mode = "character", length = 5), 
                                     "1st Ranked Driverness mediator" = vector(mode = "character", length = 5), 
                                     "1st order connections" = vector(mode = "character", length = 5), 
                                     "Commonalities with 1st module" = vector(mode = "character", length = 5), 
                                     "Commonalities with 2nd module" = vector(mode = "character", length = 5), 
                                     "Commonalities with 3rd module" = vector(mode = "character", length = 5), 
                                     "Commonalities with 4th module" = vector(mode = "character", length = 5), 
                                     "Commonalities with 5th module" = vector(mode = "character", length = 5))
  
  icgc_driver_mediator_list <- list()
  
  for(j in 1:5) {
    icgc_driver_mediator_tbl[j,1] <- names(icgc_driver_net_gene_tbl)[j]
    
    tmp_first_ranked_mediator <- icgc_driver_net_gene_tbl[[j]] %>%
      slice_max(mediator_score, n = 1) %>%
      select(gene) %>%
      unlist() %>%
      unname()
    icgc_driver_mediator_tbl[j,2] <- tmp_first_ranked_mediator
    
    tmp_first_ranked_mediator_Neighbors <- igraph::neighbors(graph = icgc_driver_net[[j]], v = tmp_first_ranked_mediator, mode = "all") %>% as_ids()
    icgc_driver_mediator_tbl[j,3] <- paste0(tmp_first_ranked_mediator_Neighbors, collapse = ", ")
    
    tmp_first_ranked_mediator_moduleNeighbors <- vector(mode = "character")
    
    for(i in 1:5) {
      tmp_module <- icgc_driver_net_modules_flt[[j]][[names(icgc_driver_net_modules_flt_scores[[j]])[i]]]
      first_orders_tmp = tmp_module$gene[which(tmp_module$gene %in% tmp_first_ranked_mediator_Neighbors)]
      icgc_driver_mediator_tbl[j,i+3] <- paste0(first_orders_tmp, collapse = ", ")
      
      if(length(first_orders_tmp) > 0) {
        tmp_first_ranked_mediator_moduleNeighbors <- append(tmp_first_ranked_mediator_moduleNeighbors, first_orders_tmp)
      }
    }
    icgc_driver_mediator_list <- append(icgc_driver_mediator_list, list(tmp_first_ranked_mediator_moduleNeighbors))
  }
  
  readr::write_delim(x = icgc_driver_mediator_tbl,
                     file = "Results/Evaluations/ICGC/icgc_driver_mediator_tbl.csv", 
                     delim = ",")
  
  names(icgc_driver_mediator_list) <- names(icgc_driver_net_gene_tbl)
  
  ###########
  
  #### Visualization of results of Breast-AdenoCA
  brca_icgc_driver_mediator_neighbors_net <- subgraph(graph = icgc_driver_net$`Breast-AdenoCA`, 
                                                      vids = which(as_ids(V(icgc_driver_net$`Breast-AdenoCA`)) %in%  c(icgc_driver_mediator_tbl[icgc_driver_mediator_tbl[,1] == "Breast-AdenoCA",2], 
                                                                                                                       icgc_driver_mediator_list$`Breast-AdenoCA`)))
  
  brca_icgc_driver_mediator_neighbors_net_tbl <- igraph::as_data_frame(brca_icgc_driver_mediator_neighbors_net)
  readr::write_delim(x = brca_icgc_driver_mediator_neighbors_net_tbl,
                     file = "Results/Evaluations/ICGC/brca_driver_mediator_neighbors_net_tbl.txt", 
                     delim = "\t")
  
##############################################################################
##############################################################################
  
  ## Save the gene data and modules as supplementary files
  
  tmp_icgc_net_gene_tbl <- icgc_driver_net_gene_tbl
  
  for(i in 1:5) {
    colnames(tmp_icgc_net_gene_tbl[[i]]) <- c("Gene", "IVI", "Mean_deletriousness_score", "Driver_score", "Mean_neighborhood_weight", "Driverness_mediator_score")
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_icgc_net_gene_tbl, file = "Results/Evaluations/ICGC/icgc_driver_net_gene_tbl.xlsx")
  
  rm(tmp_icgc_net_gene_tbl)
  
  ############################
  
  ### save the modules
  
  tmp_icgc_net_modules_flt <- icgc_driver_net_modules_flt
  
  tmp_icgc_net_modules_flt_scores <- icgc_driver_net_modules_flt_scores
  
  tmp_icgc_net_modules_flt <- lapply(tmp_icgc_net_modules_flt, function(i) {
    
    tmp_data <- i
    tmp_data <- tmp_data[match(names(tmp_icgc_net_modules_flt_scores[[which(sapply(tmp_icgc_net_modules_flt_scores, length) == length(tmp_data))]]),
                               names(i))]
    
    tmp_data
    
  })
  
  for(i in 1:5) {
    tmp_icgc_net_modules_flt[[i]] <-
      lapply(1:length(tmp_icgc_net_modules_flt[[i]]), function(j) {
        tmp_data <- tmp_icgc_net_modules_flt[[i]][[j]]
        colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
        tmp_data[tmp_data == 1] <- "Yes"
        tmp_data[tmp_data == 0] <- "No"
        tmp_data$Module <- paste0("Module number ", j)
        tmp_data <- data.frame(Module = unique(tmp_data$Module),
                               Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                               Mean_score = unique(tmp_data$Mean_score))
        tmp_data
      })
    tmp_icgc_net_modules_flt[[i]] <- do.call(rbind, tmp_icgc_net_modules_flt[[i]])
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_icgc_net_modules_flt, file = "Results/Evaluations/ICGC/icgc_driver_net_modules_flt.xlsx")
  
  rm(tmp_icgc_net_modules_flt, tmp_icgc_net_modules_flt_scores)

  
##############################################################################
  
#=============================================================================
#
#    Code chunk 14: Test the pipeline based on permutation ICGC data
#
#=============================================================================
  
  # Test the pipeline based on icgc_permutation (shuffling only PPI data)
  
  ## Prepare cor data
  
  ### Setting snv cor to the ones derived from original data
  onlyPPI_icgc_permut_snv_flt_cor <- icgc_snv_data_flt_cor
  
  ### Setting expr cor to the ones derived from NON-icgc_permutated data
  onlyPPI_icgc_permut_rnaseq_expr_cor <- icgc_tx_cor_list

  # Merge the similarity tables
  onlyPPI_icgc_permut_driver_marker_net <- onlyPPI_icgc_permut_snv_flt_cor
  
  # Add the co-expression data to the tables
  
  for(i in 1:length(onlyPPI_icgc_permut_driver_marker_net)) {
    
    genes <- onlyPPI_icgc_permut_driver_marker_net[[i]][,c(1,2)] %>% 
      unlist() %>% 
      unname() %>% 
      unique()
    
    first_level_coex_index <- c(which(onlyPPI_icgc_permut_rnaseq_expr_cor[[i]]$from %in% genes), 
                                which(onlyPPI_icgc_permut_rnaseq_expr_cor[[i]]$to %in% genes)) %>% unique()
    
    first_level_coex_genes <- c(genes, 
                                onlyPPI_icgc_permut_rnaseq_expr_cor[[i]]$from[first_level_coex_index],
                                onlyPPI_icgc_permut_rnaseq_expr_cor[[i]]$to[first_level_coex_index]) %>% unique()
    
    first_two_levels_coex_index <- c(which(onlyPPI_icgc_permut_rnaseq_expr_cor[[i]]$from %in% first_level_coex_genes), 
                                     which(onlyPPI_icgc_permut_rnaseq_expr_cor[[i]]$to %in% first_level_coex_genes)) %>% unique()
    
    onlyPPI_icgc_permut_driver_marker_net[[i]] <- rbind(onlyPPI_icgc_permut_driver_marker_net[[i]], onlyPPI_icgc_permut_rnaseq_expr_cor[[i]][first_two_levels_coex_index,])
    
  }
  
  # Add PPI data to the driver/marker net
  
  # Map gene names to stringDB ids for PPI analysis (STRING  11.5)
  
  ## get the gene names in each risk net
  onlyPPI_icgc_permut_driver_marker_net_genes <- lapply(onlyPPI_icgc_permut_driver_marker_net, function(i) {
    data.frame(gene_symbol = i[,c(1,2)] %>% 
                 unlist() %>% 
                 unname() %>% 
                 unique())
  })
  
  ## getSTRINGdb for human
  string_db <- STRINGdb$new(species=9606, version = "11.5",
                            score_threshold=900, # very high confidence (>0.9)
                            input_directory=""
  )
  
  onlyPPI_icgc_permut_driver_marker_net_genes <- lapply(onlyPPI_icgc_permut_driver_marker_net_genes, function(i) {
    
    string_db$map(my_data_frame = i, 
                  my_data_frame_id_col_names = "gene_symbol", 
                  takeFirst = TRUE, removeUnmappedRows = TRUE)
    
  })
  
  # Get the PPIs
  onlyPPI_icgc_permut_driver_marker_net_ppi <- lapply(onlyPPI_icgc_permut_driver_marker_net_genes, function(i) {
    
    string_db$get_interactions(string_ids = i$STRING_id)
    
  })
  
  # Convert stringdb ids to gene symbols
  onlyPPI_icgc_permut_driver_marker_net_ppi <- lapply(1:length(onlyPPI_icgc_permut_driver_marker_net_ppi), function(i) {
    tmp_row_index <- match(onlyPPI_icgc_permut_driver_marker_net_ppi[[i]]$from, onlyPPI_icgc_permut_driver_marker_net_genes[[i]]$STRING_id)
    onlyPPI_icgc_permut_driver_marker_net_ppi[[i]]$from <- onlyPPI_icgc_permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_row_index]
    
    tmp_column_index <- match(onlyPPI_icgc_permut_driver_marker_net_ppi[[i]]$to, onlyPPI_icgc_permut_driver_marker_net_genes[[i]]$STRING_id)
    onlyPPI_icgc_permut_driver_marker_net_ppi[[i]]$to <- onlyPPI_icgc_permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_column_index]
    
    tmp_tbl <- onlyPPI_icgc_permut_driver_marker_net_ppi[[i]][,c(1,2)]
    tmp_tbl$type <- "PPI"
    
    ## Remove duplicates
    tmp_tbl <- tmp_tbl[-which(duplicated(paste(tmp_tbl$from,
                                               tmp_tbl$to, 
                                               sep = "_"))),]
    
    ## Shuffle the data
    set.seed(1234)
    tmp_tbl$from <- tmp_tbl$from[sample(nrow(tmp_tbl))]
    
    return(tmp_tbl)
  })
  
  names(onlyPPI_icgc_permut_driver_marker_net_ppi) <- names(onlyPPI_icgc_permut_driver_marker_net_genes)
  
  # Combine PPI with the original driver/marker net
  onlyPPI_icgc_permut_driver_marker_net <- lapply(1:length(onlyPPI_icgc_permut_driver_marker_net), function(i) {
    rbind(onlyPPI_icgc_permut_driver_marker_net[[i]], onlyPPI_icgc_permut_driver_marker_net_ppi[[i]])
  })
  
  names(onlyPPI_icgc_permut_driver_marker_net) <- names(onlyPPI_icgc_permut_driver_marker_net_genes)
  
  # Reconstruct the networks
  onlyPPI_icgc_permut_driver_marker_net <- lapply(onlyPPI_icgc_permut_driver_marker_net, function(i) {
    igraph::graph_from_data_frame(i, directed=FALSE)
  })
  
  ## Check the network layers
  table(igraph::E(onlyPPI_icgc_permut_driver_marker_net$`Breast-AdenoCA`)$type)
  
  ###############################
  
  # Calculate the IVI values of genes within each net
  library(influential)
  onlyPPI_icgc_permut_driver_marker_net_ivi <- lapply(onlyPPI_icgc_permut_driver_marker_net, function(i) ivi(i, verbose = T))
  
  ###############################
  
  # Calculate the primitive gene driver scores
  onlyPPI_icgc_permut_driver_marker_net_gene_tbl <- lapply(onlyPPI_icgc_permut_driver_marker_net_ivi, function(i) {
    data.frame(gene = names(i), ivi = i)
  })
  
  ## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer
  
  onlyPPI_icgc_permut_snv_mean_deletriousness <- lapply(icgc_snv_data_flt_withCADD_summarized, function(i) {
    colMeans(i)
  })
  
  ## Generate the gene table for each net and calculate the primitive scores
  onlyPPI_icgc_permut_driver_marker_net_gene_tbl <- lapply(1:length(onlyPPI_icgc_permut_driver_marker_net_gene_tbl), function(i) {
    tmp_tbl <- onlyPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]
    tmp_tbl$mean_deletriousness <- 0
    match_index <- which(names(onlyPPI_icgc_permut_snv_mean_deletriousness[[i]]) %in% rownames(tmp_tbl))
    tmp_tbl[names(onlyPPI_icgc_permut_snv_mean_deletriousness[[i]])[match_index], "mean_deletriousness"] <- 
      onlyPPI_icgc_permut_snv_mean_deletriousness[[i]][match_index]
    
    tmp_tbl$primitive_driver_score <- (tmp_tbl$ivi)*(tmp_tbl$mean_deletriousness) # node weights
    
    return(tmp_tbl)
    
  })
  
  names(onlyPPI_icgc_permut_driver_marker_net_gene_tbl) <- names(onlyPPI_icgc_permut_snv_mean_deletriousness)
  
  # add gene scores
  onlyPPI_icgc_permut_driver_marker_net <- lapply(1:length(onlyPPI_icgc_permut_driver_marker_net), function(i) {
    set_vertex_attr(graph = onlyPPI_icgc_permut_driver_marker_net[[i]], name = "score", value = onlyPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score)
  })
  
  names(onlyPPI_icgc_permut_driver_marker_net) <- names(onlyPPI_icgc_permut_snv_mean_deletriousness)
  
  # Convert the un-weighted to a node-weighted network
  onlyPPI_icgc_permut_driver_marker_net <- lapply(1:length(onlyPPI_icgc_permut_driver_marker_net), function(i) {
    set_vertex_attr(graph = onlyPPI_icgc_permut_driver_marker_net[[i]], name = "weight", value = (onlyPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]$ivi)^(onlyPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness))
  })
  
  names(onlyPPI_icgc_permut_driver_marker_net) <- names(onlyPPI_icgc_permut_snv_mean_deletriousness)
  
  ## Check the network node weights
  summary(igraph::V(onlyPPI_icgc_permut_driver_marker_net$`Breast-AdenoCA`)$weight)
  
  # Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
  onlyPPI_icgc_permut_driver_marker_net_gene_tbl <- lapply(1:length(onlyPPI_icgc_permut_driver_marker_net), function(i) {
    tmp_tbl <- onlyPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]
    tmp_tbl$mean_neighborhood_weight <- sapply(1:nrow(tmp_tbl), function(j) {
      igraph::neighbors(graph = onlyPPI_icgc_permut_driver_marker_net[[i]], 
                        v = V(onlyPPI_icgc_permut_driver_marker_net[[i]])[j])$score %>% mean()
    })
    
    tmp_tbl$mediator_score <- tmp_tbl$primitive_driver_score * tmp_tbl$mean_neighborhood_weight
    
    return(tmp_tbl)
  })
  
  names(onlyPPI_icgc_permut_driver_marker_net_gene_tbl) <- names(onlyPPI_icgc_permut_snv_mean_deletriousness)
  
  # Detect communities/modules/clusters
  onlyPPI_icgc_permut_driver_marker_net_modules <- lapply(onlyPPI_icgc_permut_driver_marker_net, function(i) {
    set.seed(3847)
    igraph::cluster_leiden(
      graph = i,
      objective_function = "CPM",
      weights = NULL,
      resolution_parameter = 0.5,
      beta = 0.05,
      initial_membership = NULL,
      n_iterations = 100,
      vertex_weights = igraph::V(i)$weight
    )
  })
  
  sizes(onlyPPI_icgc_permut_driver_marker_net_modules$`Breast-AdenoCA`) %>% as.vector() %>% summary()
  length(onlyPPI_icgc_permut_driver_marker_net_modules$`Breast-AdenoCA`)
  
  ## Inspect and filter modules
  onlyPPI_icgc_permut_driver_marker_net_modules_flt <- lapply(1:length(onlyPPI_icgc_permut_driver_marker_net_modules), function(i) {
    
    ### Create a table of modules and their genes
    tmp_tbl <- data.frame(
      module = membership(onlyPPI_icgc_permut_driver_marker_net_modules[[i]]) %>% as.integer(),
      gene = onlyPPI_icgc_permut_driver_marker_net_modules[[i]]$names
    )
    
    ### Add if the gene is mutated or not (could be a risk gene or not)
    tmp_tbl$mutated <- 0
    tmp_tbl$mutated[which(onlyPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness > 0)] <- 1
    
    ### Separate modules
    tmp_tbl <- lapply(unique(tmp_tbl$module), function(j){
      subset(tmp_tbl, module == j)
    })
    
    names(tmp_tbl) <- sapply(tmp_tbl, function(m) {
      unique(m$module)
    })
    
    ### Remove modules that have less than 2 Risk genes or their size is less than 4
    tmp_tbl <- lapply(tmp_tbl, function(k) {
      if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
        k
      } else {
        NULL
      }
    })
    
    ### Remove NULL modules
    tmp_tbl[which(sapply(tmp_tbl, is.null))] <- NULL
    
    ### Calculate the module scores
    tmp_tbl <- lapply(tmp_tbl, function(l) {
      tmp_tbl4score <- cbind(l, mean_score = 0)
      tmp_tbl4score$mean_score <- mean(onlyPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score[which(onlyPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]$gene %in% l$gene)])
      tmp_tbl4score
    })
    
    tmp_tbl
  })
  
  names(onlyPPI_icgc_permut_driver_marker_net_modules_flt) <- names(onlyPPI_icgc_permut_driver_marker_net_modules)
  
  sapply(onlyPPI_icgc_permut_driver_marker_net_modules_flt, length)
  sapply(driver_marker_net_modules_flt, length)
  
  ## Sort the mean scores of all modules of each cancer type
  onlyPPI_icgc_permut_driver_marker_net_modules_flt_scores <-
    lapply(onlyPPI_icgc_permut_driver_marker_net_modules_flt, function(i) {
      tmp_tbl <-
        sapply(i, function(j) {unique(j$mean_score)})
      tmp_tbl <- rev(sort(tmp_tbl))
      tmp_tbl
    })
  
  ### View the number of first ranked modules
  sapply(onlyPPI_icgc_permut_driver_marker_net_modules_flt_scores, function(i) {i[1]}) %>% as.data.frame() %>% View()
  
  onlyPPI_icgc_permut_driver_marker_net_modules_flt$`Breast-AdenoCA`$`2114` %>% View()
  
  ########
  
  # Top 20 gene 
  onlyPPI_icgc_permut_driver_marker_net_gene_tbl_top20 <-
    lapply(onlyPPI_icgc_permut_driver_marker_net_gene_tbl, function(i) {
      slice_max(i, n = 20, order_by = primitive_driver_score) %>% select(gene)
    })
  
  onlyPPI_icgc_permut_driver_marker_net_gene_tbl_top20 <- do.call(cbind, onlyPPI_icgc_permut_driver_marker_net_gene_tbl_top20)
  colnames(onlyPPI_icgc_permut_driver_marker_net_gene_tbl_top20) <- names(onlyPPI_icgc_permut_driver_marker_net_gene_tbl)
  
  ###########################################
  ###########################################
  ###########################################
  
  ## Save the gene data and modules as supplementary files
  
  tmp_onlyPPI_icgc_permut_driver_marker_net_gene_tbl <- onlyPPI_icgc_permut_driver_marker_net_gene_tbl

  for(i in 1:5) {
    colnames(tmp_onlyPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]) <- c("Gene", "IVI", "Mean_deletriousness_score", "Driver_score", "Mean_neighborhood_weight", "Driverness_mediator_score")
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_onlyPPI_icgc_permut_driver_marker_net_gene_tbl, file = "Results/Evaluations/icgc_permutation/onlyPPI_icgc_permut_driver_marker_net_gene_tbl.xlsx")
  
  rm(tmp_onlyPPI_icgc_permut_driver_marker_net_gene_tbl)
  
  ############################
  
  ### save the modules
  
  tmp_onlyPPI_icgc_permut_driver_marker_net_modules_flt <- onlyPPI_icgc_permut_driver_marker_net_modules_flt

  tmp_onlyPPI_icgc_permut_driver_marker_net_modules_flt_scores <- onlyPPI_icgc_permut_driver_marker_net_modules_flt_scores
  
  tmp_onlyPPI_icgc_permut_driver_marker_net_modules_flt <- lapply(tmp_onlyPPI_icgc_permut_driver_marker_net_modules_flt, function(i) {
    
    tmp_data <- i
    tmp_data <- tmp_data[rev(order(sapply(tmp_data, function(j) unique(j$mean_score))))]
    
    tmp_data
    
  })
  
  for(i in 1:5) {
    tmp_onlyPPI_icgc_permut_driver_marker_net_modules_flt[[i]] <-
      lapply(1:length(tmp_onlyPPI_icgc_permut_driver_marker_net_modules_flt[[i]]), function(j) {
        tmp_data <- tmp_onlyPPI_icgc_permut_driver_marker_net_modules_flt[[i]][[j]]
        colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
        tmp_data[tmp_data == 1] <- "Yes"
        tmp_data[tmp_data == 0] <- "No"
        tmp_data$Module <- paste0("Module number ", j)
        tmp_data <- data.frame(Module = unique(tmp_data$Module),
                               Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                               Mean_score = unique(tmp_data$Mean_score))
        tmp_data
      })
    tmp_onlyPPI_icgc_permut_driver_marker_net_modules_flt[[i]] <- do.call(rbind, tmp_onlyPPI_icgc_permut_driver_marker_net_modules_flt[[i]])
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_onlyPPI_icgc_permut_driver_marker_net_modules_flt, file = "Results/Evaluations/icgc_permutation/onlyPPI_icgc_permut_driver_marker_net_modules_flt.xlsx")
  
  rm(tmp_onlyPPI_icgc_permut_driver_marker_net_modules_flt, tmp_onlyPPI_icgc_permut_driver_marker_net_modules_flt_scores)
  
  ###############################################################################################################
  ###############################################################################################################
  ###############################################################################################################
  
  # Test the pipeline based on icgc_permutation (shuffling both deleteriousness and expression data)
  
  ## Prepare somatic snv data
  
  icgc_permut_snv_flt_withCADD_summarized <- icgc_snv_data_flt_withCADD_summarized
  
  ## Shuffle the data
  icgc_permut_snv_flt_withCADD_summarized <- lapply(icgc_permut_snv_flt_withCADD_summarized, function(i) {
    tmp_colnames <- colnames(i)
    
    ## shuffle
    tmp_data <- apply(i, 1, function(j) {
      set.seed(1234)
      j[sample(length(j))]
    })
    tmp_data <- as.data.frame(t(tmp_data))
    colnames(tmp_data) <- tmp_colnames
    
    ## re-shuffle
    tmp_data <- apply(i, 1, function(j) {
      set.seed(5678)
      j[sample(length(j))]
    })
    tmp_data <- as.data.frame(t(tmp_data))
    colnames(tmp_data) <- tmp_colnames
    
    tmp_data
  })
  
  ## Perform correlation analysis
  
  icgc_permut_snv_flt_cor <- lapply(icgc_permut_snv_flt_withCADD_summarized, function(i){
    tmp_cor <- fcor(i)
    tmp_cor <- subset(tmp_cor, mr < 20)
    tmp_cor
  })
  
  ##################
  
  ## Prepare expression data
  icgc_permut_rnaseq_expr_list <- icgc_tx_data
  
  ## Shuffle the data
  icgc_permut_rnaseq_expr_list <- lapply(icgc_permut_rnaseq_expr_list, function(i) {
    tmp_colnames <- colnames(i)
    
    ## shuffle
    tmp_data <- apply(i, 1, function(j) {
      set.seed(1234)
      j[sample(length(j))]
    })
    tmp_data <- as.data.frame(t(tmp_data))
    colnames(tmp_data) <- tmp_colnames
    
    ## re-shuffle
    tmp_data <- apply(i, 1, function(j) {
      set.seed(5678)
      j[sample(length(j))]
    })
    tmp_data <- as.data.frame(t(tmp_data))
    colnames(tmp_data) <- tmp_colnames
    
    tmp_data
  })
  
  ## Perform correlation analysis
  icgc_permut_rnaseq_expr_cor <- lapply(icgc_permut_rnaseq_expr_list, function(i){
    
    tmp_cor <- fcor(data = i, 
                    method = "spearman", 
                    mutualRank = TRUE, pvalue = FALSE, flat = TRUE)
    
    tmp_cor <- subset(tmp_cor, mr < 20)
    tmp_cor
  })
  
  ###################
  
  # Prepare the tables for merging with other data
  icgc_permut_snv_flt_cor <- lapply(icgc_permut_snv_flt_cor, function(i) {
    tmp_tbl <- i[,c(1,2)]
    colnames(tmp_tbl) <- c("from", "to")
    tmp_tbl$type <- "co-deletriousness"
    tmp_tbl
  })
  
  icgc_permut_rnaseq_expr_cor <- lapply(icgc_permut_rnaseq_expr_cor, function(i) {
    tmp_tbl <- i[,c(1,2)]
    colnames(tmp_tbl) <- c("from", "to")
    tmp_tbl$type <- "co-expression"
    tmp_tbl
  })
  
  # Merge the similarity tables
  icgc_permut_driver_marker_net <- icgc_permut_snv_flt_cor
  
  # Make the orders the same
  icgc_permut_rnaseq_expr_cor <- icgc_permut_rnaseq_expr_cor[names(icgc_permut_snv_flt_cor)]
  
  # Add the co-expression data to the tables
  
  for(i in 1:length(icgc_permut_driver_marker_net)) {
    
    genes <- icgc_permut_driver_marker_net[[i]][,c(1,2)] %>% 
      unlist() %>% 
      unname() %>% 
      unique()
    
    first_level_coex_index <- c(which(icgc_permut_rnaseq_expr_cor[[i]]$from %in% genes), 
                                which(icgc_permut_rnaseq_expr_cor[[i]]$to %in% genes)) %>% unique()
    
    first_level_coex_genes <- c(genes, 
                                icgc_permut_rnaseq_expr_cor[[i]]$from[first_level_coex_index],
                                icgc_permut_rnaseq_expr_cor[[i]]$to[first_level_coex_index]) %>% unique()
    
    first_two_levels_coex_index <- c(which(icgc_permut_rnaseq_expr_cor[[i]]$from %in% first_level_coex_genes), 
                                     which(icgc_permut_rnaseq_expr_cor[[i]]$to %in% first_level_coex_genes)) %>% unique()
    
    icgc_permut_driver_marker_net[[i]] <- rbind(icgc_permut_driver_marker_net[[i]], icgc_permut_rnaseq_expr_cor[[i]][first_two_levels_coex_index,])
    
  }
  
  # Add PPI data to the driver/marker net
  
  # Map gene names to stringDB ids for PPI analysis (STRING  11.5)
  
  ## get the gene names in each risk net
  icgc_permut_driver_marker_net_genes <- lapply(icgc_permut_driver_marker_net, function(i) {
    data.frame(gene_symbol = i[,c(1,2)] %>% 
                 unlist() %>% 
                 unname() %>% 
                 unique())
  })
  
  ## getSTRINGdb for human
  ### WARNING: You didn't specify a version of the STRING database to use. Hence we will use STRING  11.5
  string_db <- STRINGdb$new(species=9606, version = "11.5", 
                            score_threshold=900, # very high confidence (>0.9)
                            input_directory=""
  )
  
  icgc_permut_driver_marker_net_genes <- lapply(icgc_permut_driver_marker_net_genes, function(i) {
    
    string_db$map(my_data_frame = i, 
                  my_data_frame_id_col_names = "gene_symbol", 
                  takeFirst = TRUE, removeUnmappedRows = TRUE)
    
  })
  
  # Get the PPIs
  icgc_permut_driver_marker_net_ppi <- lapply(icgc_permut_driver_marker_net_genes, function(i) {
    
    string_db$get_interactions(string_ids = i$STRING_id)
    
  })
  
  # Convert stringdb ids to gene symbols
  icgc_permut_driver_marker_net_ppi <- lapply(1:length(icgc_permut_driver_marker_net_ppi), function(i) {
    tmp_row_index <- match(icgc_permut_driver_marker_net_ppi[[i]]$from, icgc_permut_driver_marker_net_genes[[i]]$STRING_id)
    icgc_permut_driver_marker_net_ppi[[i]]$from <- icgc_permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_row_index]
    
    tmp_column_index <- match(icgc_permut_driver_marker_net_ppi[[i]]$to, icgc_permut_driver_marker_net_genes[[i]]$STRING_id)
    icgc_permut_driver_marker_net_ppi[[i]]$to <- icgc_permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_column_index]
    
    tmp_tbl <- icgc_permut_driver_marker_net_ppi[[i]][,c(1,2)]
    tmp_tbl$type <- "PPI"
    
    ## Remove duplicates
    tmp_tbl <- tmp_tbl[-which(duplicated(paste(tmp_tbl$from,
                                               tmp_tbl$to, 
                                               sep = "_"))),]
    
    return(tmp_tbl)
  })
  
  names(icgc_permut_driver_marker_net_ppi) <- names(icgc_permut_driver_marker_net_genes)
  
  # Combine PPI with the original driver/marker net
  icgc_permut_driver_marker_net <- lapply(1:length(icgc_permut_driver_marker_net), function(i) {
    rbind(icgc_permut_driver_marker_net[[i]], icgc_permut_driver_marker_net_ppi[[i]])
  })
  
  names(icgc_permut_driver_marker_net) <- names(icgc_permut_driver_marker_net_genes)
  
  # Reconstruct the networks
  icgc_permut_driver_marker_net <- lapply(icgc_permut_driver_marker_net, function(i) {
    igraph::graph_from_data_frame(i, directed=FALSE)
  })
  
  ## Check the network layers
  table(igraph::E(icgc_permut_driver_marker_net$`Breast-AdenoCA`)$type)
  
  ###############################
  
  # Calculate the IVI values of genes within each net
  library(influential)
  icgc_permut_driver_marker_net_ivi <- lapply(icgc_permut_driver_marker_net, function(i) ivi(i, verbose = T))
  
  ###############################
  
  # Calculate the primitive gene driver scores
  icgc_permut_driver_marker_net_gene_tbl <- lapply(icgc_permut_driver_marker_net_ivi, function(i) {
    data.frame(gene = names(i), ivi = i)
  })
  
  ## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer
  
  icgc_permut_snv_mean_deletriousness <- lapply(icgc_permut_snv_flt_withCADD_summarized, function(i) {
    colMeans(i)
  })
  
  ## Generate the gene table for each net and calculate the primitive scores
  icgc_permut_driver_marker_net_gene_tbl <- lapply(1:length(icgc_permut_driver_marker_net_gene_tbl), function(i) {
    tmp_tbl <- icgc_permut_driver_marker_net_gene_tbl[[i]]
    tmp_tbl$mean_deletriousness <- 0
    match_index <- which(names(icgc_permut_snv_mean_deletriousness[[i]]) %in% rownames(tmp_tbl))
    tmp_tbl[names(icgc_permut_snv_mean_deletriousness[[i]])[match_index], "mean_deletriousness"] <- 
      icgc_permut_snv_mean_deletriousness[[i]][match_index]
    
    tmp_tbl$primitive_driver_score <- (tmp_tbl$ivi)*(tmp_tbl$mean_deletriousness) # node weights
    
    return(tmp_tbl)
    
  })
  
  names(icgc_permut_driver_marker_net_gene_tbl) <- names(icgc_permut_snv_mean_deletriousness)
  
  # add gene scores
  icgc_permut_driver_marker_net <- lapply(1:length(icgc_permut_driver_marker_net), function(i) {
    set_vertex_attr(graph = icgc_permut_driver_marker_net[[i]], name = "score", value = icgc_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score)
  })
  
  names(icgc_permut_driver_marker_net) <- names(icgc_permut_snv_mean_deletriousness)
  
  # Convert the un-weighted to a node-weighted network
  icgc_permut_driver_marker_net <- lapply(1:length(icgc_permut_driver_marker_net), function(i) {
    set_vertex_attr(graph = icgc_permut_driver_marker_net[[i]], name = "weight", value = (icgc_permut_driver_marker_net_gene_tbl[[i]]$ivi)^(icgc_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness))
  })
  
  names(icgc_permut_driver_marker_net) <- names(icgc_permut_snv_mean_deletriousness)
  
  ## Check the network node weights
  summary(igraph::V(icgc_permut_driver_marker_net$`Breast-AdenoCA`)$weight)
  
  # Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
  icgc_permut_driver_marker_net_gene_tbl <- lapply(1:length(icgc_permut_driver_marker_net), function(i) {
    tmp_tbl <- icgc_permut_driver_marker_net_gene_tbl[[i]]
    tmp_tbl$mean_neighborhood_weight <- sapply(1:nrow(tmp_tbl), function(j) {
      igraph::neighbors(graph = icgc_permut_driver_marker_net[[i]], 
                        v = V(icgc_permut_driver_marker_net[[i]])[j])$score %>% mean()
    })
    
    tmp_tbl$mediator_score <- tmp_tbl$primitive_driver_score * tmp_tbl$mean_neighborhood_weight
    
    return(tmp_tbl)
  })
  
  names(icgc_permut_driver_marker_net_gene_tbl) <- names(icgc_permut_snv_mean_deletriousness)
  
  # Detect communities/modules/clusters
  icgc_permut_driver_marker_net_modules <- lapply(icgc_permut_driver_marker_net, function(i) {
    set.seed(3847)
    igraph::cluster_leiden(
      graph = i,
      objective_function = "CPM",
      weights = NULL,
      resolution_parameter = 0.5,
      beta = 0.05,
      initial_membership = NULL,
      n_iterations = 100,
      vertex_weights = igraph::V(i)$weight
    )
  })
  
  sizes(icgc_permut_driver_marker_net_modules$`Breast-AdenoCA`) %>% as.vector() %>% summary()
  length(icgc_permut_driver_marker_net_modules$`Breast-AdenoCA`)
  
  ## Inspect and filter modules
  icgc_permut_driver_marker_net_modules_flt <- lapply(1:length(icgc_permut_driver_marker_net_modules), function(i) {
    
    ### Create a table of modules and their genes
    tmp_tbl <- data.frame(
      module = membership(icgc_permut_driver_marker_net_modules[[i]]) %>% as.integer(),
      gene = icgc_permut_driver_marker_net_modules[[i]]$names
    )
    
    ### Add if the gene is mutated or not (could be a risk gene or not)
    tmp_tbl$mutated <- 0
    tmp_tbl$mutated[which(icgc_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness > 0)] <- 1
    
    ### Separate modules
    tmp_tbl <- lapply(unique(tmp_tbl$module), function(j){
      subset(tmp_tbl, module == j)
    })
    
    names(tmp_tbl) <- sapply(tmp_tbl, function(m) {
      unique(m$module)
    })
    
    ### Remove modules that have less than 2 Risk genes or their size is less than 4
    tmp_tbl <- lapply(tmp_tbl, function(k) {
      if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
        k
      } else {
        NULL
      }
    })
    
    ### Remove NULL modules
    tmp_tbl[which(sapply(tmp_tbl, is.null))] <- NULL
    
    ### Calculate the module scores
    tmp_tbl <- lapply(tmp_tbl, function(l) {
      tmp_tbl4score <- cbind(l, mean_score = 0)
      tmp_tbl4score$mean_score <- mean(icgc_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score[which(icgc_permut_driver_marker_net_gene_tbl[[i]]$gene %in% l$gene)])
      tmp_tbl4score
    })
    
    tmp_tbl
  })
  
  names(icgc_permut_driver_marker_net_modules_flt) <- names(icgc_permut_driver_marker_net_modules)
  
  ## Sort the mean scores of all modules of each cancer type
  icgc_permut_driver_marker_net_modules_flt_scores <-
    lapply(icgc_permut_driver_marker_net_modules_flt, function(i) {
      tmp_tbl <-
        sapply(i, function(j) {unique(j$mean_score)})
      tmp_tbl <- rev(sort(tmp_tbl))
      tmp_tbl
    })
  
  ### View the number of first ranked modules
  sapply(icgc_permut_driver_marker_net_modules_flt_scores, function(i) {i[1]}) %>% as.data.frame() %>% View()
  
  icgc_permut_driver_marker_net_modules_flt$`Breast-AdenoCA`$`7`
  
  ########
  
  # Top 20 gene 
  icgc_permut_driver_marker_net_gene_tbl_top20 <-
    lapply(icgc_permut_driver_marker_net_gene_tbl, function(i) {
      slice_max(i, n = 20, order_by = primitive_driver_score) %>% select(gene)
    })
  
  icgc_permut_driver_marker_net_gene_tbl_top20 <- do.call(cbind, icgc_permut_driver_marker_net_gene_tbl_top20)
  colnames(icgc_permut_driver_marker_net_gene_tbl_top20) <- names(icgc_permut_driver_marker_net_gene_tbl)
  
  ###########################################
  ###########################################
  ###########################################
  
  ## Save the gene data and modules as supplementary files
  
  tmp_icgc_permut_driver_marker_net_gene_tbl <- icgc_permut_driver_marker_net_gene_tbl

  for(i in 1:5) {
    colnames(tmp_icgc_permut_driver_marker_net_gene_tbl[[i]]) <- c("Gene", "IVI", "Mean_deletriousness_score", "Driver_score", "Mean_neighborhood_weight", "Driverness_mediator_score")
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_icgc_permut_driver_marker_net_gene_tbl, file = "Results/Evaluations/icgc_permutation/icgc_permut_driver_marker_net_gene_tbl.xlsx")
  
  rm(tmp_icgc_permut_driver_marker_net_gene_tbl)
  
  ############################
  
  ### save the modules
  
  tmp_icgc_permut_driver_marker_net_modules_flt <- icgc_permut_driver_marker_net_modules_flt

  tmp_icgc_permut_driver_marker_net_modules_flt_scores <- icgc_permut_driver_marker_net_modules_flt_scores
  
  tmp_icgc_permut_driver_marker_net_modules_flt <- lapply(tmp_icgc_permut_driver_marker_net_modules_flt, function(i) {
    
    tmp_data <- i
    tmp_data <- tmp_data[rev(order(sapply(tmp_data, function(j) unique(j$mean_score))))]
    
    tmp_data
    
  })
  
  for(i in 1:5) {
    tmp_icgc_permut_driver_marker_net_modules_flt[[i]] <-
      lapply(1:length(tmp_icgc_permut_driver_marker_net_modules_flt[[i]]), function(j) {
        tmp_data <- tmp_icgc_permut_driver_marker_net_modules_flt[[i]][[j]]
        colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
        tmp_data[tmp_data == 1] <- "Yes"
        tmp_data[tmp_data == 0] <- "No"
        tmp_data$Module <- paste0("Module number ", j)
        tmp_data <- data.frame(Module = unique(tmp_data$Module),
                               Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                               Mean_score = unique(tmp_data$Mean_score))
        tmp_data
      })
    tmp_icgc_permut_driver_marker_net_modules_flt[[i]] <- do.call(rbind, tmp_icgc_permut_driver_marker_net_modules_flt[[i]])
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_icgc_permut_driver_marker_net_modules_flt, file = "Results/Evaluations/icgc_permutation/icgc_permut_driver_marker_net_modules_flt.xlsx")
  
  rm(tmp_icgc_permut_driver_marker_net_modules_flt, tmp_icgc_permut_driver_marker_net_modules_flt_scores)
  
  ###############################################################################################################
  ###############################################################################################################
  ###############################################################################################################
  
  # Test the pipeline based on onlyTx_icgc_permutation (shuffling only expression data)
  
  ## Prepare somatic snv data
  
  onlyTx_icgc_permut_snv_flt_cor <- icgc_snv_data_flt_cor
  
  ##################
  
  ## Prepare expression data
  onlyTx_icgc_permut_rnaseq_expr_list <- icgc_tx_data
  
  ## Shuffle the data
  onlyTx_icgc_permut_rnaseq_expr_list <- lapply(onlyTx_icgc_permut_rnaseq_expr_list, function(i) {
    tmp_colnames <- colnames(i)
    
    ## shuffle
    tmp_data <- apply(i, 1, function(j) {
      set.seed(1234)
      j[sample(length(j))]
    })
    tmp_data <- as.data.frame(t(tmp_data))
    colnames(tmp_data) <- tmp_colnames
    
    ## re-shuffle
    tmp_data <- apply(i, 1, function(j) {
      set.seed(5678)
      j[sample(length(j))]
    })
    tmp_data <- as.data.frame(t(tmp_data))
    colnames(tmp_data) <- tmp_colnames
    
    tmp_data
  })
  
  ## Perform correlation analysis
  onlyTx_icgc_permut_rnaseq_expr_cor <- lapply(onlyTx_icgc_permut_rnaseq_expr_list, function(i){
    
    tmp_cor <- fcor(data = i, 
                    method = "spearman", 
                    mutualRank = TRUE, pvalue = FALSE, flat = TRUE)
    
    tmp_cor <- subset(tmp_cor, mr < 20)
    tmp_cor
  })
  
  ###################
  
  # Prepare the tables for merging with other data
  onlyTx_icgc_permut_snv_flt_cor <- lapply(onlyTx_icgc_permut_snv_flt_cor, function(i) {
    tmp_tbl <- i[,c(1,2)]
    colnames(tmp_tbl) <- c("from", "to")
    tmp_tbl$type <- "co-deletriousness"
    tmp_tbl
  })
  
  onlyTx_icgc_permut_rnaseq_expr_cor <- lapply(onlyTx_icgc_permut_rnaseq_expr_cor, function(i) {
    tmp_tbl <- i[,c(1,2)]
    colnames(tmp_tbl) <- c("from", "to")
    tmp_tbl$type <- "co-expression"
    tmp_tbl
  })
  
  # Merge the similarity tables
  onlyTx_icgc_permut_driver_marker_net <- onlyTx_icgc_permut_snv_flt_cor
  
  # Make the orders the same
  onlyTx_icgc_permut_rnaseq_expr_cor <- onlyTx_icgc_permut_rnaseq_expr_cor[names(onlyTx_icgc_permut_snv_flt_cor)]
  
  # Add the co-expression data to the tables
  
  for(i in 1:length(onlyTx_icgc_permut_driver_marker_net)) {
    
    genes <- onlyTx_icgc_permut_driver_marker_net[[i]][,c(1,2)] %>% 
      unlist() %>% 
      unname() %>% 
      unique()
    
    first_level_coex_index <- c(which(onlyTx_icgc_permut_rnaseq_expr_cor[[i]]$from %in% genes), 
                                which(onlyTx_icgc_permut_rnaseq_expr_cor[[i]]$to %in% genes)) %>% unique()
    
    first_level_coex_genes <- c(genes, 
                                onlyTx_icgc_permut_rnaseq_expr_cor[[i]]$from[first_level_coex_index],
                                onlyTx_icgc_permut_rnaseq_expr_cor[[i]]$to[first_level_coex_index]) %>% unique()
    
    first_two_levels_coex_index <- c(which(onlyTx_icgc_permut_rnaseq_expr_cor[[i]]$from %in% first_level_coex_genes), 
                                     which(onlyTx_icgc_permut_rnaseq_expr_cor[[i]]$to %in% first_level_coex_genes)) %>% unique()
    
    onlyTx_icgc_permut_driver_marker_net[[i]] <- rbind(onlyTx_icgc_permut_driver_marker_net[[i]], onlyTx_icgc_permut_rnaseq_expr_cor[[i]][first_two_levels_coex_index,])
    
  }
  
  # Add PPI data to the driver/marker net
  
  # Map gene names to stringDB ids for PPI analysis (STRING  11.5)
  
  ## get the gene names in each risk net
  onlyTx_icgc_permut_driver_marker_net_genes <- lapply(onlyTx_icgc_permut_driver_marker_net, function(i) {
    data.frame(gene_symbol = i[,c(1,2)] %>% 
                 unlist() %>% 
                 unname() %>% 
                 unique())
  })
  
  ## getSTRINGdb for human
  ### WARNING: You didn't specify a version of the STRING database to use. Hence we will use STRING  11.5
  string_db <- STRINGdb$new(species=9606, version = "11.5", 
                            score_threshold=900, # very high confidence (>0.9)
                            input_directory=""
  )
  
  onlyTx_icgc_permut_driver_marker_net_genes <- lapply(onlyTx_icgc_permut_driver_marker_net_genes, function(i) {
    
    string_db$map(my_data_frame = i, 
                  my_data_frame_id_col_names = "gene_symbol", 
                  takeFirst = TRUE, removeUnmappedRows = TRUE)
    
  })
  
  # Get the PPIs
  onlyTx_icgc_permut_driver_marker_net_ppi <- lapply(onlyTx_icgc_permut_driver_marker_net_genes, function(i) {
    
    string_db$get_interactions(string_ids = i$STRING_id)
    
  })
  
  # Convert stringdb ids to gene symbols
  onlyTx_icgc_permut_driver_marker_net_ppi <- lapply(1:length(onlyTx_icgc_permut_driver_marker_net_ppi), function(i) {
    tmp_row_index <- match(onlyTx_icgc_permut_driver_marker_net_ppi[[i]]$from, onlyTx_icgc_permut_driver_marker_net_genes[[i]]$STRING_id)
    onlyTx_icgc_permut_driver_marker_net_ppi[[i]]$from <- onlyTx_icgc_permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_row_index]
    
    tmp_column_index <- match(onlyTx_icgc_permut_driver_marker_net_ppi[[i]]$to, onlyTx_icgc_permut_driver_marker_net_genes[[i]]$STRING_id)
    onlyTx_icgc_permut_driver_marker_net_ppi[[i]]$to <- onlyTx_icgc_permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_column_index]
    
    tmp_tbl <- onlyTx_icgc_permut_driver_marker_net_ppi[[i]][,c(1,2)]
    tmp_tbl$type <- "PPI"
    
    ## Remove duplicates
    tmp_tbl <- tmp_tbl[-which(duplicated(paste(tmp_tbl$from,
                                               tmp_tbl$to, 
                                               sep = "_"))),]
    
    return(tmp_tbl)
  })
  
  names(onlyTx_icgc_permut_driver_marker_net_ppi) <- names(onlyTx_icgc_permut_driver_marker_net_genes)
  
  # Combine PPI with the original driver/marker net
  onlyTx_icgc_permut_driver_marker_net <- lapply(1:length(onlyTx_icgc_permut_driver_marker_net), function(i) {
    rbind(onlyTx_icgc_permut_driver_marker_net[[i]], onlyTx_icgc_permut_driver_marker_net_ppi[[i]])
  })
  
  names(onlyTx_icgc_permut_driver_marker_net) <- names(onlyTx_icgc_permut_driver_marker_net_genes)
  
  # Reconstruct the networks
  onlyTx_icgc_permut_driver_marker_net <- lapply(onlyTx_icgc_permut_driver_marker_net, function(i) {
    igraph::graph_from_data_frame(i, directed=FALSE)
  })
  
  ## Check the network layers
  table(igraph::E(onlyTx_icgc_permut_driver_marker_net$`Breast-AdenoCA`)$type)
  
  ###############################
  
  # Calculate the IVI values of genes within each net
  library(influential)
  onlyTx_icgc_permut_driver_marker_net_ivi <- lapply(onlyTx_icgc_permut_driver_marker_net, function(i) ivi(i, verbose = T))
  
  ###############################
  
  # Calculate the primitive gene driver scores
  onlyTx_icgc_permut_driver_marker_net_gene_tbl <- lapply(onlyTx_icgc_permut_driver_marker_net_ivi, function(i) {
    data.frame(gene = names(i), ivi = i)
  })
  
  ## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer
  
  onlyTx_icgc_permut_snv_mean_deletriousness <- lapply(icgc_snv_data_flt_withCADD_summarized, function(i) {
    colMeans(i)
  })
  
  ## Generate the gene table for each net and calculate the primitive scores
  onlyTx_icgc_permut_driver_marker_net_gene_tbl <- lapply(1:length(onlyTx_icgc_permut_driver_marker_net_gene_tbl), function(i) {
    tmp_tbl <- onlyTx_icgc_permut_driver_marker_net_gene_tbl[[i]]
    tmp_tbl$mean_deletriousness <- 0
    match_index <- which(names(onlyTx_icgc_permut_snv_mean_deletriousness[[i]]) %in% rownames(tmp_tbl))
    tmp_tbl[names(onlyTx_icgc_permut_snv_mean_deletriousness[[i]])[match_index], "mean_deletriousness"] <- 
      onlyTx_icgc_permut_snv_mean_deletriousness[[i]][match_index]
    
    tmp_tbl$primitive_driver_score <- (tmp_tbl$ivi)*(tmp_tbl$mean_deletriousness) # node weights
    
    return(tmp_tbl)
    
  })
  
  names(onlyTx_icgc_permut_driver_marker_net_gene_tbl) <- names(onlyTx_icgc_permut_snv_mean_deletriousness)
  
  # add gene scores
  onlyTx_icgc_permut_driver_marker_net <- lapply(1:length(onlyTx_icgc_permut_driver_marker_net), function(i) {
    set_vertex_attr(graph = onlyTx_icgc_permut_driver_marker_net[[i]], name = "score", value = onlyTx_icgc_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score)
  })
  
  names(onlyTx_icgc_permut_driver_marker_net) <- names(onlyTx_icgc_permut_snv_mean_deletriousness)
  
  # Convert the un-weighted to a node-weighted network
  onlyTx_icgc_permut_driver_marker_net <- lapply(1:length(onlyTx_icgc_permut_driver_marker_net), function(i) {
    set_vertex_attr(graph = onlyTx_icgc_permut_driver_marker_net[[i]], name = "weight", value = (onlyTx_icgc_permut_driver_marker_net_gene_tbl[[i]]$ivi)^(onlyTx_icgc_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness))
  })
  
  names(onlyTx_icgc_permut_driver_marker_net) <- names(onlyTx_icgc_permut_snv_mean_deletriousness)
  
  ## Check the network node weights
  summary(igraph::V(onlyTx_icgc_permut_driver_marker_net$`Breast-AdenoCA`)$weight)
  
  # Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
  onlyTx_icgc_permut_driver_marker_net_gene_tbl <- lapply(1:length(onlyTx_icgc_permut_driver_marker_net), function(i) {
    tmp_tbl <- onlyTx_icgc_permut_driver_marker_net_gene_tbl[[i]]
    tmp_tbl$mean_neighborhood_weight <- sapply(1:nrow(tmp_tbl), function(j) {
      igraph::neighbors(graph = onlyTx_icgc_permut_driver_marker_net[[i]], 
                        v = V(onlyTx_icgc_permut_driver_marker_net[[i]])[j])$score %>% mean()
    })
    
    tmp_tbl$mediator_score <- tmp_tbl$primitive_driver_score * tmp_tbl$mean_neighborhood_weight
    
    return(tmp_tbl)
  })
  
  names(onlyTx_icgc_permut_driver_marker_net_gene_tbl) <- names(onlyTx_icgc_permut_snv_mean_deletriousness)
  
  # Detect communities/modules/clusters
  onlyTx_icgc_permut_driver_marker_net_modules <- lapply(onlyTx_icgc_permut_driver_marker_net, function(i) {
    set.seed(3847)
    igraph::cluster_leiden(
      graph = i,
      objective_function = "CPM",
      weights = NULL,
      resolution_parameter = 0.5,
      beta = 0.05,
      initial_membership = NULL,
      n_iterations = 100,
      vertex_weights = igraph::V(i)$weight
    )
  })
  
  sizes(onlyTx_icgc_permut_driver_marker_net_modules$`Breast-AdenoCA`) %>% as.vector() %>% summary()
  length(onlyTx_icgc_permut_driver_marker_net_modules$`Breast-AdenoCA`)
  
  ## Inspect and filter modules
  onlyTx_icgc_permut_driver_marker_net_modules_flt <- lapply(1:length(onlyTx_icgc_permut_driver_marker_net_modules), function(i) {
    
    ### Create a table of modules and their genes
    tmp_tbl <- data.frame(
      module = membership(onlyTx_icgc_permut_driver_marker_net_modules[[i]]) %>% as.integer(),
      gene = onlyTx_icgc_permut_driver_marker_net_modules[[i]]$names
    )
    
    ### Add if the gene is mutated or not (could be a risk gene or not)
    tmp_tbl$mutated <- 0
    tmp_tbl$mutated[which(onlyTx_icgc_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness > 0)] <- 1
    
    ### Separate modules
    tmp_tbl <- lapply(unique(tmp_tbl$module), function(j){
      subset(tmp_tbl, module == j)
    })
    
    names(tmp_tbl) <- sapply(tmp_tbl, function(m) {
      unique(m$module)
    })
    
    ### Remove modules that have less than 2 Risk genes or their size is less than 4
    tmp_tbl <- lapply(tmp_tbl, function(k) {
      if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
        k
      } else {
        NULL
      }
    })
    
    ### Remove NULL modules
    tmp_tbl[which(sapply(tmp_tbl, is.null))] <- NULL
    
    ### Calculate the module scores
    tmp_tbl <- lapply(tmp_tbl, function(l) {
      tmp_tbl4score <- cbind(l, mean_score = 0)
      tmp_tbl4score$mean_score <- mean(onlyTx_icgc_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score[which(onlyTx_icgc_permut_driver_marker_net_gene_tbl[[i]]$gene %in% l$gene)])
      tmp_tbl4score
    })
    
    tmp_tbl
  })
  
  names(onlyTx_icgc_permut_driver_marker_net_modules_flt) <- names(onlyTx_icgc_permut_driver_marker_net_modules)
  
  ## Sort the mean scores of all modules of each cancer type
  onlyTx_icgc_permut_driver_marker_net_modules_flt_scores <-
    lapply(onlyTx_icgc_permut_driver_marker_net_modules_flt, function(i) {
      tmp_tbl <-
        sapply(i, function(j) {unique(j$mean_score)})
      tmp_tbl <- rev(sort(tmp_tbl))
      tmp_tbl
    })
  
  ### View the number of first ranked modules
  sapply(onlyTx_icgc_permut_driver_marker_net_modules_flt_scores, function(i) {i[1]}) %>% as.data.frame() %>% View()
  
  onlyTx_icgc_permut_driver_marker_net_modules_flt$`Breast-AdenoCA`$`1317`
  
  ########
  
  # Top 20 gene 
  onlyTx_icgc_permut_driver_marker_net_gene_tbl_top20 <-
    lapply(onlyTx_icgc_permut_driver_marker_net_gene_tbl, function(i) {
      slice_max(i, n = 20, order_by = primitive_driver_score) %>% select(gene)
    })
  
  onlyTx_icgc_permut_driver_marker_net_gene_tbl_top20 <- do.call(cbind, onlyTx_icgc_permut_driver_marker_net_gene_tbl_top20)
  colnames(onlyTx_icgc_permut_driver_marker_net_gene_tbl_top20) <- names(onlyTx_icgc_permut_driver_marker_net_gene_tbl)
  
  ###########################################
  ###########################################
  ###########################################
  
  ## Save the gene data and modules as supplementary files
  
  tmp_onlyTx_icgc_permut_driver_marker_net_gene_tbl <- onlyTx_icgc_permut_driver_marker_net_gene_tbl
  
  for(i in 1:5) {
    colnames(tmp_onlyTx_icgc_permut_driver_marker_net_gene_tbl[[i]]) <- c("Gene", "IVI", "Mean_deletriousness_score", "Driver_score", "Mean_neighborhood_weight", "Driverness_mediator_score")
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_onlyTx_icgc_permut_driver_marker_net_gene_tbl, file = "Results/Evaluations/onlyTx_icgc_permutation/onlyTx_icgc_permut_driver_marker_net_gene_tbl.xlsx")
  
  rm(tmp_onlyTx_icgc_permut_driver_marker_net_gene_tbl)
  
  ############################
  
  ### save the modules
  
  tmp_onlyTx_icgc_permut_driver_marker_net_modules_flt <- onlyTx_icgc_permut_driver_marker_net_modules_flt
  
  tmp_onlyTx_icgc_permut_driver_marker_net_modules_flt_scores <- onlyTx_icgc_permut_driver_marker_net_modules_flt_scores
  
  tmp_onlyTx_icgc_permut_driver_marker_net_modules_flt <- lapply(tmp_onlyTx_icgc_permut_driver_marker_net_modules_flt, function(i) {
    
    tmp_data <- i
    tmp_data <- tmp_data[rev(order(sapply(tmp_data, function(j) unique(j$mean_score))))]
    
    tmp_data
    
  })
  
  for(i in 1:5) {
    tmp_onlyTx_icgc_permut_driver_marker_net_modules_flt[[i]] <-
      lapply(1:length(tmp_onlyTx_icgc_permut_driver_marker_net_modules_flt[[i]]), function(j) {
        tmp_data <- tmp_onlyTx_icgc_permut_driver_marker_net_modules_flt[[i]][[j]]
        colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
        tmp_data[tmp_data == 1] <- "Yes"
        tmp_data[tmp_data == 0] <- "No"
        tmp_data$Module <- paste0("Module number ", j)
        tmp_data <- data.frame(Module = unique(tmp_data$Module),
                               Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                               Mean_score = unique(tmp_data$Mean_score))
        tmp_data
      })
    tmp_onlyTx_icgc_permut_driver_marker_net_modules_flt[[i]] <- do.call(rbind, tmp_onlyTx_icgc_permut_driver_marker_net_modules_flt[[i]])
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_onlyTx_icgc_permut_driver_marker_net_modules_flt, file = "Results/Evaluations/onlyTx_icgc_permutation/onlyTx_icgc_permut_driver_marker_net_modules_flt.xlsx")
  
  rm(tmp_onlyTx_icgc_permut_driver_marker_net_modules_flt, tmp_onlyTx_icgc_permut_driver_marker_net_modules_flt_scores)
  
  ###############################################################################################################
  ###############################################################################################################
  ###############################################################################################################
  
  # Test the pipeline based on icgc_permutation (shuffling everything) 
  
  ## Prepare cor data
  
  ### Setting snv cor to the ones derived from icgc_permutated data
  plusPPI_icgc_permut_snv_flt_cor <- icgc_permut_snv_flt_cor
  
  ### Setting expr cor to the ones derived from icgc_permutated data
  plusPPI_icgc_permut_rnaseq_expr_cor <- icgc_permut_rnaseq_expr_cor

  # Merge the similarity tables
  plusPPI_icgc_permut_driver_marker_net <- plusPPI_icgc_permut_snv_flt_cor
  
  # Add the co-expression data to the tables
  
  for(i in 1:length(plusPPI_icgc_permut_driver_marker_net)) {
    
    genes <- plusPPI_icgc_permut_driver_marker_net[[i]][,c(1,2)] %>% 
      unlist() %>% 
      unname() %>% 
      unique()
    
    first_level_coex_index <- c(which(plusPPI_icgc_permut_rnaseq_expr_cor[[i]]$from %in% genes), 
                                which(plusPPI_icgc_permut_rnaseq_expr_cor[[i]]$to %in% genes)) %>% unique()
    
    first_level_coex_genes <- c(genes, 
                                plusPPI_icgc_permut_rnaseq_expr_cor[[i]]$from[first_level_coex_index],
                                plusPPI_icgc_permut_rnaseq_expr_cor[[i]]$to[first_level_coex_index]) %>% unique()
    
    first_two_levels_coex_index <- c(which(plusPPI_icgc_permut_rnaseq_expr_cor[[i]]$from %in% first_level_coex_genes), 
                                     which(plusPPI_icgc_permut_rnaseq_expr_cor[[i]]$to %in% first_level_coex_genes)) %>% unique()
    
    plusPPI_icgc_permut_driver_marker_net[[i]] <- rbind(plusPPI_icgc_permut_driver_marker_net[[i]], plusPPI_icgc_permut_rnaseq_expr_cor[[i]][first_two_levels_coex_index,])
    
  }
  
  # icgc_permute PPI data
  plusPPI_icgc_permut_driver_marker_net_ppi <- icgc_permut_driver_marker_net_ppi
  plusPPI_icgc_permut_driver_marker_net_ppi <- lapply(plusPPI_icgc_permut_driver_marker_net_ppi, function(i) {
    tmp_data <- i
    
    ## Remove duplicates (already removed)
    
    ## Shuffle the data
    set.seed(1234)
    tmp_data$from <- tmp_data$from[sample(nrow(tmp_data))]
    tmp_data
  })
  
  # Add PPI data to the driver/marker net
  
  # Combine PPI with the original driver/marker net
  plusPPI_icgc_permut_driver_marker_net <- lapply(1:length(plusPPI_icgc_permut_driver_marker_net), function(i) {
    rbind(plusPPI_icgc_permut_driver_marker_net[[i]], plusPPI_icgc_permut_driver_marker_net_ppi[[i]])
  })
  
  names(plusPPI_icgc_permut_driver_marker_net) <- names(plusPPI_icgc_permut_driver_marker_net_ppi)
  
  # Reconstruct the networks
  plusPPI_icgc_permut_driver_marker_net <- lapply(plusPPI_icgc_permut_driver_marker_net, function(i) {
    igraph::graph_from_data_frame(i, directed=FALSE)
  })
  
  ## Check the network layers
  table(igraph::E(plusPPI_icgc_permut_driver_marker_net$`Breast-AdenoCA`)$type)
  
  ###############################
  
  # Calculate the IVI values of genes within each net
  library(influential)
  plusPPI_icgc_permut_driver_marker_net_ivi <- lapply(plusPPI_icgc_permut_driver_marker_net, function(i) ivi(i, verbose = T))
  
  ###############################
  
  # Calculate the primitive gene driver scores
  plusPPI_icgc_permut_driver_marker_net_gene_tbl <- lapply(plusPPI_icgc_permut_driver_marker_net_ivi, function(i) {
    data.frame(gene = names(i), ivi = i)
  })
  
  ## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer
  
  plusPPI_icgc_permut_snv_mean_deletriousness <- lapply(icgc_permut_snv_flt_withCADD_summarized, function(i) {
    colMeans(i)
  })
  
  ## Generate the gene table for each net and calculate the primitive scores
  plusPPI_icgc_permut_driver_marker_net_gene_tbl <- lapply(1:length(plusPPI_icgc_permut_driver_marker_net_gene_tbl), function(i) {
    tmp_tbl <- plusPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]
    tmp_tbl$mean_deletriousness <- 0
    match_index <- which(names(plusPPI_icgc_permut_snv_mean_deletriousness[[i]]) %in% rownames(tmp_tbl))
    tmp_tbl[names(plusPPI_icgc_permut_snv_mean_deletriousness[[i]])[match_index], "mean_deletriousness"] <- 
      plusPPI_icgc_permut_snv_mean_deletriousness[[i]][match_index]
    
    tmp_tbl$primitive_driver_score <- (tmp_tbl$ivi)*(tmp_tbl$mean_deletriousness) # node weights
    
    return(tmp_tbl)
    
  })
  
  names(plusPPI_icgc_permut_driver_marker_net_gene_tbl) <- names(plusPPI_icgc_permut_snv_mean_deletriousness)
  
  # add gene scores
  plusPPI_icgc_permut_driver_marker_net <- lapply(1:length(plusPPI_icgc_permut_driver_marker_net), function(i) {
    set_vertex_attr(graph = plusPPI_icgc_permut_driver_marker_net[[i]], name = "score", value = plusPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score)
  })
  
  names(plusPPI_icgc_permut_driver_marker_net) <- names(plusPPI_icgc_permut_snv_mean_deletriousness)
  
  # Convert the un-weighted to a node-weighted network
  plusPPI_icgc_permut_driver_marker_net <- lapply(1:length(plusPPI_icgc_permut_driver_marker_net), function(i) {
    set_vertex_attr(graph = plusPPI_icgc_permut_driver_marker_net[[i]], name = "weight", value = (plusPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]$ivi)^plusPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness)
  })
  
  names(plusPPI_icgc_permut_driver_marker_net) <- names(plusPPI_icgc_permut_snv_mean_deletriousness)
  
  ## Check the network node weights
  summary(igraph::V(plusPPI_icgc_permut_driver_marker_net$`Breast-AdenoCA`)$weight)
  
  # Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
  plusPPI_icgc_permut_driver_marker_net_gene_tbl <- lapply(1:length(plusPPI_icgc_permut_driver_marker_net), function(i) {
    tmp_tbl <- plusPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]
    tmp_tbl$mean_neighborhood_weight <- sapply(1:nrow(tmp_tbl), function(j) {
      igraph::neighbors(graph = plusPPI_icgc_permut_driver_marker_net[[i]], 
                        v = V(plusPPI_icgc_permut_driver_marker_net[[i]])[j])$score %>% mean()
    })
    
    tmp_tbl$mediator_score <- tmp_tbl$primitive_driver_score * tmp_tbl$mean_neighborhood_weight
    
    return(tmp_tbl)
  })
  
  names(plusPPI_icgc_permut_driver_marker_net_gene_tbl) <- names(plusPPI_icgc_permut_snv_mean_deletriousness)
  
  # Detect communities/modules/clusters
  plusPPI_icgc_permut_driver_marker_net_modules <- lapply(plusPPI_icgc_permut_driver_marker_net, function(i) {
    set.seed(3847)
    igraph::cluster_leiden(
      graph = i,
      objective_function = "CPM",
      weights = NULL,
      resolution_parameter = 0.5,
      beta = 0.05,
      initial_membership = NULL,
      n_iterations = 100,
      vertex_weights = igraph::V(i)$weight
    )
  })
  
  sizes(plusPPI_icgc_permut_driver_marker_net_modules$`Breast-AdenoCA`) %>% as.vector() %>% summary()
  length(plusPPI_icgc_permut_driver_marker_net_modules$`Breast-AdenoCA`)
  
  ## Inspect and filter modules
  plusPPI_icgc_permut_driver_marker_net_modules_flt <- lapply(1:length(plusPPI_icgc_permut_driver_marker_net_modules), function(i) {
    
    ### Create a table of modules and their genes
    tmp_tbl <- data.frame(
      module = membership(plusPPI_icgc_permut_driver_marker_net_modules[[i]]) %>% as.integer(),
      gene = plusPPI_icgc_permut_driver_marker_net_modules[[i]]$names
    )
    
    ### Add if the gene is mutated or not (could be a risk gene or not)
    tmp_tbl$mutated <- 0
    tmp_tbl$mutated[which(plusPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness > 0)] <- 1
    
    ### Separate modules
    tmp_tbl <- lapply(unique(tmp_tbl$module), function(j){
      subset(tmp_tbl, module == j)
    })
    
    names(tmp_tbl) <- sapply(tmp_tbl, function(m) {
      unique(m$module)
    })
    
    ### Remove modules that have less than 2 Risk genes or their size is less than 4
    tmp_tbl <- lapply(tmp_tbl, function(k) {
      if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
        k
      } else {
        NULL
      }
    })
    
    ### Remove NULL modules
    tmp_tbl[which(sapply(tmp_tbl, is.null))] <- NULL
    
    ### Calculate the module scores
    tmp_tbl <- lapply(tmp_tbl, function(l) {
      tmp_tbl4score <- cbind(l, mean_score = 0)
      tmp_tbl4score$mean_score <- mean(plusPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score[which(plusPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]$gene %in% l$gene)])
      tmp_tbl4score
    })
    
    tmp_tbl
  })
  
  names(plusPPI_icgc_permut_driver_marker_net_modules_flt) <- names(plusPPI_icgc_permut_driver_marker_net_modules)
  
  sapply(plusPPI_icgc_permut_driver_marker_net_modules_flt, length)

  ## Sort the mean scores of all modules of each cancer type
  plusPPI_icgc_permut_driver_marker_net_modules_flt_scores <-
    lapply(plusPPI_icgc_permut_driver_marker_net_modules_flt, function(i) {
      tmp_tbl <-
        sapply(i, function(j) {unique(j$mean_score)})
      tmp_tbl <- rev(sort(tmp_tbl))
      tmp_tbl
    })
  
  ### View the number of first ranked modules
  sapply(plusPPI_icgc_permut_driver_marker_net_modules_flt_scores, function(i) {i[1]}) %>% as.data.frame() %>% View()
  
  plusPPI_icgc_permut_driver_marker_net_modules_flt$`Breast-AdenoCA`$`1690` %>% View()
  
  ########
  
  # Top 20 gene 
  plusPPI_icgc_permut_driver_marker_net_gene_tbl_top20 <-
    lapply(plusPPI_icgc_permut_driver_marker_net_gene_tbl, function(i) {
      slice_max(i, n = 20, order_by = primitive_driver_score) %>% select(gene)
    })
  
  plusPPI_icgc_permut_driver_marker_net_gene_tbl_top20 <- do.call(cbind, plusPPI_icgc_permut_driver_marker_net_gene_tbl_top20)
  colnames(plusPPI_icgc_permut_driver_marker_net_gene_tbl_top20) <- names(plusPPI_icgc_permut_driver_marker_net_gene_tbl)
  
  ###########################################
  ###########################################
  ###########################################
  
  ## Save the gene data and modules as supplementary files
  
  tmp_plusPPI_icgc_permut_driver_marker_net_gene_tbl <- plusPPI_icgc_permut_driver_marker_net_gene_tbl

  for(i in 1:5) {
    colnames(tmp_plusPPI_icgc_permut_driver_marker_net_gene_tbl[[i]]) <- c("Gene", "IVI", "Mean_deletriousness_score", "Driver_score", "Mean_neighborhood_weight", "Driverness_mediator_score")
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_plusPPI_icgc_permut_driver_marker_net_gene_tbl, file = "Results/Evaluations/icgc_permutation/plusPPI_icgc_permut_driver_marker_net_gene_tbl.xlsx")
  
  rm(tmp_plusPPI_icgc_permut_driver_marker_net_gene_tbl)
  
  ############################
  
  ### save the modules
  
  tmp_plusPPI_icgc_permut_driver_marker_net_modules_flt <- plusPPI_icgc_permut_driver_marker_net_modules_flt

  tmp_plusPPI_icgc_permut_driver_marker_net_modules_flt_scores <- plusPPI_icgc_permut_driver_marker_net_modules_flt_scores
  
  tmp_plusPPI_icgc_permut_driver_marker_net_modules_flt <- lapply(tmp_plusPPI_icgc_permut_driver_marker_net_modules_flt, function(i) {
    
    tmp_data <- i
    tmp_data <- tmp_data[rev(order(sapply(tmp_data, function(j) unique(j$mean_score))))]
    
    tmp_data
    
  })
  
  for(i in 1:5) {
    tmp_plusPPI_icgc_permut_driver_marker_net_modules_flt[[i]] <-
      lapply(1:length(tmp_plusPPI_icgc_permut_driver_marker_net_modules_flt[[i]]), function(j) {
        tmp_data <- tmp_plusPPI_icgc_permut_driver_marker_net_modules_flt[[i]][[j]]
        colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
        tmp_data[tmp_data == 1] <- "Yes"
        tmp_data[tmp_data == 0] <- "No"
        tmp_data$Module <- paste0("Module number ", j)
        tmp_data <- data.frame(Module = unique(tmp_data$Module),
                               Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                               Mean_score = unique(tmp_data$Mean_score))
        tmp_data
      })
    tmp_plusPPI_icgc_permut_driver_marker_net_modules_flt[[i]] <- do.call(rbind, tmp_plusPPI_icgc_permut_driver_marker_net_modules_flt[[i]])
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_plusPPI_icgc_permut_driver_marker_net_modules_flt, file = "Results/Evaluations/icgc_permutation/plusPPI_icgc_permut_driver_marker_net_modules_flt.xlsx")
  
  rm(tmp_plusPPI_icgc_permut_driver_marker_net_modules_flt, tmp_plusPPI_icgc_permut_driver_marker_net_modules_flt_scores)
  
  ###############################################################################################################
  ###############################################################################################################
  ###############################################################################################################
  
  # Test the pipeline based on icgc_permutation (shuffling everything but deleteriousness data) 
  
  ## Prepare cor data
  
  ### Setting snv cor to the ones derived from original data
  noSNV_icgc_permut_snv_flt_cor <- icgc_snv_data_flt_cor
  
  ### Setting expr cor to the ones derived from icgc_permutated data
  noSNV_icgc_permut_rnaseq_expr_cor <- icgc_permut_rnaseq_expr_cor

  # Merge the similarity tables
  noSNV_icgc_permut_driver_marker_net <- noSNV_icgc_permut_snv_flt_cor
  
  # Add the co-expression data to the tables
  
  for(i in 1:length(noSNV_icgc_permut_driver_marker_net)) {
    
    genes <- noSNV_icgc_permut_driver_marker_net[[i]][,c(1,2)] %>% 
      unlist() %>% 
      unname() %>% 
      unique()
    
    first_level_coex_index <- c(which(noSNV_icgc_permut_rnaseq_expr_cor[[i]]$from %in% genes), 
                                which(noSNV_icgc_permut_rnaseq_expr_cor[[i]]$to %in% genes)) %>% unique()
    
    first_level_coex_genes <- c(genes, 
                                noSNV_icgc_permut_rnaseq_expr_cor[[i]]$from[first_level_coex_index],
                                noSNV_icgc_permut_rnaseq_expr_cor[[i]]$to[first_level_coex_index]) %>% unique()
    
    first_two_levels_coex_index <- c(which(noSNV_icgc_permut_rnaseq_expr_cor[[i]]$from %in% first_level_coex_genes), 
                                     which(noSNV_icgc_permut_rnaseq_expr_cor[[i]]$to %in% first_level_coex_genes)) %>% unique()
    
    noSNV_icgc_permut_driver_marker_net[[i]] <- rbind(noSNV_icgc_permut_driver_marker_net[[i]], noSNV_icgc_permut_rnaseq_expr_cor[[i]][first_two_levels_coex_index,])
    
  }
  
  ## get the gene names in each net
  noSNV_icgc_permut_driver_marker_net_genes <- lapply(noSNV_icgc_permut_driver_marker_net, function(i) {
    data.frame(gene_symbol = i[,c(1,2)] %>% 
                 unlist() %>% 
                 unname() %>% 
                 unique())
  })
  
  noSNV_icgc_permut_driver_marker_net_genes <- lapply(noSNV_icgc_permut_driver_marker_net_genes, function(i) {
    
    string_db$map(my_data_frame = i, 
                  my_data_frame_id_col_names = "gene_symbol", 
                  takeFirst = TRUE, removeUnmappedRows = TRUE)
    
  })
  
  # Get the PPIs
  noSNV_icgc_permut_driver_marker_net_ppi <- lapply(noSNV_icgc_permut_driver_marker_net_genes, function(i) {
    
    string_db$get_interactions(string_ids = i$STRING_id)
    
  })
  
  # Convert stringdb ids to gene symbols
  noSNV_icgc_permut_driver_marker_net_ppi <- lapply(1:length(noSNV_icgc_permut_driver_marker_net_ppi), function(i) {
    tmp_row_index <- match(noSNV_icgc_permut_driver_marker_net_ppi[[i]]$from, noSNV_icgc_permut_driver_marker_net_genes[[i]]$STRING_id)
    noSNV_icgc_permut_driver_marker_net_ppi[[i]]$from <- noSNV_icgc_permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_row_index]
    
    tmp_column_index <- match(noSNV_icgc_permut_driver_marker_net_ppi[[i]]$to, noSNV_icgc_permut_driver_marker_net_genes[[i]]$STRING_id)
    noSNV_icgc_permut_driver_marker_net_ppi[[i]]$to <- noSNV_icgc_permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_column_index]
    
    tmp_tbl <- noSNV_icgc_permut_driver_marker_net_ppi[[i]][,c(1,2)]
    tmp_tbl$type <- "PPI"
    
    ## Remove duplicates
    tmp_tbl <- tmp_tbl[-which(duplicated(paste(tmp_tbl$from,
                                               tmp_tbl$to, 
                                               sep = "_"))),]
    
    return(tmp_tbl)
  })
  
  names(noSNV_icgc_permut_driver_marker_net_ppi) <- names(noSNV_icgc_permut_driver_marker_net_genes)
  
  # Combine PPI with the original driver/marker net
  noSNV_icgc_permut_driver_marker_net <- lapply(1:length(noSNV_icgc_permut_driver_marker_net), function(i) {
    rbind(noSNV_icgc_permut_driver_marker_net[[i]], noSNV_icgc_permut_driver_marker_net_ppi[[i]])
  })
  
  names(noSNV_icgc_permut_driver_marker_net) <- names(noSNV_icgc_permut_driver_marker_net_ppi)
  
  # Reconstruct the networks
  noSNV_icgc_permut_driver_marker_net <- lapply(noSNV_icgc_permut_driver_marker_net, function(i) {
    igraph::graph_from_data_frame(i, directed=FALSE)
  })
  
  ## Check the network layers
  table(igraph::E(noSNV_icgc_permut_driver_marker_net$`Breast-AdenoCA`)$type)
  
  ###############################
  
  # Calculate the IVI values of genes within each net
  library(influential)
  noSNV_icgc_permut_driver_marker_net_ivi <- lapply(noSNV_icgc_permut_driver_marker_net, function(i) ivi(i, verbose = T))
  
  ###############################
  
  # Calculate the primitive gene driver scores
  noSNV_icgc_permut_driver_marker_net_gene_tbl <- lapply(noSNV_icgc_permut_driver_marker_net_ivi, function(i) {
    data.frame(gene = names(i), ivi = i)
  })
  
  ## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer
  
  noSNV_icgc_permut_snv_mean_deletriousness <- lapply(icgc_snv_data_flt_withCADD_summarized, function(i) {
    colMeans(i)
  })
  
  ## Generate the gene table for each net and calculate the primitive scores
  noSNV_icgc_permut_driver_marker_net_gene_tbl <- lapply(1:length(noSNV_icgc_permut_driver_marker_net_gene_tbl), function(i) {
    tmp_tbl <- noSNV_icgc_permut_driver_marker_net_gene_tbl[[i]]
    tmp_tbl$mean_deletriousness <- 0
    match_index <- which(names(noSNV_icgc_permut_snv_mean_deletriousness[[i]]) %in% rownames(tmp_tbl))
    tmp_tbl[names(noSNV_icgc_permut_snv_mean_deletriousness[[i]])[match_index], "mean_deletriousness"] <- 
      noSNV_icgc_permut_snv_mean_deletriousness[[i]][match_index]
    
    tmp_tbl$primitive_driver_score <- (tmp_tbl$ivi)*(tmp_tbl$mean_deletriousness) # node weights
    
    return(tmp_tbl)
    
  })
  
  names(noSNV_icgc_permut_driver_marker_net_gene_tbl) <- names(noSNV_icgc_permut_snv_mean_deletriousness)
  
  # add gene scores
  noSNV_icgc_permut_driver_marker_net <- lapply(1:length(noSNV_icgc_permut_driver_marker_net), function(i) {
    set_vertex_attr(graph = noSNV_icgc_permut_driver_marker_net[[i]], name = "score", value = noSNV_icgc_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score)
  })
  
  names(noSNV_icgc_permut_driver_marker_net) <- names(noSNV_icgc_permut_snv_mean_deletriousness)
  
  # Convert the un-weighted to a node-weighted network
  noSNV_icgc_permut_driver_marker_net <- lapply(1:length(noSNV_icgc_permut_driver_marker_net), function(i) {
    set_vertex_attr(graph = noSNV_icgc_permut_driver_marker_net[[i]], name = "weight", value = (noSNV_icgc_permut_driver_marker_net_gene_tbl[[i]]$ivi)^(noSNV_icgc_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness))
  })
  
  names(noSNV_icgc_permut_driver_marker_net) <- names(noSNV_icgc_permut_snv_mean_deletriousness)
  
  ## Check the network node weights
  summary(igraph::V(noSNV_icgc_permut_driver_marker_net$`Breast-AdenoCA`)$weight)
  
  # Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
  noSNV_icgc_permut_driver_marker_net_gene_tbl <- lapply(1:length(noSNV_icgc_permut_driver_marker_net), function(i) {
    tmp_tbl <- noSNV_icgc_permut_driver_marker_net_gene_tbl[[i]]
    tmp_tbl$mean_neighborhood_weight <- sapply(1:nrow(tmp_tbl), function(j) {
      igraph::neighbors(graph = noSNV_icgc_permut_driver_marker_net[[i]], 
                        v = V(noSNV_icgc_permut_driver_marker_net[[i]])[j])$score %>% mean()
    })
    
    tmp_tbl$mediator_score <- tmp_tbl$primitive_driver_score * tmp_tbl$mean_neighborhood_weight
    
    return(tmp_tbl)
  })
  
  names(noSNV_icgc_permut_driver_marker_net_gene_tbl) <- names(noSNV_icgc_permut_snv_mean_deletriousness)
  
  # Detect communities/modules/clusters
  noSNV_icgc_permut_driver_marker_net_modules <- lapply(noSNV_icgc_permut_driver_marker_net, function(i) {
    set.seed(3847)
    igraph::cluster_leiden(
      graph = i,
      objective_function = "CPM",
      weights = NULL,
      resolution_parameter = 0.5,
      beta = 0.05,
      initial_membership = NULL,
      n_iterations = 100,
      vertex_weights = igraph::V(i)$weight
    )
  })
  
  sizes(noSNV_icgc_permut_driver_marker_net_modules$`Breast-AdenoCA`) %>% as.vector() %>% summary()
  length(noSNV_icgc_permut_driver_marker_net_modules$`Breast-AdenoCA`)
  
  ## Inspect and filter modules
  noSNV_icgc_permut_driver_marker_net_modules_flt <- lapply(1:length(noSNV_icgc_permut_driver_marker_net_modules), function(i) {
    
    ### Create a table of modules and their genes
    tmp_tbl <- data.frame(
      module = membership(noSNV_icgc_permut_driver_marker_net_modules[[i]]) %>% as.integer(),
      gene = noSNV_icgc_permut_driver_marker_net_modules[[i]]$names
    )
    
    ### Add if the gene is mutated or not (could be a risk gene or not)
    tmp_tbl$mutated <- 0
    tmp_tbl$mutated[which(noSNV_icgc_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness > 0)] <- 1
    
    ### Separate modules
    tmp_tbl <- lapply(unique(tmp_tbl$module), function(j){
      subset(tmp_tbl, module == j)
    })
    
    names(tmp_tbl) <- sapply(tmp_tbl, function(m) {
      unique(m$module)
    })
    
    ### Remove modules that have less than 2 Risk genes or their size is less than 4
    tmp_tbl <- lapply(tmp_tbl, function(k) {
      if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
        k
      } else {
        NULL
      }
    })
    
    ### Remove NULL modules
    tmp_tbl[which(sapply(tmp_tbl, is.null))] <- NULL
    
    ### Calculate the module scores
    tmp_tbl <- lapply(tmp_tbl, function(l) {
      tmp_tbl4score <- cbind(l, mean_score = 0)
      tmp_tbl4score$mean_score <- mean(noSNV_icgc_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score[which(noSNV_icgc_permut_driver_marker_net_gene_tbl[[i]]$gene %in% l$gene)])
      tmp_tbl4score
    })
    
    tmp_tbl
  })
  
  names(noSNV_icgc_permut_driver_marker_net_modules_flt) <- names(noSNV_icgc_permut_driver_marker_net_modules)
  
  sapply(noSNV_icgc_permut_driver_marker_net_modules_flt, length)

  ## Sort the mean scores of all modules of each cancer type
  noSNV_icgc_permut_driver_marker_net_modules_flt_scores <-
    lapply(noSNV_icgc_permut_driver_marker_net_modules_flt, function(i) {
      tmp_tbl <-
        sapply(i, function(j) {unique(j$mean_score)})
      tmp_tbl <- rev(sort(tmp_tbl))
      tmp_tbl
    })
  
  ### View the number of first ranked modules
  sapply(noSNV_icgc_permut_driver_marker_net_modules_flt_scores, function(i) {i[1]}) %>% as.data.frame() %>% View()
  
  noSNV_icgc_permut_driver_marker_net_modules_flt$`Breast-AdenoCA`$`767` %>% View()
  
  ########
  
  # Top 20 gene 
  noSNV_icgc_permut_driver_marker_net_gene_tbl_top20 <-
    lapply(noSNV_icgc_permut_driver_marker_net_gene_tbl, function(i) {
      slice_max(i, n = 20, order_by = primitive_driver_score) %>% select(gene)
    })
  
  noSNV_icgc_permut_driver_marker_net_gene_tbl_top20 <- do.call(cbind, noSNV_icgc_permut_driver_marker_net_gene_tbl_top20)
  colnames(noSNV_icgc_permut_driver_marker_net_gene_tbl_top20) <- names(noSNV_icgc_permut_driver_marker_net_gene_tbl)
  
  ###########################################
  ###########################################
  ###########################################
  
  ## Save the gene data and modules as supplementary files
  
  tmp_noSNV_icgc_permut_driver_marker_net_gene_tbl <- noSNV_icgc_permut_driver_marker_net_gene_tbl

  for(i in 1:5) {
    colnames(tmp_noSNV_icgc_permut_driver_marker_net_gene_tbl[[i]]) <- c("Gene", "IVI", "Mean_deletriousness_score", "Driver_score", "Mean_neighborhood_weight", "Driverness_mediator_score")
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_noSNV_icgc_permut_driver_marker_net_gene_tbl, file = "Results/Evaluations/icgc_permutation/noSNV_icgc_permut_driver_marker_net_gene_tbl.xlsx")
  
  rm(tmp_noSNV_icgc_permut_driver_marker_net_gene_tbl)
  
  ############################
  
  ### save the modules
  
  tmp_noSNV_icgc_permut_driver_marker_net_modules_flt <- noSNV_icgc_permut_driver_marker_net_modules_flt

  tmp_noSNV_icgc_permut_driver_marker_net_modules_flt_scores <- noSNV_icgc_permut_driver_marker_net_modules_flt_scores
  
  tmp_noSNV_icgc_permut_driver_marker_net_modules_flt <- lapply(tmp_noSNV_icgc_permut_driver_marker_net_modules_flt, function(i) {
    
    tmp_data <- i
    tmp_data <- tmp_data[rev(order(sapply(tmp_data, function(j) unique(j$mean_score))))]
    
    tmp_data
    
  })
  
  for(i in 1:5) {
    tmp_noSNV_icgc_permut_driver_marker_net_modules_flt[[i]] <-
      lapply(1:length(tmp_noSNV_icgc_permut_driver_marker_net_modules_flt[[i]]), function(j) {
        tmp_data <- tmp_noSNV_icgc_permut_driver_marker_net_modules_flt[[i]][[j]]
        colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
        tmp_data[tmp_data == 1] <- "Yes"
        tmp_data[tmp_data == 0] <- "No"
        tmp_data$Module <- paste0("Module number ", j)
        tmp_data <- data.frame(Module = unique(tmp_data$Module),
                               Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                               Mean_score = unique(tmp_data$Mean_score))
        tmp_data
      })
    tmp_noSNV_icgc_permut_driver_marker_net_modules_flt[[i]] <- do.call(rbind, tmp_noSNV_icgc_permut_driver_marker_net_modules_flt[[i]])
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_noSNV_icgc_permut_driver_marker_net_modules_flt, file = "Results/Evaluations/icgc_permutation/noSNV_icgc_permut_driver_marker_net_modules_flt.xlsx")
  
  rm(tmp_noSNV_icgc_permut_driver_marker_net_modules_flt, tmp_noSNV_icgc_permut_driver_marker_net_modules_flt_scores)
  
  ###############################################################################################################
  ###############################################################################################################
  ###############################################################################################################
  
  # Test the pipeline based on icgc_permutation (shuffling only deleteriousness data) 
  
  ## Prepare cor data
  
  ### Setting snv cor to the ones derived from icgc_permutated data
  onlySNV_icgc_permut_snv_flt_cor <- icgc_permut_snv_flt_cor
  
  ### Setting expr cor to the ones derived from NON-icgc_permutated data
  onlySNV_icgc_permut_rnaseq_expr_cor <- icgc_tx_cor_list

  # Merge the similarity tables
  onlySNV_icgc_permut_driver_marker_net <- onlySNV_icgc_permut_snv_flt_cor
  
  # Add the co-expression data to the tables
  
  for(i in 1:length(onlySNV_icgc_permut_driver_marker_net)) {
    
    genes <- onlySNV_icgc_permut_driver_marker_net[[i]][,c(1,2)] %>% 
      unlist() %>% 
      unname() %>% 
      unique()
    
    first_level_coex_index <- c(which(onlySNV_icgc_permut_rnaseq_expr_cor[[i]]$from %in% genes), 
                                which(onlySNV_icgc_permut_rnaseq_expr_cor[[i]]$to %in% genes)) %>% unique()
    
    first_level_coex_genes <- c(genes, 
                                onlySNV_icgc_permut_rnaseq_expr_cor[[i]]$from[first_level_coex_index],
                                onlySNV_icgc_permut_rnaseq_expr_cor[[i]]$to[first_level_coex_index]) %>% unique()
    
    first_two_levels_coex_index <- c(which(onlySNV_icgc_permut_rnaseq_expr_cor[[i]]$from %in% first_level_coex_genes), 
                                     which(onlySNV_icgc_permut_rnaseq_expr_cor[[i]]$to %in% first_level_coex_genes)) %>% unique()
    
    onlySNV_icgc_permut_driver_marker_net[[i]] <- rbind(onlySNV_icgc_permut_driver_marker_net[[i]], onlySNV_icgc_permut_rnaseq_expr_cor[[i]][first_two_levels_coex_index,])
    
  }
  
  # Add PPI data to the driver/marker net
  
  # Map gene names to stringDB ids for PPI analysis (STRING  11.5)
  
  ## get the gene names in each risk net
  onlySNV_icgc_permut_driver_marker_net_genes <- lapply(onlySNV_icgc_permut_driver_marker_net, function(i) {
    data.frame(gene_symbol = i[,c(1,2)] %>% 
                 unlist() %>% 
                 unname() %>% 
                 unique())
  })
  
  ## getSTRINGdb for human
  string_db <- STRINGdb$new(species=9606, version = "11.5", 
                            score_threshold=900, # very high confidence (>0.9)
                            input_directory=""
  )
  
  onlySNV_icgc_permut_driver_marker_net_genes <- lapply(onlySNV_icgc_permut_driver_marker_net_genes, function(i) {
    
    string_db$map(my_data_frame = i, 
                  my_data_frame_id_col_names = "gene_symbol", 
                  takeFirst = TRUE, removeUnmappedRows = TRUE)
    
  })
  
  # Get the PPIs
  onlySNV_icgc_permut_driver_marker_net_ppi <- lapply(onlySNV_icgc_permut_driver_marker_net_genes, function(i) {
    
    string_db$get_interactions(string_ids = i$STRING_id)
    
  })
  
  # Convert stringdb ids to gene symbols
  onlySNV_icgc_permut_driver_marker_net_ppi <- lapply(1:length(onlySNV_icgc_permut_driver_marker_net_ppi), function(i) {
    tmp_row_index <- match(onlySNV_icgc_permut_driver_marker_net_ppi[[i]]$from, onlySNV_icgc_permut_driver_marker_net_genes[[i]]$STRING_id)
    onlySNV_icgc_permut_driver_marker_net_ppi[[i]]$from <- onlySNV_icgc_permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_row_index]
    
    tmp_column_index <- match(onlySNV_icgc_permut_driver_marker_net_ppi[[i]]$to, onlySNV_icgc_permut_driver_marker_net_genes[[i]]$STRING_id)
    onlySNV_icgc_permut_driver_marker_net_ppi[[i]]$to <- onlySNV_icgc_permut_driver_marker_net_genes[[i]]$gene_symbol[tmp_column_index]
    
    tmp_tbl <- onlySNV_icgc_permut_driver_marker_net_ppi[[i]][,c(1,2)]
    tmp_tbl$type <- "PPI"
    
    ## Remove duplicates
    tmp_tbl <- tmp_tbl[-which(duplicated(paste(tmp_tbl$from,
                                               tmp_tbl$to, 
                                               sep = "_"))),]
    
    return(tmp_tbl)
  })
  
  names(onlySNV_icgc_permut_driver_marker_net_ppi) <- names(onlySNV_icgc_permut_driver_marker_net_genes)
  
  # Combine PPI with the original driver/marker net
  onlySNV_icgc_permut_driver_marker_net <- lapply(1:length(onlySNV_icgc_permut_driver_marker_net), function(i) {
    rbind(onlySNV_icgc_permut_driver_marker_net[[i]], onlySNV_icgc_permut_driver_marker_net_ppi[[i]])
  })
  
  names(onlySNV_icgc_permut_driver_marker_net) <- names(onlySNV_icgc_permut_driver_marker_net_genes)
  
  # Reconstruct the networks
  onlySNV_icgc_permut_driver_marker_net <- lapply(onlySNV_icgc_permut_driver_marker_net, function(i) {
    igraph::graph_from_data_frame(i, directed=FALSE)
  })
  
  ## Check the network layers
  table(igraph::E(onlySNV_icgc_permut_driver_marker_net$`Breast-AdenoCA`)$type)
  
  ###############################
  
  # Calculate the IVI values of genes within each net
  library(influential)
  onlySNV_icgc_permut_driver_marker_net_ivi <- lapply(onlySNV_icgc_permut_driver_marker_net, function(i) ivi(i, verbose = T))
  
  ###############################
  
  # Calculate the primitive gene driver scores
  onlySNV_icgc_permut_driver_marker_net_gene_tbl <- lapply(onlySNV_icgc_permut_driver_marker_net_ivi, function(i) {
    data.frame(gene = names(i), ivi = i)
  })
  
  ## Calculate the mean deleteriousness (mean of probabilities (calibrated cadd_phred)) of each gene in each cancer
  
  onlySNV_icgc_permut_snv_mean_deletriousness <- lapply(icgc_permut_snv_flt_withCADD_summarized, function(i) {
    colMeans(i)
  })
  
  ## Generate the gene table for each net and calculate the primitive scores
  onlySNV_icgc_permut_driver_marker_net_gene_tbl <- lapply(1:length(onlySNV_icgc_permut_driver_marker_net_gene_tbl), function(i) {
    tmp_tbl <- onlySNV_icgc_permut_driver_marker_net_gene_tbl[[i]]
    tmp_tbl$mean_deletriousness <- 0
    match_index <- which(names(onlySNV_icgc_permut_snv_mean_deletriousness[[i]]) %in% rownames(tmp_tbl))
    tmp_tbl[names(onlySNV_icgc_permut_snv_mean_deletriousness[[i]])[match_index], "mean_deletriousness"] <- 
      onlySNV_icgc_permut_snv_mean_deletriousness[[i]][match_index]
    
    tmp_tbl$primitive_driver_score <- (tmp_tbl$ivi)*(tmp_tbl$mean_deletriousness) # node weights
    
    return(tmp_tbl)
    
  })
  
  names(onlySNV_icgc_permut_driver_marker_net_gene_tbl) <- names(onlySNV_icgc_permut_snv_mean_deletriousness)
  
  # add gene scores
  onlySNV_icgc_permut_driver_marker_net <- lapply(1:length(onlySNV_icgc_permut_driver_marker_net), function(i) {
    set_vertex_attr(graph = onlySNV_icgc_permut_driver_marker_net[[i]], name = "score", value = onlySNV_icgc_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score)
  })
  
  names(onlySNV_icgc_permut_driver_marker_net) <- names(onlySNV_icgc_permut_snv_mean_deletriousness)
  
  # Convert the un-weighted to a node-weighted network
  onlySNV_icgc_permut_driver_marker_net <- lapply(1:length(onlySNV_icgc_permut_driver_marker_net), function(i) {
    set_vertex_attr(graph = onlySNV_icgc_permut_driver_marker_net[[i]], name = "weight", value = (onlySNV_icgc_permut_driver_marker_net_gene_tbl[[i]]$ivi)^(onlySNV_icgc_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness))
  })
  
  names(onlySNV_icgc_permut_driver_marker_net) <- names(onlySNV_icgc_permut_snv_mean_deletriousness)
  
  ## Check the network node weights
  summary(igraph::V(onlySNV_icgc_permut_driver_marker_net$`Breast-AdenoCA`)$weight)
  
  # Calculate the mediator scores by multiplying the primitive scores to the mean neighborhood weight of nodes
  onlySNV_icgc_permut_driver_marker_net_gene_tbl <- lapply(1:length(onlySNV_icgc_permut_driver_marker_net), function(i) {
    tmp_tbl <- onlySNV_icgc_permut_driver_marker_net_gene_tbl[[i]]
    tmp_tbl$mean_neighborhood_weight <- sapply(1:nrow(tmp_tbl), function(j) {
      igraph::neighbors(graph = onlySNV_icgc_permut_driver_marker_net[[i]], 
                        v = V(onlySNV_icgc_permut_driver_marker_net[[i]])[j])$score %>% mean()
    })
    
    tmp_tbl$mediator_score <- tmp_tbl$primitive_driver_score * tmp_tbl$mean_neighborhood_weight
    
    return(tmp_tbl)
  })
  
  names(onlySNV_icgc_permut_driver_marker_net_gene_tbl) <- names(onlySNV_icgc_permut_snv_mean_deletriousness)
  
  # Detect communities/modules/clusters
  onlySNV_icgc_permut_driver_marker_net_modules <- lapply(onlySNV_icgc_permut_driver_marker_net, function(i) {
    set.seed(3847)
    igraph::cluster_leiden(
      graph = i,
      objective_function = "CPM",
      weights = NULL,
      resolution_parameter = 0.5,
      beta = 0.05,
      initial_membership = NULL,
      n_iterations = 100,
      vertex_weights = igraph::V(i)$weight
    )
  })
  
  sizes(onlySNV_icgc_permut_driver_marker_net_modules$`Breast-AdenoCA`) %>% as.vector() %>% summary()
  length(onlySNV_icgc_permut_driver_marker_net_modules$`Breast-AdenoCA`)
  
  ## Inspect and filter modules
  onlySNV_icgc_permut_driver_marker_net_modules_flt <- lapply(1:length(onlySNV_icgc_permut_driver_marker_net_modules), function(i) {
    
    ### Create a table of modules and their genes
    tmp_tbl <- data.frame(
      module = membership(onlySNV_icgc_permut_driver_marker_net_modules[[i]]) %>% as.integer(),
      gene = onlySNV_icgc_permut_driver_marker_net_modules[[i]]$names
    )
    
    ### Add if the gene is mutated or not (could be a risk gene or not)
    tmp_tbl$mutated <- 0
    tmp_tbl$mutated[which(onlySNV_icgc_permut_driver_marker_net_gene_tbl[[i]]$mean_deletriousness > 0)] <- 1
    
    ### Separate modules
    tmp_tbl <- lapply(unique(tmp_tbl$module), function(j){
      subset(tmp_tbl, module == j)
    })
    
    names(tmp_tbl) <- sapply(tmp_tbl, function(m) {
      unique(m$module)
    })
    
    ### Remove modules that have less than 2 Risk genes or their size is less than 4
    tmp_tbl <- lapply(tmp_tbl, function(k) {
      if(sum(k$mutated == 1) >=2 & nrow(k) >= 4) {
        k
      } else {
        NULL
      }
    })
    
    ### Remove NULL modules
    tmp_tbl[which(sapply(tmp_tbl, is.null))] <- NULL
    
    ### Calculate the module scores
    tmp_tbl <- lapply(tmp_tbl, function(l) {
      tmp_tbl4score <- cbind(l, mean_score = 0)
      tmp_tbl4score$mean_score <- mean(onlySNV_icgc_permut_driver_marker_net_gene_tbl[[i]]$primitive_driver_score[which(onlySNV_icgc_permut_driver_marker_net_gene_tbl[[i]]$gene %in% l$gene)])
      tmp_tbl4score
    })
    
    tmp_tbl
  })
  
  names(onlySNV_icgc_permut_driver_marker_net_modules_flt) <- names(onlySNV_icgc_permut_driver_marker_net_modules)
  
  sapply(onlySNV_icgc_permut_driver_marker_net_modules_flt, length)

  ## Sort the mean scores of all modules of each cancer type
  onlySNV_icgc_permut_driver_marker_net_modules_flt_scores <-
    lapply(onlySNV_icgc_permut_driver_marker_net_modules_flt, function(i) {
      tmp_tbl <-
        sapply(i, function(j) {unique(j$mean_score)})
      tmp_tbl <- rev(sort(tmp_tbl))
      tmp_tbl
    })
  
  ### View the number of first ranked modules
  sapply(onlySNV_icgc_permut_driver_marker_net_modules_flt_scores, function(i) {i[1]}) %>% as.data.frame() %>% View()
  
  onlySNV_icgc_permut_driver_marker_net_modules_flt$`Breast-AdenoCA`$`1488` %>% View()
  
  ########
  
  # Top 20 gene 
  onlySNV_icgc_permut_driver_marker_net_gene_tbl_top20 <-
    lapply(onlySNV_icgc_permut_driver_marker_net_gene_tbl, function(i) {
      slice_max(i, n = 20, order_by = primitive_driver_score) %>% select(gene)
    })
  
  onlySNV_icgc_permut_driver_marker_net_gene_tbl_top20 <- do.call(cbind, onlySNV_icgc_permut_driver_marker_net_gene_tbl_top20)
  colnames(onlySNV_icgc_permut_driver_marker_net_gene_tbl_top20) <- names(onlySNV_icgc_permut_driver_marker_net_gene_tbl)
  
  ###########################################
  ###########################################
  ###########################################
  
  ## Save the gene data and modules as supplementary files
  
  tmp_onlySNV_icgc_permut_driver_marker_net_gene_tbl <- onlySNV_icgc_permut_driver_marker_net_gene_tbl

  for(i in 1:5) {
    colnames(tmp_onlySNV_icgc_permut_driver_marker_net_gene_tbl[[i]]) <- c("Gene", "IVI", "Mean_deletriousness_score", "Driver_score", "Mean_neighborhood_weight", "Driverness_mediator_score")
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_onlySNV_icgc_permut_driver_marker_net_gene_tbl, file = "Results/Evaluations/icgc_permutation/onlySNV_icgc_permut_driver_marker_net_gene_tbl.xlsx")
  
  rm(tmp_onlySNV_icgc_permut_driver_marker_net_gene_tbl)
  
  ############################
  
  ### save the modules
  
  tmp_onlySNV_icgc_permut_driver_marker_net_modules_flt <- onlySNV_icgc_permut_driver_marker_net_modules_flt

  tmp_onlySNV_icgc_permut_driver_marker_net_modules_flt_scores <- onlySNV_icgc_permut_driver_marker_net_modules_flt_scores
  
  tmp_onlySNV_icgc_permut_driver_marker_net_modules_flt <- lapply(tmp_onlySNV_icgc_permut_driver_marker_net_modules_flt, function(i) {
    
    tmp_data <- i
    tmp_data <- tmp_data[rev(order(sapply(tmp_data, function(j) unique(j$mean_score))))]
    
    tmp_data
    
  })
  
  for(i in 1:5) {
    tmp_onlySNV_icgc_permut_driver_marker_net_modules_flt[[i]] <-
      lapply(1:length(tmp_onlySNV_icgc_permut_driver_marker_net_modules_flt[[i]]), function(j) {
        tmp_data <- tmp_onlySNV_icgc_permut_driver_marker_net_modules_flt[[i]][[j]]
        colnames(tmp_data) <- c("Module", "Genes", "Mutated", "Mean_score")
        tmp_data[tmp_data == 1] <- "Yes"
        tmp_data[tmp_data == 0] <- "No"
        tmp_data$Module <- paste0("Module number ", j)
        tmp_data <- data.frame(Module = unique(tmp_data$Module),
                               Genes =  paste0(paste(tmp_data$Genes, " (Mutated: ", tmp_data$Mutated, ")", sep = ""), collapse = ", "),
                               Mean_score = unique(tmp_data$Mean_score))
        tmp_data
      })
    tmp_onlySNV_icgc_permut_driver_marker_net_modules_flt[[i]] <- do.call(rbind, tmp_onlySNV_icgc_permut_driver_marker_net_modules_flt[[i]])
  }
  
  library(openxlsx)
  write.xlsx(x = tmp_onlySNV_icgc_permut_driver_marker_net_modules_flt, file = "Results/Evaluations/icgc_permutation/onlySNV_icgc_permut_driver_marker_net_modules_flt.xlsx")
  
  rm(tmp_onlySNV_icgc_permut_driver_marker_net_modules_flt, tmp_onlySNV_icgc_permut_driver_marker_net_modules_flt_scores)
  
  #############################################################################################
  
  # Visualize icgc_permutation results
  
  ########################################
  # define a function for the evaluation of intersections
  
  intersect_combn <- function(named_list, m = 2, vis = TRUE, 
                              label_font_size = 10, cell_font_size = 10,
                              show_column_labels = FALSE, show_row_labels = FALSE, 
                              annot_colors = rainbow(length(named_list)),
                              column_title = NULL, row_title = NULL, title_font_size = 10) { 
    
    library(magrittr)
    library(tibble)
    library(tidyr)
    library(stats)
    library(ComplexHeatmap)
    library(stringr)
    
    # named_list should be a named list, for example:
    ## list("List 1" = c(1:5), "List 2" = c(3:6))
    
    # m is the number of elements to choose
    
    comb_list = combn(seq(length(named_list)), m)
    
    if(vis & m == 2) {
      self_comb = lapply(1:length(named_list), rep, 2) %>% data.frame()
      comb_list <- cbind(comb_list, self_comb)
    }
    
    intersects =
      lapply(1:ncol(comb_list), function(i) {
        
        tmp_intersect =
          base::Reduce(intersect, named_list[comb_list[,i]])
        
        tmp_union =
          base::Reduce(union, named_list[comb_list[,i]])
        
        tmp_jaccard_sim_index = 
          length(tmp_intersect)/length(tmp_union)
        
        tmp_jaccard_dissim_index = 1 - tmp_jaccard_sim_index
        
        tmp_names_list = lapply(1:m, function(j) {
          names(named_list[comb_list[j,i]])
        })
    
        tmp_names_list = setNames(tmp_names_list, paste("Vector", c(1:m), sep = ""))

        intersect_tbl = cbind(data.frame(tmp_names_list),
                              data.frame(Intersection_length = length(tmp_intersect),
                                         Intersects = paste0(tmp_intersect, collapse = ", "),
                                         jaccard_sim_index = tmp_jaccard_sim_index,
                                         jaccard_dissim_index = tmp_jaccard_dissim_index))
        
        intersect_tbl
      })
    
    # merging tables of different combinations
    intersects <- do.call(rbind, intersects)
    
    # generating similarity matrices
    if(vis & m == 2) {
      
      ht_opt$TITLE_PADDING = unit(c(8.5, 8.5), "points")
      
      ## preparing the sim matrix for vis
      jaccard_sim_matrix <- intersects[,c(1,2,5)]
      jaccard_sim_matrix <- jaccard_sim_matrix %>% 
        pivot_wider(names_from = Vector2, values_from = jaccard_sim_index, values_fill = 0) %>% 
        column_to_rownames(var = "Vector1") %>% as.matrix()
      jaccard_sim_matrix <- jaccard_sim_matrix[,rownames(jaccard_sim_matrix)]
      
      ## Preparing the annotations
      har <- HeatmapAnnotation(`Color Legend` = stringr::str_wrap(colnames(jaccard_sim_matrix), width = 20),
                               which = "row",show_annotation_name = FALSE,simple_anno_size = unit(0.3, "cm"), 
                               annotation_legend_param = list(labels_gp = gpar(fontsize = label_font_size), title_gp = gpar(fontsize = title_font_size, fontface = "bold")),
                               col = list(`Color Legend` = setNames(annot_colors, stringr::str_wrap(colnames(jaccard_sim_matrix), width = 20))))
      ha <- HeatmapAnnotation(`Color Legend` = stringr::str_wrap(colnames(jaccard_sim_matrix), width = 20),
                              annotation_legend_param = list(labels_gp = gpar(fontsize = label_font_size), title_gp = gpar(fontsize = title_font_size, fontface = "bold")),
                              which = "column",show_legend = FALSE,show_annotation_name = FALSE, simple_anno_size = unit(0.3, "cm"),
                              col = list(`Color Legend` = setNames(annot_colors, stringr::str_wrap(colnames(jaccard_sim_matrix), width = 20))))
      
      jaccard_sim_plot <- ComplexHeatmap::Heatmap(jaccard_sim_matrix, rect_gp = gpar(type = "none"),
                                                  column_names_side = "top",
                                                  column_names_rot = 45,
                                                  column_names_gp = gpar(fontsize = label_font_size),    # Set column names font size
                                                  column_labels = stringr::str_wrap(colnames(jaccard_sim_matrix), width = 20),
                                                  row_labels = stringr::str_wrap(rownames(jaccard_sim_matrix), width = 20),
                                                  show_column_names = show_column_labels,
                                                  show_row_names = show_row_labels,
                                                  row_names_gp = gpar(fontsize = label_font_size),       # Set row names font size
                                                  column_title = column_title, 
                                                  column_title_gp = gpar(fill = "lightgrey", col = "black", border = "darkgrey", fontsize = title_font_size, fontface = "bold"),
                                                  row_title = row_title, 
                                                  row_title_gp = gpar(fill = "lightgrey", col = "black", border = "darkgrey", fontsize = title_font_size,  fontface = "bold"),
                                                  col = colorRamp2(c(0, 0.5, 1), c("lightblue", "#EEEEEE", "firebrick1")),
                                                  heatmap_legend_param = list(labels_gp = gpar(fontsize = label_font_size), title_gp = gpar(fontsize = title_font_size, fontface = "bold")),
                                                  cluster_rows = FALSE, cluster_columns = FALSE, 
                                                  name = "Jaccard\nSimilarity",
                                                  top_annotation=ha,
                                                  right_annotation=har,
                                                  cell_fun = function(j, i, x, y, w, h, fill) {
                                                    if(j >= i) {
                                                      grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "darkgrey"))
                                                      grid.text(sprintf("%.1f", jaccard_sim_matrix[i, j]), x, y, gp = gpar(fontsize = cell_font_size))
                                                    }
                                                  })
      
      ####################
      
      ## preparing the dissim matrix for vis
      jaccard_dissim_matrix <- intersects[,c(1,2,6)]
      jaccard_dissim_matrix <- jaccard_dissim_matrix %>% 
        pivot_wider(names_from = Vector2, values_from = jaccard_dissim_index, values_fill = 0) %>% 
        column_to_rownames(var = "Vector1") %>% as.matrix()
      jaccard_dissim_matrix <- jaccard_dissim_matrix[,rownames(jaccard_dissim_matrix)]
    
    jaccard_dissim_plot <- ComplexHeatmap::Heatmap(jaccard_dissim_matrix, rect_gp = gpar(type = "none"),
                                                   column_names_side = "top",column_names_rot = 45,
                                                   column_title = column_title, 
                                                   column_title_gp = gpar(fill = "lightgrey", col = "black", border = "darkgrey", fontsize = title_font_size, fontface = "bold"),
                                                   show_column_names = show_column_labels,
                                                   show_row_names = show_row_labels,
                                                   row_title = row_title, 
                                                   row_title_gp = gpar(fill = "lightgrey", col = "black", border = "darkgrey", fontsize = title_font_size, fontface = "bold"),
                                                   column_labels = stringr::str_wrap(colnames(jaccard_dissim_matrix), width = 20),
                                                   row_labels = stringr::str_wrap(rownames(jaccard_dissim_matrix), width = 20),
                                                   col = colorRamp2(c(0, 0.5, 1), c("greenyellow", "#EEEEEE", "mediumpurple1" )),
                                                   cluster_rows = FALSE, cluster_columns = FALSE,
                                                   top_annotation=ha,
                                                   right_annotation=har,
                                                   heatmap_legend_param = list(labels_gp = gpar(fontsize = label_font_size), title_gp = gpar(fontsize = title_font_size, fontface = "bold")),
                                                   name = "Jaccard\nDissimilarity",
                                                   row_names_gp = gpar(fontsize = label_font_size),       # Set row names font size
                                                   column_names_gp = gpar(fontsize = label_font_size),    # Set column names font size
                                                   cell_fun = function(j, i, x, y, w, h, fill) {
                                                     if(j >= i) {
                                                       grid.rect(x, y, w, h, gp = gpar(fill = fill, col = "darkgrey"))
                                                       grid.text(sprintf("%.1f", jaccard_dissim_matrix[i, j]), x, y, gp = gpar(fontsize = cell_font_size))
                                                     }
                                                   })
    
    
    
    # removing self intersections
    intersects <- intersects[-c((nrow(intersects) - (length(named_list) - 1)):nrow(intersects)),]
    
    return(list(intersects_tbl = intersects,
                jaccard_sim_matrix = jaccard_sim_matrix,
                jaccard_sim_plot = jaccard_sim_plot,
                jaccard_dissim_matrix = jaccard_dissim_matrix,
                jaccard_dissim_plot = jaccard_dissim_plot))
    } else {
      return(intersects)
    }
  }
  
  ########################################
  
  ## for the top 20 genes
  
  # onlySNV_icgc_permut_driver_marker_net_gene_tbl_top20
  # onlySNV_icgc_permut_driver_marker_net_modules_flt
  
  icgc_permut_eval_list <- lapply(1:5, function(i) {
    
    original = icgc_driver_net_gene_tbl_top20[,i]
    onlySNV = onlySNV_icgc_permut_driver_marker_net_gene_tbl_top20[,i]
    onlyPPI = onlyPPI_icgc_permut_driver_marker_net_gene_tbl_top20[,i]
    onlyTx = onlyTx_icgc_permut_driver_marker_net_gene_tbl_top20[,i]
    SNV_RNA = icgc_permut_driver_marker_net_gene_tbl_top20[,i]
    SNV_RNA_PPI = plusPPI_icgc_permut_driver_marker_net_gene_tbl_top20[,i]
    RNA_PPI = noSNV_icgc_permut_driver_marker_net_gene_tbl_top20[,i]
    
    intersect_combn(named_list = list("Original data" = original,
                                      "Permuted mutations" = onlySNV,
                                      "Permuted PPI" = onlyPPI,
                                      "Permuted RNA-Seq" = onlyTx,
                                      "Permuted mutations + RNA-Seq" = SNV_RNA,
                                      "Permuted mutations + RNA-Seq + PPIs" = SNV_RNA_PPI,
                                      "Permuted RNA-Seq + PPIs" = RNA_PPI), 
                    annot_colors = c("#8F1D1EFF", "#DD75D3FF", "#02734AFF", "#FFC737FF", "#3E52A1FF", "#000000FF", "#81e74a"),
                    label_font_size = 5, cell_font_size = 5, title_font_size = 6,
                    column_title =  colnames(onlySNV_icgc_permut_driver_marker_net_gene_tbl_top20)[i])
  })
  
  names(icgc_permut_eval_list) <- colnames(onlySNV_icgc_permut_driver_marker_net_gene_tbl_top20)
  
  icgc_permut_eval_tbl <- do.call(rbind, lapply(icgc_permut_eval_list, function(i) i$intersects_tbl))
  icgc_permut_eval_tbl$Cancer_type <- rep(names(icgc_permut_eval_list), each = 21)
  
  colnames(icgc_permut_eval_tbl)[c(1,2)] <- paste("Origin_of_gene_set_", c(1,2), sep = "")
  
  icgc_permut_eval_list$`Breast-AdenoCA`$jaccard_sim_plot
  
  # save plots
  Reduce(`+`, lapply(icgc_permut_eval_list, function(i) i$jaccard_sim_plot))
  
  #########################
  
  ## for the first ranked modules
  
  icgc_permut_modules_eval_list <-
    lapply(1:5, function(i) {
      
      ## prepare the first ranked modules
      original_modules <- icgc_driver_net_modules_flt[[i]]
      original_module_scores <- icgc_driver_net_modules_flt_scores[[i]]
      original_first_five_modules = sapply(1, function(j) {
        tmp_module =
          original_modules[[names(original_module_scores)[j]]]
        
        tmp_module$gene
      }) %>% unlist()
      ###########
      onlySNV_modules <- onlySNV_icgc_permut_driver_marker_net_modules_flt[[i]]
      onlySNV_module_scores <- onlySNV_icgc_permut_driver_marker_net_modules_flt_scores[[i]]
      onlySNV_first_five_modules = sapply(1, function(j) {
        tmp_module =
          onlySNV_modules[[names(onlySNV_module_scores)[j]]]
        
        tmp_module$gene
      }) %>% unlist()
      ###########
      onlyPPI_modules <- onlyPPI_icgc_permut_driver_marker_net_modules_flt[[i]]
      onlyPPI_module_scores <- onlyPPI_icgc_permut_driver_marker_net_modules_flt_scores[[i]]
      onlyPPI_first_five_modules = sapply(1, function(j) {
        tmp_module =
          onlyPPI_modules[[names(onlyPPI_module_scores)[j]]]
        
        tmp_module$gene
      }) %>% unlist()
      ###########
      onlyTx_modules <- onlyTx_icgc_permut_driver_marker_net_modules_flt[[i]]
      onlyTx_module_scores <- onlyTx_icgc_permut_driver_marker_net_modules_flt_scores[[i]]
      onlyTx_first_five_modules = sapply(1, function(j) {
        tmp_module =
          onlyTx_modules[[names(onlyTx_module_scores)[j]]]
        
        tmp_module$gene
      }) %>% unlist()
      ###########
      SNV_RNA_modules <- icgc_permut_driver_marker_net_modules_flt[[i]]
      SNV_RNA_module_scores <- icgc_permut_driver_marker_net_modules_flt_scores[[i]]
      SNV_RNA_first_five_modules = sapply(1, function(j) {
        tmp_module =
          SNV_RNA_modules[[names(SNV_RNA_module_scores)[j]]]
        
        tmp_module$gene
      }) %>% unlist()
      ###########
      SNV_RNA_PPI_modules <- plusPPI_icgc_permut_driver_marker_net_modules_flt[[i]]
      SNV_RNA_PPI_module_scores <- plusPPI_icgc_permut_driver_marker_net_modules_flt_scores[[i]]
      SNV_RNA_PPI_first_five_modules = sapply(1, function(j) {
        tmp_module =
          SNV_RNA_PPI_modules[[names(SNV_RNA_PPI_module_scores)[j]]]
        
        tmp_module$gene
      }) %>% unlist()
      ###########
      RNA_PPI_modules <- noSNV_icgc_permut_driver_marker_net_modules_flt[[i]]
      RNA_PPI_module_scores <- noSNV_icgc_permut_driver_marker_net_modules_flt_scores[[i]]
      RNA_PPI_first_five_modules = sapply(1, function(j) {
        tmp_module =
          RNA_PPI_modules[[names(RNA_PPI_module_scores)[j]]]
        
        tmp_module$gene
      }) %>% unlist()
      ###########
      
      tmp_intersect_tbl =
        intersect_combn(named_list = list("Original data" = original_first_five_modules,
                                          "Permuted mutations" = onlySNV_first_five_modules,
                                          "Permuted PPI" = onlyPPI_first_five_modules,
                                          "Permuted RNA-Seq" = onlyTx_first_five_modules,
                                          "Permuted mutations + RNA-Seq" = SNV_RNA_first_five_modules,
                                          "Permuted mutations + RNA-Seq + PPIs" = SNV_RNA_PPI_first_five_modules,
                                          "Permuted RNA-Seq + PPIs" = RNA_PPI_first_five_modules),
                        label_font_size = 8, row_title = "First ranked module")
      
      tmp_intersect_tbl
    })
  
  names(icgc_permut_modules_eval_list) <- names(icgc_permut_eval_list)
  
  icgc_permut_modules_eval_tbl <- do.call(rbind, lapply(icgc_permut_modules_eval_list, function(i) i$intersects_tbl))
  icgc_permut_modules_eval_tbl$Cancer_type <- rep(names(icgc_permut_modules_eval_list), each = 21)
  
  colnames(icgc_permut_modules_eval_tbl)[c(1,2)] <- paste("Origin_of_gene_set_", c(1,2), sep = "")
  
  # save plots
  Reduce(`+`, lapply(icgc_permut_modules_eval_list, function(i) i$jaccard_sim_plot))
  
  ## save tables of results
  openxlsx::write.xlsx(x = list("Top 20 genes" = icgc_permut_eval_tbl, "First ranked module" = icgc_permut_modules_eval_tbl), 
                       file = "Results/Evaluations/icgc_permutation/icgc_permut_eval_tbl.xlsx")
  

  ## Further evaluate the results
  for(i in 1:6) {
    
    cat(paste0(icgc_permut_eval_list$`Breast-AdenoCA`$intersects_tbl$Vector1[i], ": ",
               mean(sapply(lapply(icgc_permut_eval_list, function(i) i$intersects_tbl), function(j) {
                 j[i,3]
               })), sep = "\n"
    ))
  }
  ##############
  
  for(i in 1:6) {
    
    cat(paste0(icgc_permut_modules_eval_list$`Breast-AdenoCA`$intersects_tbl$Vector1[i], ": ",
               mean(sapply(lapply(icgc_permut_modules_eval_list, function(i) i$intersects_tbl), function(j) {
                 j[i,3]
               })), sep = "\n"
    ))
  }
  
  ###########################
  
  ## Chord diagram visualization
  
  library(circlize)
  
  par(mfrow = c(2,5))
  
  icgc_permut_eval_list_filtered <- lapply(icgc_permut_eval_list, function(i) {
    i %>% filter(Vector2 != "icgc_permuted RNA-Seq+PPIs")
  })
  
  icgc_permut_modules_eval_list_filtered <- lapply(icgc_permut_modules_eval_list, function(i) {
    i %>% filter(Vector2 != "icgc_permuted RNA-Seq+PPIs")
  })
  
  # Visualize top 20s
  for(i in 1:5) {
    
    set.seed(1234)
    chordDiagram(x = icgc_permut_eval_list_filtered,annotationTrackHeight = mm_h(c(1, 1)),
                 annotationTrack = c("grid"),
                 grid.col = setNames(c("#8F1D1EFF", "#DD75D3FF", "#02734AFF", "#FFC737FF", "#3E52A1FF", "#000000FF", "orange"), 
                                     union(icgc_permut_eval_list[[i]][,1], icgc_permut_eval_list[[i]][,2])),
                 link.border = rep("black", 5))
    
    # circos.xaxis(h = 0.5, sector.index = "Original data", 
    #              labels = rev(icgc_permut_eval_list[[i]]$Intersection_length[1:4]),
    #              major.tick = F,
    #              minor.ticks = 0, 
    #              labels.cex = 1)
    
    title(names(icgc_permut_eval_list)[i])
    
    circos.clear()
    
  }
  
  ##################
  
  # Visualize the first ranked module
  ## alaki <- icgc_permut_modules_eval_list
  ### alaki <- alaki[c(1:6, 10, 7:9)] # and use alaki instead of icgc_permut_modules_eval_list as the last one has intersection but not the three networks before that
  
  for (i in c(1:5)) {
    
    set.seed(1234)
    chordDiagram(x = icgc_permut_modules_eval_list_filtered[[i]],annotationTrackHeight = mm_h(c(1, 1)),
                 annotationTrack = c("grid"),
                 grid.col = setNames(c("#8F1D1EFF", "#DD75D3FF", "#02734AFF", "#FFC737FF", "#3E52A1FF", "#000000FF", "orange"), 
                                     union(icgc_permut_modules_eval_list[[i]][,1], icgc_permut_modules_eval_list[[i]][,2])),
                 link.border = rep("black", 5))
    
    # circos.xaxis(h = 0.5,
    #              labels = icgc_permut_modules_eval_list[[i]]$Intersection_length[icgc_permut_modules_eval_list[[i]]$Intersection_length > 0],
    #              major.tick = F,
    #              minor.ticks = 0, 
    #              labels.cex = 1)
    cat(paste("Cancer ", names(icgc_permut_modules_eval_list_filtered)[i], " is done\n"))
    
    circos.clear()
  }
  
  dev.off()
  
  par(mfrow = c(1,1))
  
  ## draw the legend
  library(ComplexHeatmap)
  
  Legend(at = union(icgc_permut_eval_list_filtered[[i]][,1], icgc_permut_eval_list_filtered[[i]][,2]), 
         type = "lines", 
         legend_gp = gpar(col = c("#8F1D1EFF", "#DD75D3FF", "#02734AFF", "#FFC737FF", "#3E52A1FF", "#000000FF"), 
                          lwd = 6), 
         title_position = "topleft", 
         title = "icgc_permutation level") %>% draw()
  
  dev.off()
  
#=============================================================================
#
#    Code chunk 15: Dependency and essentiality analysis
#
#=============================================================================
  
  BiocManager::install("depmap")
  
  
  library("dplyr")
  library("ggplot2")
  library("viridis")
  library("tibble")
  library("gridExtra")
  library("stringr")
  library("depmap")
  library("ExperimentHub")
  
  ###########################################################################################################
  
  ## create ExperimentHub query object
  eh <- ExperimentHub()
  depmap_resource <- query(eh, "depmap")
  depmap_resource <- data.frame(id = depmap_resource$ah_id,
                                title = depmap_resource$title)
  
  ###########################################################################################################
  
  # Loading the data (latest versions)
  ## The rnai dataset contains the combined genetic dependency data for RNAi - induced gene knockdown for select genes and cancer cell lines
  ### rnai and crispr. These datasets contain a score expressing how vital a particular gene is in terms of how lethal the knockout/knockdown of that gene 
  ### is on a target cell line. For example, a highly negative dependency score implies that a cell line is highly dependent on that gene.
  
  # TPM <- eh[["EH7556"]] # TPM_22Q2
  copyNumber <- eh[["EH7555"]] # copyNumber_22Q2
  mutationCalls <- eh[["EH7557"]] # mutationCalls_22Q2
  rnai <- eh[["EH3080"]] # rnai_19Q3
  metadata <- eh[["EH7558"]] # metadata_22Q2
  crispr <- eh[["EH7554"]] # crispr_22Q2
  drug_sensitivity <- eh[["EH7530"]] # drug_sensitivity_21Q2
  
  ## filter the data
  
  ### remove na/multiple dependency data
  rnai <- rnai[-which(is.na(rnai$dependency)),]
  crispr <- crispr[-which(is.na(crispr$dependency)),]
  drug_sensitivity <- drug_sensitivity[-which(is.na(drug_sensitivity$dependency)),]
  
  rnai <- rnai[-stringr::str_detect(rnai$gene, "&"), ]

  ###########################################################################################################
  
  # Add the tissue type to the data frames
  rnai_tissue_types <- sapply(unique(rnai$cell_line), function(i) {
    string <- i
    tmp_tissue_type <- 
      unlist(str_split(string = string, pattern = "_"))
    if(length(tmp_tissue_type) > 1) {
      tmp_tissue_type <- tmp_tissue_type[-1]
    }
    tmp_tissue_type <- tmp_tissue_type %>% 
      paste0(collapse = "_")
    tmp_tissue_type
  })  %>% unique()
  
  if(any(rnai_tissue_types == "") | any(rnai_tissue_types == "NA")) {
    rnai_tissue_types <- rnai_tissue_types[-which(c(rnai_tissue_types == "") | c(rnai_tissue_types == "NA"))]
  }
  
  rnai$tissue_type <- ""
  lapply(rnai_tissue_types, function(i) {
    rnai$tissue_type[str_detect(string = rnai$cell_line, pattern = i)] <<- i
  })
    
  ##############
  
  crispr_tissue_types <- sapply(unique(crispr$cell_line), function(i) {
    string <- i
    tmp_tissue_type <- 
      unlist(str_split(string = string, pattern = "_"))
    if(length(tmp_tissue_type) > 1) {
      tmp_tissue_type <- tmp_tissue_type[-1]
    }
    tmp_tissue_type <- tmp_tissue_type %>% 
      paste0(collapse = "_")
    tmp_tissue_type
  })  %>% unique()
  
  if(any(crispr_tissue_types == "") | any(crispr_tissue_types == "NA")) {
    crispr_tissue_types <- crispr_tissue_types[-which(c(crispr_tissue_types == "") | c(crispr_tissue_types == "NA"))]
  }
  
  crispr$tissue_type <- ""
  lapply(crispr_tissue_types, function(i) {
    crispr$tissue_type[str_detect(string = crispr$cell_line, pattern = i)] <<- i
  })
  
  ##############
  
  drug_sensitivity_tissue_types <- sapply(unique(drug_sensitivity$cell_line), function(i) {
    string <- i
    tmp_tissue_type <- 
      unlist(str_split(string = string, pattern = "_"))
    if(length(tmp_tissue_type) > 1) {
      tmp_tissue_type <- tmp_tissue_type[-1]
    }
    tmp_tissue_type <- tmp_tissue_type %>% 
      paste0(collapse = "_")
    tmp_tissue_type
  })  %>% unique()
  
  if(any(drug_sensitivity_tissue_types == "") | any(drug_sensitivity_tissue_types == "NA")) {
    drug_sensitivity_tissue_types <- drug_sensitivity_tissue_types[-which(c(drug_sensitivity_tissue_types == "") | c(drug_sensitivity_tissue_types == "NA"))]
  }
  
  drug_sensitivity$tissue_type <- ""
  lapply(drug_sensitivity_tissue_types, function(i) {
    drug_sensitivity$tissue_type[str_detect(string = drug_sensitivity$cell_line, pattern = i)] <<- i
  })
  
  ##############
  
  copyNumber_tissue_types <- sapply(unique(copyNumber$cell_line), function(i) {
    string <- i
    tmp_tissue_type <- 
      unlist(str_split(string = string, pattern = "_"))
    if(length(tmp_tissue_type) > 1) {
      tmp_tissue_type <- tmp_tissue_type[-1]
    }
    tmp_tissue_type <- tmp_tissue_type %>% 
      paste0(collapse = "_")
    tmp_tissue_type
  })  %>% unique()
  
  if(any(copyNumber_tissue_types == "") | any(copyNumber_tissue_types == "NA")) {
    copyNumber_tissue_types <- copyNumber_tissue_types[-which(c(copyNumber_tissue_types == "") | c(copyNumber_tissue_types == "NA"))]
  }
  
  copyNumber$tissue_type <- ""
  lapply(copyNumber_tissue_types, function(i) {
    copyNumber$tissue_type[str_detect(string = copyNumber$cell_line, pattern = i)] <<- i
  })
  
  ###########################################################################################################
  
  ## Add the lineage to the datasets
  
  rnai <- metadata %>% 
    select(depmap_id, lineage) %>% 
    right_join(rnai, by = "depmap_id")
  
  rnai$lineage[is.na(rnai$lineage)] <- rnai$tissue_type[is.na(rnai$lineage)]
  rnai$lineage[rnai$lineage == "unknown"] <- rnai$tissue_type[rnai$lineage == "unknown"]
  
  rnai$lineage <- str_replace(string = rnai$lineage, pattern = "^STOMACH$", replacement = "gastric")
  rnai$lineage <- str_replace(string = rnai$lineage, pattern = "^PANCREAS$", replacement = "pancreas")
  rnai$lineage <- str_replace(string = rnai$lineage, pattern = "^BREAST$", replacement = "breast")
  rnai$lineage <- str_replace(string = rnai$lineage, pattern = "^GASTROINTESTINAL_TRACT$", replacement = "gastrointestinal_tract")
  rnai$lineage <- str_replace(string = rnai$lineage, pattern = "^HAEMATOPOIETIC_AND_LYMPHOID_TISSUE$", replacement = "blood")
  

  ##############
  
  crispr <- metadata %>% 
    select(depmap_id, lineage) %>% 
    right_join(crispr, by = "depmap_id")
  
  crispr$lineage[is.na(crispr$lineage)] <- crispr$tissue_type[is.na(crispr$lineage)]
  crispr$lineage[crispr$lineage == "unknown"] <- crispr$tissue_type[crispr$lineage == "unknown"]
  
  crispr$lineage <- str_replace(string = crispr$lineage, pattern = "^STOMACH$", replacement = "gastric")
  crispr$lineage <- str_replace(string = crispr$lineage, pattern = "^PANCREAS$", replacement = "pancreas")
  crispr$lineage <- str_replace(string = crispr$lineage, pattern = "^BREAST$", replacement = "breast")
  crispr$lineage <- str_replace(string = crispr$lineage, pattern = "^GASTROINTESTINAL_TRACT$", replacement = "gastrointestinal_tract")
  crispr$lineage <- str_replace(string = crispr$lineage, pattern = "^HAEMATOPOIETIC_AND_LYMPHOID_TISSUE$", replacement = "blood")
  
  ##############
  
  drug_sensitivity <- metadata %>% 
    select(depmap_id, lineage) %>% 
    right_join(drug_sensitivity, by = "depmap_id")
  
  drug_sensitivity$lineage[is.na(drug_sensitivity$lineage)] <- drug_sensitivity$tissue_type[is.na(drug_sensitivity$lineage)]
  drug_sensitivity$lineage[drug_sensitivity$lineage == "unknown"] <- drug_sensitivity$tissue_type[drug_sensitivity$lineage == "unknown"]
  
  drug_sensitivity$lineage <- str_replace(string = drug_sensitivity$lineage, pattern = "^STOMACH$", replacement = "gastric")
  drug_sensitivity$lineage <- str_replace(string = drug_sensitivity$lineage, pattern = "^PANCREAS$", replacement = "pancreas")
  drug_sensitivity$lineage <- str_replace(string = drug_sensitivity$lineage, pattern = "^BREAST$", replacement = "breast")
  drug_sensitivity$lineage <- str_replace(string = drug_sensitivity$lineage, pattern = "^GASTROINTESTINAL_TRACT$", replacement = "gastrointestinal_tract")
  drug_sensitivity$lineage <- str_replace(string = drug_sensitivity$lineage, pattern = "^HAEMATOPOIETIC_AND_LYMPHOID_TISSUE$", replacement = "blood")
  
  ##############
  
  copyNumber <- metadata %>% 
    select(depmap_id, lineage) %>% 
    right_join(copyNumber, by = "depmap_id")
  
  copyNumber$lineage[is.na(copyNumber$lineage)] <- copyNumber$tissue_type[is.na(copyNumber$lineage)]
  copyNumber$lineage[copyNumber$lineage == "unknown"] <- copyNumber$tissue_type[copyNumber$lineage == "unknown"]
  
  copyNumber$lineage <- str_replace(string = copyNumber$lineage, pattern = "^STOMACH$", replacement = "gastric")
  copyNumber$lineage <- str_replace(string = copyNumber$lineage, pattern = "^PANCREAS$", replacement = "pancreas")
  copyNumber$lineage <- str_replace(string = copyNumber$lineage, pattern = "^BREAST$", replacement = "breast")
  copyNumber$lineage <- str_replace(string = copyNumber$lineage, pattern = "^GASTROINTESTINAL_TRACT$", replacement = "gastrointestinal_tract")
  copyNumber$lineage <- str_replace(string = copyNumber$lineage, pattern = "^HAEMATOPOIETIC_AND_LYMPHOID_TISSUE$", replacement = "blood")
  copyNumber$lineage <- str_replace(string = copyNumber$lineage, pattern = "^CENTRAL_NERVOUS_SYSTEM$", replacement = "central_nervous_system")
  copyNumber$lineage <- str_replace(string = copyNumber$lineage, pattern = "^UPPER_AERODIGESTIVE_TRACT$", replacement = "esophagus")
  
  copyNumber <- copyNumber[-which(copyNumber$lineage == "MATCHED_NORMAL_TISSUE"),]
  
  ##############
  
  mutationCalls <- metadata %>% 
    select(depmap_id, lineage) %>% 
    right_join(mutationCalls, by = "depmap_id")
  
  ###########################################################################################################
  ###########################################################################################################
  
  # Evaluate the dependency of cancer cell lines to top drivers of ICGC datasets
  
  ## Evaluation of the top 20 drivers
  icgc_top20_Drivers <- lapply(icgc_driver_net_gene_tbl, function(i) {
    i %>% slice_max(order_by = primitive_driver_score, n = 20) 
  })
  
  depmap_icgc_dict <- list("Breast-AdenoCA" = "breast",
                           "Liver-HCC" = "liver",
                           "Myeloid-AML" = "blood",
                           "Skin-Melanoma" = "skin",
                           "Stomach-AdenoCA" = "gastric")
  
  ########################################
  
  ### Evaluation based on rnai
  
  #### top 20 drivers of icgc
  icgc_top20_Drivers_rnai_dep <- lapply(names(depmap_icgc_dict), function(i) {
    
    tmp_data <-
    rnai %>% filter(lineage == depmap_icgc_dict[[i]])
    
    tmp_mean_dep <- mean(tmp_data$dependency)
    
    final_data <- sapply(icgc_top20_Drivers[[i]]$gene, function(j) {
      tmp_data %>% filter(gene_name == j) %>% select(dependency) %>% unlist() %>% mean()
    })
    
    names(final_data) <- icgc_top20_Drivers[[i]]$gene
    
    list(top_20_drivers = final_data,
         general_mean = tmp_mean_dep)
  })
  
  names(icgc_top20_Drivers_rnai_dep) <- names(depmap_icgc_dict)
  
  lapply(names(icgc_top20_Drivers_rnai_dep), function(i) {
    icgc_top20_Drivers_rnai_dep[[i]]$top_20_drivers[is.na(icgc_top20_Drivers_rnai_dep[[i]]$top_20_drivers)] <<- 0
  })
  
  ##############
  
  #### 20 sample icgc genes
  icgc_top20_sample_rnai_dep <- lapply(names(depmap_icgc_dict), function(i) {
    
    tmp_data <-
      rnai %>% filter(lineage == depmap_icgc_dict[[i]])
    
    tmp_mean_dep <- mean(tmp_data$dependency)
    
    tmp_sample_genes <- icgc_driver_net_gene_tbl[[i]]$gene[which(icgc_driver_net_gene_tbl[[i]]$gene %in% tmp_data$gene_name)]
    
    set.seed(9876)
    tmp_sample_genes <- sample(tmp_sample_genes, size = 2000)
    
    final_data <- sapply(tmp_sample_genes, function(j) {
      tmp_data %>% filter(gene_name == j) %>% select(dependency) %>% unlist() %>% mean()
    })
    
    names(final_data) <- tmp_sample_genes
    
    list(top_20_sample = final_data,
         general_mean = tmp_mean_dep)
  })
  
  names(icgc_top20_sample_rnai_dep) <- names(depmap_icgc_dict)
  
  ##############
  ##### Function to work with ggplot
  is_outlier <- function(x) {
    return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.5) + 1.5 * IQR(x))
  }
  
  ##############
  
  #### Creating a unified dataframe
  icgc_rnai_dep_df <- 
    sapply(icgc_top20_Drivers_rnai_dep, function(i) names(i$top_20_drivers)) %>% 
    as.data.frame() %>% 
    rbind((sapply(icgc_top20_sample_rnai_dep, function(i) names(i$top_20_sample)) %>% 
             as.data.frame())) %>% 
    pivot_longer(everything(), 
                 names_to = "Cancer Type", 
                 values_to = "Gene",
                 cols_vary = "slowest") %>% 
    cbind(sapply(icgc_top20_Drivers_rnai_dep, function(i) i$top_20_drivers) %>% 
            as.data.frame() %>% 
            rbind((sapply(icgc_top20_sample_rnai_dep, function(i) i$top_20_sample) %>% 
                     as.data.frame())) %>% 
            pivot_longer(everything(), 
                         names_to = "Cancer Type", 
                         values_to = "Dependency",
                         cols_vary = "slowest") %>% 
            select(Dependency)) %>% 
    cbind("Category" = rep(c(rep("Top 20 Drivers", 20), rep("2000 Random Genes", 2000)), 
                                             length(depmap_icgc_dict))) %>% 
    group_by(`Cancer Type`, Category) %>% 
    mutate(is_outlier = ifelse(is_outlier(-Dependency), Gene, as.numeric(NA))) %>% 
    cbind("1st-Ranked Module Dependency" = 1) %>% 
    ungroup()
  
  icgc_rnai_dep_df$`Gene Class` <- "Unknown"
  
  
  sapply(1:length(adult_drivers_list), function(i) {
    
    tmp_index <- seq((((i - 1)*2020)+1), i*2020)

    ## define known drivers
    known_driver_index <- which(icgc_rnai_dep_df$Gene[tmp_index] %in% adult_drivers_list[[i]])
    known_driver_index <- ((i - 1)*2020) + known_driver_index
    icgc_rnai_dep_df$`Gene Class`[known_driver_index] <<- "Known Driver"
    
    ## Add top drivers that are known but have a positive dependency to the list of outliers
    known_driver_low_depend_index <- which(icgc_rnai_dep_df$Dependency[tmp_index] > 0)
    known_driver_low_depend_index <- ((i - 1)*2020) + known_driver_low_depend_index
    icgc_rnai_dep_df$is_outlier[known_driver_low_depend_index] <<- icgc_rnai_dep_df$Gene[known_driver_low_depend_index]
    
    ## Add the average of the first ranked module genes of each cancer type
    tmp_mean <- mean(rnai$dependency[which(rnai$lineage == depmap_icgc_dict[[i]] & rnai$gene_name %in% na.omit(icgc_net_1st_modules[,i]))])
    icgc_rnai_dep_df$`1st-Ranked Module Dependency`[tmp_index] <<- tmp_mean
    
  })
  
  icgc_rnai_dep_df$is_outlier[icgc_rnai_dep_df$Category == "2000 Random Genes"] <- NA
  icgc_rnai_dep_df$`Dependency Score` <- icgc_rnai_dep_df$Dependency
  icgc_rnai_dep_df$`Dependency Score`[icgc_rnai_dep_df$Category == "2000 Random Genes"] <- NA
  
  ##############
  
  #### Visualization
  library(ggrepel)
  library(ggrain)
  library(extrafont)
  extrafont::font_import()
  # Restart Rstudio
  extrafont::loadfonts()
  
  icgc_rnai_dep_plot <-
  ggplot(icgc_rnai_dep_df, 
         aes(x = `Cancer Type`, y = -Dependency, 
             fill = Category, color = Category)) +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    geom_rain(boxplot.args.pos = rlang::list2(width = 0.3, position = position_dodge(width = 0.3, preserve = "single")),
              boxplot.args = rlang::list2(lwd=0.1, alpha = 0.5, outlier.shape = NA),
              violin.args = rlang::list2(lwd=0.1, alpha = 0.5),
              point.args.pos = rlang::list2(position = position_jitter(width = 0.1, height = 0, seed = 4321)),
              point.args = rlang::list2(size = 0.6, alpha = 0.5, aes(shape = `Gene Class`, x = `Cancer Type`, y = -`Dependency Score`)),
              alpha = 0.5) + 
    scale_fill_manual(values=c("#003f5c", "#e74c3c")) +
    scale_color_manual(values=c("#003f5c", "darkred")) +
    geom_text_repel(aes(label=is_outlier), na.rm=TRUE, 
                    size = 1.8,
                    family = "Arial",
                     # color = "black",
                     # arrow = arrow(length = unit(0.015, "npc")),
                    # direction = "y", 
                    nudge_x = 0.5,
                    force_pull = 0.3,
                    force = 1,
                    segment.color = NA,
                    # position = position_jitterdodge(),
                     max.overlaps = 20, alpha = 1) +
    geom_point(data = icgc_rnai_dep_df, shape = 15, size = 0.9,  alpha = 0.1, color = "#06f32b", 
               aes(x = `Cancer Type`, y = -`1st-Ranked Module Dependency`, group = Category)) +
    coord_flip() +
    theme_classic() +
    theme(text = element_text(size = 7))
  
  icgc_rnai_dep_plot
  
  ggsave(filename = "Results/DepMap/icgc_rnai_dep_plot.pdf", 
         plot = icgc_rnai_dep_plot, device = "pdf", width = 18, height = 8.4, units = "cm")
  
  
  ###########################################################################################################
  ###########################################################################################################
  
  # Evaluate the dependency of cancer cell lines to top risk genes of pediatric cancer datasets
  
  ## Evaluation of the top 20 Risk Genes
  ped_top20_risk <- lapply(risk_net_gene_tbl[-c(1,9)], function(i) {
    i %>% slice_max(order_by = primitive_risk_score, n = 20) 
  })
  
  ## Defining the dictionary of our pediatric cancer types and depmap lineages
  depmap_ped_dict <- list("DMG" = "central_nervous_system",
                           "HGG" = "central_nervous_system",
                           "Leukaemia" = "blood",
                           "NBL" = "peripheral_nervous_system",
                           "Sarcoma" = "soft_tissue",
                           "CNS embryonal" = "central_nervous_system",
                           "EPD" = "central_nervous_system",
                           "Lymphoma" = "lymphocyte",
                           "MB" = "central_nervous_system",
                           "Rhabdoid" = "kidney")
  
  ## Defining the dictionary of our pediatric cancer types and the known pediatric drivers
  ped_driver_list <- list("DMG" = hgg_driver_genes,
                          "HGG" = hgg_driver_genes,
                          "Leukaemia" = leukaemia_driver_genes,
                          "NBL" = nbl_driver_genes,
                          "Sarcoma" = sarcoma_driver_genes,
                          "CNS embryonal" = embryonal_cns_driver_genes,
                          "EPD" = epd_driver_genes,
                          "Lymphoma" = lymphoma_driver_genes,
                          "MB" = mb_driver_genes,
                          "Rhabdoid" = rhabdoid_driver_genes)
  
  ########################################
  
  ### Evaluation based on rnai
  
  #### top 20 risk genes of pediatric cancer
  ped_top20_risk_rnai_dep <- lapply(names(depmap_ped_dict), function(i) {
    
    tmp_data <-
      rnai %>% filter(lineage == depmap_ped_dict[[i]])
    
    tmp_mean_dep <- mean(tmp_data$dependency)
    
    final_data <- sapply(ped_top20_risk[[i]]$gene, function(j) {
      tmp_data %>% filter(gene_name == j) %>% select(dependency) %>% unlist() %>% mean()
    })
    
    names(final_data) <- ped_top20_risk[[i]]$gene
    
    list(top_20_risk = final_data,
         general_mean = tmp_mean_dep)
  })
  
  names(ped_top20_risk_rnai_dep) <- names(depmap_ped_dict)
  
  lapply(names(ped_top20_risk_rnai_dep), function(i) {
    ped_top20_risk_rnai_dep[[i]]$top_20_risk[is.na(ped_top20_risk_rnai_dep[[i]]$top_20_risk)] <<- 0
  })
  
  ##############
  
  #### 20 sample ped genes
  ped_top20_sample_rnai_dep <- lapply(names(depmap_ped_dict), function(i) {
    
    tmp_data <-
      rnai %>% filter(lineage == depmap_ped_dict[[i]])
    
    tmp_mean_dep <- mean(tmp_data$dependency)
    
    tmp_sample_genes <- risk_net_gene_tbl[-c(1,9)][[i]]$gene[which(risk_net_gene_tbl[-c(1,9)][[i]]$gene %in% tmp_data$gene_name)]
    
    set.seed(4321)
    tmp_sample_genes <- sample(tmp_sample_genes, size = 2000)
    
    final_data <- sapply(tmp_sample_genes, function(j) {
      tmp_data %>% filter(gene_name == j) %>% select(dependency) %>% unlist() %>% mean()
    })
    
    names(final_data) <- tmp_sample_genes
    
    list(top_20_sample = final_data,
         general_mean = tmp_mean_dep)
  })
  
  names(ped_top20_sample_rnai_dep) <- names(depmap_ped_dict)
  
  ##############
  
  #### Creating a unified dataframe
  ped_rnai_risk_dep_df <- 
    sapply(ped_top20_risk_rnai_dep, function(i) names(i$top_20_risk)) %>% 
    as.data.frame() %>% 
    rbind((sapply(ped_top20_sample_rnai_dep, function(i) names(i$top_20_sample)) %>% 
             as.data.frame())) %>% 
    pivot_longer(everything(), 
                 names_to = "Cancer Type", 
                 values_to = "Gene",
                 cols_vary = "slowest") %>% 
    cbind(sapply(ped_top20_risk_rnai_dep, function(i) i$top_20_risk) %>% 
            as.data.frame() %>% 
            rbind((sapply(ped_top20_sample_rnai_dep, function(i) i$top_20_sample) %>% 
                     as.data.frame())) %>% 
            pivot_longer(everything(), 
                         names_to = "Cancer Type", 
                         values_to = "Dependency",
                         cols_vary = "slowest") %>% 
            select(Dependency)) %>% 
    cbind("Category" = rep(c(rep("Top 20 Risk Genes", 20), rep("2000 Random Genes", 2000)), 
                           length(depmap_ped_dict))) %>% 
    group_by(`Cancer Type`, Category) %>% 
    mutate(is_outlier = ifelse(is_outlier(-Dependency), Gene, as.character(NA))) %>% 
    cbind("1st-Ranked Module Dependency" = 1)
  
  ped_rnai_risk_dep_df$`Gene Class` <- "Unknown"
  
  
  sapply(1:length(ped_driver_list), function(i) {
    
    tmp_index <- seq((((i - 1)*2020)+1), i*2020)
    
    ## define known drivers
    known_driver_index <- which(ped_rnai_risk_dep_df$Gene[tmp_index] %in% ped_driver_list[[i]])
    known_driver_index <- ((i - 1)*2020) + known_driver_index
    ped_rnai_risk_dep_df$`Gene Class`[known_driver_index] <<- "Known Driver"
    
    ## Add top drivers that are known but have a positive dependency, to the list of outliers
    known_driver_low_depend_index <- which(ped_rnai_risk_dep_df$Dependency[tmp_index] > 0)
    known_driver_low_depend_index <- ((i - 1)*2020) + known_driver_low_depend_index
    ped_rnai_risk_dep_df$is_outlier[known_driver_low_depend_index] <<- ped_rnai_risk_dep_df$Gene[known_driver_low_depend_index]
    
    ## Add the average of the first ranked module genes of each cancer type
    tmp_mean <- mean(rnai$dependency[which(rnai$lineage == depmap_ped_dict[[i]] & rnai$gene_name %in% na.omit(risk_net_1st_modules[,-c(1,9)][,i]))])
    ped_rnai_risk_dep_df$`1st-Ranked Module Dependency`[tmp_index] <<- tmp_mean
    
  })
  
  ped_rnai_risk_dep_df$is_outlier[ped_rnai_risk_dep_df$Category == "2000 Random Genes"] <- NA
  ped_rnai_risk_dep_df$`Dependency Score` <- ped_rnai_risk_dep_df$Dependency
  ped_rnai_risk_dep_df$`Dependency Score`[ped_rnai_risk_dep_df$Category == "2000 Random Genes"] <- NA
  
  ##############
  
  #### Visualization
  library(ggrepel)
  library(ggrain)
  library(extrafont)
  extrafont::font_import()
  # Restart Rstudio
  extrafont::loadfonts()
  
  ped_rnai_risk_dep_plot <-
    ggplot(ped_rnai_risk_dep_df, 
           aes(x = `Cancer Type`, y = -Dependency, 
               fill = Category, color = Category)) +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    geom_rain(boxplot.args.pos = rlang::list2(width = 0.3, position = position_dodge(width = 0.3, preserve = "single")),
              boxplot.args = rlang::list2(lwd=0.1, alpha = 0.5, outlier.shape = NA),
              violin.args = rlang::list2(lwd=0.1, alpha = 0.5),
              point.args.pos = rlang::list2(position = position_jitter(width = 0.1, height = 0, seed = 4321)),
              point.args = rlang::list2(size = 0.6, alpha = 0.5, aes(shape = `Gene Class`, x = `Cancer Type`, y = -`Dependency Score`)),
              alpha = 0.5) + 
    scale_fill_manual(values=c("#003f5c", "#e74c3c")) +
    scale_color_manual(values=c("#003f5c", "darkred")) +
    geom_text_repel(aes(label=is_outlier), na.rm=TRUE, 
                    size = 1.8,
                    family = "Arial",
                    # color = "black",
                    # arrow = arrow(length = unit(0.015, "npc")),
                    # direction = "y", 
                    nudge_x = 0.5,
                    force_pull = 0.3,
                    force = 1,
                    segment.color = NA,
                    # position = position_jitterdodge(),
                    max.overlaps = 20, alpha = 1) +
    geom_point(data = ped_rnai_risk_dep_df, shape = 15, size = 0.9,  alpha = 0.1, color = "#06f32b", 
               aes(x = `Cancer Type`, y = -`1st-Ranked Module Dependency`, group = Category)) +
    coord_flip() +
    theme_classic() +
    theme(text = element_text(size = 7))
  
  ped_rnai_risk_dep_plot
  
  ggsave(filename = "Results/DepMap/ped_rnai_risk_dep_plot.pdf", 
         plot = ped_rnai_risk_dep_plot, device = "pdf", width = 18, height = 13.6, units = "cm")
  
  ###########################################################################################################
  ###########################################################################################################
  
  # Evaluate the dependency of cancer cell lines to top Driver Genes of pediatric cancer datasets
  
  ## Evaluation of the top 20 Driver Genes
  ped_top20_driver <- lapply(driver_marker_net_gene_tbl[names(ped_top20_risk_rnai_dep)], function(i) {
    i %>% slice_max(order_by = primitive_driver_score, n = 20) 
  })
  
  ########################################
  
  ### Evaluation based on rnai
  
  #### top 20 Driver Genes of pediatric cancer
  ped_top20_driver_rnai_dep <- lapply(names(depmap_ped_dict), function(i) {
    
    tmp_data <-
      rnai %>% filter(lineage == depmap_ped_dict[[i]])
    
    tmp_mean_dep <- mean(tmp_data$dependency)
    
    final_data <- sapply(ped_top20_driver[[i]]$gene, function(j) {
      tmp_data %>% filter(gene_name == j) %>% select(dependency) %>% unlist() %>% mean()
    })
    
    names(final_data) <- ped_top20_driver[[i]]$gene
    
    list(top_20_driver = final_data,
         general_mean = tmp_mean_dep)
  })
  
  names(ped_top20_driver_rnai_dep) <- names(depmap_ped_dict)
  
  lapply(names(ped_top20_driver_rnai_dep), function(i) {
    ped_top20_driver_rnai_dep[[i]]$top_20_driver[is.na(ped_top20_driver_rnai_dep[[i]]$top_20_driver)] <<- 0
  })
  
  ##############
  
  #### 20 sample ped genes
  ped_top20_driver_sample_rnai_dep <- lapply(names(depmap_ped_dict), function(i) {
    
    tmp_data <-
      rnai %>% filter(lineage == depmap_ped_dict[[i]])
    
    tmp_mean_dep <- mean(tmp_data$dependency)
    
    tmp_sample_genes <- driver_marker_net_gene_tbl[names(ped_top20_risk_rnai_dep)][[i]]$gene[which(driver_marker_net_gene_tbl[names(ped_top20_risk_rnai_dep)][[i]]$gene %in% tmp_data$gene_name)]
    
    set.seed(4321)
    tmp_sample_genes <- sample(tmp_sample_genes, size = 2000)
    
    final_data <- sapply(tmp_sample_genes, function(j) {
      tmp_data %>% filter(gene_name == j) %>% select(dependency) %>% unlist() %>% mean()
    })
    
    names(final_data) <- tmp_sample_genes
    
    list(top_20_sample = final_data,
         general_mean = tmp_mean_dep)
  })
  
  names(ped_top20_driver_sample_rnai_dep) <- names(depmap_ped_dict)
  
  ##############
  
  #### Creating a unified dataframe
  ped_rnai_driver_dep_df <- 
    sapply(ped_top20_driver_rnai_dep, function(i) names(i$top_20_driver)) %>% 
    as.data.frame() %>% 
    rbind((sapply(ped_top20_driver_sample_rnai_dep, function(i) names(i$top_20_sample)) %>% 
             as.data.frame())) %>% 
    pivot_longer(everything(), 
                 names_to = "Cancer Type", 
                 values_to = "Gene",
                 cols_vary = "slowest") %>% 
    cbind(sapply(ped_top20_driver_rnai_dep, function(i) i$top_20_driver) %>% 
            as.data.frame() %>% 
            rbind((sapply(ped_top20_driver_sample_rnai_dep, function(i) i$top_20_sample) %>% 
                     as.data.frame())) %>% 
            pivot_longer(everything(), 
                         names_to = "Cancer Type", 
                         values_to = "Dependency",
                         cols_vary = "slowest") %>% 
            select(Dependency)) %>% 
    cbind("Category" = rep(c(rep("Top 20 Drivers", 20), rep("2000 Random Genes", 2000)), 
                           length(depmap_ped_dict))) %>% 
    group_by(`Cancer Type`, Category) %>% 
    mutate(is_outlier = ifelse(is_outlier(-Dependency), Gene, as.character(NA))) %>% 
    cbind("1st-Ranked Module Dependency" = 1)
  
  ped_rnai_driver_dep_df$`Gene Class` <- "Unknown"
  
  
  sapply(1:length(ped_driver_list), function(i) {
    
    tmp_index <- seq((((i - 1)*2020)+1), i*2020)
    
    ## define known drivers
    known_driver_index <- which(ped_rnai_driver_dep_df$Gene[tmp_index] %in% ped_driver_list[[i]])
    known_driver_index <- ((i - 1)*2020) + known_driver_index
    ped_rnai_driver_dep_df$`Gene Class`[known_driver_index] <<- "Known Driver"
    
    ## Add top drivers that are known but have a positive dependency, to the list of outliers
    known_driver_low_depend_index <- which(ped_rnai_driver_dep_df$Dependency[tmp_index] > 0)
    known_driver_low_depend_index <- ((i - 1)*2020) + known_driver_low_depend_index
    ped_rnai_driver_dep_df$is_outlier[known_driver_low_depend_index] <<- ped_rnai_driver_dep_df$Gene[known_driver_low_depend_index]
    
    ## Add the average of the first ranked module genes of each cancer type
    tmp_mean <- mean(rnai$dependency[which(rnai$lineage == depmap_ped_dict[[i]] & rnai$gene_name %in% na.omit(driver_net_1st_modules[,names(ped_top20_risk_rnai_dep)][,i]))])
    ped_rnai_driver_dep_df$`1st-Ranked Module Dependency`[tmp_index] <<- tmp_mean
    
  })
  
  ped_rnai_driver_dep_df$is_outlier[ped_rnai_driver_dep_df$Category == "2000 Random Genes"] <- NA
  ped_rnai_driver_dep_df$`Dependency Score` <- ped_rnai_driver_dep_df$Dependency
  ped_rnai_driver_dep_df$`Dependency Score`[ped_rnai_driver_dep_df$Category == "2000 Random Genes"] <- NA
  
  ##############
  
  #### Visualization
  
  ped_rnai_driver_dep_plot <-
    ggplot(ped_rnai_driver_dep_df, 
           aes(x = `Cancer Type`, y = -Dependency, 
               fill = Category, color = Category)) +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    geom_rain(boxplot.args.pos = rlang::list2(width = 0.3, position = position_dodge(width = 0.3, preserve = "single")),
              boxplot.args = rlang::list2(lwd=0.1, alpha = 0.5, outlier.shape = NA),
              violin.args = rlang::list2(lwd=0.1, alpha = 0.5),
              point.args.pos = rlang::list2(position = position_jitter(width = 0.1, height = 0, seed = 4321)),
              point.args = rlang::list2(size = 0.6, alpha = 0.5, aes(shape = `Gene Class`, x = `Cancer Type`, y = -`Dependency Score`)),
              alpha = 0.5) + 
    scale_fill_manual(values=c("#003f5c", "#e74c3c")) +
    scale_color_manual(values=c("#003f5c", "darkred")) +
    geom_text_repel(aes(label=is_outlier), na.rm=TRUE, 
                    size = 1.8,
                    family = "Arial",
                    # color = "black",
                    # arrow = arrow(length = unit(0.015, "npc")),
                    # direction = "y", 
                    nudge_x = 0.5,
                    force_pull = 0.3,
                    force = 1,
                    segment.color = NA,
                    # position = position_jitterdodge(),
                    max.overlaps = 20, alpha = 1) +
    geom_point(data = ped_rnai_driver_dep_df, shape = 15, size = 0.9,  alpha = 0.1, color = "#06f32b", 
               aes(x = `Cancer Type`, y = -`1st-Ranked Module Dependency`, group = Category)) +
    coord_flip() +
    theme_classic() +
    theme(text = element_text(size = 7))
  
  ped_rnai_driver_dep_plot
  
  ggsave(filename = "Results/DepMap/ped_rnai_driver_dep_plot.pdf", 
         plot = ped_rnai_driver_dep_plot, device = "pdf", width = 18, height = 13, units = "cm")
  
  ###########################################################################################################
  ###########################################################################################################
  
  # Evaluate the dependency of cancer cell lines to top 20 Risk and Driver Genes of pediatric pan-cancer datasets
  
  pan_cancer_driver_top_20
  pan_cancer_risk_top_20
  
  ## Creating a combined table of top 20 pan-cancer risk and driver genes
  ped_pan_cancer_top_20 <- tibble("Gene" = c(pan_cancer_driver_top_20, pan_cancer_risk_top_20),
                                  "Category" = rep(c("Driver", "Risk"), each = 20))
  
  ped_pan_cancer_top_20_dup_genes <- ped_pan_cancer_top_20$Gene[which(duplicated(ped_pan_cancer_top_20$Gene))]
  ped_pan_cancer_top_20 <- ped_pan_cancer_top_20[-which(duplicated(ped_pan_cancer_top_20$Gene)),]
  ped_pan_cancer_top_20$Category[which(ped_pan_cancer_top_20$Gene %in% ped_pan_cancer_top_20_dup_genes)] <- "Risk and Driver"
  
  ########################################
  
  ### Evaluation based on rnai
  
  pan_cancer_ped_top20_rnai_dep <-
  rnai %>% 
    filter(gene_name %in% ped_pan_cancer_top_20$Gene) %>% 
    group_by(lineage, gene_name) %>% 
    summarize("Dependency" = mean(dependency))
  
  colnames(pan_cancer_ped_top20_rnai_dep) <- c("Lineage", "Gene", "Dependency")
  pan_cancer_ped_top20_rnai_dep$Category <- ped_pan_cancer_top_20$Category[match(pan_cancer_ped_top20_rnai_dep$Gene, ped_pan_cancer_top_20$Gene)]
  
  ##############
  
  ## Define known drivers
  
  pan_cancer_ped_top20_rnai_dep$`Gene Class` <- "Unknown"
  pan_cancer_ped_top20_rnai_dep$`Gene Class`[which(pan_cancer_ped_top20_rnai_dep$Gene %in% pan_cancer_ped_drivers)] <- "Known Driver"
  
  
  ##############
  
  #### Visualization
  
  pan_cancer_ped_top20_rnai_dep_plot_color_legend_order <- 
    pan_cancer_ped_top20_rnai_dep %>% ungroup() %>% select(Gene, `Gene Class`) %>% distinct(Gene, .keep_all = T) %>% arrange(`Gene Class`)
  
  pan_cancer_ped_top20_rnai_dep_plot <-
    ggplot(pan_cancer_ped_top20_rnai_dep, 
           aes(x = Lineage, y = -Dependency, 
               color = Gene, shape = Category)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_point(aes(fill = Gene, stroke = ifelse(`Gene Class` == "Unknown", 0.3, NA)), color = 'black', size = 1.5, alpha = 0.5) +
    scale_shape_manual(values = c(22,23,21)) + # for the points to be able to have fill color their shape should be between 21-25
    scale_fill_manual(values = setNames(fishualize::fish(n = length(unique(pan_cancer_ped_top20_rnai_dep$Gene))), nm = unique(pan_cancer_ped_top20_rnai_dep$Gene)),
                      breaks = pan_cancer_ped_top20_rnai_dep_plot_color_legend_order %>% select(Gene) %>% unlist() %>% unname(),
                      guide = "legend") +
    scale_color_manual(values = setNames(fishualize::fish(n = length(unique(pan_cancer_ped_top20_rnai_dep$Gene))), nm = unique(pan_cancer_ped_top20_rnai_dep$Gene)), 
                       breaks = pan_cancer_ped_top20_rnai_dep_plot_color_legend_order %>% select(Gene) %>% unlist() %>% unname(),
                       guide = "legend") +
    coord_flip() +
    labs(y = "-Dependency") +
    theme_classic() +
    theme(text = element_text(size = 7))
  
  pan_cancer_ped_top20_rnai_dep_plot

  ggsave(filename = "Results/DepMap/pan_cancer_ped_top20_rnai_dep_plot.pdf", 
         plot = pan_cancer_ped_top20_rnai_dep_plot, device = "pdf", width = 18, height = 10, units = "cm")
  
  #############
  
  pan_cancer_ped_top20_rnai_dep_legend_plot <-
    ggplot(pan_cancer_ped_top20_rnai_dep, 
           aes(x = Lineage, y = -Dependency, 
               color = Gene, shape = Category)) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
    geom_point(aes(fill = Gene, stroke = ifelse(`Gene Class` == "Unknown", 0.3, NA)), size = 1.5, alpha = 0.5) +
    scale_shape_manual(values = c(22,23,21)) + # for the points to be able to have fill color their shape should be between 21-25
    scale_fill_manual(values = setNames(fishualize::fish(n = length(unique(pan_cancer_ped_top20_rnai_dep$Gene))), nm = unique(pan_cancer_ped_top20_rnai_dep$Gene)),
                      breaks = pan_cancer_ped_top20_rnai_dep_plot_color_legend_order %>% select(Gene) %>% unlist() %>% unname(),
                      guide = "legend") +
    scale_color_manual(values = setNames(fishualize::fish(n = length(unique(pan_cancer_ped_top20_rnai_dep$Gene))), nm = unique(pan_cancer_ped_top20_rnai_dep$Gene)), 
                       breaks = pan_cancer_ped_top20_rnai_dep_plot_color_legend_order %>% select(Gene) %>% unlist() %>% unname(),
                       guide = "legend") +
    coord_flip() +
    labs(y = "-Dependency") +
    theme_classic() +
    theme(text = element_text(size = 7))
  
  pan_cancer_ped_top20_rnai_dep_legend_plot
  
  ggsave(filename = "Results/DepMap/pan_cancer_ped_top20_rnai_dep_legend_plot.pdf", 
         plot = pan_cancer_ped_top20_rnai_dep_legend_plot, device = "pdf", width = 18, height = 10, units = "cm")
  
  ###########################################################################################################
  ###########################################################################################################
  
  # Find the Drug sensitivity for the top genes of ICGC genes with Dependency < -0.1
  
  ## Creating a unified dataframe
  icgc_top_candidate_drug_sensitivity <- 
    sapply(icgc_top20_Drivers_rnai_dep, function(i) names(i$top_20_drivers)) %>% 
    as.data.frame() %>% 
    pivot_longer(everything(), 
                 names_to = "Cancer Type", 
                 values_to = "Gene",
                 cols_vary = "slowest") %>% 
    cbind(sapply(icgc_top20_Drivers_rnai_dep, function(i) i$top_20_drivers) %>% 
            as.data.frame() %>% 
            pivot_longer(everything(), 
                         names_to = "Cancer Type", 
                         values_to = "Dependency",
                         cols_vary = "slowest") %>% 
            select(Dependency)) %>% 
    group_by(`Cancer Type`) %>% 
    mutate(is_outlier = ifelse(is_outlier(-Dependency), Gene, as.numeric(NA))) %>% 
      cbind("Gene Class" = "Unknown") 
    
    sapply(1:length(adult_drivers_list), function(i) {
      
      tmp_index <- seq((((i - 1)*20)+1), i*20)
      
      ## define known drivers
      known_driver_index <- which(icgc_top_candidate_drug_sensitivity$Gene[tmp_index] %in% adult_drivers_list[[i]])
      known_driver_index <- ((i - 1)*20) + known_driver_index
      icgc_top_candidate_drug_sensitivity$`Gene Class`[known_driver_index] <<- "Known Driver"
      
    })
    
    ### Filter the data
    icgc_top_candidate_drug_sensitivity <- icgc_top_candidate_drug_sensitivity %>% 
      filter(Dependency <= -0.1)
    
    ###############
    
    icgc_top_candidate_drug_sensitivity <- cbind(icgc_top_candidate_drug_sensitivity,
                                                 sapply(1:nrow(icgc_top_candidate_drug_sensitivity), function(i) {
                                                   drug_sensitivity %>% 
                                                     filter(target %in% icgc_top_candidate_drug_sensitivity$Gene[i] & lineage == depmap_icgc_dict[[icgc_top_candidate_drug_sensitivity$`Cancer Type`[i]]]) %>% 
                                                     group_by('Compound Name' = name) %>% 
                                                     summarise("Mean Dependency" = mean(dependency))
                                                 }) %>% t() %>% 
                                                   as.data.frame() %>% 
                                                   apply(MARGIN = 1, function(i) {
                                                       tmp_compound = unlist(i[1])
                                                       tmp_mean_dep = unlist(i[2])
                                                       tmp_compound <- tmp_compound[which(tmp_mean_dep < 0)]
                                                       tmp_mean_dep <- tmp_mean_dep[which(tmp_mean_dep < 0)]
                                                       c(paste0(tmp_compound, collapse = ", "),
                                                         paste0(tmp_mean_dep, collapse = ", "))
                                                   }) %>% t() %>% 
                                                   as.data.frame()
                                                 )
    
    colnames(icgc_top_candidate_drug_sensitivity)[c(6,7)] <- c('Compound Name', 'Mean Dependency')

    icgc_top_candidate_drug_sensitivity <- icgc_top_candidate_drug_sensitivity %>% 
      filter(`Mean Dependency` != "")
    
    icgc_top_candidate_drug_sensitivity$Group <- "Adulthood Cancer Driver"
    
    ###########################################################################################################
    ###########################################################################################################
    
    # Find the Drug sensitivity for the top genes of pediatric risk genes with Dependency < -0.1
    
    ## Creating a unified dataframe
    ped_risk_top_candidate_drug_sensitivity <- 
      sapply(ped_top20_risk_rnai_dep, function(i) names(i$top_20_risk)) %>% 
      as.data.frame() %>% 
      pivot_longer(everything(), 
                   names_to = "Cancer Type", 
                   values_to = "Gene",
                   cols_vary = "slowest") %>% 
      cbind(sapply(ped_top20_risk_rnai_dep, function(i) i$top_20_risk) %>% 
              as.data.frame() %>% 
              pivot_longer(everything(), 
                           names_to = "Cancer Type", 
                           values_to = "Dependency",
                           cols_vary = "slowest") %>% 
              select(Dependency)) %>% 
      group_by(`Cancer Type`) %>% 
      mutate(is_outlier = ifelse(is_outlier(-Dependency), Gene, as.character(NA))) %>% 
      cbind("Gene Class" = "Unknown") 
    
    sapply(1:length(ped_driver_list), function(i) {
      
      tmp_index <- seq((((i - 1)*20)+1), i*20)
      
      ## define known drivers
      known_driver_index <- which(ped_risk_top_candidate_drug_sensitivity$Gene[tmp_index] %in% ped_driver_list[[i]])
      known_driver_index <- ((i - 1)*20) + known_driver_index
      ped_risk_top_candidate_drug_sensitivity$`Gene Class`[known_driver_index] <<- "Known Driver"
      
    })
    
    ### Filter the data
    ped_risk_top_candidate_drug_sensitivity <- ped_risk_top_candidate_drug_sensitivity %>% 
      filter(Dependency <= -0.1)
    
    ###############
    
    ped_risk_top_candidate_drug_sensitivity <- cbind(ped_risk_top_candidate_drug_sensitivity,
                                                     sapply(1:nrow(ped_risk_top_candidate_drug_sensitivity), function(i) {
                                                       drug_sensitivity %>% 
                                                         filter(target %in% ped_risk_top_candidate_drug_sensitivity$Gene[i] & lineage == depmap_ped_dict[[ped_risk_top_candidate_drug_sensitivity$`Cancer Type`[i]]]) %>% 
                                                         group_by('Compound Name' = name) %>% 
                                                         summarise("Mean Dependency" = mean(dependency))
                                                     }) %>% t() %>% as.data.frame() %>% 
                                                       apply(MARGIN = 1, function(i) {
                                                         tmp_compound = unlist(i[1])
                                                         tmp_mean_dep = unlist(i[2])
                                                         tmp_compound <- tmp_compound[which(tmp_mean_dep < 0)]
                                                         tmp_mean_dep <- tmp_mean_dep[which(tmp_mean_dep < 0)]
                                                         c(paste0(tmp_compound, collapse = ", "),
                                                           paste0(tmp_mean_dep, collapse = ", "))
                                                       }) %>% t() %>% 
                                                       as.data.frame()
                                                     )
    
    
    
    colnames(ped_risk_top_candidate_drug_sensitivity)[c(6,7)] <- c('Compound Name', 'Mean Dependency')
    
    ped_risk_top_candidate_drug_sensitivity <- ped_risk_top_candidate_drug_sensitivity %>% 
      filter(`Mean Dependency` != "")
    
    ped_risk_top_candidate_drug_sensitivity$Group <- "Pediatric Cancer Risk Gene"
    
    ###########################################################################################################
    ###########################################################################################################
    
    # Find the Drug sensitivity for the top genes of pediatric driver genes with Dependency < -0.1
    
    ## Creating a unified dataframe
    ped_driver_top_candidate_drug_sensitivity <- 
      sapply(ped_top20_driver_rnai_dep, function(i) names(i$top_20_driver)) %>% 
      as.data.frame() %>% 
      pivot_longer(everything(), 
                   names_to = "Cancer Type", 
                   values_to = "Gene",
                   cols_vary = "slowest") %>% 
      cbind(sapply(ped_top20_driver_rnai_dep, function(i) i$top_20_driver) %>% 
              as.data.frame() %>% 
              pivot_longer(everything(), 
                           names_to = "Cancer Type", 
                           values_to = "Dependency",
                           cols_vary = "slowest") %>% 
              select(Dependency)) %>% 
      group_by(`Cancer Type`) %>% 
      mutate(is_outlier = ifelse(is_outlier(-Dependency), Gene, as.character(NA))) %>% 
      cbind("Gene Class" = "Unknown") 
    
    sapply(1:length(ped_driver_list), function(i) {
      
      tmp_index <- seq((((i - 1)*20)+1), i*20)
      
      ## define known drivers
      known_driver_index <- which(ped_driver_top_candidate_drug_sensitivity$Gene[tmp_index] %in% ped_driver_list[[i]])
      known_driver_index <- ((i - 1)*20) + known_driver_index
      ped_driver_top_candidate_drug_sensitivity$`Gene Class`[known_driver_index] <<- "Known Driver"
      
    })
    
    ### Filter the data
    ped_driver_top_candidate_drug_sensitivity <- ped_driver_top_candidate_drug_sensitivity %>% 
      filter(Dependency <= -0.1)
    
    ###############
    
    ped_driver_top_candidate_drug_sensitivity <- cbind(ped_driver_top_candidate_drug_sensitivity,
                                                       sapply(1:nrow(ped_driver_top_candidate_drug_sensitivity), function(i) {
                                                         drug_sensitivity %>% 
                                                           filter(target %in% ped_driver_top_candidate_drug_sensitivity$Gene[i] & lineage == depmap_ped_dict[[ped_driver_top_candidate_drug_sensitivity$`Cancer Type`[i]]]) %>% 
                                                           group_by('Compound Name' = name) %>% 
                                                           summarise("Mean Dependency" = mean(dependency))
                                                       }) %>% t() %>% as.data.frame() %>% 
                                                         apply(MARGIN = 1, function(i) {
                                                           tmp_compound = unlist(i[1])
                                                           tmp_mean_dep = unlist(i[2])
                                                           tmp_compound <- tmp_compound[which(tmp_mean_dep < 0)]
                                                           tmp_mean_dep <- tmp_mean_dep[which(tmp_mean_dep < 0)]
                                                           c(paste0(tmp_compound, collapse = ", "),
                                                             paste0(tmp_mean_dep, collapse = ", "))
                                                         }) %>% t() %>% 
                                                         as.data.frame()
    )
    
    
    
    colnames(ped_driver_top_candidate_drug_sensitivity)[c(6,7)] <- c('Compound Name', 'Mean Dependency')
    
    ped_driver_top_candidate_drug_sensitivity <- ped_driver_top_candidate_drug_sensitivity %>% 
      filter(`Mean Dependency` != "")
    
    ped_driver_top_candidate_drug_sensitivity$Group <- "Pediatric Cancer Driver"
    
    ###########################################################################################################
    ###########################################################################################################
    
    # Creating a combined drug sensitivity table for network visualization
    
    top_candidate_drug_sensitivity_tbl <- rbind(icgc_top_candidate_drug_sensitivity,
                                                ped_risk_top_candidate_drug_sensitivity,
                                                ped_driver_top_candidate_drug_sensitivity)
    
    top_candidate_drug_sensitivity_net <- top_candidate_drug_sensitivity_tbl %>% 
      select("Gene", "Cancer Type", "Dependency", "Gene Class", "Group") %>% 
      dplyr::mutate(Dependency = -Dependency) %>% 
      rename(all_of(c("Source" = "Gene", "Target" = "Cancer Type", "-(Mean Dependency)" = "Dependency"))) %>% 
      mutate("Node Type" = "") %>% 
      rbind(
    lapply(1:nrow(top_candidate_drug_sensitivity_tbl), function(i) {
      tmp_drug_gene <- tibble(Source = (str_split(top_candidate_drug_sensitivity_tbl$`Compound Name`[i], pattern = ", ") %>% unlist()),
                              Target = top_candidate_drug_sensitivity_tbl$Gene[i],
                              `-(Mean Dependency)` = -1*as.numeric((str_split(top_candidate_drug_sensitivity_tbl$`Mean Dependency`[i], pattern = ", ") %>% unlist())),
                              `Gene Class` = "",
                              Group = "") %>% 
        unique() %>% 
        mutate("Node Type" = "") %>% 
        ungroup()   %>% 
        group_by(Source, Target) %>%
        mutate(`-(Mean Dependency)` = mean(as.numeric(`-(Mean Dependency)`))) %>%
        distinct() %>%
        ungroup()
      
      tmp_drug_gene
      
    }) %>% do.call(rbind, .)
      ) %>% 
      unique() %>%
      ungroup()
    
    ########################################
    
    ## Add the cancer node type
    sapply(1:length(unique(top_candidate_drug_sensitivity_tbl$`Cancer Type`)), function(i) {
      top_candidate_drug_sensitivity_net <<- top_candidate_drug_sensitivity_net %>% 
        rbind(c(unique(top_candidate_drug_sensitivity_tbl$`Cancer Type`)[i],
                rep("", ncol(top_candidate_drug_sensitivity_net) - 2), 
                "Cancer"))
    })
    
    #############
    
    ## Add the target node type
    sapply(1:length(unique(top_candidate_drug_sensitivity_tbl$Gene)), function(i) {
      top_candidate_drug_sensitivity_net <<- top_candidate_drug_sensitivity_net %>% 
        rbind(c(unique(top_candidate_drug_sensitivity_tbl$Gene)[i],
                rep("", ncol(top_candidate_drug_sensitivity_net) - 2), 
                "Target"))
    })
    
    #############
    
    ## Add the target node type
    sapply(1:length(unique(unlist(str_split(top_candidate_drug_sensitivity_tbl$`Compound Name`, pattern = ", ")))), function(i) {
      top_candidate_drug_sensitivity_net <<- top_candidate_drug_sensitivity_net %>% 
        rbind(c(unique(unlist(str_split(top_candidate_drug_sensitivity_tbl$`Compound Name`, pattern = ", ")))[i],
                rep("", ncol(top_candidate_drug_sensitivity_net) - 2), 
                "Compound"))
    })
    
    #############
    
    top_candidate_drug_sensitivity_net <- top_candidate_drug_sensitivity_net %>% 
      group_by(Source, Target) %>%
      mutate(`-(Mean Dependency)` = mean(as.numeric(`-(Mean Dependency)`))) %>%
      distinct() %>%
      ungroup()
    
    ########################################
    
    top_candidate_drug_sensitivity_net$`Node Type`[is.na(top_candidate_drug_sensitivity_net$`Node Type`)] <- ""
    
    write_csv(x = top_candidate_drug_sensitivity_net, 
              file = "Results/DepMap/top_candidate_drug_sensitivity_net.csv")
