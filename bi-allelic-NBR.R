library(optparse)
library(tidyverse)
library(MASS)
library(reshape2)
library(readr)


parser <- OptionParser()
parser <- add_option(parser, '--dna', action = 'store', type = 'character')
parser <- add_option(parser, '--mrna', action = 'store', type = 'character')
parser <- add_option(parser, '--annotation', action = 'store', type = 'character')
parser <- add_option(parser, '--min_count', action = 'store', type = 'numeric',default=15)
parser <- add_option(parser, '--barcode_allele', action = 'store', type = 'character')
parser <- add_option(parser, '--output', action = 'store', type = 'character')

parser <- parse_args(parser)
min_count=parser$min_count

message("Reading Barcode Allele File")
barcode_allele <- read_tsv(parser$barcode_allele, col_type = 'ccci') 
var_list <- barcode_allele %>% group_by(mutation) %>% count(type) %>% dcast(mutation ~ type,value.var = "n") %>% filter(`1` >= min_count & `2` >= min_count) %>% pull(mutation)
barcode_allele_filtered <- barcode_allele %>% filter(mutation %in% var_list)

message("Reading mRNA counts")
mrna_file <- unlist(strsplit(parser$mrna, ','))
mrna_df <- do.call("rbind", lapply(mrna_file, function(f) read_delim(f, delim = ' ', col_names = c('count', 'start_aln_UMI'), col_types = 'ic')))
mrna_df <- mrna_df %>% filter(start_aln_UMI %in% barcode_allele_filtered$start_aln_UMI) %>% group_by(start_aln_UMI) %>% summarise(count = sum(count))

message("Reading DNA file")                                   
dna_df <- read_delim(parser$dna, delim = ' ', col_names = c('count', 'start_aln_UMI'), col_types = 'ic')
dna_df <- dna_df %>% filter(start_aln_UMI %in% barcode_allele_filtered$start_aln_UMI) %>% filter(count > 0)
                                                                      
dna_mrna_df <- mrna_df %>% right_join(dna_df, by = 'start_aln_UMI', suffix = c('_mrna', '_dna'))%>% replace(is.na(.), 0)

result_df= data.frame(event = character(), 
                   ref_count = integer(), 
                   alt_count = integer(),
                   ref_expr = integer(),
                   alt_expr = integer(),
                   alt_effect = double(), 
                   alt_freq = double(), 
                   pvalue = double(),
                   z = double())

message("Running NBR")
for (event in var_list){
    # get the ref/alt barcode 
    ref_bar = barcode_allele %>% filter(mutation == event, type == 1) %>% pull(start_aln_UMI)
    alt_bar = barcode_allele %>% filter(mutation == event, type == 2) %>% pull(start_aln_UMI)

    # get the barcode plasmid abundance (SNP_CL)
    # get the barcode expression level at the experiment timepoint (e.g. SNPR_1_24, SNPR_2_24, SNPR_3_25)
    ref_df = dna_mrna_df %>% filter(start_aln_UMI %in% ref_bar) 
    alt_df = dna_mrna_df %>% filter(start_aln_UMI %in% alt_bar) 
    # number of barcode supporting each type of the SNP
    ref_count = nrow(ref_df)
    alt_count = nrow(alt_df)

    ref_df$type = 'ref'
    alt_df$type = 'alt'

    # merge expression data together
    analysis_df =rbind(ref_df, alt_df)

    # compute alternative allele frequency
    ref_norm <- ref_df %>% mutate(norm = count_mrna/count_dna) %>% pull(norm) %>% mean()
    alt_norm <- alt_df%>% mutate(norm = count_mrna/count_dna) %>% pull(norm) %>% mean()
    af = alt_norm/(alt_norm + ref_norm)

    # run nb regression
    if (ref_count*alt_count>0 & any(analysis_df$count_dna > 0)){
        tryCatch(
            {
                fit <- glm.nb(count_mrna ~ type + offset(log(count_dna)), data = analysis_df)
                coef_table = summary(fit)$coefficients
                effect = -coef_table['typeref', 'Estimate']
                p = coef_table['typeref', 'Pr(>|z|)']
                z = -coef_table['typeref', 'z value']
            }, 
            warning = function(w) {
                effect = p = z = NA
            }
        )
    }else{
        effect = p = z = NA
    }

    # track results
    result_row = list(event = event, 
                      ref_count = ref_count, 
                      alt_count = alt_count, 
                      ref_expr = sum(ref_df$count_mrna),
                      alt_expr = sum(alt_df$count_mrna),
                      alt_effect = effect, 
                      alt_freq = af,
                      pvalue = p, z = z )
    result_df <- result_df %>% bind_rows(result_row) 

}

result_df <- result_df %>% filter(ref_count >= min_count, alt_count >= min_count)
message("NBR done!")
### prep for annotation.
### snp index is created chromosomeid(withoutchr):position 1-based.
temp=str_split(string=result_df$event,pattern=";",n=4,simplify = T)[,1:2]
result_df$snpindex = str_c(temp[,1] %>% gsub("chr","",.),temp[,2],sep=":")
if (!is.null(parser$annotation) | parser$annotation != "noannot"){
    message("Running Annotation")
    f <- function(x,pos) subset(x, snpindex %in% result_df$snpindex)
    snp_info=readr::read_delim_chunked(file=parser$annotation, 
        callback=DataFrameCallback$new(f),chunk_size = 10000,
        col_names=c("snpindex","SNP"),
        delim="\t",
        progress=show_progress(),
        col_types= list(snpindex=col_character(),
        SNP=col_character())
        )
    # read snp information ###
    #snp_info = read_tsv(parser$annotation, col_name = c('SNP', 'EVENT', 'G5', 'ORIGINAL'), col_type = 'cccc')
    #result_df %>% left_join(snp_info, by = c('event'='EVENT'))%>% relocate(SNP))

    right_join(result_df, snp_info,by="snpindex") %>% dplyr::select(-snpindex) %>% dplyr::relocate(SNP) -> result_df
}
### multiple-hypothesis correction
if (nrow(result_df) ==0){
    print("No NBR found. Reduce the min count threshold or check input files!")
}else{
result_df$padj=p.adjust(result_df$pvalue,n=nrow(result_df),method="fdr")
result_df %>% write_tsv(.,path = parser$output)
}

message("Analysis Done!")