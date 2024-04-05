#Load libraries and constants
#install.packages("dplyr", repos = "http://cran.us.r-project.org")
#install.packages("parallel", repos = "http://cran.us.r-project.org")

options(scipen=1000, stringsAsFactors=F)
library(readr)
library(dplyr)
library(bedr)
library(PedSV)
library(argparse)
library(parallel)
PedSV::load.constants("all")

install.packages("tibble", version = "3.2.1", repos = "http://cran.us.r-project.org")
install.packages("stringr", repos = "http://cran.us.r-project.org")
install.packages("sqldf", version = "0.4-11", repos = "http://cran.us.r-project.org")

#Parse command line arguments and options
parser <- ArgumentParser(description="assign SVs to CWAS categories")

parser$add_argument("--SV_BED", metavar=".tsv", required=TRUE,
                    help=" ")

parser$add_argument("--protein_coding_gene_list", metavar=".txt", required=TRUE,
                    help=" ")

parser$add_argument("--cosmic_and_germline_CPGs_gene_list", metavar=".txt", required=TRUE,
                    help=" ")

parser$add_argument("--lof_constrained_gene_list", metavar=".txt", required=TRUE,
                    help=" ")

parser$add_argument("--expressed_in_adrenal_gland_gene_list", metavar=".txt", required=TRUE,
                    help=" ")

parser$add_argument("--expressed_in_muscle_skeletal_gene_list", metavar=".txt", required=TRUE,
                    help=" ")

parser$add_argument("--unique_cwas_categories_listed", metavar=".txt", required=TRUE,
                    help=" ")

args <- parser$parse_args()

#Load in bed file

BED_file_col_names <- read_tsv(args$SV_BED) %>% colnames()

BED <- read.table(args$SV_BED, col.names = BED_file_col_names) %>% filter(FILTER == "PASS") %>% rename(chrom = "X.chrom")

#Define gene sets (match naming of CWAS genic categories)

protein_coding <- read_tsv(args$protein_coding_gene_list) %>% dplyr::select(value) %>% unlist() %>% unname()

lof_constrained <- read_tsv(args$lof_constrained_gene_list) %>% dplyr::select(value) %>% unlist() %>% unname()

expressed_in_adrenal_gland <- read_tsv(args$expressed_in_adrenal_gland_gene_list) %>% dplyr::select(value) %>% unlist() %>% unname()

expressed_in_muscle_skeletal <- read_tsv(args$expressed_in_muscle_skeletal_gene_list) %>% dplyr::select(value) %>% unlist() %>% unname()

cosmic_and_germline_CPGs <- read_tsv(args$cosmic_and_germline_CPGs_gene_list) %>% dplyr::select(value) %>% unlist() %>% unname()

#Unique CWAS categories listed

noncoding_CWAS_enumerated_categories_unique_with_row_names <- read_tsv(args$unique_cwas_categories_listed) %>% tibble::column_to_rownames("unique_category_id")

#Filter SV callset

svtypes_for_inclusion <- c("DUP","DEL","INS","CPX","INS:ME","INS:ME:ALU","INS:ME:LINE1","INS:ME:SVA","INV","CTX")

BED_with_call_rate <- BED %>% mutate(call_rate = N_BI_GENOS/max(BED$N_BI_GENOS, na.rm = T))

BED_filtered_for_noncoding_CWAS <- BED_with_call_rate %>% filter(svtype %in% svtypes_for_inclusion & is.na(PREDICTED_BREAKEND_EXONIC) & is.na(PREDICTED_COPY_GAIN) & is.na(PREDICTED_DUP_PARTIAL) & is.na(PREDICTED_INTRAGENIC_EXON_DUP) & is.na(PREDICTED_LOF) & is.na(PREDICTED_MSV_EXON_OVERLAP) & is.na(PREDICTED_PARTIAL_EXON_DUP) & is.na(PREDICTED_TSS_DUP) & FILTER == "PASS" & !(chrom %in% c("chrX","chrY")))

BED_filtered_for_noncoding_CWAS_duplicated_columns_removed <- BED_filtered_for_noncoding_CWAS %>% dplyr::select(-c(END,SVTYPE))

BED_filtered_for_noncoding_CWAS_coordinate_lookup_table <- BED_filtered_for_noncoding_CWAS_duplicated_columns_removed %>% dplyr::select(name,`chrom`,start,end) %>% tibble::column_to_rownames("name")

#Define function to return structural variant names by

return_structural_variant_by_category_function <- function(category_name){
  
  sv_type <- noncoding_CWAS_enumerated_categories_unique_with_row_names[category_name,"sv_type"]
  frequency <- noncoding_CWAS_enumerated_categories_unique_with_row_names[category_name,"frequency"]
  functional_intersection <- noncoding_CWAS_enumerated_categories_unique_with_row_names[category_name,"functional_intersection"]
  functional_category <- noncoding_CWAS_enumerated_categories_unique_with_row_names[category_name,"functional_category"]
  genic_relationship <- noncoding_CWAS_enumerated_categories_unique_with_row_names[category_name,"genic_relationship"]
  gene_group <- noncoding_CWAS_enumerated_categories_unique_with_row_names[category_name,"gene_group"]
  expression <- noncoding_CWAS_enumerated_categories_unique_with_row_names[category_name,"expression"]
  constraint <- noncoding_CWAS_enumerated_categories_unique_with_row_names[category_name,"constraint"]
  
  gene_group_evaluated <- eval(as.symbol(gene_group))
  expression_evaluated <- if(!expression == "ANY"){eval(as.symbol(expression))}else{eval(as.symbol(gene_group))}
  constraint_evaluated <- if(!constraint == "ANY"){eval(as.symbol(constraint))}else{eval(as.symbol(gene_group))}
  gene_group_evaluated_filtered_for_expression <- Reduce(intersect, list(gene_group_evaluated,expression_evaluated,constraint_evaluated))
  
  if(length(gene_group_evaluated_filtered_for_expression) == 0){
    
    list_of_SVs_meeting_sv_type_frequency_functional_genic_criteria <- c()
    
  }
  
  else{
  
  gene_group_only_collapsed <- stringr::str_c(stringr::str_c("^",gene_group_evaluated_filtered_for_expression,"$"), collapse = "|")
  gene_group_start_collapsed <- stringr::str_c(stringr::str_c("^",gene_group_evaluated_filtered_for_expression,"(?=,)"), collapse = "|")
  gene_group_middle_collapsed <- stringr::str_c(stringr::str_c("(?<=,)",gene_group_evaluated_filtered_for_expression,"(?=,)"), collapse = "|")
  gene_group_end_collapsed <- stringr::str_c(stringr::str_c("(?<=,)",gene_group_evaluated_filtered_for_expression,"$"), collapse = "|")
  
  #SVs meeting sv_type criteria
  
  if(sv_type == "ANY"){
    
    SVs_meeting_sv_type_criteria <- BED_filtered_for_noncoding_CWAS_duplicated_columns_removed
    
  }
  else if(sv_type == "INS_ALL"){
    
    SVs_meeting_sv_type_criteria <- BED_filtered_for_noncoding_CWAS_duplicated_columns_removed[BED_filtered_for_noncoding_CWAS_duplicated_columns_removed$svtype %in% c("INS","INS:ME:SVA","INS:ME:ALU","INS:ME","INS:ME:LINE1"),]
    
  }
  else if(sv_type == "CPX_or_INV"){
    
    SVs_meeting_sv_type_criteria <- BED_filtered_for_noncoding_CWAS_duplicated_columns_removed[BED_filtered_for_noncoding_CWAS_duplicated_columns_removed$svtype %in% c("CPX","INV"),]
    
  }
  else{
    
    SVs_meeting_sv_type_criteria <- BED_filtered_for_noncoding_CWAS_duplicated_columns_removed[BED_filtered_for_noncoding_CWAS_duplicated_columns_removed$svtype == sv_type,]
  }
  
  
  #SVs meeting frequency criteria
  
  if(frequency == "RARE"){
    
    SVs_meeting_sv_type_frequency_criteria <- SVs_meeting_sv_type_criteria[SVs_meeting_sv_type_criteria$POPMAX_AF < .01 & (SVs_meeting_sv_type_criteria$gnomad_v3.1_sv_POPMAX_AF < .01 | is.na(SVs_meeting_sv_type_criteria$gnomad_v3.1_sv_POPMAX_AF)),]
    
  }
  else if(frequency == "SINGLETON"){
    
    SVs_meeting_sv_type_frequency_criteria <- SVs_meeting_sv_type_criteria[SVs_meeting_sv_type_criteria$AC == 1 & (SVs_meeting_sv_type_criteria$gnomad_v3.1_sv_POPMAX_AF < .01 | is.na(SVs_meeting_sv_type_criteria$gnomad_v3.1_sv_POPMAX_AF)),]
    
  }
  
  
  #SVs meeting functional criteria
  
  if(functional_intersection == "PREDICTED_NONCODING_BREAKPOINT"){
    
    SVs_meeting_sv_type_frequency_functional_criteria <- SVs_meeting_sv_type_frequency_criteria[!is.na(SVs_meeting_sv_type_frequency_criteria$PREDICTED_NONCODING_BREAKPOINT) & stringr::str_detect(SVs_meeting_sv_type_frequency_criteria$PREDICTED_NONCODING_BREAKPOINT, paste0("^",functional_category,"$|^",functional_category,"(?=,)|(?<=,)",functional_category,"(?=,)|(?<=,)",functional_category,"$")),]
    
  }
  else if(functional_intersection == "PREDICTED_NONCODING_SPAN"){
    
    SVs_meeting_sv_type_frequency_functional_criteria <- SVs_meeting_sv_type_frequency_criteria[!is.na(SVs_meeting_sv_type_frequency_criteria$PREDICTED_NONCODING_SPAN) & stringr::str_detect(SVs_meeting_sv_type_frequency_criteria$PREDICTED_NONCODING_SPAN, paste0("^",functional_category,"$|^",functional_category,"(?=,)|(?<=,)",functional_category,"(?=,)|(?<=,)",functional_category,"$")),]
    
  }
  else if(functional_intersection == "ANY"){
    
    SVs_meeting_sv_type_frequency_functional_criteria <- SVs_meeting_sv_type_frequency_criteria[(!is.na(SVs_meeting_sv_type_frequency_criteria$PREDICTED_NONCODING_BREAKPOINT) & stringr::str_detect(SVs_meeting_sv_type_frequency_criteria$PREDICTED_NONCODING_BREAKPOINT, paste0("^",functional_category,"$|^",functional_category,"(?=,)|(?<=,)",functional_category,"(?=,)|(?<=,)",functional_category,"$"))) | (!is.na(SVs_meeting_sv_type_frequency_criteria$PREDICTED_NONCODING_SPAN) & stringr::str_detect(SVs_meeting_sv_type_frequency_criteria$PREDICTED_NONCODING_SPAN, paste0("^",functional_category,"$|^",functional_category,"(?=,)|(?<=,)",functional_category,"(?=,)|(?<=,)",functional_category,"$"))),]
    
  }
  
  
  #SVs meeting genic criteria
  
  if(genic_relationship == "PREDICTED_INTERGENIC"){
    
    SVs_meeting_sv_type_frequency_functional_genic_criteria <- SVs_meeting_sv_type_frequency_functional_criteria[(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_INTERGENIC == "True" & !is.na(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_NEAREST_TSS) & stringr::str_detect(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_NEAREST_TSS,paste0(gene_group_only_collapsed,"|",gene_group_start_collapsed,"|",gene_group_middle_collapsed,"|",gene_group_end_collapsed))),]
    
  }
  else if(genic_relationship == "PREDICTED_INTRONIC"){
    
    SVs_meeting_sv_type_frequency_functional_genic_criteria <- SVs_meeting_sv_type_frequency_functional_criteria[!is.na(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_INTRONIC) & stringr::str_detect(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_INTRONIC,paste0(gene_group_only_collapsed,"|",gene_group_start_collapsed,"|",gene_group_middle_collapsed,"|",gene_group_end_collapsed)),]
    
  }
  else if(genic_relationship == "PREDICTED_PROMOTER"){
    
    SVs_meeting_sv_type_frequency_functional_genic_criteria <- SVs_meeting_sv_type_frequency_functional_criteria[!is.na(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_PROMOTER) & stringr::str_detect(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_PROMOTER,paste0(gene_group_only_collapsed,"|",gene_group_start_collapsed,"|",gene_group_middle_collapsed,"|",gene_group_end_collapsed)),]
    
  }
  else if(genic_relationship == "PREDICTED_UTR"){
    
    SVs_meeting_sv_type_frequency_functional_genic_criteria <- SVs_meeting_sv_type_frequency_functional_criteria[!is.na(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_UTR) & stringr::str_detect(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_UTR,paste0(gene_group_only_collapsed,"|",gene_group_start_collapsed,"|",gene_group_middle_collapsed,"|",gene_group_end_collapsed)),]
    
  }
  else if(genic_relationship == "ANY"){
    
    SVs_meeting_sv_type_frequency_functional_genic_criteria <- SVs_meeting_sv_type_frequency_functional_criteria[(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_INTERGENIC == "True" & !is.na(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_NEAREST_TSS) & stringr::str_detect(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_NEAREST_TSS,paste0(gene_group_only_collapsed,"|",gene_group_start_collapsed,"|",gene_group_middle_collapsed,"|",gene_group_end_collapsed))) | (!is.na(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_INTRONIC) & stringr::str_detect(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_INTRONIC,paste0(gene_group_only_collapsed,"|",gene_group_start_collapsed,"|",gene_group_middle_collapsed,"|",gene_group_end_collapsed))) | (!is.na(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_PROMOTER) & stringr::str_detect(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_PROMOTER,paste0(gene_group_only_collapsed,"|",gene_group_start_collapsed,"|",gene_group_middle_collapsed,"|",gene_group_end_collapsed))) | (!is.na(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_UTR) & stringr::str_detect(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_UTR,paste0(gene_group_only_collapsed,"|",gene_group_start_collapsed,"|",gene_group_middle_collapsed,"|",gene_group_end_collapsed))) | (!is.na(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_INV_SPAN) & stringr::str_detect(SVs_meeting_sv_type_frequency_functional_criteria$PREDICTED_INV_SPAN,paste0(gene_group_only_collapsed,"|",gene_group_start_collapsed,"|",gene_group_middle_collapsed,"|",gene_group_end_collapsed))),]
    
  }
  
  list_of_SVs_meeting_sv_type_frequency_functional_genic_criteria <- SVs_meeting_sv_type_frequency_functional_genic_criteria$name  
  
  }
  
  #Return
  
  return(BED_filtered_for_noncoding_CWAS_coordinate_lookup_table[list_of_SVs_meeting_sv_type_frequency_functional_genic_criteria,] %>% tibble::rownames_to_column("SV") %>% mutate(category = category_name))
  
}

#Identify SVs falling into each category

SVs_in_each_category_list <- mclapply(row.names(noncoding_CWAS_enumerated_categories_unique_with_row_names), return_structural_variant_by_category_function)

SVs_in_each_category_list_combined <- bind_rows(SVs_in_each_category_list)

write_tsv(SVs_in_each_category_list_combined,"SVs_in_each_category_list_combined.txt",col_names = FALSE)
