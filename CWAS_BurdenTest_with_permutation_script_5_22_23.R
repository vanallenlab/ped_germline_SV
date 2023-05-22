
#Load libraries and constants
#install.packages("dplyr", repos = "http://cran.us.r-project.org")

options(scipen=1000, stringsAsFactors=F)
library(readr)
library(dplyr)
library(bedr)
library(PedSV)
library(argparse)
PedSV::load.constants("all")

#Parse command line arguments and options
parser <- ArgumentParser(description="Plot SV counts per sample")
parser$add_argument("--ad_matrix", metavar=".tsv", required=TRUE,
                    help=" ")
parser$add_argument("--gatk_sv_pediatric_cancers_combined_cohort_metadata", metavar=".txt", required=TRUE,
                    help=" ")
parser$add_argument("--SVs_in_each_category_list_combined", metavar=".tsv", required=TRUE,
                    help=" ")
parser$add_argument("--case_control_BED", metavar=".tsv", required=TRUE,
                    help=" ")
parser$add_argument("--disease", metavar="string", required=TRUE,
                    help=" ")
parser$add_argument("--study_phase", metavar="string", required=TRUE,
                    help=" ")
parser$add_argument("--category_list", metavar=".tsv", required=TRUE,
                    help=" ")
parser$add_argument("--number_of_permutations", metavar="string", required=TRUE,
                    help=" ")
args <- parser$parse_args()

#Load arguments and variables

ad_matrix <- args$ad_matrix
gatk_sv_pediatric_cancers_combined_cohort_metadata <- load.sample.metadata(args$gatk_sv_pediatric_cancers_combined_cohort_metadata)
SVs_in_each_category_list_combined <- read_tsv(args$SVs_in_each_category_list_combined)
case_control_BED <- load.sv.bed(args$case_control_BED)
disease <- args$disease
study_ph <- args$study_phase
category_list <- read_tsv(args$category_list, col_names = FALSE) %>% unlist() %>% unname()
number_of_permutations <- args$number_of_permutations %>% as.integer()

gatk_sv_pediatric_cancers_combined_cohort_metadata_filtered_for_study_phase <- gatk_sv_pediatric_cancers_combined_cohort_metadata %>% dplyr::filter(study_phase == study_ph)

count <- 0
glm_results <- list()
glm_results_permuted_z_score_matrix <- matrix(nrow = length(category_list), ncol = number_of_permutations, dimnames = list(category_list,1:number_of_permutations))

for(category_name in category_list){

SVs_of_interest <- SVs_in_each_category_list_combined[SVs_in_each_category_list_combined$category == category_name,"SV"] %>% unlist() %>% unname()

subsetted_case_control_BED <- case_control_BED[SVs_of_interest,]

SV_count_vector <- query.ad.from.sv.bed(ad_matrix, subsetted_case_control_BED, action="sum")

sample_ids <- get.eligible.samples(gatk_sv_pediatric_cancers_combined_cohort_metadata_filtered_for_study_phase,disease)
  
case_control_vector <- get.phenotype.vector(sample_ids$cases,sample_ids$controls)

SV_count_vector_subsetted <- SV_count_vector[names(case_control_vector)][!is.na(SV_count_vector[names(case_control_vector)])]

if(SV_count_vector_subsetted %>% length() > 0 & SV_count_vector_subsetted[SV_count_vector_subsetted==1] %>% length() > 0){
  
  SV_counts_cases <- SV_count_vector_subsetted[sample_ids$cases][!is.na(SV_count_vector_subsetted[sample_ids$cases])] %>% sum()
  
  SV_counts_cases_max <- SV_count_vector_subsetted[sample_ids$cases][!is.na(SV_count_vector_subsetted[sample_ids$cases])] %>% max()
  
  number_of_cases_with_zero_SVs <- length(SV_count_vector_subsetted[sample_ids$cases][!is.na(SV_count_vector_subsetted[sample_ids$cases]) & SV_count_vector_subsetted[sample_ids$cases] == 0])
  
  total_cases <- SV_count_vector_subsetted[sample_ids$cases][!is.na(SV_count_vector_subsetted[sample_ids$cases])] %>% length()
  
  SV_counts_controls <- SV_count_vector_subsetted[sample_ids$controls][!is.na(SV_count_vector_subsetted[sample_ids$controls])] %>% sum()
  
  SV_counts_controls_max <- SV_count_vector_subsetted[sample_ids$controls][!is.na(SV_count_vector_subsetted[sample_ids$controls])] %>% max()
  
  number_of_controls_with_zero_SVs <- length(SV_count_vector_subsetted[sample_ids$controls][!is.na(SV_count_vector_subsetted[sample_ids$controls]) & SV_count_vector_subsetted[sample_ids$controls] == 0])
  
  total_controls <- SV_count_vector_subsetted[sample_ids$controls][!is.na(SV_count_vector_subsetted[sample_ids$controls])] %>% length()
  
  number_of_unique_SVs <- length(SVs_of_interest)
  
  glm_results[[category_name]] <- c(pedsv.glm(gatk_sv_pediatric_cancers_combined_cohort_metadata_filtered_for_study_phase,SV_count_vector_subsetted,case_control_vector,family = binomial(), use.N.pcs = 5),SV_counts_cases = SV_counts_cases, SV_counts_cases_max = SV_counts_cases_max, number_of_cases_with_zero_SVs = number_of_cases_with_zero_SVs, total_cases = total_cases, SV_counts_controls = SV_counts_controls, SV_counts_controls_max = SV_counts_controls_max, number_of_controls_with_zero_SVs = number_of_controls_with_zero_SVs, total_controls = total_controls, number_of_unique_SVs = number_of_unique_SVs, category = category_name)
  
  #Permutated glm_results
  
  for(iteration in 1:number_of_permutations){
  
  set.seed(iteration)
  
  update_names <- names(SV_count_vector_subsetted) %>% sample()
  
  names(SV_count_vector_subsetted) <- update_names
  
  glm_results_permuted_z_score_matrix[[category_name,iteration]] <- pedsv.glm(gatk_sv_pediatric_cancers_combined_cohort_metadata_filtered_for_study_phase,SV_count_vector_subsetted,case_control_vector,family = binomial(), use.N.pcs = 5)[[3]]
  
}

}

count <- count + 1
print(count)

}

write_tsv(bind_rows(glm_results),"glm_results.txt", col_names = FALSE)

write.table(glm_results_permuted_z_score_matrix,"glm_results_permuted_z_score_matrix.txt", row.names=TRUE, col.names = FALSE,sep = "\t")
