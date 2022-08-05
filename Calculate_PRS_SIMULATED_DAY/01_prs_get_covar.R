#' Obtain some covariates to run PRSICE2 and Biglasso just as an example
#' Will obtain
#'  1. 6 PC using plink
#'  2. sex extracted from Fam files
# for data processing
library(vroom)
library(glue)
library(yaml)
library(tidyverse)


# obtain input parameters
param_file <- "./config_files/demo_param_own.yml"
params <- read_yaml(param_file)[["all_params"]]


outdir  <- file.path(params$output_rt,params$covar_output)
system(command = glue("mkdir -p {outdir}"))

genotype_file  <- params$train_data

# get covariate

fam_file  <- vroom(file  = glue("{genotype_file}.fam"),
                  col_names=c("FID","IID","MIDaa","MID","sex","pheno"),delim=" ")

meta_information_file  <- fam_file %>% select(FID,IID,sex)

# obtain principle components
plink_prune_cmd  <- paste0(
  glue("{params$plink2} ",
  glue("--bfile {genotype_file} "),
  glue("--indep-pairwise 200 50 0.25 "),
  glue("--out {file.path(outdir,'plink_eur_prune')}"))
)
system(command = plink_prune_cmd,wait=T)

plink_pc_cmd  <- paste0(
  glue("{params$plink2} "),
  glue("--bfile {genotype_file} "),
  glue("--extract {file.path(outdir,'plink_eur_prune.prune.in')} "),
  glue("--pca 6 "),
  glue("--out {file.path(outdir,'plink_train_eur_pc')}")
)

system(command = plink_pc_cmd,wait=T)


pc_data  <- vroom(glue("{file.path(outdir,'plink_train_eur_pc')}.eigenvec"),delim="\t",
  col_names=T) %>%
  dplyr::rename(FID = `#FID`)


full_covar_data  <- inner_join(x = pc_data, y = meta_information_file,by=c("FID","IID"))

meta_output  <- file.path(outdir,'final_meta_information.tsv')
vroom_write(x=full_covar_data,file=meta_output,delim="\t",col_names=T,quote='none')

