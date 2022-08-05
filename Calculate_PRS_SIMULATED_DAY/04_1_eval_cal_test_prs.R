#' Compute PRS using the selected SNPs and weights in the testing set
#' 
# for data processing
library(vroom)
library(glue)
library(yaml)
library(tidyverse)

# for lasso regression
library(biglasso)
library(bigsnpr)

# obtain input parameters
param_file <- "./config_files/demo_param_own.yml"
params <- read_yaml(param_file)[["all_params"]]

# output dir
eval_outdir  <- file.path(params$output_rt,params$eval_outdir)
system(command = glue("mkdir -p {eval_outdir}"))

# input data
lasso_outdir  <- file.path(params$output_rt,params$lasso_outdir)
lasso_plink_snp_file  <- file.path(lasso_outdir,"lasso_plink_inputs.tsv")

prsice2_outdir  <- file.path(params$output_rt,params$prsice2_out)
prsice2_plink_snp_file  <- file.path(prsice2_outdir,"plink_prsice2_snps.tsv")

test_geno  <- params$test_data

# calculate prsice2 prs in test data using plink2
test_prsice2_out  <- file.path(eval_outdir,'prsice2_test_prs')
plink_prsice2_cmd  <- paste0(
	glue("{params$plink2} --bfile {test_geno} "), # <PATH_TEST_STUFF>/TEST_EUR
	glue("--memory 15000 "),
	glue("--score {prsice2_plink_snp_file} cols=+scoresums list-variants --out {test_prsice2_out}")
)
system(command = plink_prsice2_cmd,wait=T)


# calculate lasso prs in test data using plink2
test_lasso_out  <- file.path(eval_outdir,'lasso_test_prs')
plink_lasso_cmd  <- paste0(
	glue("{params$plink2} --bfile {test_geno} "),
	glue("--memory 15000 "),
	glue("--score {lasso_plink_snp_file} cols=+scoresums list-variants --out {test_lasso_out}")
)
system(command = plink_lasso_cmd,wait=T)


