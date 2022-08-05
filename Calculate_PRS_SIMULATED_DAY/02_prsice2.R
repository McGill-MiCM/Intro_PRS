#' Will run PRSICE2.
# for data processing
library(vroom)
library(glue)
library(yaml)


# obtain input parameters
param_file <- "./config_files/demo_param_own.yml"
params <- read_yaml(param_file)[["all_params"]]

outdir  <- file.path(params$output_rt,params$prsice2_out)
system(command = glue("mkdir -p {outdir}"))

# prsice2 inputs
genotype_file  <- params$train_data
covar_file  <- file.path(params$output_rt,params$covar_output,'final_meta_information.tsv')
pheno  <- glue("{params$train_data}_HEIGHT.tsv")

prsice2_output  <- file.path(outdir,"prsice2_output")

prsice2_cmd  <- paste0(
  glue("{params$prsice2} "),
  glue("--base {params$sumstat} "),# Simulated_Height_gwas.fastGWA.gz
  glue("--target {genotype_file} "), # TRAIN_EUR
  glue("--binary-target F "),
  glue("--pheno {pheno} "), # TRAIN_EUR_HEIGHT.tsv
  glue("--pheno-col Height "),
  glue("--cov {covar_file} "), # final_meta_information.tsv
  glue("--cov-factor sex "),
  glue("--stat BETA "),
  glue("--a1 A1 --a2 A2 --beta --chr CHR --bp BP --pvalue P --snp SNP "),
  glue("--print-snp "),
  glue("--score sum "),
  glue("--out {prsice2_output}") # <NAME>
)

system(command = prsice2_cmd,wait=T)

sumstats  <- vroom(params$sumstat,col_names=T,delim="\t")
prsice_snp_lists  <- vroom(file.path(outdir,"prsice2_output.snp"),col_names=T,delim="\t")
prs_snps  <- sumstats %>% filter(SNP %in% prsice_snp_lists$SNP) %>%
	select(SNP,A1,BETA)

vroom_write(x = prs_snps,file=file.path(outdir,"plink_prsice2_snps.tsv"),
	col_names=F, quote='none',delim="\t")
