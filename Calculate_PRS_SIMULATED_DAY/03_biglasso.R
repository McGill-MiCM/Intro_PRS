#' Will run biglasso
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

prsice_outdir  <- file.path(params$output_rt,params$prsice2_out)
lasso_outdir  <- file.path(params$output_rt,params$lasso_outdir)
system(command = glue("mkdir -p {lasso_outdir}"))

genotype_file  <- params$train_data
bigsnp_bk_file  <- file.path(lasso_outdir,"TRAIN_EUR_bk")
if (file.exists(glue("{bigsnp_bk_file}.rds"))){
	unlink(paste0(bigsnp_bk_file,c(".rds",".bk")))
	
}
snp_readBed(bedfile = glue("{genotype_file}.bed"),backingfile = bigsnp_bk_file)
bigsnp_obj  <- snp_attach(glue("{bigsnp_bk_file}.rds"))
genotype_mat  <- bigsnp_obj$genotypes
map  <- bigsnp_obj$map

prsice_snps  <- vroom(file.path(prsice_outdir,"plink_prsice2_snps.tsv"),
 	col_names=c("SNP","A1","BETA"),
	delim="\t")

folds  <- 2
# due to limited sample size, need to  reduce the predictors number so can do something

predictors  <- intersect(map$marker.ID,prsice_snps$SNP)

genotype_mat_col_idx  <- na.omit(match(predictors,map$marker.ID))
subsetted_geno_matrix  <- genotype_mat[,genotype_mat_col_idx]

# obtain response variable
response  <- vroom(glue("{genotype_file}_HEIGHT.tsv"),col_names=T,delim="\t")

# obtain sex and pc covariates (can choose what if any covariates to add)
covar_file  <- file.path(params$output_rt,params$covar_output,'final_meta_information.tsv')
covar_info  <- vroom(covar_file,delim="\t",col_names=T)


# order covariates the same order as design matrix (response)
covar_info  <- covar_info[match(response$FID,covar_info$FID),]

subsetted_geno_matrix  <- cbind(subsetted_geno_matrix,as.matrix(covar_info %>% select(-c(FID,IID))))
col_names  <- c(predictors,colnames(covar_info %>% select(-c(FID,IID))))

if (file.exists(file.path(lasso_outdir,'design_mat_bk'))){
	unlink(c(file.path(lasso_outdir,c('design_mat_bk',"design_mat_desc"))))
}
subsetted_geno_matrix_bm  <- as.big.matrix(
	x = scale(subsetted_geno_matrix),
	type = 'double',
	backingfile = 'design_mat_bk',
	descriptorfile = 'design_mat_desc',
	backingpath = lasso_outdir)
options(bigmemory.allow.dimnames=TRUE)
colnames(subsetted_geno_matrix_bm)  <- col_names

res  <- cv.biglasso(
	X = subsetted_geno_matrix_bm,
	y = as.numeric(response$Height),
	verbose = TRUE,nfolds=5,
	family='gaussian')

nonzero_predictors  <- as.data.frame(
	matrix(
		data = c(names(coef(res)[which(coef(res)!=0),]),
						coef(res)[which(coef(res)!=0),]),
		ncol = 2, dimnames=list(c(),c("predictor","beta"))))

print(nrow(nonzero_predictors))

sumstats  <- vroom(params$sumstat,col_names=T,delim="\t")
lasso_plink_snps  <- sumstats %>% filter(SNP %in% nonzero_predictors$predictor) %>% select(SNP,A1)
lasso_plink_snps <- inner_join(x = nonzero_predictors,y = lasso_plink_snps,by=c("predictor" = "SNP")) %>%
	select(predictor,A1,beta)
vroom_write(lasso_plink_snps,file = file.path(lasso_outdir,'lasso_plink_inputs.tsv'),col_names=F,delim="\t",quote='none')







