
# created 19 April 2018

# use to check mtag beta, se and z

# usage:
# f_check_mtag_beta(file_mtag='mtag/mtag_trait_fdr/out_mtag_trait_fdr_trait_1.txt',
# file_raw='out_gwas/glaucoma/glaucoma_GWAS_result_all_20180322',
# snp_mtag='snpid',snp_raw='SNP',a1_raw='ALT',a2_raw='REF',
# beta_raw='OR',se_raw='SE',OR_raw=TRUE)



# file_mtag='mtag/mtag_trait_fdr/out_mtag_trait_fdr_trait_1.txt';
# file_raw='out_gwas/glaucoma/glaucoma_GWAS_result_all_20180322'
# snp_mtag='snpid' snp_raw='SNP' a1_mtag='a1' a2_mtag='a2' a1_raw='ALT'
# a2_raw='REF' beta_raw='OR' se_raw='SE' z_raw='z' OR=TRUE scale_ratio;
# variance of phenotype

library(data.table)

f_check_mtag_beta <- function(file_mtag = NULL, file_raw = NULL, snp_mtag = "snpid", 
    snp_raw = NULL, a1_raw = NULL, a2_raw = NULL, beta_raw = NULL, se_raw = NULL, 
    OR_raw = FALSE, scale_ratio = 1) {
    
    
    cat(paste0("Read in files: \n"))
    df_mtag <- fread(file_mtag)
    print(df_mtag)
    df_raw <- fread(file_raw)
    print(df_raw)
    
    stopifnot(all(c(snp_raw, beta_raw, se_raw, a1_raw, a2_raw) %in% names(df_raw)))
    
    stopifnot(all(c(snp_mtag, "z", "a1", "a2", "mtag_beta", "mtag_se", "mtag_z") %in% 
        names(df_mtag)))
    
    setnames(df_mtag, c("a1", "a2", "z", "n", "freq"), paste0(c("a1", "a2", 
        "z", "n", "freq"), "_in_mtag"))
    
    
    df_raw[, `:=`(c(a1_raw, a2_raw), lapply(.SD, toupper)), .SDcols = c(a1_raw, 
        a2_raw)]
    
    df_merge <- merge(df_mtag, df_raw, by.x = snp_mtag, by.y = snp_raw, suffixes = c("_mtag", 
        "_raw"))
    
    
    cat(paste0("\nMerged file names: \n"))
    print(head(df_merge))
    cat(paste0("\nCheck the allele: \n"))
    print(table(df_merge$a1_in_mtag == unlist(df_merge[, a1_raw, with = FALSE]), 
        df_merge$a2_in_mtag == unlist(df_merge[, a2_raw, with = FALSE])))
    if (OR_raw) {
        df_merge[, `:=`(beta_gwas_convert, lapply(.SD, log)), .SDcols = beta_raw]
        cat(paste0("\nConvert OR to beta: \n"))
    } else {
        df_merge[, `:=`(beta_gwas_convert, .SD), .SDcols = beta_raw]
    }
    
    cat("\n")
    print(head(df_merge))
    print(dim(df_merge))
    df_merge[, `:=`(beta_gwas_convert, ifelse(df_merge$a1_in_mtag == unlist(df_merge[, 
        a1_raw, with = FALSE]), beta_gwas_convert, -beta_gwas_convert))]
    
    
    print(df_merge[, lm(beta_gwas_convert ~ I(mtag_beta * scale_ratio) - 
        1)])
    print(df_merge[, lm(get(se_raw) ~ I(mtag_se * scale_ratio) - 1)])
    print(df_merge[, lm(z_in_mtag ~ mtag_z - 1)])
    
    cat("regression in sig SNP.*********************************************\n")
    print(df_merge[mtag_pval < 1e-06, lm(beta_gwas_convert ~ I(mtag_beta * 
        scale_ratio) - 1)])
    print(df_merge[mtag_pval < 1e-06, lm(get(se_raw) ~ I(mtag_se * scale_ratio) - 
        1)])
    print(df_merge[mtag_pval < 1e-06, lm(z ~ mtag_z - 1)])
}
