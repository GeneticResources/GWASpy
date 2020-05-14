
# function to fetch halopreg and regu data and gwas catlog

library(data.table)
library(haploR)
library(stringr)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)
print(args)

file_path <- args[1]
col_SNP <- args[2]
col_P <- args[3]

trait <- args[4]

dir <- args[5]


cat("*******************************************************************\n")

print(paste0("file_path: ", file_path))

print(paste0("col_SNP: ", col_SNP))

print(paste0("col_P: ", col_P))

print(paste0("trait: ", trait))
print(paste0("dir: ", dir))


cat("*******************************************************************\n")



# file_path<-'/working/lab_stuartma/xikunH/UKB/OCT/mtag/mtag_VCDR_adjustDD_VCDR_IGGC_2trait/trait_summary/cojo_summary_table_mtag_VCDR_adjustDD_VCDR_IGGC_2trait_20180322.csv'
# file_path<-'/working/lab_stuartma/xikunH/UKB/OCT/mtag/mtag_ukbPOAG_VCDR_adjustDD_IOP_4trait/trait_summary/cojo_summary_table_mtag_ukbPOAG_VCDR_adjustDD_IOP_4trait_20180322.csv'
# col_SNP<-'SNP' col_P<-'p'
# trait<-'mtag_ukbPOAG_VCDR_adjustDD_IOP_4trait'
# dir<-'/working/lab_stuartma/xikunH/UKB/OCT/mtag/'


setwd(dir)

df <- fread(file_path)

# filter sig P values
df <- df[get(col_P) < 5e-08, ]




SNP_list <- df[, get(col_SNP)]
SNP_list

# haploreg result set_config(use_proxy(url='127.0.0.1',port=0000))

df_hap <- queryHaploreg(query = SNP_list, verbose = TRUE, ldThresh = 0.8, 
    timeout = 1e+07)
df_hap <- data.table(df_hap)
df_hap

df_hap[] <- lapply(df_hap, as.character)

df_hap[, `:=`(r2, as.numeric(r2))]
# df_hap<-df_hap[order(SNP,-r2)]


# df_hap<-df_hap %>% group_by(query_snp_rsid) %>% arrange(desc(r2),chr,
# pos_hg38) %>% ungroup()

df_hap <- df_hap %>% arrange(chr, query_snp_rsid, desc(r2))

df_hap <- data.table(df_hap)


# dd<-df_hap[query_snp_rsid==rsID,]
dd <- copy(df_hap)

dest_dir <- "/working/lab_stuartma/xikunH/reference/gwas_catlog/20180523"


source("/working/lab_stuartma/xikunH/UKB/OCT/script/f_gwas_catlog.R")
df_gwas_catlog <- f_read_gwas_catlog(dest_dir = dest_dir, file = "catalog-associations_ontology-annotated.tsv", 
    verbose = TRUE)
names(df_gwas_catlog)

df_gwas_catlog[, `:=`(GWAS_report, paste0(PUBMEDID, "_", MAPPED_TRAIT, "_", 
    MAPPED_GENE, "_", `P-VALUE`))]

# df_gwas_catlog<-df_gwas_catlog[,c('DATE ADDED TO CATALOG','FIRST
# AUTHOR','DATE','JOURNAL','STUDY','DISEASE/TRAIT','MAPPED_TRAIT','CHR_ID','CHR_POS','REPORTED
# GENE(S)','MAPPED_GENE','STRONGEST SNP-RISK
# ALLELE','SNPS','CONTEXT','RISK ALLELE FREQUENCY','P-VALUE')]


df_gwas_merge <- merge(dd, df_gwas_catlog, by.x = "rsID", by.y = "SNPS", 
    all.x = TRUE)


df_gwas_merge[, `:=`(GWAS_report, paste(GWAS_report, collapse = "; \n")), 
    by = query_snp_rsid]

df_gwas_merge[, `:=`(GWAS_report, str_remove_all(GWAS_report, "NA; |\n"))]

names(df_gwas_merge)

df_gwas_merge[, length(unique(query_snp_rsid))]






# fwrite(dd,paste0('mtag/',trait,'/COJO_haploreg_',trait,'.csv'))



# Regulome result

df_Regu <- queryRegulome(query = SNP_list, verbose = TRUE, timeout = 1e+05)

df_Regu <- df_Regu$res.table
df_Regu[] <- lapply(df_Regu, as.character)
df_Regu <- as.data.table(df_Regu)
df_Regu[, `:=`(hits_Regulome = hits, score_Regulome = score)]
df_Regu <- df_Regu[, list(rsid, hits_Regulome, score_Regulome)]


df_gwas_merge <- merge(df_gwas_merge, df_Regu, by.x = "rsID", by.y = "rsid", 
    all = TRUE)

df_gwas_merge <- merge(df, df_gwas_merge, by.x = col_SNP, by.y = "rsID", 
    all = TRUE)

names(df_gwas_merge)

df_gwas_merge[, length(unique(query_snp_rsid))]

dd <- df_gwas_merge[!is.na(Chr), ]


dd <- dd[SNP == query_snp_rsid]
dd <- dd[(!duplicated(query_snp_rsid))]



fwrite(dd, paste0(dir, "/", trait, "_annotation_with_haploreg_gwas_catlog.csv"))

fwrite(df_gwas_merge, paste0(dir, "/", trait, "_annotation_with_haploreg_gwas_catlog_and_r2_0.8_SNPs.csv"))

