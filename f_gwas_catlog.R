

# function to download latest gwas catlog datasets
f_download_gwas_catlog <- function(dest_dir = NULL, var_date = format(Sys.time(), 
    "%Y%m%d"), allow_dir_exist = FALSE, download_all = FALSE) {
    if (is.null(dest_dir)) {
        stop("Please specify dest_dir!")
    } else if (dir.exists(dest_dir) && allow_dir_exist) {
        message(paste0("Download to dir: ", dest_dir))
    } else if (dir.exists(dest_dir) && !allow_dir_exist) {
        stop("dest_dir exist, please specify a new dir or allow download to exist dir")
    } else {
        dir.create(dest_dir)
        message(paste0("Download to dir: ", dest_dir))
    }
    
    # file_name<-RCurl::getURL(url,verbose=FALSE,ftp.use.epsv=FALSE,
    # dirlistonly = TRUE) file_name<-strsplit(file_name,'\n')
    # file_name<-file_name[[1]]
    
    url <- "ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
    file_info <- RCurl::getURL(url)
    file_info <- matrix(strsplit(file_info, " +|\\n")[[1]], ncol = 9, byrow = TRUE)
    file_info <- data.frame(file_info, stringsAsFactors = FALSE)
    cat("Lastest gwas catlog in ftp:\n")
    print(file_info)
    cat("\n")
    
    file_name <- file_info$X9
    if (download_all) {
        file_name <- file_name
    } else {
        file_to_download <- c("gwas-catalog-ancestry.tsv", "gwas-catalog-studies_ontology-annotated.tsv", 
            "gwas-catalog-associations_ontology-annotated.tsv", "gwas-catalog-studies.tsv", 
            "gwas_diagram.svg")
        
        if (all(file_to_download %in% file_name)) {
            file_name <- file_to_download
        } else {
            warning(paste0("File names were changed in ftp, download all files!"))
            file_name <- file_name
        }
        
    }
    cat("Download files:", "\n")
    cat(paste(c(file_name, collapse = ""), "\n"))
    
    for (var_file in file_name) {
        file_name_each <- var_file[[1]]
        out_file_each <- paste0(dest_dir, "/", var_date, "_", file_name_each)
        cat(paste0(url, file_name_each, "\n\n"))
        download.file(paste0(url, file_name_each), out_file_each)
    }
    
}




# to read in gwas catlog results
f_read_gwas_catlog <- function(dest_dir = NULL, file = NULL, verbose = FALSE) {
    
    url <- "ftp://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
    file_info <- RCurl::getURL(url)
    file_info <- matrix(strsplit(file_info, " +|\\n")[[1]], ncol = 9, byrow = TRUE)
    file_info <- data.frame(file_info, stringsAsFactors = FALSE)
    cat("Lastest gwas catlog in ftp:\n")
    print(file_info)
    cat("\n")
    
    
    file_list <- list.files(dest_dir)
    file_grep <- grep(file, file_list, value = TRUE)
    file_in <- paste0(dest_dir, "/", file_grep)
    
    if (length(file_grep) == 1 && file.exists(file_in)) {
        cat(paste0("Read in file: ", file_in, "\n"))
        df <- data.table::fread(file_in, sep = "\t", fill = TRUE)
        df[, `:=`(c("DATE ADDED TO CATALOG", "DATE"), list(as.Date(`DATE ADDED TO CATALOG`), 
            as.Date(DATE)))]
        
        df <- df[order(DATE, decreasing = TRUE)]
        
        if (verbose) {
            
            dd <- df[, list(`DATE ADDED TO CATALOG`, PUBMEDID, `FIRST AUTHOR`, 
                DATE, JOURNAL, `DISEASE/TRAIT`)]
            dd <- dd[!duplicated(dd)]
            
            dd <- dd[order(`DATE ADDED TO CATALOG`, decreasing = TRUE)]
            
            cat("********************* 20 Studies recently added to catalog. *******************************\n\n")
            print(head(dd, 20))
            
            dd <- dd[order(DATE, decreasing = TRUE)]
            
            cat("\n\n********************* 20 Studies recently published. *******************************\n")
            print(head(dd, 20))
        }
        
        
        
        return(df)
    } else {
        cat(paste(c("Read in file: ", file_grep, collapse = ""), "\n"))
        stop(paste0("Cannot read in file: ", file))
    }
}




f_queryHaploreg <- function(SNP_list = NULL, verbose = TRUE, ldThresh = 0.8, 
    timeout = 10000, save_file = NULL, ...) {
    
    cat(paste0("Start Haploreg:\n"))
    result <- vector("list", length(SNP_list))
    
    result <- lapply(SNP_list, function(snp) {
        tryCatch({
            haploR::queryHaploreg(query = snp, verbose = verbose, ldThresh = ldThresh, 
                timeout = timeout, ...)
        }, error = function(e) {
            message(paste0("Cannot get results for SNP: ", snp, " (not in 1000G?)\n"))
        })
    })
    result <- plyr::ldply(result, data.frame)
    result[] <- lapply(result, as.character)
    result$r2 <- as.numeric(result$r2)
    result$chr <- as.numeric(result$chr)
    result$pos_hg38 <- as.numeric(result$pos_hg38)
    result <- data.frame(query_snp_rsid = result$query_snp_rsid, result[, 
        names(result) != "query_snp_rsid"], stringsAsFactors = FALSE)
    
    result_Haploreg <- plyr::arrange(result, chr, query_snp_rsid, plyr::desc(is_query_snp), 
        plyr::desc(r2))
    
    cat("\n")
    cat(paste0("Start Regulome:\n"))
    
    
    result <- vector("list", length(SNP_list))
    result <- lapply(SNP_list, function(snp) {
        tryCatch({
            dd <- haploR::queryRegulome(query = snp, verbose = verbose, timeout = timeout, 
                check_bad_snps = TRUE)
            # print(NROW(dd[[1]])==0)
            if (NROW(dd[[1]]) == 0) 
                stop("failed!")
            dd[[1]]
        }, error = function(e) {
            message(paste0("Cannot get results for SNP: ", snp, " (not in 1000G?)\n"))
        })
    })
    
    result <- plyr::ldply(result, data.frame)
    result[] <- lapply(result, as.character)
    
    result$hits_Regulome <- result$hits
    result$score_Regulome <- result$score
    result_Regu <- result[, c("rsid", "hits_Regulome", "score_Regulome")]
    result <- merge(result_Haploreg, result_Regu, by.x = "rsID", by.y = "rsid", 
        all = TRUE)
    if (!is.null(save_file)) {
        fwrite(result, save_file)
    }
    return(result)
}



f_Haploreg_gaws_catlog <- function(SNP_list = NULL, gwas_catlog_dir = NULL, 
    gwas_catlog_file = NULL, ldThresh = 0.8) {
    
    # df_hap<-f_queryHaploreg(SNP_list = SNP_list,ldThresh=ldThresh)
    # df_gwascatlog<-f_read_gwas_catlog(dest_dir =
    # gwas_catlog_dir,file=gwas_catlog_file,verbose = TRUE)
    
    df_hap <- f_queryHaploreg(SNP_list = SNP_list)
    
    # fwrite(df_hap,'out_gwas/DISC_rank_chunk_mean/query_from_haploreg.csv',na='NA',quote
    # = FALSE)
    
    df_gwascatlog <- f_read_gwas_catlog(gwas_catlog_dir, gwas_catlog_file, 
        verbose = TRUE)
    
    
    
    df_hap <- df_hap[!is.na(df_hap$rsID), ]
    # df_proxy <- lookupHub::LDlink_ldproxy(SNP_list = df$SNP,pop =
    # 'EUR',ncpus = 2) gplots::venn(list(df_hap$rsID,df_proxy$proxy_SNP))
    
    
    df_merge <- merge(df_gwascatlog, df_hap, by.x = "SNPS", by.y = "rsID", 
        all.y = TRUE)
    
    
    df_merge <- as.data.table(df_merge)
    
    df_merge[!is.na(PUBMEDID), `:=`(GWAS_report_proxy_SNPS, paste0(PUBMEDID, 
        "_", `DISEASE/TRAIT`, "_", MAPPED_GENE, "_", SNPS, "_Rsq", r2, "_", 
        `P-VALUE`))]
    
    
    df_merge[, `:=`(GWAS_report_proxy_SNPS, ifelse(is.na(`P-VALUE (TEXT)`), 
        GWAS_report_proxy_SNPS, paste0(GWAS_report_proxy_SNPS, `P-VALUE (TEXT)`)))]
    
    # dd<-df_merge[!is.na(GWAS_report_proxy_SNPS),
    # .SD[duplicated(GWAS_report_proxy_SNPS)], by = query_snp_rsid]
    
    
    
    df_merge[, `:=`(GWAS_report_proxy_SNPS, paste(GWAS_report_proxy_SNPS, 
        collapse = " || ")), by = query_snp_rsid]
    
    df_merge[, `:=`(GWAS_report_proxy_SNPS, stringr::str_remove_all(GWAS_report_proxy_SNPS, 
        "(NA \\|\\| )|(NA)|\n"))]
    
    # dd<-df_merge[1:10,]
    
    df_merge[, `:=`(GWAS_report_gene, "")]
    
    # var_i<-1
    for (var_i in 1:NROW(df_merge)) {
        print(var_i)
        # print(paste0(df_merge$RefSeq_name[var_i],'_',df_merge$GENCODE_name[var_i]))
        dd <- df_gwascatlog[grepl(df_merge$RefSeq_name[var_i], `REPORTED GENE(S)`) | 
            grepl(df_merge$RefSeq_name[var_i], MAPPED_GENE) | grepl(df_merge$GENCODE_name[var_i], 
            `REPORTED GENE(S)`) | grepl(df_merge$GENCODE_name[var_i], MAPPED_GENE), 
            ]
        # print(dd)
        dd[, `:=`(GWAS_report_gene, paste0(PUBMEDID, "_", `DISEASE/TRAIT`, 
            "_", SNPS, "_", `REPORTED GENE(S)`, "_or_", MAPPED_GENE))]
        dd <- dd$GWAS_report_gene
        dd <- paste(unique(dd), sep = "; ", collapse = " || ")
        df_merge[var_i, `:=`(GWAS_report_gene, dd)]
        
    }
    
    
    dd <- df_merge[is_query_snp == 1]
    # table(duplicated(dd$query_snp_rsid))
    
    dd <- dd[order(chr, pos_hg38), ]
    dd <- dd[!duplicated(SNPS), list(SNPS, chr, pos_hg38, RefSeq_name, GWAS_report_proxy_SNPS, 
        GWAS_report_gene)]
    return(list(result_all = df_merge, result_query = dd))
    
}





f_Haploreg_gaws_catlog_lookup <- function(SNP_list = NULL, gwas_catlog_dir = NULL, 
    gwas_catlog_file = NULL, ldThresh = 0.8) {
    
    
    df_hap <- lookupHub::HaploReg(SNP_list = SNP_list)
    
    # fwrite(df_hap,'out_gwas/DISC_rank_chunk_mean/query_from_haploreg.csv',na='NA',quote
    # = FALSE)
    
    df_gwascatlog <- f_read_gwas_catlog(gwas_catlog_dir, gwas_catlog_file, 
        verbose = TRUE)
    
    
    
    df_hap <- df_hap[!is.na(df_hap$rsID), ]
    # df_proxy <- lookupHub::LDlink_ldproxy(SNP_list = df$SNP,pop =
    # 'EUR',ncpus = 2) gplots::venn(list(df_hap$rsID,df_proxy$proxy_SNP))
    
    
    df_merge <- merge(df_gwascatlog, df_hap, by.x = "SNPS", by.y = "rsID", 
        all.y = TRUE)
    
    
    df_merge <- as.data.table(df_merge)
    
    df_merge[!is.na(PUBMEDID), `:=`(GWAS_report_proxy_SNPS, paste0(PUBMEDID, 
        "_", `DISEASE/TRAIT`, "_", MAPPED_GENE, "_", SNPS, "_Rsq", r2, "_", 
        `P-VALUE`))]
    
    
    df_merge[, `:=`(GWAS_report_proxy_SNPS, ifelse(is.na(`P-VALUE (TEXT)`), 
        GWAS_report_proxy_SNPS, paste0(GWAS_report_proxy_SNPS, `P-VALUE (TEXT)`)))]
    
    # dd<-df_merge[!is.na(GWAS_report_proxy_SNPS),
    # .SD[duplicated(GWAS_report_proxy_SNPS)], by = query_snp_rsid]
    
    
    
    df_merge[, `:=`(GWAS_report_proxy_SNPS, paste(GWAS_report_proxy_SNPS, 
        collapse = " || ")), by = query_snp_rsid]
    
    df_merge[, `:=`(GWAS_report_proxy_SNPS, stringr::str_remove_all(GWAS_report_proxy_SNPS, 
        "(NA \\|\\| )|(NA)|\n"))]
    
    # dd<-df_merge[1:10,]
    
    df_merge[, `:=`(GWAS_report_gene, "")]
    
    # var_i<-1
    for (var_i in 1:NROW(df_merge)) {
        print(var_i)
        # print(paste0(df_merge$RefSeq_name[var_i],'_',df_merge$GENCODE_name[var_i]))
        dd <- df_gwascatlog[grepl(df_merge$RefSeq_name[var_i], `REPORTED GENE(S)`) | 
            grepl(df_merge$RefSeq_name[var_i], MAPPED_GENE) | grepl(df_merge$GENCODE_name[var_i], 
            `REPORTED GENE(S)`) | grepl(df_merge$GENCODE_name[var_i], MAPPED_GENE), 
            ]
        # print(dd)
        dd[, `:=`(GWAS_report_gene, paste0(PUBMEDID, "_", `DISEASE/TRAIT`, 
            "_", SNPS, "_", `REPORTED GENE(S)`, "_or_", MAPPED_GENE))]
        dd <- dd$GWAS_report_gene
        dd <- paste(unique(dd), sep = "; ", collapse = " || ")
        df_merge[var_i, `:=`(GWAS_report_gene, dd)]
        
    }
    
    
    dd <- df_merge[is_query_snp == 1]
    # table(duplicated(dd$query_snp_rsid))
    
    dd <- dd[order(chr, pos_hg38), ]
    dd <- dd[!duplicated(SNPS), list(SNPS, chr, pos_hg38, RefSeq_name, GWAS_report_proxy_SNPS, 
        GWAS_report_gene)]
    return(list(result_all = df_merge, result_query = dd))
    
}


