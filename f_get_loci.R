

# data = df_sig; SNP = 'SNP' P='P_BOLT_LMM'; locus_name = 'locus'; chr =
# 'CHR.x'; bp = 'BP.x' window = 4e5; p_threshold = 5e-8



# test two methods, second one is better, but check high LD region





# based on roll the window
f_get_loci_window <- function(data = NULL, SNP = "SNP", chr = "CHR", bp = "BP", 
    P = "P", locus_name = "locus", window = 5e+05, p_threshold = 5e-08) {
    
    df <- copy(data)
    
    if (locus_name %in% names(df)) {
        warnings("Already have `locus_name` column, will overwriten.")
    }
    
    
    setnames(df, c(SNP, chr, bp, P), c("SNP_loci_XYZ", "CHR_loci_XYZ", "BP_loci_XYZ", 
        "P_loci_XYZ"))
    
    if (any(duplicated(df$SNP_loci_XYZ))) 
        stop("Duplicated SNPs")
    
    if (NROW(df[P_loci_XYZ >= p_threshold]) == 0) {
    } else {
        cat(paste0(NROW(df[P_loci_XYZ >= p_threshold]), " SNPs removed with P values >= ", 
            p_threshold, ".\n"))
        df <- df[P_loci_XYZ < p_threshold]
    }
    
    df <- df[order(CHR_loci_XYZ, BP_loci_XYZ), ]
    
    chr_loop <- unique(df$CHR_loci_XYZ)
    df[, `:=`(c(locus_name), 0)]  # default locus is 0
    
    for (i_chr in chr_loop) {
        # print(i_chr)
        df_chr <- df[CHR_loci_XYZ == i_chr]
        while (NROW(df_chr) > 0) {
            
            df_chr_each <- df_chr[, .SD[which.min(P_loci_XYZ)]][1, ]  # if P value is the same?
            df[SNP_loci_XYZ == df_chr_each$SNP_loci_XYZ, `:=`(c(locus_name), 
                1)]
            
            # print(paste0(df_chr_each$CHR_loci_XYZ,'-',df_chr_each$BP_loci_XYZ))
            
            df_chr <- df_chr[!(BP_loci_XYZ >= df_chr_each$BP_loci_XYZ - window & 
                BP_loci_XYZ <= df_chr_each$BP_loci_XYZ + window), ]
            
        }
    }
    cat(paste0(sum(unlist(df[, c(locus_name), with = FALSE]))), " locus.\n")
    setnames(df, c("SNP_loci_XYZ", "CHR_loci_XYZ", "BP_loci_XYZ", "P_loci_XYZ"), 
        c(SNP, chr, bp, P))
    
    return(df)
}





f_get_loci <- function(data = NULL, SNP = "SNP", chr = "CHR", bp = "BP", 
    P = "P", locus_name = "locus", window = 5e+05, threshold = 5e-8) {
    df <- copy(data)
    
    if (locus_name %in% names(df)) {
        warnings("Already have `locus_name` column, will overwriten.")
    }
    setnames(df, c(SNP, chr, bp, P), c("SNP_loci_XYZ", "CHR_loci_XYZ", "BP_loci_XYZ", 
        "P_loci_XYZ"))
    if (any(duplicated(df$SNP_loci_XYZ))) 
        stop("Duplicated SNPs")
    
    if (NROW(df[P_loci_XYZ >= threshold]) == 0) {
    } else {
        cat(paste0(NROW(df[P_loci_XYZ >= threshold]), " SNPs removed with P values >= 5e-8.\n"))
        df <- df[P_loci_XYZ < threshold]
    }
    
    
    
    df <- df[order(CHR_loci_XYZ, BP_loci_XYZ), ]
    
    df[, `:=`(BP_distance, list(BP_loci_XYZ * 1 - shift(BP_loci_XYZ))), by = "CHR_loci_XYZ"]
    
    df[, `:=`(BP_window, as.integer(BP_distance < window))]
    df[, `:=`(BP_window, ifelse(is.na(BP_distance), 0, BP_window))]
    
    
    v_count <- 0
    
    for (i_each in 1:NROW(df)) {
        if (df[i_each, BP_window] == 0) {
            v_count <- v_count + 1
            df[i_each, `:=`(locus_index, v_count)]
        } else if (df[i_each, BP_window] == 1) {
            df[i_each, `:=`(locus_index, v_count)]
        }
    }
    
    df_loci <- df[, .SD[which.min(P_loci_XYZ)], by = locus_index]
    
    cat(paste0(max(df_loci$locus_index)), " locus.\n")
    
    df[, `:=`(c(locus_name), ifelse(SNP_loci_XYZ %in% df_loci$SNP_loci_XYZ, 
        1, 0))]
    
    df[, `:=`(c("BP_window"), NULL)]
    
    setnames(df, c("SNP_loci_XYZ", "CHR_loci_XYZ", "BP_loci_XYZ", "P_loci_XYZ"), 
        c(SNP, chr, bp, P))
    
    return(df)
}




f_get_gene <- function(data = NULL, SNP = "SNP", chr = "CHR", bp = "BP", 
    gene_col = "Gene", allow_dup = FALSE) {
    df <- copy(data)
    
    # put as a data in package later
    df_hg19 <- fread("/working/lab_stuartma/xikunH/UKB/OCT/data/glist-hg19")
    setnames(df_hg19, c("CHR", "Start", "End", "Gene"))
    
    df_hg19[, `:=`(CHR, ifelse(CHR == "Y", 24, CHR))]
    df_hg19[, `:=`(CHR, ifelse(CHR %in% c("X", "XY"), 23, CHR))]
    df_hg19[, `:=`(CHR, as.integer(CHR))]
    
    
    if (gene_col %in% names(df)) {
        warnings("Already have `gene_col` column, will overwriten.")
    }
    
    setnames(df, c(SNP, chr, bp), c("SNP_loci_XYZ", "CHR_loci_XYZ", "BP_loci_XYZ"))
    
    if (any(duplicated(df$SNP_loci_XYZ)) & (!allow_dup)) 
        stop("Duplicated SNPs")
    
    df[, `:=`(CHR_loci_XYZ, ifelse(CHR_loci_XYZ == "Y", 24, CHR_loci_XYZ))]
    df[, `:=`(CHR_loci_XYZ, ifelse(CHR_loci_XYZ %in% c("X", "XY"), 23, CHR_loci_XYZ))]
    df[, `:=`(CHR_loci_XYZ, as.integer(CHR_loci_XYZ))]
    
    df <- df[order(CHR_loci_XYZ, BP_loci_XYZ), ]
    df[, `:=`((gene_col), " ")]
    
    for (i_each in 1:NROW(df)) {
        i_chr <- df[i_each, CHR_loci_XYZ]
        i_bp <- df[i_each, BP_loci_XYZ]
        
        df_gene <- df_hg19[CHR == i_chr]
        df_gene[, `:=`(start_dis, i_bp - Start)]
        df_gene[, `:=`(end_dis, i_bp - End)]
        
        df_res <- df_gene[start_dis > 0 & end_dis < 0, ]
        
        
        if (NROW(df_res) > 0) {
            df_res <- paste(df_res$Gene, collapse = "; ")
        } else {
            df_res_start <- df_gene[end_dis > 0][which.min(end_dis)]
            df_res_end <- df_gene[start_dis < 0][which.max(start_dis)]
            
            df_res <- paste(c(df_res_start$Gene, df_res_end$Gene), collapse = " - ")
        }
        print(paste0(i_each, " ", df_res))
        df[i_each, `:=`((gene_col), df_res)]
    }
    
    setnames(df, c("SNP_loci_XYZ", "CHR_loci_XYZ", "BP_loci_XYZ"), c(SNP, 
        chr, bp))
    
    return(df)
}






# chr = 'CHR'; bp = 'BP' p = 'P_BOLT_LMM' snp = 'SNP' data <- df filter_p
# = 0.001


# library(ggplot2) p<-f_manhattan(df,p='P_BOLT_LMM',highlight_data =
# highlight_data,show_text_known = FALSE) p
# ggsave(paste0('result/figure/test.pdf'),width =15,height =7)


f_manhattan <- function(data = NULL, chr = "CHR", bp = "BP", p = "P", snp = "SNP", 
    filter_p = 0.001, highlight_data = NULL, color = if (is.null(highlight_data)) RColorBrewer::brewer.pal(3, 
        "Set1")[1:2] else RColorBrewer::brewer.pal(8, "Pastel2")[c(1, 8)], 
    suggestive_line = NULL, genomewide_line = -log10(5e-08), show_text_known = TRUE, 
    col_point_highlight_novel = "red", col_text_highlight_novel = "black", 
    col_point_highlight_known = "purple", col_text_highlight_known = "purple", 
    seed_number = 1234) {
    
    
    if (!is.data.table(data)) {
        data.table::setDT(data)
    }
    
    
    f_chr_as_integer <- function(data = NULL, chr = "CHR") {
        df <- copy(data)
        setnames(df, c(chr), c("CHR_XYZ"))
        warning(paste0("The class of `CHR` is ", class(df$CHR_XYZ), ", will be converted to integer."))
        cat("\nBefore convert `CHR`\n")
        print(df[, table(CHR_XYZ, useNA = "ifany")])
        cat("\nAfter convert `CHR`\n")
        
        df[, `:=`(CHR_XYZ, as.character(CHR_XYZ))]
        df[, `:=`(CHR_XYZ, stringr::str_remove_all(CHR_XYZ, "chr"))]
        df[, `:=`(CHR_XYZ, ifelse(CHR_XYZ %in% c("X", "XY"), "23", CHR_XYZ))]  ## or change to 25
        df[, `:=`(CHR_XYZ, ifelse(CHR_XYZ %in% c("Y"), "24", CHR_XYZ))]
        df[, `:=`(CHR_XYZ, as.integer(CHR_XYZ))]
        print(df[, table(CHR_XYZ, useNA = "ifany")])
        setnames(df, c("CHR_XYZ"), c(chr))
        return(df)
    }
    
    
    df <- copy(data)
    df <- df[, c(snp, chr, bp, p), with = FALSE]
    
    setnames(df, c(snp, chr, bp, p), c("SNP", "CHR", "BP", "P"))
    
    df <- df[, list(SNP, CHR, BP, P)]
    
    # check chr class
    if (!(class(df$CHR) == "integer" || class(df$CHR) == "numeric")) {
        df <- f_chr_as_integer(data = df, chr = "CHR")
    }
    
    if (class(df$BP) != "numeric") {
        df[, `:=`(BP, as.numeric(as.character(BP)))]
    }
    
    if (class(df$P) != "numeric") {
        df[, `:=`(P, as.numeric(as.character(P)))]
    }
    
    if (filter_p < 1 && filter_p > 0) {
        df <- df[P < filter_p]
    }
    v_n <- NROW(df)
    cat("\n Filter SNPs with P values < ", filter_p, ". ", v_n, "SNPs used in manhattan plot. \n")
    
    df <- df[complete.cases(df)]
    if (NROW(df) != v_n) {
        warning(paste0(v_n - NROW(df), " SNPs with missing values for CHR, BP, or P value are removed."))
    }
    df <- df[order(CHR, BP)]
    
    df[, `:=`(logp, -log10(P))]
    
    df[, `:=`(pos, NA_real_)]
    df[, `:=`(index, NA_integer_)]
    
    v_count <- 0L
    for (i_each in unique(df$CHR)) {
        v_count <- v_count + 1
        df[CHR == i_each, `:=`(index, v_count)]
    }
    
    nchr <- length(unique(df$CHR))
    
    if (nchr == 1) {
        df[, `:=`(pos, BP)]
        ticks = floor(length(df$pos)/2) + 1
        xlabel = paste("Chromosome", unique(df$CHR), "position")
        labs = ticks
    } else {
        lastbase = 0
        ticks = NULL
        for (i_each in unique(df$index)) {
            if (i_each == 1) {
                df[index == i_each, `:=`(pos, BP)]
            } else {
                lastbase <- lastbase + unlist(df[index == i_each - 1, list(BP)][.N, 
                  ])
                df[index == i_each, `:=`(pos, BP + lastbase)]
            }
            ticks <- c(ticks, (min(df[index == i_each, pos]) + max(df[index == 
                i_each, pos]))/2 + 1)
        }
        xlabel = "Chromosome"
        labs <- unique(df$CHR)
    }
    
    set.seed(seed_number)
    
    
    
    p_ggplot2 <- ggplot(data = df, aes(pos, logp, colour = factor(CHR))) + 
        guides(colour = FALSE) + geom_point() + labs(y = expression(bold(-log["10"](italic(p))))) + 
        scale_x_continuous(name = "Chromosome", breaks = ticks, labels = labs) + 
        scale_colour_manual(values = rep(color, 12))
    
    if (!is.null(suggestive_line)) {
        p_ggplot2 <- p_ggplot2 + geom_hline(yintercept = suggestive_line, 
            colour = "blue", size = 0.1)
    }
    
    if (!is.null(genomewide_line)) {
        p_ggplot2 <- p_ggplot2 + geom_hline(yintercept = genomewide_line, 
            colour = "red", size = 0.1)
    }
    if (!is.null(highlight_data)) {
        
        df_highlight <- as.data.table(highlight_data)
        if (!all(c("SNP", "CHR", "BP") %in% names(df_highlight))) {
            stop("highlight_data should at least have three colums: SNP, CHR, BP.")
        }
        
        if (!("novel" %in% names(df_highlight))) {
            df_highlight[, `:=`(novel, 1)]
        }
        
        if (!("gene" %in% names(df_highlight))) {
            df_highlight[, `:=`(gene, SNP)]
        }
        
        if (!(class(df_highlight$CHR) == "integer" || class(df_highlight$CHR) == 
            "numeric")) {
            df_highlight <- f_chr_as_integer(data = df_highlight, chr = "CHR")
        }
        
        df_select <- df[SNP %in% df_highlight$SNP, ]
        df_highlight <- merge(df_select, df_highlight, by = "SNP", suffixes = c("", 
            "_post"))
        cat("\n Dimension of hightlight data after merge with coordinate data.\n")
        print(dim(df_highlight))
        cat("\nThe first head rows:\n")
        print(head(df_highlight))
        cat("\nNovel loci:\n")
        print(table(df_highlight$novel))
        
        df_highlight[, `:=`(col_point, ifelse(novel == 1, col_point_highlight_novel, 
            col_point_highlight_known))]
        df_highlight[, `:=`(col_text, ifelse(novel == 1, col_text_highlight_novel, 
            col_text_highlight_known))]
        
        df_highlight_novel <- df_highlight[novel == 1, ]
        
        p_ggplot2 <- p_ggplot2 + geom_point(data = df_highlight, aes(pos, 
            logp), color = df_highlight$col_point)
        
        if (show_text_known) {
            p_ggplot2 <- p_ggplot2 + ggrepel::geom_text_repel(data = df_highlight, 
                aes(x = pos, y = logp), label = df_highlight$gene, cex = 3.5, 
                color = df_highlight$col_text, fontface = "bold", point.padding = unit(0.1, 
                  "lines"), segment.color = "grey50")
        } else {
            p_ggplot2 <- p_ggplot2 + ggrepel::geom_text_repel(data = df_highlight_novel, 
                aes(x = pos, y = logp), label = df_highlight_novel$gene, 
                cex = 4, color = df_highlight_novel$col_text, fontface = "bold", 
                point.padding = unit(0.1, "lines"), segment.color = "grey50")
        }
    }
    
    p_ggplot2 <- p_ggplot2 + theme(legend.position = "none", panel.background = element_blank(), 
        axis.line.x = element_line(colour = "black", size = 1.1), axis.line.y = element_line(colour = "black", 
            size = 1.1), axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"), 
        text = element_text(face = "bold", size = 14))
    return(p_ggplot2)
}


