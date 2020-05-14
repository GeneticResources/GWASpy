
# based on qqman, use ggolot2 style and avoid overlap for annotations.


library(calibrate)


RColorBrewer::display.brewer.all()
# var_col <- c(RColorBrewer::brewer.pal(8, "Paired")[c( 1)],
#              RColorBrewer::brewer.pal(8, "Pastel2")[c( 8)])
# var_col <- RColorBrewer::brewer.pal(8, "Pastel2")[c(1, 8)]

f_textplot <- function(x, y, words, cex = 1, new = TRUE, show.lines = TRUE, 
    col = "red", ...) {
    if (new) 
        plot(x, y, type = "n", ...)
    lay <- wordcloud::wordlayout(x, y, words, cex, rotate90 = FALSE, tstep = 0.1, 
        rstep = 0.1, ...)
    if (show.lines) {
        for (i in 1:length(x)) {
            xl <- lay[i, 1]
            yl <- lay[i, 2]
            w <- lay[i, 3]
            h <- lay[i, 4]
            if (x[i] < xl || x[i] > xl + w || y[i] < yl || y[i] > yl + h) {
                points(x[i], y[i], pch = 16, col = col, cex = 0.1)
                nx <- xl + 0.5 * w
                ny <- yl + 0.5 * h
                lines(c(x[i], nx), c(y[i], ny + 1.5), col = RColorBrewer::brewer.pal(8, 
                  "Set1")[3])
            }
        }
    }
    text(lay[, 1] + 0.5 * lay[, 3] + 1, lay[, 2] + 0.5 * lay[, 4] + 2, words, 
        cex = cex, col = col, srt = 60, ...)
}



library(ggplot2)

f_manhattan <- function(x, chr = "CHR", bp = "BP", p = "P", snp = "SNP", 
    col = c("gray10", "gray60"), chrlabs = NULL, suggestiveline = -log10(1e-05), 
    genomewideline = -log10(5e-08), highlight = NULL, highlight_known = NULL, 
    highlight_ggplot2 = NULL, text_known = TRUE, logp = TRUE, loglogp=TRUE, annotatePval = NULL, 
    annotateTop = TRUE, var_hight_col = RColorBrewer::brewer.pal(12, "Paired")[6], 
    var_hight_known_col = RColorBrewer::brewer.pal(12, "Paired")[10], only_position = FALSE, 
    ...) {
    CHR = BP = P = index = NULL
    if (!(chr %in% names(x))) 
        stop(paste("Column", chr, "not found!"))
    if (!(bp %in% names(x))) 
        stop(paste("Column", bp, "not found!"))
    if (!(p %in% names(x))) 
        stop(paste("Column", p, "not found!"))
    if (!(snp %in% names(x))) 
        warning(paste("No SNP column found. OK unless you're trying to highlight."))
    if (!is.numeric(x[[chr]])) 
        stop(paste(chr, "column should be numeric. Do you have 'X', 'Y', 'MT', etc? If so change to numbers and try again."))
    if (!is.numeric(x[[bp]])) 
        stop(paste(bp, "column should be numeric."))
    if (!is.numeric(x[[p]])) 
        stop(paste(p, "column should be numeric."))
    d = data.frame(CHR = x[[chr]], BP = x[[bp]], P = x[[p]])
    if (!is.null(x[[snp]])) 
        d = transform(d, SNP = x[[snp]])
    d <- subset(d, (is.numeric(CHR) & is.numeric(BP) & is.numeric(P)))
    d <- d[order(d$CHR, d$BP), ]
    # if (loglogp) {
    #     d$logp <- log10(-log10(d$P))
    # } 
    
    if(logp) { 
        d$logp <- -log10(d$P)
    } else {
        d$logp <- d$P
    }
    d$pos = NA
    d$index = NA
    ind = 0
    for (i in unique(d$CHR)) {
        ind = ind + 1
        d[d$CHR == i, ]$index = ind
    }
    nchr = length(unique(d$CHR))
    if (nchr == 1) {
        d$pos = d$BP
        ticks = floor(length(d$pos))/2 + 1
        xlabel = paste("Chromosome", unique(d$CHR), "position")
        labs = ticks
    } else {
        lastbase = 0
        ticks = NULL
        for (i in unique(d$index)) {
            if (i == 1) {
                d[d$index == i, ]$pos = d[d$index == i, ]$BP
            } else {
                lastbase = lastbase + tail(subset(d, index == i - 1)$BP, 
                  1)
                d[d$index == i, ]$pos = d[d$index == i, ]$BP + lastbase
            }
            ticks = c(ticks, (min(d[d$index == i, ]$pos) + max(d[d$index == 
                i, ]$pos))/2 + 1)
        }
        xlabel = "Chromosome"
        labs <- unique(d$CHR)
    }
    xmax = ceiling(max(d$pos) * 1.03)
    xmin = floor(max(d$pos) * -0.03)
    def_args <- list(xaxt = "n", bty = "n", xaxs = "i", yaxs = "i", las = 1, 
        pch = 20, xlim = c(xmin, xmax), ylim = c(0, ceiling(max(d$logp))), 
        xlab = xlabel, ylab = expression(-log[10](italic(p))))
    dotargs <- list(...)
    do.call("plot", c(NA, dotargs, def_args[!names(def_args) %in% names(dotargs)]))
    if (!is.null(chrlabs)) {
        if (is.character(chrlabs)) {
            if (length(chrlabs) == length(labs)) {
                labs <- chrlabs
            } else {
                warning("You're trying to specify chromosome labels but the number of labels != number of chromosomes.")
            }
        } else {
            warning("If you're trying to specify chromosome labels, chrlabs must be a character vector")
        }
    }
    if (nchr == 1) {
        axis(1, ...)
    } else {
        axis(1, at = ticks, labels = labs, ...)
    }
    col = rep(col, max(d$CHR))
    if (nchr == 1) {
        with(d, points(pos, logp, pch = 20, col = col[1], ...))
    } else {
        icol = 1
        for (i in unique(d$index)) {
            with(d[d$index == unique(d$index)[i], ], points(pos, logp, col = col[icol], 
                pch = 20, ...))
            icol = icol + 1
        }
    }
    if (suggestiveline) 
        abline(h = suggestiveline, col = "blue")
    if (genomewideline) 
        abline(h = genomewideline, col = "red")
    if (!is.null(highlight)) {
        print(names(highlight))
        print(names(d))
        if (any(!(highlight$SNP %in% d$SNP))) 
            warning("You're trying to highlight SNPs that don't exist in your results.")
        d.highlight <- d[which(d$SNP %in% highlight$SNP), ]
        d.highlight <- base::merge(d.highlight, highlight, by = "SNP", suffixes = c("", 
            "_post"))
        print(dim(d.highlight))
        
        par(xpd = TRUE)
        
        with(d.highlight, f_textplot(pos, -log10(P), d.highlight$gene_hight, 
            cex = 0.7, new = FALSE, col = "black"), ...)
        with(d.highlight, points(pos, logp, col = var_hight_col, pch = 20, 
            cex = 0.7, ...))
    }
    
    
    if (!is.null(highlight_known)) {
        print(names(highlight_known))
        print(names(d))
        if (any(!(highlight_known$SNP %in% d$SNP))) 
            warning("You're trying to highlight_known SNPs that don't exist in your results.")
        d.highlight_known <- d[which(d$SNP %in% highlight_known$SNP), ]
        d.highlight_known <- base::merge(d.highlight_known, highlight_known, 
            by = "SNP", suffixes = c("", "_post"))
        print(dim(d.highlight_known))
        
        # print(d.highlight_known)
        
        
        par(xpd = TRUE)
        
        if (text_known) {
            with(d.highlight_known, f_textplot(pos, -log10(P), d.highlight_known$gene_hight, 
                cex = 0.5, new = FALSE, col = var_hight_known_col), ...)
        }
        
        
        with(d.highlight_known, points(pos, logp, col = var_hight_known_col, 
            pch = 20, cex = 0.7, ...))
    }
    
    print(head(d))
    
    if (only_position) {
        return(list(data = d, ticks = ticks))
    }
    
    if (!is.null(highlight_ggplot2)) {
        d.highlight <- d[which(d$SNP %in% highlight_ggplot2$SNP), ]
        d.highlight <- base::merge(d.highlight, highlight_ggplot2, by = "SNP", 
            suffixes = c("", "_post"))
        print(dim(d.highlight))
       print(d.highlight)
        print(table(d.highlight$novel))
        
        # d.highlight$col_novel<-ifelse(d.highlight$novel==1,var_hight_col,var_hight_known_col)
        
        d.highlight$col_novel <- ifelse(d.highlight$novel == 1, "red", var_hight_known_col)
        d.highlight$col_text <- ifelse(d.highlight$novel == 1, "black", var_hight_known_col)
        
        d.highlight_novel <- d.highlight[d.highlight$novel == 1, ]
        
        p_ggplot2 <- ggplot(data = d, aes(pos, logp, colour = factor(CHR))) + 
            guides(colour = FALSE) + geom_point() + 
            # labs(y = expression(bold("Z"))) + 
            # labs(y = expression(bold(-log["10"](italic(p))))) + 
            labs(y = expression(bold(log["10"](-log["10"](italic(p)))))) + 
            scale_x_continuous(name = "Chromosome", breaks = ticks, labels = (unique(d$CHR))) + 
            scale_colour_manual(values = rep(var_col, 12)) + geom_hline(yintercept = genomewideline, 
            colour = "red", size = 0.1) + geom_point(data = d.highlight, 
            aes(pos, logp), color = d.highlight$col_novel)
        if (text_known) {
            p_ggplot2 <- p_ggplot2 + ggrepel::geom_text_repel(data = d.highlight, 
                aes(x = pos, y = logp), label = d.highlight$gene, cex = 3.5, 
                color = d.highlight$col_text, fontface = "bold", point.padding = unit(0.1, 
                  "lines"), segment.color = "grey50")
        } else {
            p_ggplot2 <- p_ggplot2 + ggrepel::geom_text_repel(data = d.highlight_novel, 
                aes(x = pos, y = logp), label = d.highlight_novel$gene, 
                cex = 6, color = d.highlight_novel$col_text, fontface = "bold", 
                point.padding = unit(0.1, "lines"), segment.color = "grey50")
        }
        
        
        p_ggplot2 <- p_ggplot2 + theme(legend.position = "none", panel.background = element_blank(), 
            axis.line.x = element_line(colour = "black", size = 1.1), axis.line.y = element_line(colour = "black", 
                size = 1.1), axis.text.x = element_text(colour = "black"), 
            axis.text.y = element_text(colour = "black"), text = element_text(face = "bold", 
                size = 14))
        
    }
    
    
    if (!is.null(annotatePval)) {
        topHits = subset(d, P <= annotatePval)
        par(xpd = TRUE)
        if (annotateTop == FALSE) {
            with(subset(d, P <= annotatePval), textplot(pos, -log10(P), offset = 0.3, 
                labs = topHits$SNP, cex = 0.2), ...)
        } else {
            topHits <- topHits[order(topHits$P), ]
            topSNPs <- NULL
            for (i in unique(topHits$CHR)) {
                chrSNPs <- topHits[topHits$CHR == i, ]
                topSNPs <- rbind(topSNPs, chrSNPs[1, ])
            }
            textxy(topSNPs$pos, -log10(topSNPs$P), offset = 0.625, labs = topSNPs$SNP, 
                cex = 0.5, ...)
        }
    }
    par(xpd = FALSE)
    return(list(data = d, ticks = ticks, p_ggplot2 = p_ggplot2))
}
