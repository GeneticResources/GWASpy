library(ggplot2)


# example f_xyplot(data=df, beta_x = 'mtag_beta_four_trait', beta_y =
# 'Effect_glaucoma_ANZRAG_flip', se_x = 'mtag_se_four_trait', se_y =
# 'StdErr_glaucoma_ANZRAG', aes_col = 'mtag_color_group', aes_col_value
# =var_col, lab_x = 'Effect Size on MTAG glaucoma', lab_y = 'Effect Size
# on ANZRAG glaucoma', show_errbar = FALSE)

f_xyplot <- function(data = NULL, beta_x = NULL, beta_y = NULL, se_x = NULL, 
    se_y = NULL, aes_col = NULL, aes_col_value = NULL, 
    aes_shape = NULL, aes_shape_value = NULL, 
    lab_x = "Effect Size on X", 
    lab_y = "Effect Size on Y", show_errbar = TRUE, scale_errbar = 1.96, 
    lab_snp = NULL) {
    data <- as.data.table(data)
    
    data[, `:=`(c("ymin_errbar", "ymax_errbar", "xmin_errbar", "xmax_errbar"), 
        list(get(beta_y) - get(se_y) * scale_errbar, get(beta_y) + get(se_y) * 
            scale_errbar, get(beta_x) - get(se_x) * scale_errbar, get(beta_x) + 
            get(se_x) * scale_errbar))]
    
    print(data[1:5, c(1:10)])
    if (is.null(aes_col_value)) {
        aes_col_value <- RColorBrewer::brewer.pal(9, "Set1")
    }
    
    res_cor <- cor.test(~data[, get(beta_x)] + data[, get(beta_y)])
    print(res_cor)
    print(summary(lm(data[, get(beta_y)] ~ data[, get(beta_x)])))
    
    p <- ggplot(data = copy(data), aes_string(x = beta_x, y = beta_y))
    
    if (show_errbar) {
        
        p <- p + geom_errorbar(aes(ymin = ymin_errbar, ymax = ymax_errbar), 
            colour = "grey")
        
        p <- p + geom_errorbarh(aes(xmin = xmin_errbar, xmax = xmax_errbar), 
            colour = "grey")
        # print(p)
    }
    
    
    p <- p + geom_point(aes_string(color = aes_col,shape = aes_shape)) + 
        scale_color_manual(values = aes_col_value) + 
        scale_shape_manual(values = aes_shape_value) +
        geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, 
        linetype = "dashed") + labs(x = lab_x, y = lab_y, title = "") + theme(plot.title = element_text(hjust = 0.5)) + 
        stat_smooth(color = "red", size = 0.5, se = TRUE, method = "lm", 
            fill = "grey", aes(weight = 1/(data[, get(se_y)]^2)))
    if (!is.null(lab_snp)) {
        p <- p + ggrepel::geom_text_repel(label = lab_snp, cex = 2)
    }
    # scale_x_continuous(limits = c(-0.25, 0.35))+ scale_y_continuous(limits
    # = c(-0.35, 0.35))+
    p <- p + theme(legend.justification = c(0, 0), legend.position = c(0.6, 
        0.05), legend.title = element_blank(), legend.text = element_text(size = 10), 
        legend.background = element_blank(), legend.key = element_blank(), 
        panel.background = element_blank(), legend.box.background = element_rect(colour = "black"), 
        axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), 
        text = element_text(face = "bold", size = 12))
    return(p)
    
}



f_xyplot_shape <- function(data = NULL, beta_x = NULL, beta_y = NULL, se_x = NULL, 
                     se_y = NULL, aes_shape = NULL, aes_shape_value = NULL, lab_x = "Effect Size on X", 
                     lab_y = "Effect Size on Y", show_errbar = TRUE, scale_errbar = 1.96, 
                     lab_snp = NULL) {
    data <- as.data.table(data)
    
    data[, `:=`(c("ymin_errbar", "ymax_errbar", "xmin_errbar", "xmax_errbar"), 
                list(get(beta_y) - get(se_y) * scale_errbar, get(beta_y) + get(se_y) * 
                         scale_errbar, get(beta_x) - get(se_x) * scale_errbar, get(beta_x) + 
                         get(se_x) * scale_errbar))]
    
    print(data[1:5, c(1:10)])
    if (is.null(aes_shape_value)) {
        aes_col_value <- c(1:9)
    }
    
    res_cor <- cor.test(~data[, get(beta_x)] + data[, get(beta_y)])
    print(res_cor)
    print(summary(lm(data[, get(beta_y)] ~ data[, get(beta_x)])))
    
    p <- ggplot(data = copy(data), aes_string(x = beta_x, y = beta_y))
    
    if (show_errbar) {
        
        p <- p + geom_errorbar(aes(ymin = ymin_errbar, ymax = ymax_errbar), 
                               colour = "grey")
        
        p <- p + geom_errorbarh(aes(xmin = xmin_errbar, xmax = xmax_errbar), 
                                colour = "grey")
        # print(p)
    }
    
    
    p <- p + geom_point(aes_string(shape = aes_shape)) + scale_shape_manual(values = aes_shape_value) + 
        geom_hline(yintercept = 0, linetype = "dashed") + geom_vline(xintercept = 0, 
                                                                     linetype = "dashed") + labs(x = lab_x, y = lab_y, title = "") + theme(plot.title = element_text(hjust = 0.5)) + 
        stat_smooth(color = "red", size = 0.5, se = TRUE, method = "lm", 
                    fill = "grey", aes(weight = 1/(data[, get(se_y)]^2)))
    if (!is.null(lab_snp)) {
        p <- p + ggrepel::geom_text_repel(label = lab_snp, cex = 2)
    }
    # scale_x_continuous(limits = c(-0.25, 0.35))+ scale_y_continuous(limits
    # = c(-0.35, 0.35))+
    p <- p + theme(legend.justification = c(0, 0), legend.position = c(0.6, 
                                                                       0.05), legend.title = element_blank(), legend.text = element_text(size = 10), 
                   legend.background = element_blank(), legend.key = element_blank(), 
                   panel.background = element_blank(), legend.box.background = element_rect(colour = "black"), 
                   axis.line.x = element_line(colour = "black"), axis.line.y = element_line(colour = "black"), 
                   text = element_text(face = "bold", size = 12))
    return(p)
    
}
