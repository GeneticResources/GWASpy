
library(ggplot2)

PRS_cut_interval <- function(var_df, number = 3, ref = "01") {
    var_df <- cut_number(var_df, number)
    var_df <- factor(var_df, levels = paste0(levels(var_df)), labels = sprintf("%02.f", 
        1:number))
    var_df <- relevel(var_df, ref = ref)
    return(var_df)
}
