format_sprint<-function(x) {
    x<-ifelse(abs(x)<0.01,sprintf(fmt="%.1E",x),sprintf(fmt="%.2f",x))
    x
}