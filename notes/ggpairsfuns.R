trim_gg <- function(gg,hack_spaces=TRUE) {
    n <- gg$nrow
    gg$nrow <- gg$ncol <- n-1
    v <- 1:n^2
    gg$plots <- gg$plots[v>n & v%%n!=0]
    gg$xAxisLabels <- gg$xAxisLabels[-n]
    gg$yAxisLabels <- gg$yAxisLabels[-1]
    ## ggpairs tries to parse labels, so this fails if they
    ## have spaces in them ... maybe there's a better way?
    if (hack_spaces) {
        gg$xAxisLabels <- gsub("_"," ",gg$xAxisLabels)
        gg$yAxisLabels <- gsub("_"," ",gg$yAxisLabels)
    }
    return(gg)
}

tweak_legends_gg <- function(gg,pos=1.5) {
    n <- gg$nrow
    for (i in 1:n) {
        for (j in 1:n) {
            inner <- getPlot(gg,i,j)
            if (i==1 & j==1) {
                inner <- inner + ggplot2::theme(legend.position=c(pos,0.5)) +
                    guides(colour = guide_legend(override.aes = list(size=6)))
            } else {
                inner <- inner + ggplot2::theme(legend.position="none")
            }
            gg <- putPlot(gg,inner,i,j)
        }
    }
    return(gg)
}
tweak_colours_gg <- function(gg,sc=scale_color_viridis()) {
    n <- gg$nrow
    for (i in 1:n) {
        for (j in 1:n) {
            inner <- getPlot(gg,i,j)
            inner <- inner+sc
            gg <- putPlot(gg,inner,i,j)
        }
    }
    return(gg)
}
