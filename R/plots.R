# GGPLOT functions #############################################################

#' Multiplot for ggplot2
#'
#' Multiplot for ggplot2 by Winston Chang from
#' \href{http://www.cookbook-r.com/}{Cookbook for R}
#'
#' @param ... ggplot2. Input plot objects
#' @param plotlist list. (optional) Input ggplots as a list
#' @param cols numeric. Number of columns to distribute plots
#' @param layout matrix. Optional layout to order plots
#'
#' @return NULL. plots input ggplots
#' @export
#'
#' @author Winston Chang
multiplot <- function(..., plotlist = NULL, cols = 1, layout = NULL) {
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    # If layout is NULL, then use 'cols' to determine layout
    if(is.null(layout)){
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    if(numPlots==1){
        print(plots[[1]])
    }else{
        # Set up the page
        grid::grid.newpage()
        grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow(layout), ncol(layout))))
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            print(plots[[i]], vp = grid::viewport(layout.pos.row = matchidx$row,
                                                  layout.pos.col = matchidx$col))
        }
    }
}

#' Interval heatmap
#'
#' Plots the mean value of a value Z when divided into groups defined by
#' intervals at variables X and Y
#'
#' @param x numeric. X dimension
#' @param y numeric. Y dimension
#' @param z logical or numeric. Z dimension to be summarized
#' @param x_n numeric. Number of groups for x dimension to be cut into
#' @param y_n numeric. Number of groups for y dimension to be cut into
#' @param nCores numeric. Number of cores to use to process data
#' @param dig.lab numeric. Number of decimal positions to place in labels
#' @param outputMatrix logical. Set to TRUE to output instead the matrix
#' @param dispNum logical. Set to FALSE for not showing the numeric labels at
#' each cell
#' @param na.rm logical. Indicating if NA values should be removed to calculate
#' mean Z. Default: TRUE.
#' @param FUN character. Name of function to be used over the values that fall
#' inside the combinations of x and y intervals. e.g. "mean" or "sum".
#'
#' @return matrix / pheatmap
#' @export
#'
#' @examples
#' intervalHeatmap(x = iris$Sepal.Length, y = iris$Sepal.Width,
#'                 z = rep(1, nrow(iris)), x_n = 5, y_n = 5, FUN = "sum")
intervalHeatmap <- function(x, y, z, x_n, y_n, nCores = 1, dig.lab = 1,
                            outputMatrix = FALSE, dispNum = TRUE, na.rm = TRUE, FUN = "mean"){
    tmpX_fct <- ggplot2::cut_interval(x, x_n, dig.lab = dig.lab)
    tmpY_fct <- ggplot2::cut_interval(y, y_n, dig.lab = dig.lab)
    tmpO <- parallel::mclapply(mc.cores = nCores, rev(levels(tmpY_fct)), function(yF){
        unlist(lapply(levels(tmpX_fct), function(xF){
            FUN <- get(FUN)
            FUN(z[which(tmpX_fct == xF & tmpY_fct == yF)], na.rm = na.rm)
        }))
    }) %>%
        do.call(what = "rbind") %>%
        magrittr::set_rownames(rev(levels(tmpY_fct))) %>%
        magrittr::set_colnames(levels(tmpX_fct))
    if(outputMatrix){return(tmpO)}
    pheatmap::pheatmap(tmpO, cluster_cols = FALSE, cluster_rows = FALSE,
                       display_numbers = dispNum, number_format = "%.2f")
}

#' ggBoxplot
#'
#' ggplot2 boxplot shortcut
#'
#' @param xFactor factor. Factor vector specifying the groups
#' @param yNumeric numeric. Numeric vector specifying the values per observation
#' @param outLCol character. Color of outliers. NA for hidding outliers.
#' @param outLAlpha numeric. Value from 0 to 1 to set transparency of outliers.
#'
#' @return ggplot
#' @export
#'
#' @examples
#' ggBoxplot(iris$Species, iris$Petal.Length)
ggBoxplot <- function(xFactor, yNumeric, outLCol = NA, outLAlpha = 0.5){
    if(length(xFactor) != length(yNumeric)){stop("Length of vectors is not the same.")}
    ggplot2::ggplot(data = data.frame(x = xFactor, y = yNumeric),
                    ggplot2::aes(x = .data$x, y = .data$y)) +
        ggplot2::geom_boxplot(outlier.colour = outLCol, outlier.size = 1, outlier.alpha = outLAlpha) +
        ggplot2::theme_minimal() + EnvStats::stat_n_text(size = 3, angle = 90) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
}


#' Angle ggplot x labels
#'
#' @param angle Angle used to tilt x labels (default = 90)
#'
#' @return ggplot
#' @export
#'
#' @example
#' ggBoxplot(iris$Species, iris$Petal.Length) + ggAngleXLabs(45)
ggAngleXLabs <- function(angle = 90){
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = angle, vjust = 0.5, hjust=1))
}


# # Boxplot shortcut function for matrices
# ggBoxplot_matrix <- function(matrix, title = "Title", xlab = "x", ylab = "y", outLCol = NA){
#     ggplot(data= reshape2::melt(as.data.frame(matrix)), aes(variable, value)) +
#         geom_boxplot(outlier.colour= outLCol, outlier.size = 1) + xlab(xlab) + ylab(ylab) +
#         ggtitle(title) + theme_classic() +  stat_n_text(size = 3, angle = 90) +
#         theme(axis.text.x = element_text(angle = 90, hjust = 1))
# }
#
# # Pie Chart Nucleotide frequency
# ggPieNucFreq <- function(nucFreq, labSize = 5){
#     palette(brewer.pal(9, "Set1")) #RColorBrewer
#     tmpDF <- data.frame(nucs = names(nucFreq), Percent = nucFreq, stringsAsFactors = F)
#     tmpDF <- tmpDF[order(tmpDF[,2], decreasing = T),]
#     tmpDF <- data.frame(rbind(tmpDF[1:5,], c("All Others", sum(tmpDF[-1:-5,2]))))
#     tmpDF[,2] <- as.numeric(tmpDF[,2]) / sum(as.numeric(tmpDF[,2]))
#     tmpDF[,1] <- factor(tmpDF[,1], levels = tmpDF[,1])
#     ggPie <- ggplot(tmpDF, aes(x="", y=Percent, fill=nucs)) +
#         geom_bar(width = 1, stat = "identity") + coord_polar(theta = "y",start = 0,direction = 1) +
#         geom_text(aes(label = round(Percent,2)), size= labSize, position = position_stack(vjust = 0.5)) +
#         theme(axis.text.x =element_blank()) + theme_classic()
#     return(ggPie)
# }
#
# # GGplot alternative to pairs function (additionally it fits linear models to all pair-wise comparisons)
# ggPairs <- function(DF, alpha = 1){
#     iCol <- colnames(DF)
#     matD <- combinations(n = length(iCol), r = 2, v = 1:length(iCol))
#     ggSC <- lapply(1:nrow(matD), function(x){
#         tmpL <- lm(DF[,matD[x,2]] ~ DF[,matD[x,1]])
#         if(tmpL$coefficients[1]>=0){
#             linModEq = paste("y = x *",tmpL$coefficients[2] %>% signif(2), "+", tmpL$coefficients[1]  %>% signif(2))
#         }else if(tmpL$coefficients[1]<0){linModEq = paste("y = x *", signif(tmpL$coefficients[2],2), "-",
#                                                            tmpL$coefficients[1]  %>% signif(2) %>% abs)}
#         tmpC <- cor(DF[,matD[x,1]], DF[,matD[x,2]], use = "p") %>% round(4)
#         tmpP <- cor.test(DF[,matD[x,1]], DF[,matD[x,2]], use = "p")$p.value %>% signif(4)
#         tmpC2 <- cor(DF[,matD[x,1]], DF[,matD[x,2]], use = "p", method = "spearman") %>% round(4)
#         tmpP2 <- cor.test(DF[,matD[x,1]], DF[,matD[x,2]], use = "p", method = "spearman")$p.value %>% signif(4)
#         ggplot(DF, aes(x= DF[,matD[x,1]], y= DF[,matD[x,2]])) +
#             geom_point(alpha = alpha, shape = 16) +
#             geom_smooth(method = lm, se=FALSE, fullrange= TRUE, aes(group=1)) +
#             geom_abline(intercept = 0, slope = 1, colour = "gray") +
#             theme_classic() + xlab(iCol[matD[x,1]]) + ylab(iCol[matD[x,2]]) +
#             ggtitle(paste("R =", tmpC, "p = ", tmpP, "\nrho =", tmpC2, "p =", tmpP2),
#                     subtitle = linModEq) +
#             coord_cartesian(ylim = range(DF, na.rm = T), xlim = range(DF, na.rm = T))
#     })
#     ggLabs <- lapply(iCol, function(x){
#         df <- data.frame(x = 1, y = 1, text = x)
#         ggO <- ggplot(df, aes(x, y)) +
#             geom_text(aes(label = text), size = 5) + theme_classic() +
#             theme(panel.border = element_rect(colour = 1, fill = NA), axis.line = element_line())+
#             theme(axis.title.x=element_blank(),
#                   axis.text.x=element_blank(),
#                   axis.ticks.x=element_blank()) +
#             theme(axis.title.y=element_blank(),
#                   axis.text.y=element_blank(),
#                   axis.ticks.y=element_blank())
#         return(ggO)
#     })
#     ggCf <- lapply(1:nrow(matD), function(x){
#         return()
#     })
#     lOut <- matrix(NA, ncol = ncol(DF), nrow = ncol(DF))
#     for(i in 1:nrow(matD)){lOut[matD[i,2], matD[i,1]] <- i}
#     for(i in 1:length(iCol)){lOut[i, i] <- length(ggSC) + i}
#     for(i in 1:nrow(matD)){lOut[matD[i,1], matD[i,2]] <- length(ggSC) + length(iCol) + i}
#     multiplot(plotlist = c(ggSC, ggLabs), layout = lOut)
# }
#
# # Simple Barplot function
# ggBarplot <- function(x, ci = NA, title = NULL, subt = NULL, xLab = "Names", yLab = "Values"){
#     if(is.null(names(x))){names(x) <- 1:length(x)}
#     df <- data.frame(names = factor(names(x), levels = names(x)), value = x, CI = ci)
#     outGG <- ggplot(data=df, aes(x=names, y=value)) +
#         geom_bar(stat="identity") + theme_classic() +
#         geom_errorbar(aes(ymin=value-CI, ymax=value+CI), width=.2, position=position_dodge(.9)) +
#         theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(title, subt) +
#         ylab(yLab) + xlab(xLab)
#     return(outGG)
# }
#
# # Scatterplot with linear model fitted line and calculates correlation
# ggScattLinePlot <- function(x, y, title = "", xLab = "", yLab = "", alpha = 1){
#     tmpC <- cor(x, y, use = "p") %>% round(4)
#     tmpP <- cor.test(x, y, use = "p")$p.value %>% signif(3)
#     tmpC2 <- cor(x, y, use = "p", method = "spearman") %>% round(4)
#     tmpP2 <- cor.test(x, y, use = "p", method = "spearman")$p.value %>% signif(3)
#     tmpDF <- data.frame(var1 = x, var2 = y)
#     ggSCLINE <-  ggplot(tmpDF, aes(x = var1, y = var2)) + geom_point(alpha = alpha) +
#         geom_smooth(method = lm, se=FALSE, fullrange= TRUE, aes(group=1)) +
#         ggtitle(title, paste("R =", tmpC, "p = ", tmpP, "\nrho =", tmpC2, "p =", tmpP2)) +
#         ylab(yLab) + xlab(xLab) + theme_classic()
#     return(ggSCLINE)
# }
#
# # ggplot heatmap
# ggHeatmap <- function(x, y, countLogTransform = T, tendencyLine = F, nBins = 100){
#     tmpC <- cor(x, y, use = "p") %>% round(4)
#     tmpP <- cor.test(x, y, use = "p")$p.value %>% signif(3)
#     tmpC2 <- cor(x, y, use = "p", method = "spearman") %>% round(4)
#     tmpP2 <- cor.test(x, y, use = "p", method = "spearman")$p.value %>% signif(3)
#     tmpDF <- data.frame(var1 = x, var2 = y)
#     tmpGG <- ggplot(tmpDF, aes(x = var1, y = var2)) + geom_bin2d(bins = nBins) +
#         theme_minimal() + ggtitle("", paste("R =", tmpC, "p = ", tmpP, "\nrho =", tmpC2, "p =", tmpP2)) +
#         theme(legend.position = "bottom")
#     if(countLogTransform == T){
#         tmpGG <- tmpGG + scale_fill_gradientn(trans = "log", colours = rev(brewer.pal(9, "Spectral")))
#     }else{
#         tmpGG <- tmpGG + scale_fill_gradientn(colours = rev(brewer.pal(9, "Spectral")))
#     }
#     if(tendencyLine == T){
#         tmpGG <- tmpGG + geom_smooth(method = lm, se=FALSE, fullrange= TRUE, aes(group=1))
#     }
#     return(tmpGG)
# }
#
# # ggboxplot with variables cut into ordered categories by Interval
# ggBoxPlot_cutInterval <- function(x, y, nCut_x){
#     tmpDF <- data.frame(var1 = x, var2 = y)
#     tmpDF$cut_x <- cut_interval(tmpDF$var1, nCut_x)
#     if(mean(is.na(tmpDF$cut_x > 0))){
#       levels(tmpDF$cut_x) <- c("NA", levels(tmpDF$cut_x))
#       tmpDF$cut_x[is.na(tmpDF$cut_x)] <- "NA"
#     }
#     ggplot(tmpDF, aes(y= var2, x = cut_x)) + geom_boxplot(outlier.colour = NA) + theme_classic() +
#         stat_n_text(size = 3, angle = 90) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# }
#
# # ggboxplot with variables cut into ordered categories by number of data points
# ggBoxPlot_cutNumber <- function(x, y, nCut_x){
#     tmpDF <- data.frame(var1 = x, var2 = y)
#     tmpDF$cut_x <- cut_number(tmpDF$var1, nCut_x)
#     ggplot(tmpDF, aes(y= var2, x = cut_x)) + geom_boxplot(outlier.colour = NA) + theme_classic() +
#         stat_n_text(size = 3, angle = 90) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
# }
#
