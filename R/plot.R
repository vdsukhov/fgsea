#' Plots table of enrichment graphs using ggplot and gridExtra
#' @param pathways Pathways to plot table, as in `fgsea` function
#' @param stats Gene-level stats, as in `fgsea` function
#' @param fgseaRes Table with fgsea results
#' @param gseaParam GSEA-like parameter. Adjusts displayed statistic values,
#'      values closer to 0 flatten plots. Default = 1, value of 0.5 is a good
#'      choice too.
#' @param colwidths Vector of five elements corresponding to column width for
#'      grid.arrange
#' @import ggplot2
#' @import gridExtra
#' @import grid
#' @export
#' @examples
#' data(examplePathways)
#' data(exampleRanks)
#' \dontrun{
#' fgseaRes <- fgsea(examplePathways, exampleRanks, nperm=1000,
#'                   minSize=15, maxSize=100)
#' topPathways <- fgseaRes[head(order(pval), n=15)][order(NES), pathway]
#' plotGseaTable(examplePathways[topPathways], exampleRanks,
#'               fgseaRes, gseaParam=0.5)
#' }
plotGseaTable <- function(pathways, stats, fgseaRes,
                          gseaParam=1,
                          colwidths=c(5, 3, 0.8, 1.2, 1.2)) {

    rnk <- rank(-stats)
    ord <- order(rnk)

    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj) ^ gseaParam)
    statsAdj <- statsAdj / max(abs(statsAdj))

    pathways <- lapply(pathways, function(p) {
        if (is.character(p)) {
            unname(as.vector(na.omit(match(p, names(statsAdj)))))
        } else {
            rnk[p]
        }
    })
    tt <- melt(pathways, as.is=TRUE)
    colnames(tt) <- c("rank", "pathway")

    ps <- lapply(names(pathways), function(pn) {
        p <- pathways[[pn]]
        annotation <- fgseaRes[match(pn, fgseaRes$pathway), ]
        list(
            textGrob(pn, just="right", x=unit(0.95, "npc")),
            ggplot() +
                geom_segment(aes(x=p, xend=p,
                                 y=0, yend=statsAdj[p])) +
                scale_x_continuous(limits=c(0, length(statsAdj)),
                                   expand=c(0, 0)) +
                scale_y_continuous(limits=c(-1, 1),
                                   expand=c(0, 0)) +
                xlab(NULL) + ylab(NULL) +
                theme(panel.background = element_blank(),
                      axis.line=element_blank(),
                      axis.text=element_blank(),
                      axis.ticks=element_blank(),
                      panel.grid = element_blank(),
                      axis.title=element_blank(),
                      axis.line = element_blank(),
                      panel.margin = element_blank(),
                      plot.margin = rep(unit(0,"null"),4),
                      panel.margin = rep(unit(0,"null"),4)
                ),
            textGrob(sprintf("%.2f", annotation$NES)),
            textGrob(sprintf("%.1e", annotation$pval)),
            textGrob(sprintf("%.1e", annotation$padj))
            )
    })


    rankPlot <-
        ggplot() +
        geom_blank() +
        scale_x_continuous(limits=c(0, length(statsAdj)),
                           expand=c(0, 0)) +
        scale_y_continuous(limits=c(-1, 1),
                           expand=c(0, 0)) +
        xlab(NULL) + ylab(NULL) +
        theme(panel.background = element_blank(),
              axis.line=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid = element_blank(),
              axis.title=element_blank(),
              axis.line = element_blank(),
              panel.margin = element_blank(),
              plot.margin = unit(c(0,0,0.5,0), "npc"),
              panel.margin = unit(c(0,0,0,0), "npc")
        )

    grid.arrange(grobs=c(
            lapply(c("Pathway", "Gene ranks", "NES", "pval", "padj"), textGrob),
            unlist(ps, recursive = F),
            list(nullGrob(),
                 rankPlot
                 )),
        ncol=5, widths=colwidths)
}