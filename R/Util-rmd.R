options(stringsAsFactors = FALSE)
library(data.table)
library(tidyverse)

`%&%` <- function(a, b) paste0(a, b)

sigmoid <- function(x) 1/(1 + exp(-x))

logit <- function(x) log(x) - log(1 - x)

num.sci <- function(x) format(x, digits = 2, scientific = TRUE)

num.round <- function(x, d=2) round(x, digits = d)

num.int <- function(x) format(x, big.mark = ',')

.gsub.rm <- function(x, pat) gsub(x, pattern = pat, replacement = '')

.unlist <- function(...) unlist(..., use.names = FALSE)

.select <- function(...) select(...) %>% .unlist()

log.msg <- function(...) {
    ss = as.character(date())
    cat(sprintf('[%s] ', ss), sprintf(...), '\n', file = stderr(), sep = '')
}

################################################################

.cor.str <- function(x, y, ...) {
  .ret = cor.test(x, y, ...)
  num.round(.ret$estimate) %&%
    " (p=" %&% num.sci(.ret$p.value) %&% ")"
}


################################################################
.mkdir <- function(...) dir.create(..., recursive=TRUE, showWarnings=FALSE)

.gg.save <- function(filename, ..., cat.link = TRUE) {
    .mkdir(dirname(filename))
    if(file.exists(filename)) {
        log.msg('File already exits: %s', filename)
    } else {
        ggplot2::ggsave(..., filename = filename,
                        limitsize = FALSE,
                        units = 'in',
                        dpi = 300,
                        useDingbats = FALSE)
    }
    if(cat.link) {
        cat("\n[PDF](",filename,")\n\n",sep="")
    }
}

.gg.plot <- function(...) {
    ggplot2::ggplot(...) +
        ggplot2::theme_classic() +
        lemon::coord_capped_cart(left = 'both', bottom = 'both') +
        ggplot2::theme(plot.background = element_blank(),
                       plot.margin = unit(c(0,.5,0,.5), 'lines'),
                       plot.title = element_text(size = 10),
                       panel.background = element_blank(),
                       strip.background = element_blank(),
                       strip.text = element_text(size=4),
                       legend.background = element_blank(),
                       legend.text = element_text(size = 6),
                       legend.title = element_text(size = 6),
                       axis.title = element_text(size = 8),
                       legend.key.width = unit(.2, 'lines'),
                       legend.key.height = unit(.5, 'lines'),
                       legend.key.size = unit(.2, 'lines'),
                       axis.line = element_line(color = 'gray20', size = .2),
                       axis.ticks = element_line(color = 'gray20', size = .2),
                       axis.text = element_text(size = 4))
}

.remove <- function(...) gsub(..., replacement = '')

match.widths.grob <- function(g.list) {

    max.width = g.list[[1]]$widths[2:7]

    for(j in 2:length(g.list)) {
        max.width = grid::unit.pmax(max.width, g.list[[j]]$widths[2:7])
    }

    for(j in 1:length(g.list)) {
        g.list[[j]]$widths[2:7] = as.list(max.width)
    }
    return(g.list)
}

match.widths <- function(p.list) {
    g.list = lapply(p.list, ggplotGrob)
    return(match.widths.grob(g.list))
}

grid.vcat <- function(p.list, ...) {
    g.list = match.widths(p.list)
    ret = gridExtra::grid.arrange(grobs = g.list, ncol = 1, ...)
    return(ret)
}

match.heights.grob <- function(g.list, stretch = TRUE)  {
    max.height = g.list[[1]]$heights[2:7]

    if(stretch) {
        for(j in 2:length(g.list)) {
            max.height = grid::unit.pmax(max.height, g.list[[j]]$heights[2:7])
        }
    }

    for(j in 1:length(g.list)) {
        g.list[[j]]$heights[2:7] = as.list(max.height)
    }

    return(g.list)
}

match.heights <- function(p.list, stretch = FALSE) {
    g.list = lapply(p.list, ggplotGrob)
    return(match.heights.grob(g.list, stretch))
}

grid.hcat <- function(p.list, ...) {
    g.list = match.heights(p.list, stretch = TRUE)
    ret = gridExtra::grid.arrange(grobs = g.list, nrow = 1, ...)
    return(ret)
}

################################################################
match.2by2 <- function(p1, p2, p3, p4) {

    g1 = ggplotGrob(p1)
    g2 = ggplotGrob(p2)
    g3 = ggplotGrob(p3)
    g4 = ggplotGrob(p4)

    gg = match.widths.grob(list(g1, g3))
    g1 = gg[[1]]
    g3 = gg[[2]]

    gg = match.widths.grob(list(g2, g4))
    g2 = gg[[1]]
    g4 = gg[[2]]

    gg = match.heights.grob(list(g1, g2))
    g1 = gg[[1]]
    g2 = gg[[2]]

    gg = match.heights.grob(list(g3, g4))
    g3 = gg[[1]]
    g4 = gg[[2]]

    list(g1, g2, g3, g4)
}

################################################################
#' @param mat
row.order <- function(mat) {
    require(cba)
    require(proxy)

    if(nrow(mat) < 3) {
        return(1:nrow(mat))
    }

    D = suppressMessages(proxy::dist(mat, method <- function(a,b) 1 - cor(a,b, method = 'spearman')))
    D[!is.finite(D)] = 0
    h.out = hclust(D)
    o.out = cba::order.optimal(D, h.out$merge)
    return(o.out$order)
}

#' @param pair.ta
#' @param .ro
#' @param ret.tab
col.order <- function(pair.tab, .ro, ret.tab = FALSE) {

    M = pair.tab %>%
        dplyr::select(row, col, weight) %>%
        dplyr::mutate(row = factor(row, .ro)) %>%
        tidyr::spread(key = col, value = weight, fill = 0)

    co = order(apply(M[, -1], 2, which.max), decreasing = TRUE)
    .co = colnames(M)[-1][co]
    if(ret.tab) {
        ret = pair.tab %>%
            mutate(row = factor(row, .ro)) %>% 
            mutate(col = factor(col, .co))
    } else {
        ret = .co
    }
    return(ret)
}

#' @param pair.tab
#' @param ret.tab
order.pair <- function(pair.tab, ret.tab=FALSE) {

    require(tidyr)
    require(dplyr)
    
    .tab = pair.tab %>% dplyr::select(row, col, weight)

    M = .tab %>% tidyr::spread(key = col, value = weight, fill = 0)
    rr = M[, 1] %>% unlist(use.names = FALSE)
    cc = colnames(M)[-1] %>% unlist(use.names = FALSE)

    ## log.msg('Built the Mat: %d x %d', nrow(M), ncol(M))
    ro = row.order(M %>% dplyr::select(-row) %>% as.matrix())

    ## log.msg('Sort the rows: %d', length(ro))
    co = order(apply(M[ro, -1], 2, which.max), decreasing = TRUE)

    ## co = row.order(t(M %>% dplyr::select(-row) %>% as.matrix()))
    ## log.msg('Sort the columns: %d', length(co))

    if(ret.tab){
        ret = pair.tab %>%
            mutate(row = factor(row, rr[ro])) %>%
            mutate(col = factor(col, cc[co]))
    } else {
        ret = list(rows = rr[ro], cols = cc[co], M = M)
    }

    return(ret)
}

################################################################
#' Read gene set information
read.geneset <- function(ensg) {

    ANNOT.FILE = "../data/msigdb.txt.gz"

    annot.tab = fread(ANNOT.FILE)[, transcript_start := as.integer(transcript_start)]
    annot.tab = annot.tab[, transcript_end := as.integer(transcript_end)]
    annot.tab = annot.tab[ensembl_gene_id %in% ensg]

    gs.tab = annot.tab[, .(ensembl_gene_id, hgnc_symbol, gs_subcat, gs_name)][, val := 1]

    C.mat = dcast(gs.tab, ensembl_gene_id + hgnc_symbol ~ gs_name,
                  fun.aggregate = length, value.var = "val")

    list(C = C.mat, gs = gs.tab)
}

################################################################
#' Read gene ontology and coding genes
read.ontology <- function() {

    .temp.file <- ".ensembl.ontology.rdata"

    if(file.exists(.temp.file)) {
        load(.temp.file)

    } else {

        ensembl = biomaRt::useMart(biomart='ENSEMBL_MART_ENSEMBL',
                                   host='uswest.ensembl.org',
                                   path='/biomart/martservice',
                                   dataset='hsapiens_gene_ensembl')

        ensembl.hs = biomaRt::useDataset('hsapiens_gene_ensembl',mart=ensembl)

        .attr = c('ensembl_gene_id',
                  'hgnc_symbol',
                  'chromosome_name',
                  'transcription_start_site',
                  'transcript_start',
                  'transcript_end',
                  'description',
                  'percentage_gene_gc_content')

        .temp = biomaRt::getBM(attributes=.attr,
                               filters=c('biotype'),
                               values=c('protein_coding'),
                               mart=ensembl.hs,
                               useCache = FALSE)

        genes.desc = .temp %>%
            dplyr::select(hgnc_symbol, description) %>%
            unique()

        .temp <- as.data.table(.temp)

        coding.genes = .temp[,
                             .(hgnc_symbol = paste(unique(hgnc_symbol), collapse='|'),
                               transcription_start_site = mean(transcription_start_site),
                               transcript_start = min(transcript_start),
                               transcript_end = max(transcript_end)),
                             by = .(chromosome_name, ensembl_gene_id)]

        ensg.tot = unique(coding.genes$ensembl_gene_id)
        go.attr = c('ensembl_gene_id', 'go_id',
                    'name_1006', 'namespace_1003')

        genes.go = biomaRt::getBM(filters = 'ensembl_gene_id',
                                  values = ensg.tot,
                                  attributes = go.attr,
                                  mart = ensembl.hs,
                                  useCache = FALSE)
        .ensembl.ontology <- 
            list(coding = coding.genes,
                 go = genes.go,
                 desc = genes.desc)

        save(.ensembl.ontology, file=.temp.file)
    }

    return(.ensembl.ontology)
}
