library(dplyr)
library(tidyr)
library(ggplot2)
library(grid)
library(MASS)

# Manhattan plot of SNVs with their corresonding (-log) p-values 
# associated with gene expression.
# results, snps, chooseGene, gtf, width
plot_manhattan_gene <- function(res, snps, whichGene, gtf, w = 2){
  
    if(whichGene==''){
        whichGene = res$gene_name[1] 
    }
    manh = res[res$gene_name==whichGene, ]
    manh = manh %>% left_join(snps, by='snpid')
    manh = manh [,c("chr", "pos", "pValue", "snpid")]
  
    geneForPlot = gtf[gtf$gene_name==whichGene, 
                      c('gene_name','seqid','start','end')] %>% 
                  group_by(gene_name) %>% 
                  summarise(start = min(start), 
                            end = max(end), 
                            chr = dplyr::first(seqid))
  
    w = (geneForPlot$end-geneForPlot$start) * w  # Plot width
  # Take gene neighboors (w form start/end)
    g = gtf %>% filter(seqid == geneForPlot$chr) %>% 
        dplyr::select(gene_name, gene_id, type, seqid, start, end, strand) %>% 
        filter(end > max(0,(geneForPlot$start - w))) %>% 
        filter(start < (geneForPlot$end + w)) 
  
    exons = g[g$type == 'exon',]
    g = g[g$type == 'gene',]
  
  # Find overlapping genes and sort them using count
    gg <- GRanges(seqnames = g$seqid,
                  ranges = IRanges(start = g$start, end=g$end), 
                  strand = g$strand, 
                  gene_name = g$gene_name)
  
    overlapps = as.data.frame(findOverlaps(gg,gg,ignore.strand=T))
    gg$count =  overlapps %>% 
                filter(queryHits >= subjectHits) %>%
                group_by(queryHits) %>%
                summarise(count=n()) %>%
                dplyr::select(count)
    
    g = g %>% left_join(data.frame(gene_name = gg$gene_name,
                                   count = gg$count), 
                        by='gene_name') # add count to genes
    g = exons %>% 
        right_join(g[,c('gene_name','count')], by='gene_name') %>% 
        dplyr::union(g)  # add gene count to exon and then make union of genes and exons
  
    maxP = max(-log10(manh$pValue)) # size of Y axis
    gls = as.integer(maxP * 0.03) + 1  # size of genes, exons below the points
  
    igv = data.frame(gene_name = g$gene_name, 
                     type = g$type, 
                     xmin = g$start, 
                     xmax = g$end, 
                     strand = g$strand, 
                     count = as.integer(g$count))
    igv = igv[order(igv$type),] # to plot exons over
    igv$where <- 1
    igv$ymax = - gls*igv$where - gls * igv$count
    igv$ymin = - gls*(igv$where+1) - gls * igv$count
    igv$xmin = pmax(igv$xmin, geneForPlot$start - w)
    igv$xmax = pmin(igv$xmax, geneForPlot$end + w)
    igvText = data.frame(gene_name = igv$gene_name, 
                         xpos = (igv$xmin+igv$xmax)/2, 
                         ymin=igv$ymin, 
                         ymax=igv$ymax, 
                         count = igv$count, 
                         type=igv$type)
    igvText = igvText %>%
                filter(type=='gene') %>% 
                group_by(count) %>% 
                mutate(sync = (cumsum(ymin)/dplyr::first(ymin)) %% 4 + 1)
    
    igvText$ypos = (igvText$ymin * (igvText$sync) + 
                    igvText$ymax * (7 - igvText$sync))/7
  
  #### Create plots
  
  p1 <- ggplot(manh, aes(x=pos, y=-log10(pValue))) +
    ggtitle(paste(whichGene, "associated SNVs")) +
    geom_point(color="steelblue", fill="steelblue", size=3) +
    geom_point(data = manh[which.min(manh$pValue),], 
               aes(x=pos, y=-log10(pValue)),
               color="orange", 
               size=5) +
    scale_y_continuous(breaks = seq(from=0,
                                    to=as.integer(maxP)+1, 
                                    by = as.integer(maxP/5)), 
                       limits = c(0, as.integer(maxP)+5)) +
    scale_x_continuous(limits = c(max(0,geneForPlot$start - w), 
                                  geneForPlot$end + w)) +
    ylab("-log(pValue)") +
    theme_classic() +
    theme(plot.title = element_text(size=20, hjust = 0.5),
          text=element_text(size=17),
          axis.text=element_text(size=12),
          axis.title.x=element_blank(),
          axis.title.y=element_text(vjust=2.5, hjust = 0.5, size=17),
          axis.text.y=element_text(size=12),
          plot.margin=unit(c(5.5,5.5,0,5.5), "pt"))
  
  p2 <- ggplot() + geom_rect(data = igv[igv$type=='gene',],
                             aes(xmax = xmax, 
                                 xmin = xmin, 
                                 ymax = ymax, 
                                 ymin = ymin, 
                                 fill = type),
                             inherit.aes = FALSE, 
                             show.legend = TRUE) +
    geom_rect(data = igv[igv$type=='exon',], 
              aes(xmax = xmax, 
                  xmin = xmin, 
                  ymax = ymax, 
                  ymin = ymin, 
                  fill = type), 
              inherit.aes = FALSE, 
              show.legend = TRUE) +
    geom_text(data = igvText, 
              aes(x=xpos, y=ypos, label=gene_name), 
              size = 3, 
              inherit.aes = TRUE, 
              check_overlap = TRUE) +
    geom_segment(data = filter(igv, type=='gene', strand=='+'), 
                 mapping = aes(x = xmin, 
                               y = (ymax+6*ymin)/7, 
                               xend = xmax, 
                               yend = (ymax+6*ymin)/7), 
                 arrow = arrow(length = unit(0.2, "cm")), 
                 size = 0.5, 
                 color = "blue", 
                 inherit.aes = T) +
    geom_segment(data = filter(igv, type=='gene', strand=='-'), 
                 mapping = aes(xend = xmin, 
                               yend = (ymax+6*ymin)/7, 
                               x = xmax, 
                               y = (ymax+6*ymin)/7), 
                 arrow = arrow(length = unit(0.2, "cm")), 
                 size = 0.5, 
                 color = "blue", 
                 inherit.aes = T) +
    scale_x_continuous(limits = c(max(0,geneForPlot$start - w), 
                                  geneForPlot$end + w)) +
    xlab("Genomic position (bp)") +
    theme(panel.border = element_rect(colour = "black", fill=NA),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          size = 0.5, linetype = "solid"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.y=element_text(vjust = 2, 
                                    hjust = 0.7, 
                                    size = 9, 
                                    color = 'white'),
          axis.text.y=element_text(size = 9, 
                                   color = 'white'),
          axis.ticks.y=element_blank(),
          legend.position="bottom",
          plot.margin=unit(c(5.5,5.5,5.5,5.5), "pt"),
          axis.title.x = element_text(vjust = 0, 
                                      hjust = 0.7, 
                                      size = 17),
          axis.text.x = element_text(size = 12),
          legend.title=element_text(size = 12),
          legend.text=element_text(size = 12))
  
  maxOverlaps = min(max(igv$count),6)   
  pp <- rbind(ggplotGrob(p1), ggplotGrob(p2), size = "last")
  idpanels <- unique(pp$layout[grepl("panel",pp$layout$name), "t"])
  pp$heights[idpanels][1] = unit(x=6,units = 'null')
  pp$heights[idpanels][2] = unit(x=maxOverlaps,units = 'null')
  grid.draw(pp)
}
