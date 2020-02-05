brary(dplyr)
library(ggplot2)
library(seqminer)
library(EnsDb.Hsapiens.v75)
library(ggbio)
library(patchwork)
library(stringr)
library(biomaRt)
library(dplyr)
library(Seurat)
options(stringsAsFactors=FALSE)

plot_region_pileups = function(chrom, # chromosome
                                start, # start coordinate
                                end, # end coordinate
                                bam_bed_files, # list of tabix indexed BED files with chrom,start,end,cell
                                grouping, # dataframe with cell,group,total_reads
                                ensdb=EnsDb.Hsapiens.v75, # EnsDb to allow display of genes
                                x_steps=500) { # limits number of actual to plot data for

  # Validate the grouping dataframe
  if (! all(c('cell', 'group', 'total_reads') %in% colnames(grouping))) {
    stop('Grouping DF must have cell, grouping, and total_reads columns.')
  }

  
  if (! is.vector(bam_bed_files)) {
      bam_bed_files = c(bam_bed_files)
  }

  # Get per-base-pair coverage over region for each cell
  per_base_coverage = lapply(bam_bed_files, function(bam_bed_file) {

    print(paste0('Processing: ', bam_bed_file))

    # Get reads for the requested region of BED file (using tabix)
    range_string = paste0(chrom, ':', start, '-', end)
    reads = seqminer::tabix.read.table(tabixRange=range_string, tabixFile=bam_bed_file)
    colnames(reads) = c('chrom', 'start', 'stop', 'cell')

    # Get count of cells per group for normalization
    grouping = grouping %>%
        group_by(group) %>%
        mutate(cells_in_group = n()) %>%
        ungroup()

      # Restrict to cell barcodes specified by user
    reads = dplyr::filter(reads, cell %in% grouping$cell)
    reads = dplyr::inner_join(reads, grouping, by='cell')

    # Expand reads so have one line per base
    per_base_coverage = dplyr::bind_rows(lapply(1:nrow(reads), function(i) {
        interval = reads[i, 'start']:reads[i, 'stop']
        return(data.frame(position=interval,
                        value=1,
                        cell=reads[i, 'cell'],
                        group=reads[i, 'group'],
                        scaling_factor=log10(reads[i, 'total_reads']),
                        cells_in_group=reads[i, 'cells_in_group']))
    }))

    return(per_base_coverage)
  })

  # Combine over all files
  per_base_coverage = dplyr::bind_rows(per_base_coverage)

  # Apply scaling factor and account for different cell counts
  # Also apply windowed smoothing
  per_base_coverage = per_base_coverage %>%
    group_by(position, group) %>%
    summarize(total=sum(value / scaling_factor) / mean(cells_in_group)) %>%
    group_by(group) %>%
    arrange(position) %>%
    mutate(total_smoothed=total) %>%
    #mutate(total_smoothed=zoo::rollapply(total, smoothing_window, mean, align='center',fill=NA)) %>%
    ungroup()

  # Finally, restrict to x_steps points for plotting to speed things up
  total_range = end - start
  stepsize = ceiling(total_range / x_steps)
  positions_to_show = seq(start, end, by=stepsize)
  per_base_coverage = subset(per_base_coverage, position %in% positions_to_show)

  # Make plots faceted by group
  max_value = max(per_base_coverage$total_smoothed, na.rm = TRUE)

  plot_object = ggplot() +
    geom_bar(data=per_base_coverage, aes(position, total_smoothed), stat='identity', color='black', fill='#d3d3d3', alpha=0.3) +
    facet_wrap(~group, ncol=1, strip.position="right") +
    ylim(0, max_value) +
    xlab('Position (bp)') +
    ylab(paste0('Normalized Signal (', smoothing_window, ' bp window)')) +
    theme_classic()

  # Add gene track using ensembl ID
  gr <- GRanges(seqnames = str_replace(chrom, 'chr', ''), IRanges(start, end), strand = "*")
  filters = AnnotationFilterList(GRangesFilter(gr), GenebiotypeFilter('protein_coding'))

  genes = autoplot(ensdb, filters, names.expr = "gene_name")

  genes_plot = genes@ggplot +
    xlim(start, end) +
    theme_classic()

  # Now combine the two plots into one (using patchwork)
  return(plot_object / (genes_plot + xlim(start, end)))
}

gene_body_coords = function(gene_name, mart) {

    coords = getBM(attributes=c("chromosome_name", "start_position", "end_position"),
        filters="hgnc_symbol", values=gene_name, mart=mart)

    if (nrow(coords) == 0) {
        stop(paste0('No coordinates found for gene: ', gene_name))
    }

    return(list('chrom'=paste0('chr', coords$chromosome_name[1]), 'start'=coords$start_position[1], 'end'=coords$end_position[1]))
}

seurat_to_grouping_table = function(seurat_obj, cluster_column=NULL) {
    grouping_table = seurat_obj@meta.data
    grouping_table$cell = rownames(grouping_table)

    if (!is.null(cluster_column)) {
        grouping_table$group = grouping_table[, cluster_column]
    } else {
        grouping_table$group = Idents(seurat_obj)
    }

    grouping_table$total_reads = grouping_table$nCount_RNA
    return(grouping_table %>% dplyr::select(cell, group, total_reads))
}
