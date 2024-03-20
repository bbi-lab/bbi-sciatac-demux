#!/usr/bin/env Rscript

#
# TODO
#   o  specify individual PCR wells (understand better)
#   o  bad well info
#   o  tag i7 switch
#   o  make summary table here? No.
#   o  build pcr_plate_list
#   o  build include_norm
#   o  well_df_out
#   o  barcode_list
#   o  build bad_well_list
#   o  pcr_combo_list
#   o  lig_combo_list
#   o  include_norm

# variable                       value source
# --------                       ------------
# lane_list                      Illumina run directory (barcode reader)
# plate_list                     sample sheet
# pcr combo list                 sample sheet
# include norm                   boolean sample sheet
# bad_wells_barcodes             barcode_list in generate_html.R from na in counts
# bad_wells                      bad_well_list in generate_html.R from na in counts
# pcr_well_info                  well-wise only  see bbi-dmux generate_html.R
#

# test of PCR row/col vs wells
# tst.R --p7_rows=0 --p5_cols=0 --p7_wells='D3 E3' --p5_wells='B2 B7'
# args$p7_rows:  0 
# args$p5_cols:  0 
# args$p7_wells:  D3 E3 
# args$p5_wells:  B2 B7 
# p7_rows:  D03 E03 
# p5_cols:  B02 B07 
# well_wise:  TRUE 
# well_df:
# 'data.frame':   2 obs. of  2 variables:
#  $ p5: chr  "B02" "B07"
#  $ p7: chr  "D03" "E03"
#
# tst.R --p7_rows="E F" --p5_cols="7 8"
# args$p7_rows:  E F 
# args$p5_cols:  7 8 
# args$p7_wells:  0 
# args$p5_wells:  0 
# p7_rows:  E F 
# p5_cols:  7 8 
# well_wise:  FALSE 


suppressPackageStartupMessages({
  library(ggplot2)
  library(argparse)
  library(jsonlite)
  library(shiny)
  library(stringr)
})


parser = argparse::ArgumentParser(description='Script to generate demux dashboard.')
parser$add_argument('--input_folder', required=TRUE, help='Input folder, typically <run_process_dir>/demux_out/fastqs_barcode.')
parser$add_argument('--samplesheet', required=TRUE, help='Samplesheet file name.')
parser$add_argument('--project_name', required=TRUE, help='Name of the project.')
parser$add_argument('--image_output_folder', required=TRUE, help='Path to image output folder.')
args = parser$parse_args()

# the following arguments are no longer used because the samplesheet
# JSON file has the information
# parser$add_argument('--p7_rows', required=TRUE, help='p7 rows')
# parser$add_argument('--p5_cols', required=TRUE, help='p5 cols')
# parser$add_argument('--p7_wells', required=TRUE, help='p7 wells')
# parser$add_argument('--p5_wells', required=TRUE, help='p5 wells')


# call ggplot to make u-titer plate plots
ggplot_data <- function(data=data, plot_type=c('tag','pcr','normalized','no_sample'), normalization_factor = NULL, sample_type=NULL)
{
  plot <- ggplot(data=data) + theme_bw() + labs(x = '', y = '')
      
  if(plot_type == 'tag') {
    plot <- plot + geom_point(mapping = aes(as.factor(Var1), Var2, fill = ReadCount), shape = 21, size = 10) +
            scale_fill_gradient(low = 'white', high = 'blue')
  } else
  if(plot_type == 'pcr') {
    plot <- plot + geom_point(mapping = aes(as.factor(Var1), Var2, fill = ReadCount, color = outlier, stroke = outlier), shape = 21, size = 10) +
            scale_fill_gradient(low = 'white', high = 'blue') +
            scale_color_manual(values = c('black', 'red'), guide = FALSE) +
            scale_discrete_manual(aesthetics = 'stroke', values = c(0.5,1), guide = FALSE)
  } else
  if(plot_type == 'normalized') {
    plot <- plot + geom_point(mapping = aes(as.factor(Var1), Var2, fill = log2(ReadCount/normalization_factor)), shape = 21, size = 10) +
            scale_fill_gradient2(name = 'log2 norm value', low = 'red', mid = 'white', high = 'blue')
  } else
  if(plot_type == 'no_sample') {
    plot <- plot + geom_point(aes(as.factor(Var1), Var2), shape = 21, size = 10) +
        geom_text(aes(x = 6.5, y = 'D', label = paste0('No ', sample_type, ' detected for this plate')))
  }
  return(plot)
}


#
# Infer run lanes from csv files.
# Return a list of lanes, e.g., ('L001', ...).
#
get_run_lanes <- function(input_folder)
{
  file_list <- list.files(input_folder, pattern = 'index_counts.csv$')
  grxp <- '(RUN001_|.index_counts.csv)'
  lane_list <- gsub(grxp, '', file_list, perl = TRUE)
  return(lane_list)
}


# read 'index' file of read counts by tag/pcr well
read_index_file <- function(index_filename)
{
  index_data <- read.csv(index_filename, header = TRUE, stringsAsFactors = FALSE)
  return(index_data)
}


# make tagi7 or tagi5 tagmentation plots
make_tag_plots <- function(lane, tag_data, sent_wells, barn_wells, tag_name = c('tagi7', 'tagi5'))
{
  # set up u-titer plate data frame
  rows <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
  cols <- 1:12
  plate <- data.frame(expand.grid(cols, rows))
  levels(plate$Var2) <- rev(levels(plate$Var2))
 
  # parse out plate and well strings 
  tag_data$plate <- stringr::str_split_fixed(tag_data$plate_well, '-', 2)[,1]
  tag_data$well  <- stringr::str_split_fixed(tag_data$plate_well, '-', 2)[,2]
  tag_data$rows <- substring(tag_data$well, first = 1, last = 1)
  tag_data$cols <- suppressWarnings(as.numeric(substring(tag_data$well, first = 2, last = nchar(tag_data$well))))
 
  # identify bad wells, which are wells with (fill in here) 
  bad_well <- tag_data[is.na(tag_data$rows) | is.na(tag_data$cols) | !tag_data$rows %in% rows | !tag_data$cols %in% cols,]

  # identify wells with valid data (?)
  # Note: I don't understand the is.na(tag_data$cols).
  tag_data <- tag_data[!(is.na(tag_data$rows) | is.na(tag_data$cols) | !tag_data$rows %in% rows | !tag_data$cols %in% cols),]

  # obvious, no?
  for(p in unique(tag_data$plate)) {
    # choose data for plate 'p'
    plate_df <- tag_data[tag_data$plate == p,]

    # calculate plate sentinel normalization factor, if relevant    
    sent_norm <- NULL
    if (!is.null(sent_wells)) {
      if (any(sent_wells %in% tag_data$plate_well)) {
        sent_norm <- mean(tag_data[tag_data$plate_well %in% sent_wells,]$ReadCount)
      }
    }
    
    # calculate plate barnyard normalization factor, if relevant    
    barn_norm <- NULL
    if (!is.null(barn_wells)) {
      if (any(barn_wells %in% tag_data$plate_well)) {
        barn_norm <- mean(tag_data[tag_data$plate_well %in% barn_wells,]$ReadCount)
      }
    }
   
    # Map counts to plates/wells.
    data <- merge(plate, 
        plate_df[,c('rows', 'cols', 'ReadCount')],
        by.x=c('Var1', 'Var2'),
        by.y=c('cols', 'rows'), all.x=T)
    data$ReadCount[is.na(data$ReadCount)] <- 0
  
#   diagnostic
    fn <- paste0('lig_',tag_name,'_lane_',lane,'_plate_',p,'_data_frame.txt')
    write.table(data,file=fn)
     
    # make plate plot with sentinel normalization (or none) 
    if (!is.null(sent_norm)) {
#      file_name <- paste0(args$image_output_folder,'/', lane, '_', p, '.', tag_name, '_plate_sent_norm.png')
#      png(file = file_name, width = 6, height = 4, res = 200, units = 'in')
#      print(ggplot_data(data = data, plot_type = 'normalized', normalization_factor = sent_norm))
#      dev.off()

      file_name <- paste0(args$image_output_folder,'/', lane, '_', p, '.', tag_name, '_plate_sent_norm.png')
      ggp_obj <- ggplot_data(data = data, plot_type = 'normalized', normalization_factor = sent_norm)
      ggsave(filename=file_name, plot=ggp_obj, device='png', width=6, height=4, dpi=200, units='in')

    } else {
#      file_name <- paste0(args$image_output_folder,'/', lane, '_', p, '.', tag_name, '_plate_sent_norm.png')
#      png(file = file_name, width = 5.5, height = 4, res = 200, units = 'in')
#      print(ggplot_data(data = data, plot_type = 'no_sample', sample_type = 'sentinel'))
#      dev.off()

      file_name <- paste0(args$image_output_folder,'/', lane, '_', p, '.', tag_name, '_plate_sent_norm.png')
      ggp_obj <- ggplot_data(data = data, plot_type = 'no_sample', sample_type = 'sentinel')
      ggsave(filename=file_name, plot=ggp_obj, device='png', width=5.5, height=4, dpi=200, units='in')

    }    
    
    # make plate plot with sentinel normalization (or none) 
    if (!is.null(barn_norm)) {
#      file_name <- paste0(args$image_output_folder,'/', lane, '_', p, '.', tag_name, '_plate_barn_norm.png')
#      png(file = file_name, width = 6, height = 4, res = 200, units = 'in')
#      print(ggplot_data(data = data, plot_type = 'normalized', normalization_factor = barn_norm))
#      dev.off()

      file_name <- paste0(args$image_output_folder,'/', lane, '_', p, '.', tag_name, '_plate_barn_norm.png')
      ggp_obj <- ggplot_data(data = data, plot_type = 'normalized', normalization_factor = barn_norm)
      ggsave(filename=file_name, plot=ggp_obj, device='png', width=6, height=4, dpi=200, units='in')

    } else {
#      file_name <- paste0(args$image_output_folder,'/', lane, '_', p, '.', tag_name, '_plate_barn_norm.png')
#      png(file = file_name, width = 5.5, height = 4, res = 200, units = 'in')
#      print(ggplot_data(data = data, plot_type = 'no_sample', sample_type = 'barnyard'))
#      dev.off()

      file_name <- paste0(args$image_output_folder,'/', lane, '_', p, '.', tag_name, '_plate_barn_norm.png')
      ggp_obj <- ggplot_data(data = data, plot_type = 'no_sample', sample_type = 'barnyard')
      ggsave(filename=file_name, plot=ggp_obj, device='png', width=5.5, height=4, dpi=200, units='in')

    }

    # make plate plot (no normalization)
#    file_name <- paste0(args$image_output_folder,'/', lane, '_', p, '.', tag_name, '_plate.png')
#    png(file = file_name, width = 6, height = 4, res = 200, units = 'in')
#    print(ggplot_data(data = data, plot_type = 'tag'))
#    dev.off()

    file_name <- paste0(args$image_output_folder,'/', lane, '_', p, '.', tag_name, '_plate.png')
    ggp_obj <- ggplot_data(data = data, plot_type = 'tag')
    ggsave(filename=file_name, plot=ggp_obj, device='png', width=6, height=4, dpi=200, units='in')

  }
}


# make i7 and i5 tagmentation plot sets
make_tag_plot_sets <- function(args, lane_list)
{
  for (lane in lane_list) {
    # read lane index file
    file_name <- paste0(args$input_folder, '/RUN001_', lane, '.index_counts.csv')
    index_data <- read.csv(file_name, header = TRUE, stringsAsFactors = FALSE)
    
    # identify wells with sentinel or barnyard samples
    if(length(index_data$sample_name_tagi5 == 'Sentinel') > 0) {
      sent_wells <- as.character(index_data[index_data$sample_name_tagi5 == 'Sentinel',]$i5_well)
    } else {
      sent_wells <- NULL
    }
    if(length(index_data$sample_name_tagi5 == 'Barnyard') > 0) {
      barn_wells <- as.character(index_data[index_data$sample_name_tagi5 == 'Barnyard',]$i5_well)
    } else {
      barn_wells <- NULL
    }
    
    # tagi5 plots
    tag_data <- index_data[,c('i5_well', 'tagi5_count')]
    colnames(tag_data) <- c('plate_well', 'ReadCount')
    make_tag_plots(lane, tag_data, sent_wells, barn_wells, 'tag_i5')

    # tagi7 plots
    tag_data <- index_data[,c('i7_well', 'tagi7_count')]
    colnames(tag_data) <- c('plate_well', 'ReadCount')
    make_tag_plots(lane, tag_data, sent_wells, barn_wells, 'tag_i7') 
  }
}


make_pcr_plot_sets <- function(args, lane_list) {
  # get PCR plate information from samplesheet json file 
  sample_data <- read_json(args$samplesheet)
  pcr_format <- sample_data$pcr_format

  # We can handle only row-column PCR format
  # at this moment.
  if( pcr_format != 'row_col' )
    return( 0 )

  # set up u-titer plate data frame
  rows <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H')
  cols <- 1:12
  plate <- data.frame(expand.grid(cols, rows))
  levels(plate$Var2) <- rev(levels(plate$Var2))
  well_fix <- list('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12')

  # get PCR plate information from samplesheet json file 
  sample_data <- read_json(args$samplesheet)
  pcr_format <- sample_data$pcr_format
  p7_rows <- unlist(sample_data$p7_row_list)
  p5_cols <- unlist(sample_data$p5_col_list)
  well_wise <- FALSE

  # make sure that p7 row specs are upper-case
  for( irow in seq(length(p7_rows))) {
    p7_rows[irow] <- toupper(p7_rows[irow])
  }

  # the commented out code below is mostly obsolete because the
  # information is in the JSON file
  # set up PCR plate rows and columns used
  # also...I have not attempted to understand the meaning of p7_rows[1] == 'none' or
  # p5_cols[1] == 'none' so the code below is unlikely to function correctly in those
  # cases.
#  if (args$p7_rows != 0) {
#    if (args$p7_rows != 'none') {
#      p7_rows <- unlist(stringr::str_split(args$p7_rows, ' '))
#    } else {
#      p7_rows <- 'none'
#    } 
#    if (args$p5_cols != 'none') {
#      p5_cols <- unlist(stringr::str_split(args$p5_cols, ' '))
#    } else {
#      p5_cols <- 'none'
#    }
#    well_wise <- FALSE
#  } else {
#    if (args$p7_wells != 'none') {
#      p7_rows <- unlist(stringr::str_split(args$p7_wells, ' '))
#      p7_rows <- paste0(substr(p7_rows, start=1, stop=1), well_fix[as.numeric(substr(p7_rows, start = 2, stop=3))])
#    } else {
#      p7_rows <- 'none'
#    }
#    if (args$p5_wells != 'none') {
#      p5_cols <- unlist(stringr::str_split(args$p5_wells, ' '))
#      p5_cols <- paste0(substr(p5_cols, start=1, stop=1), well_fix[as.numeric(substr(p5_cols, start = 2, stop=3))])
#    } else {
#      p5_cols <- 'none'
#    }
#    well_df <- data.frame(p5 = p5_cols, p7 = p7_rows)
#    well_wise <- TRUE
#  }

  for (lane in lane_list) {
    # read lane data file
    file_name <- paste0(args$input_folder, '/RUN001_', lane, '.pcr_pair_counts.csv')
    pcr_data <- read.csv(file = file_name, header = TRUE, stringsAsFactors = FALSE)

    # parse out i7 and i5 well information
    # (edit out plate specification: there is only one PCR plate)
    pcr_data <- pcr_data[,c('pcri5_well','pcri7_well','pcr_pair_count')]
    pcr_data[,'pcri5_well'] <- sub('P1-', '', pcr_data[,'pcri5_well'])
    pcr_data[,'pcri7_well'] <- sub('P1-', '', pcr_data[,'pcri7_well'])
    colnames(pcr_data) <- c('V1','V2','V3')

    # parse out i5 well row and column strings    
    if(p5_cols[1] != 'none') {
      pcr_data$p5_row <- substring(pcr_data$V1, first = 1, last = 1)
      pcr_data$p5_col <- as.numeric(substring(pcr_data$V1, first = 2,
              last = nchar(pcr_data$V1)))
    }
    
    # parse out i7 well row and column strings    
    if (p7_rows[1] != 'none') {
      pcr_data$p7_row <- substring(pcr_data$V2, first = 1, last = 1)
      pcr_data$p7_col <- as.numeric(substring(pcr_data$V2, first = 2,
              last = nchar(pcr_data$V2)))
    }

    # ?    
    if (p7_rows[1] == 'none' || p5_cols[1] == 'none') {
      pcr_plate_list <- c()
      rel_barc <- ifelse(p7_rows == 'none', 'p5', 'p7')
      if(p7_rows[1] == 'none') {
        pcr_data$rows <- substring(pcr_data$V1, first = 1, last = 1)
        pcr_data$cols <- as.numeric(substring(pcr_data$V1, first = 2,
                last = nchar(pcr_data$V1)))
      } else {
        pcr_data$rows <- substring(pcr_data$V2, first = 1, last = 1)
        pcr_data$cols <- as.numeric(substring(pcr_data$V2, first = 2,
                last = nchar(pcr_data$V2)))
      }
      p5_df <- data.frame(rows = pcr_data$rows,
          cols = pcr_data$cols, ReadCount = pcr_data$V3)
      data <- merge(plate, p5_df, by.x=c('Var1', 'Var2'),
          by.y=c('cols', 'rows'), all.x=T)
      data$ReadCount[is.na(data$ReadCount)] <- 0
     
      # identify outlier wells (read count < 0.01 of median count) 
      data$outlier <- data$ReadCount < 0.01 * median(data$ReadCount)

      pcr_plate_list <- c(pcr_plate_list, paste0(rel_barc))
#      file_name <- paste0(args$image_output_folder,'/', lane,'_',rel_barc, '.pcr_plate.png')
#      png(file = paste0(args$image_output_folder,'/', lane,'_',rel_barc, '.pcr_plate.png'), width = 6, height = 4, res = 200, units = 'in')
#      print(ggplot_data(data = data, plot_type = 'pcr'))
#      dev.off()

      file_name <- paste0(args$image_output_folder,'/', lane,'_',rel_barc, '.pcr_plate.png')
      ggp_obj <- ggplot_data(data = data, plot_type = 'pcr')
      ggsave(filename=file_name, plot=ggp_obj, device='png', width=6, height=4, dpi=200, units='in')

      well_df_lane <- data.frame()
    } else if (!well_wise) {
      pcr_plate_list <- c()
      for (i in  1:length(p7_rows)) {
        sub <- subset(pcr_data, p7_row == p7_rows[i] & p5_col == as.numeric(p5_cols[i]))
        p5_df <- data.frame(rows = sub$p5_row,
                            cols = sub$p7_col,
                            ReadCount = sub$V3)
        data <- merge(plate, p5_df,
                      by.x=c('Var1', 'Var2'),
                      by.y=c('cols', 'rows'),
                      all.x=T)

#        diagnostic
        fn <- paste0('pcr_lane_',lane,'_p7_row_',p7_rows[i],'_p5_col_',p5_cols[i],'_data_frame.txt')
        write.table(data,file=fn)
     
        data$ReadCount[is.na(data$ReadCount)] <- 0
        
        data$outlier <- data$ReadCount < 0.01 * median(data$ReadCount)
        pcr_plate_list <- c(pcr_plate_list, paste0(p7_rows[i], p5_cols[i]))
#        file_name <- paste0(args$image_output_folder,'/', lane,'_', p7_rows[i], p5_cols[i], '.pcr_plate.png')
#        png(file = paste0(args$image_output_folder,'/', lane,'_', p7_rows[i], p5_cols[i], '.pcr_plate.png'), width = 6, height = 4, res = 200, units = 'in')
#        print(ggplot_data(data = data, plot_type = 'pcr'))
#        dev.off()

        file_name <- paste0(args$image_output_folder,'/', lane,'_', p7_rows[i], p5_cols[i], '.pcr_plate.png')
        ggp_obj <- ggplot_data(data = data, plot_type = 'pcr')
        ggsave(filename=file_name, plot=ggp_obj, device='png', width=6, height=4, dpi=200, units='in')

        well_df_lane <- data.frame()
      }
      
    } else {
      pcr_plate_list <- c()
      well_df_lane <- merge(well_df, pcr_data[,c('V1', 'V2', 'V3')], by.x=c('p5', 'p7'), by.y=c('V1', 'V2'))
      names(well_df_lane)[3] <- gsub('L00', 'Lane ', lane)
    }
  }
}


make_demux_image_dir <- function(args) {
  output_dir <- args$image_output_folder
  if (!dir.exists(output_dir)) {
    dir.create(output_dir,recursive=TRUE)
  }
}


make_plots <- function(args) {
  lane_list <- get_run_lanes(args$input_folder)
  lane_names <- gsub('L00', 'Lane ', lane_list)
  lane_nums <- gsub('L00', '', lane_list)

  make_demux_image_dir(args)  
  make_tag_plot_sets(args, lane_list)
  make_pcr_plot_sets(args, lane_list)
}

make_plots(args)
