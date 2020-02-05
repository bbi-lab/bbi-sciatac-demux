library(Matrix)
library(argparse)
library(Seurat)
library(stringr)

thisFile <- function() {
        cmdArgs <- commandArgs(trailingOnly = FALSE)
        needle <- "--file="
        match <- grep(needle, cmdArgs)
        if (length(match) > 0) {
                # Rscript
                return(normalizePath(sub(needle, "", cmdArgs[match])))
        } else {
                # 'source'd via R console
                return(normalizePath(sys.frames()[[1]]$ofile))
        }
}

source(paste0(dirname(thisFile()), '/', 'r_helper_functions/io_functions.R'))

parser = argparse::ArgumentParser(description='Script that takes a set RDS files containing sparse matrices and merges them together into one. Supports RDS/mtx input and output.')
parser$add_argument('--matrices', nargs='+', required=TRUE, help='List of RDS files, each of which contains a sparse matrix.')
parser$add_argument('--output', '-o', required=TRUE, help='Output RDS file for the combined matrix.')
args = parser$parse_args()

# Read in data
if (stringr::str_detect(args$matrices[[1]], '[.]rds$')) {
	matrices = lapply(args$matrices, readRDS)
} else if (stringr::str_detect(args$matrices[[1]], '[.]mtx')) {
	matrices = lapply(args$matrices, load_mtx_file)
}

# Merge
# NOTE that if Seurat ever removes RowMergeSparseMatrices, I also have it in schelpers
if (length(args$matrices) < 2) {
	merged_data = matrices[[1]]
} else if (length(matrices) == 2) {
	merged_data = Seurat:::RowMergeSparseMatrices(matrices[[1]], matrices[[2]])
} else {
	merged_data = matrices[[1]]
	for (i in 2:length(matrices)) {
		merged_data <<- Seurat:::RowMergeSparseMatrices(merged_data, matrices[[i]])
	}
}

print(paste0('Writing merged matrix with ', nrow(merged_data), ' rows and ', ncol(merged_data), ' columns.'))
if (stringr::str_detect(args$output, '[.]rds$')) {
	saveRDS(merged_data, file=args$output)
} else if (stringr::str_detect(args$output, '[.]mtx')) {
	write_mtx_file(merged_data, file=args$output)
}