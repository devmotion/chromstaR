

#' Read bed-file into GRanges
#'
#' This is a simple convenience function to read a bed(.gz)-file into a \code{\link{GRanges}} object. The bed-file is expected to have at least the following three fields: \code{chromosome, start, end}.
#'
#' @param bedfile Filename of the bed or bed.gz file.
#' @param skip Number of lines to skip at the beginning.
#' @return A \code{\link{GRanges}} object with the contents of the bed-file.
#' @author Aaron Taudt, David Widmann
#' @importFrom utils count.fields read.table
#' @export
#'
#'@examples
#'## Get an example BED file
#'bedfile <- system.file("extdata", "liver-H3K4me3-BN-male-bio2-tech1.bed.gz",
#'                        package="chromstaRData")
#'## Import the file and skip the first 10 lines
#'data <- readBedFile(bedfile, skip=10)
#'
readBedFile <- function(bedfile, skip=0) {
    ncols <- max(utils::count.fields(bedfile, skip=skip))
    if ( ncols < 3 )
        stop("Not a correct BED3 file format.")

    if ( ncols < 6 ) {
        data <- utils::read.table(bedfile, colClasses=c("character", rep("numeric", 2), rep("NULL", ncols-3)), skip=skip)
        names(data) <- c("chrom", "chromStart", "chromEnd")
    } else {
        data <- utils::read.table(bedfile, colClasses=c("character", rep("numeric", 2),
                                                        rep("NULL", 2), "character",
                                                        rep("NULL", ncols-6)), skip=skip)
        names(data) <- c("chrom", "chromStart", "chromEnd", "strand")

        ## adjust strand information
        data$strand <- sub("^\\.$", "*", data$strand)
    }

    ## convert to GRanges object
    gr <- GenomicRanges::GRanges(seqnames=data$chrom,
                                 ranges=IRanges(start=data$chromStart+1,     # Convert from 0-based half open to 1-based closed
                                                end=data$chromEnd),
                                 strand=data$strand)

    return(gr)
}

#' Write GRanges to bed3-file
#'
#' This is a simple convenience function to write a \code{\link{GRanges}} object to a bed(.gz)-file. The bed-file will have the following three fields: \code{chromosome, start, end}.
#'
#' @param gr A \code{\link{GRanges}} object.
#' @param filename The name of the file that will be written. The ending ".bed.gz". Any existing file will be overwritten.
#' @param append Whether or not to append to an existing file.
#' @param chromosome.format Desired format of the chromosomes. Either 'NCBI' for (1,2,3 ...) or 'UCSC' for (chr1,chr2,chr3 ...).
#' @return \code{NULL}
#' @importFrom utils write.table
#' @author Aaron Taudt, David Widmann
#' @export
exportBed3File <- function(gr, filename, append=FALSE, chromosome.format='NCBI') {
    filename <- paste0(filename,".bed.gz")
    if (append) {
        filename.gz <- gzfile(filename, 'a')
        message('Appending to file ', filename)
    } else {
        filename.gz <- gzfile(filename, 'w')
        message('Writing to file ', filename)

        # Write first line to file
        cat("", file=filename.gz)
    }

    # Collect genomic data
    data <- data.frame("chrom" = as.factor(GenomeInfoDb::mapSeqlevels(seqnames=as.character(seqnames(gr)), style=chromosome.format)),
                       "chromStart" = start(gr) - 1, # Convert from 1-based closed to 0-based half open
                       "chromEnd" = end(gr))

    utils::write.table(format(data, scientific=FALSE, trim=TRUE), file=filename.gz,
                       append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

    close(filename.gz)
}
