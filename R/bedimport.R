

#' Read bed-file into GRanges
#'
#' This is a simple convenience function to read a bed(.gz)-file into a \code{\link{GRanges}} object. The bed-file is expected to have the following fields: \code{chromosome, start, end, name, score, strand}.
#'
#' @param bedfile Filename of the bed or bed.gz file.
#' @param skip Number of lines to skip at the beginning.
#' @param chromosome.format Desired format of the chromosomes. Either 'NCBI' for (1,2,3 ...) or 'UCSC' for (chr1,chr2,chr3 ...).
#' @return A \code{\link{GRanges}} object with the contents of the bed-file.
#' @author Aaron Taudt
#' @importFrom utils read.table
#' @export
#'
#'@examples
#'## Get an example BED file
#'bedfile <- system.file("extdata", "liver-H3K4me3-BN-male-bio2-tech1.bed.gz",
#'                        package="chromstaRData")
#'## Import the file and skip the first 10 lines
#'data <- readCustomBedFile(bedfile, skip=10)
#'
readCustomBedFile <- function(bedfile, skip=0, chromosome.format='NCBI') {

    # File with reads, determine classes first for faster import (0-based)
    classes <- c('character','numeric','numeric','character','integer','character')
    data <- utils::read.table(bedfile, colClasses=classes, skip=skip)
    # GRanges compatible strand information
    data[,6] <- sub('.','*',data[,6])
    # Adjust chromosome format
    data[,1] <- sub('^chr', '', data[,1])
    if (chromosome.format=='UCSC') {
        data[,1] <- paste0('chr', data[,1])
    }
    # Convert to GRanges object
    gr <- GenomicRanges::GRanges(seqnames=data[,1],
                                ranges=IRanges(start=data[,2]+1, end=data[,3]),     # Convert from 0-based half open to 1-based closed
                                strand=data[,6],
                                name=data[,4],
                                score=data[,5])

    return(gr)

}

readBed3File <- function(bedfile, skip=0, chromosome.format='NCBI') {
    ncols <- max(count.fields(bedfile, skip=skip))
    if ( ncols < 3 )
        stop("Not a correct BED3 file format.")

    data <- utils::read.table(bedfile, colClasses=c("character", rep("numeric", 2), rep("NULL", ncols-3)), skip=skip)
    names(data) <- c("chrom", "chromStart", "chromEnd")

    # Adjust chromosome format
    data[,"chrom"] <- GenomeInfoDb::MapSeqlevels(seqnames=data[,"chrom"], style=chromosome.format)

    # convert to GRanges object
    gr <- with(data, GenomicRanges::GRanges(seqnames=chrom,
                                            ranges=IRanges(start=chromStart+1,     # Convert from 0-based half open to 1-based closed
                                                           end=chromEnd)))

    return(gr)
}

readBed6File <- function(bedfile, skip=0, chromosome.format='NCBI') {
    ncols <- max(count.fields(bedfile, skip=skip))
    if ( ncols < 6 )
        stop("Not a correct BED6 file format.")

    data <- utils::read.table(bedfile, colClasses=c("character", rep("numeric", 2),
                                                    "character", "integer", "character",
                                                    rep("NULL", ncols-6)),
                              skip=skip)
    names(data) <- c("chrom", "chromStart", "chromEnd", "name", "score", "strand")

    # Adjust chromosome format
    data[,"chrom"] <- GenomeInfoDb::MapSeqlevels(seqnames=data[,"chrom"], style=chromosome.format)

    # Adjust strand information if provided
    data[,"strand"] <- sub('.', '*', data[,"strand"])

    # convert to GRanges object
    gr <- with(data, GenomicRanges::GRanges(seqnames=chrom,
                                            ranges=IRanges(start=chromStart+1,     # Convert from 0-based half open to 1-based closed
                                                           end=chromEnd),
                                            name=name,
                                            score=score,
                                            strand=strand))

    return(gr)
}

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
    data <- data.frame("chrom" = as.factor(GenomeInfoDb::MapSeqlevels(seqnames=seqnames(x), style=chromosome.format)),
                       "chromStart" = start(gr) - 1, # Convert from 1-based closed to 0-based half open
                       "chromEnd" = end(gr))

    utils::write.table(format(data, scientific=FALSE, trim=TRUE), file=filename.gz,
                       append=TRUE, row.names=FALSE, col.names=FALSE, quote=FALSE, sep='\t')

    close(filename.gz)
}
