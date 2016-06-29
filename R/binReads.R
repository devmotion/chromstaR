

#' Convert aligned reads from various file formats into read counts in equidistant bins
#'
#' Convert aligned reads in .bam or .bed(.gz) format into read counts in equidistant windows.
#'
#' Convert aligned reads from .bam or .bed(.gz) files into read counts in equidistant windows (bins). This function uses \code{\link[GenomicRanges]{countOverlaps}} to calculate the read counts, or alternatively \code{\link[bamsignals]{bamProfile}} if option \code{use.bamsignals} is set (only effective for .bam files).
#'
#' @aliases binning
#' @param file A file with aligned reads. Alternatively a \code{\link{GRanges}} with aligned reads if format is set to 'GRanges'.
#' @param experiment.table An \code{\link{experiment.table}} containing the supplied \code{file}. This is necessary to uniquely identify the file in later steps of the workflow. Set to \code{NULL} if you don't have it (not recommended).
#' @inheritParams readBamFileAsGRanges
#' @inheritParams readBedFileAsGRanges
#' @param blacklist A \code{\link{GRanges}} or a bed(.gz) file with blacklisted regions. Bins falling into those regions will be discarded.
#' @param binsizes An integer vector specifying the bin sizes to use.
#' @param bins A named \code{list} with \code{\link{GRanges}} containing precalculated bins produced by \code{\link{fixedWidthBins}} or \code{\link{variableWidthBins}}. Names must correspond to the binsize.
#' @param reads.per.bin Approximate number of desired reads per bin. The bin size will be selected accordingly.
#' @param variable.width.reference A BAM file that is used as reference to produce variable width bins. See \code{\link{variableWidthBins}} for details.
#' @param chromosomes If only a subset of the chromosomes should be binned, specify them here.
#' @param use.bamsignals If \code{TRUE} the \pkg{\link[bamsignals]{bamsignals}} package is used for parsing of BAM files. This gives tremendous speed advantage for only one binsize but linearly increases for multiple binsizes, while \code{use.bamsignals=FALSE} has a binsize dependent runtime and might be faster if many binsizes are calculated.
#' @return If only one bin size was specified for option \code{binsizes}, the function returns a single \code{\link{GRanges}} object with meta data column 'counts' that contains the read count. If multiple \code{binsizes} were specified , the function returns a named \code{list()} of \link{GRanges} objects.
#' @importFrom Rsamtools BamFile indexBam
#' @importFrom bamsignals bamCount
#' @export
#'
#'@examples
#'## Get an example BAM file with ChIP-seq reads
#'file <- system.file("extdata", "euratrans",
#'                    "lv-H3K27me3-BN-male-bio2-tech1.bam",
#'                     package="chromstaRData")
#'## Bin the file into bin size 1000bp
#'data(rn4_chrominfo)
#'data(experiment_table)
#'binned <- binReads(file, experiment.table=experiment_table,
#'                   assembly=rn4_chrominfo, binsizes=1000,
#'                   chromosomes='chr12')
#'print(binned)
#'
binReads <- function(file, experiment.table=NULL, assembly, bamindex=file, chromosomes=NULL, pairedEndReads=FALSE, min.mapq=10, remove.duplicate.reads=TRUE, max.fragment.width=1000, blacklist=NULL, binsizes=1000, reads.per.bin=NULL, bins=NULL, variable.width.reference=NULL, use.bamsignals=TRUE) {

    ## Determine format
    if (is.character(file)) {
        file.clean <- sub('\\.gz$','', file)
        format <- rev(strsplit(file.clean, '\\.')[[1]])[1]
    } else if (class(file)=='GRanges') {
        format <- 'GRanges'
    }

    if (format=='bed') {
        temp <- assembly # trigger error if not defined
    }

    ## Variables for bamsignals
    paired.end <- ifelse(pairedEndReads, 'filter', 'ignore')
    
    ## Create INFO object as row from the experiment.table
    if (!is.null(experiment.table)) {
        check.experiment.table(experiment.table)
        info <- experiment.table[basename(as.character(experiment.table$file))==basename(file),]
        ID <- paste0(info$mark, '-', info$condition, '-rep', info$rep)
        info$ID <- ID
        if (pairedEndReads != info$pairedEndReads) {
            warning("Option 'pairedEndReads' overwritten by entry in 'experiment.table'.")
        }
    } else {
        info <- NULL
        ID <- ifelse(format=='GRanges', 'GRanges', basename(file))
    }

    ### Read in the data
    data <- NULL
    if (format == "bed") {
        ## BED (0-based)
        data <- readBedFileAsGRanges(file, assembly=assembly, chromosomes=chromosomes, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, max.fragment.width=max.fragment.width)
        chrom.lengths <- seqlengths(data)
    } else if (format == "bam") {
        ## BAM (1-based)
        if (use.bamsignals) {
            ## Check if bamindex exists
            bamindex.raw <- sub('\\.bai$', '', bamindex)
            bamindex <- paste0(bamindex.raw,'.bai')
            if (!file.exists(bamindex)) {
                ptm <- startTimedMessage("Making bam-index file ...")
                bamindex.own <- Rsamtools::indexBam(file)
                # warning("Couldn't find BAM index-file ",bamindex,". Creating our own file ",bamindex.own," instead.")
                bamindex <- bamindex.own
                stopTimedMessage(ptm)
            }
            ptm <- startTimedMessage(paste0("Reading header from ", file, " ..."))
            chrom.lengths <- GenomeInfoDb::seqlengths(Rsamtools::BamFile(file))
            stopTimedMessage(ptm)
        } else {
            data <- readBamFileAsGRanges(file, bamindex, chromosomes=chromosomes, pairedEndReads=pairedEndReads, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, max.fragment.width=max.fragment.width)
            chrom.lengths <- seqlengths(data)
        }
    } else if (format == "GRanges") {
        ## GRanges (1-based)
        data <- file
        chrom.lengths <- seqlengths(data)
    }

    ## Select chromosomes to bin
    if (is.null(chromosomes)) {
        chromosomes <- names(chrom.lengths)
    }
    chroms2use <- intersect(chromosomes, names(chrom.lengths))
  	## Stop if none of the specified chromosomes exist
  	if (length(chroms2use)==0) {
  		chrstring <- paste0(chromosomes, collapse=', ')
  		stop('Could not find length information for any of the specified chromosomes: ', chrstring)
  	}
    skipped.chroms <- setdiff(chromosomes, chroms2use)
    if (length(skipped.chroms)>0) {
        warning("Could not find chromosomes ", paste0(skipped.chroms, collapse=', '), ".")
    }

    ## Select only desired chromosomes
    if (!is.null(data)) {
        ptm <- startTimedMessage("Subsetting specified chromosomes ...")
        data <- data[seqnames(data) %in% chroms2use]
        data <- keepSeqlevels(data, as.character(unique(seqnames(data))))
        ## Drop seqlevels where seqlength is NA
        na.seqlevels <- seqlevels(data)[is.na(seqlengths(data))]
        data <- data[seqnames(data) %in% seqlevels(data)[!is.na(seqlengths(data))]]
        data <- keepSeqlevels(data, as.character(unique(seqnames(data))))
        if (length(na.seqlevels) > 0) {
            warning("Dropped seqlevels because no length information was available: ", paste0(na.seqlevels, collapse=', '))
        }
        stopTimedMessage(ptm)
    }
  
    ### Make variable width bins ###
    if (!is.null(variable.width.reference)) {
        message("Making variable width bins:")
        if (is.character(variable.width.reference)) {
            variable.width.reference.clean <- sub('\\.gz$','', variable.width.reference)
            vformat <- rev(strsplit(variable.width.reference.clean, '\\.')[[1]])[1]
        } else if (class(variable.width.reference)=='GRanges') {
            vformat <- 'GRanges'
        }
        if (vformat == 'bam') {
            refreads <- readBamFileAsGRanges(variable.width.reference, bamindex=variable.width.reference, chromosomes=chroms2use, pairedEndReads=pairedEndReads, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, max.fragment.width=max.fragment.width)
        } else if (vformat == 'bed') {
            refreads <- readBedFileAsGRanges(variable.width.reference, assembly=assembly, chromosomes=chroms2use, remove.duplicate.reads=remove.duplicate.reads, min.mapq=min.mapq, max.fragment.width=max.fragment.width)
        }
        bins.binsize <- variableWidthBins(refreads, binsizes=binsizes, chromosomes=chroms2use)
        message("Finished making variable width bins.")
    }
        
    ### Make fixed width bins ###
    if (is.null(variable.width.reference)) {
        bins.binsize <- fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=chroms2use, binsizes=binsizes)
    }
  
    ### Make reads.per.bin bins ###
    bins.rpb <- NULL
    if (!is.null(data)) {
        numcountsperbp <- length(data) / sum(as.numeric(seqlengths(data)))
        binsizes.rpb <- round(reads.per.bin / numcountsperbp, -2)
        bins.rpb <- fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=chroms2use, binsizes=binsizes.rpb)
    } else if (use.bamsignals & !is.null(reads.per.bin)) {
        ptm <- startTimedMessage("Parsing bamfile to determine binsize for reads.per.bin option ...")
        bins.helper <- suppressMessages( fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=chroms2use, binsizes=1e6)[[1]] )
        counts <- bamsignals::bamCount(file, bins.helper, mapqual=min.mapq, paired.end=paired.end, tlenFilter=c(0, max.fragment.width), verbose=FALSE)
        stopTimedMessage(ptm)
        numcountsperbp <- sum(as.numeric(counts)) / sum(as.numeric(chrom.lengths[chroms2use]))
        binsizes.rpb <- round(reads.per.bin / numcountsperbp, -2)
        bins.rpb <- fixedWidthBins(chrom.lengths=chrom.lengths, chromosomes=chroms2use, binsizes=binsizes.rpb)
    }
    
    ### Combine in bins.list ###
    bins.list <- c(bins, bins.binsize, bins.rpb)
 
    ### Loop over all binsizes ###
    if (!use.bamsignals) {
        ptm <- startTimedMessage("Splitting into strands ...")
        data.plus <- data[strand(data)=='+']
        data.minus <- data[strand(data)=='-']
        data.star <- data[strand(data)=='*']
        ptm <- stopTimedMessage(ptm)
    }
    for (ibinsize in 1:length(bins.list)) {
        binsize <- as.numeric(names(bins.list)[ibinsize])

        ## Filter blacklisted bins
        bins <- blacklistGRanges(bins.list[[ibinsize]], blacklist)

        if (format == 'bam' & use.bamsignals) {
            ptm <- startTimedMessage("Counting overlaps for binsize ", binsize, " ...")
            bins$counts <- bamsignals::bamCount(file, bins, mapqual=min.mapq, paired.end=paired.end, tlenFilter=c(0, max.fragment.width), verbose=FALSE)
            stopTimedMessage(ptm)    
        } else {
            readsperbin <- round(length(data) / sum(as.numeric(seqlengths(data))) * binsize, 2)
            ptm <- startTimedMessage("Counting overlaps for binsize ",binsize," with on average ",readsperbin," reads per bin ...")
            scounts <- suppressWarnings( GenomicRanges::countOverlaps(bins, data.star) )
            mcounts <- suppressWarnings( GenomicRanges::countOverlaps(bins, data.minus) )
            pcounts <- suppressWarnings( GenomicRanges::countOverlaps(bins, data.plus) )
            counts <- mcounts + pcounts + scounts
            countmatrix <- matrix(c(counts,mcounts,pcounts), ncol=3)
            colnames(countmatrix) <- c('counts','mcounts','pcounts')
            mcols(bins)$counts <- countmatrix[,'counts']
            stopTimedMessage(ptm)
        }

        attr(bins, 'info') <- info
        attr(bins, 'min.mapq') <- min.mapq
        bins.list[[ibinsize]] <- bins    
    } ### end loop binsizes ###

    if (length(bins.list) == 1) {
        return(bins.list[[1]])
    } else {
        return(bins.list)
    }
}

#' Remove blacklist from GRanges
#'
#' This is a simple convenience function to remove regions from a bed(.gz)-file from a \code{\link{GRanges}} object. The bed-file is expected to have at least the following three fields: \code{chromosome, start, end}.
#'
#' @param gr A \code{\link{GRanges}} object.
#' @param blacklist Filename of the bed or bed.gz file which contains the blacklisted regions.
#' @return A \code{\link{GRanges}} object with the regions not overlapping any regions of the blacklist.
#' @importFrom utils write.table
#' @author Aaron Taudt, David Widmann
blacklistGRanges <- function(gr, blacklist=NULL) {
    if (is.null(blacklist))
        return(gr)

    ## Input checks
    if ( !(is.character(blacklist) | class(blacklist)=='GRanges') )
        stop("'blacklist' has to be either a bed(.gz) file or a GRanges object")

    ptm <- startTimedMessage("Filtering blacklisted regions ...")

    if (is.character(blacklist))
        blacklist <- readBedFile(blacklist, skip=0)

    # Convert both 'gr' and 'blacklist' to the same chromsome format
    chromosome.format <- GenomeInfoDb::seqlevelsStyle(gr)
    seqnames(blacklist) <- GenomeInfoDb::mapSeqlevels(seqnames=as.character(seqnames(gr)),
                                                      style=chromosome.format)

    overlaps <- findOverlaps(gr, blacklist)
    idx <- setdiff(1:length(gr), S4Vectors::queryHits(overlaps))
    gr <- gr[idx]
    stopTimedMessage(ptm)

    return(gr)
}


