% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/readCustomBedFile.R
\name{readCustomBedFile}
\alias{readCustomBedFile}
\title{Read bed-file into GRanges}
\usage{
readCustomBedFile(bedfile, col.names = c("chromosome", "start", "end", "name",
  "score", "strand"), col.classes = NULL, skip = 0,
  chromosome.format = "NCBI")
}
\arguments{
\item{bedfile}{Filename of the bed or bed.gz file.}

\item{col.names}{A character vector giving the names of the columns in the \code{bedfile}. Must contain at least \code{c('chromosome','start','end')}.}

\item{col.classes}{A character vector giving the classes of the columns in \code{bedfile}. Speeds up the import.}

\item{skip}{Number of lines to skip at the beginning.}

\item{chromosome.format}{Desired format of the chromosomes. Either 'NCBI' for (1,2,3 ...) or 'UCSC' for (chr1,chr2,chr3 ...).}
}
\value{
A \code{\link{GRanges}} object with the contents of the bed-file.
}
\description{
This is a simple convenience function to read a bed(.gz)-file into a \code{\link{GRanges}} object. The bed-file is expected to have the following fields: \code{chromosome, start, end, name, score, strand}.
}
\examples{
## Get an example BED file
bedfile <- system.file("extdata", "liver-H3K4me3-BN-male-bio2-tech1.bed.gz",
                       package="chromstaRData")
## Import the file and skip the first 10 lines
data <- readCustomBedFile(bedfile, skip=10)

}
\author{
Aaron Taudt
}

