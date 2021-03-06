% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combinatorialStates.R
\name{combinatorialStates}
\alias{combinatorialStates}
\title{Get the (decimal) combinatorial states of a list of univariate HMM models}
\usage{
combinatorialStates(hmm.list, binary = FALSE)
}
\arguments{
\item{hmm.list}{A list of models generated by \code{\link{callPeaksUnivariate}}, e.g. 'list(model1,model2,...)'.}

\item{binary}{If \code{TRUE}, a matrix of binary instead of decimal states will be returned.}
}
\value{
Output is a vector of integers representing the combinatorial state of each bin.
}
\description{
Get the combinatorial states of a list of models generated by \code{\link{callPeaksUnivariate}}. The function returns the decimal combinatorial states for each bin (see details for an explanation of combinatorial state).
}
\details{
For a given model, each genomic bin can be either called 'unmodified' or 'modified', depending on the posterior probabilities estimated by the Baum-Welch. Thus, a list of models defines a binary combinatorial state for each bin. This binary combinatorial state can be expressed as a decimal number.
Example: We have 4 histone modifications, and we run the univariate HMM for each of them. Then we use a false discovery rate of 0.5 to call each bin either 'unmodified' or 'modified'. The resulting binary combinatorial states can then be converted to decimal representation. The following table illustrates this:

\tabular{crrrrr}{
bin \tab modification state \tab \tab \tab \tab decimal state\cr
    \tab model1 \tab model2 \tab model3 \tab model4 \tab     \cr
1   \tab      0 \tab      0 \tab      1 \tab      0 \tab 2   \cr
2   \tab      0 \tab      0 \tab      0 \tab      0 \tab 0   \cr
3   \tab      0 \tab      1 \tab      1 \tab      0 \tab 6   \cr
4   \tab      0 \tab      1 \tab      1 \tab      1 \tab 7   \cr
}
}
\examples{
# Get example BAM files for 3 different marks in hypertensive rat (SHR)
file.path <- system.file("extdata","euratrans", package='chromstaRData')
files <- list.files(file.path, full.names=TRUE, pattern='SHR.*bam$')[c(1,4,6)]
# Bin the data
data(rn4_chrominfo)
binned.data <- list()
for (file in files) {
 binned.data[[basename(file)]] <- binReads(file, binsizes=1000, stepsizes=500,
                                           assembly=rn4_chrominfo, chromosomes='chr12')
}
# Obtain the univariate fits
models <- list()
for (i1 in 1:length(binned.data)) {
 models[[i1]] <- callPeaksUnivariate(binned.data[[i1]], max.time=60, eps=1)
}
## Get the decimal representation of the combinatorial state of this combination of models
states <- chromstaR:::combinatorialStates(models, binary=FALSE)
## Show number of each state
table(states)

}
\author{
Aaron Taudt
}
\seealso{
\code{\link{dec2bin}}, \code{\link{bin2dec}}
}

