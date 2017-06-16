\name{thetaWC.pair}
\alias{thetaWC.pair}
\title{Weir and Cockerham's theta adapted for pairwise Fst.}
\description{
This function estimates Fst between population pairs based on Weir and Cockerham's theta (Weir & Cockerham 1984) adapted for pairwise comparison from a GENEPOP data object (Rousset 2008). Missing genotype values in the GENEPOP file ("0000" or "000000") are simply ignored.
}
\usage{
thetaWC.pair(popdata)
}
\arguments{
  \item{popdata}{Population data object created by read.genepop function from a GENEPOP file.}
}
\details{
Weir and Cockerham (1984) derived an unbiased estimator of a coancestry coefficient (theta) based on a random effect model. It expresses the extent of genetic heterogeneity within the population. The second stage common approach is to investigate the detailed pattern of the population structure, based on a measure of genetic difference between pairs of subpopulations (demes). We call this by pairwise Fst. This function follows the formula of Weir and Cockerham's theta with the sample size \emph{r} = 2. Given the pair, our finite sample correction multiplies \emph{a} of Weir & Cockerham's theta by (\emph{r} - 1) / \emph{r} (equation 2 in p.1359 of Weir & Cockerham 1984).
}
\value{
Matrix of estimated pairwise Fst by theta with finite sample correction.
}
\references{
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.\cr\cr
Weir BS, Cockerham CC (1984) Estimating F-statistics for the analysis of population structure. \emph{Evolution}, 38, 1358-1370.
}
\author{Reiichiro Nakamichi, Hirohisa Kishino, Shuichi Kitada}
\seealso{
\code{\link[FinePop]{read.genepop}}
}
\examples{
# Example of GENEPOP file
data(jsmackerel)
cat(jsmackerel$MS.genepop, file="JSM_MS_genepop.txt", sep="\n")
cat(jsmackerel$popname, file="JSM_popname.txt", sep=" ")

# Data load
# Prepare your GENEPOP file and population name file in the working directory
# (Here, these files were provided as "JSM_MS_genepop.txt" and "JSM_popname.txt".)
popdata <- read.genepop(genepop="JSM_MS_genepop.txt", popname="JSM_popname.txt")

# theta estimation
result.theta.pair <- thetaWC.pair(popdata)
write.csv(result.theta.pair, "result_thetaWCpair.csv", na="")
print(as.dist(result.theta.pair))
}
