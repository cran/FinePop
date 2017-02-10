\name{EBDJ}
\alias{EBDJ}
\title{Empirical Bayes estimator of Jost's D}
\description{
This function estimates pairwise D (Jost 2008) among subpopulations using empirical Bayes method (Kitada et al. 2007). This function accepts two types of data object, GENEPOP data (Rousset 2008) and allele (haplotype) frequency data (Kitada et al. 2007). Missing genotype values in the GENEPOP file ("0000" or "000000") are simply ignored.
}
\usage{
EBDJ(popdata)
}
\arguments{
  \item{popdata}{Genotype data object of populations created by read.genepop function from a GENEPOP file. Allele (haplotype) frequency data object created by read.frequency function from a frequency format file also is acceptable.}
}
\details{
Frequency format file is a plain text file containing allele (haplotype) count data. This format is mainly for mitochondrial DNA (mtDNA) haplotype frequency data, however nuclear DNA (nDNA) data also is applicable. In the data object created by read.frequency function, "number of samples" means haplotype count. Therefore, it equals the number of individuals in mtDNA data, however it is the twice of the number of individuals in nDNA data. First part of the frequency format file is the number of subpopulations, second part is the number of loci, and latter parts are [\emph{population} x \emph{allele}] matrices of the observed allele (haplotype) counts at each locus. Two examples of frequency format files are attached in this package. See \code{\link[FinePop]{data.JSM}}.
}
\value{
Matrix of estimated pairwise Jost's D.
}
\references{
Jost L (2008) Gst and its relatives do not measure differentiation. \emph{Molecular Ecology}, 17, 4015-4026.\cr\cr
Kitada S, Kitakado T, Kishino H (2007) Empirical Bayes inference of pairwise FST and its distribution in the genome. \emph{Genetics}, 177, 861-873.\cr\cr
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.
}
\author{Reiichiro Nakamichi, Hirohisa Kishino, Shuichi Kitada}
\seealso{
\code{\link[FinePop]{read.genepop}}, \code{\link[FinePop]{read.frequency}}
}
\examples{
# Example of GENEPOP file
data(data.JSM)
cat(data.JSM$MS.genepop, file="JSM_MS_genepop.txt", sep="\n")
cat(data.JSM$popname, file="JSM_popname.txt", sep=" ")

# Data load
# Prepare your GENEPOP file and population name file in the working directory
# (Here, these files were provided as "JSM_MS_genepop.txt" and "JSM_popname.txt".)
popdata <- read.genepop(genepop="JSM_MS_genepop.txt", popname="JSM_popname.txt")

# Jost's D estimation
result.EBDJ <- EBDJ(popdata)
write.csv(result.EBDJ, "result_EBDJ.csv", na="")
print(as.dist(result.EBDJ))
}
