\name{GstN}
\alias{GstN}
\title{Nei's Gst.}
\description{
This function estimates pairwise Gst among subpopulations (Nei 1973) from a GENEPOP data object (Rousset 2008). Missing genotype values in the GENEPOP file ("0000" or "000000") are simply ignored.
}
\usage{
GstN(popdata)
}
\arguments{
  \item{popdata}{Population data object created by read.genepop function from a GENEPOP file.}
}
\value{
Matrix of estimated pairwise Gst.}
\references{
Nei M (1973) Analysis of Gene Diversity in Subdivided Populations. \emph{Proc. Nat. Acad. Sci.}, 70, 3321-3323.\cr\cr
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.
}
\author{Reiichiro Nakamichi, Hirohisa Kishino, Shuichi Kitada}
\seealso{
\code{\link[FinePop]{read.genepop}}
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

# Gst estimation
result.gstN <- GstN(popdata)
write.csv(result.gstN, "result_GstN.csv", na="")
print(as.dist(result.gstN))
}
