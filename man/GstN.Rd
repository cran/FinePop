\name{GstN}
\alias{GstN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nei's Gst.
}
\description{
This function estimates pairwise Gst among subpopulations (Nei 1973).
}
\usage{
GstN(popdata)
}
\arguments{
  \item{popdata}{Population data object created by read.genepop function from a GENEPOP file (Rousset 2008).}
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
# data
data(data.JSM)
cat(data.JSM$genepop, file="JSM_genepop.txt", sep="\n")
cat(data.JSM$popname, file="JSM_popname.txt", sep=" ")
popdata <- read.genepop(genepop="JSM_genepop.txt", popname="JSM_popname.txt")

# Gst estimation
result.gst <- GstN(popdata)
write.csv(result.gst, "result_Gst.csv", na="")
}
