\name{GstH}
\alias{GstH}
\title{Hedrick's G'st}
\description{
This function estimates pairwise G'st (Hedrick 2005) among subpopulations from a GENEPOP data object (Rousset 2008). Missing genotype values in the GENEPOP file ("0000" or "000000") are simply ignored.
}
\usage{
GstH(popdata)
}
\arguments{
  \item{popdata}{Population data object created by read.genepop function from a GENEPOP file.}
}
\value{
Matrix of estimated pairwise Hedrick's G'st.
}
\references{
Hedrick P (2005) A standardized genetic differentiation measure. \emph{Evolution}, 59, 1633-1638.\cr\cr
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.
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

# Hedrick's G'st estimation
result.GstH <- GstH(popdata)
write.csv(result.GstH, "result_GstH.csv", na="")
print(as.dist(result.GstH))
}
