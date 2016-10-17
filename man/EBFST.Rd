\name{EBFST}
\alias{EBFST}
\title{Empirical Bayes estimator of Fst.}
\description{
This function estimates global/pairwise Fst among subpopulations using empirical Bayes method (Kitada et al. 2007). Preciseness of estimated pairwise Fst is evaluated by bootstrap method.}
\usage{
EBFST(popdata)
}
\arguments{
  \item{popdata}{Population data object created by read.genepop function from a GENEPOP file (Rousset 2008).}
}
\value{
global:
  \item{theta}{Rate of gene flow.}
  \item{fst}{Estimated global Fst.}
pairwise:
  \item{fst}{Estimated pairwise Fst.}
  \item{fst.boot}{Bootstrap mean of estimated Fst.}
  \item{fst.boot.sd}{Bootstrap standard deviation of estimated Fst.}
}
\references{
Kitada S, Kitakado T, Kishino H (2007) Empirical Bayes inference of pairwise FST and its distribution in the genome. \emph{Genetics}, 177, 861-873.\cr\cr
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.
}
\author{Reiichiro Nakamichi, Hirohisa Kishino, Shuichi Kitada}
\seealso{
\code{\link[FinePop]{read.genepop}}, \code{\link[stats]{as.dist}}, \code{\link[stats]{as.dendrogram}},
\code{\link[stats]{hclust}}, \code{\link[stats]{cmdscale}}, \code{\link[ape]{nj}}
}
\examples{
# data
data(data.JSM)
cat(data.JSM$genepop, file="JSM_genepop.txt", sep="\n")
cat(data.JSM$popname, file="JSM_popname.txt", sep=" ")
popdata <- read.genepop(genepop="JSM_genepop.txt", popname="JSM_popname.txt")

# Fst estimation
result.eb <- EBFST(popdata)
ebfst <- result.eb$pairwise$fst
write.csv(ebfst, "result_EBFST.csv", na="")

# dendrogram
ebfst.d <- as.dist(ebfst)
ebfst.hc <- hclust(ebfst.d,method="average")
plot(as.dendrogram(ebfst.hc), xlab="",ylab="",main="", las=1)

# MDS plot
mds <- cmdscale(ebfst.d)
plot(mds, type="n", xlab="",ylab="")
text(mds[,1],mds[,2], popdata$pop_names)

# NJ tree
library(ape)
ebfst.nj <- nj(ebfst.d)
plot(ebfst.nj,type="u",main="",sub="")
}
