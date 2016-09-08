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
Kitada S, Kitakado T, Kishino H (2007) Empirical Bayes inference of pairwise FST and its distribution in the genome. Genetics, 177, 861-873.\cr\cr
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. Mol. Ecol. Resources 8, 103-106.
}
\author{Reiichiro Nakamichi, Hirohisa Kishino, Shuichi Kitada}
\seealso{
\code{\link[FinePop]{read.genepop}}, \code{\link[stats]{as.dist}}, \code{\link[stats]{as.dendrogram}},
\code{\link[stats]{hclust}}, \code{\link[stats]{cmdscale}}, \code{\link[ape]{nj}}
}
\examples{
# load data
data(data.JSM)
cat(data.JSM$genepop, file="JSM_genepop.txt", sep="\n")
popdata <- read.genepop("JSM_genepop.txt")
name.pop <- data.JSM$popname
num.pop <- popdata$npop

# Fst estimation
result <- EBFST(popdata)
fst <- result$pairwise$fst
dimnames(fst) <- list(name.pop[-1], name.pop[-num.pop])
write.csv(fst, "result_EBFST.csv", na="")

# dendrogram
fst.d <- as.dist(fst)
fst.hc <- hclust(fst.d,method="average")
plot(as.dendrogram(fst.hc), xlab="",ylab="",main="", las=1)

# MDS plot
mds <- cmdscale(fst.d)
plot(mds, type="n", xlab="",ylab="")
text(mds[,1],mds[,2], name.pop)

# NJ tree
library(ape)
fst.nj <- nj(fst.d)
plot(fst.nj,type="u",main="",sub="")
}
