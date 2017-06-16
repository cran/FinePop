\name{EBFST}
\alias{EBFST}
\title{Empirical Bayes estimator of Fst.}
\description{
This function estimates global/pairwise Fst among subpopulations using empirical Bayes method (Kitada et al. 2007, 2017). Preciseness of estimated pairwise Fst is evaluated by bootstrap method. This function accepts two types of data object, GENEPOP data (Rousset 2008) and allele (haplotype) frequency data (Kitada et al. 2007). Missing genotype values in the GENEPOP file ("0000" or "000000") are simply ignored.
}
\usage{
EBFST(popdata)
}
\arguments{
  \item{popdata}{Genotype data object of populations created by read.genepop function from a GENEPOP file. Allele (haplotype) frequency data object created by read.frequency function from a frequency format file also is acceptable.}
}
\details{
Frequency format file is a plain text file containing allele (haplotype) count data. This format is mainly for mitochondrial DNA (mtDNA) haplotype frequency data, however nuclear DNA (nDNA) data also is applicable. In the data object created by read.frequency function, "number of samples" means haplotype count. Therefore, it equals the number of individuals in mtDNA data, however it is the twice of the number of individuals in nDNA data. First part of the frequency format file is the number of subpopulations, second part is the number of loci, and latter parts are [\emph{population} x \emph{allele}] matrices of the observed allele (haplotype) counts at each locus. Two examples of frequency format files are attached in this package. See \code{\link[FinePop]{jsmackerel}}.
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
Kitada S, Nakamichi R, Kishino H (2017) The empirical Bayes estimators of fine-scale population structure in high gene flow species. \emph{Mol. Ecol. Resources}, DOI: 10.1111/1755-0998.12663 \cr\cr
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.
}
\author{Reiichiro Nakamichi, Hirohisa Kishino, Shuichi Kitada}
\seealso{
\code{\link[FinePop]{read.genepop}}, \code{\link[FinePop]{read.frequency}},
\code{\link[stats]{as.dist}}, \code{\link[stats]{as.dendrogram}},
\code{\link[stats]{hclust}}, \code{\link[stats]{cmdscale}}, \code{\link[ape]{nj}}
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

# Fst estimation
result.eb <- EBFST(popdata)
ebfst <- result.eb$pairwise$fst
write.csv(ebfst, "result_EBFST.csv", na="")
ebfst.d <- as.dist(ebfst)
print(ebfst.d)

# dendrogram
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
