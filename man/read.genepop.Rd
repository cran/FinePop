\name{read.genepop}
\alias{read.genepop}
\title{Create a genetic data object of populations from a GENEPOP format file.}
\description{
This function reads a GENEPOP format file (Rousset 2008) and parse it into an R data object. This data object provides a summary of genotype/haplotype of each sample, allele frequency in each population, and marker status. This data object is used in downstream analysis of this package. This function is a "lite" version of readGenepop function in diveRsity package (Keenan 2015).
}
\usage{
read.genepop(genepop)
}
\arguments{
  \item{genepop}{A character value specifying the name of the GENEPOP file to be analyzed.}
}
\value{
  \item{npops}{Number of subpopulations.}
  \item{nloci}{Number of loci.}
  \item{pop_sizes}{Number of samples in each subpopulation.}
  \item{loci_names}{Names of loci.}
  \item{all_alleles}{A list of alleles at each locus.}
  \item{nalleles}{Number of alleles at each locus.}
  \item{indtyp}{Number of genotyped samples in each subpopulation at each locus.}
  \item{ind_names}{Names of samples in each subpopulation.}
  \item{pop_alleles}{Genotypes of each sample at each locus in haploid designation.}
  \item{pop_list}{Genotypes of each sample at each locus in diploid designation.}
  \item{obs_allele_num}{Allele counts of each alleles at each locus in each subpopulation.}
  \item{allele_freq}{Allele frequencies of each alleles at each locus in each subpopulation.}
  \item{call_rate}{Rate of genotyped samples at each locus.}
}
\references{
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. Mol. Ecol. Resources 8, 103-106.\cr
Keenan K (2015) diveRsity: A Comprehensive, General Purpose Population Genetics Analysis Package. R package version 1.9.89. \url{http://CRAN.R-project.org/package=diveRsity}
}
\author{Reiichiro Nakamichi}
\seealso{
\code{\link[diveRsity]{readGenepop}}
}
\examples{
data(data.JSM)
cat(data.JSM$genepop, file="JSM_genepop.txt", sep="\n")
popdata <- read.genepop("JSM_genepop.txt")
}
