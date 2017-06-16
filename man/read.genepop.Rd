\name{read.genepop}
\alias{read.genepop}
\title{Create a genotype data object of populations from a GENEPOP format file.}
\description{
This function reads a GENEPOP format file (Rousset 2008) and parse it into an R data object. This data object provides a summary of genotype/haplotype of each sample, allele frequency in each population, and marker status. This data object is used in downstream analysis of this package. This function is a "lite" version of readGenepop function in diveRsity package (Keenan 2015).
}
\usage{
read.genepop(genepop, popname = NULL)
}
\arguments{
  \item{genepop}{A character value specifying the name of the GENEPOP file to be analyzed.}
  \item{popname}{A character value specifying the name of the plain text file containing the names of subpopulations to be analyzed. This text file must not contain other than subpopulation names. The names must be separated by spaces, tabs or line breaks. If this argument is omitted, serial numbers will be assigned as subpopulation names.}
}
\value{
  \item{npops}{Number of subpopulations.}
  \item{pop_sizes}{Number of samples in each subpopulation.}
  \item{pop_names}{Names of subpopulations.}
  \item{nloci}{Number of loci.}
  \item{loci_names}{Names of loci.}
  \item{all_alleles}{A list of alleles at each locus.}
  \item{nalleles}{Number of alleles at each locus.}
  \item{indtyp}{Number of genotyped samples in each subpopulation at each locus.}
  \item{ind_names}{Names of samples in each subpopulation.}
  \item{pop_alleles}{Genotypes of each sample at each locus in haploid designation.}
  \item{pop_list}{Genotypes of each sample at each locus in diploid designation.}
  \item{obs_allele_num}{Observed allele counts at each locus in each subpopulation.}
  \item{allele_freq}{Observed allele frequencies at each locus in each subpopulation.}
  \item{call_rate}{Rate of genotyped samples at each locus.}
}
\references{
Keenan K (2015) diveRsity: A Comprehensive, General Purpose Population Genetics Analysis Package. R package version 1.9.89. \url{https://CRAN.R-project.org/package=diveRsity}\cr\cr
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.
}
\author{Reiichiro Nakamichi}
\seealso{
\code{\link[diveRsity]{readGenepop}}
}
\examples{
# Example of GENEPOP file
data(jsmackerel)
cat(jsmackerel$MS.genepop, file="JSM_MS_genepop.txt", sep="\n")
cat(jsmackerel$popname, file="JSM_popname.txt", sep=" ")

# Read GENEPOP file with subpopulation names
# Prepare your GENEPOP file and population name file in the working directory
# (Here, these files were provided as "JSM_MS_genepop.txt" and "JSM_popname.txt".)
popdata <- read.genepop(genepop="JSM_MS_genepop.txt", popname="JSM_popname.txt")

# Read GENEPOP file without subpopulation names
popdata.noname <- read.genepop(genepop="JSM_MS_genepop.txt")
}
