\name{read.frequency}
\alias{read.frequency}
\title{Create an allele (haplotype) frequency data object of populations from a frequency format file.}
\description{
This function reads a frequency format file (Kitada et al. 2007) and parse it into an R data object. This data object provides a summary of allele (haplotype) frequency in each population and marker status. This data object is used by EBFST function of this package.
}
\usage{
read.frequency(frequency, popname = NULL)
}
\arguments{
  \item{frequency}{A character value specifying the name of the frequency format file to be analyzed.}
  \item{popname}{A character value specifying the name of the plain text file containing the names of subpopulations to be analyzed. This text file must not contain other than subpopulation names. The names must be separated by spaces, tabs or line breaks. If this argument is omitted, serial numbers will be assigned as subpopulation names.}
}
\details{
Frequency format file is a plain text file containing allele (haplotype) count data. This format is mainly for a mitochondrial DNA (mtDNA) haplotype frequency data, however nuclear DNA (nDNA) data also is applicable. In the data object created by read.frequency function, "number of samples" means haplotype count. Therefore, it equals the number of individuals in mtDNA data, however it is the twice of the number of individuals in nDNA data. First part of the frequency format file is the number of subpopulations, second part is the number of loci, and latter parts are [\emph{population} x \emph{allele}] matrices of the observed allele (haplotype) counts at each locus. Two examples of frequency format files are attached in this package. See \code{\link[FinePop]{data.JSM}}.
}
\value{
  \item{npops}{Number of subpopulations.}
  \item{pop_sizes}{Number of samples in each subpopulation.}
  \item{pop_names}{Names of subpopulations.}
  \item{nloci}{Number of loci.}
  \item{loci_names}{Names of loci.}
  \item{all_alleles}{A list of alleles (haplotypes) at each locus.}
  \item{nalleles}{Number of alleles (haplotypes) at each locus.}
  \item{indtyp}{Number of genotyped samples in each subpopulation at each locus.}
  \item{obs_allele_num}{Observed allele (haplotype) counts at each locus in each subpopulation.}
  \item{allele_freq}{Observed allele (haplotype) frequencies at each locus in each subpopulation.}
  \item{call_rate}{Rate of genotyped samples at each locus.}
}
\references{
Kitada S, Kitakado T, Kishino H (2007) Empirical Bayes inference of pairwise FST and its distribution in the genome. \emph{Genetics}, 177, 861-873.
}
\author{Reiichiro Nakamichi, Hirohisa Kishino, Shuichi Kitada}
\seealso{
\code{\link[FinePop]{data.JSM}}, \code{\link[FinePop]{EBFST}}
}
\examples{
# data
data(data.JSM)
cat(data.JSM$mtDNA.freq, file="JSM_mtDNA_freq.txt", sep="\n")
cat(data.JSM$popname, file="JSM_popname.txt", sep=" ")

# read frequency format file with subpopulation names
popdata.mt <- read.frequency(frequency="JSM_mtDNA_freq.txt", popname="JSM_popname.txt")

# read frequency file without subpopulation names
popdata.mt.noname <- read.frequency(frequency="JSM_mtDNA_freq.txt")

# Fst estimation
result.eb.mt <- EBFST(popdata.mt)
ebfst.mt <- result.eb.mt$pairwise$fst
write.csv(ebfst.mt, "result_EBFST_mtDNA.csv", na="")
}
