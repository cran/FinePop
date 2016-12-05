\name{data.JSM}
\alias{data.JSM}
\docType{data}
\title{An example dataset of GENEPOP format and frequency format.}
\description{An example of a genetic data for a Japanese Spanish mackerel population (Nakajima et al. 2014). It contains genotypic information of 5 microsatellite markers and mtDNA D-loop region from 8 subpopulations of 715 individuals. GENEPOP format (Rousset 2008) and frequency format (Kitada et al. 2007) text files are available. Name list of subpopulations also is attached.}
\usage{data("data.JSM")}
\format{
 $ MS.genepop: Genotypic information of 5 microsatellites in GENEPOP format text data.\cr
 $ MS.freq: Allele frequency of 5 microsatellites  in frequency format text data.\cr
 $ mtDNA.freq: Haplotype frequency of mtDNA D-loop region in frequency format text data.\cr
 $ popname: Names of subpopulations.
}
\details{
Frequency format file is a plain text file containing allele (haplotype) count data. This format is mainly for mitochondrial DNA (mtDNA) haplotype frequency data, however nuclear DNA (nDNA) data also is applicable. In the data object created by read.frequency function, "number of samples" means haplotype count. Therefore, it equals the number of individuals in mtDNA data, however it is the twice of the number of individuals in nDNA data. First part of the frequency format file is the number of subpopulations, second part is the number of loci, and latter parts are [\emph{population} x \emph{allele}] matrices of the observed allele (haplotype) counts at each locus. Two examples of frequency format files are attached in this package.
}
\references{
Nakajima K et al. (2014) Genetic effects of marine stock enhancement: a case study based on the highly piscivorous Japanese Spanish mackerel. \emph{Canadian Journal of Fisheries and Aquatic Sciences}, 71, 301-314.\cr\cr
Kitada S, Kitakado T, Kishino H (2007) Empirical Bayes inference of pairwise FST and its distribution in the genome. \emph{Genetics}, 177, 861-873.\cr\cr
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.
}
\seealso{
\code{\link[FinePop]{read.genepop}}, \code{\link[FinePop]{read.frequency}}
}
\examples{
data(data.JSM)
cat(data.JSM$MS.genepop, file="JSM_MS_genepop.txt", sep="\n")
cat(data.JSM$MS.freq, file="JSM_MS_freq.txt", sep="\n")
cat(data.JSM$mtDNA.freq, file="JSM_mtDNA_freq.txt", sep="\n")
cat(data.JSM$popname, file="JSM_popname.txt", sep=" ")

# See four text files in your working directory.
#  JSM_MS_genepop.txt  : GENEPOP format file of microsatellite data
#  JSM_MS_freq.txt     : frequency format file of microsatellite data
#  JSM_mtDNA_freq.txt  : frequency format file of mtDNA D-loop region data
#  JSM_popname.txt     : plain text file of subpopulation names
}
\keyword{datasets}
