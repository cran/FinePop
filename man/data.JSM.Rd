\name{data.JSM}
\alias{data.JSM}
\docType{data}
\title{An example dataset of GENEPOP format.}
\description{An example of a GENEPOP (Rousset 2008) format text for a Japanese Spanish mackerel population (Nakajima et al. 2014). It contains genotypic information of 5 microsatellite markers from 9 subpopulations of 791 individuals. Name list of subpopulations also is attached.}
\usage{data("data.JSM")}
\format{
 $ genepop: GENEPOP format text data.\cr
 $ popname: Names of subpopulations.
}
\references{
Nakajima K et al. (2014) Genetic effects of marine stock enhancement: a case study based on the highly piscivorous Japanese Spanish mackerel. Canadian Journal of Fisheries and Aquatic Sciences. 71, 301-314.\cr\cr
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. Mol. Ecol. Resources 8, 103-106.
}
\examples{
data(data.JSM)
cat(data.JSM$genepop, file="JSM_genepop.txt", sep="\n")
print(data.JSM$popname)
}
\keyword{datasets}
