\name{clip.genepop.maf}
\alias{clip.genepop.maf}
\title{Remove markers with too high major allele frequency from a GENEPOP file.}
\description{
This function reads a GENEPOP file (Rousset 2008), remove markers with too high major allele frequencies, and write a GENEPOP file of clipped data.
}
\usage{
clip.genepop.maf(infile, outfile, major.af)
}
\arguments{
  \item{infile}{A character value specifying the name of the GENEPOP file to be clipped.}
  \item{outfile}{A character value specifying the name of the clipped GENEPOP file.}
  \item{major.af}{A numeric value specifying the threshold of major allele frequency for marker removal. Markers with major allele frequencies higher than this value will be removed. This value must be between 0 and 1.}
}
\references{
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. Mol. Ecol. Resources 8, 103-106.
}
\author{Reiichiro Nakamichi}
\examples{
data(data.JSM)
cat(data.JSM$genepop, file="JSM_genepop.txt", sep="\n")
clip.genepop.maf("JSM_genepop.txt", "JSM_genepop_clipped.txt", 0.5)
}
