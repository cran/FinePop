\name{clip.genepop.name}
\alias{clip.genepop.name}
\title{Remove designated markers from a GENEPOP file.}
\description{This function reads a GENEPOP file (Rousset 2008), remove designated markers, and write a GENEPOP file of clipped data.}
\usage{
clip.genepop.name(infile, outfile, remove.list)
}
\arguments{
  \item{infile}{A character value specifying the name of the GENEPOP file to be clipped.}
  \item{outfile}{A character value specifying the name of the clipped GENEPOP file.}
  \item{remove.list}{A character value or vector specifying the names of the markers to be removed. The names must be included in the target GENEPOP file.}
}
\references{
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. Mol. Ecol. Resources 8, 103-106.
}
\author{Reiichiro Nakamichi}
\examples{
data(data.JSM)
cat(data.JSM$genepop, file="JSM_genepop.txt", sep="\n")
clip.genepop.name("JSM_genepop.txt", "JSM_genepop_clipped.txt", c("Sni21","Sni26"))
}
