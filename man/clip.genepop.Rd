\name{clip.genepop}
\alias{clip.genepop}
\title{Remove designated markers from a GENEPOP file.}
\description{
This function reads a GENEPOP file (Rousset 2008), remove designated markers, and write a GENEPOP file of clipped data. The user can directly designate the names of the markers to be removed. The user also can set the filtering threshold of major allele frequency. 
}
\usage{
clip.genepop(infile, outfile, remove.list = NULL, major.af = NULL)
}
\arguments{
  \item{infile}{A character value specifying the name of the GENEPOP file to be clipped.}
  \item{outfile}{A character value specifying the name of the clipped GENEPOP file.}
  \item{remove.list}{A character value or vector specifying the names of the markers to be removed. The names must be included in the target GENEPOP file.}
  \item{major.af}{A numeric value specifying the threshold of major allele frequency for marker removal. Markers with major allele frequencies higher than this value will be removed. This value must be between 0 and 1.}
}
\references{
Rousset F (2008) Genepop'007: a complete reimplementation of the Genepop software for Windows and Linux. \emph{Mol. Ecol. Resources}, 8, 103-106.
}
\author{Reiichiro Nakamichi}
\examples{
# Example of GENEPOP file
data(jsmackerel)
cat(jsmackerel$MS.genepop, file="JSM_MS_genepop.txt", sep="\n")

# remove markers designated by their names
clip.genepop("JSM_MS_genepop.txt", "JSM_genepop_clipped_name.txt",
             remove.list=c("Sni21","Sni26"))

# remove markers with high major allele frequencies (in this example, > 0.5)
clip.genepop("JSM_MS_genepop.txt", "JSM_genepop_clipped_maf.txt",
             major.af=0.5)

# remove markers both by their names and by major allele frequencies
clip.genepop("JSM_MS_genepop.txt", "JSM_genepop_clipped_both.txt",
             remove.list=c("Sni21","Sni26"), major.af=0.5)
}
