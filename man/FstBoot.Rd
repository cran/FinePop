\name{FstBoot}
\alias{FstBoot}
\title{Bootstrap sampler of Fst}
\description{
This function provides bootstrapped estimators of Fst to evaluate the environmental effects on the genetic diversity. See Details.
}
\usage{
FstBoot(popdata, fst.method = "EBFST", bsrep = 100, log.bs = F, locus = F)
}
\arguments{
  \item{popdata}{Genotype data object of populations created by read.genepop function from a GENEPOP file.}
  \item{fst.method}{A character value specifying the Fst estimation method to be used. Currently, "EBFST", "EBGstH", "EBDJ", "GstN", "GstNC", "GstH", "DJ" and "thetaWC.pair" are available.}
  \item{bsrep}{A positive integer value specifying the trial times of bootstrapping.}
  \item{log.bs}{A logical value specifying whether the bootstrapped data of each trial should be saved. If TRUE, GENEPOP format files named "gtdata_bsXXX.txt" (XXX=trial number) are saved in the working directory.}
  \item{locus}{A Logical argument indicating whether locus-specific Fst values should be calculated.}
}
\details{
	FinePop provides a method for regression analyses of the pairwise Fst values against geographical distance and the differences to examine the effect of environmental variables on population differentiation (Kitada et. al 2017). First, \code{\link[FinePop]{FstBoot}} function resamples locations with replacement, and then, we also resample the member individuals with replacement from the sampled populations. It calculates pairwise Fst for each bootstrap sample. Second, \code{\link[FinePop]{FstEnv}} function estimates regression coefficients (\code{\link[stats]{lm}} function) of for the Fst values for each iteration. It then computes the standard deviation of the regression coefficients, Z-values and P-values of each regression coefficient. All possible model combinations for the environmental explanatory variables were examined, including their interactions. The best fit model with the minimum information criterion (TIC, Takeuchi 1976, Burnham & Anderson 2002) is selected. Performance for detecting environmental effects on population structuring is evaluated by the R2 value.
}
\value{
\item{bs.pop.list}{List of subpopulations in bootstrapped data}
\item{bs.fst.list}{List of genome-wide pairwise Fst matrices for bootstrapped data.}
\item{org.fst}{Genome-wide pairwise Fst matrix for original data.}
\item{bs.fst.list.locus}{List of locus-specific pairwise Fst matrices for bootstrapped data. (If locus = TRUE)}
\item{org.fst.locus}{Locus-specific pairwise Fst matrix for original data. (If locus = TRUE)}
}
\references{
Burnham KP, Anderson DR (2002) \emph{Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach.} Springer, New York.\cr\cr
Kitada S, Nakamichi R, Kishino H (2017) The empirical Bayes estimators of fine-scale population structure in high gene flow species. \emph{Mol. Ecol. Resources}, DOI: 10.1111/1755-0998.12663. \cr\cr
Takeuchi K (1976) Distribution of information statistics and criteria for adequacy of Models. \emph{Mathematical Science}, 153, 12-18 (in Japanese).
}
\author{Reiichiro Nakamichi, Hirohisa Kishino, Shuichi Kitada}
\seealso{
\code{\link[FinePop]{read.genepop}}, \code{\link[FinePop]{FstEnv}},
\code{\link[FinePop]{EBFST}}, \code{\link[FinePop]{EBGstH}},\code{\link[FinePop]{EBDJ}},
\code{\link[FinePop]{GstN}}, \code{\link[FinePop]{GstNC}}, \code{\link[FinePop]{GstH}},
\code{\link[FinePop]{DJ}}, \code{\link[FinePop]{thetaWC.pair}},
\code{\link[FinePop]{herring}}
}
\examples{
# Example of genotypic and environmental dataset
data(herring)

# Data bootstrapping and Fst estimation
# fstbs <- FstBoot(herring$popdata)

# Effects of environmental factors on genetic differentiation
# fstenv <- FstEnv(fstbs, herring$environment, herring$distance)

# Since these calculations are too heavy, pre-caluculated results are included in this dataset.
fstbs <- herring$fst.bootstrap
fstenv <- herring$fst.env
summary(fstenv)
}
