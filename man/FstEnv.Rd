\name{FstEnv}
\alias{FstEnv}
\title{Regression analysis of environmental factors on genetic differentiation}
\description{
This function provides linear regression analysis for Fst against environmental factors to evaluate the environmental effects on the genetic diversity. See Details.
}
\usage{
FstEnv(fst.bs, environment, distance = NULL)
}
\arguments{
  \item{fst.bs}{Bootstrap samples of pairwise Fst matrices provided by FstBoot function from a GENEPOP file.}
  \item{environment}{A table object of environmental factors. Rows are subpopulations, and columns are environmental factors. Names of subpopulations (row names) must be same as those in fst.bs.}
  \item{distance}{A square matrix of distance among subpopulations (omittable). Names of subpopulations (row/column names) must be same as those in fst.bs.}
}
\details{
	FinePop provides a method for regression analyses of the pairwise Fst values against geographical distance and the differences to examine the effect of environmental variables on population differentiation (Kitada et. al 2017). First, \code{\link[FinePop]{FstBoot}} function resamples locations with replacement, and then, we also resample the member individuals with replacement from the sampled populations. It calculates pairwise Fst for each bootstrap sample. Second, \code{\link[FinePop]{FstEnv}} function estimates regression coefficients (\code{\link[stats]{lm}} function) of for the Fst values for each iteration. It then computes the standard deviation of the regression coefficients, Z-values and P-values of each regression coefficient. All possible model combinations for the environmental explanatory variables were examined, including their interactions. The best fit model with the minimum Takeuchi information criterion (TIC, Takeuchi 1976, Burnham & Anderson 2002) is selected. Performance for detecting environmental effects on population structuring is evaluated by the R2 value.
}
\value{
A list of regression result:
  \item{model}{Evaluated model of environmental factors on genetic differentiation.}
  \item{coefficients}{Estimated coefficient, standard deviation, Z value and p value of each factor.}
  \item{TIC}{Takeuchi information criterion.}
  \item{R2}{coefficient of determination.}
}
\references{
Burnham KP, Anderson DR (2002) \emph{Model Selection and Multimodel Inference: A Practical Information-Theoretic Approach.} Springer, New York.\cr\cr
Kitada S, Nakamichi R, Kishino H (2017) The empirical Bayes estimators of fine-scale population structure in high gene flow species. \emph{Mol. Ecol. Resources}, DOI: 10.1111/1755-0998.12663. \cr\cr
Takeuchi K (1976) Distribution of information statistics and criteria for adequacy of Models. \emph{Mathematical Science}, 153, 12-18 (in Japanese).
}
\author{Reiichiro Nakamichi, Hirohisa Kishino, Shuichi Kitada}
\seealso{
\code{\link[FinePop]{FstBoot}}, \code{\link[FinePop]{read.genepop}},
\code{\link[stats]{lm}}, \code{\link[FinePop]{herring}}
}
\examples{
# Example of genotypic and environmental dataset
data(herring)

# Data bootstrapping and Fst estimation
# fstbs <- FstBoot(herring$popdata)

# Effects of environmental factors on genetic differentiation
# fstenv <- FstEnv(fstbs, herring$environment, herring$distance)

# Since these calculations are too heavy, pre-calculated results are included in this dataset.
fstbs <- herring$fst.bootstrap
fstenv <- herring$fst.env
summary(fstenv)
}