% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utree.R
\name{utree}
\alias{utree}
\alias{utree.formula}
\alias{print.utree}
\alias{nodeprune.utree}
\title{Fitting uplift trees.}
\usage{
\method{utree}{formula}(formula, data, na.action, classLevel = NULL,
  treatLevel = NULL, control = utree_control(...), ...)

\method{print}{utree}(x, ...)

\method{nodeprune}{utree}(x, ...)
}
\arguments{
\item{formula}{A model formula of the form y ~ x1 + ....+ xn + trt(), where
the left-hand side corresponds to the observed response, the right-hand side
corresponds to the predictors, and 'trt' is the special expression to mark
the treatment term. At the moment, \code{utree} only handles binary responses.}

\item{data}{A data frame in which to interpret the variables named in the
formula.}

\item{na.action}{A missing-data filter function.}

\item{classLevel}{A character string for the class of interest. Defaults to
the last level of the factor.}

\item{treatLevel}{A character string for the treatment level of interest.
Defaults to the last level of the treatment factor.}

\item{control}{A list with control parameters, see \code{\link{utree_control}}.}

\item{\dots}{Arguments passed to \code{\link{utree_control}}.}

\item{x}{An object of class \code{"utree"}}
}
\value{
An object of class \code{"utree"}.
}
\description{
\code{utree} implements recursive partitioning for uplift modeling.
}
\details{
Roughly, the algorithm works as follows:\enumerate{ \item For each terminal node in the tree
we test the global null hypothesis of no interaction effect between the treatment indicator and any of
the covariates. Stop if this hypothesis cannot be rejected. Otherwise, select the input variable
with strongest interaction effect. The interaction effect is measured by a p-value corresponding
to an asymptotic or permutation test (Strasser and Weber, 1999) for the partial null hypothesis of independence
between each covariate and a transformed response. Specifically, the response is
transformed so the impact of the covariate on the response has a causal interpretation
for the treatment effect (see details in Guelman et al. 2015)
\item Implement a binary split in the selected input variable.
\item Recursively repeate the two steps above.}

Function \code{nodeprune} is not yet implemented for \code{utree} objects.
}
\examples{

set.seed(1)
df <- sim_uplift(n = 1000, p = 50, response = "binary")
form <- create_uplift_formula(x = names(df)[-c(1:3)], y = "y", trt = "T")
fit <- utree(form, data = df, maxdepth = 3)
fit
}
\references{
Guelman, L., Guillen, M., and Perez-Marin A.M. (2015). "A decision support
framework to implement optimal personalized marketing interventions." Decision
Support Systems, Vol. 72, pp. 24--32.

Hothorn, T., Hornik, K. and Zeileis, A. (2006). "Unbiased recursive partitioning:
A conditional inference framework". Journal of Computational and Graphical Statistics,
15(3): 651--674.

Rzepakowski, Piotr and Jaroszewicz, Szymon. (2011). "Decision trees for uplift modeling
with single and multiple treatments". Knowledge and Information Systems, 32(2) 303--327.

Strasser, H. and Weber, C. (1999). "On the asymptotic theory of permutation statistics".
Mathematical Methods of Statistics, 8: 220--250.

Su, X., Tsai, C.-L., Wang, H., Nickerson, D. M. and Li, B. (2009). "Subgroup Analysis via Recursive Partitioning".
Journal of Machine Learning Research 10, 141--158.
}
\seealso{
\code{\link{plot.utree}}.
}
\author{
Leo Guelman \email{leo.guelman@gmail.com}
}
