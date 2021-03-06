% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inspect_balance.R
\name{inspect_balance}
\alias{inspect_balance}
\alias{inspect_balance.formula}
\alias{print.inspect_balance}
\alias{summary.inspect_balance}
\title{Inspect balance of covariates.}
\usage{
\method{inspect_balance}{formula}(formula, data, method = "dev",
  nPerm = NULL, midpval = TRUE, na.rm = FALSE, treatLevel = NULL, ...)

\method{print}{inspect_balance}(x, ...)

\method{summary}{inspect_balance}(object, ...)
}
\arguments{
\item{formula}{A formula containing an indicator of treatment assignment on
the left hand side and covariates on the right. The treatment indicator
should be numeric with 0/1 values.}

\item{data}{A data frame in which to interpret the variables named in the
formula.}

\item{method}{The method used to compute a p-value associated with the balance
test. See details.}

\item{nPerm}{The number of random permutations of the treatment assignment.
Only applicable with methods \code{"pdev"} and \code{"paic"}.}

\item{midpval}{Should the mid p-value be used?}

\item{na.rm}{Should observations with NAs on any variables named in the RHS of
formula be removed from the covariate balance summary table?}

\item{treatLevel}{A character string for the treatment level of interest. By
default, the treatment is coerced to a factor and the last level is used as
the \code{treatLevel}. This argument is only relevant for calculating the
standardized bias of covariates.}

\item{\dots}{Additional arguments passed to the various methods. Specifically,
for methods \code{"dev"} and \code{"pdev"}, arguments are passed to
\code{stats::glm}. For \code{"paic"}, arguments are passed to
\code{brglm::brglm}, and for \code{"hansen"} they are passed to
\code{RItools::xBalance}.}

\item{x}{A \code{inspect_balance} object.}

\item{object}{A \code{inspect_balance} object.}
}
\value{
An object of class \code{inspect_balance}, which is a list with the
 following components \itemize{ \item \code{fit} The fitted model object
 (\code{NULL} for \code{method ="hansen"}). \item \code{pvalue} The p-value
 of the test. \item \code{nObs} The number of observations used by the
 procedure. \item \code{cbs} The covariate balance summary table. \item
 \code{pdata} The underlying data used in \code{ggplot.inspect_balance}.
 \item \code{treatLevel} The treatment level of interest. \item \code{yLabel}
 The name of the treatment indicator. \item \code{call} The call to
 \code{inspect_balance}.}
}
\description{
\code{inspect_balance} calculates standardized differences for each covariate
between two treatment levels, and tests for conditional independence between
the treatment and the covariates.
}
\details{
In randomized experiments, the assignment of subjects to treatment and control
groups is independent of their baseline covariates. As the sample size grows,
random assignment tends to balance covariates, in the sense that both groups
have similar distributions of covariates. Following Rosenbaum and Rubin
(1985), we define the standardized bias on a covariate as

\deqn{\frac{\bar{x}_t-\bar{x}_c}{\sqrt{\frac{s_t^2 + s_c^2}{2}}}}

where \eqn{\bar{x_t}} and \eqn{\bar{x_c}} represent the sample means of a
covariate in the treated and control groups, respectively, and \eqn{s_t^2} and
\eqn{s_c^2} reresent their sample variances.

Another way to think about balance is that covariates \eqn{X} should have no
predictive power for treatment assignment \eqn{Z}. That is, \eqn{Prob(Z|X) =
Prob(Z)}. Logistic regression is well suited for this task. If \code{method =
"dev"} (default), we follow the approach suggested by Imai (2005). First
regress treatment assignment \eqn{Z} on the covariates \eqn{X} and a constant,
then on a constant alone, and then compare the two fits using a standard
asymptotic likelihood-ratio test. This test is likely to perform poorly (i.e.,
high Type I error rates) in small samples (see Hansen, 2008). If \code{method
= "pdev"}, we compute a permutation distribution of the likelihood ratio
statistic between the two models and compare it to the observe test statistic
to obtain a p-value. Models are fitted using standard logistic regression. If
\code{method = "paic"}, the test statistic is given by the difference in AIC
between the two models. A permutation distribution of this test statistic is
computed and compared to its observed value to obtain a p-value for the test.
Models are fitted using penalized likelihood using Jeffreys prior (Firth,
1993). Finally, if \code{method = "hansen"}, p-values are computed using
\code{RItools::xBalance} (Hansen, 2008).

We note that balance tests of this kind are subject to criticism, since
balance is a characteristic of the sample, not some hypothetical population
(see Ho et al., 2007).
}
\examples{
set.seed(343)
df <- sim_uplift(n = 200, p = 50, response = "binary")
df$T <- ifelse(df$T == 1, 1, 0)
ib <- inspect_balance(T~ X1 + X2 + X3, data = df, method ="pdev", nPerm = 500)
ib
summary(ib)
}
\references{
Firth, D. (1993). "Bias reduction of maximum likelihood estimates". Biometrika
80, pp.27-38

Hansen, B.B. and Bowers, J. (2008)."Covariate Balance in Simple, Stratified
and Clustered Comparative Studies". Statistical Science, 23, pp.219--236.

Ho, D., Kosuke I., King, G. and Stuart, E. (2007). "Matching as Nonparametric
Preprocessing for Reducing Model Dependence in Parametric Causal Inference".
Political Analysis 15, pp.199--236.

Kosuke, I. (2005). "Do Get-Out-The-Vote Calls Reduce Turnout? The Importance
of Statistical Methods for Field Experiments". American Political Science
Review, Vol. 99, No. 2 (May), pp. 283--300.

Rosenbaum, P.R. and Rubin, D.B. (1985). "Constructing a control group using
multivariate matched sampling methods that incorporate the propensity score".
The American Statistician, 39, pp.33--38.
}
\seealso{
\code{\link{ggplot.inspect_balance}}.
}
\author{
Leo Guelman \email{leo.guelman@gmail.com}
}
