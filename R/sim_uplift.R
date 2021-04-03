#'Uplift simulations.
#'
#'Numerical simulations for uplift modeling, as described in Tian et al. (2014).
#'
#'For the gaussian case, \code{sim_uplift} simulates data according to the
#'following specification:
#'
#'
#'\deqn{y = (\beta_0 + \sum_{j=1}^p \beta_{j}X_{j})^2 + (\gamma_0 + \sum_{j=1}^p
#'\gamma_{j}X_{j} + 0.8 X_1 X_2) T + \sigma_{0}\epsilon}
#'
#'where the covariates \eqn{(X_{1}, \ldots, X_{p})} follow a mean zero
#'multivariate normal distribution with a compound symmetric variance-covariance
#'matrix, \eqn{(1-\rho)\mathbf{I}_{p} +\rho \mathbf{1}^{'}\mathbf{1}},
#'\eqn{\beta_0} = \code{beta.par}^-1, \eqn{\beta_j} =  (2 * \code{beta.par})^-1,
#'\eqn{\gamma_0 = 0.4}, \eqn{\gamma_j = (0.8, -0.8, 0.8, -0.8, 0, \ldots, 0)},
#'\eqn{T=[-1,1]} is the treatment indicator gerated with equal probability at
#'random, \eqn{\epsilon} is \eqn{N(0,1)}, and \eqn{\sigma_{0}} = \code{sigma0}.
#'
#'For the binary case,
#'
#'\deqn{y = I((\beta_0 + \sum_{j=1}^p \beta_{j}X_{j})^2 + (\gamma_0 +
#'\sum_{j=1}^p \gamma_{j}X_{j} + 0.8 X_1 X_2) T + \sigma_{0}\epsilon \ge 0)}
#'
#'For further details, see Tian et al. (2014).
#'
#'@name sim_uplift
#'
#'@aliases sim_uplift
#'
#'@export sim_uplift
#'
#'@param n The number of observations.
#'@param p The number of predictors.
#'@param rho The correlation between predictors.
#'@param beta.par Size of main effects. See details.
#'@param sigma0 Multiplier of error term. See details.
#'@param response The type of response distribution. Possible values are
#'  \code{"gaussian"} and \code{"binary"}.
#'
#'@return A data frame including the response variable (\code{y}), the treatment
#'  indicator (\code{T}), the "true" uplift effect (\code{trueUplift}), and the
#'  predictors (\code{X}).
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@references
#'
#'Tian, L., Alizadeh, A., Gentles, A. and Tibshirani, R. (2014). "A simple
#'method for detecting interactions between a treatment and a large number of
#'covariates." Journal of the American Statistical Association, 109:508, pp.
#'1517--1532.
#'
#'
#' @examples
#'
#'set.seed(324)
#'df1 <- sim_uplift(p = 30, response = "binary")
#'str(df1)
#'df2 <- sim_uplift(n = 10000, p = 20)
#'str(df2)


sim_uplift <- function(n = 1000, p = 20, rho = 0.2,
                       beta.par = sqrt(6), sigma0 =  sqrt(2),
                       response = "gaussian") {
  if (p <= 10)
    stop("uplift: the number predictors must be greater than 10")
  if (rho < 0 | rho > 1)
    stop("uplift: rho must be between 0 and 1")
  if (sigma0 < 0)
    stop("uplift: sigma0 must be equal or greater than 0")
  validresponse <- charmatch(response,
                             c("gaussian", "binary"))
  if (any(is.na(validresponse)))
    stop("uplift: not recognized 'response' argument")
  cov.mat <- matrix(rho, p, p)
  diag(cov.mat) <- 1
  x <- cbind(1, mv_rnorm(n = n, p = p, mu = rep(0, p), s = cov.mat))
  beta <- c(beta.par^-1, 0, 0, rep((2 * beta.par)^-1, 8), rep(0, p-10))
  gamma <- c(.4, .8, -.8, .8, -.8, rep(0, p-4))
  T0 <- sample(c(0, 1), n, TRUE)
  T <- 2 * T0 - 1
  t1 <- (x %*% beta)^2
  t2 <- ((x %*% gamma) + .8 * x[, 2] * x[, 3])
  y <- t1 + t2 * T + sigma0 * rnorm(n)
  if (response == "gaussian") {
    trueUplift <- 1.6 * (.5 + x[, 2] - x[, 3] + x[, 4] - x[, 5] + x[, 2] * x[, 3])
  } else if (response == "binary") {
    y <- as.factor(1 * (y >= 0))
    trueUplift <-  pnorm((-t1 + t2)/sigma0) - pnorm((-t1 - t2)/sigma0)
  }
  df <- data.frame(y, T, trueUplift, x[, -1])
  df
}

mv_rnorm <- function(n, p, mu, s) {
  eS <- eigen(s, symmetric = TRUE)
  ev <- eS$values
  X <- matrix(rnorm(p * n), n)
  X <- drop(mu) + eS$vectors %*% diag(sqrt(pmax(ev, 0)), p) %*% t(X)
  t(X)
}
