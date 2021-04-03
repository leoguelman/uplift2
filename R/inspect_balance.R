inspect_balance <- function(x, ...) UseMethod("inspect_balance")

inspect_balance.default <- function(x, ...) stop("uplift: 'x' should be a formula object")

#'Inspect balance of covariates.
#'
#'\code{inspect_balance} calculates standardized differences for each covariate
#'between two treatment levels, and tests for conditional independence between
#'the treatment and the covariates.
#'
#'In randomized experiments, the assignment of subjects to treatment and control
#'groups is independent of their baseline covariates. As the sample size grows,
#'random assignment tends to balance covariates, in the sense that both groups
#'have similar distributions of covariates. Following Rosenbaum and Rubin
#'(1985), we define the standardized bias on a covariate as
#'
#'\deqn{\frac{\bar{x}_t-\bar{x}_c}{\sqrt{\frac{s_t^2 + s_c^2}{2}}}}
#'
#'where \eqn{\bar{x_t}} and \eqn{\bar{x_c}} represent the sample means of a
#'covariate in the treated and control groups, respectively, and \eqn{s_t^2} and
#'\eqn{s_c^2} reresent their sample variances.
#'
#'Another way to think about balance is that covariates \eqn{X} should have no
#'predictive power for treatment assignment \eqn{Z}. That is, \eqn{Prob(Z|X) =
#'Prob(Z)}. Logistic regression is well suited for this task. If \code{method =
#'"dev"} (default), we follow the approach suggested by Imai (2005). First
#'regress treatment assignment \eqn{Z} on the covariates \eqn{X} and a constant,
#'then on a constant alone, and then compare the two fits using a standard
#'asymptotic likelihood-ratio test. This test is likely to perform poorly (i.e.,
#'high Type I error rates) in small samples (see Hansen, 2008). If \code{method
#'= "pdev"}, we compute a permutation distribution of the likelihood ratio
#'statistic between the two models and compare it to the observe test statistic
#'to obtain a p-value. Models are fitted using standard logistic regression. If
#'\code{method = "paic"}, the test statistic is given by the difference in AIC
#'between the two models. A permutation distribution of this test statistic is
#'computed and compared to its observed value to obtain a p-value for the test.
#'Models are fitted using penalized likelihood using Jeffreys prior (Firth,
#'1993). Finally, if \code{method = "hansen"}, p-values are computed using
#'\code{RItools::xBalance} (Hansen, 2008).
#'
#'We note that balance tests of this kind are subject to criticism, since
#'balance is a characteristic of the sample, not some hypothetical population
#'(see Ho et al., 2007).
#'
#'@name inspect_balance
#'
#'@aliases inspect_balance
#'
#'@export inspect_balance inspect_balance.default inspect_balance.formula
#'
#'@param formula A formula containing an indicator of treatment assignment on
#'  the left hand side and covariates on the right. The treatment indicator
#'  should be numeric with 0/1 values.
#'@param data A data frame in which to interpret the variables named in the
#'  formula.
#'@param method The method used to compute a p-value associated with the balance
#'  test. See details.
#'@param nPerm The number of random permutations of the treatment assignment.
#'  Only applicable with methods \code{"pdev"} and \code{"paic"}.
#'@param midpval Should the mid p-value be used?
#'@param na.rm Should observations with NAs on any variables named in the RHS of
#'  formula be removed from the covariate balance summary table?
#'@param treatLevel A character string for the treatment level of interest. By
#'  default, the treatment is coerced to a factor and the last level is used as
#'  the \code{treatLevel}. This argument is only relevant for calculating the
#'  standardized bias of covariates.
#'@param \dots Additional arguments passed to the various methods. Specifically,
#'  for methods \code{"dev"} and \code{"pdev"}, arguments are passed to
#'  \code{stats::glm}. For \code{"paic"}, arguments are passed to
#'  \code{brglm::brglm}, and for \code{"hansen"} they are passed to
#'  \code{RItools::xBalance}.
#'@param x A \code{inspect_balance} object.
#'@param object A \code{inspect_balance} object.
#'
#'@return An object of class \code{inspect_balance}, which is a list with the
#'  following components \itemize{ \item \code{fit} The fitted model object
#'  (\code{NULL} for \code{method ="hansen"}). \item \code{pvalue} The p-value
#'  of the test. \item \code{nObs} The number of observations used by the
#'  procedure. \item \code{cbs} The covariate balance summary table. \item
#'  \code{pdata} The underlying data used in \code{ggplot.inspect_balance}.
#'  \item \code{treatLevel} The treatment level of interest. \item \code{yLabel}
#'  The name of the treatment indicator. \item \code{call} The call to
#'  \code{inspect_balance}.}
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@seealso \code{\link{ggplot.inspect_balance}}.
#'
#'@references
#'
#'Firth, D. (1993). "Bias reduction of maximum likelihood estimates". Biometrika
#'80, pp.27-38
#'
#'Hansen, B.B. and Bowers, J. (2008)."Covariate Balance in Simple, Stratified
#'and Clustered Comparative Studies". Statistical Science, 23, pp.219--236.
#'
#'Ho, D., Kosuke I., King, G. and Stuart, E. (2007). "Matching as Nonparametric
#'Preprocessing for Reducing Model Dependence in Parametric Causal Inference".
#'Political Analysis 15, pp.199--236.
#'
#'Kosuke, I. (2005). "Do Get-Out-The-Vote Calls Reduce Turnout? The Importance
#'of Statistical Methods for Field Experiments". American Political Science
#'Review, Vol. 99, No. 2 (May), pp. 283--300.
#'
#'Rosenbaum, P.R. and Rubin, D.B. (1985). "Constructing a control group using
#'multivariate matched sampling methods that incorporate the propensity score".
#'The American Statistician, 39, pp.33--38.
#'
#'
#'@examples
#'set.seed(343)
#'df <- sim_uplift(n = 200, p = 50, response = "binary")
#'df$T <- ifelse(df$T == 1, 1, 0)
#'ib <- inspect_balance(T~ X1 + X2 + X3, data = df, method ="pdev", nPerm = 500)
#'ib
#'summary(ib)

inspect_balance.formula <- function(formula,
                                    data,
                                    method = "dev",
                                    nPerm = NULL,
                                    midpval = TRUE,
                                    na.rm = FALSE,
                                    treatLevel = NULL,
                                    ...)
{
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  if (!na.rm) mf$na.action <- na.pass
  mf[[1]] <- as.name("model.frame")
  mf <- eval.parent(mf)
  terms <- attr(mf, "terms")
  if (!attr(terms, "response") > 0)
    stop("uplift: formula must specify a treatment variable as the response.")
  y <- model.response(mf)
  if (!is.numeric(y) || !all(unique(y) %in% c(0,1)))
    stop("uplift: LHS of formula should be numeric with 0/1 values.")
  if (any(is.na(y))) stop('uplift: NAs on LHS of formula not allowed.')
  yLabel <- deparse(terms[[2]])
  ocontrasts <- options(contrasts = c(unordered = "contr_identity", ordered = "contr_identity"))
  on.exit(options(ocontrasts))
  X <- model.matrix(terms, mf, contrasts)
  intercept <- which(colnames(X) == "(Intercept)")
  if (length(intercept > 0)) X <- X[, -intercept]
  y <- as.factor(y)
  yl <- levels(y)
  ynl <- length(yl)
  if (is.null(treatLevel)) treatLevel <- yl[ynl]
  treatInd <- which(yl == treatLevel)
  if (treatInd == 1) y <- relevel(y, ref = yl[2])
  df <- data.frame(X, y)
  names(df)[ncol(df)] <- yLabel
  nObs <- dplyr::summarize(dplyr::group_by_(df, yLabel), nObs=n())
  df <- df %>% reshape2::melt(id.var = yLabel) %>%
    dplyr::group_by_(yLabel, "variable") %>%
    dplyr::summarize(mean = mean(value, na.rm = na.rm), sd = sd(value, na.rm = na.rm))
  df <- as.data.frame(data.table::dcast.data.table(data.table::setDT(df),
                                                   as.formula(paste("variable ~", yLabel)), value.var = c("mean", "sd")))
  df <- dplyr::mutate(df,
                      sd.bias = (df[,3] - df[, 2]) / ((df[,4]^2 + df[,5]^2)/2),
                      var.ratio = df[,5]^2 / df[,4]^2)

  if (method %in% c("pdev", "paic")) {
    if (is.null(nPerm)) stop("uplift: argument \"nPerm\" must be supplied with method: ", method)
    respInd <- which(names(data) == formula[[2]])
  } else {
    if (!is.null(nPerm)) warning("uplift: argument \"nPerm\" ignored with method: ", method)
  }
  if (method == "hansen") {
    fit0 <- NULL
    pvalue <- as.numeric((RItools::xBalance(fmla = formula,
                                            data = data, report = "all", ...))$overall[3])
  } else if (method == "dev") {
    fit0 <- glm(formula = formula, data = data, family = binomial, ...)
    pvalue <- 1-pchisq(abs(fit0$null.deviance - fit0$deviance),
                       abs(fit0$df.null - fit0$df.residual))
  } else {
    if (method == "pdev") {
      fit0 <- glm(formula = formula, data = data, family = binomial, ...)
      stat0 <- abs(fit0$null.deviance - fit0$deviance)
    } else if (method == "paic") {
      fit0 <- brglm::brglm(formula = formula, data = data, family = binomial, ...)
      #(fit0$null.deviance + 2) = aic from null model
      stat0 <- abs(fit0$null.deviance + 2 - fit0$aic)
    }
    stat <- vapply(1:nPerm, dot_every(100,
                                      function(i) get_perm_stat(formula = formula, data = data,
                                                                id = respInd, method = method, ...)),
                   numeric(1))
    mult <- if (midpval) .5 else 1
    lower <- sum(stat0 > stat) / (nPerm + 1)
    equal <- sum(stat0 == stat) / (nPerm + 1)
    pvalue <- lower + mult * equal
  }

  lsOut <- list(fit = fit0,
                pvalue = pvalue,
                nObs = nObs,
                cbs = df, # covaiates balance summary
                pdata = mf, # data for plots
                method = method,
                treatLevel = treatLevel,
                yLabel = yLabel,
                call = call)
  class(lsOut) <- "inspect_balance"
  lsOut
}

get_perm_stat <- function(formula, data, id, method, ...) {
  data[, id] <- sample(data[, id])
  if (method == "pdev") {
    fitPerm <- glm(formula = formula, data = data, family = binomial, ...)
    stat <- abs(fitPerm$null.deviance - fitPerm$deviance)
  } else if (method == "paic") {
    fitPerm <- brglm::brglm(formula = formula, data = data, family = binomial, ...)
    stat <- abs(fitPerm$null.deviance + 2 - fitPerm$aic)
  }
  stat
}

dot_every <- function(m, f) {
  i <- 1
  function(...) {
    if (i %% m == 0) cat(".")
    i <<- i + 1
    f(...)
  }
}

contr_identity <- function(n, contrasts) {
  contr.treatment(n, contrasts = FALSE)
}

#'@rdname inspect_balance
#'@method print inspect_balance
#'@export print.inspect_balance

print.inspect_balance <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("inspect_balance", "\n", sep="")
  cat("\n")
  cat("Number of observations used: ", sum(x$nObs$nObs), "\n", sep="")
  cat("Treatment indicator: ", x$yLabel, "\n", sep="")
  cat("treatLevel: ", x$treatLevel, "\n", sep="")
  cat("Test method: ", x$method, "\n", sep="")
  cat("p-value: ", x$pvalue, "\n", sep="")
  invisible(x)
}

#' @rdname inspect_balance
#' @method summary inspect_balance
#' @export summary.inspect_balance print.summary.inspect_balance

summary.inspect_balance <- function(object, ...) {
  res <- list (call = object$call,
               cbs = object$cbs)
  class(res) <- "summary.inspect_balance"
  return(res)
}

print.summary.inspect_balance <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$cbs)
  invisible(x)
}

#'Plot a `inspect_balance' object.
#'
#'ggplot method for class \code{inspect_balance}.
#'
#'@param x An object of class \code{inspect_balance}.
#'@param i.var A numeric vector of indices of the variables to plot. The
#'  variables should be indexed in the same order that they appear in the
#'  initial inspect_balance formula. The default is to plot all variables.
#'@param n.type The type of plot for numeric variables. Boxplots are generated
#'  by default. Alternatively, use \code{n.type = "density"} for density plots,
#'  or \code{n.type = "qqplot"} for qq-plots.
#'@param f.type The type of plot for categorical variables. 100-percent stacked
#'  columns are generated by default. The alternative option is \code{"counts"},
#'  which show stacked columns of counts.
#'@param legend.position The position of legends \code{("left", "right",
#'  "bottom", "top")}. The default is \code{"right"}.
#'@param nrows Number of rows for plots.
#'@param ncols Number of columns for plots.
#'
#'@export ggplot.inspect_balance
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@examples
#'set.seed(343)
#'df <- sim_uplift(n = 200, p = 50, response = "binary")
#'df$T <- ifelse(df$T == 1, 1, 0)
#'ib <- inspect_balance(T~ X1 + X2 + X3, data = df, method ="pdev", nPerm = 500)
#'ggplot(ib)
#'ggplot(ib, n.type = "density")
#'ggplot(ib, i.var = c(1,2), ncols =2, n.type ="density")


ggplot.inspect_balance <-  function(x,
                                    i.var = NULL,
                                    n.type = "boxplot",
                                    f.type = "percent",
                                    legend.position = NULL,
                                    nrows = NULL,
                                    ncols = NULL)
{

  if (!inherits(x, "inspect_balance"))
    stop("uplift: x must be a 'inspect_balance' class object")
  n.validType <- charmatch(n.type,
                           c("boxplot", "density", "qqplot"))
  f.validType <- charmatch(f.type,
                           c("percent", "counts"))
  if (any(is.na(n.validType)))
    stop("uplift: not recognized plot 'n.type' argument")
  if (any(is.na(f.validType)))
    stop("uplift: not recognized plot 'f.type' argument")

  pdata <- x$pdata
  pdataC <- attributes(attr(pdata, "terms"))$dataClasses[-1]
  pdataL <- attributes(attr(pdata, "terms"))$term.labels
  yLabel <- x$yLabel
  pdata[, 1] <- as.factor(pdata[, 1])

  nVars <- length(pdataL)

  if (is.null(i.var)) {
    i.var <- 1:nVars
  } else if (any(i.var > nVars)) {
    stop("uplift: \"i.var\" must be between 1 and ",nVars)
  }
  if (is.null(legend.position)) legend.position <- "right"

  plots <- vector("list", length(i.var))
  for (i in seq_along(i.var)) {

    if (pdataC[i.var[i]] == "numeric") {
      if (n.type == "boxplot") {
        plots[[i]] <- ggplot(pdata, aes_string(x=yLabel, y =  pdataL[i.var[i]], fill=yLabel)) +
          geom_boxplot(width = 0.3) + guides(fill=FALSE) +
          theme(legend.position=legend.position)
      } else if (n.type == "density") {
        plots[[i]] <- ggplot(pdata, aes_string(x=pdataL[i.var[i]], fill=yLabel)) +
          geom_density(alpha=0.25) +
          theme(legend.position=legend.position)
      } else {
        ### get qq-plot data
        qplotd <- as.data.frame(stats::qqplot(x = pdata[pdata[yLabel] == 0, i.var[i]+1],
                                              y = pdata[pdata[yLabel] == 1, i.var[i]+1],
                                              plot.it = FALSE))
        limits <- range(qplotd$x, qplotd$y)
        plots[[i]] <- ggplot(qplotd) + geom_point(aes(x=x, y=y)) + xlab(paste0(yLabel, " = 0")) +
          scale_x_continuous(limits = limits) +
          scale_y_continuous(limits = limits) +
          geom_abline(aes(intercept=0, slope=1), col = "red") +
          ylab(paste0(yLabel, " = 1")) + ggtitle(pdataL[i.var[i]])
      }
    } else {
      if (f.type == "percent") {
        plots[[i]] <- ggplot(pdata, aes_string(yLabel)) +
          geom_bar(aes_string(fill = pdataL[i.var[i]]), position = 'fill')  +
          scale_y_continuous(labels = percent_format()) + ylab("percent") +
          theme(legend.position=legend.position)
      } else {
        plots[[i]] <- ggplot(pdata, aes_string(yLabel)) +
          geom_bar(aes_string(fill = pdataL[i.var[i]]))  +
          theme(legend.position=legend.position)
      }
    }
  }

  nPlots <- length(plots)
  if (nPlots <= 2) {
    if (is.null(nrows)) nrows <- 1
    if (is.null(ncols)) ncols <- 1
  } else {
    if (is.null(nrows)) nrows <- 2
    if (is.null(ncols)) ncols <- 2
  }
  if (nPlots <= nrows * ncols) {
    gridExtra::grid.arrange(grobs = plots, ncol = ncols, nrow = nrows)
  } else {
    gridExtra::marrangeGrob(grobs = plots, ncol = ncols, nrow = nrows)
  }
}

