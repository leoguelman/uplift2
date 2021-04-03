niv <- function(x, ...) UseMethod("niv")

niv.default <- function(x, ...) stop("uplift: 'x' should be a formula object")

#'Net Information Value.
#'
#'\code{niv} computes the net information value for each uplift predictor. This
#'can be a helpful exploratory tool to (preliminary) determine the predictive
#'power of each variable for uplift.
#'
#'Given a binary response variable \eqn{y \in (0,1)}, the information value
#'(Siddiqi, 2006) from a predictor \eqn{x} is given by
#'
#'\deqn{IV = \sum_{i=1}^{G} \left (P(x=i|y=1) - P(x=i|y=0) \right) \times WOE_i}
#'
#'where \eqn{G} is the number of groups created from a numeric predictor or
#'levels from a categorical predictor, and \eqn{WOE_i = ln
#'(\frac{P(x=i|y=1)}{P(x=i|y=0)})}.
#'
#'To avoid an undefined WOE, an adjustment factor A is used. Specifically,
#'\eqn{WOE_i = ln(\frac{(N(x=i|y=1)+A)/(N(y=1))}{(N(x=i|y=0)+A)/(N(y=0))})},
#'where \eqn{N} represents observation counts.
#'
#'The net information value (NIV) proposed by Larsen (2009) is a natural
#'extension of the IV for the case of uplift. It is computed as
#'
#'\deqn{NIV = \sum_{i=1}^{G}(P(x=i|y=1)^{T} \times P(x=i|y=0)^{C}  -
#'P(x=i|y=0)^{T} \times P(x=i|y=1)^{C}) \times NWOE_i} where \eqn{NWOE_i =
#'WOE_i^{T} - WOE_i^{C}}, and \eqn{T} and \eqn{C} refer to treatment and control
#'groups, respectively.
#'
#'The adjusted net information value (ANIV) is computed as follows
#'
#'\enumerate{ \item Draw B bootstrap samples from the training data and compute
#'the NIV for each variable in each sample. \item Compute the mean of the NIV
#'(\eqn{NIV_{mean}}) and sd of the NIV (\eqn{NIV_{sd}}) for each variable over
#'the \eqn{B} replications. \item The ANIV for a given variable is computed by
#'subtracting a penalty term from the mean NIV. Specifically, \eqn{ANIV =
#'NIV_{mean} - \frac{NIV_{sd}}{\sqrt{B}}}.}
#'
#'@name niv
#'
#'@aliases niv
#'
#'@export niv niv.default niv.formula
#'
#'@param formula A model formula of the form y ~ x1 + ....+ xn + trt(), where
#'  the left-hand side corresponds to the observed response, the right-hand side
#'  corresponds to the predictors, and 'trt' is the special expression to mark
#'  the treatment term. If the treatment term is not a factor, it is converted to one.
#'  \code{niv} only handles response variables of class factor.
#'@param data A data frame in which to interpret the variables named in the
#'  formula.
#'@param subset Expression indicating which subset of the rows of data should be
#'  included. All observations are included by default.
#'@param na.action A missing-data filter function. Defaults to \code{na.pass}.
#'@param nBins The number of bins created from numeric predictors. The bins are
#'  created based on sample quantiles, with a default value of 10 (deciles).
#'@param continuous Specifies the threshold for when bins should be created from
#'  numeric predictors. If there are less or equal than n (i.e.,
#'  \code{continuous = n}) unique values in the numeric predictor, it is
#'  coverted to a factor without binning. The default is \code{continuous = 4}.
#'@param B The number of bootstraps.
#'@param woeAdj The adjustment factor used to avoid an undefined WOE. The value
#'  should be between [0, 1]. By default \code{woeAdj = 0.5}. See details.
#'@param parallel If \code{TRUE}, computations are performed in parallel,
#'  otherwise they are done sequentially.
#'@param nCore The number of cores used. Default is: number of available
#'  cores-1.
#'@param digitsB Number of digits used in formatting the breaks in numeric
#'  predictors.
#'@param classLevel A character string for the class of interest. Defaults to
#'  the last level of the factor.
#'@param treatLevel A character string for the treatment level of interest.
#'  Defaults to the last level of the treatment factor.
#'@param x A \code{niv} object.
#'@param object A \code{niv} object.
#'
#'@return An object of class \code{niv}, which is a list with the following
#'  components (among others passed to the S3 methods): \itemize{ \item
#'  \code{nwoeData} A list of data frames, one for each variable. The columns
#'  represent: \itemize{ \item \code{y00} the number of non-event records
#'  (response != \code{classLevel}) in the control group (treatment !=
#'  \code{treatLevel}). \item \code{y10} the number of event records (response
#'  == \code{classLevel}) in the control group (treatment != \code{treatLevel}).
#'  \item \code{y01} the number of non-event records in the treatment group
#'  (treatment == \code{treatLevel}). \item \code{y11} the number of event
#'  records in the treatment group. \item \code{py00} proportion of non-event
#'  records in the control group. \item \code{py10} proportion of event records
#'  in the control group. \item \code{py01} proportion of non-event records in
#'  the treatment group. \item \code{py11} proportion of event records in the
#'  treatment group. \item \code{woe0} the control group weight-of-evidence.
#'  \item \code{woe1} the treatment group weight-of-evidence. \item \code{nwoe}
#'  the net weight-of-evidence. \item \code{niv} the net information value. }
#'  The values above are computed based on the entire data. \item \code{nivData}
#'  A data frame with the following columns: niv (the average net information
#'  value for each variable over all bootstrap samples), the penalty term, and
#'  the adjusted net information value.}
#'
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@seealso \code{\link{ggplot.niv}}.
#'
#'@references
#'
#'Larsen, K. (2009). Net lift models. In: M2009 - 12th Annual SAS Data Mining
#'Conference.
#'
#'Siddiqi, N. (2006). Credit Risk Scorecards: Developing and Implementing
#'Intelligent Credit Scoring. Wiley, Hoboken, NJ.
#'
#'@examples
#'set.seed(1)
#'df <- sim_uplift(n = 1000, p = 20, response = "binary")
#'f <- create_uplift_formula(names(df)[-c(1:3)], "y", "T")
#'netInf <- niv(f, data = df, B=10, parallel = FALSE)
#'head(netInf$nivData)

niv.formula <- function(formula,
                        data,
                        subset,
                        na.action = na.pass,
                        nBins = 10,
                        continuous = 4,
                        B = 10,
                        woeAdj = 0.5,
                        parallel = TRUE,
                        nCore = NULL,
                        digitsB = NULL,
                        classLevel = NULL,
                        treatLevel = NULL
)
{
  call <- match.call()
  if (!inherits(formula, "formula"))
    stop("uplift: method is only for formula objects.")
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("formula", "data", "subset", "na.action"),
                names(mf), 0L)
  mf <- mf[c(1L, args)]
  mf$drop.unused.levels <- TRUE
  if (missing(na.action)) mf$na.action <- na.pass
  mf[[1L]] <- quote(stats::model.frame)
  special <- "trt"
  mt <- if (missing(data)) terms(formula, special) else
    terms(formula, special, data = data)
  mf$formula <- mt
  mf <- eval(mf, parent.frame())
  terms <- attr(mf, "terms")
  attr(terms, "intercept") <- 0
  trtVar <- attr(terms, "specials")$trt
  if (!is.null(trtVar))
    .tt. <- mf[, trtVar] else
      stop("uplift: formula must include a treatment variable.")
  y <- model.response(mf, "any")
  yFactor <- is.factor(y)
  if (yFactor) {
    if(is.null(classLevel)) classLevel <- levels(y)[nlevels(y)]
    y <- ifelse(y == classLevel, 1, 0)
  } else stop("uplift: response must be a factor.")
  attrTerms <- attributes(terms)
  yLabel <- names(attrTerms$dataClasses[1])
  XLabels <- attrTerms$term.labels[-(trtVar-1)]
  X <- mf[XLabels]
  ttLabel <- attrTerms$term.labels[trtVar-1]
  ttLabel <- sub("\\).*", "", sub(".*\\(", "", ttLabel))
  .tt.f <- factor(.tt.)
  ttLevels <- levels(.tt.f)
  ttNLevels <- length(ttLevels)
  if (is.null(treatLevel)) treatLevel <- ttLevels[ttNLevels]
  .tt.mod <- ifelse(.tt.f == treatLevel, 1L, 0L)
  dataClasses <- attrTerms$dataClasses[-c(1, trtVar)]
  numVars <- dataClasses == "numeric"
  if(is.null(digitsB)) digitsB <- 6
  if (woeAdj < 0 | woeAdj > 1)
    stop("uplift: woeAdj must be between 0 and 1")
  if (any(numVars)) {
    X[numVars] <- data.frame(lapply(X[numVars],
                                    create_bins_numeric,
                                    nBins = nBins, continuous = continuous,
                                    excludeNA = FALSE, digitsB = digitsB))
  }
  df <- data.frame(y, .tt.mod = .tt.mod, X)
  varNames <- names(df)[3:ncol(df)]
  nwoeData <- lapply(varNames, function(i) create_nwoe_table(df, i, woeAdj = woeAdj))
  dfb <- create_bootstrap(df$y, times = B)
  if (parallel) {
    if (is.null(nCore)) {
      nCore <- parallel::detectCores()-1
    }
    doParallel::registerDoParallel(nCore)
    nivStore <- foreach::foreach(i=1:B, .combine=cbind, .export = 'create_nwoe_table', .packages = 'dplyr') %dopar% {
      dfbi <- df[dfb[, i], ]
      vapply(varNames,
             function(j) create_nwoe_table(dfbi, j, woeAdj = woeAdj, nivOnly = TRUE), numeric(1))
    }
    doParallel::stopImplicitCluster()
  } else {
    nivStore <- foreach::foreach(i=1:B, .combine=cbind) %do% {
      dfbi <- df[dfb[, i], ]
      vapply(varNames,
             function(j) create_nwoe_table(dfbi, j, woeAdj = woeAdj, nivOnly = TRUE), numeric(1))
    }
  }
  nivData <- t(apply(nivStore, 1, sum_funs, B = B))
  row.names(nivData) <- varNames
  nivData <- as.data.frame(nivData)
  nivData$adj.niv <- nivData$niv - nivData$penalty
  ordInd <- order(nivData$adj.niv, decreasing = TRUE)
  nivData <-  nivData[ordInd, ]
  lsOut <- list(nwoeData = nwoeData,
                nivData = nivData,
                nObs = nrow(df),
                nVars = length(varNames),
                classLevel = classLevel,
                treatLevel = treatLevel,
                responseLabel = yLabel,
                treatmentLabel = ttLabel,
                B = B,
                varNames = varNames,
                ordInd = ordInd,
                parallel = parallel,
                nCore = nCore,
                call = call)
  class(lsOut) <- "niv"
  lsOut
}

sum_funs <- function(x, B) c(niv = mean(x, na.rm = TRUE),
                             penalty = sd(x, na.rm = TRUE) / sqrt(B))


create_bins_numeric <- function(x, nBins, continuous, excludeNA, digitsB) {
  nas <- is.na(x)
  xul <- length(unique(x[!nas]))
  if (xul > continuous) {
    xc <-  cut(x, breaks = unique(quantile(x, seq(0, 1, 1/nBins), na.rm = TRUE)),
               include.lowest = TRUE, dig.lab = digitsB)
  } else {
    xc <- factor(x)
  }
  if (!excludeNA && any(nas)) {
    xc <- factor(xc, exclude = NULL)
    levels(xc)[length(levels(xc))] <- "Missing"
  }
  xc
}

create_nwoe_table <- function(data, x, woeAdj, nivOnly = FALSE) {

  dfOut <- data %>%
    dplyr::group_by_(x) %>%
    dplyr::summarise(y00 = sum(y == 0L & .tt.mod == 0L),
                     y10 = sum(y == 1L & .tt.mod == 0L),
                     y01 = sum(y == 0L & .tt.mod == 1L),
                     y11 = sum(y == 1L & .tt.mod == 1L)) %>%
    dplyr::mutate(py00 = y00 / sum(y00),
                  py10 = y10 / sum(y10),
                  py01 = y01 / sum(y01),
                  py11 = y11 / sum(y11),
                  woe0 = ifelse(y00 + y10 == 0L, 0L,
                                ifelse(y00 == 0L | y10 == 0L, log(((y10 + woeAdj) / sum(y10)) / ((y00 + woeAdj) / sum(y00))),
                                       log(py10/py00))),
                  woe1 = ifelse(y01 + y11 == 0L, 0L,
                                ifelse(y01 == 0L | y11 == 0L, log(((y11 + woeAdj) / sum(y11)) / ((y01 + woeAdj) / sum(y01))),
                                       log(py11/py01))),
                  nwoe = (woe1 - woe0),
                  niv = (py11 * py00 - py01 * py10) * nwoe)
  if (nivOnly) {
    sum(dfOut$niv)
  } else {
    as.data.frame(dfOut)
  }
}

#'@rdname niv
#'@method print niv
#'@export print.niv

print.niv <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("niv", "\n", sep="")
  cat("\n")
  cat("Number of observations used: ", x$nObs, "\n", sep="")
  cat("Number of variables: ", x$nVars, "\n", sep="")
  cat("Number of bootstraps: ", x$B, "\n", sep="")
  cat("Response variable: ", x$responseLabel, "\n", sep="")
  cat("Treatment indicator: ", x$treatmentLabel, "\n", sep="")
  cat("classLevel: ", x$classLevel, "\n", sep="")
  cat("treatLevel: ", x$treatLevel, "\n", sep="")
  if (x$parallel)  cat("Number of cores: ", x$nCore, "\n", sep="")

  invisible(x)
}

#' @rdname niv
#' @method summary niv
#' @export summary.niv print.summary.niv

summary.niv <- function(object, ...) {
  res <- list (call = object$call,
               nivData = object$nivData)
  class(res) <- "summary.niv"
  return(res)
}

print.summary.niv<- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$nivData)
  invisible(x)
}


### ggplot method

#'Net weight of evidence plots from a `niv' object.
#'
#'ggplot method for class \code{niv}.
#'
#'@param x An object of class \code{niv}.
#'@param i.var A numeric vector of indices of the variables to plot. The
#'  variables should be indexed in the same order that they appear in the
#'  initial niv formula. The default is to plot all variables.
#'@param same.limits Should limits for the y-axis be the same in all plots?
#'@param col The fill colour of the bars.
#'@param xlab A character string giving the title for the x axis (the length
#'  should be the same as the number of plots).
#'@param ylab A character string of length 1 giving the title for the y axis.
#'@param title The main title for the plot (the length should be the same as the
#'  number of plots).
#'@param nrows Number of rows for plots.
#'@param ncols Number of columns for plots.
#'
#'@export ggplot.niv
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@examples
#'set.seed(1)
#'df <- sim_uplift(n = 1000, p = 20, response = "binary")
#'netInf <- niv(y ~ trt(T) + X1 + X2 + X9 + X10, data = df,
#'              B=10, parallel = FALSE, digitsB = 1)
#'ggplot(netInf)


ggplot.niv <- function(x,
                       i.var = NULL,
                       same.limits = FALSE,
                       col = NULL,
                       xlab = NULL,
                       ylab = NULL,
                       title = NULL,
                       nrows = NULL,
                       ncols = NULL
)
{

  if (!inherits(x, "niv"))
    stop("uplift: x must be a 'niv' class object")
  nwoeData <- x$nwoeData
  nivData <- x$nivData
  ordInd <- x$ordInd
  nivData$rank <- 1:nrow(nivData)
  nivData <- nivData[order(ordInd), ]
  varNames <- x$varNames
  nVars <- length(varNames)
  if (is.null(i.var)) {
    i.var <- 1:nVars
  } else if (any(i.var > nVars)) {
    stop("uplift: \"i.var\" must be between 1 and ",nVars)
  }
  if (is.null(col)) col <- "blue"
  if (is.null(ylab)) ylab <- "NWOE"
  if (is.null(xlab)) {
    xlabs <- varNames[i.var]
  } else xlabs <- xlab
  if (same.limits) ylim <- range(unlist(lapply(nwoeData[i.var], function(x) x$nwoe)))
  nPlots <- length(i.var)
  plotList <- vector("list", nPlots)
  titles <- character(nPlots)
  for (i in seq_along(i.var)) {

    if (is.null(title)) {
      titles[i] <- paste0("adj.niv = ", round(nivData$adj.niv[i.var[i]], 4),
                          " | Rank = ", nivData$rank[i.var[i]])
    } else titles <- title

    plotList[[i]] <- ggplot(data=nwoeData[[i.var[i]]], aes_string(x=varNames[i.var[i]], y="nwoe")) +
      geom_bar(position = "dodge", stat = "identity", fill = col) +
      theme(axis.text.x =element_text(angle = 90, hjust = 1, vjust=0.5)) +
      xlab(xlabs[i]) + ylab(ylab) + ggtitle(titles[i])
    if (same.limits) {
      plotList[[i]] <- plotList[[i]] + ylim(ylim)
    }
  }
  if (nPlots <= 2) {
    if (is.null(nrows)) nrows <- 1L
    if (is.null(ncols)) ncols <- 1L
  } else {
    if (is.null(nrows)) nrows <- 2L
    if (is.null(ncols)) ncols <- 2L
  }
  if (nPlots <= nrows * ncols) {
    gridExtra::grid.arrange(grobs = plotList, ncol = ncols, nrow = nrows)
  } else {
    gridExtra::marrangeGrob(grobs = plotList, ncol = ncols, nrow = nrows)
  }
}

