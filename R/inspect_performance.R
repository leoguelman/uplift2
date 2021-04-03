inspect_performance <- function(x, ...) UseMethod("inspect_performance")

inspect_performance.default <- function(x, ...) stop("uplift: 'x' should be a formula object")

#'Inspect performance from a fitted uplift model.
#'
#'\code{inspect_performance} returns various performance measures from a fitted
#'uplift model, including predicted uplift versus observed response under each
#'treatment, and the Qini coefficient. The results can be used to seamlessly
#'produce Qini curves and calibration plots.
#'
#'The Qini curve (Radcliffe, 2007) is a two-dimensional depiction of model
#'performance for uplift models. It represents a natural extension of the Gains
#'curve (Blattberg et al., 2008, p. 319) for uplift models. It is constructed as
#'follows: (i) observations are sorted by the uplift predictions in descending
#'order, (ii) the cumulative number of responses are computed for each
#'treatment and expressed as a percentage of the total number of observations
#'within each treatment, (iii) the "net lift" is computed as the difference
#'between the values obtained in (ii) for the reference treatment
#'(\code{treatLevel}) and the alternative treatment. The Qini curve is a plot
#'of the net lift versus the corresponding fraction of observations in the data.
#'
#'The interpretation of the Qini curve is as follows: on the x-axis we show the
#'fraction of subjects in the population in which the treatment is performed,
#'and on the y-axis we show the difference in the success rate between the
#'reference treatment and the alternative treatment.
#'
#'A benchmark for a given uplift model can be represented by the strategy of
#'randomly selecting subjects to perform the treatment. This is represented in the
#'figure by the diagonal line. For example, if we perform the treatment on 30 percent of the
#'population, we expect to obtain 30 percent of the net lift relative to
#'performing the action on the entire population.
#'
#'The Qini coefficient is a single estimate of model performance ranging from +100
#'(best possible model), to -100 (worst possible model). A value of zero represents a performance
#'equivalent to that of a random model. The Qini coefficient is computed as
#'
#'\deqn{(AUQC_m - AUQC_r) / (AUQC_o - AUQC_r)}
#'
#'where \eqn{AUQC_m}, \eqn{AUQC_r}, and \eqn{AUQC_o} represent the area under
#'the Qini curve for the fitted uplift model, the random model, and the optimal
#'model, respectively.
#'
#'If \code{method = "quantile"} (the default), 'n' bins (\code{nBins}) are
#'created based on quantiles of the distribution of the predicted values (in an
#'attempt to have the same number of observations in each bin). If \code{method
#'= "bucket"}, bins are created by dividing the predicted values into equally
#'spaced intervals based on the difference between the minimum and maximum
#'values. Unlike \code{method = "quantile"}, the number of observations in each
#'group is typically unequal. If \code{method = "user"}, bins are created
#'according to user-specified breaks (\code{userBreaks}).
#'
#'In some cases, it may not be feasible to obtain the number of bins requested
#'(\code{nBins}) (e.g., due to many ties in the predicted values). The function
#'returns the effective number of bins created for each model
#'(\code{actualBins}).
#'
#'@name inspect_performance
#'
#'@aliases inspect_performance
#'
#'@export inspect_performance inspect_performance.default inspect_performance.formula
#'
#'@param formula A model formula of the form y ~ pred_1 + ....+ pred_n + trt(),
#'  where the left-hand side corresponds to the observed response, and the
#'  right-hand side corresponds to the predicted values from an arbitrary number
#'  of uplift models, and 'trt' is the special expression to mark the treatment
#'  term. If the treatment term is not a factor, it is converted to one.
#'@param data A data frame in which to interpret the variables named in the
#'  formula.
#'@param subset Expression indicating which subset of the rows of data should be
#'  included. All observations are included by default.
#'@param na.action A missing-data filter function.
#'@param method Possible values are \code{"quantile"} (default) if you want to
#'  create intervals from the predicted uplift values with approximately the
#'  same number of observations in each group, \code{"bucket"} if you want to
#'  divide the predicted response values into equally spaced intervals, or
#'  \code{"user"} to create intervals from user-specified breaks (see details
#'  below).
#'@param nBins The number of bins to create.
#'@param qini Return the Qini coefficient for each model?
#'@param qini.nBins The number of cutoffs in the predicted values used for
#'  constructing the qini curves. Cutoffs are created based on quantiles of the
#'  distribution of the predicted values. By default, cutoffs are given by unique
#'  values of the predictions.
#'@param userBreaks A user-specified numeric vector of breaks for the predicted
#'  uplift values from which to create bins. It is required when \code{method =
#'  "user"}, and ignored otherwise.
#'@param classLevel A character string for the class of interest. Only
#'  applicable when the response is a factor. Defaults to the last level of the
#'  factor.
#'@param treatLevel A character string for the treatment level of interest.
#'  Defaults to the last level of the treatment factor.
#'@param x A \code{inspect_performance} object.
#'@param object A \code{inspect_performance} object.
#'@param \dots Additional arguments for the S3 methods.
#'
#'
#'@return An object of class \code{"inspect_performance"}, which is a list with
#'  the following components: \itemize{ \item \code{data} The data for
#'  underlying calibration plots \item \code{method} The method used to create
#'  the bins \item \code{models} The labels of the variables supplied in the
#'  right-hand side of the model formula \item \code{nBins} The number of bins
#'  requested \item \code{actualnBins} The effective number of bins created for
#'  each model \item \code{yFactor} Is response a factor? \item
#'  \code{classLevel} The class of interest \item \code{treatLevel} The
#'  treatment level of interest \item \code{qiniCall} Whether the qini
#'  coefficient was requested in the function call \item \code{qiniData} The
#'  data for plotting the qini curves \item \code{qiniC} The qini coefficient
#'  \item \code{treatInd} The position of the treatment level of interest \item
#'  \code{call} The original call to \code{inspect_performance} }
#'
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@seealso \code{\link{ggplot.inspect_performance}}.
#'
#'@references
#'
#'Blattberg, R. C., Do, K. B., and Scott, N. A. (2008). "Database Marketing:
#'Analyzing and Managing Customers". Springer Science+Business Media, New York,
#'NY.
#'
#'Radcliffe, N (2007). "Using control groups to target on predicted lift:
#'Building and assessing uplift models." Direct Marketing Analytics Journal, An
#'Annual Publication from the Direct Marketing Association Analytics Council:
#'pp. 14--21.
#'
#'
#' @examples
#'
#'set.seed(324)
#'df_train <- sim_uplift(p = 30, response = "binary")
#'df_test <- sim_uplift(n = 10000, p = 30, response = "binary")
# summary(df_test$trueUplift)
#'fit_t1 <- glm(as.formula(paste('y ~', paste('X', 1:30, sep = '', collapse = "+"))),
#'              family = "binomial", data = df_train, subset = T==1)
#'fit_t0 <- glm(as.formula(paste('y ~', paste('X', 1:30, sep = '', collapse = "+"))),
#'              family = "binomial", data = df_train, subset = T==-1)
#'uplift_score <- predict(fit_t1, df_test, type = "response") -
#'                predict(fit_t0, df_test, type = "response")
#'df_test$uplift_score <- uplift_score
#'res <- inspect_performance(y ~ uplift_score + trueUplift + trt(T),
#'                           data = df_test, qini = TRUE)
#'res
#'summary(res)

inspect_performance.formula <- function(formula,
                                        data,
                                        subset,
                                        na.action = na.fail,
                                        method = "quantile",
                                        nBins = 10,
                                        qini = FALSE,
                                        qini.nBins = NULL,
                                        userBreaks = NULL,
                                        classLevel = NULL,
                                        treatLevel = NULL)



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
  if (missing(na.action)) {
    mf$na.action <- na.fail
  } else if (na.action == "na.pass") stop("uplift: passing NAs is not allowed.")
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
    .tt. <- mf[, trtVar, drop = TRUE] else
      stop("uplift: formula must include a treatment variable.")
  y <- model.response(mf, "any")
  yFactor <- is.factor(y)
  if (yFactor) {
    if(is.null(classLevel)) classLevel <- levels(y)[nlevels(y)]
    y <- ifelse(y == classLevel, 1, 0)
  } else {
    if (!is.null(classLevel))
      warning("uplift: classLevel is ignored when response is not a factor.")
    classLevel = "NA"}
  attrTerms <- attributes(terms)
  yLabel <- names(attrTerms$dataClasses[1])
  XLabels <- attrTerms$term.labels[-(trtVar-1)]
  X <- mf[XLabels]
  ttLabel <- attrTerms$term.labels[trtVar-1]
  ttLabel <- sub("\\).*", "", sub(".*\\(", "", ttLabel))
  XClasses <- attrTerms$dataClasses[-c(1, trtVar)]
  if (!all(unique(XClasses) == "numeric"))
    stop("uplift: model predictions must be numeric.")
  validMethod <- c("quantile", "bucket", "user")
  methodOk <- charmatch(method, validMethod, -1)
  if (any(methodOk == -1)) {
    stop(paste("uplift: Invalid option for method:", paste(method[methodOk == -1], collapse=", ")))
  }
  if (method == "user") {
    if (is.null(userBreaks)) {
      stop("uplift: 'userBreaks must be supplied when method = 'user'.") } else {
        nBinsNew <- length(unique(userBreaks)) - 1
        if (nBins != nBinsNew) {
          nBins <- nBinsNew
          warning("uplift: nBins has been reset to ", nBinsNew, ", as implied by the userBreaks.")
        }
      }
  }
  if (!yFactor && length(unique(y)) < 5)
    warning("uplift:", sQuote(yLabel), " has less than 5 unique values.
            Consider treating as factor instead of numeric.")
  df <- data.frame(y, .tt. = factor(.tt.), X)
  dfm <- reshape2::melt(df, id.vars = c("y", ".tt."), variable.name = "models", value.name = "pred")
  levels(dfm$models) <- XLabels # fix issue when predictors are named y or .tt.
  df2 <- plyr::ddply(dfm, plyr::.(models), create_bins, method, nBins, userBreaks)
  dfS2 <- dplyr::summarize(dplyr::group_by(df2, models, bin, .tt.),
                           obsCount = n(),
                           meanObservedResponse = round(mean(y), 6))
  dfS2W1 <- reshape2::dcast(dfS2, models + bin ~ .tt., value.var = "meanObservedResponse")
  dfS2W2 <- reshape2::dcast(dfS2, models + bin ~ .tt., value.var = "obsCount", fill = 0) # As NAs missing zero records
  dfS22 <- dplyr::summarize(dplyr::group_by(df2, models, bin),
                            meanPredictedResponse = round(mean(pred), 6))
  dfOut <- data.frame(dfS22,
                      dfS2W2[, -c(1:2)],
                      dfS2W1[, -c(1:2)])
  dfOut <- dplyr::arrange(dfOut, models, meanPredictedResponse)
  ttLevels <- levels(dfm$.tt.)
  ttNLevels <- length(ttLevels)
  if (is.null(treatLevel)) treatLevel <- ttLevels[ttNLevels]
  treatInd <- which(ttLevels == treatLevel)
  names(dfOut)[4:(4+ttNLevels-1)] <- paste("N_", ttLabel, "=", ttLevels, sep = "")
  names(dfOut)[(4+ttNLevels):(4+2 * ttNLevels-1)] <- paste("meanObservedResponse_", ttLabel, "=", ttLevels, sep = "")
  if (qini) {
    qiniData  <- plyr::ddply(dfm, plyr::.(models), get_qini_data,
                             yLabel,
                             classLevel,
                             treatLevel,
                             ttLabel,
                             treatInd,
                             qini.nBins,
                             yFactor,
                             optimal = FALSE)
    auqcData <- as.data.frame(dplyr::summarize(dplyr::group_by(qiniData, models),
                                               auqc = compute_auqc(N, cumPctnN, netLift)))
    auqcR <- as.data.frame(dplyr::summarize(dplyr::group_by(qiniData, models),
                                            auqc = compute_auqc(N, cumPctnN, gainsRandom)))
    if (yFactor) {
      qiniOptData <- dplyr::mutate(df[, c(1,2)],
                                   models = ".optScore.",
                                   pred = ifelse(y == 1 & .tt. == treatLevel, 4,
                                                 ifelse(y == 0 & .tt. != treatLevel, 3,
                                                        ifelse(y == 0 & .tt. == treatLevel, 2,
                                                               ifelse(y == 1 & .tt. != treatLevel, 1, NA)))))
    } else {
      qiniOptData <- dplyr::mutate(rbind(df[, c(1,2)][df$.tt. == treatLevel, ] %>% dplyr::arrange(desc(y)),
                                         df[, c(1,2)][df$.tt. != treatLevel, ] %>% dplyr::arrange(y)),
                                   models = ".optScore.",
                                   pred = nrow(df):1)
    }
    qiniOptData <- plyr::ddply(qiniOptData, plyr::.(models), get_qini_data,
                               yLabel,
                               classLevel,
                               treatLevel,
                               ttLabel,
                               treatInd,
                               qini.nBins,
                               yFactor,
                               optimal = TRUE)
    auqcOpt <- as.data.frame(dplyr::summarize(dplyr::group_by(qiniOptData, models),
                                              auqc = compute_auqc(N, cumPctnN, netLift)))
    qiniC <- data.frame(models = auqcData[, 1],
                        qiniCoef = 100 * (auqcData[, 2] - auqcR[, 2]) / (auqcOpt[, 2]  - auqcR[, 2]))
  } else {
    qiniData <- NULL
    qiniC <- NULL
  }

  lsOut <- list(data = dfOut,
                method = method,
                models = XLabels,
                nBins = nBins,
                actualnBins = tapply(dfOut$bin, dfOut$models, length),
                yFactor = yFactor,
                classLevel = classLevel,
                treatLevel = treatLevel,
                qiniCall = qini,
                qiniData = qiniData,
                qiniC = qiniC,
                treatInd = treatInd,
                call = call)

  class(lsOut) <- "inspect_performance"
  lsOut

}


trt <- function(x) x

create_bins <- function(x, method, nBins, userBreaks)
{
  Breaks <- switch(method,
                   quantile = unique(quantile(x$pred, probs = seq(0, 1, 1/nBins))),
                   bucket = seq(min(x$pred), max(x$pred), (max(x$pred) - min(x$pred)) / nBins),
                   user = unique(userBreaks))
  dfB <-  data.frame(y = x$y,
                     .tt. = x$.tt.,
                     pred = x$pred,
                     bin = cut(x$pred, breaks = Breaks, include.lowest = TRUE))
  dfB
}

get_qini_data <- function(x,
                          yLabel,
                          classLevel,
                          treatLevel,
                          ttLabel,
                          treatInd,
                          qini.nBins,
                          yFactor,
                          optimal)
{

  nObs <- nrow(x)
  ttNperLevel <- tapply(x$y, x$.tt., length)
  if (!optimal) {
    if (is.null(qini.nBins)) {
      cutoffs <- unique(x[, "pred"])
    } else {
      cutoffs <- unique(quantile(x[, "pred"], probs = seq(0, 1, 1/qini.nBins)))
    }
    x$bin <- cut(x$pred, breaks = cutoffs, labels = FALSE, include.lowest = TRUE)
  } else {
    x$bin <- as.factor(x$pred)
  }
  x <- dplyr::mutate(x,
                     y11 = ifelse(.tt. == treatLevel, y, 0),
                     y12 = ifelse(.tt. != treatLevel, y, 0))
  x <- dplyr::summarize(dplyr::group_by(x, models, bin),
                        N = n(),
                        y11 = sum(y11),
                        y12 = sum(y12))
  x <- dplyr::arrange(x, desc(bin))
  x <- dplyr::mutate(x,
                     cumN = cumsum(N),
                     cumPctnN = cumN / nObs,
                     cumy11 = cumsum(y11),
                     cumy12 = cumsum(y12),
                     gainsL1 = cumy11 / ttNperLevel[treatInd],
                     gainsL2 = cumy12 / sum(ttNperLevel[-treatInd]),
                     netLift = gainsL1 - gainsL2,
                     gainsRandom = netLift[length(netLift)] * cumPctnN)
  label1 <- paste(paste(yLabel, ifelse(yFactor, classLevel, ""), sep = ""), "_",
                  paste(ttLabel, treatLevel, sep = ""), sep ="")
  label2 <- paste(paste(yLabel, ifelse(yFactor, classLevel, ""), sep = ""), "_",
                  paste(ttLabel, "Other", sep = ""), sep ="")
  colnames(x)[c(4, 5, 8, 9, 10, 11)] <- c(paste("sum_", label1, sep =""),
                                          paste("sum_", label2, sep =""),
                                          paste("cumSum_", label1, sep =""),
                                          paste("cumSum_", label2, sep =""),
                                          paste("gains_", label1, sep =""),
                                          paste("gains_", label2, sep =""))
  x <- as.data.frame(x)
  x
}

compute_auqc <- function(N, cumPctnN, netLiftVar) {
  finiteBool <- is.finite(cumPctnN) & is.finite(netLiftVar)
  xVal <- cumPctnN[finiteBool]
  yVal <- netLiftVar[finiteBool]
  nObs <- sum(N)
  # Extrapolate curve to x ~ zero
  if (1/nObs < xVal[1]) {
    xVal <- c(1/nObs, xVal)
    yVal <- c(1/nObs * yVal[1] / xVal[1], yVal)
  }
  K <- length(xVal)
  if (K < 2)
    stop("uplift: Not enough distinct predictions to compute
         the qini coefficient for at least one of the models.")
  auqc <- sum(0.5 * (xVal[2:K] - xVal[1:K-1]) * (yVal[2:K] + yVal[1:K-1]))
  auqc
}


#' @rdname inspect_performance
#' @method print inspect_performance
#' @export print.inspect_performance

print.inspect_performance <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("inspect_performance", "\n", sep="")
  cat("Response type: ", ifelse(x$yFactor, "factor", "numeric"), "\n", sep="")
  cat("Models:", paste(unique(x$models), collapse = ", "), "\n")
  cat("Method: ", x$method, "\n", sep="")
  cat("nBins: ", x$nBins, "\n", sep="")
  cat("Actual Number of Bins: ",
      paste(names(x$actualnBins), "=", x$actualnBins, collapse = ", "), "\n")
  cat("classLevel: ", x$classLevel, "\n", sep="")
  cat("treatLevel: ", x$treatLevel, "\n", sep="")
  cat("Return qini: ", x$qiniCall, "\n", sep="")
  invisible(x)
}

#' @rdname inspect_performance
#' @method summary inspect_performance
#' @export summary.inspect_performance print.summary.inspect_performance

summary.inspect_performance <- function(object, ...) {
  res <- list (call = object$call,
               data = object$data,
               qini = object$qini)
  class(res) <- "summary.inspect_performance"
  return(res)
}

print.summary.inspect_performance <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
  cat("\n")
  print(x$data)
  if (!is.null(x$qini)) {
    cat("\nqini coefficient(s):\n")
    print(x$qini)
  }
  invisible(x)
}

### ggplot method

#'Qini and Calibration plots from a \code{inspect_performance} object.
#'
#'ggplot method for class `inspect_performance'
#'
#'@param x An object of class \code{"inspect_performance"}.
#'@param type The type of plot. Possible values are \code{qini} for Qini curves
#'  (default) and \code{calib} for calibration plots.
#'@param legend.position The position of legends ("left", "right", "bottom",
#'  "top").
#'@param fillCol Panel's background color.
#'@param facets Lay out panels in a grid?
#'@param pointSize For \code{type = "calib"}, make size of points in plot
#'  proportional to the number of observations?
#'@param diagCol For \code{type = "calib"}, the color of the diagonal line.
#'@param xlim,ylim Numeric vectors of length 2, giving the x and y coordinates
#'  ranges.
#'@param xBreaks,yBreaks Points at which x, y gridlines appear.
#'@param xlab,ylab Title for the x, y axes.
#'@param title The main title for the plot.
#'
#'@export ggplot.inspect_performance
#'
#'@import ggplot2
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@examples
#'
#'set.seed(324)
#'df_train <- sim_uplift(p = 30, response = "binary")
#'df_test <- sim_uplift(n = 10000, p = 30, response = "binary")
# summary(df_test$trueUplift)
#'fit_t1 <- glm(as.formula(paste('y ~', paste('X', 1:30, sep = '', collapse = "+"))),
#'              family = "binomial", data = df_train, subset = T==1)
#'fit_t0 <- glm(as.formula(paste('y ~', paste('X', 1:30, sep = '', collapse = "+"))),
#'              family = "binomial", data = df_train, subset = T==-1)
#'uplift_score <- predict(fit_t1, df_test, type = "response") -
#'                predict(fit_t0, df_test, type = "response")
#'df_test$uplift_score <- uplift_score
#'res <- inspect_performance(y ~ uplift_score + trueUplift + trt(T),
#'                           data = df_test, qini = TRUE)
#'res
#'summary(res)
#'ggplot(res)
#'ggplot(res, type = "calib", pointSize = TRUE)

ggplot.inspect_performance <- function(x,
                                       type = "qini",
                                       legend.position = "top",
                                       fillCol = "white",
                                       facets = FALSE,
                                       pointSize = FALSE,
                                       diagCol = "grey",
                                       xlim = NULL,
                                       ylim = NULL,
                                       xBreaks = NULL,
                                       yBreaks = NULL,
                                       xlab = NULL,
                                       ylab = NULL,
                                       title = NULL)
{

  if (!inherits(x, "inspect_performance"))
    stop("uplift: x must be a 'inspect_performance' class object")

  if (type == "qini" && !x$qiniCall)
    stop("uplift: The qini option must be set to 'TRUE' in the call to
         'inspect_performance' to enable qini plots.")

  validType <- charmatch(type, c("qini", "calib"))
  if (any(is.na(validType)))
    stop("uplift: not a valid plot type.")

  if (type == "qini") {

    qini_plot(x = x,
              legend.position = legend.position,
              fillCol = fillCol,
              facets = facets,
              xlim = xlim,
              ylim = ylim,
              xBreaks = xBreaks,
              yBreaks = yBreaks,
              xlab = xlab,
              ylab = ylab,
              title)

  } else if (type == "calib") {

    calib_plot(x = x,
               legend.position = legend.position,
               fillCol = fillCol,
               facets = facets,
               pointSize = pointSize,
               diagCol = diagCol,
               xlim = xlim,
               ylim = ylim,
               xBreaks = xBreaks,
               yBreaks = yBreaks,
               xlab = xlab,
               ylab = ylab,
               title)
  }
}


qini_plot <- function(x,
                      legend.position,
                      fillCol,
                      facets,
                      xlim,
                      ylim,
                      xBreaks,
                      yBreaks,
                      xlab,
                      ylab,
                      title)
{

  if (is.null(xlab)) xlab <- "Proportion of population targeted"
  if (is.null(ylab)) ylab <- "Net lift"
  if (is.null(title)) title <- "Qini plot"
  qiniData <- x$qiniData
  obsPerModel <- dplyr::summarize(dplyr::group_by(qiniData, models),
                                  N = n())
  qiniData <- rbind.data.frame(qiniData[, c(1, 7, 12)],
                               data.frame(models = "Random",
                                          cumPctnN = qiniData[1:as.integer(obsPerModel[1, 2]), 7],
                                          netLift = qiniData[1:as.integer(obsPerModel[1, 2]), 13]))
  yLimits <- range(qiniData$netLift)
  xLimits <- range(qiniData$cumPctnN)
  if (is.null(xlim)) xlim <- xLimits
  if (is.null(ylim)) ylim <- yLimits
  xAxisBreaks <- pretty(xLimits, x$nBins)
  yAxisBreaks <- pretty(yLimits, x$nBins)
  if (is.null(xBreaks)) xBreaks <- xAxisBreaks
  if (is.null(yBreaks)) yBreaks <- yAxisBreaks
  if (facets) {
    plotOut <- ggplot(data = qiniData,
                      aes(x = cumPctnN, y = netLift)) +
      facet_grid(.~ models)
  } else {
    plotOut <- ggplot(data = qiniData,
                      aes(x = cumPctnN, y = netLift,
                          colour = models))
  }
  plotOut <- plotOut +
    geom_line() +
    scale_x_continuous(limits = c(xlim[1], xlim[2]), breaks = xBreaks) +
    scale_y_continuous(limits = c(ylim[1], ylim[2]), breaks = yBreaks) +
    xlab(xlab) + ylab(ylab) + ggtitle(title) +
    theme(legend.position=legend.position,
          legend.key =  element_rect(fill = fillCol),
          panel.background = element_rect(fill = fillCol, colour = 'black'))
  plotOut

}


calib_plot <- function(x,
                       legend.position,
                       fillCol,
                       facets,
                       pointSize,
                       diagCol,
                       xlim,
                       ylim,
                       xBreaks,
                       yBreaks,
                       xlab,
                       ylab,
                       title)
{

  if (is.null(xlab)) xlab <- "Mean predicted uplift"
  if (is.null(ylab)) ylab <- "Mean observed uplift"
  if (is.null(title)) title <- "Calibration plot"
  xdata <- x$data
  refTreatName1 <- names(xdata)[6:ncol(xdata)][x$treatInd]
  refTreatName2 <- names(xdata)[6:ncol(xdata)][-x$treatInd][1] #if more than 1 treatment, take 1 level as the ref.
  xdata <- dplyr::mutate(xdata,
                         obsCount = xdata[, 4] + xdata[, 5],
                         obsUplift = xdata[, refTreatName1] - xdata[,  refTreatName2])
  limits <- range(xdata$meanPredictedResponse, xdata$obsUplift)
  if (is.null(xlim)) xlim <- limits
  if (is.null(ylim)) ylim <- limits
  axisBreaks <- pretty(limits, x$nBins)
  if (is.null(xBreaks)) xBreaks <- axisBreaks
  if (is.null(yBreaks)) yBreaks <- axisBreaks
  if (facets) {
    plotOut <- ggplot(data = xdata,
                      aes(x = meanPredictedResponse, y = obsUplift)) +
      facet_grid(.~ models)
  } else {
    plotOut <- ggplot(data = xdata,
                      aes(x = meanPredictedResponse, y = obsUplift,
                          colour = models))
  }
  plotOut <- plotOut +
    geom_line() +
    geom_abline(intercept = 0, slope = 1, col = diagCol) +
    scale_x_continuous(limits = c(xlim[1], xlim[2]), breaks = xBreaks) +
    scale_y_continuous(limits = c(ylim[1], ylim[2]), breaks = yBreaks) +
    xlab(xlab) + ylab(ylab) + ggtitle(title) +
    theme(legend.position=legend.position,
          legend.key = element_rect(fill = fillCol),
          panel.background = element_rect(fill = fillCol, colour = 'black'))
  if (pointSize) {
    plotOut <- plotOut + geom_point(aes(size=obsCount)) } else {
      plotOut <- plotOut + geom_point()
    }

  plotOut

}

