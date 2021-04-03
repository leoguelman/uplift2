mom <- function(x, ...) UseMethod("mom")

mom.default <- function(x, ...) stop("uplift: 'x' should be a formula object")

#'Modified outcome method for uplift modeling.
#'
#'\code{mom} transforms the response variable in a way that is relevant for
#'subsequent uplift modeling. It handles continuous (uncensored) and categorical
#'responses. A model fitted to this transformed response has a causal
#'interpretation for the treatment effect conditional on the covariates.
#'
#'Let \eqn{T \in [-1,1]} be a binary treatment indicator with \eqn{T=1} being
#'the treatment level of interest (i.e., the treatment group). Also, let \eqn{y}
#'be a response variable. If the response is a factor, the transformed response
#'is set to 1 if \eqn{T=1} and \eqn{y=1}, or if \eqn{T=-1} and \eqn{y=0}
#'(assuming the \code{classLevel} of interest for \eqn{y} is 1). Otherwise, the
#'transformed response is set to 0. Under the specific case in which
#'\eqn{Prob(T=1) = Prob(T=-1) = 1/2}, it is easy to show that
#'
#'\deqn{2 * Prob(z=1|X) - 1 = Prob(y=1|T = 1, X) - Prob(y=1|T = -1, X)}
#'(Jaskowski and Jaroszewicz, 2012), where y, z, and X denote the original
#'response variable, the transformed response, and the covariates, respectively.
#'
#'If the response is numeric, it is transformed as \eqn{z = 2 * (y - \bar{y}) *
#'T}. A model fitted to \eqn{z} effectively estimates \eqn{E[y|T = 1, X] - E[y|T
#'= -1, X]} (Tian et al., 2014).
#'
#'The argument \code{sampling} can be used to obtain a balanced treatment
#'distribution. Specifically, if \code{sampling = "oversample"}, observations
#'from the treatment minority class are duplicated (by sampling with
#'replacement), so that the resulting data frame has exactly the same number of
#'observations under each treatment level. Alternatively, if \code{sampling =
#'"undersample"}, observations from the treatment majority class are dropped (by
#'sampling without replacement), so that the resulting data frame has exactly
#'the same number of observations under each treatment level. If \code{sampling
#'= "none"}, no sampling is done. Lastly, if \code{sampling = "weights"}, the
#'returned data frame includes a weight variable that equals (1 - \eqn{\pi}) for
#'T = \code{treatLevel}  and \eqn{\pi} otheriwse, where \eqn{\pi = Prob(T =
#'treatLevel)}. The weight variable can be subsequently used to perform
#'case-weighted regression/classification on the transformed response.
#'
#'@name mom
#'
#'@aliases mom
#'
#'@export mom mom.default mom.formula
#'
#'@param formula A model formula of the form y ~ x1 + ....+ xn + trt(), where
#'  the left-hand side corresponds to the observed response, the right-hand side
#'  corresponds to the predictors, and 'trt' is the special expression to mark
#'  the treatment term. If the treatment term is not a factor, it is converted
#'  to one.
#'@param data A data frame in which to interpret the variables named in the
#'  formula.
#'@param subset Expression indicating which subset of the rows of data should be
#'  included. All observations are included by default.
#'@param na.action A missing-data filter function. Defaults to \code{na.omit}.
#'@param sampling The sampling method used to balance the treatment variable.
#'  See details.
#'@param newRespName The name for the transformed response variable.
#'@param classLevel A character string for the class of interest. Only
#'  applicable when the response is a factor. Defaults to the last level of the
#'  factor.
#'@param treatLevel A character string for the treatment level of interest.
#'  Defaults to the last level of the treatment factor.
#'@param x A \code{mom} object.
#'@param \dots Additional arguments for the S3 methods.
#'
#'
#'@return An object of class \code{"mom"}, which is a list with the following
#'  components (among others passed to the S3 methods): \itemize{\item
#'  \code{data} The data set including the original response variable, the
#'  treatment indicator, the transformed response, the predictors, and
#'  (optionally) a weight variable. \item\code{call} The original call to
#'  \code{mom}. }
#'
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@references
#'
#'Guelman, L., Guillen, M., and Perez-Marin A.M. (2015). "A decision support
#'framework to implement optimal personalized marketing interventions." Decision
#'Support Systems, Vol. 72, pp. 24--32.
#'
#'Jaskowski, M. and Jaroszewicz, S. (2012)."Uplift modeling for clinical trial
#'data". In ICML 2012 Workshop on Machine Learning for Clinical Data Analysis,
#'Edinburgh, Scotland.
#'
#'Tian, L., Alizadeh, A., Gentles, A. and Tibshirani, R. (2014). "A simple
#'method for detecting interactions between a treatment and a large number of
#'covariates." Journal of the American Statistical Association, 109:508, pp.
#'1517--1532,
#'
#'@examples
#'
#'set.seed(324)
#'df_train <- sim_uplift(p = 15, response = "binary")
#'df_train_mcm <- mom(y ~  X1 + X2 + X3 + trt(T),
#'                    data = df_train, sampling = "undersample")

mom.formula <- function(formula,
                        data,
                        subset,
                        na.action,
                        sampling =  "none",
                        newRespName = "z",
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
  attrTerms <- attributes(terms)
  yLabel <- names(attrTerms$dataClasses[1])
  XLabels <- attrTerms$term.labels[-(trtVar-1)]
  XClasses <- attrTerms$dataClasses[-c(1, trtVar)]
  X <- mf[XLabels]
  ttLabel <- attrTerms$term.labels[trtVar-1]
  ttLabel <- sub("\\).*", "", sub(".*\\(", "", ttLabel))
  .tt.f <- factor(.tt.)
  ttLevels <- levels(.tt.f)
  ttNLevels <- length(ttLevels)
  if (is.null(treatLevel)) treatLevel <- ttLevels[ttNLevels]
  .tt.mod <- ifelse(.tt.f == treatLevel, 1L, -1L)
  validSampling <- charmatch(sampling,
                             c("none", "undersample", "oversample", "weights"))
  if (any(is.na(validSampling)))
    stop("uplift: not recognized 'sampling' argument")
  s <- balance_sample(.tt.mod, sampling, 1L)
  weights <- s$weights
  inbag <- s$inbag
  yFactor <- is.factor(y)
  if (!yFactor && length(unique(y)) == 2)
    warning("uplift:", sQuote(yLabel), " has 2 unique values.
            Consider treating as factor instead of numeric.")
  if (yFactor) {
    if (is.null(classLevel)) classLevel <- levels(y)[nlevels(y)]
    z <- factor(ifelse((y == classLevel & .tt.f == treatLevel) |
                         (y != classLevel & .tt.f != treatLevel), 1, 0), levels = c("0", "1"))
  } else {
    if (!is.null(classLevel))
      warning("uplift: classLevel is ignored when response is not a factor.")
    classLevel = "NA"
    mean.y <- mean(y, na.rm = TRUE)
    z <- ifelse(.tt.f == treatLevel, 2 * (y - mean.y), -2 * (y - mean.y))
  }
  df.temp <- setNames(data.frame(y, .tt.f, z), c(yLabel, ttLabel, newRespName))
  if (sampling == "weights") {
    res <- data.frame(df.temp, X, weights = weights)
  } else {
    res <- data.frame(df.temp, X)[inbag, ]
  }

  lsOut <- list(data = res,
                yFactor = yFactor,
                classLevel = classLevel,
                treatLevel = treatLevel,
                newRespName = newRespName,
                XLabels = XLabels,
                XClasses =  XClasses,
                yLabel = yLabel,
                ttLabel = ttLabel,
                sampling = sampling,
                call = call)

  class(lsOut) <- "mom"
  lsOut
}

#' @rdname mom
#' @method print mom
#' @export print.mom

print.mom <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("mom", "\n", sep="")
  cat("Response type: ", ifelse(x$yFactor, "factor", "numeric"), "\n", sep="")
  cat("classLevel: ", x$classLevel, "\n", sep="")
  cat("treatLevel: ", x$treatLevel, "\n", sep="")
  cat("newRespName: ", x$newRespName, "\n", sep="")
  cat("sampling: ", x$sampling, "\n", sep="")
  invisible(x)
}
