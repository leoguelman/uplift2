uplift_glm <- function(x, ...) UseMethod("uplift_glm")

uplift_glm.default <- function(x, ...) stop("uplift: 'x' should be a formula object")

#'Fitting Uplift Generalized Linear Models.
#'
#'\code{uplift_glm} fits Uplift Generalized Linear Models, optionally with lasso
#'or elasticnet regularization.
#'
#'The function follows the method for uplift modeling proposed by Tian et al.
#'(2014). This method consists in modifying the covariates in a simple way, and
#'then fitting an appropriate regression model using the modified covariates and
#'no main effects. See Tian et al. (2014) for details.
#'
#'The argument \code{sampling} can be used to obtain a balanced treatment
#'distribution. Specifically, if \code{sampling = "oversample"}, observations
#'from the treatment minority class are duplicated (by sampling with
#'replacement), so that the data frame used in model fitting has exactly the
#'same number of observations under each treatment level. Alternatively, if
#'\code{sampling = "undersample"}, observations from the treatment majority
#'class are dropped (by sampling without replacement), so that the data frame
#'used in model fitting has exactly the same number of observations under each
#'treatment level. If \code{sampling = "none"}, no sampling is done. Lastly, if
#'\code{sampling = "weights"}, the returned data frame includes a weight
#'variable that equals (1 - \eqn{\pi}) for T = \code{treatLevel}  and \eqn{\pi}
#'otherwise, where \eqn{\pi = Prob(T = treatLevel)}. These weights are
#'subsequently used as case weights in the fitting process.
#'
#'@name uplift_glm
#'
#'@aliases uplift_glm
#'
#'@export uplift_glm uplift_glm.default uplift_glm.formula
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
#'@param na.action A missing-data filter function.
#'@param family Response type. For \code{family = "gaussian"} (default), the
#'  response must be presented as numeric. For \code{family = "binomial"}, the
#'  response must be a factor with two levels. If the response is numeric,  it
#'  will be coerced into a factor. For \code{family = "cox"}, the response must
#'  be a survival object, as returned by \code{survival::Surv}.
#'@param method The method used for model fitting. If \code{method = "glm"}
#'  (default), the model is fitted on the modified covariates (see details)
#'  using \code{stats::glm}. If \code{method = "glmStepAIC"}, the model is first
#'  fitted using \code{stats::glm} and then this is passed to
#'  \code{MASS::stepAIC} for AIC stepwise selection. Alternatively, for
#'  \code{method = "glmnet"} and \code{method = "cv.glmnet"}, models are fitted
#'  using \code{glmnet::glmnet} and \code{glmnet::cv.glmnet}, respectively.
#'@param sampling The sampling method used to balance the treatment variable.
#'  See details.
#'@param treatLevel A character string for the treatment level of interest.
#'  Defaults to the last level of the treatment factor.
#'@param Anova If \code{TRUE}, the analysis-of-variance table is returned using
#'  the function \code{car::Anova}. It does not apply to \code{method =
#'  "cv.glmnet"} or \code{method = "glmnet"}.
#'@param \dots Additional arguments passed to the regression method selected in
#'  \code{method}.
#'@param x A \code{uplift_glm} object.
#'
#'
#'@return An object of class \code{"uplift_glm"}, which is a list with the
#'  following components, in addition to the ones returned by the specific
#'  fitting method:
#'
#'  \itemize{ \item \code{call} The calling expression \item \code{na.action}
#'  Information returned by \code{model.frame} on the special handling of NAs.
#'  \item \code{xlevels} The levels of predictors of class factor. \item
#'  \code{Family} The \code{family} used. \item \code{method} The \code{method}
#'  used. \item \code{sampling} The \code{sampling} method used. \item
#'  \code{dataClasses} The data classes of predictors. \item \code{treatLevel}
#'  The reference treatment level. \item \code{ttReLabel} The label of the
#'  transformed treatment indicator. \item \code{modForm} The model formula.
#'  \item \code{modData} The data frame used in model fitting. \item
#'  \code{inbag} The index of of which observations were used for fitting. \item
#'  \code{weightVector} The vector of weights used for fitting.}
#'
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
#'@examples
#'
#'set.seed(1)
#'df_train <- sim_uplift(p = 50, response = "binary")
#'df_test<- sim_uplift(p = 50, n = 10000, response = "binary")
#'form <- as.formula(paste('y ~', 'trt(T) +',
#'        paste('X', 1:(ncol(df_train)-3), sep = '', collapse = "+")))
#'fit1 <- uplift_glm(form,
#'                   family = "binomial",
#'                   method = "glm",
#'                   data = df_train)
#'fit1
#'fit2 <- uplift_glm(form,
#'                   family = "binomial",
#'                   method = "glmStepAIC",
#'                   data = df_train)
#'fit2
#'fit3 <- uplift_glm(form,
#'                   family = "binomial",
#'                   method = "cv.glmnet",
#'                   data = df_train)
#'lambda.opt <- fit3$lambda.min
#'fit3 <- uplift_glm(form,
#'                   family = "binomial",
#'                   method = "glmnet",
#'                   data = df_train)
#'upliftPred1 <- predict(fit1, df_test)
#'upliftPred2 <- predict(fit2, df_test)
#'upliftPred3 <- predict(fit3, df_test, s=lambda.opt)
#'df_eval<- data.frame(upliftPred1 = upliftPred1,
#'                     upliftPred2 = upliftPred2,
#'                     upliftPred3 = upliftPred3,
#'                     y = df_test$y,
#'                     T = df_test$T)
#'res <- inspect_performance(y ~ upliftPred1 + upliftPred2 + upliftPred3 + trt(T),
#'                           data = df_eval, qini = TRUE)
#'res
#'summary(res)
#'ggplot(res)
#'res$qiniC

uplift_glm.formula <- function(formula,
                               data,
                               subset,
                               na.action,
                               family = "gaussian",
                               method = "glm",
                               sampling =  "weights",
                               treatLevel = NULL,
                               Anova = FALSE,
                               ...)
{
  call <- match.call()
  if (!inherits(formula, "formula"))
    stop("uplift: method is only for formula objects.")
  if(!(sampling %in% c("none", "undersample", "oversample", "weights")))
    stop("uplift: sampling must be either \"none\", \"undersample\", \"oversample\", \"weights\".")
  if(!(family %in% c("gaussian", "binomial", "cox")))
    stop("uplift: family must be either \"gaussian\", \"binomial\", \"cox\".")
  if(!(method %in% c("glm", "glmStepAIC", "cv.glmnet", "glmnet")))
    stop("uplift: method must be either \"glm\", \"glmStepAIC\", \"cv.glmnet\", \"glmnet\".")

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
    .tt. <- mf[, trtVar] else
      stop("uplift: formula must include a treatment variable.")
  .tt.f <- factor(.tt.)
  ttLevels <- levels(.tt.f)
  ttNLevels <- length(ttLevels)
  if (is.null(treatLevel)) treatLevel <- ttLevels[ttNLevels]
  .tt.mod <- ifelse(.tt.f == treatLevel, .5, -.5)
  s <- balance_sample(.tt.mod, sampling, .5)
  weights <- s$weights
  inbag <- s$inbag
  attrTerms <- attributes(terms)
  ttLabel <- attrTerms$term.labels[trtVar-1]
  ttLabel <- sub("\\).*", "", sub(".*\\(", "", ttLabel))
  ttReLabel <- paste("T_", ttLabel, sep = "")
  yLabel <- names(attrTerms$dataClasses[1])
  XLabels <- attrTerms$term.labels[-(trtVar-1)]
  dataClasses <- attrTerms$dataClasses[-c(1, trtVar)]
  Terms <- terms(reformulate(XLabels))
  modData <- data.frame(mf[inbag, -trtVar, drop = FALSE], weights = weights,
                        .Int = 1L, .tt.mod[inbag])
  modData <- setNames(modData, c(colnames(modData)[-ncol(modData)], ttReLabel))
  modForm <- as.formula(paste(yLabel, '~  -1 +',
                              paste(c(".Int", XLabels), ":", ttReLabel, sep = ' ', collapse = "+")))

  if (method == "glm" && family == "cox") {
    method <- "coxph"
  } else if (method == "glmStepAIC") {
    # intercept should not be eliminated from the stepAIC procedure
    lScope <- paste("~.Int:", ttReLabel, sep = '')
  }

  res <- switch(method,
                glm = withCallingHandlers(glm(modForm,
                                              family = family,
                                              data = modData,
                                              weights = modData$weights,
                                              na.action = na.action,
                                              ...), warning = remove_glmbinom_warn),

                glmStepAIC = withCallingHandlers(MASS::stepAIC(withCallingHandlers(glm(modForm,
                                                                                       family = family,
                                                                                       data = modData,
                                                                                       weights = modData$weights,
                                                                                       na.action = na.action,
                                                                                       ...), warning = remove_glmbinom_warn),
                                                               trace = FALSE, scope = list(lower = eval(parse(text = lScope))),
                                                               ...), warning = remove_glmbinom_warn),

                glmnet = glm_net(modForm,
                                 family = family,
                                 data = modData,
                                 which = "glmnet",
                                 weights = modData$weights,
                                 intercept = FALSE,
                                 na.action = na.action,
                                 ...),

                cv.glmnet = glm_net(modForm,
                                    family = family,
                                    data = modData,
                                    which = "cv.glmnet",
                                    weights = modData$weights,
                                    intercept = FALSE,
                                    na.action = na.action,
                                    ...),

                coxph = survival::coxph(modForm,
                                        data = modData,
                                        weights = modData$weights,
                                        na.action = na.action,
                                        ...)
  )
  if (Anova) {
    if (method %in% c("glm", "glmStepAIC", "coxph")) {
      res$Anova <- withCallingHandlers(car::Anova(res), warning = remove_glmbinom_warn)
    } else {
      res$Anova <- NULL
      warning("uplift: \"Anova\" is not applicable to method ", method)
    }
  }

  res$call <- call
  res$na.action <- attr(mf, "na.action")
  #res$Terms <- Terms
  #res$contrasts <- attr(X, "contrast")
  res$xlevels <- .getXlevels(Terms, mf)
  res$Family <- family
  res$method <- method
  res$sampling <- sampling
  res$dataClasses <- dataClasses
  res$treatLevel <- treatLevel
  res$ttReLabel <- ttReLabel
  res$modForm <- modForm
  res$modData <- modData
  res$inbag <- inbag
  res$weightVector <- weights
  class(res) <- c("uplift_glm", class(res))
  res
}

#######################################
# Print method
#######################################

#'@rdname uplift_glm
#'@method print uplift_glm
#'@export print.uplift_glm

print.uplift_glm <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("uplift_glm", "\n", sep="")
  cat("\n")
  cat("family: ", x$Family, "\n", sep="")
  cat("method: ", x$method, "\n", sep="")
  cat("sampling: ", x$sampling, "\n", sep="")
  cat("Number of observations used: ", nrow(x$modData), "\n", sep="")
  cat("treatLevel: ", x$treatLevel, "\n", sep="")
  invisible(x)
}

#######################################
# Formula method for glmnet
#######################################


glm_net <- function(formula, data, subset, na.action, which = "glmnet", ...) {
  call <- match.call()
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("formula", "data", "subset", "na.action"),
                names(mf), 0)
  mf <- mf[c(1, args)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval.parent(mf)
  terms <- attr(mf, "terms")
  y <- model.response(mf)
  X <- model.matrix(terms, mf, contrasts)
  intercept <- which(colnames(X) == "(Intercept)")
  if (length(intercept > 0)) X <- X[, -intercept]
  switch(which,
         glmnet = withCallingHandlers(glmnet::glmnet(X, y, ...), warning = remove_coxglmnet_warn),
         cv.glmnet = withCallingHandlers(glmnet::cv.glmnet(X, y, ...), warning = remove_coxglmnet_warn)
  )
}

#######################################
# Remove specific warnings
#######################################

# removes warning produced by glmnet when fitting a cox model with option interept = FALSE
remove_coxglmnet_warn <- function(x) {
  if(any(grepl("Cox model has no intercept", x)))
    invokeRestart("muffleWarning")
}


# Need a better way to handle the warning message below, essentially due to the fact that
# model.matrix may create different design matrices simply by swapping the interaction
# terms in the model formula. See example below.

#set.seed(1)
#x1 <- rnorm(100)
#f1 <- factor(sample(letters[1:3], 100, replace = TRUE))
#trt <- sample(c(-1,1), 100, replace = TRUE)
#y <- factor(sample(c(0,1), 100, T))
#df <- data.frame(y=y, x1=x1, f1=f1, trt=trt)
#fit1 <- glm(y ~ x1:trt + f1:trt, data = df, family = binomial)
#coef(fit1)
#fit2 <- glm(y ~ f1:trt + x1:trt, data = df, family = binomial)
#coef(fit2)
#all.equal(fitted(fit1), fitted(fit2))

# removes warning in rank-deficient fit due to not including main effects but only interaction effects in the model formula in glm
remove_glmpred_warn <- function(x) {
  if(any(grepl("prediction from a rank-deficient fit may be misleading", x)))
    invokeRestart("muffleWarning")
}
# remove warning when fitting binomial glm with weights

remove_glmbinom_warn <- function(x) {
  if(any(grepl("non-integer #successes in a binomial glm!", x)))
    invokeRestart("muffleWarning")
}


#######################################
# Prediction for "coxnet" objects
#######################################

predict_coxnet <- function(object,
                           newx,
                           s = NULL,
                           type = c("link", "response"),
                           offset,
                           ...) {
  type <- match.arg(type)
  nbeta <- object$beta
  if (!is.null(s)) {
    vnames <- dimnames(nbeta)[[1]]
    dimnames(nbeta) <- list(NULL, NULL)
    lambda <- object$lambda
    lamlist <- lambda.interp(lambda,s)
    nbeta <- nbeta[,lamlist$left, drop = FALSE] * lamlist$frac + nbeta[, lamlist$right, drop = FALSE] * (1-lamlist$frac)
    dimnames(nbeta) <- list(vnames, paste0(seq_along(s)))
  }
  nfit <- as.matrix(newx %*% nbeta)
  if (object$offset) {
    if (missing(offset))
      stop("uplift: No offset provided for prediction, yet used in fit of glmnet", call.=FALSE)
    nfit <- nfit + array(offset, dim = dim(nfit))
  }
  switch(type,
         response = exp(nfit),
         link = nfit
  )
}

#######################################
# predict method
#######################################

#'Predict method for uplift_glm fits.
#'
#'Obtains predictions from a fitted \code{uplift_glm} object.
#'
#'@param object A fitted object inheriting from \code{uplift_glm}.
#'@param newdata An optional set of data to predict on. If \code{NULL}, then the
#'  original data are used.
#'@param type The type of predictions required. Only \code{"uplift"} predictions
#'  are allowed.
#'@param na.action The method for handling missing data.
#'@param \dots Additional arguments to be passed to other methods.
#'
#'@return A numeric vector of predictions for all methods, except when the model
#'  was fitted using \code{method = "glmnet"}, in which case the returned object
#'  is a matrix of predictions for the entire sequence of the penalty parameter
#'  used to create the model.
#'
#'@export predict.uplift_glm
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@seealso \code{\link{uplift_glm}}


predict.uplift_glm <- function(object, newdata = NULL,  type = "uplift", na.action = na.omit, ...) {

  if (!inherits(object, "uplift_glm"))
    stop("uplift: object must be of class 'uplift_glm'.")
  if(type != "uplift")
    stop("uplift: only \"uplift\" type is allowed.")
  if(!is.null(newdata)) {
    newdata <- as.data.frame(newdata)
  } else newdata <- object$modData
  .tt.Ind <- which(colnames(newdata) == object$ttReLabel)
  if (length(.tt.Ind > 0)) newdata <- newdata[, -.tt.Ind]
  newdata$.tt. <- 1 # predict using non-transformed predictors
  newdata$.Int <- 0 # exclude intercept from predictions
  ttLabelInd <- which(colnames(newdata) == ".tt.")
  colnames(newdata)[ttLabelInd] <- object$ttReLabel
  if (object$method %in% c("glmnet", "cv.glmnet")) {
    Terms <- delete.response(terms(object$modForm))
    mf <- model.frame(Terms, newdata, na.action = na.action, xlev = object$xlevels)
    if (!is.null(object$dataClasses)) .checkMFClasses(object$dataClasses, mf)
    newdata <- model.matrix(Terms, mf, contrasts)
  }
  if (object$method %in% c("glm", "glmStepAIC")) {
    lpred <- withCallingHandlers(stats::predict.glm(object = object, newdata = newdata, type = "link",
                                                    na.action = na.action, ...),
                                 warning = remove_glmpred_warn)
  } else if (object$method == "coxph") {
    lpred <- survival:::predict.coxph(object = object, newdata = newdata, type = "lp",
                                      na.action = na.action, ...)
  } else if (object$method == "glmnet" && object$Family != "cox") {
    lpred <- glmnet::predict.glmnet(object = object, newx = newdata, type = "link", ...)
    if (dim(lpred)[2] == 1L) lpred <- as.vector(lpred)
  } else if (object$method == "cv.glmnet" && object$Family != "cox") {
    lpred <- as.vector(glmnet::predict.cv.glmnet(object = object, newx = newdata, type = "link", ...))
  } else if (object$method %in% c("glmnet", "cv.glmnet") && object$Family == "cox") {
    lpred <- predict_coxnet(object = object, newx = newdata, type = "link", ...)
  }
  # Follow Tian et al. (JASA, 2014, 109:508, 1517-1532)
  if (object$Family == "binomial") {
    resp <- (exp(lpred/2) - 1) / (exp(lpred/2) + 1)
  } else if (object$Family == "cox") {
    resp <- exp(-lpred)
  } else resp <- lpred
  resp
}


