#'Data splitting functions for uplift.
#'
#'\code{create_data_partition} creates one or more random data partitions into
#'training and test sets. \code{create_kfolds} splits the data into k-folds (or
#'groups) with approximatly the same number of observations.
#'\code{create_bootstrap} draws bootstrap samples.
#'
#'If \code{y} is a factor, sampling is done within the levels of \code{y} in an
#'attempt to balance the class distributions between the partitions. If \code{y}
#'is numeric, groups are first created based on the quantiles of its
#'distribution and then sampling is done within these groups.
#'
#'If \code{trt} is supplied, the data partitions are stratified by the treatment
#'variable.
#'
#'Notice that in addition to \code{create_bootstrap}, bootstrap samples can also
#'by created using \code{create_data_partition} with \code{p = 1} and
#'\code{replace = TRUE}.
#'
#'
#'@name create_data_partition
#'
#'@aliases create_data_partition
#'
#'@export create_data_partition
#'
#'@param y An atomic vector.
#'@param trt An optional treatment variable.
#'@param p The proportion of training observations.
#'@param times The number of partitions to create.
#'@param groups For numeric y, the number of breaks in the quantiles.
#'@param replace Should sampling be done with replacement?
#'@param k The number of folds.
#'
#'@return \code{create_data_partition} and \code{create_bootstrap} return a
#'  matrix of row position integers corresponding to the training set and to the
#'  bootstrap sample, respectively. \code{create_kfolds} returns a matrix with
#'  the row integers corresponding to the folds.
#'
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@examples
#'
#'set.seed(545)
#'r <- factor(sample(c(0,1), 1000, replace = TRUE))
#'t <- factor(sample(c(0,1), 1000, replace = TRUE))
#'df <- data.frame(r, t)
#'trainIndex <- create_data_partition(df$r, df$t)
#'dfTrain <- df[trainIndex, ]
#'dfTest <- df[-trainIndex, ]
#'table(df$r, df$t)
#'table(dfTrain$r, dfTrain$t)
#'table(dfTest$r, dfTest$t)
#'# Create k-folds
#'head(create_kfolds(r, t, times = 5))
#'
#'# Create 10 bootstrap samples
#'set.seed(1)
#'x <- rnorm(100)
#'xb <- create_bootstrap(x, times = 10)


create_data_partition <- function(y,
                                  trt = NULL,
                                  p = 0.5,
                                  times = 1,
                                  groups = 5,
                                  replace = FALSE) {

  if (!is.atomic(y)) stop("uplift: y must be an atomic vector.")
  if (p <= 0 | p > 1) stop("uplift: p must be in the (0,1] interval.")
  if (times < 1) stop("uplift: times must be equal or greater than 1.")
  if (!is.null(trt)) {
    trtName <- deparse(substitute(trt))
    trt <- as.factor(trt)
  } else {
    if (!missing(trt))
      stop("uplift: ", deparse(substitute(trt)), " does not exist.")
  }
  yName <- deparse(substitute(y))
  if (is.numeric(y) && length(unique(y)) < 5)
    warning("uplift:", sQuote(yName), " has less than 5 unique values.
            Consider treating as factor instead of numeric.")
  n <- length(y)
  if (is.factor(y)) {
    if (is.null(trt)) {
      ySplit <- split(1:n, y)
    } else {
      ySplit <- split(1:n, list(y, trt))
    }
  } else if (is.numeric(y)) {
    if (groups < 2) stop("uplift: Number of groups must be at least 2.")
    yQuants <- cut(y, breaks = unique(quantile(y, probs = seq(0, 1, 1 / groups))),
                   include.lowest = TRUE, labels = FALSE)
    if (is.null(trt)) {
      ySplit <- split(1:n, yQuants)
    } else {
      ySplit <- split(1:n, list(yQuants, trt))
    }
  }

  nPerSplit <- lapply(ySplit, length)
  nClasses <- length(ySplit)
  mPicks <- lapply(nPerSplit, function(i) round(i * p))
  m0Picks <- vapply(mPicks, function(x) x == 0, logical(1))
  if (is.factor(y) && any(m0Picks)) {
    warning("uplift: Some levels of ", sQuote(yName)," and/or ", sQuote(trtName),
            " are not represented in the training samples.")
  } else if (is.numeric(y) && any(m0Picks)) {
    warning("uplift: Some levels of ", sQuote(trtName),
            " are not represented in the training samples.")
  }
  indexMatrix <- matrix(nrow = sum(unlist(mPicks)), ncol = times)
  for (i in 1:times) {
    indexMatrix[, i] <- unlist(lapply(1:nClasses, function(k) sample(ySplit[[k]],
                                                                     mPicks[[k]], replace = replace)))
  }
  colnames(indexMatrix) <- paste("trainIndex", 1:ncol(indexMatrix), sep = '')
  indexMatrix
}


#' @rdname create_data_partition
#' @export create_kfolds

create_kfolds <- function(y,
                          trt = NULL,
                          k = 10,
                          times = 1,
                          groups = 5) {

  if (!is.atomic(y)) stop("uplift: y must be an atomic vector.")
  if (k < 1) stop("uplift: k must be equal or greater than 1.")
  if (times < 1) stop("uplift: times must be equal or greater than 1.")
  if (!is.null(trt)) {
    trtName <- deparse(substitute(trt))
    trt <- as.factor(trt)
  } else {
    if (!missing(trt))
      stop("uplift: ", deparse(substitute(trt)), " does not exist.")
  }
  yName <- deparse(substitute(y))
  if (is.numeric(y) && length(unique(y)) < 5)
    warning("uplift:", sQuote(yName), " has less than 5 unique values.
            Consider treating as factor instead of numeric.")

  n <- length(y)
  if (is.factor(y)) {
    if (is.null(trt)) {
      ySplit <- split(1:n, y)
    } else {
      ySplit <- split(1:n, list(y, trt))
    }
  } else if (is.numeric(y)) {
    if (groups < 2) stop("uplift: Number of groups must be at least 2.")
    yQuants <- cut(y, breaks = unique(quantile(y, probs = seq(0, 1, 1 / groups))),
                   include.lowest = TRUE, labels = FALSE)
    if (is.null(trt)) {
      ySplit <- split(1:n, yQuants)
    } else {
      ySplit <- split(1:n, list(yQuants, trt))
    }
  }

  nPerSplit <- lapply(ySplit, length)
  nClasses <- length(ySplit)
  foldMatrix <- matrix(nrow = n, ncol = times)
  for (i in 1:times) {
    foldsInSplit <- lapply(1:nClasses, function(j) if (nPerSplit[[j]] > 0)
      cut(sample(1:nPerSplit[[j]], nPerSplit[[j]], replace = FALSE),
          breaks = k,
          include.lowest = TRUE,
          labels = FALSE))
    folds <- unlist(foldsInSplit)[order(unlist(ySplit))]

    if (is.null(trt)) {
      na <- any(table(y, folds) == 0)
    } else {
      na <- any(table(y, folds, trt) == 0)
    }
    if (na)
      warning("uplift: Some levels of ", sQuote(yName)," and/or ", sQuote(trtName),
              " are not represented in at least one of the ", sQuote(k)," folds,
              for at least one of the ", sQuote(times), " created partitions.")
    foldMatrix[, i] <- folds
  }
  colnames(foldMatrix) <- paste("Resample", 1:ncol(foldMatrix), sep = '')
  foldMatrix
}


#' @rdname create_data_partition
#' @export create_bootstrap

create_bootstrap <- function(y, times = 10) {

  if (!is.atomic(y)) stop("uplift: y must be an atomic vector.")
  if (times < 1) stop("uplift: times must be equal or greater than 1.")

  n <- length(y)
  bootMatrix <- vapply(1:times, function(i) sample(1:n, n, replace = TRUE), numeric(n))
  colnames(bootMatrix) <- paste("bootstrap", 1:ncol(bootMatrix), sep = '')
  bootMatrix
}


