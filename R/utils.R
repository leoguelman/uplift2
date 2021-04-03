### Make misc. functions available
`%>%` <- magrittr::`%>%`
`%dopar%` <- foreach::`%dopar%`
`%do%` <- foreach::`%do%`

### Mark treatment term

#'Mark treatment term.
#'
#'An identity function, simply used to mark the treatment term in the model formula.
#'
#'@name trt
#'
#'@aliases trt
#'
#'@export trt
#'
#'@param x The treatment indicator variable.
#'
#'@return Same as \code{x}.
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}

trt <- function(x) x

### Create bins


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


### Balance treatment sampling

balance_sample <- function(T, sampling, treatLevelNum) {

  n <- length(T)
  ttSplit <- split(1:n, list(T))
  nPerSplit <- unlist(lapply(ttSplit, length))
  max.tt.l <- which.max(nPerSplit)
  min.tt.l <- which.min(nPerSplit)
  maxPicks <- nPerSplit[max.tt.l]
  minPicks <- nPerSplit[min.tt.l]

  if (sampling == "undersample") {
    inbag <- c(ttSplit[[min.tt.l]], sample(ttSplit[[max.tt.l]],
                                           minPicks, replace = FALSE))
    weights <- rep(1, length(inbag))
  } else if (sampling == "oversample") {
    inbag <- c(ttSplit[[max.tt.l]], sample(ttSplit[[min.tt.l]],
                                           maxPicks, replace = TRUE))
    weights <- rep(1, length(inbag))
  } else if (sampling == "weights") {
    inbag <-  1L:n
    pPerSplit <- nPerSplit / n
    weights <- ifelse(T == treatLevelNum, pPerSplit[1], pPerSplit[2])
  } else {
    inbag <- 1L:n
    weights <- rep(1, n)
  }
  list(inbag = inbag, weights = weights)
}



#######################################
# unit L2-norm standardization
# follows lars package
#######################################

std_unit_l2norm <- function(X, eps = .Machine$double.eps) {
  np <- dim(X)
  n <- np[1]
  p <- np[2]
  one <- rep(1, n)
  normX <- sqrt(drop(one %*% (X^2)))
  nosignal <- normX/sqrt(n) < eps
  if(any(nosignal)) normX[nosignal] <- eps*sqrt(n)
  names(normX) <- NULL
  X <- scale(X, FALSE, normX)
}

#######################################
# An alternative contrasts function for unordered factors
# Ensures symmetric treatment of all levels of a factor
# (penalized package)
#######################################

contr.none <- function(n, contrasts) {
  if (length(n) == 1)
    contr.treatment(n, contrasts = n<=2)
  else
    contr.treatment(n, contrasts = length(unique(n))<=2)
}

#######################################
# An alternative contrasts function for ordered factors
# Ensures use of a difference penalty for such factors
# (penalized package)
#######################################

contr.diff <- function (n, contrasts = TRUE)
{
  if (is.numeric(n) && length(n) == 1) {
    if (n > 1)
      levs <- 1:n
    else stop("not enough degrees of freedom to define contrasts")
  }
  else {
    levs <- n
    n <- length(n)
  }
  contr <- array(0, c(n, n), list(levs, paste(">=", levs, sep="")))
  contr[outer(1:n,1:n, ">=")] <- 1
  if (n < 2)
    stop(gettextf("contrasts not defined for %d degrees of freedom",
                  n - 1), domain = NA)
  if (contrasts)
    contr <- contr[, -1, drop = FALSE]
  contr
}
