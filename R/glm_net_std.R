glm_net_std <- function(formula, data, subset, na.action, which = "glmnet", penCts = TRUE, stdx = TRUE, ...) {
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
  if (penCts) {
    ocontrasts <- options(contrasts = c(unordered = "contr.none", ordered = "contr.diff"))
    on.exit(options(ocontrasts))
  }
  X <- model.matrix(terms, mf, contrasts)
  browser()
  if (stdx) X <- std_unit_l2norm(X)
  intercept <- which(colnames(X) == "(Intercept)")
  if (length(intercept > 0)) X <- X[, -intercept]
  switch(which,
         glmnet = glmnet::glmnet(X, y, ...),
         cv.glmnet = glmnet::cv.glmnet(X, y, ...))
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
