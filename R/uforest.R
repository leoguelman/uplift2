uforest <- function(x, ...) UseMethod("uforest")

uforest.default <- function(x, ...) stop("uplift: 'x' should be a formula object")

#'Control for uplift random forest.
#'
#'Various parameters that control aspects of the \code{uforest} fit.
#'
#'@name uforest_control
#'
#'@aliases uforest_control
#'
#'@export uforest_control
#'
#'@param ntree Number of uplift trees to fit.
#'@param minsplit The minimum number of observations in a node in order to be
#'  considered for splitting.
#'@param minbucket.t The minimum number of treatment observations in any
#'  terminal <leaf> node. The \code{treatLevel} can be used to determine the
#'  treatment level of interest.
#'@param minbucket.c The minimum number of control observations in any terminal
#'  <leaf> node.
#'@param var.select.criterion The criterion used to select the variable for
#'  splitting. At the moment, only \code{"pvalue"} is accepted. The variable
#'  with minimum pvalue is selected for splitting.
#'@param var.select.test The conditional null distribution of the test
#'  statistic. This is passed to the \code{distribution} argument in
#'  \code{coin::independence_test}. For example, for an approximative (Monte
#'  Carlo) reference distribution with B Monte Carlo replicates, use
#'  \code{approximate(B=999)}.
#'@param alpha The maximum acceptable pvalue required in order to make a split.
#'@param bonferroni Apply bonferroni adjustment to pvalue?
#'@param balance.sample The sampling method used to balance the treatment
#'  variable. This attempts to have an equal representation of each treatment
#'  before implementing the independence test described in
#'  \code{var.select.test}. The options are \code{"undersample"} (default),
#'  \code{"oversample"}, \code{"none"}. See the argument \code{sampling} in
#'  \code{\link{mom}} for details.
#'@param split.criterion The split criteria used at each node of each tree;
#'  Possible values are: \code{"uplift"} (default), \code{"kld"}
#'  (Kullback-Leibler divergence), \code{"ed"} (Euclidean divergence),
#'  \code{"l1d"} (L1-norm divergence). See details in Guelman et al. (2015).
#'@param maxdepth Maximum depth of the tree. The default \code{maxdepth = Inf}
#'  means that no restrictions are applied to tree sizes.
#'@param mtry Number of input variables randomly sampled as candidates at each
#'  node. The default is \eqn{\sqrt{p}}, where \eqn{p} represents the number of
#'  covariates.
#'@param parallel If \code{TRUE}, computations are performed in parallel,
#'  otherwise they are done sequentially.
#'@param nCore The number of cores used. Default is: number of available
#'  cores-1.
#'
#'@return A list.
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'

uforest_control <- function(ntree = 100, minsplit = 40L, minbucket.t = 20L, minbucket.c = 20L,
                            var.select.criterion = "pvalue", var.select.test = "asymptotic()",
                            alpha = 0.05, bonferroni = FALSE, balance.sample = "undersample",
                            split.criterion = "uplift", maxdepth = Inf, mtry = NULL, parallel = TRUE,
                            nCore = NULL)
{
  list(ntree = ntree, minsplit = minsplit, minbucket.t = minbucket.t, minbucket.c = minbucket.c,
       var.select.criterion = var.select.criterion, var.select.test = var.select.test,
       alpha = alpha, bonferroni = bonferroni,  balance.sample = balance.sample,
       split.criterion = split.criterion, maxdepth = maxdepth, mtry = mtry, parallel = parallel,
       nCore = nCore)
}

#'Fitting uplift random forest.
#'
#'\code{uforest} implements uplift random forests.
#'
#'\code{uforest} builds a sequence of de-correlated uplift trees (see
#'\code{\link{utree}}) fitted on bootstrap samples of the training data.
#'Additionally, the best split at each node is selected among a subset of
#'predictors randomly selected at that node. See Guelman et al. (2015) for
#'details.
#'
#'@name uforest
#'
#'@aliases uforest
#'
#'@export uforest uforest.default uforest.formula
#'
#'@param formula A model formula of the form y ~ x1 + ....+ xn + trt(), where
#'  the left-hand side corresponds to the observed response, the right-hand side
#'  corresponds to the predictors, and 'trt' is the special expression to mark
#'  the treatment term. At the moment, \code{uforest} only handles binary
#'  responses.
#'@param data A data frame in which to interpret the variables named in the
#'  formula.
#'@param na.action A missing-data filter function.
#'@param classLevel A character string for the class of interest. Defaults to
#'  the last level of the factor.
#'@param treatLevel A character string for the treatment level of interest.
#'  Defaults to the last level of the treatment factor.
#'@param control A list with control parameters, see \code{\link{uforest_control}}.
#'@param \dots Arguments passed to \code{\link{uforest_control}}.
#'@param x An object of class \code{"uforest"}
#'
#'@return An object of class \code{"uforest"}.
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@references
#'
#'Guelman, L., Guillen, M., and Perez-Marin A.M. (2015). "A decision support
#'framework to implement optimal personalized marketing interventions." Decision
#'Support Systems, Vol. 72, pp. 24--32.
#'
#'Hothorn, T., Hornik, K. and Zeileis, A. (2006). "Unbiased recursive
#'partitioning: A conditional inference framework". Journal of Computational and
#'Graphical Statistics, 15(3): 651--674.
#'
#'Rzepakowski, Piotr and Jaroszewicz, Szymon. (2011). "Decision trees for uplift
#'modeling with single and multiple treatments". Knowledge and Information
#'Systems, 32(2) 303--327.
#'
#'Strasser, H. and Weber, C. (1999). "On the asymptotic theory of permutation
#'statistics". Mathematical Methods of Statistics, 8: 220--250.
#'
#'Su, X., Tsai, C.-L., Wang, H., Nickerson, D. M. and Li, B. (2009). "Subgroup
#'Analysis via Recursive Partitioning". Journal of Machine Learning Research 10,
#'141--158.
#'
#'@examples
#'
#'set.seed(1)
#'df <- sim_uplift(n = 1000, p = 50, response = "binary")
#'form <- create_uplift_formula(x = names(df)[-c(1:3)], y = "y", trt = "T")
#'fit <- uforest(form, data = df, maxdepth = 3, ntree = 10, nCore = 2)
#'fit
#'t1 <- fit$forest[[1]] # see structure of first tree
#'plot(t1,  main = "first tree...", gp = grid::gpar(cex = 0.5))

uforest.formula <- function(formula,
                            data,
                            na.action,
                            classLevel = NULL,
                            treatLevel = NULL,
                            control = uforest_control(...),
                            ...
)
{

  call <- match.call()

  stopifnot(control$minsplit > 2L, control$minbucket.c > 1L, control$minbucket.t > 1L)

  if (control$ntree < 1)
    stop("uplift: ntree must be >=1")

  if (!(control$var.select.criterion %in% "pvalue"))
    stop("uplift: currently only var.select.criterion \"pvalue\" is accepted.")

  if (!(is.character(control$var.select.test)))
    stop("uplift: var.select.test must be a character string.")

  if (!(control$split.criterion %in% c("uplift", "kld", "ed", "l1d")))
    stop("uplift: split.criterion must be one of \"uplift\", \"kld\", \"ed\", or \"l1d\".")

  if (control$alpha < 0 | control$alpha > 1)
    stop("uplift: alpha must be between 0 and 1")

  if (!(is.logical(control$bonferroni)))
    stop("uplift: bonferroni must be a logical string.")

  if (!(control$balance.sample %in% c("none", "undersample", "oversample")))
    stop("uplift: sampling must be either \"none\", \"undersample\", \"oversample\".")

  if (!(is.logical(control$parallel)))
    stop("uplift: parallel must be logical.")

  momObj <- mom(formula = formula, data = data, na.action = na.action, sampling = "none",
                classLevel = classLevel, treatLevel = treatLevel)

  df <- momObj$data
  yn <- momObj$newRespName
  yon <- momObj$yLabel
  xn <- momObj$XLabels
  tn <- momObj$ttLabel
  yc <- momObj$yFactor
  classLevel <- momObj$classLevel
  treatLevel <- momObj$treatLevel
  XClasses <- momObj$XClasses
  weights <- rep(1L, nrow(df))
  if (!yc)
    stop("uplift: uforest can only handle binary response variables of class factor.")

  if (is.null(control$mtry)) {
    mtry <- floor(sqrt(length(xn)))
  } else mtry <- control$mtry

  ### pre-sorting numeric predictors
  df.s <- lapply(1:length(xn), function(i) sort_fun(df, i))
  names(df.s) <- xn
  df.s <- Filter(Negate(is.null), df.s)

  ### create bootstrap indices
  dfb <- create_bootstrap(df[, yn], times = control$ntree)

  if (control$parallel) {
    if (is.null(control$nCore)) {
      nCore <- parallel::detectCores()-1
    } else nCore <- control$nCore
    doParallel::registerDoParallel(nCore)
    outForest <- foreach::foreach(i=1:control$ntree,
                         .export = c('grow_tree_f', 'find_split_f'),
                         .packages = c('dplyr', 'coin', 'lazyeval', 'partykit')) %dopar% {
                           b.in <- dfb[, i]
                           dfbi <- df[b.in, ]
                           nodes <- grow_tree_f(id = 1L, xn, yn, dfbi, weights, minsplit = control$minsplit, maxdepth = control$maxdepth,
                                                mtry = mtry, alpha = control$alpha, balance.sample = control$balance.sample,
                                                bonferroni = control$bonferroni,split.criterion = control$split.criterion,
                                                var.select.criterion = control$var.select.criterion, var.select.test = control$var.select.test,
                                                minbucket.c = control$minbucket.c, minbucket.t = control$minbucket.t,
                                                yon, tn, classLevel, treatLevel, XClasses, yc, df.s, b.in)
                           allids <- partykit::nodeids(nodes)
                           tids <- partykit::nodeids(nodes, terminal = TRUE)
                           ntids <- allids[!allids %in% tids]
                           node.stats <- partykit::nodeapply(nodes, ids = allids, FUN = function(n) (n$info)$node.stats)
                           fitted <- partykit::fitted_node(nodes, data = data[xn])
                           t.node.stats <- partykit::nodeapply(nodes, ids = tids, FUN = function(n) (n$info)$node.stats)
                           t.node.stats <- data.frame(node = tids, uplift =  vapply(t.node.stats, function(x) x[, 1], numeric(1)))
                           resp <- (dplyr::inner_join(data.frame(node = fitted), t.node.stats, by = 'node'))$uplift
                           varids <- unlist(partykit::nodeapply(nodes, ids = ntids, FUN = function(x) partykit::varid_split(partykit::split_node(x))))
                           split.val <- unlist(partykit::nodeapply(nodes, ids = ntids, FUN = function(n) (n$info)$split.val))
                           split.val <- data.frame(varnames = xn[varids], varids = varids, split.val = split.val, stringsAsFactors = FALSE)
                           outTree <- partykit::party(node = nodes,
                                                      data = data[xn],
                                                      fitted = data.frame("(fitted)" = fitted,
                                                                          "(response)" = resp,
                                                                          check.names = FALSE),
                                                      terms = terms(formula),
                                                      info = list(node.stats = node.stats,
                                                                  split.val = split.val,
                                                      yon = yon,
                                                      tn = tn,
                                                      classLevel = classLevel,
                                                      treatLevel = treatLevel,
                                                      call = call))
                           outTree <- partykit::as.constparty(outTree)
                           class(outTree) <- append("utree", class(outTree))
                           outTree
                         }
    doParallel::stopImplicitCluster()
  } else {
    outForest <- foreach::foreach(i=1:control$ntree) %do% {
      b.in <- dfb[, i]
      dfbi <- df[b.in, ]
      nodes <- grow_tree_f(id = 1L, xn, yn, dfbi, weights, minsplit = control$minsplit, maxdepth = control$maxdepth,
                           mtry = mtry, alpha = control$alpha, balance.sample = control$balance.sample,
                           bonferroni = control$bonferroni,split.criterion = control$split.criterion,
                           var.select.criterion = control$var.select.criterion, var.select.test = control$var.select.test,
                           minbucket.c = control$minbucket.c, minbucket.t = control$minbucket.t,
                           yon, tn, classLevel, treatLevel, XClasses, yc, df.s, b.in)
      allids <- partykit::nodeids(nodes)
      tids <- partykit::nodeids(nodes, terminal = TRUE)
      ntids <- allids[!allids %in% tids]
      node.stats <- partykit::nodeapply(nodes, ids = allids, FUN = function(n) (n$info)$node.stats)
      fitted <- partykit::fitted_node(nodes, data = data[xn])
      t.node.stats <- partykit::nodeapply(nodes, ids = tids, FUN = function(n) (n$info)$node.stats)
      t.node.stats <- data.frame(node = tids, uplift =  vapply(t.node.stats, function(x) x[, 1], numeric(1)))
      resp <- (dplyr::inner_join(data.frame(node = fitted), t.node.stats, by = 'node'))$uplift
      varids <- unlist(partykit::nodeapply(nodes, ids = ntids, FUN = function(x) partykit::varid_split(partykit::split_node(x))))
      split.val <- unlist(partykit::nodeapply(nodes, ids = ntids, FUN = function(n) (n$info)$split.val))
      split.val <- data.frame(varnames = xn[varids], varids = varids, split.val = split.val, stringsAsFactors = FALSE)

      outTree <- partykit::party(node = nodes,
                                 data = data[xn],
                                 fitted = data.frame("(fitted)" = fitted,
                                                     "(response)" = resp,
                                                     check.names = FALSE),
                                 terms = terms(formula),
                                 info = list(node.stats = node.stats,
                                             split.val = split.val,
                                             yon = yon,
                                             tn = tn,
                                             classLevel = classLevel,
                                             treatLevel = treatLevel,
                                             call = call))
      outTree <- partykit::as.constparty(outTree)
      class(outTree) <- append("utree", class(outTree))
      outTree


    }
  }
  lsOut <- list(forest = outForest,
                ntree = control$ntree,
                mtry = mtry,
                responseLabel = yon,
                treatmentLabel = tn,
                classLevel = classLevel,
                treatLevel = treatLevel,
                nCore = ifelse(control$parallel, nCore, NA),
                call = call)
  class(lsOut) <- "uforest"
  lsOut
}

grow_tree_f <- function(id, xn, yn, data, weights, minsplit, maxdepth, mtry, alpha, balance.sample,
                        bonferroni, split.criterion, var.select.criterion, var.select.test, minbucket.c, minbucket.t,
                        yon, tn, classLevel, treatLevel, XClasses, yc, df.s, b.in, cenv = NULL) {

  node.stats <- get_node_uplift(data[weights == 1L, ], yon, tn, classLevel, treatLevel)

  if (is.null(cenv)) {
    cenv <- new.env()
    depth <- 0
  } else {
    depth <- get("depth", envir = cenv)
    if (depth >= maxdepth)
      return(partykit::partynode(id = id, info = list(split.val = NULL, node.stats = node.stats)))
  }

  tt <- ifelse(data[weights == 1L, tn] == treatLevel, 1L, -1L)
  inbag <- (balance_sample(tt, balance.sample, 1L))$inbag

  if (sum(weights) < minsplit | (yc && length(levels(data[weights == 1L, yn][inbag][drop=TRUE])) < 2))
    return(partykit::partynode(id = id, info = list(split.val = NULL, node.stats = node.stats)))

  sp <- find_split_f(xn, yn, data, weights, inbag, df.s, b.in,mtry, alpha, balance.sample,
                     bonferroni, split.criterion, var.select.criterion,
                     var.select.test, minbucket.c, minbucket.t,
                     yon, tn, classLevel, treatLevel, XClasses)

  if (is.null(sp)) return(partykit::partynode(id = id, info = list(split.val = NULL, node.stats = node.stats)))
  kidids <- partykit::kidids_split(sp, data = data[xn])
  ## set up all daugther nodes
  kids <- vector(mode = "list", length = max(kidids, na.rm = TRUE))
  for (kidid in 1L:length(kids)) {
    ## select observations for current node
    w <- weights
    w[kidids != kidid] <- 0
    assign("depth", depth + 1, envir = cenv)
    ## get next node id
    if (kidid > 1) {
      myid <- max(partykit::nodeids(kids[[kidid - 1]]))
    } else {
      myid <- id
    }
    ## start recursion on this daugther node
    kids[[kidid]] <- grow_tree_f(id = as.integer(myid + 1), xn, yn, data, w, minsplit, maxdepth, mtry, alpha,
                                 balance.sample, bonferroni, split.criterion, var.select.criterion, var.select.test,
                                 minbucket.c, minbucket.t, yon, tn, classLevel, treatLevel, XClasses, yc, df.s, b.in,
                                 cenv = cenv)
  }
  ## return nodes
  partykit::partynode(id = as.integer(id), split = sp, kids = kids,
                      info = list(split.val = partykit::info_split(sp)$split.val, node.stats = node.stats))
}


find_split_f <- function(xn,
                         yn,
                         data,
                         weights,
                         inbag,
                         df.s,
                         b.in,
                         mtry,
                         alpha,
                         balance.sample,
                         bonferroni,
                         split.criterion,
                         var.select.criterion,
                         var.select.test,
                         minbucket.c,
                         minbucket.t,
                         yon,
                         tn,
                         classLevel,
                         treatLevel,
                         XClasses)


{
  data <- data[weights == 1L, ]
  if (mtry < Inf) {
    inputs <- !(logical(length(xn)))
    mtry <- min(length(xn), mtry)
    ### length(xn) == 1 will lead to sample.int instead of sample; see partykit::ctree
    resample <- function(x, ...) x[sample.int(length(x), ...)]
    s <- resample(which(inputs), mtry)
    inp <- logical(length(inputs))
    inp[s] <- TRUE
  } else {
    inp <- !(logical(length(xn)))
  }
  xstat <- vapply(1:length(xn), function(i)
    select_split_var(yn, xn[i], inp[i], data[inbag, ], XClasses[i], var.select.criterion, var.select.test), numeric(1))
  if (all(is.na(xstat))) return(NULL)
  if (bonferroni) xstat <- stats::p.adjust(xstat, method = "bonferroni") # xstat * length(xn)
  names(xstat) <- xn
  split.var <- as.integer(which_is_min(xstat))
  min.xstat <- as.numeric(xstat[split.var])
  split.var.n <- xn[split.var]
  if (min.xstat > alpha) return(NULL)
  var.type <- XClasses[split.var] #class(df[, 1])

  if (var.type == "numeric") {
    nodeObsInd <- b.in[which(weights == 1L)]
    df <- dplyr::distinct(dplyr::inner_join(df.s[[split.var.n]],
                                            dplyr::data_frame(index = b.in), by = 'index'), index, .keep_all = TRUE)
    df <- dplyr::inner_join(df, dplyr::data_frame(index = nodeObsInd), by = 'index')[, 1:3]
  } else {
      df <- data[c(split.var.n, yon, tn)]
  }
  if (var.type == "numeric") {
    sindex <- NULL
    df <- dplyr::group_by_(df, split.var.n) %>%
      dplyr::summarize_(y11 = lazyeval::interp(~sum(v1 == classLevel & v2 == treatLevel),
                                               v1 = as.name(yon), v2 = as.name(tn)),
                        y10 = lazyeval::interp(~sum(v1 == classLevel & v2 != treatLevel),
                                               v1 = as.name(yon), v2 = as.name(tn)),
                        n0 = lazyeval::interp(~sum(v2 != treatLevel), v2 = as.name(tn)),
                        n1 = lazyeval::interp(~sum(v2 == treatLevel), v2 = as.name(tn)))
  } else {
    sbreaks <- NULL
    lev <- levels(df[, 1][drop = TRUE])
    lev.all <- levels(df[, 1])
    if (length(lev) == 2) {
      sindex <- ifelse(lev.all %in% lev, cumsum(lev.all %in% lev), NA)
    } else {
      df <- dplyr::group_by_(df, split.var.n) %>%
        dplyr::summarize_(y11 = lazyeval::interp(~sum(v1 == classLevel & v2 == treatLevel),
                                                 v1 = as.name(yon), v2 = as.name(tn)),
                          y10 = lazyeval::interp(~sum(v1 == classLevel & v2 != treatLevel),
                                                 v1 = as.name(yon), v2 = as.name(tn)),
                          n1 = lazyeval::interp(~sum(v2 == treatLevel), v2 = as.name(tn)),
                          n0 = lazyeval::interp(~sum(v2 != treatLevel), v2 = as.name(tn))) %>%
        dplyr::mutate(ate = y11/n1 - y10/n0) %>% #y11 is always the reference treatment
        dplyr::arrange(ate)
    }
  }

  if (var.type == "numeric" || length(lev) > 2) {

    df <- dplyr::mutate(df,
                        n0.l = cumsum(n0),
                        n0.r = sum(n0) - n0.l,
                        n1.l = cumsum(n1),
                        n1.r = sum(n1) - n1.l,
                        py11.l = cumsum(y11) / n1.l,
                        py10.l = cumsum(y10) / n0.l,
                        py11.r = (sum(y11) - cumsum(y11)) / n1.r,
                        py10.r = (sum(y10) - cumsum(y10)) / n0.r,
                        splitOK = ifelse(n0.l > minbucket.c & n0.r > minbucket.c & n1.l > minbucket.t &
                                           n1.r > minbucket.t, 1L, 0L))

    split.val <- get_split_val(n0=df$n0, n1=df$n1, n0.l=df$n0.l, n0.r=df$n0.r,
                               n1.l=df$n1.l, n1.r=df$n1.r, y11=df$y11, y10=df$y10,
                               py11.l=df$py11.l, py10.l=df$py10.l, py11.r=df$py11.r,
                               py10.r=df$py10.r, sc = split.criterion)

    split.val[is.na(split.val) | is.infinite(split.val)] <- 0
    if (all(split.val * df$splitOK <= 0)) return(NULL)
    split.val.pos <- which_is_max(split.val * df$splitOK)

    if (var.type == "numeric") {
      sbreaks <- as.numeric(df[split.val.pos, xn[split.var]])
    } else {
      left.lev <- as.character((as.data.frame(df[, xn[split.var]]))[, 1])[1L:split.val.pos]
      sindex <- !(levels(df[[xn[split.var]]]) %in% left.lev)
      sindex[!(levels(df[[xn[split.var]]]) %in% lev)] <- NA_integer_
      sindex <- sindex - min(sindex, na.rm = TRUE) + 1L
    }
  }

  partykit::partysplit(varid = split.var,
                       breaks = sbreaks,
                       index = sindex,
                       info = list(split.val = split.val[split.val.pos]))
}


#'@rdname uforest
#'@method print uforest
#'@export print.uforest

print.uforest <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  cat("Uplift Random Forest:", "\n", sep="")
  cat("\n")
  cat("Number of trees: ", x$ntree, "\n", sep="")
  cat("Number of randomly selected inputs: ", x$mtry, "\n", sep="")
  cat("Average tree depth: ", mean(vapply(x$forest, depth, numeric(1))), "\n", sep="")
  cat("Response variable: ", x$responseLabel, "\n", sep="")
  cat("Treatment indicator: ", x$treatmentLabel, "\n", sep="")
  cat("classLevel: ", x$classLevel, "\n", sep="")
  cat("treatLevel: ", x$treatLevel, "\n", sep="")
  if (!is.null(x$nCore))  cat("Number of cores: ", x$nCore, "\n", sep="")
  cat("\n")
  cat("\n")
  invisible(x)
}

#'Predict method for uforest fits.
#'
#'Obtains predictions from a fitted \code{uforest} object.
#'
#'@param object A fitted object inheriting from \code{uforest}.
#'@param newdata A data frame to predict on.
#'@param ntrees Number of trees used in prediction. All fitted trees are
#'  used by default.
#'@param type The type of predictions required. Only uplift predictions are
#'  allowed.
#'@param agg.fun The function used to combine the predictions from the
#'  individual trees. Must be functions that operate on a vector and produce a
#'  single value, as \code{"mean"} and \code{"median"} do. It may be a
#'  user-defined function.
#'@param \dots Additional arguments passed to \code{agg.fun}.
#'
#'@return A numeric vector of predictions.
#'
#'@export predict.uforest
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@seealso \code{\link{uforest}}.

predict.uforest <- function(object, newdata, ntrees = NULL, type = "uplift", agg.fun = "mean", ...) {
  if (type != "uplift") {
    stop("uplift: only \"uplift\" predictions are available.")
  } else type <- "response"
  fst <- object$forest
  if (is.null(ntrees)) {
    ntrees <- length(fst)
  } else if (ntrees > length(fst)) stop("uplift: Number of trees cannot exceeded number fitted.")
  predMatrix <- vapply(1:ntrees, function(i)
                       partykit:::predict.party(object = fst[[i]], newdata = newdata, type = type),
                       numeric(nrow(newdata)))
  apply(predMatrix, 1, eval(parse(text = agg.fun)), ...)
}

#'@rdname var_importance
#'@method var_importance uforest
#'@export var_importance.uforest

var_importance.uforest <- function(x, type = "I", valid.data = NULL, error.fun = "sel")
{
  if (!inherits(x, "uforest"))
    stop("uplift: x must be a 'uforest' class object")
  fst <- x$forest
  ntrees <- length(fst)
  imp <- lapply(1:ntrees,
                function(i) var_importance.utree(fst[[i]], type = type, valid.data = valid.data, error.fun = error.fun))
  imp <- do.call("rbind", imp)
  imp <- as.data.frame(dplyr::group_by(imp, varnames) %>% dplyr::summarize(importance = mean(importance, na.rm = TRUE)) %>%
                         dplyr::arrange(desc(importance)))
  imp
}


