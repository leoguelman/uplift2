utree <- function(x, ...) UseMethod("utree")

utree.default <- function(x, ...) stop("uplift: 'x' should be a formula object")

#'Control for uplift trees.
#'
#'Various parameters that control aspects of the \code{utree} fit.
#'
#'@name utree_control
#'
#'@aliases utree_control
#'
#'@export utree_control
#'
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
#'  possible values are: \code{"uplift"} (default), \code{"kld"}
#'  (Kullback-Leibler divergence), \code{"ed"} (Euclidean divergence), or
#'  \code{"l1d"} (L1-norm divergence). See details in Guelman et al. (2015).
#'@param maxdepth Maximum depth of the tree. The default \code{maxdepth = Inf}
#'  means that no restrictions are applied to tree sizes.
#'@param mtry Number of input variables randomly sampled as candidates at each
#'  node. The default \code{mtry = Inf} means that no random selection takes
#'  place.
#'
#'@return A list.
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'

utree_control <- function(minsplit = 40L, minbucket.t = 20L, minbucket.c = 20L,
                          var.select.criterion = "pvalue", var.select.test = "asymptotic()",
                          alpha = 0.05, bonferroni = FALSE, balance.sample = "undersample",
                          split.criterion = "uplift", maxdepth = Inf, mtry = Inf)
{
  list(minsplit = minsplit, minbucket.t = minbucket.t, minbucket.c = minbucket.c,
       var.select.criterion = var.select.criterion, var.select.test = var.select.test,
       alpha = alpha, bonferroni = bonferroni,  balance.sample = balance.sample,
       split.criterion = split.criterion, maxdepth = maxdepth, mtry = mtry)
}

#'Fitting uplift trees.
#'
#'\code{utree} implements recursive partitioning for uplift modeling.
#'
#'Roughly, the algorithm works as follows:\enumerate{ \item For each terminal node in the tree
#'we test the global null hypothesis of no interaction effect between the treatment indicator and any of
#'the covariates. Stop if this hypothesis cannot be rejected. Otherwise, select the input variable
#'with strongest interaction effect. The interaction effect is measured by a p-value corresponding
#'to an asymptotic or permutation test (Strasser and Weber, 1999) for the partial null hypothesis of independence
#'between each covariate and a transformed response. Specifically, the response is
#'transformed so the impact of the covariate on the response has a causal interpretation
#'for the treatment effect (see details in Guelman et al. 2015)
#'\item Implement a binary split in the selected input variable.
#'\item Recursively repeate the two steps above.}
#'
#'Function \code{nodeprune} is not yet implemented for \code{utree} objects.
#'
#'@name utree
#'
#'@aliases utree
#'
#'@export utree utree.default utree.formula
#'
#'@import partykit
#'
#'@param formula A model formula of the form y ~ x1 + ....+ xn + trt(), where
#'  the left-hand side corresponds to the observed response, the right-hand side
#'  corresponds to the predictors, and 'trt' is the special expression to mark
#'  the treatment term. At the moment, \code{utree} only handles binary responses.
#'@param data A data frame in which to interpret the variables named in the
#'  formula.
#'@param na.action A missing-data filter function.
#'@param classLevel A character string for the class of interest. Defaults to
#'  the last level of the factor.
#'@param treatLevel A character string for the treatment level of interest.
#'  Defaults to the last level of the treatment factor.
#'@param control A list with control parameters, see \code{\link{utree_control}}.
#'@param \dots Arguments passed to \code{\link{utree_control}}.
#'@param x An object of class \code{"utree"}
#'
#'@return An object of class \code{"utree"}.
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@seealso \code{\link{plot.utree}}.
#'
#'@references
#'
#'Guelman, L., Guillen, M., and Perez-Marin A.M. (2015). "A decision support
#'framework to implement optimal personalized marketing interventions." Decision
#'Support Systems, Vol. 72, pp. 24--32.
#'
#'Hothorn, T., Hornik, K. and Zeileis, A. (2006). "Unbiased recursive partitioning:
#' A conditional inference framework". Journal of Computational and Graphical Statistics,
#' 15(3): 651--674.
#'
#'Rzepakowski, Piotr and Jaroszewicz, Szymon. (2011). "Decision trees for uplift modeling
#'with single and multiple treatments". Knowledge and Information Systems, 32(2) 303--327.
#'
#'Strasser, H. and Weber, C. (1999). "On the asymptotic theory of permutation statistics".
#'Mathematical Methods of Statistics, 8: 220--250.
#'
#'Su, X., Tsai, C.-L., Wang, H., Nickerson, D. M. and Li, B. (2009). "Subgroup Analysis via Recursive Partitioning".
#'Journal of Machine Learning Research 10, 141--158.
#'
#'@examples
#'
#'set.seed(1)
#'df <- sim_uplift(n = 1000, p = 50, response = "binary")
#'form <- create_uplift_formula(x = names(df)[-c(1:3)], y = "y", trt = "T")
#'fit <- utree(form, data = df, maxdepth = 3)
#'fit

utree.formula <- function(formula,
                          data,
                          na.action,
                          classLevel = NULL,
                          treatLevel = NULL,
                          control = utree_control(...),
                          ...
)
{

  call <- match.call()

  stopifnot(control$minsplit > 2L, control$minbucket.c > 1L, control$minbucket.t > 1L)

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
    stop("uplift: utree can only handle binary response variables of class factor.")

  ### pre-sorting numeric predictors
  df.s <- lapply(1:length(xn), function(i) sort_fun(df, i))
  names(df.s) <- xn
  df.s <- Filter(Negate(is.null), df.s)

  nodes <- grow_tree(id = 1L, xn, yn, df, weights, minsplit = control$minsplit, maxdepth = control$maxdepth,
                     mtry = control$mtry, alpha = control$alpha, balance.sample = control$balance.sample,
                     bonferroni = control$bonferroni,split.criterion = control$split.criterion,
                     var.select.criterion = control$var.select.criterion, var.select.test = control$var.select.test,
                     minbucket.c = control$minbucket.c, minbucket.t = control$minbucket.t,
                     yon, tn, classLevel, treatLevel, XClasses, yc, df.s)

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

grow_tree <- function(id, xn, yn, data, weights, minsplit, maxdepth, mtry, alpha, balance.sample,
                      bonferroni, split.criterion, var.select.criterion, var.select.test, minbucket.c, minbucket.t,
                      yon, tn, classLevel, treatLevel, XClasses, yc, df.s, cenv = NULL) {

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

  sp <- find_split(xn, yn, data, weights, inbag, df.s, mtry, alpha, balance.sample,
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
    kids[[kidid]] <- grow_tree(id = as.integer(myid + 1), xn, yn, data, w, minsplit, maxdepth, mtry, alpha,
                               balance.sample, bonferroni, split.criterion, var.select.criterion, var.select.test,
                               minbucket.c, minbucket.t, yon, tn, classLevel, treatLevel, XClasses, yc, df.s,
                               cenv = cenv)
  }
  ## return nodes
  partykit::partynode(id = as.integer(id), split = sp, kids = kids,
                      info = list(split.val = partykit::info_split(sp)$split.val, node.stats = node.stats))
}


find_split <- function(xn,
                       yn,
                       data,
                       weights,
                       inbag,
                       df.s,
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
    nodeObsInd <- which(weights == 1L)
    #df <- df.s[[split.var.n]][df.s[[split.var.n]]$index %in% nodeObsInd, ][, 1:3]
    df <- dplyr::inner_join(df.s[[split.var.n]], dplyr::data_frame(nodeObsInd),
                            by = c("index" = "nodeObsInd"))[, 1:3]
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

select_split_var <- function(y,
                             x,
                             inp,
                             data,
                             var.class,
                             var.select.criterion,
                             var.select.test)
{
  if (var.class == "factor" && length(levels(data[, x][drop=TRUE])) < 2) return(NA)
  if (var.class == "numeric" && length(unique(data[, x])) < 2) return(NA)
  if (!inp) return(NA)
  switch(var.select.criterion,
         pvalue = as.numeric(coin::pvalue(coin::independence_test(new("IndependenceProblem", y = data[y], x = data[x]),
                                                                  distribution = eval(parse(text = paste0("coin::", var.select.test, sep =""))))))
  )
}


get_split_val <- function(n0, n1, n0.l, n0.r, n1.l, n1.r, y11, y10, py11.l, py10.l, py11.r, py10.r, sc = NULL)
{

  if (sc == "uplift") {
    sse <-  n0.l * py10.l * (1 - py10.l) +
      n1.l * py11.l * (1 - py11.l) +
      n0.r * py10.r * (1 - py10.r) +
      n1.r * py11.r * (1 - py11.r)
    split.val <-  ((sum(n0 + n1) - 4) * ((py11.l - py10.l) - (py11.r - py10.r))^2) /
      (sse * (1/n0.l + 1/n0.r + 1/n1.l + 1/n1.r))
  } else if (sc %in% c("kld", "ed", "l1d")) {
    .sc <- get(sc, mode = "function", envir = parent.frame())
    d.node <- .sc(sum(y11)/sum(n1), sum(y10)/sum(n0))
    d.l <- mapply(.sc, py11.l, py10.l)
    d.r <- mapply(.sc, py11.r, py10.r)
    d.lr <- ((n0.l + n1.l)/sum(n0 + n1)) * d.l + ((n0.r + n1.r)/sum(n0 + n1)) * d.r
    d.gain <- d.lr - d.node
    T.norm <-  mapply(.sc, n1.l/sum(n1), n0.l/sum(n0))
    if (sc == "kld") {
      T <- ep(sum(n1)/sum(n0+n1))
      T1 <- vapply(n1.l/sum(n1), ep, numeric(1))
      T0 <- vapply(n0.l/sum(n0), ep, numeric(1))
    } else {
      T <- gini(sum(n1)/sum(n0+n1))
      T1 <- vapply(n1.l/sum(n1), gini, numeric(1))
      T0 <- vapply(n0.l/sum(n0), gini, numeric(1))
    }
    norm <- T.norm * T + T1 * sum(n1)/sum(n0+n1) + T0 * sum(n0)/sum(n0+n1) + 0.5
    split.val <- d.gain/norm
  }
  split.val
}


### See "Elements of Information Theory" (Cover)

ep <- function(p, base = 2) {
  if (is.na(p)) return(NA)
  if (p==0 | p==1) return(0L)
  - 1 * (p * log(p, base = base) +  (1-p) * log((1-p), base = base))
}

kld <- function(p, q, base = 2) {
  if (is.na(p) | is.na(q)) return(NA)
  if (p==0 | p==1) return(0L)
  else if (p==0 && q ==0 | p==1 && q==1) return(0L)
  else if (p>0 && (q==0 | q==1)) q <- 1E-11
  p * log(p/q, base = base) + (1-p) * log((1-p)/(1-q), base = base)

}

ed <- function(p, q) {
  if (is.na(p) | is.na(q)) return(NA)
  (p-q)^2
}

l1d <- function(p, q) {
  if (is.na(p) | is.na(q)) return(NA)
  abs(p-q)
}

gini <- function(p) {
  2 * p * (1 - p)
}

# Follows nnet::which.is.max

which_is_min <- function(x, na.rm = TRUE)
{
  y <- seq_along(x)[x == min(x, na.rm = na.rm)]
  if (na.rm) y <- y[!is.na(y)]
  if (length(y) > 1L)
    sample(y, 1L)
  else y
}

which_is_max <- function(x, na.rm = TRUE)
{
  y <- seq_along(x)[x == max(x, na.rm = na.rm)]
  if (na.rm) y <- y[!is.na(y)]
  if (length(y) > 1L)
    sample(y, 1L)
  else y
}

get_node_uplift <- function(data, yon, tn, classLevel, treatLevel)
{
  dplyr::summarize_(data,
                    y11 = lazyeval::interp(~sum(v1 == classLevel & v2 == treatLevel),
                                           v1 = as.name(yon), v2 = as.name(tn)),
                    y10 = lazyeval::interp(~sum(v1 == classLevel & v2 != treatLevel),
                                           v1 = as.name(yon), v2 = as.name(tn)),
                    n.treat = lazyeval::interp(~sum(v2 == treatLevel), v2 = as.name(tn)),
                    n.ctrl = lazyeval::interp(~sum(v2 != treatLevel), v2 = as.name(tn))) %>%
    dplyr::mutate(uplift = y11/n.treat - y10/n.ctrl) %>%
    dplyr::select(c(5L, 3L, 4L))
}


get_pred_error <- function(data, yon, tn, classLevel, treatLevel, error.fun)
{
  ddf <-  dplyr::group_by(data, data[, 4], data[, 3]) %>%
    dplyr::summarize_(y11 = lazyeval::interp(~sum(v1 == classLevel & v2 == treatLevel),
                                             v1 = as.name(yon), v2 = as.name(tn)),
                      y10 = lazyeval::interp(~sum(v1 == classLevel & v2 != treatLevel),
                                             v1 = as.name(yon), v2 = as.name(tn)),
                      n.treat = lazyeval::interp(~sum(v2 == treatLevel), v2 = as.name(tn)),
                      n.ctrl = lazyeval::interp(~sum(v2 != treatLevel), v2 = as.name(tn))) %>%
    dplyr::mutate(actual.u = y11/n.treat - y10/n.ctrl)
  err <- switch(error.fun,
                ### error weighted by obs counts
                sel = sum(((ddf[, 2] - ddf[, 7])^2) * (ddf[, 5] + ddf[, 6]), na.rm = TRUE),
                abs = sum((abs(ddf[, 2] - ddf[, 7])) * (ddf[, 5] + ddf[, 6]), na.rm = TRUE))
  err
}

sort_fun <- function(data, x)
{
  if (is.numeric(data[, x+3])) {
    d <- data[, c(x+3, 1, 2)]
    index <- order(d[, 1]) # the index gives original position of the obs
    d <- d[index, ]
    d$index <- index
    d
  } else NULL
}

#'@rdname utree
#'@method print utree
#'@export print.utree

print.utree <- function(x, ...)
{
  cat("Call:\n")
  print((x$`info`)$call)
  cat("\n")
  cat("Uplift tree structure:", "\n", sep="")
  cat("\n")
  node.stats <- ((x$`info`)$node.stats)
  node.uplift <- vapply(node.stats, function(x) x[, 1], numeric(1))
  node.n.treat <- vapply(node.stats, function(x) x[, 2], numeric(1))
  node.n.ctrl <- vapply(node.stats, function(x) x[, 3], numeric(1))
  partykit:::print.constparty(x, header = FALSE, FUN = function(...)
    paste0("Leaf info: uplift =",
           format(node.uplift, digits = getOption("digits") - 4),
           "; n.treat = ", node.n.treat,
           "; n.ctrl = ", node.n.ctrl))
  invisible(x)
}

#'Predict method for utree fits.
#'
#'Obtains predictions from a fitted \code{utree} object.
#'
#'@param object A fitted object inheriting from \code{utree}.
#'@param newdata A data frame to predict on.
#'@param type The type of predictions required. The default is \code{"uplift"}
#'  for uplift predictions Alternatively, \code{"node"}  returns an integer
#'  vector of terminal node identifiers.
#'@param \dots Additional arguments passed to \code{partykit::predict.party}.
#'
#'@return A numeric vector of predictions.
#'
#'@export predict.utree
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@seealso \code{\link{utree}}.

predict.utree <- function(object, newdata, type = "uplift", ...) {
  if (!(type %in% c("uplift", "node"))) stop("uplift: type must be either \"uplift\" or \"node\".")
  if (type == "uplift") type <- "response"
  partykit:::predict.party(object = object, newdata = newdata, type = type, ...)

}

#'@rdname utree
#'@method nodeprune utree
#'@export nodeprune.utree

nodeprune.utree <- function(x, ...) stop("uplift: 'nodeprune' not yet implemented for utree objects")

#'Visualization of uplift trees.
#'
#'plot method for class \code{utree}.
#'
#'@param x An object of class \code{utree}.
#'@param \dots Arguments passed to \code{partykit::plot.constparty}.
#'
#'@export plot.utree
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@examples
#'set.seed(1)
#'df <- sim_uplift(n = 1000, p = 50, response = "binary")
#'form <- create_uplift_formula(x = names(df)[-c(1:3)], y = "y", trt = "T")
#'fit <- utree(form, data = df, maxdepth = 3)
#'plot(fit, main = "uplift tree", gp = gpar(cex = 0.5))

plot.utree <- function(x, ...) {

  xl <- as.list(x$node)
  for(i in 1:length(xl)) {
    xl[[i]]$info$node.stats <- paste(
      format(names(xl[[i]]$info$node.stats)), "=",
      format(xl[[i]]$info$node.stats, digits = getOption("digits") - 4))
  }
  x$node <- partykit::as.partynode(xl)

  partykit:::plot.constparty(x = x,
                             terminal_panel = partykit::node_terminal,
                             tp_args = list(FUN = function(node) node$node.stats), ...)
}


var_importance <- function(x, ...) UseMethod("var_importance")

#'Variable importance for uplift trees and uplift random forest.
#'
#'This is the extractor function for variable importance measures as produced by
#'\code{utree} and \code{uforest}.
#'
#'For type I, the measure of importance given to a predictor is the sum of the
#'values given by the split-criterion produced over all internal nodes for which
#'it was chosen as the splitting variable. For uplift random forest, this
#'relative influence measure is naturally extended by averaging the importance
#'for each variable over the collection of trees. For type II, variable
#'importance is measured based on an independent validation sample, with the aim
#'of quantifying the prediction strength of each variable. This is achieved by
#'first measuring the prediction accuracy on this validation sample.
#'Subsequently, the values for the jth variable are randomly permuted, and the
#'accuracy again computed. The decrease in accuracy as a result of this
#'permutation is the importance attributed to the jth variable.The accuracy is
#'measured by the squared-error or absolute error between the predicted and true
#'uplift on each terminal node of the tree.
#'
#'
#'@name var_importance
#'
#'@aliases var_importance
#'
#'@export var_importance var_importance.utree
#'
#'@param x An object of class \code{"utree"} or \code{"uforest"}
#'@param type Either \code{"I"} or \code{"II"}, specifying the type of
#'  importance measure. See details.
#'@param valid.data For \code{type = "II"}, importance is measured based on a
#'  validation data frame, which must be provided.
#'@param error.fun The prediction error used to compute variable importance when
#'  \code{type = "II"}. Possible values are \code{"sel"} for squared-error loss
#'  (default), or \code{"abs"} for absolute loss. See details.
#'
#'@return A data frame with the variable importance.
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@examples
#'set.seed(1)
#'df <- sim_uplift(n = 1000, p = 50, response = "binary")
#'form <- create_uplift_formula(x = names(df)[-c(1:3)], y = "y", trt = "T")
#'fit <- utree(form, data = df, maxdepth = 3)
#'var_importance(fit)

var_importance.utree <- function(x, type = "I", valid.data = NULL, error.fun = "sel")
{
  if (!inherits(x, "utree"))
    stop("uplift: x must be a 'utree' class object")
  if (!(type %in% c("I", "II")))
    stop("uplift: type must be either \"I\", or \"II\".")
  if (!error.fun %in% c("sel", "abs"))
    stop("uplift: error.fun must be either \"sel\", or \"abs\".")
  if (type == "I") {
    sv <- (x$`info`)$split.val
    imp <- dplyr::group_by(sv, varnames) %>%
           dplyr::summarize(importance = sum(split.val)) %>%
           dplyr::arrange(desc(importance))
    imp$importance <- 100 * imp$importance / sum(imp$importance)
    imp <- as.data.frame(imp)
  } else {
    if (is.null(valid.data)) stop("uplift: valid.data must be provided with type \"II\".")
    varnames <- unique((x$`info`)$split.val[, 1])
    yon <- (x$`info`)$yon
    tn <- (x$`info`)$tn
    classLevel <- (x$`info`)$classLevel
    treatLevel <- (x$`info`)$treatLevel
    newnames <- names(valid.data)
    var.check <- varnames %in% newnames
    if (!(all(var.check)))
      stop("uplift: variables ", paste0(varnames[!var.check], collapse=", "), " missing in valid.data.")
    if (!(any(newnames %in% yon)))
      stop("uplift: valid.data must contain the response", yon, ".")
    if (!(any(newnames %in% tn)))
      stop("uplift: valid.data must contain the treatment indicator", tn, ".")
    valid.data[, tn] <- as.factor(valid.data[, tn])
    valid.data0 <- data.frame(valid.data[c(yon, tn)], pred.u = predict(x, valid.data),
                              pred.n = predict(x, valid.data, type = "node"))
    err0 <- get_pred_error(valid.data0, yon, tn, classLevel, treatLevel, error.fun = error.fun)
    nVars <- length(varnames)
    err <- numeric(nVars)
    for (i in 1:nVars) {
      valid.data.p <- valid.data[c(yon, tn, varnames)]
      valid.data.p[, varnames[i]] <- sample(valid.data.p[, varnames[i]])
      valid.data.p <- data.frame(valid.data.p[c(yon, tn)], pred.u = predict(x, valid.data.p),
                                 pred.n =  predict(x, valid.data.p, type = "node"))
      err[i] <- get_pred_error(valid.data.p, yon, tn, classLevel, treatLevel, error.fun = error.fun)
    }
    imp <- data.frame(varnames, importance = err/err0)
    imp$importance <- 100 * imp$importance / sum(imp$importance)
    imp <- dplyr::arrange(imp, desc(importance))
  }
  imp
}
