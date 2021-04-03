score_profile <- function(x, ...) UseMethod("score_profile")

score_profile.default <- function(x, ...) stop("uplift: 'x' should be a formula object")

#'Profile deciles from a fitted model.
#'
#'This function can be used to profile the deciles from a fitted model. Given a
#'vector of numeric scores (fitted values) and predictors, it computes basic
#'summary statistics for each predictor by score quantile.
#'
#'This function ranks the variable supplied in the left-hand side of the model
#'formula and classifies it into groups with approximately the same number of observations.
#'It subsequently calls the function \code{tables::tabular} to compute the average
#'of each numeric predictor, and the distribution of each factor within each
#'group.
#'
#'@name score_profile
#'
#'@aliases score_profile
#'
#'@export score_profile score_profile.default score_profile.formula
#'
#'@import tables
#'
#'@param formula A formula expression of the form score ~ predictors, where the
#'  score represents the predictions from a fitted model.
#'@param data A data frame in which to interpret the variables named in the
#'  formula.
#'@param groups Number of groups of equal observations in which to partition the
#'  data set to show results. The default value is 10 (deciles).
#'@param statistic Functions that operate on a vector and produce a single
#'  value, as \code{mean} and \code{sd} do. It may be a user-defined function.
#'  To request several statistics, use the \code{+} operator. For example,
#'  \code{statistic = "mean + min + max"}. This argument only applies to numeric
#'  variables when \code{categorize = FALSE}. Factors are always shows as
#'  percentages within each group.
#'@param direction Possible values are \code{"D"} or \code{"I"}, for group
#'  number labels which are decreasing or increasing with the model score,
#'  respectively.
#'@param categorize Should numeric predictors be categorized at their quantiles?
#'@param nBins The number of bins created for numeric variables. The bins are
#'  created based on quantiles, with a default value of 4 (quartiles). Only
#'  applicable when \code{categorize=TRUE}.
#'@param continuous When \code{categorize=TRUE}, it specifies the threshold for
#'  when a numeric variable should be categorized at their quantiles, or at
#'  their unique values. When there are at least \code{continuous} unique
#'  values, bins are created based on quantiles. Otherwise, the variables is
#'  converted to factor with levels being equal to the variable's unique values.
#'@param digitsN Number of decimal places to show for numeric predictors.
#'@param digitsF Number of decimal places to show for factor predictors.
#'@param digitsB Number of digits used in formatting the breaks
#'@param groupVar A character string with the variable name in the data which
#'  holds the grouped predictions. If this argument is not null, groups of
#'  predictions are not created based on their quantiles but already declared
#'  from the named variable supplied to this argument.
#'@param excludeNA Should the results exclude observations with missing values
#'  in any of the variables named in the formula?
#'@param LaTex Should the function output LaTex code?
#'@param x A \code{score_profile} object.
#'@param \dots Additional arguments for the S3 methods.
#'
#'@return An object of class \code{score_profile}, which is a list with the
#'  following components: \itemize{ \item \code{data} The data frame containing
#'  the data used for plotting. \item \code{Table} An object of class
#'  \code{tabular} See \code{?tables::tabular} for details.}
#'
#'@author Leo Guelman \email{leo.guelman@@rbc.com}
#'
#'@seealso \code{\link{ggplot.score_profile}}.
#'
#'@examples
#'
#' ### Simulate some data
#' set.seed(123)
#' x1 <- rnorm(1000)
#' x2 <- rnorm(1000)
#' f1 <- sample(c(0, 1), 1000, replace = TRUE)
#' z <- 1 + 2 * x1 + 3 * x2  + f1
#' pr <- 1 / (1 + exp( -z))
#' y <- rbinom(1000, 1, pr)
#' df <- data.frame(y = y, x1 = x1, x2 = x2, f1 = factor(f1))
#' ### Fit model and get fitted values
#' Fitted <- fitted(glm(y ~ x1 + x2 + f1, data = df, family = "binomial"))
#' ### Profile deciles
#' score_profile(Fitted ~ x1 + x2 + f1, data = df, direction = "I")


score_profile.formula <- function(formula,
                                  data,
                                  groups = 10,
                                  statistic = "mean",
                                  direction = "D",
                                  categorize = TRUE,
                                  nBins = 4,
                                  continuous = 4,
                                  digitsN = NULL,
                                  digitsF = NULL,
                                  digitsB = NULL,
                                  groupVar = NULL,
                                  excludeNA = FALSE,
                                  LaTex = FALSE) {


  call <- match.call()
  if (!inherits(formula, "formula"))
    stop("uplift: method is only for formula objects")
  validDir <- charmatch(direction,
                        c("I", "D"))
  if (any(is.na(validDir)))
    stop("uplift: direction must be either 'I' or 'D'")
  mf <- match.call(expand.dots = FALSE)
  args <- match(c("formula", "data"),
                names(mf), 0L)
  mf <- mf[c(1L, args)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  dataCL <- attributes(mt)$dataClasses[-1] # exclude response
  if (!all(unique(dataCL) %in% c("numeric", "factor", "ordered")))
    stop("uplift: variable classes in formula must be either numeric, integer, factor or ordered")
  numVars <- unname(which(dataCL ==  "numeric"))
  facVars <- unname(which(dataCL %in% c("factor", "ordered")))
  termLabels <- attributes(mt)$term.labels
  numVarNames <- termLabels[numVars]
  facVarNames <- termLabels[facVars]
  respVarName <- names(attributes(mt)$dataClasses[1])
  attr(mt, "intercept") <- 0
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  if (!is.numeric(Y))
    stop("uplift: the left-hand side of the model formula must be a numeric vector")
  if (is.null(digitsN)) digitsN <- 2
  if (is.null(digitsF)) digitsF <- 1
  if(is.null(digitsB)) digitsB <- 6
  if (is.null(groupVar)) {
    if (direction == "D")
      Yval <- (-Y) else Yval <- Y
      Breaks <- unique(quantile(Yval, probs = seq(0, 1, 1 / groups)))
      BreaksLen <- length(Breaks)
      if (BreaksLen - 1 < groups)
        warning("uplift: Table created with ", (BreaksLen-1),
                " groups instead of ", groups)
      Group <- cut(Yval, breaks = Breaks, labels = 1:(BreaksLen-1), include.lowest = TRUE)
      respVarNameGrp <- paste(respVarName, "group", sep ="_")
  } else {
    Group <- as.factor(data[, groupVar])
    respVarNameGrp <- groupVar
  }
  if (categorize && length(numVars) != 0L)
    mf[numVarNames] <- data.frame(lapply(mf[numVarNames], create_bins_numeric,
                                         nBins = nBins, continuous = continuous,
                                         excludeNA = excludeNA, digitsB = digitsB))
  dframe <- data.frame(mf, Group)
  names(dframe)[ncol(dframe)] <- respVarNameGrp
  if (excludeNA) dframe <- na.omit(dframe)
  if (!categorize && length(numVars) != 0L) {
    t1 <- paste("+", paste(numVarNames, collapse = "+"))
  } else {
    t1 <- ""}
  if (categorize) {
    facVars <- c(facVars, numVars)
    facVarNames <- c(facVarNames, numVarNames)
    numVarNames <- NULL
  }
  if (length(facVars) != 0L) {
    t2 <- paste("+ Format(sprintf('%.", digitsF, "f')) * ((",
                paste("tables::Factor(", facVarNames, ")", sep = "",  collapse = " + "),
                ") * ",
                "(Pctn. = Percent('col')))", sep = "")} else {
                  t2 <- ""}
  tabForm <- paste("tables::tabular((n=1) + Format(sprintf('%.", digitsN, "f'))", " * ((",
                   respVarName,
                   t1, ") * (", statistic, "))", t2,
                   " ~ Justify(r) * (", respVarNameGrp," + 1),",
                   " data = dframe)", sep ="")
  if (LaTex) tabForm <- paste("latex(", tabForm, ")", sep ="")
  Table <- eval(parse(text = tabForm))
  lsOut <- list(data = dframe,
                Table = Table,
                call = call)
  class(lsOut) <- append("score_profile", class(lsOut))
  lsOut
}

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

#' @rdname score_profile
#' @method is score_profile
#' @export is.score_profile

is.score_profile <- function(x) inherits(x, "score_profile")


#' @rdname score_profile
#' @method print score_profile
#' @export print.score_profile

print.score_profile <- function(x, ...)
{
  print(x$Table)
}

#'Plot a \code{score_profile} object.
#'
#'ggplot method for class `score_profile'
#'
#'@param x An object of class \code{score_profile}.
#'@param i.var A numeric vector of indices of the variables to plot. The
#'  variables should be indexed in the same order that they appear in the
#'  initial inspect_balance formula. The default is to plot all variables.
#'@param n.type The type of plot for numeric variables. Boxplots are generated
#'  by default. For alternative summary statistics, use \code{summary}.
#'@param f.type The type of plot for categorical variables. 100 percent stacked
#'  columns are generated by default. The alternative option is \code{counts},
#'  which show stacked columns of counts.
#'@param statistic Functions that operate on a vector and produce a single
#'  value, as \code{mean} and \code{sd} do. It may be a user-defined function.
#'@param geom The geometric object to display the data. Argument is passed to
#'  \code{ggplot2::stat_summary} when \code{type = summary}. For \code{type =
#'  summary}, \code{geom = "point"}.
#'@param legend.position The position of legends \code{("left", "right",
#'  "bottom", "top")}. The default is \code{"right"}.
#'@param col Color of \code{geom}. The default is \code{"red"}.
#'@param size Size of \code{geom}. The default is \code{2}.
#'@param nrows Number of rows for plots.
#'@param ncols Number of columns for plots.
#'@param xlab Title for the x.
#'@param \dots Additional arguments passed to \code{ggplot2::stat_summary}.
#'
#'@export ggplot.score_profile
#'
#'@author Leo Guelman \email{leo.guelman@@rbc.com}
#'
#'@examples
#'
#'set.seed(123)
#'N <- 10000
#'eps <- rnorm(N)
#'age <- round(rnorm(N, 50, 10))
#'income <- rnorm(N, 60000, 10000) + 200 * age
#'gender <- gl(2, N/2,  labels = c("F", "M"))
#'insurance <- gl(4, N/4,  labels = c("HOME", "AUTO", "LIFE", "HEALTH"))
#'z <- 1e-01 + 0.1 * age - 1e-04 * income + 0.3 * (gender == "F") + eps
#'pr <- 1 / (1 + exp( -z))
#'purchase <- rbinom(N, 1, pr)
#'df <- data.frame(purchase, age, income, gender, insurance)
#'### Fit glm
#'pred <- fitted(glm(purchase ~ age + income + gender,
#'                   data = df, family = "binomial"))
#'profileForm <- pred ~ age + income + gender + insurance
#'prof1 <- score_profile(profileForm, data = df)
#'prof1
#'ggplot(prof1)

ggplot.score_profile <- function(x,
                                 i.var = NULL,
                                 n.type = "boxplot",
                                 f.type = "percent",
                                 statistic = "mean",
                                 geom = "point",
                                 legend.position = NULL,
                                 col = NULL,
                                 size = NULL,
                                 nrows = NULL,
                                 ncols = NULL,
                                 xlab = NULL,
                                 ...)
{

  if (!inherits(x, "score_profile"))
    stop("uplift: x must be a 'score_profile' class object")
  n.validType <- charmatch(n.type,
                           c("boxplot", "summary"))
  f.validType <- charmatch(f.type,
                           c("percent", "counts"))
  if (any(is.na(n.validType)))
    stop("uplift: not recognized plot 'n.type' argument")
  if (any(is.na(f.validType)))
    stop("uplift: not recognized plot 'f.type' argument")
  if (is.null(size)) size <- 2
  if (is.null(col)) col <- "red"
  if (is.null(legend.position)) legend.position <- "right"
  data <- x$data
  respVarNameGrp <- names(data)[ncol(data)]
  if (is.null(xlab)) xlab <- respVarNameGrp
  ylab <- paste("Percent within", xlab, sep = " ")
  if (is.null(xlab)) xlab <- respVarNameGrp
  varC <- vapply(data[-c(1L, ncol(data))], class, character(1))
  varL <- names(varC)
  nVars <- length(varL)
  if (is.null(i.var)) {
    i.var <- 1:nVars
  } else if (any(i.var > nVars)) {
    stop("uplift: \"i.var\" must be between 1 and ",nVars)
  }
  plots <- vector("list", length(i.var))
  for (i in seq_along(i.var)) {
    if (varC[i.var[i]] == "numeric") {
      if (n.type == "boxplot") {
        plots[[i]] <- ggplot(data, mapping = aes_string(x=respVarNameGrp, y=varL[i.var[i]])) +
          geom_boxplot() +
          stat_summary(fun.y = statistic, geom = "point", col = col, size = size, ...) + xlab(xlab)
      } else {
        plots[[i]] <- ggplot(data, mapping = aes_string(x=respVarNameGrp, y=varL[i.var[i]], group=1)) +
          stat_summary(fun.y = statistic, geom = geom, col = col, size = size, ...) + xlab(xlab)
      }
    } else {
      if (f.type == "percent") {
        plots[[i]] <- ggplot(data, aes_string(respVarNameGrp)) + geom_bar(aes_string(fill = varL[i.var[i]]), position = 'fill')  +
          scale_y_continuous(labels = scales::percent_format()) + xlab(xlab) + ylab(ylab) +
          theme(legend.position=legend.position)
      } else {
        plots[[i]] <- ggplot(data, aes_string(respVarNameGrp)) + geom_bar(aes_string(fill = varL[i.var[i]]))  +
          xlab(xlab) +
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


