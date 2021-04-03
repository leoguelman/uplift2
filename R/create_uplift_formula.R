#'Create an uplift formula.
#'
#'\code{create_uplift_formula} is a helper function to create a formula object
#'as required by most functions in the uplift package.
#'
#'@name create_uplift_formula
#'
#'@aliases create_uplift_formula
#'
#'@export create_uplift_formula
#'
#'@param x A character vector of predictor names.
#'@param y A character vector of the response name.
#'@param trt A character vector of the treatment name.
#'@param env The environment associated with the result, the calling environment
#'  by default.
#'
#'@return A formula object.
#'
#'@author Leo Guelman \email{leo.guelman@@gmail.com}
#'
#'@examples
#'set.seed(1)
#'df <- sim_uplift(n = 100, p = 20, response = "binary")
#'f <- create_uplift_formula(names(df)[-c(1:3)], "y", "T")
#'class(f)
#'environment(f) # the callling environment

create_uplift_formula <- function(x, y, trt, env = parent.frame()) {

  if (missing(x)) stop("uplift: argument \"x\" is missing.")
  if (missing(y)) stop("uplift: argument \"y\" is missing.")
  if (missing(trt)) stop("uplift: argument \"trt\" is missing.")
  if (!(is.character(y) && length(y) == 1)) stop("uplift: argument \"y\" should be a character of length 1.")
  if (!(is.character(trt) && length(trt) == 1)) stop("uplift: argument \"trt\" should be a character of length 1.")
  if (!(is.character(x))) stop("uplift: argument \"x\" should be a character vector.")

  mform <- as.formula(paste(y, '~', 'trt(', trt, ')+',
                            paste(x, sep = '', collapse = "+")), env = env)

}
