#' @import purrr
#' @import stats
#' @import Rcpp
#' @import RcppArmadillo
#' @importFrom magrittr %>%
#' @details
#' Linear Regression with Little Bag of Bootstraps
"_PACKAGE"


## quiets concerns of R CMD check re: the .'s that appear in pipelines
# from https://github.com/jennybc/googlesheets/blob/master/R/googlesheets.R
utils::globalVariables(c("."))


#' @export
blblm <- function(formula, data, m = 10, B = 5000, num_core = 1) {
  # browser()
  data_list <- split_data(data, m)
  mechine_core <- detectCores()
  if (num_core == 1) {
    estimates <- map(
      data_list,
      ~ lm_each_subsample(formula = formula, data = ., n = nrow(data), B = B))
    res <- list(estimates = estimates, formula = formula)
  } else {
    # detect num of core
    if ((num_core >= mechine_core) | (num_core <= 0)) {
      run_core <- mechine_core
    } else {
      run_core <- num_core
    }
    print(glue::glue("all core: {mechine_core}, ues {run_core} cores in this function"))
    # start parallel
    nrow_data <<- nrow(data)
    B <<- B
    formula_in_this <<- formula
    cl <- makeCluster(run_core)
    clusterExport(cl, c("lm_each_subsample","lm_each_boot", "lm1", "blbcoef", "blbsigma",
                        "B", "formula_in_this", "nrow_data"))
    estimates <- pblapply(X = data_list,
                          FUN = function(temp_data){
                            lm_each_subsample(formula = formula_in_this,
                                              data = temp_data,
                                              n = nrow_data,
                                              B = B)}, cl = cl)

    stopCluster(cl)
    res <- list(estimates = estimates, formula = formula)
  }

  class(res) <- "blblm"
  invisible(res)
}



#' split data into m parts of approximated equal sizes
#' @export
split_data <- function(data, m) {
  idx <- sample.int(m, nrow(data), replace = TRUE)
  data %>% split(idx)
}


#' compute the estimates
#' @export
lm_each_subsample <- function(formula, data, n, B) {
  replicate(B, lm_each_boot(formula, data, n), simplify = FALSE)
}


#' compute the regression estimates for a blb dataset
#' @export
lm_each_boot <- function(formula, data, n) {
  freqs <- rmultinom(1, n, rep(1, nrow(data)))
  lm1(formula, data, freqs)
}


#' estimate the regression estimates based on given the number of repetitions
#' @export
lm1 <- function(formula, data, freqs) {
  # drop the original closure of formula,
  # otherwise the formula will pick a wront variable from the global scope.
  environment(formula) <- environment()
  fit <- lm(formula, data, weights = freqs)
  list(coef = blbcoef(fit), sigma = blbsigma(fit))
}


#' compute the coefficients from fit
#' @export
blbcoef <- function(fit) {
  coef(fit)
}


#' compute sigma from fit
#' @export
blbsigma <- function(fit) {
  p <- fit$rank
  y <- model.extract(fit$model, "response")
  e <- fitted(fit) - y
  w <- fit$weights
  sqrt(sum(w * (e^2)) / (sum(w) - p))
}


#' @export
#' @method print blblm
print.blblm <- function(x, ...) {
  cat("blblm model:", capture.output(x$formula))
  cat("\n")
}


#' @export
#' @method sigma blblm
sigma.blblm <- function(object, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  sigma <- mean(map_dbl(est, ~ mean(map_dbl(., "sigma"))))
  if (confidence) {
    alpha <- 1 - 0.95
    limits <- est %>%
      map_mean(~ quantile(map_dbl(., "sigma"), c(alpha / 2, 1 - alpha / 2))) %>%
      set_names(NULL)
    return(c(sigma = sigma, lwr = limits[1], upr = limits[2]))
  } else {
    return(sigma)
  }
}

#' @export
#' @method coef blblm
coef.blblm <- function(object, ...) {
  est <- object$estimates
  map_mean(est, ~ map_cbind(., "coef") %>% rowMeans())
}


#' @export
#' @method confint blblm
confint.blblm <- function(object, parm = NULL, level = 0.95, ...) {
  if (is.null(parm)) {
    parm <- attr(terms(fit$formula), "term.labels")
  }
  alpha <- 1 - level
  est <- object$estimates
  out <- map_rbind(parm, function(p) {
    map_mean(est, ~ map_dbl(., list("coef", p)) %>% quantile(c(alpha / 2, 1 - alpha / 2)))
  })
  if (is.vector(out)) {
    out <- as.matrix(t(out))
  }
  dimnames(out)[[1]] <- parm
  out
}

#' @export
#' @method predict blblm
predict.blblm <- function(object, new_data, confidence = FALSE, level = 0.95, ...) {
  est <- object$estimates
  X <- model.matrix(reformulate(attr(terms(object$formula), "term.labels")), new_data)
  if (confidence) {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>%
      apply(1, mean_lwr_upr, level = level) %>%
      t())
  } else {
    map_mean(est, ~ map_cbind(., ~ X %*% .$coef) %>% rowMeans())
  }
}

#' @export
mean_lwr_upr <- function(x, level = 0.95) {
  alpha <- 1 - level
  c(fit = mean(x), quantile(x, c(alpha / 2, 1 - alpha / 2)) %>% set_names(c("lwr", "upr")))
}

#' @export
map_mean <- function(.x, .f, ...) {
  (map(.x, .f, ...) %>% reduce(`+`)) / length(.x)
}

#' @export
map_cbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(cbind)
}

#' @export
map_rbind <- function(.x, .f, ...) {
  map(.x, .f, ...) %>% reduce(rbind)
}


################################################################################
# cpp part



#' @export
my_lm_with_cpp <- function(formula, data) {
  Rcpp::sourceCpp(code = '
                  #include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
Rcpp::List fastLm(const arma::mat& X, const arma::colvec& y) {
  int n = X.n_rows, k = X.n_cols;

  arma::colvec coef = arma::solve(X, y);    // fit model y ~ X
  arma::colvec res  = y - X*coef;           // residuals

  // std.errors of coefficients
  double s2 = std::inner_product(res.begin(), res.end(), res.begin(), 0.0)/(n - k);

  arma::colvec std_err = arma::sqrt(s2 * arma::diagvec(arma::pinv(arma::trans(X)*X)));

  return Rcpp::List::create(Rcpp::Named("coefficients") = coef,
                            Rcpp::Named("stderr")       = std_err,
                            Rcpp::Named("df.residual")  = n - k);
}
                  ')
  my_formula <- formula

  y_and_x <- unlist(strsplit(as.character(my_formula), split = '~'))

  y_name <- y_and_x[2]
  x_name <- y_and_x[3]
  data <- mtcars
  y_value <- data[ , y_name, drop = TRUE]
  x_matrix <- model.matrix(as.formula(paste0("~", x_name)), data = data)

  my_fit <- fastLm(x_matrix, y_value)

  my_fit$"fitted.values" <- x_matrix %*% my_fit$coefficients
  my_fit$"residuals" <- y_value - c(my_fit$"fitted.values")
  my_fit$"x_part" <- as.formula(paste0("~", x_name))
  return(my_fit)

}



#' @export
predict_my_lm_with_cpp <- function(object, new_data){
  x_part <- model.matrix(object$x_part, new_data)
  predict_value <- x_part %*% object[["coefficients"]]
  return(predict_value)
}

#' @export
confit.my_lm_with_cpp <- function(object, confidence = FALSE, level = 0.95) {
  estimate_value <- object$coefficients
  stderror_value <- object$stderr
  rownames(estimate_value) <- colnames(object$x_matrix)
  rownames(stderror_value) <- colnames(object$x_matrix)

  if (confidence) {

    output_matrix <- matrix(data = 0, nrow = length(estimate_value), ncol = 4)
    output_matrix[, 1] <- estimate_value
    output_matrix[, 2] <- stderror_value
    output_matrix[, 3] <- estimate_value + stderror_value * qnorm(level)
    output_matrix[, 4] <- estimate_value - stderror_value * qnorm(level)
    colnames(output_matrix) <- c("fit_value", "stderr", "uwr", "lwr")
    rownames(output_matrix) <- colnames(object$x_matrix)

  } else {
    output_matrix <- cbind(estimate_value, stderror_value)
    colnames(output_matrix) <- c("fit", "stderr")
  }

  return(output_matrix)
}


#' @export
hello <- function() {"this is written by cpp"}




#' @export
lm_in_one_file <- function(temp_file_path, formula) {
  # browser()
  data_in_r <- read.csv(temp_file_path)
  lm_temp <- lm(formula, data = data_in_r)
  return(lm_temp)
}



#' @export
cal_csv <- function(file_path, formula) {
  # each function
  # browser()
  all_lm_result <- list()
  cat("start\n")

  for (j_ in seq_along(file_path)) {
    all_lm_result[[j_]] <- lm_in_one_file(file_path[j_], formula)
    cat(j_, "\n")
  }

  cat("done\n")
  return(all_lm_result)
}





