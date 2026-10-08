#' Print a CaPri fit
#'
#' @param x An object of class \code{capri}.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{x}.
#' @export
print.capri <- function(x, ...) {
  cat("CaPri fit\n")
  cat("---------\n")
  cat("Input type:", x$input_type, "\n")
  cat("Posterior probability of causal effect q_delta:", round(x$q_delta, 4), "\n")
  cat("Estimated causal effect hat_gamma:", round(x$hat_gamma, 6), "\n")
  cat("Estimated pleiotropic effect mu_betaY:", round(x$mu_betaY, 6), "\n")
  cat("Number of instruments:", length(x$pip_c), "\n")
  cat("Number of iterations:", x$n_iter, "\n")

  invisible(x)
}

#' Summarize a CaPri fit
#'
#' @param object An object of class \code{capri}.
#' @param ... Additional arguments.
#'
#' @return Invisibly returns \code{object}.
#' @export
summary.capri <- function(object, ...) {
  print(object)

  cat("\nPleiotropy posterior probabilities:\n")
  print(summary(object$pip_c))

  cat("\nTop pleiotropic instruments by posterior probability:\n")
  ord <- order(object$pip_c, decreasing = TRUE)
  top <- utils::head(ord, min(10, length(ord)))

  out <- data.frame(
    instrument = top,
    pip_c = object$pip_c[top]
  )

  print(out, row.names = FALSE)

  invisible(object)
}

#' Extract coefficients from a CaPri fit
#'
#' @param object An object of class \code{capri}.
#' @param ... Additional arguments.
#'
#' @return A named numeric vector.
#' @export
coef.capri <- function(object, ...) {
  c(
    hat_gamma = object$hat_gamma,
    q_delta = object$q_delta,
    mu_betaY = object$mu_betaY
  )
}

#' Plot a CaPri fit
#'
#' @param x An object of class \code{capri}.
#' @param ... Additional arguments passed to \code{graphics::plot}.
#'
#' @return Invisibly returns \code{x}.
#' @export
plot.capri <- function(x, ...) {
  graphics::plot(
    x$pip_c,
    type = "h",
    ylim = c(0, 1),
    xlab = "Instrument index",
    ylab = "Posterior probability of pleiotropy",
    main = "CaPri pleiotropy probabilities",
    ...
  )

  invisible(x)
}
