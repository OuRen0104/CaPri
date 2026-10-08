#' Fit CaPri using summary statistics
#'
#' Fits the CaPri Mendelian randomization model using SNP-level summary
#' statistics. This summary-statistics version assumes independent variants,
#' corresponding to an identity LD matrix by default.
#'
#' @param muY Numeric vector of SNP-outcome summary association estimates.
#' @param betaX Numeric vector of SNP-exposure summary association estimates.
#' @param n Sample size used for the summary statistics.
#' @param allele_mode Character string specifying allele orientation.
#' If "as_is", summary statistics are used without modification.
#' If "minor", SNP effects are flipped so that effects correspond to
#' the minor allele direction.
#' @param eaf Numeric vector of effect allele frequencies. Required when
#' \code{allele_mode="minor"}.
#' @param pi_delta Prior probability of a causal effect.
#' @param pi_c Prior probability that an instrument is pleiotropic.
#' @param tau_gamma Prior precision for the causal effect.
#' @param tau_beta Prior precision for the pleiotropic effect.
#' @param sigma_Y Residual standard deviation of the outcome model.
#' @param max_iter Maximum number of variational inference iterations.
#' @param eps_theta Parameter convergence threshold.
#' @param patience Number of consecutive stable iterations required for convergence.
#' @param verbose Whether to print progress during variational inference.
#'
#' @return An object of class \code{capri}.
#'
#' @export
fit_capri_sumstats <- function(
    muY,
    betaX,
    n,
    allele_mode = c("as_is", "minor"),
    eaf = NULL,
    pi_delta = 0.01,
    pi_c = 0.5,
    tau_gamma = 1.0,
    tau_beta = 0,
    sigma_Y = 1.0,
    max_iter = 2000,
    eps_theta = 1e-6,
    patience = 5,
    verbose = 1
) {

  allele_mode <- match.arg(allele_mode)

  muY <- as.numeric(muY)
  betaX <- as.numeric(betaX)
  n <- as.numeric(n)

  if (length(muY) != length(betaX)) {
    stop("muY and betaX must have the same length.")
  }

  if (allele_mode == "minor") {
    if (is.null(eaf)) {
      stop(
        "eaf must be provided when allele_mode='minor'."
      )
    }

    eaf <- as.numeric(eaf)

    if (length(eaf) != length(betaX)) {
      stop(
        "eaf must have the same length as betaX."
      )
    }
    if (any(!is.finite(eaf)) ||
        any(eaf <= 0) ||
        any(eaf >= 1)) {
      stop(
        "eaf must contain values between 0 and 1."
      )
    }

    # Effect allele frequency > 0.5
    # means effect allele is the major allele.
    # Flip beta so effects correspond to minor allele.

    flip_idx <- eaf > 0.5

    betaX[flip_idx] <- -betaX[flip_idx]
    muY[flip_idx] <- -muY[flip_idx]
  }


  if (length(n) != 1 || 
      !is.finite(n) || 
      n <= 0) {
    stop("n must be a positive scalar.")
  }
  if (pi_delta <= 0 || pi_delta >= 1) {
    stop("pi_delta must be between 0 and 1.")
  }
  if (pi_c <= 0 || pi_c >= 1) {
    stop("pi_c must be between 0 and 1.")
  }
  if (tau_gamma <= 0) {
    stop("tau_gamma must be positive.")
  }


  K <- length(muY)
  eps <- 1e-16
  sigma2 <- as.numeric(sigma_Y^2)

  sigmoid <- function(x) {
    1 / (1 + exp(-x))
  }


  logit_pi_delta <- 
    log(pi_delta + eps) -
    log(1 - pi_delta + eps)


  logit_pi_c <- 
    log(pi_c + eps) -
    log(1 - pi_c + eps)



  xTx_hat <- n * sum(betaX^2)
  yTx_hat <- n * sum(betaX * muY)


  q_delta <- 0.5
  mu_gamma <- 0.0
  tau_gamma_star <- tau_gamma + 1.0


  p_c <- rep(pi_c, K)


  mu_beta <- 0.0
  tau_beta_star <- tau_beta + 1.0


  good_count <- 0L
  last <- NULL


  for (it in seq_len(max_iter)) {
    Eg <- q_delta * mu_gamma
    mu_r_beta <- muY - Eg * betaX
    c_hat <- as.numeric(p_c)
    GcTGc <- n * sum(c_hat)
    GcTr <- n * sum(c_hat * mu_r_beta)
    tau_beta_star <- tau_beta + GcTGc / sigma2
    mu_beta <- as.numeric(
      GcTr /
        (sigma2 * tau_beta_star)
    )



    var_beta <- 1 /
      (tau_beta_star + eps)
    Eb2 <- mu_beta^2 + var_beta

    for (k in seq_len(K)) {
      mu_r_k <- muY[k] -
        Eg * betaX[k]

      diff <- (
        mu_beta * n * mu_r_k -
          0.5 * Eb2 * n
      ) / sigma2


      p_c[k] <- as.numeric(
        sigmoid(
          logit_pi_c + diff
        )
      )

    }


    c_hat <- as.numeric(p_c)

    mu_r_gamma <-
      muY -
      mu_beta * c_hat

    Xtr_gamma <-
      n * sum(betaX * mu_r_gamma)

    tau_gamma_star <-
      tau_gamma +
      xTx_hat / sigma2

    mu_gamma <-
      as.numeric(
        Xtr_gamma /
          (sigma2 * tau_gamma_star)
      )

    logit_delta <-
      logit_pi_delta +
      0.5 *
      (
        log(
          (tau_gamma + eps) /
          (tau_gamma_star + eps)
        ) +
          tau_gamma_star * mu_gamma^2
      )

    q_delta <-
      as.numeric(
        sigmoid(logit_delta)
      )

    if (is.null(last)) {
      dtheta <- Inf
    } else {


      dtheta <- max(
        abs(q_delta - last$q_delta),
        abs(mu_gamma - last$mu_gamma),
        abs(mu_beta - last$mu_beta),
        max(abs(p_c - last$p_c)),
        abs(tau_beta_star - last$tau_beta_star),
        abs(tau_gamma_star - last$tau_gamma_star),
        abs(
          q_delta * mu_gamma -
          last$q_delta * last$mu_gamma
        )
      )

    }

    last <- list(
      q_delta = q_delta,
      mu_gamma = mu_gamma,
      mu_beta = mu_beta,
      p_c = p_c,
      tau_beta_star = tau_beta_star,
      tau_gamma_star = tau_gamma_star
    )

    if (dtheta < eps_theta) {
      good_count <- good_count + 1L
    } else {
      good_count <- 0L
    }


    if (verbose &&
        (it == 1 || it %% 10 == 0)) {


      cat(
        sprintf(
          paste0(
            "[sumstats iter %4d] ",
            "dtheta=%.2e  ",
            "good=%d/%d  ",
            "qdelta=%.4f  ",
            "hat_gamma=%.6f  ",
            "mu_betaY=%.6f  ",
            "mean(p)=%.4f\n"
          ),
          it,
          dtheta,
          good_count,
          patience,
          q_delta,
          q_delta * mu_gamma,
          mu_beta,
          mean(p_c)
        )
      )
    }



    if (good_count >= patience) {
      if(verbose){
        cat(
          sprintf(
            paste0(
              "Summary-statistics VI converged ",
              "at iter %d.\n"
            ),
            it
          )
        )
      }
      break
    }
  }


  out <- list(
    q_delta = q_delta,
    mu_gamma = mu_gamma,
    hat_gamma =
      q_delta * mu_gamma,
    mu_betaY = mu_beta,
    pip_c = p_c,
    tau_betaY_star =
      tau_beta_star,
    tau_gamma_star =
      tau_gamma_star,
    betaX = betaX,
    xTx_hat = xTx_hat,
    yTx_hat = yTx_hat,
    n_iter = it,
    input_type = "summary",
    allele_mode = allele_mode,
    call = match.call()

  )


  class(out) <- "capri"


  out

}
