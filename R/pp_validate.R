# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2016 Trustees of Columbia University
# Copyright (C) 2005 Samantha Cook
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

#' Model validation via simulation
#' 
#' The \code{pp_validate} function is based on the methods described in
#' Cook, Gelman, and Rubin (2006) for validating software developed to fit
#' particular Bayesian models. Here we take the perspective that models 
#' themselves are software and thus it is useful to apply this validation 
#' approach to individual models.
#' 
#' @export
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @param nreps The number of replications to be performed.
#' @param seed A seed passed to Stan to use when refitting the model.
#' @param ... Arguments passed to \code{\link{geom_point}} to control the
#'   appearance of the plot.
#'   
#' @details 
#' We repeat \code{nreps} times the process of simulating parameters and data 
#' from the model and refitting the model to this simulated data. For each of
#' the \code{nreps} replications we do the following:
#' \enumerate{
#' \item Refit the model but \emph{without} conditioning on the data (setting 
#' \code{prior_PD=TRUE}), obtaining draws \eqn{\theta^{true}}{\theta_true} 
#' from the \emph{prior} distribution of the model parameters.
#' \item Given \eqn{\theta^{true}}{\theta_true}, simulate data \eqn{y^\ast}{y*} 
#' from the \emph{prior} predictive distribution (calling 
#' \code{\link{posterior_predict}} on the fitted model object obtained in step
#' 1).
#' \item Fit the model to the simulated outcome \eqn{y^\ast}{y*}, obtaining 
#' parameters \eqn{\theta^{post}}{\theta_post}.
#' }
#' For any individual parameter, the quantile of the "true" parameter value with
#' respect to its posterior distribution \emph{should} be uniformly distributed.
#' The validation procedure entails looking for deviations from uniformity by 
#' computing statistics for a test that the quantiles are uniformly distributed.
#' The absolute values of the computed  test statistics are plotted for batches 
#' of parameters (e.g., non-varying coefficients are grouped into a batch called
#' "beta", parameters that vary by group level are in batches named for the 
#' grouping variable, etc.). See Cook, Gelman, and Rubin (2006) for more details
#' on the validation procedure.
#' 
#' @return A ggplot object that can be further customized using the 
#'   \pkg{ggplot2} package.
#' 
#' @references 
#' Cook, S., Gelman, A., and Rubin, D. 
#' (2006). Validation of software for Bayesian models using posterior quantiles.
#' \emph{Journal of Computational and Graphical Statistics}. 15(3), 675--692.
#' 
#' @seealso 
#' \code{\link{pp_check}} for graphical posterior predictive checks and 
#' \code{\link{posterior_predict}} to draw from the posterior predictive 
#' distribution.
#' 
#' 
#' @examples 
#' \dontrun{
#' pp_validate(example_model)
#' }
#' 
pp_validate <- function(object, nreps = 20, seed = 12345, ...) {
  # based on Samantha Cook's BayesValidate::validate
  quant <- function(draws) {
    n <- length(draws)
    rank_theta <- c(1:n)[order(draws) == 1] - 1
    quants <- (rank_theta + 0.5) / n
    return(quants)
  }
  
  if (!used.sampling(object))
    STOP_sampling_only("pp_validate")
  if (nreps < 2)
    stop("'nreps' must be at least 2.")
  
  dims <- object$stanfit@par_dims[c("alpha", "beta", "b", "dispersion")]
  dims <- dims[!sapply(dims, is.null)]
  dims <- sapply(dims, prod)
  dims <- dims[dims > 0]
  if ("b" %in% names(dims)) {
    mark <- which(names(dims) == "b")
    vals <- sapply(ranef(object), function(x) length(as.matrix(x)))
    dims <- append(dims, values = vals, after = mark)
    dims <- dims[-mark]
  }
  batches <- dims
  params_batch <- names(dims)
  num_batches <- length(batches)
  num_params <- sum(dims)
  batch_ind <- rep(0, num_batches + 1)
  plot_batch <- rep(1, batches[1])
  for (i in 1:num_batches) 
    batch_ind[i+1] <- batch_ind[i] + batches[i]
  for (i in 2:num_batches) 
    plot_batch <- c(plot_batch, rep(i, batches[i]))
  quantile_theta <- matrix(NA_real_, nrow = nreps, ncol = num_params + num_batches)
  for (reps in 1:nreps) {
    post <- suppressWarnings(update(object, prior_PD = TRUE, seed = seed,
                                    warmup = 1000, iter = 1000 + 2, chains = 1))
    theta_true <- as.matrix(post)[1,, drop = FALSE]
    data_rep <- posterior_predict(post)[1, ]
    mf <- model.frame(object)
    if (NCOL(mf[, 1]) == 2) { # binomial models
      mf[, 1] <- c(data_rep)
      colnames(mf)[1] <- colnames(get_y(object))[1]
    } else {
      mf[, 1] <- c(data_rep)
    }
    theta_draws <- as.matrix(update(object, data = mf, seed = seed))
    if (!is.null(batches)){
      for (i in 1:num_batches) {
        if (batches[i] > 1) {
          sel <- (batch_ind[i]+1):batch_ind[(i+1)]
          theta_draws <-  cbind(theta_draws, 
                                apply(theta_draws[, sel], 1, mean))
          theta_true <- c(theta_true, mean(theta_true[sel]))
        } else {
          theta_draws <- cbind(theta_draws, theta_draws[, (batch_ind[i]+1)])
          theta_true <- c(theta_true, theta_true[(batch_ind[i]+1)])
        }
      }
    }
    theta_draws <- rbind(theta_true, theta_draws)
    quantile_theta[reps, ] <- apply(theta_draws, 2, quant)
  }
  quantile_trans <- (apply(quantile_theta, 2, qnorm))^2
  q_trans <- apply(quantile_trans, 2, sum) 
  p_vals <- pchisq(q_trans, df = nreps, lower.tail = FALSE)
  z_stats <- abs(qnorm(p_vals))
  if (is.null(batches)) {
    adj_min_p <- num_params * min(p_vals)
  } else {
    z_batch <- z_stats[(num_params + 1):length(p_vals)]
    p_batch <- p_vals[(num_params + 1):length(p_vals)]
    adj_min_p <- num_batches * min(p_batch)
  }
  
  upper_lim <- max(max(z_stats + 1), 3.5)
  plotdata <- data.frame(x = z_batch, y = params_batch)
  defaults <- list(shape = 21, fill = .PP_FILL, color = "black", 
                   size = 2.5, alpha = 1)
  geom_args <- set_geom_args(defaults, ...)
  ggplot(plotdata, aes_string(x = "x", y = "y")) + 
    geom_segment(aes_string(x = "0", xend = "x", y = "y", yend = "y")) +
    do.call("geom_point", geom_args) + 
    scale_x_continuous(limits = c(0, upper_lim), expand = c(0, 0)) + 
    labs(y = NULL, x = expression("Absolute " * z[theta] * " Statistics")) + 
    pp_check_theme(no_y = FALSE) + 
    theme(panel.grid.major.x = element_line(size = 0.1, color = "gray"))
}  

# If we wanted to return the actual pvals: 
#
# if (is.null(batches)){
#   return(list(p_vals = p_vals, adj_min_p = adj_min_p))
# } else {
#   if (length(z_batch) == num_params)
#     return(list(p_batch = p_batch, adj_min_p = adj_min_p)) 
#   else 
#     return(list(p_vals = p_vals[1:num_params], p_batch = p_batch, 
#                 adj_min_p = adj_min_p))
# }
