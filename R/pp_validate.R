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
#' @param seed A seed passed to Stan to use when refitting the model to the
#'   simulated data.
#' @param ... Arguments passed to \code{\link{geom_point}} to control the
#'   appearance of the plot.
#'   
#' @details 
#' We repeat \code{nreps} times the process of simulating data and parameteres 
#' from the model and refitting the model to this simulated data. 
#' For each of the \code{nreps} replications we 
#' \enumerate{
#' \item Draw parameter values from the \emph{prior} predictive distribution.
#' \item Simulate the outcome \eqn{y} using these parameters.
#' \item Fit the model to the simulated outcome.
#' \item Check the posterior quantiles for the parameters in the model fit
#' to the simulated outcome. 
#' }
#' The posterior quantiles for these parameters \emph{should} follow a uniform
#' distribution, so we look for deviations from uniformity by computing
#' statistics for a test that the quantiles are uniformly distributed. The
#' absolute values of the computed z-statistics are plotted for batches of 
#' parameters (e.g., non-varying slope parameters are grouped into a batch 
#' called "beta").
#' 
#' @return A ggplot object that can be further customized using the 
#'   \pkg{ggplot2} package.
#' 
#' @references 
#' Cook, S., Gelman, A., and Rubin, D. 
#' (2006). Validation of software for Bayesian models using posterior quantiles.
#' \emph{Journal of Computational and Graphical Statistics}. 15(3), 675--692.
#' 
#' @examples 
#' \dontrun{
#' pp_validate(example_model, nreps = 5)
#' }
#' 
pp_validate <- function(object, nreps, seed = 12345, ...) {
  # based on Samantha Cook's BayesValidate::validate
  quant <- function(draws) {
    n <- length(draws)
    rank.theta <- c(1:n)[order(draws) == 1] - 1
    quantile.theta <- (rank.theta +.5) / n
    return(quantile.theta)
  }
  
  if (missing(nreps))
    stop("'nreps' must be specified")
  
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
  n.batch <- dims
  params.batch <- names(dims)
  num.batches <- length(n.batch)
  n.param <- sum(dims)
  batch.ind <- rep(0, num.batches + 1)
  plot.batch <- rep(1, n.batch[1])
  for (i in 1:num.batches) 
    batch.ind[i+1] <- batch.ind[i] + n.batch[i]
  for (i in 2:num.batches) 
    plot.batch <- c(plot.batch, rep(i, n.batch[i]))
  quantile.theta <- matrix(NA_real_, nrow = nreps, ncol = n.param + num.batches)
  for (reps in 1:nreps) {
    post <- suppressWarnings(update(object, prior_PD = TRUE, 
                 warmup = 1000, iter = 1000 + 2, chains = 1, seed = seed))
    theta.true <- as.matrix(post)[1,, drop = FALSE]
    data.rep <- posterior_predict(post)[1, ]
    mf <- model.frame(object)
    if (NCOL(mf[, 1]) == 2) { # binomial models
      mf[, 1] <- c(data.rep)
      colnames(mf)[1] <- colnames(get_y(object))[1]
    } else {
      mf[, 1] <- c(data.rep)
    }
    theta.draws <- as.matrix(update(object, data = mf, seed = seed))
    if (!is.null(n.batch)){
      for (i in 1:num.batches) {
        if (n.batch[i] > 1) {
          sel <- (batch.ind[i]+1):batch.ind[(i+1)]
          theta.draws <-  cbind(theta.draws, 
                                apply(theta.draws[, sel], 1, mean))
          theta.true <- c(theta.true, mean(theta.true[sel]))
        } else {
          theta.draws <- cbind(theta.draws, theta.draws[, (batch.ind[i]+1)])
          theta.true <- c(theta.true, theta.true[(batch.ind[i]+1)])
        }
      }
    }
    theta.draws <- rbind(theta.true, theta.draws)
    quantile.theta[reps, ] <- apply(theta.draws, 2, quant)
  }
  quantile.trans <- (apply(quantile.theta, 2, qnorm))^2
  q.trans <- apply(quantile.trans, 2, sum) 
  p.vals <- pchisq(q.trans, df = nreps, lower.tail = FALSE)
  z.stats <- abs(qnorm(p.vals))
  if (is.null(n.batch)) {
    adj.min.p <- n.param * min(p.vals)
  } else {
    z.batch <- z.stats[(n.param + 1):length(p.vals)]
    p.batch <- p.vals[(n.param + 1):length(p.vals)]
    adj.min.p <- num.batches * min(p.batch)
  }
  
  upper.lim <- max(max(z.stats + 1), 3.5)
  plotdata <- data.frame(x = z.batch, y = params.batch)
  defaults <- list(shape = 21, fill = .PP_FILL, color = "black", 
                   size = 2.5, alpha = 1)
  geom_args <- set_geom_args(defaults, ...)
  ggplot(plotdata, aes(x, y)) + 
    geom_segment(aes(x = 0, xend = x, y = y, yend = y)) +
    do.call("geom_point", geom_args) + 
    scale_x_continuous(limits = c(0, upper.lim), expand = c(0, 0)) + 
    labs(y = NULL, x = expression("Absolute " * z[theta] * " Statistics")) + 
    pp_check_theme(no_y = FALSE) + 
    theme(panel.grid.major.x = element_line(size = 0.1, color = "gray"))
  
  
  # if (is.null(n.batch)){
  #   return(list(p.vals = p.vals, adj.min.p = adj.min.p))
  # } else {
  #   if (length(z.batch) == n.param)
  #     return(list(p.batch = p.batch, adj.min.p = adj.min.p)) 
  #   else 
  #     return(list(p.vals = p.vals[1:n.param], p.batch = p.batch, 
  #                 adj.min.p = adj.min.p))
  # }
}
