#' Plots
#' 
#' All models fit using \code{algorithm='sampling'} are compatible with a 
#' variety of plotting functions. Each function returns at least one
#' \code{\link[ggplot2]{ggplot}} object that can be customized further using the
#' \pkg{ggplot2} package.
#' 
#' @name plots
#' 
#' @section Plotting functions:
#' \describe{
#' \item{Posterior predictive checks}{\code{\link{ppcheck}}}
#' \item{Posterior intervals and point estimates}{\code{\link[rstan]{stan_plot}}}
#' \item{Traceplots}{\code{\link[rstan]{stan_trace}}}
#' \item{Histograms}{\code{\link[rstan]{stan_hist}}}
#' \item{Kernel density estimates}{\code{\link[rstan]{stan_dens}}}
#' \item{Scatterplots}{\code{\link[rstan]{stan_scat}}}
#' \item{Diagnostics for Hamiltonian Monte Carlo and the No-U-Turn Sampler}{\code{\link[rstan]{stan_diag}}}
#' \item{Rhat}{\code{\link[rstan]{stan_rhat}}}
#' \item{Ratio of effective sample size to total posterior sample size}{\code{\link[rstan]{stan_ess}}}
#' \item{Ratio of Monte Carlo standard error to posterior standard deviation}{\code{\link[rstan]{stan_mcse}}}
#' \item{Autocorrelation}{\code{\link[rstan]{stan_ac}}}
#' }
#' 
#' @seealso The \code{\link[=shinystan]{ShinyStan}} graphical user interface.
#' 
#' @examples 
#' \dontrun{
#' data("clouds", package = "HSAUR3")
#' f <- rainfall ~ seeding * (sne + cloudcover + prewetness + echomotion) + time
#' fit <- stan_lm(f, data = clouds, prior = R2(location = 0.25), 
#'                cores = 4, seed = 12345)
#'                
#' # stan_plot: posterior intervals and point estimates
#' stan_plot(fit, ci_level = 0.8)
#' stan_plot(fit, pars = c("prewetness", "echomotionstationary"), show_density = TRUE)               
#' 
#' # posterior predictive checks (see ?ppcheck for more details and examples)
#' ppcheck(fit, check = "distributions")
#' ppcheck(fit, check = "distributions", overlay = TRUE)
#' ppcheck(fit, check = "residuals")
#' ppcheck(fit, check = "test", test = sd)
#' 
#' # traceplot
#' (trace <- stan_trace(fit, pars = "(Intercept)"))
#' trace + scale_color_discrete()
#' 
#' # distributions 
#' stan_hist(fit, fill = "skyblue") + ggtitle("Example Plot")
#' stan_dens(fit, pars = c("sne", "cloudcover"), separate_chains = TRUE, alpha = 0.1)
#' 
#' # scatterplot
#' stan_scat(fit, pars = c("sne", "cloudcover"))
#' }
NULL