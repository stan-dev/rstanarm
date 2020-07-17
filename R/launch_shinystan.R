# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016, 2017 Trustees of Columbia University
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

#' Using the ShinyStan GUI with rstanarm models
#' 
#' The ShinyStan interface provides visual and numerical summaries of model
#' parameters and convergence diagnostics.
#' 
#' @aliases launch_shinystan
#' @export
#' 
#' @templateVar stanregArg object
#' @template args-stanreg-object
#' @inheritParams shinystan::launch_shinystan
#' @param ppd Should \pkg{rstanarm} draw from the posterior predictive 
#'   distribution before launching ShinyStan? The default is \code{TRUE}, 
#'   although for very large objects it can be convenient to set it to 
#'   \code{FALSE} as drawing from the posterior predictive distribution can be 
#'   time consuming. If \code{ppd} is \code{TRUE} then graphical posterior 
#'   predictive checks are available when ShinyStan is launched.
#' @param seed Passed to \link[=pp_check]{pp_check} if 
#'   \code{ppd} is \code{TRUE}.
#' @param model_name,note Optional arguments passed to
#'   \code{\link[shinystan]{as.shinystan}}.
#'   
#' @details The \code{\link[shinystan]{launch_shinystan}} function will accept a
#'   \code{\link[=stanreg-objects]{stanreg}} object as input. Currently, almost 
#'   any model fit using one of \pkg{rstanarm}'s model-fitting functions can be 
#'   used with ShinyStan. The only exception is that ShinyStan does not 
#'   currently support \pkg{rstanarm} models fit using 
#'   \code{algorithm='optimizing'}. See the 
#'   \pkg{\link[=shinystan-package]{shinystan}} package documentation for more 
#'   information.
#'   
#' @section Faster launch times:
#' For some \pkg{rstanarm} models ShinyStan may take a very long time to launch.
#' If this is the case with one of your models you may be able to speed up
#' \code{launch_shinystan} in one of several ways:
#' \describe{
#'   \item{Prevent ShinyStan from preparing graphical posterior predictive
#'   checks:}{
#'   When used with a \code{\link[=stanreg-objects]{stanreg}} object 
#'   (\pkg{rstanarm} model object) ShinyStan will draw from the posterior 
#'   predictive distribution and prepare graphical posterior predictive checks 
#'   before launching. That way when you go to the PPcheck page the plots are 
#'   immediately available. This can be time consuming for models fit to very
#'   large datasets and you can prevent this behavior by creating a shinystan
#'   object before calling \code{launch_shinystan}. To do this use 
#'   \code{\link[shinystan]{as.shinystan}} with optional argument \code{ppd} set
#'   to \code{FALSE} (see the Examples section below). When you then launch
#'   ShinyStan and go to the PPcheck page the plots will no longer be 
#'   automatically generated and you will be presented with the standard
#'   interface requiring you to first specify the appropriate \eqn{y} and
#'   \eqn{yrep}, which can be done for many but not all \pkg{rstanarm} models.
#'   }
#'   \item{Use a shinystan object:}{
#'   Even if you don't want to prevent ShinyStan from preparing graphical
#'   posterior predictive checks, first creating a shinystan object using
#'   \code{\link[shinystan]{as.shinystan}} can reduce \emph{future} launch
#'   times. That is, \code{launch_shinystan(sso)} will be faster than
#'   \code{launch_shinystan(fit)}, where \code{sso} is a shinystan object and
#'   \code{fit} is a stanreg object. It still may take some time for 
#'   \code{as.shinystan} to create \code{sso} initially, but each time you
#'   subsequently call \code{launch_shinystan(sso)} it will reuse \code{sso}
#'   instead of internally creating a shinystan object every time. See the
#'   Examples section below.}
#' }
#' 
#' @template reference-bayesvis
#' @template reference-muth
#'   
#' @examples
#' \dontrun{
#' if (!exists("example_model")) example(example_model) 
#' 
#' # Launch the ShinyStan app without saving the resulting shinystan object
#' if (interactive()) launch_shinystan(example_model)
#' 
#' # Launch the ShinyStan app (saving resulting shinystan object as sso)
#' if (interactive()) sso <- launch_shinystan(example_model)
#' 
#' # First create shinystan object then call launch_shinystan
#' sso <- shinystan::as.shinystan(example_model)
#' if (interactive()) launch_shinystan(sso)
#' 
#' # Prevent ShinyStan from preparing graphical posterior predictive checks that
#' # can be time consuming. example_model is small enough that it won't matter
#' # much here but in general this can help speed up launch_shinystan
#' sso <- shinystan::as.shinystan(example_model, ppd = FALSE)
#' if (interactive()) launch_shinystan(sso)
#' }
#' 
launch_shinystan.stanreg <-
  function(object,
           ppd = TRUE, 
           seed = 1234, 
           model_name = NULL, 
           note = NULL, 
           rstudio = getOption("shinystan.rstudio"), 
           ...) {
    sso <-
      shinystan::as.shinystan(
        object,
        ppd = ppd,
        seed = seed,
        model_name = model_name,
        note = note
      )
    shinystan::launch_shinystan(sso, rstudio = rstudio, ...)
  }
