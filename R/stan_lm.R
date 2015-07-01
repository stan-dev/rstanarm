# This file is part of rstanarm.
# Copyright 2013 Stan Development Team
# rstanarm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# rstanarm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with rstanarm.  If not, see <http://www.gnu.org/licenses/>.

#' @rdname stan_glm
#' @export
#' 
stan_lm <- function(formula, data, weights, subset,
                    na.action = NULL, start = NULL, offset = NULL, 
                    model = TRUE, x = FALSE, y = TRUE, contrasts = NULL,
                    prior = normal(), prior.for.intercept = normal(),
                    prior.options = prior_options(), ...) {
  
  if ("family" %in% names(list(...)))
    warning("extra argument ", sQuote("family"), 
            " is disregarded. Use stan_glm.", domain = NA)
  
  mf <- call <- match.call()
  mf[["family"]] <- gaussian(link = "identity")
  mf[[1L]] <- as.name("stan_glm")
  fit <- eval(mf, parent.frame())
  fit$call <- call
  fit
}



