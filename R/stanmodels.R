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

loadModule("stan_fit4bernoulli_mod", TRUE)
loadModule("stan_fit4binomial_mod", TRUE)
loadModule("stan_fit4continuous_mod", TRUE)
loadModule("stan_fit4count_mod", TRUE)
loadModule("stan_fit4lm_mod", TRUE)
loadModule("stan_fit4polr_mod", TRUE)

MODELS_HOME <- "exec"
if (!file.exists(MODELS_HOME)) MODELS_HOME <- sub("R$", "exec", getwd())

#' @importFrom methods new
make_stanfit <- function(f) {
  model_cppname <- sub("\\.stan$", "", basename(f))
  program <- c(readLines(file.path(MODELS_HOME, "functions.txt")), 
               readLines(f))
  program <- paste(program, collapse = "\n")
  stanfit <- rstan::stanc(model_code = program, model_name = model_cppname, 
                          obfuscate_model_name = FALSE)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name, 
                            model_cppcode = stanfit$cppcode)
  return(do.call(new, args = c(stanfit[-(1:3)], Class = "stanmodel", 
                 mk_cppmodule = function(x) get(paste0("model_", model_cppname)))))
}

stanfit_bernoulli <- make_stanfit(file.path(MODELS_HOME, "bernoulli.stan"))
stanfit_binomial <- make_stanfit(file.path(MODELS_HOME, "binomial.stan"))
stanfit_continuous  <- make_stanfit(file.path(MODELS_HOME, "continuous.stan"))
stanfit_count <- make_stanfit(file.path(MODELS_HOME, "count.stan"))
stanfit_lm <- make_stanfit(file.path(MODELS_HOME, "lm.stan"))
stanfit_polr <- make_stanfit(file.path(MODELS_HOME, "polr.stan"))
