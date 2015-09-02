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

# if you change a .stan file, source() stanmodels.R when the working 
# directory is the root of rstanarm/ in order to update the .rda file 
# and reduce Build & Reload time

MODELS_HOME <- "exec"
if (!file.exists(MODELS_HOME)) MODELS_HOME <- sub("R$", "exec", getwd())
  
stanfit_lm <- rstan::stan_model(file.path(MODELS_HOME, "lm.stan"),
                                model_name = "Linear Regression",
                                auto_write = interactive(), 
                                obfuscate_model_name = FALSE)
stanfit_continuous <- rstan::stan_model(file.path(MODELS_HOME, "continuous.stan"), 
                                        model_name = "Continuous GLM",
                                        auto_write = interactive(),
                                        obfuscate_model_name = FALSE)
stanfit_bernoulli <- rstan::stan_model(file.path(MODELS_HOME, "bernoulli.stan"), 
                                      model_name = "Bernoulli GLM",
                                      auto_write = interactive(),
                                      obfuscate_model_name = FALSE)
stanfit_binomial <- rstan::stan_model(file.path(MODELS_HOME, "binomial.stan"), 
                                      model_name = "Binomial GLM",
                                      auto_write = interactive(),
                                      obfuscate_model_name = FALSE)
stanfit_count <- rstan::stan_model(file.path(MODELS_HOME, "count.stan"), 
                                   model_name = "Count GLM",
                                   auto_write = interactive(),
                                   obfuscate_model_name = FALSE)
stanfit_polr <- rstan::stan_model(file.path(MODELS_HOME, "polr.stan"), 
                                   model_name = "Proportional Odds GLM",
                                   auto_write = interactive(),
                                  obfuscate_model_name = FALSE)
