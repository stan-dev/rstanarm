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

MODELS_HOME <- file.path(dirname(system.file(package = "rstanarm")), 
                         "rstanarm", "exec")
stanfit_lm <- rstan::stan_model(file.path(MODELS_HOME, "lm.stan"),
                                model_name = "Linear Regression",
                                auto_write = TRUE)

stanfit_gaussian <- rstan::stan_model(file.path(MODELS_HOME, "gaussian2.stan"), 
                                      model_name = "Gaussian GLM",
                                      auto_write = TRUE)

stanfit_discrete <- rstan::stan_model(file.path(MODELS_HOME, "discrete2.stan"), 
                                      model_name = "Discrete GLM",
                                      auto_write = TRUE)
