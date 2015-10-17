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

stanfit_bernoulli <- rstan::stanc(file.path(MODELS_HOME, "bernoulli.stan"),
                                  obfuscate_model_name = FALSE)
stanfit_bernoulli$model_cpp <- list(model_cppname = stanfit_bernoulli$model_name, 
                                    model_cppcode = stanfit_bernoulli$cppcode)
stanfit_bernoulli <- do.call(new, args = c(stanfit_bernoulli[-(1:3)], Class = "stanmodel", 
                                           mk_cppmodule = function(x) model_bernoulli))

stanfit_binomial <- rstan::stanc(file.path(MODELS_HOME, "binomial.stan"),
                                 obfuscate_model_name = FALSE)
stanfit_binomial$model_cpp <- list(model_cppname = stanfit_binomial$model_name, 
                                   model_cppcode = stanfit_binomial$cppcode)
stanfit_binomial <- do.call(new, args = c(stanfit_binomial[-(1:3)], Class = "stanmodel", 
                                          mk_cppmodule = function(x) model_binomial))

stanfit_continuous <- rstan::stanc(file.path(MODELS_HOME, "continuous.stan"), 
                                   obfuscate_model_name = FALSE)
stanfit_continuous$model_cpp <- list(model_cppname = stanfit_continuous$model_name, 
                                     model_cppcode = stanfit_continuous$cppcode)
stanfit_continuous <- do.call(new, args = c(stanfit_continuous[-(1:3)], Class = "stanmodel", 
                                            mk_cppmodule = function(x) model_continuous))

stanfit_count <- rstan::stanc(file.path(MODELS_HOME, "count.stan"), 
                              obfuscate_model_name = FALSE)
stanfit_count$model_cpp <- list(model_cppname = stanfit_count$model_name, 
                                model_cppcode = stanfit_count$cppcode)
stanfit_count <- do.call(new, args = c(stanfit_count[-(1:3)], Class = "stanmodel", 
                                            mk_cppmodule = function(x) model_count))

stanfit_lm <- rstan::stanc(file.path(MODELS_HOME, "lm.stan"), 
                           obfuscate_model_name = FALSE)
stanfit_lm$model_cpp <- list(model_cppname = stanfit_lm$model_name, 
                             model_cppcode = stanfit_lm$cppcode)
stanfit_lm <- do.call(new, args = c(stanfit_lm[-(1:3)], Class = "stanmodel", 
                                    mk_cppmodule = function(x) model_lm))

stanfit_polr <- rstan::stanc(file.path(MODELS_HOME, "polr.stan"), 
                           obfuscate_model_name = FALSE)
stanfit_polr$model_cpp <- list(model_cppname = stanfit_polr$model_name, 
                               model_cppcode = stanfit_polr$cppcode)
stanfit_polr <- do.call(new, args = c(stanfit_polr[-(1:3)], Class = "stanmodel", 
                                    mk_cppmodule = function(x) model_polr))
