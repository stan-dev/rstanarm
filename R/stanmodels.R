# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015 Trustees of Columbia University
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

MODELS_HOME <- "exec"
if (!file.exists(MODELS_HOME)) MODELS_HOME <- sub("R$", "exec", getwd())

make_stanmodel <- function(f) { # nocov start
  model_cppname <- sub("\\.stan$", "", basename(f))
  stanfit <- rstan::stanc_builder(f)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name, 
                            model_cppcode = stanfit$cppcode)
  return(do.call(methods::new, args = c(stanfit[-(1:3)], Class = "stanmodel", 
                 mk_cppmodule = function(x) get(paste0("model_", model_cppname)))))
} # nocov end

stan_files <- dir(MODELS_HOME, pattern = "stan$", full.names = TRUE)
stanmodels <- sapply(stan_files, make_stanmodel)
names(stanmodels) <- sub("\\.stan$", "", basename(names(stanmodels)))
