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

# This file is only intended to be used during the installation process
# nocov start
MODELS_HOME <- "exec"
if (!file.exists(MODELS_HOME)) MODELS_HOME <- sub("R$", "exec", getwd())

stan_files <- dir(MODELS_HOME, pattern = "stan$", full.names = TRUE)
stanmodels <- sapply(stan_files, function(f) {
  model_cppname <- sub("\\.stan$", "", basename(f))
  isystem <- system.file("chunks", package = methods::getPackageName(environment(), FALSE))
  if (!file.exists(file.path(isystem, "common_functions.stan")))
    isystem <- file.path("inst", "chunks")
  if (!file.exists(file.path(isystem, "common_functions.stan")))
    isystem <- file.path("..", "inst", "chunks")
  stanfit <- rstan::stanc_builder(f, isystem, allow_undefined = TRUE)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name, 
                            model_cppcode = stanfit$cppcode)
  return(do.call(methods::new, args = c(stanfit[-(1:3)], Class = "stanmodel", 
                 mk_cppmodule = function(x) get(paste0("model_", model_cppname)))))
  }
)
names(stanmodels) <- sub("\\.stan$", "", basename(names(stanmodels)))
rm(MODELS_HOME)
# nocov end
