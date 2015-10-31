# This file is part of rstanarm.
# Copyright 2015 Jonah Gabry and Benjamin Goodrich
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

MODELS_HOME <- "exec"
if (!file.exists(MODELS_HOME)) MODELS_HOME <- sub("R$", "exec", getwd())

make_stanmodel <- function(f) {
  model_cppname <- sub("\\.stan$", "", basename(f))
  stanfit <- rstan::stanc_builder(f)
  stanfit$model_cpp <- list(model_cppname = stanfit$model_name, 
                            model_cppcode = stanfit$cppcode)
  return(do.call(methods::new, args = c(stanfit[-(1:3)], Class = "stanmodel", 
                 mk_cppmodule = function(x) get(paste0("model_", model_cppname)))))
}

stan_files <- dir(MODELS_HOME, pattern = "stan$", full.names = TRUE)
stanmodels <- sapply(stan_files, make_stanmodel)
names(stanmodels) <- sub("\\.stan$", "", basename(names(stanmodels)))
