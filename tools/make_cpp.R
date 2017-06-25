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
options(warn = 3L)
stan_files <- dir("exec", pattern = "stan$", full.names = TRUE)
include_files <- dir("src", pattern = "hpp$")
cat(readLines(file.path("inst", "chunks", "license.stan")),
  "#ifndef MODELS_HPP", "#define MODELS_HPP",
  "#define STAN__SERVICES__COMMAND_HPP", "#include <rstan/rstaninc.hpp>",
  if (length(include_files)) paste0('#include "', include_files, '"'),
  sapply(stan_files, FUN = function(f) {
    cppcode <- rstan::stanc_builder(f, allow_undefined = TRUE,
                 isystem = file.path("inst", "chunks"))$cppcode
    cppcode <- gsub("typedef.*stan_model.*;", "", cppcode, perl = TRUE)
    return(cppcode)
  }), "#endif", file = file.path("src", "include", "models.hpp"), 
  sep = "\n", append = FALSE)

options("useFancyQuotes" = FALSE)

sapply(sub(".stan", "", basename(stan_files), fixed = TRUE), function(f) {
  Rcpp::exposeClass(class = paste0("model_", f),
                    constructors = list(c("SEXP", "SEXP", "SEXP")), fields = character(),
                    methods = c("call_sampler", 
                                "param_names", "param_names_oi", "param_fnames_oi", 
                                "param_dims",  "param_dims_oi", "update_param_oi", "param_oi_tidx", 
                                "grad_log_prob", "log_prob", 
                                "unconstrain_pars", "constrain_pars", "num_pars_unconstrained", 
                                "unconstrained_param_names", "constrained_param_names"), 
                    file = paste0(f, "Module.cc"), header = '#include "include/models.hpp"', 
                    module = paste0("stan_fit4", f, "_mod"), 
                    CppClass = paste0("rstan::stan_fit<model_", f, "_namespace::model_", f,
                                      ", boost::random::ecuyer1988> "),
                    Rfile = FALSE)
  return(invisible(NULL))
})
