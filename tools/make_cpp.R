stan_files <- dir("exec", pattern = "stan$", full.names = TRUE)
cat("#ifndef MODELS_HPP", "#define MODELS_HPP",
  "#define STAN__SERVICES__COMMAND_HPP", "#include <rstan/rstaninc.hpp>",
  sapply(stan_files, FUN = function(f) {
    cppcode <- rstan::stanc_builder(f)$cppcode
    cppcode <- gsub("typedef.*stan_model.*;", "", cppcode, perl = TRUE)
    return(cppcode)
  }), "#endif", file = file.path("src", "include", "models.hpp"), sep = "\n", append = FALSE)

options("useFancyQuotes" = FALSE)

sapply(sub(".stan", "", basename(stan_files), fixed = TRUE), function(f) {
  Rcpp::exposeClass(class = paste0("model_", f),
                    constructors = list(c("SEXP", "SEXP")), fields = character(),
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
