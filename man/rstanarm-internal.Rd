\name{internal-rstanarm}
\title{Internal objects of the package}
\description{
There are not any but this entry is used to silently generate the other 
documentation pages.
\Sexpr[results=hide,stage=build]{roxygen2::roxygenize()}
\Sexpr[results=hide,stage=build]{
  MODELS_HOME <- "exec"
  stanfit_lm <- rstan::stan_model(file.path(MODELS_HOME, "lm.stan"),
                                  model_name = "Linear Regression")

  stanfit_gaussian <- rstan::stan_model(file.path(MODELS_HOME, "gaussian2.stan"), 
                                        model_name = "Gaussian GLM")

  stanfit_discrete <- rstan::stan_model(file.path(MODELS_HOME, "discrete2.stan"), 
                                        model_name = "Discrete GLM")
  save(stanfit_lm, stanfit_gaussian, stanfit_discrete, 
       file = file.path("R", "sysdata.rda"), compress = TRUE)
}
}
\keyword{internal}
