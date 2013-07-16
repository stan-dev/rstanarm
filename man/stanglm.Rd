\name{stanglm}
\alias{stanglm}
\alias{stanglm.fit}
\alias{stanlm}
\title{
Fitting Bayesian generalized linear models via Stan
}
\description{
Full Bayesian inference for generalized linear modeling with student t prior distributions for the coefficients.
}
\usage{
stanglm(formula, family = gaussian(), data, weights, subset, na.action, start = NULL, 
        etastart, mustart, offset, control = list(), model = TRUE, method = "glm.fit", 
        x = FALSE, y = TRUE, contrasts = NULL, prior.mean = 0, prior.scale = NULL, 
        prior.df = 1, prior.mean.for.intercept = 0, prior.scale.for.intercept = NULL, 
        prior.df.for.intercept = 1, min.prior.scale = 1e-12, scaled = TRUE, 
        prior.scale.for.dispersion = 5, ...)

stanglm.fit(x, y, weights = rep(1, NROW(x)), start = NULL, etastart = NULL, 
           mustart = NULL, offset = rep(0, NROW(x)), family = gaussian(), 
           control = list(), intercept = TRUE,
           prior.mean = 0, prior.scale = NULL, prior.df = 1,
           prior.mean.for.intercept = 0, prior.scale.for.intercept = NULL,
           prior.df.for.intercept = 1, min.prior.scale = 1e-12, scaled = TRUE,
           prior.scale.for.dispersion = 5, ...)

stanlm(formula, data, subset, weights, na.action, method = "qr", 
       model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
       contrasts = NULL, offset = NULL,
       prior.mean = 0, prior.scale = NULL, prior.df = 1,
       prior.mean.for.intercept = 0, prior.scale.for.intercept = NULL,
       prior.df.for.intercept = 1, min.prior.scale = 1e-12, scaled = TRUE,
       prior.scale.for.dispersion = 5, ...)
}
\arguments{
  \item{formula}{
  Same as in \code{\link[stats]{glm}}
}
  \item{family}{
  Same as in \code{\link[stats]{glm}}
}
  \item{data}{
  Same as in \code{\link[stats]{glm}}
}
  \item{weights}{
  Same as in \code{\link[stats]{glm}}
}
  \item{subset}{
  Same as in \code{\link[stats]{glm}}
}
  \item{na.action}{
  Same as in \code{\link[stats]{glm}}
}
  \item{start}{
  Same as in \code{\link[stats]{glm}}. Also, if not \code{NULL} used as starting
  values for the MCMC. If \code{NULL} (the default), then \code{\link[rstan]{stan}}
  is initialized with \code{init = 'random'}.
}
  \item{etastart}{
  Same as in \code{\link[stats]{glm}} and ignored by \code{\link[rstan]{stan}}.
}
  \item{mustart}{
  Same as in \code{\link[stats]{glm}} and ignored by \code{\link[rstan]{stan}}.
}
  \item{offset}{
  Same as in \code{\link[stats]{glm}}
}
  \item{control}{
  Included for compatibility with \code{\link[stats]{glm}} but ignored by \code{\link[rstan]{stan}}.
}
  \item{intercept}{
  Included for compatibility with \code{\link[stats]{glm.fit}} but ignored by \code{\link[rstan]{stan}}.
}
  \item{model}{
  Included for compatibility with \code{\link[stats]{glm}} but ignored by \code{\link[rstan]{stan}}.
}
  \item{method}{
  Included for compatibility with \code{\link[stats]{glm}} but ignored by \code{\link[rstan]{stan}}.
}
  \item{x}{
  For \code{stanglm} and \code{stanlm}, this argument it included for compatibility with 
  \code{\link[stats]{glm}} but ignored by \code{\link[rstan]{stan}}. For \code{stanglm.fit},
  a \code{\link{model.matrix}} like that for \code{\link{glm.fit}}.
}
  \item{y}{
  For \code{stanglm} and \code{stanlm}, this argument it included for compatibility with 
  \code{\link[stats]{glm}} but ignored by \code{\link[rstan]{stan}}. For \code{stanglm.fit},
  a \code{\link{model.response}} like that for \code{\link{glm.fit}}.
}
  \item{qr}{
  Included for compatibility with \code{\link[stats]{lm}} but ignored by \code{\link[rstan]{stan}}.  
}
  \item{singular.ok}{
  Included for compatibility with \code{\link[stats]{lm}} but ignored by \code{\link[rstan]{stan}}.  
}  
  \item{contrasts}{
  Same as in \code{\link[stats]{glm}}.
}
  \item{prior.mean}{
  Prior mean (technically location) vector for the coefficient(s) in the model, not counting the intercept,
  which defaults to 0 and is replicated to the appropriate length.
}
  \item{prior.scale}{
  Prior scale vector for the coefficient(s) in the model, not counting the intercept, whose default
  depend on the nature of the model. See the Details section.
}
  \item{prior.df}{
  Prior degrees of freedom vector for the coefficient(s) in the model, not countingthe intercept, which
  defaults to 1 (implying a Cauchy prior) and is replicated to the appropriate length
}
  \item{prior.mean.for.intercept}{
  Prior mean (technically location) for the intercept, which defaults to 0.
}
  \item{prior.scale.for.intercept}{
  Prior scale for the intercept, which defaults to 10.
}
  \item{prior.df.for.intercept}{
  Prior degrees of freedom for the intercept, which defaults to 1 and implies a Cauchy prior.
}
  \item{min.prior.scale}{
  Minimum prior scale for the intercept and coefficients. See the Details section.
}
  \item{scaled}{
  Logical scalar, defaulting to \code{TRUE}, and if \code{TRUE} further scales the prior.scale
  by the range of the predictor if the predictor has exactly two unique values and scales
  prior.scale by twice the standard deviation of the predictor if it has more than two unique values.
}
  \item{prior.scale.for.dispersion}{
  Prior scale for the standard error of the regression in Gaussian models, which is given a 
  half-Cauchy prior truncated at zero.
}
  \item{\dots}{
  Further arguments passed to \code{\link[rstan]{stan}}.
}
}
\details{
The \code{stanglm} function is similar in syntax to \code{\link[stats]{glm}} but rather than
performing maximum likelihood estimation of generalized linear models, full Bayesian
estimation is performed via Markov Chain Monte Carlo. Thus, the Bayesian model
adds independent student t priors on the coefficients of the generalized linear
model. The \code{stanlm} function calls \code{stanglm} with \code{family = gaussian}.
}
\value{
An object of \code{\link[rstan]{stanfit-class}}.
}
\author{
Ben Goodrich
}
\seealso{
\code{\link[stats]{glm}}, \code{\link[rstan]{stan}}
}
\examples{
fit <- stanglm(RTEN ~ INTG + ORAL + WRIT, data = USJudgeRatings)
fit
}
