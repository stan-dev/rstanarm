#' Estimation algorithms available for \pkg{rstanarm} models
#' 
#' @name available-algorithms
#' 
#' @section Estimation algorithms:
#' The modeling functions in the \pkg{rstanarm} package take an \code{algorithm}
#' argument that can be one of the following:
#' \describe{
#'  \item{\strong{Sampling} (\code{algorithm="sampling"})}{
#'  Uses Markov Chain Monte Carlo (MCMC) --- in particular, Hamiltonian Monte
#'  Carlo (HMC) with a tuned but diagonal mass matrix --- to draw from the
#'  posterior distribution of the parameters. See \code{\link[rstan:stanmodel-method-sampling]{sampling}}
#'  (\pkg{rstan}) for more details. This is the slowest but most reliable of the
#'  available estimation algorithms and it is \strong{the default and
#'  recommended algorithm for statistical inference.}
#'  }
#'  \item{\strong{Mean-field} (\code{algorithm="meanfield"})}{
#'  Uses mean-field variational inference to draw from an approximation to the
#'  posterior distribution. In particular, this algorithm finds the set of
#'  independent normal distributions in the unconstrained space that --- when
#'  transformed into the constrained space --- most closely approximate the
#'  posterior distribution. Then it draws repeatedly from these independent
#'  normal distributions and transforms them into the constrained space. The
#'  entire process is much faster than HMC and yields independent draws but
#'  \strong{is not recommended for final statistical inference}. It can be
#'  useful to narrow the set of candidate models in large problems, particularly
#'  when specifying \code{QR=TRUE} in \code{\link{stan_glm}},
#'  \code{\link{stan_glmer}}, and \code{\link{stan_gamm4}}, but is \strong{only
#'  an approximation to the posterior distribution}.
#'  }
#'  \item{\strong{Full-rank} (\code{algorithm="fullrank"})}{
#'  Uses full-rank variational inference to draw from an approximation to the
#'  posterior distribution by finding the multivariate normal distribution in
#'  the unconstrained space that --- when transformed into the constrained space
#'  --- most closely approximates the posterior distribution. Then it draws
#'  repeatedly from this multivariate normal distribution and transforms the
#'  draws into the constrained space. This process is slower than meanfield
#'  variational inference but is faster than HMC. Although still an
#'  approximation to the posterior distribution and thus \strong{not recommended
#'  for final statistical inference}, the approximation is more realistic than
#'  that of mean-field variational inference because the parameters are not
#'  assumed to be independent in the unconstrained space. Nevertheless, fullrank
#'  variational inference is a more difficult optimization problem and the
#'  algorithm is more prone to non-convergence or convergence to a local
#'  optimum.
#'  }
#'  \item{\strong{Optimizing} (\code{algorithm="optimizing"})}{
#'  Finds the posterior mode using a C++ implementation of the LBGFS algorithm.
#'  See \code{\link[rstan:stanmodel-method-optimizing]{optimizing}} for more details. If there is no prior
#'  information, then this is equivalent to maximum likelihood, in which case
#'  there is no great reason to use the functions in the \pkg{rstanarm} package
#'  over the emulated functions in other packages. However, if priors are
#'  specified, then the estimates are penalized maximum likelihood estimates,
#'  which may have some redeeming value. Currently, optimization is only
#'  supported for \code{\link{stan_glm}}.
#'  }
#' }
#' 
#' @seealso \url{https://mc-stan.org/rstanarm/}
#' 
NULL
