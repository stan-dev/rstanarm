stan_betareg.fit <- function (x, y, z, weights = rep(1, NROW(x)), offset = rep(0, NROW(x)),
                              link = c("logit", "probit", "cloglog", "cauchit", "log", "loglog"), 
                              link.phi = c("log", "identity", "sqrt"), ...,
                              prior = normal(), prior_intercept = normal(),
                              prior_ops = prior_options(), prior_PD = FALSE, 
                              algorithm = c("sampling", "optimizing", "meanfield", "fullrank"),
                              adapt_delta = NULL, QR = FALSE, sparse = FALSE, Z_true) {
  
  # lots of tedious but simple stuff including standata which is a big list to pass to data {}
  # process the prior information like stan_glm.fit() does
  
  algorithm <- match.arg(algorithm)
  
  # no family argument
  famname <- "beta"
  is_beta <- is.beta(famname)
  
  # link for X variables
  link <- match.arg(link)
  supported_links <- c("logit", "probit", "cloglog", "cauchit", "log", "loglog")
  link_num <- which(supported_links == link)
  if (!length(link)) 
    stop("'link' must be one of ", paste(supported_links, collapse = ", "))
  
  # link for Z variables
  link.phi <- match.arg(link.phi)
  supported_phi_links <- c("log", "identity", "sqrt")
  link_num_phi <- which(supported_phi_links == link.phi)
  if (!length(link_num_phi)) 
    stop("'link' must be one of ", paste(supported_phi_links, collapse = ", "))
  
  # useless assignments to pass R CMD check
  has_intercept <- min_prior_scale <- prior_df <- prior_df_for_intercept <-
    prior_dist <- prior_dist_for_intercept <- prior_mean <- prior_mean_for_intercept <-
    prior_scale_for_dispersion <- scaled <- NULL

  x_stuff <- center_x(x, sparse)
  for (i in names(x_stuff)) # xtemp, xbar, has_intercept
    assign(i, x_stuff[[i]])
  nvars <- ncol(xtemp)
  
  # z_stuff <- center_x(z, sparse)
  # ztemp <- z_stuff$xtemp
  # zbar <- z_stuff$xbar
  # has_intercept_z <- z_stuff$has_intercept
  
  for (i in names(prior_ops)) # scaled, min_prior_dispersion, prior_scale_for_dispersion
    assign(i, prior_ops[[i]])
  
  ok_dists <- nlist("normal", student_t = "t", "cauchy", "hs", "hs_plus")
  ok_intercept_dists <- ok_dists[1:3]
  
  # prior distributions (handle_glm_prior() from data_block.R)
  prior_stuff <- handle_glm_prior(prior, nvars, link, default_scale = 2.5)
  for (i in names(prior_stuff)) # prior_{dist, mean, scale, df}
    assign(i, prior_stuff[[i]])
  prior_intercept_stuff <- handle_glm_prior(prior_intercept, nvars = 1, default_scale = 10,
                                            link, ok_dists = 
                                              nlist("normal", student_t = "t", "cauchy"))
  names(prior_intercept_stuff) <- paste0(names(prior_intercept_stuff), "_for_intercept")
  for (i in names(prior_intercept_stuff)) # prior_{dist, mean, scale, df}_for_intercept
    assign(i, prior_intercept_stuff[[i]])
  
  # create entries in the data block of the .stan file
  standata <- nlist(
    N = nrow(xtemp), K = ncol(xtemp), xbar = as.array(xbar), dense_X = !sparse, # TRUE, sparse = FALSE,
    X = array(xtemp, dim = c(1L, dim(xtemp))),
    nnz_X = 0L, 
    w_X = double(), 
    v_X = integer(), 
    u_X = integer(), 
    y = y, 
    prior_PD, has_intercept, family = 4L, link = link_num, 
    prior_dist, prior_mean, prior_scale = as.array(pmin(.Machine$double.xmax, prior_scale)), prior_df,
    prior_dist_for_intercept, prior_mean_for_intercept = c(prior_mean_for_intercept), 
    prior_scale_for_intercept = min(.Machine$double.xmax, prior_scale_for_intercept), 
    prior_df_for_intercept = c(prior_df_for_intercept),
    prior_scale_for_dispersion = prior_scale_for_dispersion %ORifINF% 0,
    has_weights = length(weights) > 0, weights = double(),
    has_offset = length(offset) > 0, offset = double(),
    t = 0L, 
    p = integer(), 
    l = integer(), 
    q = 0L, 
    len_theta_L = 0L, shape = double(), scale = double(), 
    len_concentration = 0L, concentration = double(),
    len_regularization = 0L, regularization = double(),
    num_non_zero = 0L, 
    w = double(), 
    v = integer(), 
    u = integer(),
    no_Z = Z_true,
    betareg_Z_dim = ncol(z),
    link_phi = link_num_phi,
    betareg_Z = z
    )
  print("*** LINK PHI ***")
  print(link.phi)
  print(link_num_phi)
  
  # call stan() to draw from posterior distribution
  stanfit <- stanmodels$continuous
  if (Z_true == 1) {
    pars <- c(if (has_intercept) "alpha", "beta", "omega", "mean_PPD")
  }
  else {
    pars <- c(if (has_intercept) "alpha", "beta", "dispersion", "mean_PPD")
  }
  
  if (algorithm == "optimizing") {
    out <- optimizing(stanfit, data = standata, draws = 1000, constrained = TRUE, ...)
    out$par <- out$par[!grepl("eta_Z", names(out$par))] # kinda sketch - might need fixing
    out$theta_tilde <- out$theta_tilde[,!grepl("eta_Z", colnames(out$theta_tilde))] # kinda sketch - might need fixing
    new_names <- names(out$par)
    mark <- grepl("^beta\\[[[:digit:]]+\\]$", new_names)
    new_names[mark] <- colnames(xtemp)
    new_names[new_names == "alpha[1]"] <- "(Intercept)"
    if (Z_true == 1) {
      mark_z <- grepl("^omega\\[[[:digit:]]+\\]$", new_names)
      new_names[mark_z] <- paste0("(phi)_", colnames(z))
    }
    else {
      new_names[new_names == "dispersion"] <- "(phi)"
    }
    names(out$par) <- new_names
    colnames(out$theta_tilde) <- new_names
    out$stanfit <- suppressMessages(sampling(stanfit, data = standata, chains = 0))
    return(out)
  }
  else {
    if (algorithm == "sampling") {
      sampling_args <- set_sampling_args(
        object = stanfit, 
        prior = prior, 
        user_dots = list(...), 
        user_adapt_delta = adapt_delta, 
        data = standata, 
        pars = pars, 
        show_messages = FALSE)
      stanfit <- do.call(sampling, sampling_args)
    }
    else if (algorithm == "meanfield") { # FIXME
      stanfit <- rstan::vb(stanfit, pars = pars, data = standata,
                           algorithm = algorithm, init = 0.001, ...)
    }
    else if (algorithm == "fullrank") { # FIXME
      stanfit <- rstan::vb(stanfit, pars = pars, data = standata,
                           algorithm = algorithm, init = 0.001, ...)
    }
    if (Z_true == 1) {
      new_names <- c(if (has_intercept) "(Intercept)", 
                     colnames(xtemp), paste0("(phi)_", colnames(z)),
                     "mean_PPD", "log-posterior")
    }
    else {
      new_names <- c(if (has_intercept) "(Intercept)", 
                     colnames(xtemp), 
                     "(phi)",  "mean_PPD", "log-posterior")
    }
    stanfit@sim$fnames_oi <- new_names
    return(stanfit)
  }
}
