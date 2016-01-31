stan_mnp <- function(formula, data = parent.frame(), choiceX = NULL, cXnames = NULL, 
                     base = NULL, ...) {
  
  if (!requireNamespace("MNP"))
    stop("the 'MNP' package must be installed to use 'stan_mnp'")
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "choiceX", "cXnames", "base"), 
             names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf[[1L]] <- quote(MNP::mnp)
  mf$n.draws <- 0L
  mf <- eval(mf, parent.frame())
  
}

stan_mnp.fit <- function(y, X, ...) {
  stopifnot(is.matrix(X), is.vector(y), 
            all(y == as.integer(y)), length(y) == nrow(X))
}

debug(stan_mnp)
res1 <- stan_mnp(choice ~ 1, choiceX = list(Surf=SurfPrice, Tide=TidePrice,
                                       Wisk=WiskPrice, EraPlus=EraPlusPrice,
                                       Solo=SoloPrice, All=AllPrice),
                 data = detergent, cXnames = "price")