stopifnot(require(rstanarm))
stopifnot(require(loo))
stopifnot(require(ggplot2))
stopifnot(require(parallel))
options(mc.cores = parallel::detectCores())

SEED <- 12345L
REFRESH <- 1000L

if (R.version$major < 3 || R.version$minor < 2.0)
  warning("This demo may not work on older versions of R due to HTTPS URLs")

ROOT <- "https://raw.githubusercontent.com/stan-dev/example-models/master/"
DATA_ENV <- new.env()
