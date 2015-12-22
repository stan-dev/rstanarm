# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015 Trustees of Columbia University
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

# tests can be run using devtools::test() or manually by loading testthat 
# package and then running the code below possibly with options(mc.cores = 4).

library(rstanarm)
SEED <- 123

threshold <- 0.03

context("stan_polr")
test_that("stan_polr returns expected result for esoph example", {
  library(MASS)
  f <- tobgp ~ agegp + alcgp
  fit <- stan_polr(f, data = esoph, prior = R2(location = 0.4, what = "median"),
                   chains = 2, iter = 400, seed = SEED)
  fit <- stan_polr(factor(tobgp == "30+") ~ agegp + alcgp, data = esoph, 
                   prior = R2(location = 0.4), method = "loglog",
                   chains = 2, iter = 400, seed = SEED)
  
  # fit <- stan_polr(f, data = esoph, prior = NULL, 
  #                  algorithm = "fullrank", seed = SEED)
  # check <- polr(f, data = esoph)
  # expect_equal(coef(fit), coef(check), threshold)
})
