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
set.seed(SEED)
CHAINS <- 2
ITER <- 100
REFRESH <- 0

SW <- suppressWarnings

plink <- function(fit, nd = NULL, sef = TRUE) 
  predict(fit, newdata = nd, type = "link", se.fit = sef)
presp <- function(fit, nd = NULL, sef = TRUE) 
  predict(fit, newdata = nd, type = "response", se.fit = sef)

context("predict")
test_that("predict ok for binomial", {
  # example from help(predict.glm)
  ldose <- rep(0:5, 2)
  numdead <- c(1, 4, 9, 13, 18, 20, 0, 2, 6, 10, 12, 16)
  sex <- factor(rep(c("M", "F"), c(6, 6)))
  SF <- cbind(numdead, numalive = 20-numdead)
  
  glmfit <- glm(SF ~ sex*ldose, family = binomial)
  stanfit <- SW(stan_glm(SF ~ sex*ldose, family = binomial, chains = CHAINS, 
                         iter = ITER, seed = SEED, refresh = REFRESH))
  stanfit_opt <- SW(update(stanfit, algorithm = "optimizing"))
  
  
  pg <- plink(glmfit)
  ps <- plink(stanfit)
  pso <- plink(stanfit_opt)
  expect_equal(pg$fit, ps$fit, tol = 0.1)
  expect_equal(pg$fit, pso$fit, tol = 0.05)
  expect_equal(pg$se.fit, ps$se.fit, tol = 0.1)
  expect_equal(pg$se.fit, pso$se.fit, tol = 0.1)
  expect_equal(presp(glmfit)[1:2], presp(stanfit_opt), tol = 0.05)
  expect_error(presp(stanfit))
  
  ld <- seq(0, 5, 0.1)
  newd <- data.frame(ldose = ld, sex = factor(rep("M", length(ld)), 
                                              levels = levels(sex)))
  pg <- plink(glmfit, newd)
  ps <- plink(stanfit, newd)
  pso <- plink(stanfit_opt, newd)
  expect_equal(pg$fit, ps$fit, tol = 0.05)
  expect_equal(pg$fit, pso$fit, tol = 0.05)
  expect_equal(pg$se.fit, ps$se.fit, tol = 0.1)
  expect_equal(pg$se.fit, pso$se.fit, tol = 0.1)
  expect_equal(presp(glmfit, newd)[1:2], presp(stanfit_opt, newd), tol = 0.1)
})

test_that("predict ok for gaussian", {
  glmfit <- glm(mpg ~ wt, data = mtcars)
  stanfit <- SW(stan_glm(mpg ~ wt, data = mtcars, chains = CHAINS,
                      iter = 2 * ITER, seed = SEED, refresh = REFRESH))
  stanfit_opt <- SW(update(stanfit, algorithm = "optimizing"))
  
  pg <- plink(glmfit)
  ps <- plink(stanfit)
  pso <- plink(stanfit_opt)
  expect_equal(pg$fit, ps$fit, tol = 0.05)
  expect_equal(pg$fit, pso$fit, tol = 0.05)
  expect_equal(pg$se.fit, ps$se.fit, tol = 0.1)
  expect_equal(pg$se.fit, pso$se.fit, tol = 0.1)
  expect_equal(presp(glmfit)[1:2], presp(stanfit_opt), tol = 0.1)
  expect_error(presp(stanfit))

  newd <- data.frame(wt = c(1,5))
  pg <- plink(glmfit, newd)
  ps <- plink(stanfit, newd)
  pso <- plink(stanfit_opt, newd)
  expect_equal(pg$fit, ps$fit, tol = 0.05)
  expect_equal(pg$fit, pso$fit, tol = 0.05)
  expect_equal(pg$se.fit, ps$se.fit, tol = 0.1)
  expect_equal(pg$se.fit, pso$se.fit, tol = 0.1)
  expect_equal(presp(glmfit, newd)[1:2], presp(stanfit_opt, newd), tol = 0.1)
})

test_that("predict ok for Poisson", {
  dat <- data.frame(counts = c(18,17,15,20,10,20,25,13,12),
                    outcome = gl(3,1,9), treatment = gl(3,3))

  glmfit <- glm(counts ~ outcome + treatment, data = dat, family = poisson())
  stanfit <- SW(stan_glm(counts ~ outcome + treatment, data = dat, family = poisson(), 
                         chains = CHAINS, iter = ITER, seed = SEED, refresh = REFRESH))
  stanfit_opt <- SW(update(stanfit, algorithm = "optimizing"))

  pg <- plink(glmfit)
  ps <- plink(stanfit)
  pso <- plink(stanfit_opt)
  expect_equal(pg$fit, ps$fit, tol = 0.05)
  expect_equal(pg$fit, pso$fit, tol = 0.05)
  expect_equal(pg$se.fit, ps$se.fit, tol = 0.1)
  expect_equal(pg$se.fit, pso$se.fit, tol = 0.1)
  expect_equal(presp(glmfit)[1:2], presp(stanfit_opt), tol = 0.1)
  expect_error(presp(stanfit))

  expect_equal(plink(stanfit, sef = FALSE), plink(glmfit, sef = FALSE), tol = 0.05)
  expect_equal(presp(stanfit, sef = FALSE), presp(glmfit, sef = FALSE), tol = 0.05)

  newd <- dat[1:2, ]
  pg <- plink(glmfit, newd)
  ps <- plink(stanfit, newd)
  pso <- plink(stanfit_opt, newd)
  expect_equal(pg$fit, ps$fit, tol = 0.05)
  expect_equal(pg$fit, pso$fit, tol = 0.05)
  expect_equal(pg$se.fit, ps$se.fit, tol = 0.1)
  expect_equal(pg$se.fit, pso$se.fit, tol = 0.1)
  expect_equal(presp(glmfit, newd)[1:2], presp(stanfit_opt, newd), tol = 0.1)
})
