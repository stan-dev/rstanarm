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

.onLoad <- function(libname, pkgname) {
  modules <- paste0("stan_fit4", names(stanmodels), "_mod")
  if (.Platform$OS.type == "windows" && .Platform$r_arch =="i386") 
    modules <- grep("jm", modules, invert = TRUE, value = TRUE)
  for (m in modules) loadModule(m, what = TRUE)
}

.onAttach <- function(...) {
  ver <- utils::packageVersion("rstanarm")
  packageStartupMessage("This is rstanarm version ", ver)
  packageStartupMessage("- See https://mc-stan.org/rstanarm/articles/priors for changes to default priors!")
  packageStartupMessage("- Default priors may change, so it's safest to specify priors, even if equivalent to the defaults.")
  packageStartupMessage("- For execution on a local, multicore CPU with excess RAM we recommend calling")
  packageStartupMessage("  options(mc.cores = parallel::detectCores())")
}

