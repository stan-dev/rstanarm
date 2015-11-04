# This file is part of rstanarm.
# Copyright 2013 Stan Development Team
# rstanarm is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# rstanarm is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with rstanarm.  If not, see <http://www.gnu.org/licenses/>.

.onLoad <- function(libname, pkgname) {
  loadRcppModules()
}

.onAttach <- function(...) {
  rstanarmLib <- dirname(system.file(package = "rstanarm"))
  pkgdesc <- utils::packageDescription("rstanarm", lib.loc = rstanarmLib)
  builddate <- gsub(';.*$', '', pkgdesc$Packaged)
  packageStartupMessage(paste("rstanarm (Version ", pkgdesc$Version, ", packaged: ", builddate, ")", sep = ""))
  packageStartupMessage("- Do not expect the default priors to remain the same in future rstanarm versions.")
  packageStartupMessage("Thus, R scripts should specify priors explicitly, even if they are just the defaults.")
  packageStartupMessage("- For execution on a local, multicore CPU with excess RAM we recommend calling")
  packageStartupMessage("options(mc.cores = parallel::detectCores())")
}

