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

stanlm <-
  function (formula, data, subset, weights, na.action, method = "qr", 
            model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE, 
            contrasts = NULL, offset = NULL,
            # above arguments from glm(), below arguments from arm:::bayesglm()            
            prior.mean = 0, prior.scale = NULL, prior.df = 1,
            prior.mean.for.intercept = 0, prior.scale.for.intercept = NULL,
            prior.df.for.intercept = 1, min.prior.scale = 1e-12, scaled = TRUE,
            prior.scale.for.dispersion = 5, ...) { # further arguments to stan()
    
    mf <- match.call(expand.dots = TRUE)
    mf[[1L]] <- as.name("stanglm")
    mf$qr <- mf$singular.ok <- NULL
    return(eval(mf, parent.frame()))
    
}


