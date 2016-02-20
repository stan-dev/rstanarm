# Part of the rstanarm package for estimating model parameters
# Copyright (C) 2015, 2016 Trustees of Columbia University
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
#
#' Datasets for rstanarm examples
#' 
#' Datasets for use in \pkg{rstanarm} examples and vignettes.
#'
#' @name rstanarm-datasets
#' @keywords internal
#' @aliases kidiq roaches wells bball1970 bball2006 mortality tumors
#' @format 
#' \describe{
#' \item{\code{bball1970}}{
#' Data on hits and at-bats from the 1970 Major League Baseball season for 18
#' players.
#' 
#' Source: Efron and Morris (1975).
#' 
#' 18 obs. of 5 variables
#' \itemize{
#' \item \code{Player} Player's last name
#' \item \code{Hits} Number of hits in the first 45 at-bats of the season
#' \item \code{AB} Number of at-bats (45 for all players)
#' \item \code{RemainingAB} Number of remaining at-bats (different for most players)
#' \item \code{RemainingHits} Number of remaining hits
#' }
#' }
#' \item{\code{bball2006}}{
#' Hits and at-bats for the entire 2006 American League season of Major League
#' Baseball.
#' 
#' Source: Carpenter (2009)
#' 
#' 302 obs. of 2 variables
#' \itemize{
#' \item \code{y} Number of hits
#' \item \code{K} Number of at-bats
#' }
#' }
#' \item{\code{kidiq}}{
#' Data from a survey of adult American women and their children 
#' (a subsample from the National Longitudinal Survey of Youth).
#' 
#' Source: Gelman and Hill (2007)
#' 
#' 434 obs. of 4 variables
#' \itemize{
#' \item \code{kid_score} Child's IQ score
#' \item \code{mom_hs} Indicator for whether the mother has a high school degree
#' \item \code{mom_iq} Mother's IQ score
#' \item \code{mom_age} Mother's age
#' }
#' }
#' \item{\code{mortality}}{
#' Surgical mortality rates in 12 hospitals performing cardiac surgery
#' in babies.
#' 
#' Source: Spiegelhalter et al. (1996).
#'
#' 12 obs. of 2 variables
#' \itemize{
#' \item \code{y} Number of deaths
#' \item \code{K} Number of surgeries
#' }
#' }
#' \item{\code{roaches}}{
#' Data on the efficacy of a pest management system at reducing the number of
#' roaches in urban apartments.
#' 
#' Source: Gelman and Hill (2007)
#' 
#' 262 obs. of 6 variables
#' \itemize{
#' \item \code{y} Number of roaches caught
#' \item \code{roach1} Pretreatment number of roaches
#' \item \code{treatment} Treatment indicator
#' \item \code{senior} Indicator for only eldery residents in building
#' \item \code{exposure2} Number of days for which the roach traps were used
#' }
#' }
#' \item{\code{tumors}}{
#' Tarone (1982) provides a data set of tumor incidence in historical
#' control groups of rats; specifically endometrial stromal polyps in
#' female lab rats of type F344.  
#' 
#' Source: Gelman and Hill (2007)
#' 
#' 71 obs. of 2 variables
#' \itemize{
#' \item \code{y} Number of rats with tumors
#' \item \code{K} Number of rats
#' }
#' }
#' \item{\code{wells}}{
#' A survey of 3200 residents in a small area of Bangladesh suffering from
#' arsenic contamination of groundwater. Respondents with elevated arsenic
#' levels in their wells had been encouraged to switch their water source to a
#' safe public or private well in the nearby area and the survey was conducted
#' several years later to learn which of the affected residents had switched
#' wells.
#' 
#' Souce: Gelman and Hill (2007)
#' 
#' 3020 obs. of 5 variables
#' \itemize{
#' \item \code{switch} Indicator for well-switching
#' \item \code{arsenic} Arsenic level in respondent's well
#' \item \code{dist} Distance (meters) from the respondent's house to the
#' nearest well with safe drinking water.
#' \item \code{association} Indicator for member(s) of household participate
#' in community organizations
#' \item \code{educ} Years of education (head of household)
#' }
#' }
#' }
#' 
#' @references 
#' Carpenter, B. (2009) Bayesian estimators for the beta-binomial model of
#' batting ability. \url{http://lingpipe-blog.com/2009/09/23/}
#' 
#' Efron, B. and Morris, C. (1975) Data analysis using Stein's estimator and its
#' generalizations. \emph{Journal of the American Statistical Association}
#' \strong{70}(350), 311--319.
#' 
#' @templateVar armRef \url{http://stat.columbia.edu/~gelman/arm/}
#' @template reference-gelman-hill
#' 
#' @references
#' Spiegelhalter, D., Thomas, A., Best, N., & Gilks, W. (1996) BUGS 0.5 
#' Examples. MRC Biostatistics Unit, Institute of Public health, Cambridge, UK.
#' 
#' Tarone, R. E. (1982) The use of historical control information in testing for
#' a trend in proportions. \emph{Biometrics} \strong{38}(1):215--220.
#' 
#' @examples 
#' fit <- stan_lm(kid_score ~ mom_hs * mom_iq, data = kidiq, 
#'                prior = R2(location = 0.30, what = "mean"),
#'                # the next line is only to make the example go fast enough
#'                chains = 2, iter = 500, seed = 12345)
#' pp_check(fit, nreps = 20)
#' 
NULL
