#' Audio 'announcements' when printing model summaries
#' 
#' The only function of this package is to augment the printing of
#' model summaries with appropriate sound effects based on whether
#' you have achieved the magical p<0.05 level of significance for
#' at least one variable (not including the intercept - we have to
#' try to be a \emph{little} bit sensible ...
#'
#' Sounds are taken from \url{http://freesounds.org}
#'
#' \itemize{
#' \item fanfare: user primordiality, \url{https://www.freesound.org/people/primordiality/sounds/78823/} (CC BY 3.0)
#' \item trombone: user kirbydx, \url{https://www.freesound.org/people/kirbydx/sounds/175409/} (CC0)
#' }
#' @examples
#' \dontrun{
#' m1 <- lm(speed~dist,cars)
#' m2 <- lm(Income~Population,data.frame(state.x77))
#' summary(m1)
#' summary(m2)
#' }
#' 
#' @name celebrate
NULL
