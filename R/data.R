#' Global Monthly and Annual Temperature Anomalies (degrees C), 1850-2015
#' (Relative to the 1961-1990 Mean) (May 2016)
#'
#' @name cdiac
#' @docType data
#' @source \url{https://ess-dive.lbl.gov/}
#' @references \url{https://ess-dive.lbl.gov/}
#' @keywords data
#' @format A data frame with 166 rows and 14 variables:
#' \describe{
#'   \item{year}{Year}
#'   \item{jan}{Anomaly for month of January}
#'   \item{feb}{Anomaly for month of February}
#'   \item{mar}{Anomaly for month of March}
#'   \item{apr}{Anomaly for month of April}
#'   \item{may}{Anomaly for month of May}
#'   \item{jun}{Anomaly for month of June}
#'   \item{jul}{Anomaly for month of July}
#'   \item{aug}{Anomaly for month of August}
#'   \item{sep}{Anomaly for month of September}
#'   \item{oct}{Anomaly for month of October}
#'   \item{nov}{Anomaly for month of November}
#'   \item{dec}{Anomaly for month of December}
#'   \item{annual}{Annual anomaly for the year}
#' }
"cdiac"

#' Direct Standardization: Population
#'
#' Randomly generated data for direct standardization example.
#' Sex was drawn from a Bernoulli distribution, and age was drawn from a uniform distribution on {10,...,60}.
#' The response was drawn from a normal distribution with a mean that depends on sex and age, and a variance of 1.
#'
#' @name dspop
#' @docType data
#' @keywords data
#' @seealso \link[cvxr]{dssamp}
#' @format A data frame with 1000 rows and 3 variables:
#' \describe{
#'  \item{y}{Response variable}
#'  \item{sex}{Sex of individual, coded male (0) and female (1)}
#'  \item{age}{Age of individual}
#' }
"dspop"

#' Direct Standardization: Sample
#'
#' A sample of \code{\link{dspop}} for direct standardization example.
#' The sample is skewed such that young males are overrepresented in comparison to the population.
#'
#' @name dssamp
#' @docType data
#' @keywords data
#' @seealso \link[cvxr]{dspop}
#' @format A data frame with 100 rows and 3 variables:
#' \describe{
#'  \item{y}{Response variable}
#'  \item{sex}{Sex of individual, coded male (0) and female (1)}
#'  \item{age}{Age of individual}
#' }
"dssamp"
