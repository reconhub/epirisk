#' Compute infection risks from past cases
#'
#' This function computes the probability of secondary infections per infected
#' cases, based on their dates of onset and isolation, and on the average
#' reproduction number (i.e. number of secondary cases per index case). Only
#' cases isolated within a defined time window are retained for the
#' analysis. NOTE: this function is experimental and developed for the response
#' to Ebola in North Kivu. Do not use it without consulting the authors first.
#'
#' @details
#'
#' If $d(.)$ is the probability mass function of the delay from onset
#' of the index case to secondary infection, and $R$ the average reproduction
#' number, the average risk of secondary infections is computed as:
#'
#' $$ risk(t) = R * \sum_{t = onset}^{isolation - 1} w(t) $$
#' 
#' @export
#'
#' @author Thibaut Jombart (thibautjombart@@gmail.com)
#'
#' @param onset: a vector of dates of onset
#'
#' @param isolation: a vector of dates of isolation / outcome (when patients are
#'   assumed to ' no longer seed infections)
#'
#' @param R the average reproduction number, i.e. number of secondary cases per
#'   infected case
#'
#' @param time_period the number of days since the isolation for cases to be
#'   included in the analysis
#'
#' @param reference_day the day to use as current date, as a `Date` object;
#'   defaults to the current day
#'
#' @param dist_delay distribution of the delay from symptom onset to secondary
#'   infections, as a `distcrete` object (see the RECON package of the same name
#'   for more information)
#' 

case_risk <- function(onset, isolation, R, time_period, dist_delay,
                      reference_day = today()) {
  
  ## check that inputs are correct
  if (!inherits(onset, "Date")) {
    msg <- sprintf("onset is not a `Date` object but a %s",
                   class(onset)[1])
    stop(msg)
  }
  if (!inherits(isolation, "Date")) {
    msg <- sprintf("isolation is not a `Date` object but a %s",
                   class(isolation)[1])
    stop(msg)
  }
  if (any(na.omit(onset > isolation))) {
    msg <- "Some dates of isolation predate the onset"
    stop(msg)
  }
  if (any(is.na(onset))) {
    msg <- "Some onset dates are missing: these cases will be ignored"
    warning(msg)
  }
  if (any(is.na(isolation))) {
    msg <- "Some isolation dates are missing: these cases will be ignored"
    warning(msg)
  }
  if (is.null(R) || !is.finite(R)) {
    msg <- "`R` needs to be a finite number"
    stop(msg)
  }
  if (R < 0) {
    msg <- "`R` needs to be a positive number"
    stop(msg)
  }
  if (!is.finite(time_period)) {
    msg <- "`time_period` needs to be a finite number"
    stop(msg)
  }
  if (time_period < 1) {
    msg <- "`time_period` needs to be greater than 0"
    stop(msg)
  }
  if (!inherits(reference_day, "Date")) {
    msg <- "`reference_day` needs to be a `Date` object"
  }

  if (!inherits(dist_delay, "distcrete")) {
    msg <- sprintf("`dist_delay` is not a `distcrete` object, but a %s",
                   class(dist_delay)[1])
    stop(msg)
  }

  ## extract the cumulative mass function
  p_delay <- dist_delay$p

  ## select cases to retain
  stop_at <- reference_day - time_period
  to_keep <- !is.na(isolation) & (isolation >= stop_at)
  if (!any(to_keep)) {
    msg <- sprintf("No cases have an isolation date in the `time_period` %s - %s",
                   format(stop_at, "%d %b %Y"),
                   format(reference_day, "%d %b %Y"))
    warning("msg")
  }


  ## compute delays
  delays <- as.integer(isolation - onset)

  ## get average number of secondary cases
  out <- p_delay(delays) * R

  ## set risks of cases outside the time window to zero
  out[to_keep] <- 0
  out

}
