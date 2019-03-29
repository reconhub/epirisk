#' Compute infection risks from past cases
#'
#' This function computes the probability of secondary infections per infected
#' cases, based on their dates of onset and isolation, and on the average
#' reproduction number (i.e. number of secondary cases per index case). Only
#' cases isolated within a defined time window are retained for the
#' analysis. NOTE: this function is experimental and developed for the response
#' to Ebola in North Kivu. Do not use it without consulting the authors first.
#'
#' 
#' @details
#'
#' If \eqn{d(.)} is the probability mass function of the delay from onset
#' of the index case to secondary infection, and $R$ the average reproduction
#' number, the average risk of secondary infections is computed as:
#'
#' \deqn{ risk(t) = R * \sum_{t = onset}^{isolation - 1} w(t) }
#' 
#' 
#' @export
#' 
#'
#' @author Thibaut Jombart (thibautjombart@@gmail.com)
#'
#' @param onset a vector of dates of onset
#'
#' @param isolation a vector of dates of isolation / outcome (when patients are
#'   assumed to ' no longer seed infections)
#'
#' @param R the average reproduction number, i.e. number of secondary cases per
#'   infected case; can be a single value, or a vector with one value per case
#'   (recycled if needed)
#'
#' @param time_period the number of days since the isolation for cases to be
#'   included in the analysis
#'
#' @param reference_day the day to use as current date, as a `Date` object;
#'   defaults to the current day
#'
#' @param p_delay function giving the cumulative mass function of the delay
#'   from symptom onset to secondary infections (in usual distributions,
#'   starting with `p` e.g. `pgeom`)
#'
#' 
#' @examples
#'
#' isolation <- today()  -  rpois(50, 10)
#' onset <- isolation - rpois(50, 7)
#' delay <- function(t) pgeom(t, prob = .2)
#' risks <- case_risk(onset, isolation, R = 1.2,
#'                    time_period = 14, p_delay = delay)
#' risks

#' df <- data.frame(id = factor(1:50), onset, isolation, risks)

#' if (require(ggplot2)) {
#'   ggplot(df, aes(y = id, yend = id)) +
#'     geom_segment(aes(x = onset, xend = isolation, color = risks),
#'                  size = 4) +
#'     geom_vline(aes(xintercept = today())) +
#'     scale_color_gradientn("Risk 2nd cases", colors = c("#80aaff", "gold", "#b3003b")) +
#'   labs(title = "Risk and infectious period",
#'        x = "Infectious period (onset -> isolation)")
#' }
#' 

case_risk <- function(onset, isolation, R, time_period, p_delay,
                      reference_day = today()) {
  
  ## check that inputs are correct
  if (!inherits(onset, "Date")) {
    msg <- sprintf("onset is not a `Date` object but a `%s`",
                   class(onset)[1])
    stop(msg)
  }
  if (!inherits(isolation, "Date")) {
    msg <- sprintf("isolation is not a `Date` object but a `%s`",
                   class(isolation)[1])
    stop(msg)
  }
  if (any(stats::na.omit(onset > isolation))) {
    msg <- "Some dates of isolation predate the onset"
    warning(msg)
  }
  if (any(is.na(onset))) {
    msg <- "Some onset dates are missing: these cases will be ignored"
    warning(msg)
  }
  if (any(is.na(isolation))) {
    msg <- "Some isolation dates are missing: these cases will be ignored"
    warning(msg)
  }
  if (is.null(R) || any(!is.finite(R))) {
    msg <- "`R` needs to be a finite number"
    stop(msg)
  }
  if (any(R < 0)) {
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

  if (!inherits(p_delay, "function")) {
    msg <- sprintf("`dist_delay` is not a `function`, but a `%s`",
                   class(p_delay)[1])
    stop(msg)
  }


  ## select cases to retain
  stop_at <- reference_day - time_period
  to_keep <- !is.na(isolation) & (isolation >= stop_at)
  if (!any(to_keep)) {
    msg <- sprintf("No cases have an isolation date in the `time_period` %s - %s",
                   format(stop_at, "%d %b %Y"),
                   format(reference_day, "%d %b %Y"))
    warning(msg)
  }


  ## compute delays
  delays <- as.integer(isolation - onset)

  ## get average number of secondary cases
  out <- p_delay(delays) * R

  ## set risks of cases outside the time window to zero
  out[!to_keep] <- 0
  out

}
