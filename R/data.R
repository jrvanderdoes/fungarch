#' SPY 10-minute closing-price curves for 2018--2020
#'
#' A data frame containing intraday SPY closing-price curves for the 2018,
#' 2019, and 2020 calendar years. Each column represents one trading day and
#' each row represents a 10-minute position during the regular US equity
#' market session. Values are the `Close` prices from the corresponding
#' 10-minute observations.
#'
#' The data retain observations from 09:30 through 15:59, remove dates whose
#' source prices are constant, and carry the most recent observed price
#' forward when a day has a missing intraday observation. Leading missing
#' values are filled with that day's first observed price.
#'
#' @format A data frame with one numeric column per trading day. Column names
#'   are dates in `MM-DD-YYYY` format.
#' @source Historical SPY data prepared from `SPY_2018.csv`, `SPY_2019.csv`,
#'   and `SPY_2020.csv` in `data-raw/`.
#' @usage data(SPY_18_20)
#' @seealso [SPY_22_24]
#' @keywords datasets
"SPY_18_20"

#' SPY 10-minute closing-price curves for 2022--2024
#'
#' A data frame containing intraday SPY closing-price curves for the 2022,
#' 2023, and 2024 calendar years. Each column represents one trading day and
#' each row represents a 10-minute position during the regular US equity
#' market session. Values are the `Close` prices from the corresponding
#' 10-minute observations.
#'
#' The data retain observations from 09:30 through 15:59, remove dates whose
#' source prices are constant, and carry the most recent observed price
#' forward when a day has a missing intraday observation. Leading missing
#' values are filled with that day's first observed price.
#'
#' @format A data frame with one numeric column per trading day. Column names
#'   are dates in `MM-DD-YYYY` format.
#' @source Historical SPY data prepared from `SPY_2022.csv`, `SPY_2023.csv`,
#'   and `SPY_2024.csv` in `data-raw/`, which were generated from the
#'   one-minute `SPY_clean_1min.csv` source.
#' @usage data(SPY_22_24)
#' @seealso [SPY_18_20]
#' @keywords datasets
"SPY_22_24"
