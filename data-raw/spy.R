library(tidyr)

# Prepare the 10-minute SPY close curves for 2018--2020.
# Run this script from the fungarch package root.

input_files <- file.path("data-raw", paste0("SPY_", 2018:2020, ".csv"))
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop("Missing input file(s): ", paste(missing_files, collapse = ", "))
}

SPY2018 <- read.csv(input_files[1], stringsAsFactors = FALSE)
SPY2019 <- read.csv(input_files[2], stringsAsFactors = FALSE)
SPY2020 <- read.csv(input_files[3], stringsAsFactors = FALSE)

spy <- rbind(SPY2018, SPY2019, SPY2020)
spy$time <- as.POSIXct(
  spy$Local.time,
  format = "%d.%m.%Y %H:%M:%OS",
  tz = "UTC"
)

if (anyNA(spy$time)) {
  stop("Some Local.time values could not be parsed.")
}

# Keep regular US market hours, including the 09:30 bar and excluding bars
# after 15:59. The package uses the Close column as the functional curve.
data <- spy[, c("time", "Close")]
hours <- as.numeric(format(data$time, "%H"))
minutes <- as.numeric(format(data$time, "%M"))
keep_market_hours <- hours >= 9 & hours <= 15 &
  !(hours == 9 & minutes < 30)
data <- data[keep_market_hours, , drop = FALSE]

# Remove dates that contain only one unique close value, which are typically
# non-trading/flat days in the historical source data.
date_key <- format(data$time, "%Y-%m-%d")
unique_close_count <- tapply(data$Close, date_key, function(x) length(unique(x)))
drop_dates <- names(unique_close_count)[unique_close_count == 1]
data <- data[!(format(data$time, "%Y-%m-%d") %in% drop_dates), , drop = FALSE]

data$clock <- format(data$time, "%H %M")
data$day <- format(data$time, "%m-%d-%Y")

data_use <- tidyr::pivot_wider(
  data[, c("Close", "clock", "day")],
  names_from = day,
  values_from = Close
)
SPY_18_20 <- as.data.frame(data_use[, -1, drop = FALSE])

# Carry the most recent available value forward within each trading-day
# column. For a leading gap, use that day's first observed close.
for (i in seq_len(ncol(SPY_18_20))) {
  missing_idx <- which(is.na(SPY_18_20[[i]]))
  if (!length(missing_idx)) next

  nonmissing <- which(!is.na(SPY_18_20[[i]]))
  if (!length(nonmissing)) {
    stop("Column ", names(SPY_18_20)[i], " contains no observed closes.")
  }

  first_value <- SPY_18_20[[i]][nonmissing[1]]
  for (j in missing_idx) {
    if (j == 1) {
      SPY_18_20[[i]][j] <- first_value
    } else {
      SPY_18_20[[i]][j] <- SPY_18_20[[i]][j - 1]
    }
  }
}

SPY_18_20 <- as.data.frame(SPY_18_20)

# Old output retained for easy rollback:
# usethis::use_data(spy, overwrite = TRUE)
usethis::use_data(SPY_18_20, overwrite = TRUE)
