# Convert SPY 1-minute data to the package's 10-minute CSV format.
#
# Put SPY_clean_1min.csv in data-raw/ before running this file. The input
# should contain a timestamp and a close-price column. Open/high/low/volume
# are optional; the output close is the last observed 1-minute close in each
# 10-minute interval.

input_file <- file.path("data-raw", "SPY_clean_1min.csv")
output_dir <- "data-raw"

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file,
       "\nPlace SPY_clean_1min.csv in data-raw/ and run this script from the package root.")
}

find_column <- function(data, candidates, required = TRUE) {
  normalized <- tolower(trimws(names(data)))
  candidate_names <- tolower(trimws(candidates))
  idx <- which(normalized %in% candidate_names)[1]

  if (required && is.na(idx)) {
    stop("Could not find a column matching: ", paste(candidates, collapse = ", "),
         "\nAvailable columns: ", paste(names(data), collapse = ", "))
  }

  if (is.na(idx)) return(NULL)
  names(data)[idx]
}

parse_timestamp <- function(x) {
  if (inherits(x, "POSIXt")) return(as.POSIXct(x, tz = "UTC"))

  # Numeric timestamps are interpreted as Unix seconds, milliseconds, or
  # microseconds based on their magnitude.
  if (is.numeric(x)) {
    scale <- ifelse(max(abs(x), na.rm = TRUE) > 1e14, 1e6,
                    ifelse(max(abs(x), na.rm = TRUE) > 1e11, 1e3, 1))
    return(as.POSIXct(x / scale, origin = "1970-01-01", tz = "UTC"))
  }

  x <- trimws(as.character(x))
  # SPY_clean_1min.csv contains exchange-local timestamps without an offset.
  # Treat those as New York time; timestamps that explicitly include Z or an
  # offset are parsed as UTC so their stated timezone is not changed.
  has_timezone <- grepl("([zZ]|[+-][0-9]{2}:?[0-9]{2})$", x)
  parse_tz <- if (any(has_timezone)) "UTC" else "America/New_York"
  parsed <- as.POSIXct(x, tz = parse_tz, tryFormats = c(
    "%Y-%m-%d %H:%M:%OS", "%Y-%m-%dT%H:%M:%OS",
    "%Y/%m/%d %H:%M:%OS", "%m/%d/%Y %H:%M:%OS",
    "%d.%m.%Y %H:%M:%OS", "%Y-%m-%d", "%m/%d/%Y"
  ))

  if (all(is.na(parsed))) {
    stop("Could not parse any timestamps. Check the timestamp format in the input file.")
  }
  parsed
}

format_local_time <- function(x) {
  # Match the existing SPY_2018.csv through SPY_2020.csv format.
  format(x, "%d.%m.%Y %H:%M:%OS3 GMT%z", tz = "America/New_York")
}

make_ten_minute_data <- function(data, timestamp_col, close_col,
                                 open_col = NULL, high_col = NULL,
                                 low_col = NULL, volume_col = NULL) {
  timestamp <- parse_timestamp(data[[timestamp_col]])
  close <- suppressWarnings(as.numeric(data[[close_col]]))

  if (anyNA(timestamp) || anyNA(close)) {
    keep <- !is.na(timestamp) & !is.na(close)
    timestamp <- timestamp[keep]
    close <- close[keep]
    data <- data[keep, , drop = FALSE]
  }

  if (!length(timestamp)) stop("No valid timestamp/close rows were found.")

  # Use elapsed 10-minute UTC buckets. The output time zone is only a display
  # choice; this avoids daylight-saving-time bucket errors.
  bucket <- as.POSIXct(floor(as.numeric(timestamp) / 600) * 600,
                       origin = "1970-01-01", tz = "UTC")
  order_idx <- order(timestamp)
  data <- data[order_idx, , drop = FALSE]
  timestamp <- timestamp[order_idx]
  close <- close[order_idx]
  bucket <- bucket[order_idx]

  groups <- split(seq_along(bucket), bucket)
  rows <- lapply(groups, function(idx) {
    out <- data.frame(
      bucket = bucket[idx[1]],
      Open = close[idx[1]],
      High = max(close[idx], na.rm = TRUE),
      Low = min(close[idx], na.rm = TRUE),
      Close = close[idx[length(idx)]],
      Volume = NA_real_,
      stringsAsFactors = FALSE
    )

    # If OHLC exists, aggregate it conventionally while retaining the
    # last observed close. Otherwise all price fields are derived from close.
    if (!is.null(open_col)) out$Open <- as.numeric(data[[open_col]][idx[1]])
    if (!is.null(high_col)) out$High <- max(as.numeric(data[[high_col]][idx]), na.rm = TRUE)
    if (!is.null(low_col)) out$Low <- min(as.numeric(data[[low_col]][idx]), na.rm = TRUE)
    if (!is.null(volume_col)) out$Volume <- sum(as.numeric(data[[volume_col]][idx]), na.rm = TRUE)
    out
  })

  result <- do.call(rbind, rows)
  result$year <- format(result$bucket, "%Y", tz = "America/New_York")
  result$`Local time` <- format_local_time(result$bucket)
  result <- result[, c("year", "Local time", "Open", "High", "Low", "Close", "Volume")]
  rownames(result) <- NULL
  result
}

raw <- read.csv(input_file, check.names = FALSE, stringsAsFactors = FALSE)

timestamp_col <- find_column(raw, c("timestamp", "datetime", "date", "time", "localtime"))
close_col <- find_column(raw, c("close", "closeprice", "adjclose"))
open_col <- find_column(raw, "open", required = FALSE)
high_col <- find_column(raw, "high", required = FALSE)
low_col <- find_column(raw, "low", required = FALSE)
volume_col <- find_column(raw, c("volume", "vol"), required = FALSE)

ten_minute <- make_ten_minute_data(
  raw, timestamp_col, close_col, open_col, high_col, low_col, volume_col
)

years <- as.integer(ten_minute$year)
ten_minute$year <- NULL

for (year in 2021:2024) {
  output_file <- file.path(output_dir, paste0("SPY_", year, ".csv"))
  if (file.exists(output_file)) {
    message("Skipping existing file: ", output_file)
    next
  }

  result <- ten_minute[years == year, , drop = FALSE]
  if (!nrow(result)) {
    warning("No rows found for ", year, "; no file written.")
    next
  }

  write.csv(result, output_file, row.names = FALSE, na = "")
  message("Wrote ", nrow(result), " rows to ", output_file)
}
