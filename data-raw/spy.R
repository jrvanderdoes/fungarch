library(tidyverse)

## code to prepare `spy` dataset goes here
SPY2018 <- read.csv("./data-raw/SPY_2018.csv")
SPY2019 <- read.csv("./data-raw/SPY_2019.csv")
SPY2020 <- read.csv("./data-raw/SPY_2020.csv")

# Create Temp File
spy <- rbind(SPY2018, SPY2019, SPY2020)
spy$time <- as.POSIXct(spy$Local.time,
                       format = "%d.%m.%Y %H:%M:%OS ",
                       tz = "UTC")

# Organize
data = spy[,c("time","Close")]
# base drop
data <- data[as.numeric(format(data$time, '%H'))>=9 &
               as.numeric(format(data$time, '%H'))<=15, ]
data <- data[!(as.numeric(format(data$time, '%H'))==9 &
                 as.numeric(format(data$time, '%M'))<30),]
drop_vals = c()
dates = unique(format(data$time,'%Y-%m-%d'))
for(i in dates){
  cat(i,', ')
  if(length(unique(data[format(data$time,'%Y-%m-%d')==i,'Close']))==1){
    drop_vals = c(drop_vals, i)
  }
}
data <- data[!(format(data$time,'%Y-%m-%d') %in% drop_vals), ]
data$clock <- format(data$time,'%H %M')
data$day <- format(data$time,'%m-%d-%Y')

data_use <- data[,c("Close","clock","day")] %>%
  pivot_wider(names_from = day, values_from = Close)
data_use <- data_use[,-1] %>% as.data.frame()

# Widen and fix
stock_wide <- data_use
for(i in 1:ncol(stock_wide)){
  # If any missing
  if(sum(is.na(stock_wide[,i]))>0){
    missing_idx <- which(is.na(stock_wide[,i]))

    nonmissing <- na.omit(stock_wide[,i])[[1]]
    first_val <- nonmissing[1]

    for(j in missing_idx){
      if(j==1) stock_wide[j,i] <- first_val
      else stock_wide[j,i] <- stock_wide[j-1,i]
    }
  }
}

# Prep data
spy <- as.data.frame(stock_wide)

usethis::use_data(spy, overwrite = TRUE)
