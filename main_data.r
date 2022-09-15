### ======================================================================== ###
library("readr")
### ======================================================================== ###


### ====================================================================== ###
### data: price and volume (No. of shares) 
ticker <- "AAPL"
data <- na.omit(read_csv(paste("data/", ticker, ".csv", sep="")))
### ====================================================================== ###
dates <- as.Date(data$Date, format="%d/%m/%Y")
prices <- data$`Adj Close`
volumes <- data$Volume
### ====================================================================== ###
### process zero volumes [volume on the denominator]
### replace with the average of near data [not optimal]
volumes[which(volumes==0)] <- NA
while(sum(is.na(volumes))>0) {
    volumes <- na.approx(volumes, rule=2)
}
### ====================================================================== ###
### individual stock: adjust to dollar volume=share volume*price
volumesDollar <- volumes*prices
### ====================================================================== ###
returns <- diff(log(prices))   ### return series
n <- length(returns)   ### No. of observations
liquidity <- abs(returns)/volumesDollar[2:(n+1)]
### ======================================================================== ###


### ======================================================================== ###
### zero observations in liquidity [due to zero returns]
### replace with average of near observations or a very small value [not optimal]
while(sum(liquidity==0)>0) {
    # liquidity[liquidity==0] <- 1e-20
    liquidity[liquidity==0] <- (liquidity[which(liquidity==0)-1]
                                +liquidity[which(liquidity==0)+1])/2
}
### ======================================================================== ###

