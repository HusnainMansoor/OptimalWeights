Gathering Data from Yahoo Finance and Generating Optimal Weights for Alpha, Beta. We then produce yearly weights


Libraries and ETF List

library(PerformanceAnalytics)
library(quantmod)
library(DEoptim)
PerformanceAnalytics: Provides tools for performance and risk analysis of financial portfolios.
quantmod: Used for quantitative financial modeling, specifically for retrieving and managing financial data.
DEoptim: Implements Differential Evolution optimization, useful for portfolio optimization problems.

Define the list of ETFs

etfs <- list(
  GOOG <- c("GOOG"),
  MSFT <- c("MSFT"),
  ...
)

You define a list of ETFs (in this case, individual stocks from various sectors).

Load ETF data

e <- new.env()
getSymbols(unlist(etfs), env = e, from = "2021-01-01")
new.env(): Creates a new environment to store the data.
getSymbols(): Downloads historical stock data for each symbol from the list, starting from January 1, 2021.
Extract adjusted monthly prices

aa <- eapply(e, function(x) to.monthly(x, name = names(x)))
port <- do.call(merge, lapply(aa, Ad))
colnames(port) <- gsub(".Adjusted", "", names(port))
port <- ROC(port, type = "discrete")
port[is.na(port)] <- 0
eapply(): Applies the to.monthly function to each dataset to convert daily prices to monthly.
merge: Combines all monthly datasets into one.
Ad: Extracts adjusted closing prices.
ROC(): Calculates the rate of change, converting prices to returns.
port[is.na(port)] <- 0: Replaces any NA values with 0.

Group ETFs by sectors
# Process each stock individually
GOOG_data <- port[, GOOG]
GOOG_data <- reclass(coredata(GOOG_data) %*% rep(1/ncol(GOOG_data), ncol(GOOG_data)), match.to = GOOG_data)
...
For each stock, it extracts the data and normalizes it by applying equal weights.

Merge all stocks

ALL <- merge(GOOG_data, MSFT_data, AAPL_data, NVDA_data, AMZN_data, ...)
Combines all individual stock datasets into a single dataframe.

Plot performance summary

charts.PerformanceSummary(ALL, geometric = FALSE, colorset = rich12equal, cex.legend = 0.45)
Plots a performance summary of the combined dataset.

Generate random weights

random_weights <- runif(ncol(ALL))
random_weights <- random_weights / sum(random_weights)
print(data.frame(Ticker = colnames(ALL), Weight = random_weights))
Generates a random weight for each stock and normalizes them so that they sum to 1.

Calculate weighted returns

WTD <- as.matrix(ALL) %*% random_weights
WTD <- reclass(WTD, match.to = ALL)
Calculates the weighted returns of the portfolio using the random weights.

Benchmark against SPY

SPY <- ROC(Ad(to.monthly(getSymbols("^GSPC", auto.assign = FALSE, from = "2021-01-01"))), type = "discrete")
SPY[is.na(SPY)] <- 0
colnames(SPY) <- "SPY"
charts.PerformanceSummary(merge(WTD, SPY), geometric = FALSE)
Downloads SPY data, calculates its returns, and plots the performance of the weighted portfolio against SPY.

Optimization
Return optimization function

toOptim = function(n) {
  if (length(n) != ncol(ALL)) {
    stop("Number of elements in 'n' must match the number of columns in ALL")
  }
  WTD <- reclass(coredata(ALL) %*% n, match.to = ALL)
  wt.penalty = 100 * (1 - sum(n))^2
  return(-colSums(WTD) + wt.penalty)
}
Defines the function to optimize returns. It includes a penalty term to ensure the weights sum to 1.

Differential Evolution optimization for returns

r <- DEoptim(toOptim, lower = rep(0.00, ncol(ALL)), upper = rep(0.25, ncol(ALL)),
             control = list(itermax = 10000), fnMap = fnmap_f)
r$optim$bestmem
sum(r$optim$bestmem)
res <- reclass(coredata(ALL) %*% as.numeric(r$optim$bestmem), match.to = ALL)
colnames(res) <- "bestRet"
charts.PerformanceSummary(merge(WTD, SPY, res), geometric = FALSE)
Runs the optimization, extracts the best weights, calculates the optimized returns, and plots the performance.

Beta optimization function

BETAoptim = function(n) {
  WTD <- reclass(coredata(ALL) %*% n, match.to = ALL)
  obj = CAPM.beta(Ra = WTD, Rb = SPY, Rf = 0)
  wt.penalty = 100 * (1 - sum(n))^2
  return(obj + wt.penalty)
}
Defines the function to optimize the portfolio beta.

Differential Evolution optimization for beta

r2 <- DEoptim(BETAoptim, lower = rep(0.00, ncol(ALL)), upper = rep(0.25, ncol(ALL)),
              control = list(itermax = 5000), fnMap = fnmap_f)
r2$optim$bestmem
sum(r2$optim$bestmem)
res2 <- reclass(coredata(ALL) %*% as.numeric(r2$optim$bestmem), match.to = ALL)
colnames(res2) <- "BETA"
charts.PerformanceSummary(merge(WTD, SPY, res, res2), geometric = FALSE)
Runs the optimization for beta, extracts the best weights, calculates the optimized beta portfolio, and plots the performance.

Yearly Analysis
Repeats the above process for each year (2021 and 2022) by:

Loading data for the specific year.
Extracting monthly prices.
Merging stock data for the year.
Defining optimization functions.
Running optimization for returns and beta.
Plotting performance summaries for each year.
