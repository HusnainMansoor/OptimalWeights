Gathering Data from Yahoo Finance and Generating Optimal Weights for Alpha and Beta
Libraries and ETF List
1.	PerformanceAnalytics: Provides tools for performance and risk analysis of financial portfolios.
2.	quantmod: Used for quantitative financial modeling, specifically for retrieving and managing financial data.
3.	DEoptim: Implements Differential Evolution optimization, useful for portfolio optimization problems.
Define the list of ETFs:
•	You define a list of ETFs, which in this case are individual stocks from various sectors (e.g., GOOG, MSFT).
Load ETF Data
•	Creating a new environment: new.env() creates a new environment to store the data.
•	Downloading historical stock data: getSymbols() downloads historical stock data for each symbol from the list, starting from January 1, 2021.
Extract Adjusted Monthly Prices
•	Convert daily to monthly prices: eapply() applies the to.monthly function to each dataset to convert daily prices to monthly.
•	Combine datasets: merge combines all monthly datasets into one.
•	Extract adjusted closing prices: Ad extracts adjusted closing prices.
•	Calculate returns: ROC() calculates the rate of change, converting prices to returns.
•	Replace NA values: port[is.na(port)] <- 0 replaces any NA values with 0.
Group ETFs by Sectors
•	Process each stock individually: For each stock, it extracts the data and normalizes it by applying equal weights.
Merge All Stocks
•	Combining datasets: Combines all individual stock datasets into a single dataframe.
Plot Performance Summary
•	Performance summary: charts.PerformanceSummary() plots a performance summary of the combined dataset.
Generate Random Weights
•	Random weight generation: Generates a random weight for each stock and normalizes them so that they sum to 1.
Calculate Weighted Returns
•	Weighted returns: Calculates the weighted returns of the portfolio using the random weights.
Benchmark Against SPY
•	Downloading SPY data: Downloads SPY data, calculates its returns, and plots the performance of the weighted portfolio against SPY.
Optimization
Return Optimization Function
•	Define function: Defines the function to optimize returns, including a penalty term to ensure the weights sum to 1.
Differential Evolution Optimization for Returns
•	Run optimization: Runs the optimization, extracts the best weights, calculates the optimized returns, and plots the performance.
Beta Optimization Function
•	Define function: Defines the function to optimize the portfolio beta, including a penalty term to ensure the weights sum to 1.
Differential Evolution Optimization for Beta
•	Run optimization: Runs the optimization for beta, extracts the best weights, calculates the optimized beta portfolio, and plots the performance.
Yearly Analysis
Repeats the above process for each year (2021 and 2022) by:
1.	Loading data for the specific year.
2.	Extracting monthly prices.
3.	Merging stock data for the year.
4.	Defining optimization functions.
5.	Running optimization for returns and beta.
6.	Plotting performance summaries for each year.

