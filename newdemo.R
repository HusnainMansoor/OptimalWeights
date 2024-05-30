library(PerformanceAnalytics)
library(quantmod)
library(DEoptim)

# Define the list of ETFs
etfs <- list(
  GOOG <- c("GOOG"),
  MSFT <- c("MSFT"),
  AAPL <- c("AAPL"),
  NVDA <- c("NVDA"),
  AMZN <- c("AMZN"),
  CRM <- c("CRM"),
  IBM <- c("IBM"),
  ADBE <- c("ADBE"),
  INTC <- c("INTC"),
  CSCO <- c("CSCO"),
  ORCL <- c("ORCL"),
  
  VZ <- c("VZ"),
  CMCSA <- c("CMCSA"),
  T <- c("T"),
  CHTR <- c("CHTR"),
  TMUS <- c("TMUS"),
  DIS <- c("DIS"),
  
  JPM <- c("JPM"),
  V <- c("V"),
  BAC <- c("BAC"),
  MA <- c("MA"),
  PYPL <- c("PYPL"),
  GS <- c("GS"),
  MS <- c("MS"),
  WFC <- c("WFC"),
  AXP <- c("AXP"),
  
  PLD <- c("PLD"),
  AMT <- c("AMT"),
  CCI <- c("CCI"),
  SPG <- c("SPG"),
  EQIX <- c("EQIX"),
  PSA <- c("PSA"),
  O <- c("O"),
  DLR <- c("DLR"),
  AVB <- c("AVB"),
  WELL <- c("WELL"),
  
  CVX <- c("CVX"),
  XOM <- c("XOM"),
  EOG <- c("EOG"),
  COP <- c("COP"),
  DOW <- c("DOW"),
  FCX <- c("FCX"),
  KMI <- c("KMI"),
  PSX <- c("PSX"),
  MPC <- c("MPC"),
  APD <- c("APD")
)


# Load ETF data
e <- new.env()
getSymbols(unlist(etfs), env = e, from = "2021-01-01")

# Extract adjusted monthly prices
aa <- eapply(e, function(x) to.monthly(x, name = names(x)))
port <- do.call(merge, lapply(aa, Ad))
colnames(port) <- gsub(".Adjusted", "", names(port))
port <- ROC(port, type = "discrete")
port[is.na(port)] <- 0

# Group ETFs by sectors
# Process each stock individually
GOOG_data <- port[, GOOG]
GOOG_data <- reclass(coredata(GOOG_data) %*% rep(1/ncol(GOOG_data), ncol(GOOG_data)), match.to = GOOG_data)

MSFT_data <- port[, MSFT]
MSFT_data <- reclass(coredata(MSFT_data) %*% rep(1/ncol(MSFT_data), ncol(MSFT_data)), match.to = MSFT_data)

AAPL_data <- port[, AAPL]
AAPL_data <- reclass(coredata(AAPL_data) %*% rep(1/ncol(AAPL_data), ncol(AAPL_data)), match.to = AAPL_data)

NVDA_data <- port[, NVDA]
NVDA_data <- reclass(coredata(NVDA_data) %*% rep(1/ncol(NVDA_data), ncol(NVDA_data)), match.to = NVDA_data)

AMZN_data <- port[, AMZN]
AMZN_data <- reclass(coredata(AMZN_data) %*% rep(1/ncol(AMZN_data), ncol(AMZN_data)), match.to = AMZN_data)

CRM_data <- port[, CRM]
CRM_data <- reclass(coredata(CRM_data) %*% rep(1/ncol(CRM_data), ncol(CRM_data)), match.to = CRM_data)

IBM_data <- port[, IBM]
IBM_data <- reclass(coredata(IBM_data) %*% rep(1/ncol(IBM_data), ncol(IBM_data)), match.to = IBM_data)

ADBE_data <- port[, ADBE]
ADBE_data <- reclass(coredata(ADBE_data) %*% rep(1/ncol(ADBE_data), ncol(ADBE_data)), match.to = ADBE_data)

INTC_data <- port[, INTC]
INTC_data <- reclass(coredata(INTC_data) %*% rep(1/ncol(INTC_data), ncol(INTC_data)), match.to = INTC_data)

CSCO_data <- port[, CSCO]
CSCO_data <- reclass(coredata(CSCO_data) %*% rep(1/ncol(CSCO_data), ncol(CSCO_data)), match.to = CSCO_data)

ORCL_data <- port[, ORCL]
ORCL_data <- reclass(coredata(ORCL_data) %*% rep(1/ncol(ORCL_data), ncol(ORCL_data)), match.to = ORCL_data)

VZ_data <- port[, VZ]
VZ_data <- reclass(coredata(VZ_data) %*% rep(1/ncol(VZ_data), ncol(VZ_data)), match.to = VZ_data)

CMCSA_data <- port[, CMCSA]
CMCSA_data <- reclass(coredata(CMCSA_data) %*% rep(1/ncol(CMCSA_data), ncol(CMCSA_data)), match.to = CMCSA_data)

T_data <- port[, T]
T_data <- reclass(coredata(T_data) %*% rep(1/ncol(T_data), ncol(T_data)), match.to = T_data)

CHTR_data <- port[, CHTR]
CHTR_data <- reclass(coredata(CHTR_data) %*% rep(1/ncol(CHTR_data), ncol(CHTR_data)), match.to = CHTR_data)

TMUS_data <- port[, TMUS]
TMUS_data <- reclass(coredata(TMUS_data) %*% rep(1/ncol(TMUS_data), ncol(TMUS_data)), match.to = TMUS_data)

DIS_data <- port[, DIS]
DIS_data <- reclass(coredata(DIS_data) %*% rep(1/ncol(DIS_data), ncol(DIS_data)), match.to = DIS_data)

JPM_data <- port[, JPM]
JPM_data <- reclass(coredata(JPM_data) %*% rep(1/ncol(JPM_data), ncol(JPM_data)), match.to = JPM_data)

V_data <- port[, V]
V_data <- reclass(coredata(V_data) %*% rep(1/ncol(V_data), ncol(V_data)), match.to = V_data)

BAC_data <- port[, BAC]
BAC_data <- reclass(coredata(BAC_data) %*% rep(1/ncol(BAC_data), ncol(BAC_data)), match.to = BAC_data)

MA_data <- port[, MA]
MA_data <- reclass(coredata(MA_data) %*% rep(1/ncol(MA_data), ncol(MA_data)), match.to = MA_data)

PYPL_data <- port[, PYPL]
PYPL_data <- reclass(coredata(PYPL_data) %*% rep(1/ncol(PYPL_data), ncol(PYPL_data)), match.to = PYPL_data)


GS_data <- port[, GS]
GS_data <- reclass(coredata(GS_data) %*% rep(1/ncol(GS_data), ncol(GS_data)), match.to = GS_data)

MS_data <- port[, MS]
MS_data <- reclass(coredata(MS_data) %*% rep(1/ncol(MS_data), ncol(MS_data)), match.to = MS_data)

WFC_data <- port[, WFC]
WFC_data <- reclass(coredata(WFC_data) %*% rep(1/ncol(WFC_data), ncol(WFC_data)), match.to = WFC_data)

AXP_data <- port[, AXP]
AXP_data <- reclass(coredata(AXP_data) %*% rep(1/ncol(AXP_data), ncol(AXP_data)), match.to = AXP_data)

PLD_data <- port[, PLD]
PLD_data <- reclass(coredata(PLD_data) %*% rep(1/ncol(PLD_data), ncol(PLD_data)), match.to = PLD_data)

AMT_data <- port[, AMT]
AMT_data <- reclass(coredata(AMT_data) %*% rep(1/ncol(AMT_data), ncol(AMT_data)), match.to = AMT_data)

CCI_data <- port[, CCI]
CCI_data <- reclass(coredata(CCI_data) %*% rep(1/ncol(CCI_data), ncol(CCI_data)), match.to = CCI_data)

SPG_data <- port[, SPG]
SPG_data <- reclass(coredata(SPG_data) %*% rep(1/ncol(SPG_data), ncol(SPG_data)), match.to = SPG_data)

EQIX_data <- port[, EQIX]
EQIX_data <- reclass(coredata(EQIX_data) %*% rep(1/ncol(EQIX_data), ncol(EQIX_data)), match.to = EQIX_data)

PSA_data <- port[, PSA]
PSA_data <- reclass(coredata(PSA_data) %*% rep(1/ncol(PSA_data), ncol(PSA_data)), match.to = PSA_data)

O_data <- port[, O]
O_data <- reclass(coredata(O_data) %*% rep(1/ncol(O_data), ncol(O_data)), match.to = O_data)

DLR_data <- port[, DLR]
DLR_data <- reclass(coredata(DLR_data) %*% rep(1/ncol(DLR_data), ncol(DLR_data)), match.to = DLR_data)

AVB_data <- port[, AVB]
AVB_data <- reclass(coredata(AVB_data) %*% rep(1/ncol(AVB_data), ncol(AVB_data)), match.to = AVB_data)

WELL_data <- port[, WELL]
WELL_data <- reclass(coredata(WELL_data) %*% rep(1/ncol(WELL_data), ncol(WELL_data)), match.to = WELL_data)

CVX_data <- port[, CVX]
CVX_data <- reclass(coredata(CVX_data) %*% rep(1/ncol(CVX_data), ncol(CVX_data)), match.to = CVX_data)

XOM_data <- port[, XOM]
XOM_data <- reclass(coredata(XOM_data) %*% rep(1/ncol(XOM_data), ncol(XOM_data)), match.to = XOM_data)

EOG_data <- port[, EOG]
EOG_data <- reclass(coredata(EOG_data) %*% rep(1/ncol(EOG_data), ncol(EOG_data)), match.to = EOG_data)

COP_data <- port[, COP]
COP_data <- reclass(coredata(COP_data) %*% rep(1/ncol(COP_data), ncol(COP_data)), match.to = COP_data)

DOW_data <- port[, DOW]
DOW_data <- reclass(coredata(DOW_data) %*% rep(1/ncol(DOW_data), ncol(DOW_data)), match.to = DOW_data)

FCX_data <- port[, FCX]
FCX_data <- reclass(coredata(FCX_data) %*% rep(1/ncol(FCX_data), ncol(FCX_data)), match.to = FCX_data)

KMI_data <- port[, KMI]
KMI_data <- reclass(coredata(KMI_data) %*% rep(1/ncol(KMI_data), ncol(KMI_data)), match.to = KMI_data)

PSX_data <- port[, PSX]
PSX_data <- reclass(coredata(PSX_data) %*% rep(1/ncol(PSX_data), ncol(PSX_data)), match.to = PSX_data)

MPC_data <- port[, MPC]
MPC_data <- reclass(coredata(MPC_data) %*% rep(1/ncol(MPC_data), ncol(MPC_data)), match.to = MPC_data)

APD_data <- port[, APD]
APD_data <- reclass(coredata(APD_data) %*% rep(1/ncol(APD_data), ncol(APD_data)), match.to = APD_data)

# Merge all stocks
ALL <- merge(GOOG_data, MSFT_data, AAPL_data, NVDA_data, AMZN_data, CRM_data, IBM_data, ADBE_data, INTC_data, CSCO_data, ORCL_data, 
             VZ_data, CMCSA_data, T_data, CHTR_data, TMUS_data, DIS_data, 
             JPM_data, V_data, BAC_data, MA_data, PYPL_data, GS_data, MS_data, WFC_data, AXP_data, 
             PLD_data, AMT_data, CCI_data, SPG_data, EQIX_data, PSA_data, O_data, DLR_data, AVB_data, WELL_data, 
             CVX_data, XOM_data, EOG_data, COP_data, DOW_data, FCX_data, KMI_data, PSX_data, MPC_data, APD_data)


# Plot performance summary
charts.PerformanceSummary(ALL, geometric = FALSE, colorset = rich12equal, cex.legend = 0.45)

# Generate 77 random weights
random_weights <- runif(ncol(ALL))
random_weights <- random_weights / sum(random_weights)  # Normalize to sum up to 1

# Display the generated random weights
print(data.frame(Ticker = colnames(ALL), Weight = random_weights))

# Calculate weighted returns
WTD <- as.matrix(ALL) %*% random_weights
WTD <- reclass(WTD, match.to = ALL)

SPY <- ROC(Ad(to.monthly(getSymbols("^GSPC", auto.assign = FALSE, from = "2021-01-01"))), type = "discrete")
SPY[is.na(SPY)] <- 0
colnames(SPY) <- "SPY"
charts.PerformanceSummary(merge(WTD,SPY), geometric = FALSE)

toOptim = function(n) {
  if (length(n) != ncol(ALL)) {
    stop("Number of elements in 'n' must match the number of columns in ALL")
  }
  WTD <- reclass(coredata(ALL) %*% n, match.to = ALL)
  wt.penalty = 100 * (1 - sum(n))^2
  return(-colSums(WTD) + wt.penalty)
}


fnmap_f <- function(x){c(round(x,2))}

r <- DEoptim(toOptim, lower = rep(0.00, ncol(ALL)),upper = rep(0.25, ncol(ALL)),
             control = list(itermax = 10000), fnMap = fnmap_f)

r$optim$bestmem
sum(r$optim$bestmem)
res <- reclass(coredata(ALL) %*% as.numeric(r$optim$bestmem), match.to = ALL)
colnames(res) <- "bestRet"
charts.PerformanceSummary(merge(WTD,SPY,res), geometric = FALSE)

BETAoptim = function(n)
{
  WTD <- reclass(coredata(ALL) %*% c(n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], n[12], n[13], n[14], n[15], n[16], n[17], n[18], n[19], n[20], n[21], n[22], n[23], n[24], n[25], n[26], n[27], n[28], n[29], n[30], n[31], n[32], n[33], n[34], n[35], n[36], n[37], n[38], n[39], n[40], n[41], n[42], n[43], n[44], n[45], n[46]), match.to = ALL)
  obj = CAPM.beta(Ra = WTD, Rb = SPY, Rf = 0) 
  wt.penalty = 100 * (1 - sum(n))^2
  return(obj + wt.penalty)
}


r2 <- DEoptim(BETAoptim, lower = rep(0.00, ncol(ALL)), upper = rep(0.25, ncol(ALL)),
              control = list(itermax = 5000), fnMap = fnmap_f)

r2$optim$bestmem
sum(r2$optim$bestmem)
res2 <- reclass(coredata(ALL) %*% as.numeric(r2$optim$bestmem), match.to = ALL)
colnames(res2) <- "BETA"
charts.PerformanceSummary(merge(WTD,SPY,res, res2), geometric = FALSE)

###############################################################################################


#Extracting for each year

#################


# Load ETF data for the specified time period
e <- new.env()
getSymbols(unlist(etfs), env = e, from = "2021-01-01", to = "2022-01-01")

# Extract adjusted monthly prices
aa <- eapply(e, function(x) to.monthly(x, name = names(x)))
port <- do.call(merge, lapply(aa, Ad))
colnames(port) <- gsub(".Adjusted", "", names(port))
port <- ROC(port, type = "discrete")
port[is.na(port)] <- 0

ALL2021 <- merge(GOOG_data, MSFT_data, AAPL_data, NVDA_data, AMZN_data, CRM_data, IBM_data, ADBE_data, INTC_data, CSCO_data, ORCL_data, 
             VZ_data, CMCSA_data, T_data, CHTR_data, TMUS_data, DIS_data, 
             JPM_data, V_data, BAC_data, MA_data, PYPL_data, GS_data, MS_data, WFC_data, AXP_data, 
             PLD_data, AMT_data, CCI_data, SPG_data, EQIX_data, PSA_data, O_data, DLR_data, AVB_data, WELL_data, 
             CVX_data, XOM_data, EOG_data, COP_data, DOW_data, FCX_data, KMI_data, PSX_data, MPC_data, APD_data)

toOptim = function(n) {
  if (length(n) != ncol(ALL2021)) {
    stop("Number of elements in 'n' must match the number of columns in ALL")
  }
  WTD <- reclass(coredata(ALL2021) %*% n, match.to = ALL2021)
  wt.penalty = 100 * (1 - sum(n))^2
  return(-colSums(WTD) + wt.penalty)
}

fnmap_f <- function(x){c(round(x,2))}

r <- DEoptim(toOptim, lower = rep(0.00, ncol(ALL2021)),upper = rep(0.25, ncol(ALL2021)),
             control = list(itermax = 50), fnMap = fnmap_f)

r$optim$bestmem
sum(r$optim$bestmem)
res <- reclass(coredata(ALL2021) %*% as.numeric(r$optim$bestmem), match.to = ALL2021)
colnames(res) <- "bestRet"
charts.PerformanceSummary(merge(WTD,SPY,res), geometric = FALSE)

BETAoptim = function(n)
{
  WTD <- reclass(coredata(ALL2021) %*% c(n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], n[12], n[13], n[14], n[15], n[16], n[17], n[18], n[19], n[20], n[21], n[22], n[23], n[24], n[25], n[26], n[27], n[28], n[29], n[30], n[31], n[32], n[33], n[34], n[35], n[36], n[37], n[38], n[39], n[40], n[41], n[42], n[43], n[44], n[45], n[46]), match.to = ALL2024)
  obj = CAPM.beta(Ra = WTD,Rb=SPY, Rf =0) 
  wt.penalty = 100 * (1 - sum(n))^2
  return(obj + wt.penalty)
}

r2 <- DEoptim(BETAoptim, lower = rep(0.00, ncol(ALL2021)), upper = rep(0.25, ncol(ALL2021)),
              control = list(itermax = 500), fnMap = fnmap_f)

r2$optim$bestmem
sum(r2$optim$bestmem)
res2 <- reclass(coredata(ALL2021) %*% as.numeric(r2$optim$bestmem), match.to = ALL2021)
colnames(res) <- "BETA"
charts.PerformanceSummary(merge(WTD,SPY,res, res2), geometric = FALSE)

########################################################################

e <- new.env()
getSymbols(unlist(etfs), env = e, from = "2022-01-01", to = "2023-01-01")

# Extract adjusted monthly prices
aa <- eapply(e, function(x) to.monthly(x, name = names(x)))
port <- do.call(merge, lapply(aa, Ad))
colnames(port) <- gsub(".Adjusted", "", names(port))
port <- ROC(port, type = "discrete")
port[is.na(port)] <- 0

ALL2022 <- merge(GOOG_data, MSFT_data, AAPL_data, NVDA_data, AMZN_data, CRM_data, IBM_data, ADBE_data, INTC_data, CSCO_data, ORCL_data, 
                 VZ_data, CMCSA_data, T_data, CHTR_data, TMUS_data, DIS_data, 
                 JPM_data, V_data, BAC_data, MA_data, PYPL_data, GS_data, MS_data, WFC_data, AXP_data, 
                 PLD_data, AMT_data, CCI_data, SPG_data, EQIX_data, PSA_data, O_data, DLR_data, AVB_data, WELL_data, 
                 CVX_data, XOM_data, EOG_data, COP_data, DOW_data, FCX_data, KMI_data, PSX_data, MPC_data, APD_data)

toOptim = function(n) {
  if (length(n) != ncol(ALL2021)) {
    stop("Number of elements in 'n' must match the number of columns in ALL")
  }
  WTD <- reclass(coredata(ALL2022) %*% n, match.to = ALL2022)
  wt.penalty = 100 * (1 - sum(n))^2
  return(-colSums(WTD) + wt.penalty)
}

fnmap_f <- function(x){c(round(x,2))}

r <- DEoptim(toOptim, lower = rep(0.00, ncol(ALL2022)),upper = rep(0.25, ncol(ALL2022)),
             control = list(itermax = 500), fnMap = fnmap_f)

r$optim$bestmem
sum(r$optim$bestmem)
res <- reclass(coredata(ALL2022) %*% as.numeric(r$optim$bestmem), match.to = ALL2022)
colnames(res) <- "bestRet"
charts.PerformanceSummary(merge(WTD,SPY,res), geometric = FALSE)

BETAoptim = function(n)
{
  WTD <- reclass(coredata(ALL2022) %*% c(n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], n[12], n[13], n[14], n[15], n[16], n[17], n[18], n[19], n[20], n[21], n[22], n[23], n[24], n[25], n[26], n[27], n[28], n[29], n[30], n[31], n[32], n[33], n[34], n[35], n[36], n[37], n[38], n[39], n[40], n[41], n[42], n[43], n[44], n[45], n[46]), match.to = ALL2024)
  obj = CAPM.beta(Ra = WTD,Rb=SPY, Rf =0) 
  wt.penalty = 100 * (1 - sum(n))^2
  return(obj + wt.penalty)
}

r2 <- DEoptim(BETAoptim, lower = rep(0.00, ncol(ALL2022)), upper = rep(0.25, ncol(ALL2022)),
              control = list(itermax = 500), fnMap = fnmap_f)

r2$optim$bestmem
sum(r2$optim$bestmem)
res2 <- reclass(coredata(ALL2022) %*% as.numeric(r2$optim$bestmem), match.to = ALL2022)
colnames(res) <- "BETA"
charts.PerformanceSummary(merge(WTD,SPY,res, res2), geometric = FALSE)

###############################################

# Load ETF data for the specified time period
e <- new.env()
getSymbols(unlist(etfs), env = e, from = "2023-01-01", to = "2023-12-01")

# Extract adjusted monthly prices
aa <- eapply(e, function(x) to.monthly(x, name = names(x)))
port <- do.call(merge, lapply(aa, Ad))
colnames(port) <- gsub(".Adjusted", "", names(port))
port <- ROC(port, type = "discrete")
port[is.na(port)] <- 0

ALL2023 <- merge(GOOG_data, MSFT_data, AAPL_data, NVDA_data, AMZN_data, CRM_data, IBM_data, ADBE_data, INTC_data, CSCO_data, ORCL_data, 
                 VZ_data, CMCSA_data, T_data, CHTR_data, TMUS_data, DIS_data, 
                 JPM_data, V_data, BAC_data, MA_data, PYPL_data, GS_data, MS_data, WFC_data, AXP_data, 
                 PLD_data, AMT_data, CCI_data, SPG_data, EQIX_data, PSA_data, O_data, DLR_data, AVB_data, WELL_data, 
                 CVX_data, XOM_data, EOG_data, COP_data, DOW_data, FCX_data, KMI_data, PSX_data, MPC_data, APD_data)

toOptim = function(n) {
  if (length(n) != ncol(ALL2023)) {
    stop("Number of elements in 'n' must match the number of columns in ALL")
  }
  WTD <- reclass(coredata(ALL2023) %*% n, match.to = ALL2023)
  wt.penalty = 100 * (1 - sum(n))^2
  return(-colSums(WTD) + wt.penalty)
}

fnmap_f <- function(x){c(round(x,2))}

r <- DEoptim(toOptim, lower = rep(0.00, ncol(ALL2023)),upper = rep(0.25, ncol(ALL2023)),
             control = list(itermax = 500), fnMap = fnmap_f)

r$optim$bestmem
sum(r$optim$bestmem)
res <- reclass(coredata(ALL2023) %*% as.numeric(r$optim$bestmem), match.to = ALL2023)
colnames(res) <- "bestRet"
charts.PerformanceSummary(merge(WTD,SPY,res), geometric = FALSE)

BETAoptim = function(n)
{
  WTD <- reclass(coredata(ALL2023) %*% c(n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], n[12], n[13], n[14], n[15], n[16], n[17], n[18], n[19], n[20], n[21], n[22], n[23], n[24], n[25], n[26], n[27], n[28], n[29], n[30], n[31], n[32], n[33], n[34], n[35], n[36], n[37], n[38], n[39], n[40], n[41], n[42], n[43], n[44], n[45], n[46]), match.to = ALL2024)
  obj = CAPM.beta(Ra = WTD,Rb=SPY, Rf =0) 
  wt.penalty = 100 * (1 - sum(n))^2
  return(obj + wt.penalty)
}

r2 <- DEoptim(BETAoptim, lower = rep(0.00, ncol(ALL2023)), upper = rep(0.25, ncol(ALL2023)),
              control = list(itermax = 500), fnMap = fnmap_f)

r2$optim$bestmem
sum(r2$optim$bestmem)
res2 <- reclass(coredata(ALL2023) %*% as.numeric(r2$optim$bestmem), match.to = ALL2023)
colnames(res) <- "BETA"
charts.PerformanceSummary(merge(WTD,SPY,res, res2), geometric = FALSE)

#################################
# Load ETF data for the specified time period
e <- new.env()
getSymbols(unlist(etfs), env = e, from = "2024-01-01", to = "2024-05-01")

# Extract adjusted monthly prices
aa <- eapply(e, function(x) to.monthly(x, name = names(x)))
port <- do.call(merge, lapply(aa, Ad))
colnames(port) <- gsub(".Adjusted", "", names(port))
port <- ROC(port, type = "discrete")
port[is.na(port)] <- 0

ALL2024 <- merge(GOOG_data, MSFT_data, AAPL_data, NVDA_data, AMZN_data, CRM_data, IBM_data, ADBE_data, INTC_data, CSCO_data, ORCL_data, 
                 VZ_data, CMCSA_data, T_data, CHTR_data, TMUS_data, DIS_data, 
                 JPM_data, V_data, BAC_data, MA_data, PYPL_data, GS_data, MS_data, WFC_data, AXP_data, 
                 PLD_data, AMT_data, CCI_data, SPG_data, EQIX_data, PSA_data, O_data, DLR_data, AVB_data, WELL_data, 
                 CVX_data, XOM_data, EOG_data, COP_data, DOW_data, FCX_data, KMI_data, PSX_data, MPC_data, APD_data)

toOptim = function(n) {
  if (length(n) != ncol(ALL2024)) {
    stop("Number of elements in 'n' must match the number of columns in ALL")
  }
  WTD <- reclass(coredata(ALL2024) %*% n, match.to = ALL2024)
  wt.penalty = 100 * (1 - sum(n))^2
  return(-colSums(WTD) + wt.penalty)
}

fnmap_f <- function(x){c(round(x,2))}

r <- DEoptim(toOptim, lower = rep(0.00, ncol(ALL2024)),upper = rep(0.25, ncol(ALL2024)),
             control = list(itermax = 500), fnMap = fnmap_f)

r$optim$bestmem
sum(r$optim$bestmem)
res <- reclass(coredata(ALL2024) %*% as.numeric(r$optim$bestmem), match.to = ALL2024)
colnames(res) <- "bestRet"
charts.PerformanceSummary(merge(WTD,SPY,res), geometric = FALSE)

BETAoptim = function(n)
{
  WTD <- reclass(coredata(ALL2024) %*% c(n[1], n[2], n[3], n[4], n[5], n[6], n[7], n[8], n[9], n[10], n[11], n[12], n[13], n[14], n[15], n[16], n[17], n[18], n[19], n[20], n[21], n[22], n[23], n[24], n[25], n[26], n[27], n[28], n[29], n[30], n[31], n[32], n[33], n[34], n[35], n[36], n[37], n[38], n[39], n[40], n[41], n[42], n[43], n[44], n[45], n[46]), match.to = ALL2024)
  obj = CAPM.beta(Ra = WTD,Rb=SPY, Rf =0) 
  wt.penalty = 100 * (1 - sum(n))^2
  return(obj + wt.penalty)
}

r2 <- DEoptim(BETAoptim, lower = rep(0.00, ncol(ALL2024)), upper = rep(0.25, ncol(ALL2024)),
              control = list(itermax = 1000), fnMap = fnmap_f)

r2$optim$bestmem
sum(r2$optim$bestmem)
res2 <- reclass(coredata(ALL2024) %*% as.numeric(r2$optim$bestmem), match.to = ALL2024)
colnames(res) <- "BETA"
charts.PerformanceSummary(merge(WTD,SPY,res, res2), geometric = FALSE)
################################################




