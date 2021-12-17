library(e1071)
library(ggplot2)
library(doParallel)

registerDoParallel(cores = 8)

# --------- Preparing Env -------------
#Working dir
setwd("~/Desktop/code/Lab_6_TDM")

#Load CSV Data
data_g5_01 <- read.csv("./data/G5_001.csv", header = TRUE)
data_g5_02 <- read.csv("./data/G5_002.csv", header = TRUE)

#Compare data length
g5_01_length <- length(data_g5_01$PAM)
g5_02_length <- length(data_g5_02$PAM)
print(g5_01_length == g5_02_length)

#G5-01 tiene 1494 muestras y G5-02 tiene 1495
#Se elimina la ultima fila de G5-02
data_g5_02 <- data_g5_02[-g5_02_length,]


#Tiempo muestreo
Ts <- 0.2
Tiempo <- seq(Ts, length(data_g5_01$VFSC)*Ts, Ts)

#Plot G5-001 - VFSC
ggplot(data_g5_01, aes(x=Tiempo, y=VFSC)) +
  geom_line(color="#60c70c", size=0.7, alpha=1, linetype=1) +
  ggtitle("Time vs VFSC - G5-001") + 
  labs(y= "VFSC (ml/sec)", x = "Time (s)")

#Plot G5-001 - PAM
ggplot(data_g5_01, aes(x=Tiempo, y=PAM)) +
  geom_line(color="#fc7303", size=0.7, alpha=1, linetype=1) +
  ggtitle("Time vs PAM - G5-001") + 
  labs(y= "PAM (mmHg)", x = "Time (s)")

#Plot G5-002 - VFSC
ggplot(data_g5_02, aes(x=Tiempo, y=VFSC)) +
  geom_line(color="#60c70c", size=0.7, alpha=1, linetype=1) +
  ggtitle("Time vs VFSC - G5-001") + 
  labs(y= "VFSC (ml/sec)", x = "Time (s)")

#Plot G5-002 - PAM
ggplot(data_g5_02, aes(x=Tiempo, y=PAM)) +
  geom_line(color="#fc7303", size=0.7, alpha=1, linetype=1) +
  ggtitle("Time vs PAM - G5-001") + 
  labs(y= "PAM (mmHg)", x = "Time (s)")

retardos_multi <- function(SignalData, lags) {
  signal.uni <- SignalData
  max.lag <- max(unlist(lags)) + 1
  indices <- 1:nrow(signal.uni)
  lag.mat <- embed(indices, max.lag)
  col.names <- list("PAMn", "VFSCn")
  columns <- NULL
  lagged.columns.names <- c()
  
  for(colname in col.names) {
    lag.order <- lags[[colname]]
    columns[[colname]] <- signal.uni[lag.mat[, 1], colname]
    
    if(!is.null(lag.order) && lag.order > 0) {
      for(i in 1:lag.order) {
        new.colname <- paste(colname, paste0("lag", i), sep=".")
        lagged.columns.names <- c(lagged.columns.names, new.colname)
        columns[[new.colname]] <- signal.uni[lag.mat[, i+1], colname]
      }
    }
  }
  folded.signal <- data.frame(columns)
  sorting <- order(lag.mat[, 1])
  folded.signal <- folded.signal[sorting, ]
  
  return=list(folded.signal = folded.signal, lagged.columns.names = lagged.columns.names)
}

#Generate parameters
grid_cost <- 2^seq(-4,12,1)
grid_nu <- seq(0.1, 0.9, 0.1)
grid_gamma <- 2^seq(-4, 12, 1)
grid_lag <- seq(1,5,1)

PAMn <- (data_g5_01$PAM - min(data_g5_01$PAM))/(max(data_g5_01$PAM) - min(data_g5_01$PAM))
VFSCn <- (data_g5_01$VFSC - min(data_g5_01$VFSC))/(max(data_g5_01$VFSC) - min(data_g5_01$VFSC))

PAMn_test <-(data_g5_02$PAM - min(data_g5_02$PAM))/(max(data_g5_02$PAM) - min(data_g5_02$PAM))
VFSCn_test <-(data_g5_02$VFSC - min(data_g5_02$VFSC))/(max(data_g5_02$VFSC) - min(data_g5_02$VFSC))

train_data <- data.frame(PAMn, VFSCn)
test_data <- data.frame(PAMn_test, VFSCn_test)
  
params <-expand.grid(lagsList = grid_lag, cost = grid_cost, nu = grid_nu, gamma = grid_gamma)

output <- (c(foreach(i = 1:nrow(params), combine = rbind, .inorder = FALSE)
             %dopar% {
               c <- params[i,]$cost
               n <- params[i,]$nu
               g <- params[i,]$gamma
               l <- params[i,]$lagsList
               
               lag <- list(PAMn=1, VFSCn=0)
               signal.train <- retardos_multi(train_data, lag)
               retDatos = signal.train$folded.signal
               
               x = subset(retDatos, select = -VFSCn)
               y = retDatos$VFSCn
               
               modelo <- svm(x, y, type = "nu-regression", kernel = "radial", 
                             cost = c, nu = n, gamma = g)
               
               
             }))
