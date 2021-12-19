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

# #Plot G5-001 - VFSC
# ggplot(data_g5_01, aes(x=Tiempo, y=VFSC)) +
#   geom_line(color="#60c70c", size=0.7, alpha=1, linetype=1) +
#   ggtitle("Time vs VFSC - G5-001") + 
#   labs(y= "VFSC (ml/sec)", x = "Time (s)")
# 
# #Plot G5-001 - PAM
# ggplot(data_g5_01, aes(x=Tiempo, y=PAM)) +
#   geom_line(color="#fc7303", size=0.7, alpha=1, linetype=1) +
#   ggtitle("Time vs PAM - G5-001") + 
#   labs(y= "PAM (mmHg)", x = "Time (s)")
# 
# #Plot G5-002 - VFSC
# ggplot(data_g5_02, aes(x=Tiempo, y=VFSC)) +
#   geom_line(color="#60c70c", size=0.7, alpha=1, linetype=1) +
#   ggtitle("Time vs VFSC - G5-001") + 
#   labs(y= "VFSC (ml/sec)", x = "Time (s)")
# 
# #Plot G5-002 - PAM
# ggplot(data_g5_02, aes(x=Tiempo, y=PAM)) +
#   geom_line(color="#fc7303", size=0.7, alpha=1, linetype=1) +
#   ggtitle("Time vs PAM - G5-001") + 
#   labs(y= "PAM (mmHg)", x = "Time (s)")

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

getModels <- function(train, test, params) {
  start_time <- Sys.time()
  output_par <- (c(foreach(i = 1:nrow(params), combine = rbind, .inorder = FALSE)
                   %dopar% {
                     c <- params[i,]$cost
                     n <- params[i,]$nu
                     g <- params[i,]$gamma
                     l <- params[i,]$lagsList
                     
                     lag <- list(PAMn=1, VFSCn=0)
                     signal.train <- retardos_multi(train, lag)
                     retDatos = signal.train$folded.signal
                     
                     x = subset(retDatos, select = -VFSCn)
                     y = retDatos$VFSCn
                     
                     modelo <- svm(x, y, type = "nu-regression", kernel = "radial", 
                                   cost = c, nu = n, gamma = g)
                     
                     
                     signal.test <- retardos_multi(test, lag)
                     retDatos.test = signal.test$folded.signal
                     
                     x_test = subset(retDatos.test, select = -VFSCn)
                     y_test = retDatos.test$VFSCn
                     
                     pred <- predict(modelo, x_test)
                     corr_pred <- cor(pred, y_test, method = "pearson")
                     
                     c(l, c, n, g, corr_pred)
                   }))
  
  end_time <- Sys.time()
  print(end_time - start_time)
  
  output <- matrix(unlist(output_par), ncol = 5, byrow = TRUE)
  best_models <- output[order(output[, 5], decreasing = TRUE),]
  return(best_models)
}

#Generate parameters
grid_cost <- 2^seq(-4,12,1)
grid_nu <- seq(0.1, 0.9, 0.1)
grid_gamma <- 2^seq(-4, 12, 1)
grid_lag <- seq(1,5,1)

#=============================================
#Models generated using G5_01
#=============================================
PAMn <- (data_g5_01$PAM - min(data_g5_01$PAM))/(max(data_g5_01$PAM) - min(data_g5_01$PAM))
VFSCn <- (data_g5_01$VFSC - min(data_g5_01$VFSC))/(max(data_g5_01$VFSC) - min(data_g5_01$VFSC))

normalized_data <- data.frame(PAMn, VFSCn)

ind <- sample(2, nrow(normalized_data), replace = TRUE,prob = c(0.5, 0.5))

train_data_A <- normalized_data[ind==1,]
train_data_B <- normalized_data[ind==2,]
  
test_data_A <- normalized_data[ind==2,]
test_data_B <- normalized_data[ind==1,]

params <-expand.grid(lagsList = grid_lag, cost = grid_cost, nu = grid_nu, gamma = grid_gamma)

#G5_01 -> A (train) - B (test)
best_models_g5_01_A <- getModels(train_data_A, test_data_A, params)
df_g5_01_A <- as.data.frame(best_models_g5_01_A)
colnames(df_g5_01_A) <- c("lag", "cost", "nu", "gamma", "corr_pred")
df_g5_01_A <- df_g5_01_A[order(df_g5_01_A$corr_pred), ]

write.csv(df_g5_01_A, "./output/best_models_g5_01_A.csv", row.names = FALSE)

#G5_01 -> B (train) - A (test)
best_models_g5_01_B <- getModels(train_data_B, test_data_B, params)
df_g5_01_B <- as.data.frame(best_models_g5_01_B)
colnames(df_g5_01_B) <- c("lag", "cost", "nu", "gamma", "corr_pred")
df_g5_01_B <- df_g5_01_B[order(df_g5_01_B$corr_pred), ]

write.csv(df_g5_01_B, "./output/best_models_g5_01_B.csv", row.names = FALSE)


#=============================================
#Models generated using G5_02
#=============================================
PAMn_g5_02 <- (data_g5_02$PAM - min(data_g5_02$PAM))/(max(data_g5_02$PAM) - min(data_g5_02$PAM))
VFSCn_g5_02 <- (data_g5_02$VFSC - min(data_g5_02$VFSC))/(max(data_g5_02$VFSC) - min(data_g5_02$VFSC))

normalized_data_g5_02 <- data.frame(PAMn, VFSCn)

ind_g5_02 <- sample(2, nrow(normalized_data_g5_02), replace = TRUE,prob = c(0.5, 0.5))

train_data_A_g5_02 <- normalized_data_g5_02[ind_g5_02==1,]
train_data_B_g5_02 <- normalized_data_g5_02[ind_g5_02==2,]

test_data_A_g5_02 <- normalized_data[ind==2,]
test_data_B_g5_02 <- normalized_data[ind==1,]

#G5_02 -> A (train) - B (test)
best_models_g5_02_A <- getModels(train_data_A_g5_02, test_data_A_g5_02, params)
df_g5_02_A <- as.data.frame(best_models_g5_02_A)
colnames(df_g5_02_A) <- c("lag", "cost", "nu", "gamma", "corr_pred")
df_g5_02_A <- df_g5_02_A[order(df_g5_02_A$corr_pred), ]

write.csv(df_g5_02_A, "./output/best_models_g5_02_A.csv", row.names = FALSE)

#G5_02 -> B (train) - A (test)
best_models_g5_02_B <- getModels(train_data_B_g5_02, test_data_B_g5_02, params)
df_g5_02_B <- as.data.frame(best_models_g5_02_B)
colnames(df_g5_02_B) <- c("lag", "cost", "nu", "gamma", "corr_pred")
df_g5_02_B <- df_g5_02_B[order(df_g5_02_B$corr_pred), ]

write.csv(df_g5_02_B, "./output/best_models_g5_02_B.csv", row.names = FALSE)

