library(forecast)

#DEFINE DATASET
dataset_y <- c(1,2,3,4,3,2,1,2,3,4,3,2,1)
print(dataset_y)

#FIT MODEL
num_predictions <- 20
algo <- Arima(dataset_y, order=c(1,0,0))

#PREDICT NEW VALUES USING ARIMA
fc <- forecast(algo,  h=num_predictions)

#Figure stuff
#plot(fc) 

###########################
# MANUAL AR MODEL
#############################

#FIT AR MODEL; ar.burg or ar.yw or ar.mle
fitted_ar <- ar.burg(dataset_y, aic = FALSE, order.max = 1, method = "ols")

#PRINT COEFFICIENT
print(fitted_ar$ar[,,1])

#FORECAST
fc <- forecast(fitted_ar,  h=num_predictions)

#PLOT
plot(fc) 