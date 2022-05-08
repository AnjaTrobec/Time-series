library(dplyr)
library(tseries)
library(LSTS)
library(forecast)
library("TTR")


#______________________________________________________________________________
#DVOMI IN VPRAŠANJA
# Zakaj harmonična regresija ni grafično ustrezna?
# Kaj je lag, kako prepoznati parametre za MA, ARMA, AR
# Graf residualov, je ustrezen?
# Kako poiskati periodo pri HR?


rm(list=ls()) 
#______________________________________________________________________________
setwd("C:/Users/aanja/OneDrive/Dokumenti/fmf/magisterij/časovne vrste/seminar/Time-series")
# 1. Narišite graf in komentirajte, ali se iz njega vidi kakšen trend ali sezonskost.
file2 <- file("A09tri.txt")
str_data2 <- readLines(file2)
data2 <- as.numeric(unlist(strsplit(str_data2, " ")))
x <- ts(data2)

# vsa opazovanja:
plot(x, main="Time series 2", ylim = c(-100,2000), type='l')

# le prvih 100:
plot(x, main="Time series 2", xlim=c(0,100), ylim = c(-100,2000), type='l')

# ne opazimo trenda, opazimo sezonskost

#______________________________________________________________________________

# 2. Odstranite morebiten trend in sezonskost z metodami, uporabljenimi pri tečaju:
# (zaporedno) diferenciranje, logaritmiranje, neposredna ocena sezonskih komponent,
# polinomski trend stopnje največ 3 ali prileganje periodične funkcije (ali kakšna kombinacija teh metod).
# Potem ko odstranite morebiten trend, narišite tudi surovi in zglajeni periodogram ter komentirajte,
# ali se vidi kakšna sezonskost in kakšna naj bi bila perioda.

# REŠEVANJE:
# logaritmiranja ne moremo uporabiti, ker imamo negativne podatke
# diferenciranje:

d <- diff(x)
plot(d, xlim = c(0,100), ylab=expression(paste(nabla,x[t]," [pp]")),type='l')

abline(h=mean(d),col="blue")
abline(h=c(mean(d)-sd(d),mean(d)+sd(d)),col="blue",lty="dotted")

#Iz grafa ne razberemo trenda. Sezonskost se opazi, zato jo bomo
#v nadaljevanju odstranili. Periodo bomo ugotovili s periodogramom.

# periode <-  periodogram(x)
# indeks <- which(periode$periodogram == max(periode$periodogram))
# per_indeks <- periode$periodogram[indeks]
# perioda1 <- 1/per_indeks
# perioda2 <- 1/periode$periodogram[215]

periodogram(d)
perioda <-  periodogram(x)$lambda[which.max(periodogram(x)$periodogram)] 
perioda2 <- periodogram(x)$periodogram[208] #ni statistično značilna!
# perioda3 <- periodogram(x)$periodogram[206]
# perioda4 <- periodogram(x)$periodogram[209]
# perioda5 <- periodogram(x)$periodogram[205]  

t <- time(x)
mod.H <- lm(x ~  I(sin(perioda*t)) + I(cos(perioda*t)) + I(sin(perioda2*t)) + I(cos(perioda2*t)))
summary(mod.H)

x.fit <- mod.H$fitted.values
x.res <- mod.H$residuals

x.fit <- ts(x.fit,start=start(x),frequency=frequency(x))
x.res <- ts(x.res,start=start(x),frequency=frequency(x))

par(mfrow=c(1,1))
plot(x, type='l', xlim=c(0,100))
points(x.fit,col="red",type="o",pch=16,cex=0.5)

#PROBLEM- bolj gost x.fit bi potrebovali, dobiva špice, ker imava premalo točk, 
# kako zgostiti število točk?

################################################################################
# DRUG PRISTOP:
library(TSA)
x.p <- periodogram(x)
#we can see some seasonality

order(x.p$spec,decreasing=TRUE) #these are heights

# Top frequencies at 216 215 223 222
#what are those frequencies?

x.p$freq[c(216,215,223,222)]

# Raw periodogram - we take frequences
# sin(2*pi*0.2500000*t) + cos(2*pi*0.2500000*t)
# sin(2*pi*0.2488426*t) + cos(2*pi*0.2488426*t)
# sin(2*pi*0.2581019*t) + cos(2*pi*0.2581019*t)
# sin(2*pi*0.2569444*t) + cos(2*pi*0.2569444*t)

#what are the periods of those?
1/x.p$freq[c(216,215,223,222)]

# Smoothed periodogram
# ====================
# Raw periodogram as a density estimator (identical numbers)

spectrum(x,
         log="no",
         demean=TRUE,
         detrend=FALSE,
         taper=0,
         lwd = 2)
spectrum(x,
         log="no",
         demean=TRUE,
         detrend=FALSE,
         taper=0,
         lwd = 2,
         xlim=c(0.22,0.28))

# we have to smooth it! By using kernels!

# Smoothing kernels
kernel(coef='daniell',m=2) #2 steps to the left and 2 steps to the right
kernel(coef='modified.daniell',m=1) #1 step left & 1 step right

# Iterating smoothing kernels
kernel(coef='daniell', m=c(1,1))
kernel(coef='modified.daniell', m=c(1,1))

# Visually choose kernel (Command spec.pgram does the same!)
spectrum(x,
         kernel=kernel(coef='modified.daniell', m=2),
         log="no",
         demean=TRUE,
         detrend=FALSE,
         taper=0,
         xlim=c(0.22,0.28))

spectrum(x,
         spans=7, # L = 2m+1 for Modified Daniell kernel (default); "rounds up" if even
         log="no",
         demean=TRUE,
         detrend=FALSE,
         taper=0,
         xlim=c(0.22,0.28))

x.sp <- spectrum(d,
                 kernel=kernel(coef='daniell', m=c(3,3)),
                 log="no",
                 demean=TRUE,
                 detrend=FALSE,
                 taper=0,
                 lwd = 1.5,
                 xlim=c(0.22,0.28))

order(x.sp$spec,decreasing = TRUE)

# Top frequencies at 216 217 215 218

#SM. Periodogram
x.sp$freq[c(216,217,215,218)]

# Period
1/x.sp$freq[c(216,217,215,218)]

# Harmonic regression
# ===================
library(TSA)
plot(x, ylab=expression(paste(nabla,x[t]," [pp]")))

t <- time(x)

mod.H1 <- lm(x ~ I(sin(2*pi*t*4)) + I(cos(2*pi*t*4)))
summary(mod.H1)

mod.H2 <- lm(x ~ I(sin(2*pi*t*4)) + I(cos(2*pi*t*4)) + I(sin(2*pi*t*3.981567)) + I(cos(2*pi*t*3.981567)))
summary(mod.H2)

# Fitted values and residuals

x.fit <- mod.H1$fitted.values
x.fit <- ts(x.fit,start=start(d),frequency=frequency(d))

points(x.fit,col="red",type="o",pch=16,cex=0.5)

r <- residuals(mod.H2)

#______________________________________________________________________________
## 3. TOČKA: Narišite graf rezidualov in komentirajte, ali so videti stacionarni. 
# Stacionarnost tudi preizkusite z uporabo ustreznih statističnih metod.

plot(x.res, main="Residuals")

summary(mod.H)$adj.r.squared

plot(x.res,
     ylab="Residuals",type='l')

abline(h=mean(x.res),col="blue")
abline(h=c(mean(x.res)-sd(x.res),mean(x.res)+sd(x.res)),col="blue",lty="dotted")


library(tseries)
# Augmented Dickey-Fuller Test za preverjanje stacionarnosti
adf.test(x.res,alternative="stationary")

# we say: ro - 1 = gamma 
# Ho : gamma = 0 -> random walk (not stationary)
# H1 : gamma < 0 -> stationary series

# this test is only appropriate if we remove seasonal and trend component
# result: we reject Ho, so series is stationary!

###############################################################################
## 4. TOČKA: Na rezidualih naredite grafikona ACF in PACF in na njuni podlagi predlagajte vsaj en model vrste AR(p) ali MA(q).

library(forecast)

#zelo slabo, to ne more biti okej?
Acf(x.res, main = "Acf residuals")  # model za to ?
Pacf(x.res, main = "Pacf residuals") # AR(3)?

diferen <- diff(x.res, lag = 1)
Pacf(diferen) 
Acf(diferen)

#PREDLOG MODELA:

###############################################################################
## 5. TOČKA
#Na podlagi Yule–Walkerjevih cenilk in kriterija AIC izberite najboljši model AR(p). 
# Primerjajte ga z najboljšim modelom ARMA(p, q) za p + q ≤ 3 po kriteriju AIC 
# (pozor: kriterij AIC je lahko definiran drugače od postopka do postopka). 
# Če je videti smiselno, pa namesto tega uporabite model GARCH.
# ==========
# ARMA model
# ==========

ar(x.res, arg = "yule–walker", aic=TRUE)

# my.arma <- stats::arima(r, order=c(1,0,1), include.mean = FALSE, method ='ML')
# 
# AIC(my.arma)
# 
# ic.values <- matrix(NA,nrow=5,ncol=5,
#                     dimnames=list(paste("AR=",0:4,sep=""),
#                                   paste("MA=",0:4,sep="")))
# 
# 
# for(i in 1:5){
#   for(j in 1:5){
#     my.arma <- stats::arima(r, order=c(i-1,0,j-1), include.mean = FALSE, method ='ML')
#     ic.values[i,j] <- AIC(my.arma)
#     remove(my.arma)
#   }
# }
# 
# which.min(ic.values) #or
# ic.values == min(ic.values)
# 
# #MA - moving average
# #AR - autoregressive function
# 
# # The optimal model is AR(1) or ARMA(1,0)
# 
# my.arma <- stats::arima(r, order=c(1,0,0), include.mean = FALSE, method ='ML')
# 
# confint(my.arma, level = 0.95)




ar(x.res, arg = "yule–walker", aic=TRUE)

#Primerjajte ga z najboljšim modelom ARMA(p, q)
#aic želimo minimizirati
arima(x.res, order=c(1,0,1)) #aic = 836.11 - to je arma (1,1)
arima(x.res, order=c(1,0,0))
arima(x.res, order=c(0,0,1)) #877.33
arima(x.res, order=c(2,0,1))
arima(x.res, order=c(1,0,2))
arima(x.res, order=c(2,0,0))
arima(x.res, order=c(0,0,2))
arima(x.res, order=c(3,0,0))
arima(x.res, order=c(0,0,3))
arima(x.res, order=c(4,0,0))


#Če je videti smiselno, pa namesto tega uporabite model GARCH.

AIC(garch(x.res, order = c(0,1))) #814,54
AIC(garch(x.res, order = c(1,1)))
AIC(garch(x.res, order = c(2,1)))
AIC(garch(x.res, order = c(1,0)))
AIC(garch(x.res, order = c(1,2)))
AIC(garch(x.res, order = c(0,3)))
AIC(garch(x.res, order = c(3,0)))


###############################################################################
## 6. TOČKA
# Izberite »optimalni« model in ocenite vse njegove parametre.
#best <- garch(diff(x.res, 1), order = c(1,1)) 
best <- arima(x.res, order=c(1,0,1)) #aic = 836.11
plot(best$residuals)
hist(best$residuals)
shapiro.test(best$residuals)
xx <- seq(-4, 4, by=0.001)
y1 <- dnorm(xx, mean=mean(best$residuals), sd=(sd(best$residuals)))
plot(density(best$residuals),col="blue", main ="Normal density funcion comparison", xlab="x")
lines(xx,y1, col="red")
legend("topleft",
       c("Residual density","Normal function"),
       col=c("blue","red"),
       lty="solid",
       bty="n")
# It does not look normal

###############################################################################
## 7. TOČKA
#Oglejte si ostanke po vašem modelu in komentirajte, ali so videti kot beli šum. 
#Primerjajte njihovo porazdelitev z normalno.

#odgovori
#Če je white noise nobena frekvenca ni dominantna. Bi mogel biti graf dokaj raven

plot(periodogram(best$residuals)$lambda, periodogram(best$residuals)$periodogram, type="l")


Box.test(x.res, lag = 1, type="Box-Pierce")
Box.test(x.res, lag = 1, type="Ljung-Box")

#p-value = 1.776e-15

#The P-Value of the Ljung-Box white noise test 
#is greater than significance level (i.e. α), so 
#we don't reject the white noise hypothesis (Ho), 
#or, simply stated; there is no statistical evidence of a serial 
#correlation, so the data can be white noise.

#Is not white noise.

###############################################################################
## 8. TOČKA
#Z uporabo izbranega modela in pod predpostavko normalnosti z R-ovo funkcijo predict 
#konstruirajte 90% napovedni interval za naslednjo vrednost. Ne pozabite vračunati tudi 
#odstranjenega trenda in sezonskosti.

prediction <- predict(best, n.ahead = 1, interval = "prediction", level = 0.90)
val1 <- prediction$pred[1]
sd1 <- prediction$se[1]
z <- 1.64

n <- 261
napoved <- 3.17670 + 2.92612*sin(perioda*n) + 2.11376 * cos(perioda*n)
naslednji <- exp(napoved + val1)
optimal_interval <- c(exp(napoved + val1 - z*sd1), exp(napoved + val1+z*sd1))
plot(append(naslednji, ts2)[210:261], type="l", main="Prediction", ylab="Value")

###############################################################################
## 9. TOČKA
# Dobljeni napovedni interval primerjajte z napovednim intervalom, 
#ki bi ga dobili, če bi naivno privzeli, da so podatki kar Gaussov 
#beli šum – pred in po odstranitvi trenda in sezonskosti.
#http://jse.amstat.org/v13n1/olsson.html


naive_gaussian <- mean(ts2)
naive_gaussian_interval <- c(naive_gaussian - z*sd(ts2), naive_gaussian + z*sd(ts2))
plot(append(ts2, naive_gaussian), type="l")

better_gaussian <- exp(napoved + mean(x.res) + (sd(x.res)^2)/2)

spodnja_meja <- exp(napoved + mean(x.res) + (sd(x.res)^2)/2 - z * sqrt((sd(x.res))^2/length(x.res) + (sd(x.res)^4)/(2*(length(x.res)-1))))
zgornja_meja <- exp(napoved + mean(x.res) + (sd(x.res)^2)/2 + z * sqrt((sd(x.res))^2/length(x.res) + (sd(x.res)^4)/(2*(length(x.res)-1))))
better_gaussian_interval <- c(spodnja_meja, zgornja_meja)

plot(append(ts2, better_gaussian), type="l")

plot(append(ts2, naslednji), type="l")
plot(append(ts2, naive_gaussian), type="l")
plot(append(ts2, better_gaussian), type="l")

