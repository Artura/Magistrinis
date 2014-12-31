rm(list = ls())
library(ggplot2)
library(zoo)
library(plyr)
library(sampling)
library(R2jags)
library(ggmcmc)
library(forecast)
library(abind)
library(xtable)
  
data <- read.csv2("duomenys.csv", header = TRUE, sep = ";", dec = ".") 
attach(data)

#----------------------------------Duomenø nuskaitymas-------------------------------------------

y1 <- c(y10_1, y10_2, y10_3, y10_4, y10_5, y10_6, y10_7, y10_8, y10_9, y10_10, y10_11, y10_12,
       y11_1, y11_2, y11_3, y11_4, y11_5, y11_6, y11_7, y11_8, y11_9, y11_10, y11_11, y11_12,
       y12_1, y12_2, y12_3, y12_4, y12_5, y12_6, y12_7, y12_8, y12_9, y12_10, y12_11, y12_12,
       y13_1, y13_2, y13_3, y13_4, y13_5, y13_6, y13_7, y13_8, y13_9, y13_10, y13_11, y13_12)
  
x1 <- c(x10_1, x10_2, x10_3, x10_4, x10_5, x10_6, x10_7, x10_8, x10_9, x10_10, x10_11, x10_12,
       x11_1, x11_2, x11_3, x11_4, x11_5, x11_6, x11_7, x11_8, x11_9, x11_10, x11_11, x11_12,
       x12_1, x12_2, x12_3, x12_4, x12_5, x12_6, x12_7, x12_8, x12_9, x12_10, x12_11, x12_12,
       x13_1, x13_2, x13_3, x13_4, x13_5, x13_6, x13_7, x13_8, x13_9, x13_10, x13_11, x13_12)

z1 <- c(z10_1, z10_2, z10_3, z10_4, z10_5, z10_6, z10_7, z10_8, z10_9, z10_10, z10_11, z10_12,
       z11_1, z11_2, z11_3, z11_4, z11_5, z11_6, z11_7, z11_8, z11_9, z11_10, z11_11, z11_12,
       z12_1, z12_2, z12_3, z12_4, z12_5, z12_6, z12_7, z12_8, z12_9, z12_10, z12_11, z12_12,
       z13_1, z13_2, z13_3, z13_4, z13_5, z13_6, z13_7, z13_8, z13_9, z13_10, z13_11, z13_12)

from <- as.Date("2010-01-01")
to <- as.Date("2013-12-01")
months <- seq.Date(from = from, to = to, by = "month")

menesiai <- c(rep(months[1], 3812), rep(months[2], 3812), rep(months[3], 3812), rep(months[4], 3812),
        rep(months[5], 3812), rep(months[6], 3812), rep(months[7], 3812), rep(months[8], 3812),
        rep(months[9], 3812), rep(months[10], 3812), rep(months[11], 3812), rep(months[12], 3812),
        rep(months[13], 3812), rep(months[14], 3812), rep(months[15], 3812), rep(months[16], 3812), 
        rep(months[17], 3812), rep(months[18], 3812), rep(months[19], 3812), rep(months[20], 3812),
        rep(months[21], 3812), rep(months[22], 3812), rep(months[23], 3812), rep(months[24], 3812),
        rep(months[25], 3812), rep(months[26], 3812), rep(months[27], 3812), rep(months[28], 3812),
        rep(months[29], 3812), rep(months[30], 3812), rep(months[31], 3812), rep(months[32], 3812),
        rep(months[33], 3812), rep(months[34], 3812), rep(months[35], 3812), rep(months[36], 3812),
        rep(months[37], 3812), rep(months[38], 3812), rep(months[39], 3812), rep(months[40], 3812),
        rep(months[41], 3812), rep(months[42], 3812), rep(months[43], 3812), rep(months[44], 3812),
        rep(months[45], 3812), rep(months[46], 3812), rep(months[47], 3812), rep(months[48], 3812))

numeris <- rep(nr, 48)
veikla <- rep(veikla2, 48)
apskritis <- rep(apskr, 48)
vieta <- rep(vietove, 48)
dataNA <- data.frame(menesiai, numeris, veikla, apskritis, vieta, y = y1, x = x1, z = z1)
dataG <- na.omit(dataNA)
dataG <- dataG[order(dataG$numeris),]

month0 <- as.factor(dataG$menesiai)
nb0 <- as.factor(dataG$numeris)
dataG$mixed <- month0:nb0

dataG$logy <- log(dataG$y)
dataG$logx <- log(dataG$x)
dataG$logz <- log(dataG$z)
is.na(dataG) <- sapply(dataG, is.infinite) 
dataGo <- na.omit(dataG)

#-------------------------------Isskirèiø tyrimas---------------------------------------

veikla45 <- dataGo[dataGo$veikla %in% sort(unique(dataGo$veikla))[1:6],]
veikla46 <- dataGo[dataGo$veikla %in% sort(unique(dataGo$veikla))[7:54],]
veikla47 <- dataGo[dataGo$veikla %in% sort(unique(dataGo$veikla))[55:90],]

veikla45B <- veikla45[veikla45$logx > 15.13,]
veikla46B <- veikla46[veikla46$logx > 15.699,]
veikla47B <- veikla47[veikla47$logx > 14.35,]
firmBAD <- c(unique(veikla45B$numeris), unique(veikla46B$numeris),
  unique(veikla47B$numeris))
dataOut <- dataGo[ ! dataGo$numeris %in% firmBAD, ] 

#-----------------------------------Imtis------------------------------------------------

kiek <- rep(0, length(unique(dataOut$veikla)))
ve <- unique(dataOut$veikla)
for(i in 1:length(ve)){
  kiek[i] <- ceiling(length(unique(dataOut$numeris[dataOut$veikla == ve[i]]))*0.10)
}

prob <- split(rep(0, sum(kiek)), rep(1:length(kiek), times = kiek))

probs <- ddply(dataOut, .(veikla), function(x) {
  x$pr <- inclusionprobabilities(x$y, kiek[ve == x$veikla[1]])
  ddply(x, .(numeris), transform, pr = sum(pr))
})
prTmp <- unique(probs[, c("numeris", "pr", "veikla")])

set.seed(12320144)
s <- na.omit(unlist(sapply(split(prTmp$pr, prTmp$veikla), UPpivotal)))
dataFinal <- dataOut[dataOut$numeris %in% prTmp$numeris[s == 1], ]

dataPan <- dataFinal[c("menesiai", "numeris", "veikla", "logy", "logx", "logz")]

tbl <- table(dataPan$numeris, !is.na(dataPan$logy))
fns <- rownames(tbl)[tbl == 48]
ns <- table(dataPan$veikla[dataPan$numeris %in% fns]) / 48
vs <- dimnames(ns)[[1]]

w1 <- array(0, c(48, max(ns), length(ns)))
for(f in fns) {
  sset <- dataPan[dataPan$numeris == f, ]
  v <- sset$veikla[1]
  ord <- which(colSums(w1[, , vs == v]^2) == 0)[1]
  w1[, ord, vs == v] <- sset$logx
}

w2 <- array(0, c(48, max(ns), length(ns)))
for(f in fns) {
  sset <- dataPan[dataPan$numeris == f, ]
  v <- sset$veikla[1]
  ord <- which(colSums(w2[, , vs == v]^2) == 0)[1]
  w2[, ord, vs == v] <- sset$logz
}

yfull <- array(0, c(48, max(ns), length(ns)))
for(f in fns) {
  sset <- dataPan[dataPan$numeris == f, ]
  v <- sset$veikla[1]
  ord <- which(colSums(yfull[, , vs == v]^2) == 0)[1]
  yfull[, ord, vs == v] <- sset$logy
}

megarray <- array(NA, c(48, 48, 10, 71))
for(i in 1:71)
  for(f in 1:max(ns)) {
    sigma <- var(yfull[, f, i][-48])
    megarray[, , f, i] <- diag(ifelse(sigma > 0, sigma, 1), 48)
}

ko <- matrix(0, 10, 71)
for(i in 1:71)
  for(f in 1:max(ns)){
    sigma <- var(yfull[, f, i][-48])
    ko[f, i]<- ifelse(sigma > 0, sigma, 1)
}

y <- yfull[1:47, , ]  
y13 <- yfull[48, ,]
yNA <- y13
yNA[which(yNA!=0)] <- NA
ybugs <- abind(y, yNA, along = 1)

dataPanFull <- dataPan[dataPan$numeris %in% fns,]
i2 <- unique(dataPanFull$menesiai) 
j2 <- vs

kiek2 <- matrix(rep(0, 3408), nrow = 48)
for(i in 1:48){
  for(j in 1:71){
    kiek2[i, j] <- length(unique(dataPanFull$numeris[dataPanFull$menesiai == i2[i] & 
  dataPanFull$veikla == j2[j]]))
  }
}

#-------------------------------------Fay-Herriot modelis--------------------------------------------

par11 <- c(paste("y[48,", 1, ",", 1:71, "]", sep = ""), paste("y[48,", 2, ",", 1:71, "]", sep = ""),
  paste("y[48,", 3, ",", 1:71, "]", sep = ""), paste("y[48,", 4, ",", 1:71, "]", sep = ""),
  paste("y[48,", 5, ",", 1:71, "]", sep = ""), paste("y[48,", 6, ",", 1:71, "]", sep = ""),
  paste("y[48,", 7, ",", 1:71, "]", sep = ""), paste("y[48,", 8, ",", 1:71, "]", sep = ""),
  paste("y[48,", 9, ",", 1:71, "]", sep = ""), paste("y[48,", 10, ",", 1:71, "]", sep = ""), 
  "beta1", "beta2", "tauv", "tau", "tauu", "rho")
ns2 <- which(ns != max(ns))
deita1 <- list(y = ybugs, z1 = w1, z2 = w2, n = as.numeric(ns), ns2 = ns2, l = length(ns2),
  kova = ko)
inits1 <- function(){
   list(beta1 = 0.97, beta2 = 0.02, sigmav = 1, sigmau = 1, sigma = 1, rho = 0.2)
}

system.time({outj13 <- jags(data = deita1, inits = inits1, parameters.to.save = par11,
  model.file = "C:/Users/Artura/Desktop/FayHerriot.txt", 
  n.chains = 3, n.iter = 100000, n.thin = 20, DIC = TRUE) 
})/60

#----------------------------------------MCMC-------------------------------------------------------

print(outj13, dig = 3)
which(outj13$BUGSoutput$summary[, 8] > 1.1) 
outj13$BUGSoutput$summary[1:7, 8]

summary(outj13$BUGSoutput$sims.list$beta1)
sd(outj13$BUGSoutput$sims.list$beta1)
summary(outj13$BUGSoutput$sims.list$beta2)
sd(outj13$BUGSoutput$sims.list$beta2)
summary(outj13$BUGSoutput$sims.list$rho)
sd(outj13$BUGSoutput$sims.list$rho)
summary(outj13$BUGSoutput$sims.list$tau)
sd(outj13$BUGSoutput$sims.list$tau)
summary(outj13$BUGSoutput$sims.list$tauu)
sd(outj13$BUGSoutput$sims.list$tauu)
summary(outj13$BUGSoutput$sims.list$tauv)
sd(outj13$BUGSoutput$sims.list$tauv)

outj11.mcmc <- as.mcmc(outj13) 
param11 <- c("beta1", "beta2", "tauv", "tau", "tauu", "rho")
mcmc11 <-  as.mcmc(list(outj11.mcmc[[1]][, param11], outj11.mcmc[[2]][,param11], 
  outj11.mcmc[[3]][,param11]))
# gelman.plot(mcmc11) 
# geweke.plot(mcmc11)
# raftery.diag(mcmc11) 
heidel.diag(mcmc11)

# ggmcmc(gg11)
P <- data.frame(Parameter = c("beta1", "beta2", "tauv", "tau", "tauu", "rho"),
  Label = c("Beta1", "Beta2", "Sigmav", "Sigma", "Sigmau", "Rho"))
gg11 <- ggs(mcmc11, par_labels = P) 

# ggh1 <- ggs_histogram(gg11) + theme_bw() + xlab("Reikðmë") + ylab("Daþnis")
ggd1 <- ggs_density(gg11) + theme_bw() + xlab("Reikðmë") + ylab("Daþnis") +
  scale_fill_discrete("Grandinë") + scale_colour_discrete("Grandinë")
ggtr1 <- ggs_traceplot(gg11) + theme_bw() + xlab("Iteracija") + ylab("Reikðmë") +
  scale_fill_discrete("Grandinë") + scale_colour_discrete("Grandinë")
ggr1 <- ggs_running(gg11) + theme_bw() + xlab("Iteracija") + ylab("Vidurkis") +
  scale_fill_discrete("Grandinë") + scale_colour_discrete("Grandinë") 
# ggs_compare_partial(gg11)
gga1 <- ggs_autocorrelation(gg11) + theme_bw() + xlab("Ankstinys") + ylab("Autokoreliacija") +
  scale_fill_discrete("Grandinë") + scale_colour_discrete("Grandinë") 
ggcor1 <- ggs_crosscorrelation(gg11, absolute_scale = TRUE) + theme_bw() + 
  theme(legend.title = element_blank())
# ggs_Rhat(gg11) 
gggew1 <- ggs_geweke(gg11) + theme_bw() + xlab("z-reikðmë") + ylab("Parametras") +
  labs(title = "Geweke diagnostika") + scale_fill_discrete("Grandinë") + 
  scale_colour_discrete("Grandinë") 
# gghpd1 <- ggs_caterpillar(gg11, horizontal = FALSE) + theme_bw() + ylab("Parametras")

#-----------------------------------Prognozë----------------------------------------------

progn1 <- outj13$BUGSoutput$sims.list$y 
mr1 <- array(c(progn1), c(7500, 71, 10))
dimnames(outj13$BUGSoutput$sims.matrix)[[2]][1:7]
proglong1 <- outj13$BUGSoutput$sims.matrix[, -1:-7]
paror1 <- colnames(proglong1)
for(i in 1:(10 * 71)) {
  indF1 <- as.numeric(gsub(".*,(.*),.*", "\\1", paror1[i]))
  indI1 <- as.numeric(gsub(".*,(.*)]", "\\1", paror1[i]))
  mr1[, indI1, indF1] <- proglong1[, i]
}

progn1tr <- aperm(mr1, c(2:3, 1))
vidurpr1 <- apply(progn1tr, 1:2, mean) 
dfpr1 <- adply(vidurpr1, 2:1, identity)
y13_tikri <- adply(y13, 1:2, identity)
dfpr1$y13 <- y13_tikri$V1

w113 <- w1[48, ,]
w113_tikri <- adply(w113, 1:2, identity)
dfpr1$w113 <- w113_tikri$V1

names(dfpr1) <- c("firmos", "veikla", "prognoze", "tikros", "pvm")

dataPan13 <- dataPan[grep("2013-12", dataPan$menesiai), ] 
data2013 <- dataPan13[dataPan13$numeris %in% fns,] 
ordered <- data2013[ with( data2013, order(veikla)) , ]

dataPro1 <- dfpr1[apply(dfpr1[,3:4], 1, function(row) all(row !=0 )),]

dataProgn1 <- data.frame(numeris = ordered$numeris,
  veikla = ordered$veikla, logy = ordered$logy, prognoze = dataPro1$prognoze,
  z = exp(ordered$logz))

dataProgn1$skirtumas <- dataProgn1$logy - dataProgn1$prognoze

#--------------------Asimetriniø Normaliøjø pajamø modelis----------------------------

deita2 <- list(y = ybugs, z1 = w1, z2 = w2, n = as.numeric(ns), ns2 = ns2, l = length(ns2),
  kova = ko, kiek2 = kiek2)

par21 <- c(paste("y[48,", 1, ",", 1:71, "]", sep = ""), paste("y[48,", 2, ",", 1:71, "]", sep = ""),
  paste("y[48,", 3, ",", 1:71, "]", sep = ""), paste("y[48,", 4, ",", 1:71, "]", sep = ""),
  paste("y[48,", 5, ",", 1:71, "]", sep = ""), paste("y[48,", 6, ",", 1:71, "]", sep = ""),
  paste("y[48,", 7, ",", 1:71, "]", sep = ""), paste("y[48,", 8, ",", 1:71, "]", sep = ""),
  paste("y[48,", 9, ",", 1:71, "]", sep = ""), paste("y[48,", 10, ",", 1:71, "]", sep = ""), 
  "beta1", "beta2", "tauv", "tau", "tauu", "rho", "lambda") 
inits20 <- function(){
   list(beta1 = 0.97, beta2 = 0.018, sigma = 28.8, sigmau = 23.26, sigmav = 25, rho = 0.22, 
     lambda = 5)
}
system.time({outj20y <- jags(data = deita2, inits = inits20, parameters.to.save = par21,
  model.file = "C:/Users/Artura/Desktop/SkewNormalY.txt", 
  n.chains = 3, n.iter = 100000, n.thin = 20, DIC = TRUE) 
})/60

#----------------------------------------MCMC-------------------------------------------------------

print(outj20y, dig = 3)
which(outj20y$BUGSoutput$summary[, 8] > 1.1) 
outj20y$BUGSoutput$summary[1:8, 8]

summary(outj20y$BUGSoutput$sims.list$beta1)
sd(outj20y$BUGSoutput$sims.list$beta1)
summary(outj20y$BUGSoutput$sims.list$beta2)
sd(outj20y$BUGSoutput$sims.list$beta2)
summary(outj20y$BUGSoutput$sims.list$rho)
sd(outj20y$BUGSoutput$sims.list$rho)
summary(outj20y$BUGSoutput$sims.list$tau)
sd(outj20y$BUGSoutput$sims.list$tau)
summary(outj20y$BUGSoutput$sims.list$tauu)
sd(outj20y$BUGSoutput$sims.list$tauu)
summary(outj20y$BUGSoutput$sims.list$tauv)
sd(outj20y$BUGSoutput$sims.list$tauv)
summary(outj20y$BUGSoutput$sims.list$lambda)
sd(outj20y$BUGSoutput$sims.list$lambda)

outj20y.mcmc <- as.mcmc(outj20y) 
param20 <- c("beta1", "beta2", "tauv", "tau", "tauu", "rho", "lambda")
mcmc20 <-  as.mcmc(list(outj20y.mcmc[[1]][, param20], outj20y.mcmc[[2]][,param20], 
  outj20y.mcmc[[3]][,param20]))
gelman.plot(mcmc20) 
geweke.plot(mcmc20)
# raftery.diag(mcmc20)
heidel.diag(mcmc20) 

P2 <- data.frame(Parameter = c("beta1", "beta2", "tauv", "tau", "tauu", "rho", "lambda"),
  Label = c("Beta1", "Beta2", "Sigmav", "Sigma", "Sigmau", "Rho", "Lambda"))
gg20 <- ggs(mcmc20, par_labels = P2)

# ggmcmc(gg20)
# ggh2 <- ggs_histogram(gg20) + theme_bw() + xlab("Reikðmë") + ylab("Daþnis")
ggd2 <- ggs_density(gg20) + theme_bw() + xlab("Reikðmë") + ylab("Daþnis") +
  scale_fill_discrete("Grandinë") + scale_colour_discrete("Grandinë")
ggtr2 <- ggs_traceplot(gg20) + theme_bw() + xlab("Iteracija") + ylab("Reikðmë") +
  scale_fill_discrete("Grandinë") + scale_colour_discrete("Grandinë")
ggr2 <- ggs_running(gg20) + theme_bw() + xlab("Iteracija") + ylab("Vidurkis") +
  scale_fill_discrete("Grandinë") + scale_colour_discrete("Grandinë")  
# ggs_compare_partial(gg31)
gga2 <- ggs_autocorrelation(gg20) + theme_bw() + xlab("Ankstinys") + ylab("Autokoreliacija") +
  scale_fill_discrete("Grandinë") + scale_colour_discrete("Grandinë") 
ggcor2 <- ggs_crosscorrelation(gg20, absolute_scale = TRUE) + theme_bw() + 
  theme(legend.title = element_blank()) 
# ggs_Rhat(gg20) 
gggew2 <- ggs_geweke(gg20) + theme_bw() + xlab("z-reikðmë") + ylab("Parametras") +
  labs(title = "Geweke diagnostika") + scale_fill_discrete("Grandinë") + 
  scale_colour_discrete("Grandinë") 
gghpd2 <- ggs_caterpillar(gg20, horizontal = FALSE) + theme_bw() + ylab("Parametras")

#-----------------------------------Prognozë----------------------------------------------

progn2 <- outj20y$BUGSoutput$sims.list$y 
mr2 <- array(c(progn2), c(7500, 71, 10))
dimnames(outj20y$BUGSoutput$sims.matrix)[[2]][1:9]
proglong2 <- outj20y$BUGSoutput$sims.matrix[, -1:-8]
paror2 <- colnames(proglong2)
for(i in 1:(10 * 71)) {
  indF2 <- as.numeric(gsub(".*,(.*),.*", "\\1", paror2[i]))
  indI2 <- as.numeric(gsub(".*,(.*)]", "\\1", paror2[i]))
  mr2[, indI2, indF2] <- proglong2[, i]
}

progn2tr <- aperm(mr2, c(2:3, 1))
vidurpr2 <- apply(progn2tr, 1:2, mean) 
dfpr2 <- adply(vidurpr2, 2:1, identity) 
dfpr2$y13 <- y13_tikri$V1
dfpr2$w113 <- w113_tikri$V1
names(dfpr2) <- c("firmos", "veikla", "prognoze", "tikros", "pvm")

dataPro2 <- dfpr2[apply(dfpr2[,3:4], 1, function(row) all(row !=0 )),]

dataProgn2 <- data.frame(numeris = ordered$numeris,
  veikla = ordered$veikla, logy = ordered$logy, prognoze = dataPro2$prognoze,
  z = exp(ordered$logz), pvm = dataPro2$pvm)

dataProgn2$skirtumas <- dataProgn2$logy - dataProgn2$prognoze

#--------------------Asimetriniø Normaliøjø atsitiktiniø efektø modelis------------------------

deita3 <- list(y = ybugs, z1 = w1, z2 = w2, n = as.numeric(ns), ns2 = ns2, l = length(ns2),
  kova = ko)

par31 <- c(paste("y[48,", 1, ",", 1:71, "]", sep = ""), paste("y[48,", 2, ",", 1:71, "]", sep = ""),
  paste("y[48,", 3, ",", 1:71, "]", sep = ""), paste("y[48,", 4, ",", 1:71, "]", sep = ""),
  paste("y[48,", 5, ",", 1:71, "]", sep = ""), paste("y[48,", 6, ",", 1:71, "]", sep = ""),
  paste("y[48,", 7, ",", 1:71, "]", sep = ""), paste("y[48,", 8, ",", 1:71, "]", sep = ""),
  paste("y[48,", 9, ",", 1:71, "]", sep = ""), paste("y[48,", 10, ",", 1:71, "]", sep = ""), 
  "beta1", "beta2", "tau", "tauu", "rho", "lambda") 
par30 <- c(paste("y[48,", 1, ",", 1:71, "]", sep = ""), paste("y[48,", 2, ",", 1:71, "]", sep = ""),
  paste("y[48,", 3, ",", 1:71, "]", sep = ""), paste("y[48,", 4, ",", 1:71, "]", sep = ""),
  paste("y[48,", 5, ",", 1:71, "]", sep = ""), paste("y[48,", 6, ",", 1:71, "]", sep = ""),
  paste("y[48,", 7, ",", 1:71, "]", sep = ""), paste("y[48,", 8, ",", 1:71, "]", sep = ""),
  paste("y[48,", 9, ",", 1:71, "]", sep = ""), paste("y[48,", 10, ",", 1:71, "]", sep = ""), 
  "beta1", "beta2", "tau", "tauu", "rho", "lambda", "sigmav") 

inits32 <- function(){
   list(beta1 = 0.97, beta2 = 0.018, sigma = 28.8, sigmau = 23.26, rho = 0.22, lambda = 14.6)
}

system.time({outj32y <- jags(data = deita3, inits = inits32, parameters.to.save = par31,
  model.file = "C:/Users/Artura/Desktop/SkewNormalV.txt", 
  n.chains = 3, n.iter = 100000, n.thin = 20, DIC = TRUE) 
})/60 

system.time({outj30y <- jags(data = deita3, inits = inits32, parameters.to.save = par30,
  model.file = "C:/Users/Artura/Desktop/SkewNormalV.txt", 
  n.chains = 3, n.iter = 100000, n.thin = 20, DIC = TRUE) 
})/60 

#----------------------------------------MCMC-------------------------------------------------------

print(outj32y, dig = 3)
which(outj32y$BUGSoutput$summary[, 8] > 1.1) 
outj30y$BUGSoutput$summary[, 8]

summary(outj32y$BUGSoutput$sims.list$beta1)
sd(outj32y$BUGSoutput$sims.list$beta1)
summary(outj32y$BUGSoutput$sims.list$beta2)
sd(outj32y$BUGSoutput$sims.list$beta2)
summary(outj32y$BUGSoutput$sims.list$rho)
sd(outj32y$BUGSoutput$sims.list$rho)
summary(outj32y$BUGSoutput$sims.list$tau)
sd(outj32y$BUGSoutput$sims.list$tau)
summary(outj32y$BUGSoutput$sims.list$tauu)
sd(outj32y$BUGSoutput$sims.list$tauu)
summary(outj32y$BUGSoutput$sims.list$lambda)
sd(outj32y$BUGSoutput$sims.list$lambda)

sdv <- rep(0, 71)
minv <- rep(0, 71)
qu1v <- rep(0, 71)
medv <- rep(0, 71)
mev <- rep(0, 71)
qu3v <- rep(0, 71)
maxv <- rep(0, 71)
for( i in 1:71 ){
  sdv[i] <- sd(outj30y$BUGSoutput$sims.list$sigmav[,i])
  minv[i] <- summary(outj30y$BUGSoutput$sims.list$sigmav[,i])[1]
  qu1v[i] <- summary(outj30y$BUGSoutput$sims.list$sigmav[,i])[2]
  medv[i] <- summary(outj30y$BUGSoutput$sims.list$sigmav[,i])[3]
  mev[i] <- summary(outj30y$BUGSoutput$sims.list$sigmav[,i])[4]
  qu3v[i] <- summary(outj30y$BUGSoutput$sims.list$sigmav[,i])[5]
  maxv[i] <- summary(outj30y$BUGSoutput$sims.list$sigmav[,i])[6]
}

rhat <- outj30y$BUGSoutput$summary[6:76, 8]

svdf <- data.frame(summary = c(sdv, minv, qu1v, medv, mev, qu3v), indeksas = rep(1:71, 6),
  Statistika = c(rep("Standartinis nuokrypis", 71), rep("Minimumas", 71),
    rep("1 kvartilis", 71), rep("Mediana", 71),
    rep("Vidurkis", 71), rep("3 kvartilis", 71)))
ggplot(data = svdf, aes(x = indeksas, y = summary, color = Statistika)) +
  geom_line() + xlab("Veikla") + ylab("") + theme_bw() + 
  scale_x_continuous(breaks = c(1:71), labels = c(as.character(heiddf$veikla))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

svdf2 <- data.frame(rhat, maxv, indeksas = 1:71)
ggplot(data = svdf2, aes(x = indeksas, y = maxv)) +
  geom_line() + xlab("Veikla") + ylab("Maksimumas") + theme_bw() + 
  scale_x_continuous(breaks = c(1:71), labels = c(as.character(heiddf$veikla))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  
ggplot(data = svdf2, aes(x = indeksas, y = rhat)) +
  geom_line() + xlab("Veikla") + ylab("Gelman ir Rubin konvergavimo statistika") + theme_bw() + 
  scale_x_continuous(breaks = c(1:71), labels = c(as.character(heiddf$veikla))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))  

densv <- data.frame(sigv = as.vector(outj30y$BUGSoutput$sims.list$sigmav),
  veikla = factor(sort(rep(1:71, 7500))))
ggplot(densv, aes(x = sigv, fill = veikla)) + geom_density(alpha = .3) + theme_bw() +
  ylab("Tankis") + xlab("Reikðmë") + scale_x_continuous(limits = c(0, 10)) + 
  theme(legend.position = "none")  

outj31.mcmc <- as.mcmc(outj32y) 
param31 <- c("beta1", "beta2", "tau", "tauu", "rho", "lambda")
mcmc31 <-  as.mcmc(list(outj31.mcmc[[1]][, param31], outj31.mcmc[[2]][,param31], 
  outj31.mcmc[[3]][,param31]))
gelman.plot(mcmc31) 
geweke.plot(mcmc31)
# raftery.diag(mcmc31) 
heidel.diag(mcmc31) 

outj30.mcmc <- as.mcmc(outj30y)
param32 <- c(paste("sigmav[", 1:71, "]", sep = ""))
mcmc32 <-  as.mcmc(list(outj30.mcmc[[1]][, param32], outj30.mcmc[[2]][,param32], 
  outj30.mcmc[[3]][,param32]))

heidel <- heidel.diag(mcmc32)
heiddf <- data.frame(pvalue1 = heidel[[1]][, 3], pvalue2 = heidel[[2]][, 3], 
  pvalue3 = heidel[[3]][, 3], veikla = unique(ordered$veikla), indeksas = 1:71)
plotdf <- data.frame(x = heiddf$indeksas, y = unlist(heiddf[, 1:3]), veikla = heiddf$veikla,
  Grandinë = factor(rep(1:3, each = nrow(heiddf))))
ggplot(data = plotdf, aes(x = x, y = y, color = Grandinë)) +
  geom_line() + xlab("Veikla") + ylab("p-reikðmë") + theme_bw() + 
  geom_hline(aes(yintercept = 0.05)) +
  scale_x_continuous(breaks = c(1:71), labels = c(as.character(heiddf$veikla))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) 

P3 <- data.frame(Parameter = c("beta1", "beta2", "tau", "tauu", "rho", "lambda"),
  Label = c("Beta1", "Beta2", "Sigma", "Sigmau", "Rho", "Lambda"))
gg31 <- ggs(mcmc31, par_labels = P3)
# ggmcmc(gg31)
# ggh2 <- ggs_histogram(gg31) + theme_bw() + xlab("Reikðmë") + ylab("Daþnis")
ggd3 <- ggs_density(gg31) + theme_bw() + xlab("Reikðmë") + ylab("Daþnis") +
  scale_fill_discrete("Grandinë") + scale_colour_discrete("Grandinë")
ggtr3 <- ggs_traceplot(gg31) + theme_bw() + xlab("Iteracija") + ylab("Reikðmë") +
  scale_fill_discrete("Grandinë") + scale_colour_discrete("Grandinë")
ggr3 <- ggs_running(gg31) + theme_bw() + xlab("Iteracija") + ylab("Vidurkis") +
  scale_fill_discrete("Grandinë") + scale_colour_discrete("Grandinë")
# ggs_compare_partial(gg31)
gga3 <- ggs_autocorrelation(gg31) + theme_bw() + xlab("Ankstinys") + ylab("Autokoreliacija") +
  scale_fill_discrete("Grandinë") + scale_colour_discrete("Grandinë") 
ggcor3 <- ggs_crosscorrelation(gg31, absolute_scale = TRUE) + theme_bw() + 
  theme(legend.title = element_blank())
# ggs_Rhat(gg31) 
gggew3 <- ggs_geweke(gg31) + theme_bw() + xlab("z-reikðmë") + ylab("Parametras")  +
  labs(title = "Geweke diagnostika") + scale_fill_discrete("Grandinë") + 
  scale_colour_discrete("Grandinë") 
# gghpd3 <- ggs_caterpillar(gg21, horizontal = FALSE) + theme_bw() + ylab("Parametras")

#-----------------------------------Prognozë----------------------------------------------

progn3 <- outj32y$BUGSoutput$sims.list$y 
mr3 <- array(c(progn3), c(7500, 71, 10))
dimnames(outj32y$BUGSoutput$sims.matrix)[[2]][1:7]
proglong3 <- outj32y$BUGSoutput$sims.matrix[, -1:-7]
paror3 <- colnames(proglong3)
for(i in 1:(10 * 71)) {
  indF3 <- as.numeric(gsub(".*,(.*),.*", "\\1", paror3[i]))
  indI3 <- as.numeric(gsub(".*,(.*)]", "\\1", paror3[i]))
  mr3[, indI3, indF3] <- proglong3[, i]
}

progn3tr <- aperm(mr3, c(2:3, 1))
vidurpr3 <- apply(progn3tr, 1:2, mean) 
dfpr3 <- adply(vidurpr3, 2:1, identity) 
dfpr3$y13 <- y13_tikri$V1
dfpr3$w113 <- w113_tikri$V1

names(dfpr3) <- c("firmos", "veikla", "prognoze", "tikros", "pvm")

dataPro3 <- dfpr3[apply(dfpr3[,3:4], 1, function(row) all(row !=0 )),]

dataProgn3 <- data.frame(numeris = ordered$numeris,
  veikla = ordered$veikla, logy = ordered$logy, prognoze = dataPro3$prognoze,
  z = exp(ordered$logz), pvm = dataPro3$pvm)

dataProgn3$skirtumas <- dataProgn3$logy - dataProgn3$prognoze

#---------------------------Visø modeliø prognoziø palyginimas------------------------------------------

accuracy(dataProgn1$logy, dataProgn1$prognoze)
accuracy(dataProgn2$logy, dataProgn2$prognoze)
accuracy(dataProgn3$logy, dataProgn3$prognoze)

prodf <- data.frame(prognoze = c(dataProgn1$prognoze, dataProgn2$prognoze, dataProgn3$prognoze),
  logy = rep(dataProgn1$logy, 3), Modelis = factor(c(rep(1, 171), rep(2, 171), rep(3, 171))))
ggplot(data = prodf, aes(x = prognoze, y = logy, color = Modelis)) + geom_point(size = 2.3) + 
  xlab("Tikrosios reikðmë s") + ylab("Prognozuotos reikðmë s") + theme_bw() + 
  geom_abline(intercept = 0, slope = 1, colour = "#666666") 

densdf <- data.frame(dens = c(dataProgn1$prognoze, dataProgn2$prognoze, dataProgn3$prognoze, 
  dataProgn1$logy), kas = factor(c(rep("1 modelis", 171), rep("2 modelis", 171), 
    rep("3 modelis", 171), rep("Tikrosios", 171))))
ggplot(densdf, aes(x = dens, fill = kas)) + geom_density(alpha=.3) + theme_bw() +
  ylab("Tankis") + xlab("Reikðmë")  + 
  theme(legend.title = element_blank()) 

pakldf <- data.frame(dens = c(dataProgn1$skirtumas, dataProgn2$skirtumas, dataProgn3$skirtumas), 
  kas = factor(c(rep("1 modelis", 171), rep("2 modelis", 171), 
    rep("3 modelis", 171))))
ggplot(pakldf, aes(x = dens, fill = kas)) + geom_density(alpha=.3) + theme_bw() +
  ylab("Tankis") + xlab("Reikðmë")  + 
  theme(legend.title = element_blank())

paklvsprdf <- data.frame(skirt = pakldf$dens, 
  progn = c(dataProgn1$prognoze, dataProgn2$prognoze, dataProgn3$prognoze), 
  Modelis = factor(c(rep(1, 171), rep(2, 171), rep(3, 171))))
ggplot(data = paklvsprdf, aes(x = progn, y = skirt, color = Modelis)) + geom_point(size = 2.3) + 
  xlab("Prognozë") + ylab("Prognozë s paklaida") + theme_bw() + geom_hline(aes(intercept = 0), 
    color = "#666666", size = 1) 

ynames <- c("y[48,1,10]", "y[48,1,11]", "y[48,1,12]", "y[48,1,13]", "y[48,1,14]", "y[48,1,15]",
  "y[48,1,16]", "y[48,1,17]", "y[48,1,18]", "y[48,1,19]", "y[48,1,1]" ,"y[48,1,20]", "y[48,1,21]",
  "y[48,1,22]", "y[48,1,23]", "y[48,1,24]", "y[48,1,25]", "y[48,1,26]", "y[48,1,27]", "y[48,1,28]",
  "y[48,1,29]", "y[48,1,2]", "y[48,1,30]", "y[48,1,31]", "y[48,1,32]", "y[48,1,33]", "y[48,1,34]",
  "y[48,1,35]", "y[48,1,36]", "y[48,1,37]", "y[48,1,38]", "y[48,1,39]", "y[48,1,3]", "y[48,1,40]",
  "y[48,1,41]", "y[48,1,42]", "y[48,1,43]", "y[48,1,44]", "y[48,1,45]", "y[48,1,46]", "y[48,1,47]",
  "y[48,1,48]", "y[48,1,49]", "y[48,1,4]", "y[48,1,50]", "y[48,1,51]", "y[48,1,52]", "y[48,1,53]",
  "y[48,1,54]", "y[48,1,55]", "y[48,1,56]", "y[48,1,57]", "y[48,1,58]", "y[48,1,59]", "y[48,1,5]",
  "y[48,1,60]", "y[48,1,61]", "y[48,1,62]", "y[48,1,63]", "y[48,1,64]", "y[48,1,65]", "y[48,1,66]",
  "y[48,1,67]", "y[48,1,68]", "y[48,1,69]",  "y[48,1,6]", "y[48,1,70]", "y[48,1,71]", "y[48,1,7]",
  "y[48,1,8]", "y[48,1,9]", "y[48,10,43]", "y[48,2,11]", "y[48,2,12]", "y[48,2,14]", "y[48,2,15]",
  "y[48,2,16]", "y[48,2,17]", "y[48,2,18]", "y[48,2,19]", "y[48,2,20]", "y[48,2,21]", "y[48,2,22]",
  "y[48,2,24]", "y[48,2,25]" , "y[48,2,26]", "y[48,2,27]", "y[48,2,28]", "y[48,2,29]", "y[48,2,2]",
  "y[48,2,31]", "y[48,2,33]", "y[48,2,34]", "y[48,2,35]", "y[48,2,36]", "y[48,2,37]", "y[48,2,38]",
  "y[48,2,39]", "y[48,2,40]", "y[48,2,41]", "y[48,2,42]", "y[48,2,43]", "y[48,2,44]", "y[48,2,49]",
  "y[48,2,4]", "y[48,2,51]", "y[48,2,53]", "y[48,2,54]", "y[48,2,55]", "y[48,2,56]", "y[48,2,57]",
  "y[48,2,59]", "y[48,2,5]", "y[48,2,60]", "y[48,2,61]", "y[48,2,62]", "y[48,2,65]", "y[48,2,68]",
  "y[48,2,7]", "y[48,2,8]", "y[48,3,11]", "y[48,3,12]", "y[48,3,15]", "y[48,3,16]", "y[48,3,17]",
  "y[48,3,18]", "y[48,3,19]", "y[48,3,22]", "y[48,3,25]", "y[48,3,26]", "y[48,3,27]", "y[48,3,28]",
  "y[48,3,29]", "y[48,3,2]", "y[48,3,35]", "y[48,3,36]", "y[48,3,37]", "y[48,3,38]", "y[48,3,39]",
  "y[48,3,41]", "y[48,3,42]", "y[48,3,43]", "y[48,3,4]", "y[48,3,55]", "y[48,3,5]", "y[48,3,60]",
  "y[48,3,61]", "y[48,3,62]", "y[48,3,65]", "y[48,3,7]", "y[48,4,12]", "y[48,4,22]", "y[48,4,28]",
  "y[48,4,29]", "y[48,4,35]", "y[48,4,37]", "y[48,4,38]", "y[48,4,42]", "y[48,4,43]", "y[48,4,4]",
  "y[48,4,61]", "y[48,5,35]", "y[48,5,38]", "y[48,5,42]", "y[48,5,43]", "y[48,5,4]", "y[48,6,43]",
  "y[48,6,4]", "y[48,7,43]", "y[48,8,43]", "y[48,9,43]")     

sd1 <- outj13$BUGSoutput$sims.matrix[, dimnames(outj13$BUGSoutput$sims.matrix)[[2]]%in%ynames]
sd1 <- sd1[, match(ynames, colnames(sd1))]
sd2 <- outj20y$BUGSoutput$sims.matrix[, dimnames(outj20y$BUGSoutput$sims.matrix)[[2]]%in%ynames]
sd2 <- sd2[, match(ynames, colnames(sd2))]
sd3 <- outj30y$BUGSoutput$sims.matrix[, dimnames(outj30y$BUGSoutput$sims.matrix)[[2]]%in%ynames]
sd3 <- sd3[, match(ynames, colnames(sd3))]

sdpr1 <- rep(0, 171)
sdpr2 <- rep(0, 171)
sdpr3 <- rep(0, 171)
for(i in 1:171){
  sdpr1[i] <- sd(sd1[,i])  
  sdpr2[i] <- sd(sd2[,i])
  sdpr3[i] <- sd(sd3[,i])
}

sddf <- data.frame(sd = c(sdpr2, sdpr3), Modelis = c(rep("2 modelis", 171), rep("3 modelis", 171)),
  index = rep(1:171, 2))
ggplot(sddf, aes(x = index, y = sd, colour = Modelis)) + geom_line() + theme_bw() + 
  ylab("Prognozë s dispersija") + xlab("Indeksas")  

#---------------------------------HPD palyginimas-----------------------------------------

p1 <- c("beta1", "beta2", "tau", "tauu", "rho")
p2 <- c("lambda", "beta1")
p3 <- c("tauv", "beta1")
p32 <- c(paste("sigmav[", 1:71, "]", sep = ""))

hpd11 <-  as.mcmc(list(outj11.mcmc[[1]][, p1], outj11.mcmc[[2]][,p1], 
  outj11.mcmc[[3]][,p1]))
hpd12 <-  as.mcmc(list(outj20y.mcmc[[1]][, p1], outj20y.mcmc[[2]][,p1], 
  outj20y.mcmc[[3]][,p1]))
hpd13 <-  as.mcmc(list(outj31.mcmc[[1]][, p1], outj31.mcmc[[2]][,p1], 
  outj31.mcmc[[3]][,p1]))

P1 <- data.frame(Parameter = c("beta1", "beta2", "tau", "tauu", "rho"), 
  Label = c("Beta1", "Beta2", "Sigma", "Sigmau", "Rho"))

ggs_caterpillar(list("1 modelis" = ggs(hpd11, par_labels = P1), 
  "2 modelis" = ggs(hpd12, par_labels = P1), 
  "3 modelis" = ggs(hpd13, par_labels = P1))) + theme_bw() + ylab("") 

hpd22 <-  as.mcmc(list(outj20y.mcmc[[1]][, p2], outj20y.mcmc[[2]][,p2], 
  outj20y.mcmc[[3]][,p2]))
hpd23 <-  as.mcmc(list(outj31.mcmc[[1]][, p2], outj31.mcmc[[2]][,p2], 
  outj31.mcmc[[3]][,p2]))
P22 <- data.frame(Parameter = c("lambda", "beta1"), Label = c("Lambda", "Beta1"))
ggs_caterpillar(list("2 modelis" = ggs(hpd22, par_labels = P22, family = "lambda") , 
  "3 modelis" = ggs(hpd23, par_labels = P22, family = "lambda"))) + theme_bw() + ylab("") 
P3 <- data.frame(Parameter = c("beta1", "beta2", "tau", "tauu", "rho"), 
  Label = c("Beta1", "Beta2", "Sigma", "Sigmau", "Rho"))

ggs_caterpillar(list("1 modelis" = ggs(hpd11, par_labels = P1), 
  "2 modelis" = ggs(hpd12, par_labels = P1), 
  "3 modelis" = ggs(hpd13, par_labels = P1))) + theme_bw() + ylab("") 
hpd31 <-  as.mcmc(list(outj11.mcmc[[1]][, p3], outj11.mcmc[[2]][,p3], 
  outj11.mcmc[[3]][,p3]))
hpd32 <-  as.mcmc(list(outj20y.mcmc[[1]][, p3], outj20y.mcmc[[2]][,p3], 
  outj20y.mcmc[[3]][,p3]))
hpd33 <-  as.mcmc(list(outj30.mcmc[[1]][, p32], outj30.mcmc[[2]][,p32], 
  outj30.mcmc[[3]][,p32]))
P32 <- data.frame(Parameter = c("tauv", "beta1"), Label = c("Sigmav", "Beta1"))
ggs_caterpillar(list("1 modelis" = ggs(hpd31, par_labels = P32, family = "tauv") , 
  "2 modelis" = ggs(hpd32, par_labels = P32, family = "tauv"))) + theme_bw() + ylab("") 

mcmc2y <-  as.mcmc(list(outj20y.mcmc[[1]][, ynames], outj20y.mcmc[[2]][,ynames], 
  outj20y.mcmc[[3]][, ynames]))
mcmc3y <- as.mcmc(list(outj30.mcmc[[1]][, ynames], outj30.mcmc[[2]][,ynames], 
  outj30.mcmc[[3]][, ynames]))

ggs_caterpillar(list("2 modelis" = ggs(mcmc2y), "3 modelis" = ggs(mcmc3y)), horizontal = F) + 
  theme_bw() + ylab("Prognozë") +
  theme(axis.text.x = element_blank(), axis.ticks = element_blank())

#-------------------------3 modelio su 1 palyginimas, kur asimetrija didþiausia-------------------

lambdai <- rep(0, 71)
for( i in 1:71 ){
  lambdai[i] <- 14.6/sqrt(deita3$n[i])
}
range(lambdai)
# 4.616925 14.600000
which(lambdai == max(lambdai))
# 1  3  6  9 10 13 23 30 32 45 46 47 48 50 52 58 63 64 66 67 69 70 71

veiklos <- data.frame(senos = ordered$veikla, naujos = dataPro1$veikla)
tikros <- veiklos$senos[veiklos$naujos %in% which(lambdai == max(lambdai))]

trys1 <- progn3tr[1,,]
trys2 <- progn3tr[3,,]
trys3 <- progn3tr[6,,]
trys4 <- progn3tr[9,,]
trys5 <- progn3tr[10,,]
trys6 <- progn3tr[13,,]
trys7 <- progn3tr[23,,]
trys8 <- progn3tr[30,,]
trys9 <- progn3tr[32,,]
trys10 <- progn3tr[45,,]
trys11 <- progn3tr[46,,]
trys12 <- progn3tr[47,,]
trys13 <- progn3tr[48,,]
trys14 <- progn3tr[50,,]
trys15 <- progn3tr[52,,]
trys16 <- progn3tr[58,,]
trys17 <- progn3tr[63,,]
trys18 <- progn3tr[64,,]
trys19 <- progn3tr[66,,]
trys20 <- progn3tr[67,,]
trys21 <- progn3tr[69,,]
trys22 <- progn3tr[70,,]
trys23 <- progn3tr[71,,]

vienas1 <- progn1tr[1,,]
vienas2 <- progn1tr[3,,]
vienas3 <- progn1tr[6,,]
vienas4 <- progn1tr[9,,]
vienas5 <- progn1tr[10,,]
vienas6 <- progn1tr[13,,]
vienas7 <- progn1tr[23,,]
vienas8 <- progn1tr[30,,]
vienas9 <- progn1tr[32,,]
vienas10 <- progn1tr[45,,]
vienas11 <- progn1tr[46,,]
vienas12 <- progn1tr[47,,]
vienas13 <- progn1tr[48,,]
vienas14 <- progn1tr[50,,]
vienas15 <- progn1tr[52,,]
vienas16 <- progn1tr[58,,]
vienas17 <- progn1tr[63,,]
vienas18 <- progn1tr[64,,]
vienas19 <- progn1tr[66,,]
vienas20 <- progn1tr[67,,]
vienas21 <- progn1tr[69,,]
vienas22 <- progn1tr[70,,]
vienas23 <- progn1tr[71,,]


plot(trys1[trys1!=0], col = 2, type = "l")
lines(vienas1[vienas1!=0])
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[1]], col = 3)
plot(density(vienas1[vienas1!=0]))
lines(density(trys1[trys1!=0]), col = 2)

plot(trys2[trys2!=0], col = 2, type = "l")
lines(vienas2[vienas2!=0])
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[2]], col = 3)
lines(progn2tr[3,,][progn2tr[3,,]!=0], col = 4)
plot(density(vienas2[vienas2!=0]))
lines(density(trys2[trys2!=0]), col = 2)

summary(trys2[trys2!=0])
summary(vienas2[vienas2!=0])
dataProgn1$logy[dataProgn1$veikla == tikros[2]]

plot(vienas3[vienas3!=0], type = "l")
lines(trys3[trys3!=0], col = 2) 
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[3]], col = 3)
lines(progn2tr[6,,][progn2tr[6,,]!=0], col = 4)

plot(exp(trys3[trys3!=0]), col = 2, type = "l")
lines(exp(vienas3[vienas3!=0])) 
abline(h = exp(dataProgn1$logy[dataProgn1$veikla == tikros[3]]), col = 3)

plot(density(vienas3[vienas3!=0]))
lines(density(trys3[trys3!=0]), col = 2)
sd(vienas3[vienas3!=0])
sd(trys3[trys3!=0])

plot(vienas4[vienas4!=0], type = "l")
lines(trys4[trys4!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[4]], col = 3)
plot(density(vienas4[vienas4!=0]))
lines(density(trys4[trys4!=0]), col = 2)

plot(vienas5[vienas5!=0], type = "l")
lines(trys5[trys5!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[5]], col = 3)
plot(density(vienas5[vienas5!=0]))
lines(density(trys5[trys5!=0]), col = 2)

plot(vienas6[vienas6!=0], type = "l")
lines(trys6[trys6!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[6]], col = 3)
plot(density(vienas5[vienas5!=0]))
lines(density(trys5[trys5!=0]), col = 2)

plot(vienas7[vienas7!=0], type = "l")
lines(trys7[trys7!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[7]], col = 3)
plot(density(vienas7[vienas7!=0]))
lines(density(trys7[trys7!=0]), col = 2)

plot(vienas8[vienas8!=0], type = "l")
lines(trys8[trys8!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[8]], col = 3)
plot(density(vienas8[vienas8!=0]))
lines(density(trys8[trys8!=0]), col = 2)

plot(vienas9[vienas9!=0], type = "l")
lines(trys9[trys9!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[9]], col = 3)
plot(density(vienas9[vienas9!=0]))
lines(density(trys9[trys9!=0]), col = 2)

plot(vienas10[vienas10!=0], type = "l")
lines(trys10[trys10!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[10]], col = 3)
plot(density(vienas10[vienas10!=0]))
lines(density(trys10[trys10!=0]), col = 2)

plot(exp(vienas10[vienas10!=0]), type = "l")
lines(exp(trys10[trys10!= 0]), col = 2)
plot(density(exp(vienas10[vienas10!=0])))
lines(density(exp(trys10[trys10!=0]), col = 2))

plot(vienas11[vienas11!=0], type = "l")
lines(trys11[trys11!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[11]], col = 3)
plot(density(vienas11[vienas11!=0]))
lines(density(trys11[trys11!=0]), col = 2)

plot(vienas12[vienas12!=0], type = "l")
lines(trys12[trys12!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[12]], col = 3)
plot(density(vienas12[vienas12!=0]))
lines(density(trys12[trys12!=0]), col = 2)

plot(vienas13[vienas13!=0], type = "l")
lines(trys13[trys13!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[13]], col = 3)
plot(density(vienas13[vienas13!=0]))
lines(density(trys13[trys13!=0]), col = 2)

plot(vienas14[vienas14!=0], type = "l")
lines(trys14[trys14!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[14]], col = 3)
plot(density(vienas14[vienas14!=0]))
lines(density(trys14[trys14!=0]), col = 2)

plot(vienas15[vienas15!=0], type = "l")
lines(trys15[trys15!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[15]], col = 3)
plot(density(vienas15[vienas15!=0]))
lines(density(trys15[trys15!=0]), col = 2)

plot(vienas16[vienas16!=0], type = "l")
lines(trys16[trys16!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[16]], col = 3)
plot(density(vienas16[vienas16!=0]))
lines(density(trys16[trys16!=0]), col = 2)
plot(density(exp(vienas16[vienas16!=0])))
lines(density(exp(trys16[trys16!=0])), col = 2)

plot(vienas17[vienas17!=0], type = "l")
lines(trys17[trys17!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[17]], col = 3)
plot(density(vienas17[vienas17!=0]))
lines(density(trys17[trys17!=0]), col = 2)

plot(vienas18[vienas18!=0], type = "l")
lines(trys18[trys18!= 0], col = 2)
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[18]], col = 3)
plot(density(vienas16[vienas18!=0]))
lines(density(trys16[trys18!=0]), col = 2)

plot(trys19[trys19!=0], col = 2, type = "l")
lines(vienas19[vienas19!= 0])
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[19]], col = 3)
plot(density(vienas16[vienas19!=0]))
lines(density(trys16[trys19!=0]), col = 2)

plot(trys20[trys20!=0], col = 2, type = "l")
lines(vienas20[vienas20!= 0])
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[20]], col = 3)
plot(density(vienas20[vienas20!=0]))
lines(density(trys20[trys20!=0]), col = 2)

plot(trys21[trys21!=0], col = 2, type = "l")
lines(vienas21[vienas21!= 0])
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[21]], col = 3)
plot(density(vienas21[vienas21!=0]))
lines(density(trys21[trys21!=0]), col = 2)

plot(trys22[trys22!=0], col = 2, type = "l")
lines(vienas22[vienas22!= 0])
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[22]], col = 3)
plot(density(vienas22[vienas22!=0]))
lines(density(trys22[trys22!=0]), col = 2)

plot(trys23[trys23!=0], col = 2, type = "l")
lines(vienas23[vienas23!= 0])
abline(h = dataProgn1$logy[dataProgn1$veikla == tikros[23]], col = 3)
plot(density(vienas23[vienas23!=0]))
lines(density(trys23[trys23!=0]), col = 2)

sd(trys1[trys1!=0])
sd(vienas1[vienas1!=0])

sd(trys2[trys2!=0])
sd(vienas2[vienas2!=0])

sd(trys3[trys3!=0])
sd(vienas3[vienas3!=0])

sd(trys4[trys4!=0])
sd(vienas4[vienas4!=0])

sd(trys5[trys5!=0])
sd(vienas5[vienas5!=0])

sd(trys6[trys6!=0])
sd(vienas6[vienas6!=0])

sd(trys7[trys7!=0])
sd(vienas7[vienas7!=0])

sd(trys8[trys8!=0])
sd(vienas8[vienas8!=0]) 

sd(trys9[trys9!=0])
sd(vienas9[vienas9!=0]) 

sd(trys10[trys10!=0])
sd(vienas10[vienas10!=0]) 

sd(trys11[trys11!=0])
sd(vienas11[vienas11!=0]) 

sd(trys12[trys12!=0])
sd(vienas12[vienas12!=0])

sd(trys13[trys13!=0])
sd(vienas13[vienas13!=0])

sd(trys14[trys14!=0])
sd(vienas14[vienas14!=0]) 

sd(trys15[trys15!=0])
sd(vienas15[vienas15!=0]) 

sd(trys16[trys16!=0]) 
sd(vienas16[vienas16!=0])

sd(trys17[trys17!=0]) 
sd(vienas17[vienas17!=0])

sd(trys10[trys18!=0])
sd(vienas10[vienas18!=0]) 

sd(trys10[trys19!=0])
sd(vienas10[vienas19!=0])

sd(trys10[trys20!=0])
sd(vienas10[vienas20!=0])

sd(trys10[trys21!=0])
sd(vienas10[vienas21!=0])

sd(trys10[trys22!=0])
sd(vienas10[vienas22!=0])

sd(trys10[trys23!=0])
sd(vienas10[vienas23!=0])

#------------------Generuotos pajamos pagal 3 modelá--------------------------------------

beta1 <- 0.9686916
beta2 <- 0.01796294
rho <- 0.2227492
lambda <- 14.38185
sigma <- 1/0.03464309
sigmau <- 1/0.04385805
sigmav <- outj30y$BUGSoutput$mean$sigmav*2

u <- array(0, dim = c(48, 10, 71))
eps <- array(0, dim = c(48, 10, 71))
pajnew <- array(0, dim = c(48, 10, 71))

set.seed(321)
for ( i in 1:71 ) {
  W0 <- rnorm(1000, 0, 1)
  W <- W0[W0 > 0][1]
  lambdai <- lambda/sqrt(deita3$n[i])
  delta <- lambdai/sqrt( 1 + lambdai^2 )
  v <- rnorm(1, sigmav[i]*delta*W, sigmav[i]^2 * ( 1 - delta^2) )
	
	for ( f in 1:deita3$n[i] ) {
	u[1, f, i] <- rnorm(1, 0, 1/sigmau^2 ) 
	
	for( m in 2:48 ) { 
	eps[m, f, i] <- rnorm(1, 0, 1/sigma^2 )
	u[m, f, i] <- rho*u[m-1, f, i] + eps[m, f, i] 
	}
	
  for( t in 1:48){
  mu <- beta1*deita3$z1[t, f, i] + beta2*deita3$z2[t, f, i] + v + u[t, f, i]
  pajnew[t, f, i] <- rnorm(1, mu, deita3$kova[f, i])
  }
  }
  }	
range(pajnew[pajnew != 0])
range(y13[y13 != 0])
plot(pajnew[pajnew != 0], type = "l")

#-----------------------------Modeliai su sugeneruotais duomenimis---------------------------

pajnew1 <- pajnew[1:47, , ]  
pajnew13 <- pajnew[48, ,]
pajnewNA <- pajnew13
pajnewNA[which(pajnewNA!=0)] <- NA
pajbugs <- abind(pajnew1, pajnewNA, along = 1)

dei1 <- list(y = pajbugs, z1 = w1, z2 = w2, n = as.numeric(ns), ns2 = ns2, l = length(ns2),
  kova = ko)

system.time({new1 <- jags(data = dei1, inits = inits1, parameters.to.save = par11,
  model.file = "C:/Users/Artura/Desktop/FayHerriot.txt", 
  n.chains = 3, n.iter = 30000, n.thin = 10, DIC = TRUE) 
})/60 

geweke.diag(new1)
print(new1)
summary(new1$BUGSoutput$sims.list$beta1)
sd(new1$BUGSoutput$sims.list$beta1)
summary(new1$BUGSoutput$sims.list$beta2)
sd(new1$BUGSoutput$sims.list$beta2)
summary(new1$BUGSoutput$sims.list$rho)
sd(new1$BUGSoutput$sims.list$rho)
summary(new1$BUGSoutput$sims.list$tau)
sd(new1$BUGSoutput$sims.list$tau)
summary(new1$BUGSoutput$sims.list$tauu)
sd(new1$BUGSoutput$sims.list$tauu)

new1$BUGSoutput$last.value

dei2 <- list(y = pajbugs, z1 = w1, z2 = w2, n = as.numeric(ns), ns2 = ns2, l = length(ns2),
  kova = ko, kiek2 = kiek2)

system.time({new2 <- jags(data = dei2, inits = inits20, parameters.to.save = par21,
  model.file = "C:/Users/Artura/Desktop/SkewNormalY.txt", 
  n.chains = 3, n.iter = 30000, n.thin = 10, DIC = TRUE) 
})/60

new2$BUGSoutput$last.value

dei3 <- list(y = pajbugs, z1 = w1, z2 = w2, n = as.numeric(ns), ns2 = ns2, l = length(ns2),
  kova = ko)

system.time({new3 <- jags(data = dei3, inits = inits32, parameters.to.save = par30,
  model.file = "C:/Users/Artura/Desktop/SkewNormalV.txt", 
  n.chains = 3, n.iter = 30000, n.thin = 10, DIC = TRUE) 
})/60 

new3$BUGSoutput$last.value

summary(new3$BUGSoutput$sims.list$beta1)
sd(new3$BUGSoutput$sims.list$beta1)
summary(new3$BUGSoutput$sims.list$beta2)
sd(new3$BUGSoutput$sims.list$beta2)
summary(new3$BUGSoutput$sims.list$rho)
sd(new3$BUGSoutput$sims.list$rho)
summary(new3$BUGSoutput$sims.list$tau)
sd(new3$BUGSoutput$sims.list$tau)
summary(new3$BUGSoutput$sims.list$tauu)
sd(new3$BUGSoutput$sims.list$tauu)
summary(new3$BUGSoutput$sims.list$lambda)
sd(new3$BUGSoutput$sims.list$lambda)

progn12 <- new1$BUGSoutput$sims.list$y 
mr12 <- array(c(progn12), c(4500, 71, 10))
dimnames(new1$BUGSoutput$sims.matrix)[[2]][1:8]
proglong12 <- new1$BUGSoutput$sims.matrix[, -1:-7]
paror12 <- colnames(proglong12)
for(i in 1:(10 * 71)) {
  indF12 <- as.numeric(gsub(".*,(.*),.*", "\\1", paror12[i]))
  indI12 <- as.numeric(gsub(".*,(.*)]", "\\1", paror12[i]))
  mr12[, indI12, indF12] <- proglong12[, i]
}

progn12tr <- aperm(mr12, c(2:3, 1))
vidurpr12 <- apply(progn12tr, 1:2, mean)
dfpr12 <- adply(vidurpr12, 2:1, identity) 
paj13_tikri <- adply(pajnew13, 1:2, identity)
dfpr12$paj13 <- paj13_tikri$V1
names(dfpr12) <- c("firmos", "veikla", "prognoze", "tikros")
dataPro12 <- dfpr12[apply(dfpr12[,3:4], 1, function(row) all(row !=0 )),]

dataProgn12 <- data.frame(numeris = ordered$numeris,
  veikla = ordered$veikla, logy = dataPro12$tikros, prognoze = dataPro12$prognoze)

progn22 <- new2$BUGSoutput$sims.list$y 
mr22 <- array(c(progn22), c(4500, 71, 10))
dimnames(new2$BUGSoutput$sims.matrix)[[2]][1:9]
proglong22 <- new2$BUGSoutput$sims.matrix[, -1:-8]
paror22 <- colnames(proglong22)
for(i in 1:(10 * 71)) {
  indF22 <- as.numeric(gsub(".*,(.*),.*", "\\1", paror22[i]))
  indI22 <- as.numeric(gsub(".*,(.*)]", "\\1", paror22[i]))
  mr22[, indI22, indF22] <- proglong22[, i]
}

progn22tr <- aperm(mr22, c(2:3, 1))
vidurpr22 <- apply(progn22tr, 1:2, mean) 
dfpr22 <- adply(vidurpr22, 2:1, identity)
dfpr22$paj13 <- paj13_tikri$V1
names(dfpr22) <- c("firmos", "veikla", "prognoze", "tikros")
dataPro22 <- dfpr22[apply(dfpr22[,3:4], 1, function(row) all(row !=0 )),]

dataProgn22 <- data.frame(numeris = ordered$numeris,
  veikla = ordered$veikla, logy = dataPro22$tikros, prognoze = dataPro22$prognoze)

progn32 <- new3$BUGSoutput$sims.list$y 
mr32 <- array(c(progn32), c(4500, 71, 10))
dimnames(new3$BUGSoutput$sims.matrix)[[2]][1:79]
proglong32 <- new3$BUGSoutput$sims.matrix[, -1:-78]
paror32 <- colnames(proglong32)
for(i in 1:(10 * 71)) {
  indF32 <- as.numeric(gsub(".*,(.*),.*", "\\1", paror32[i]))
  indI32 <- as.numeric(gsub(".*,(.*)]", "\\1", paror32[i]))
  mr32[, indI32, indF32] <- proglong32[, i]
}

progn32tr <- aperm(mr32, c(2:3, 1))
vidurpr32 <- apply(progn32tr, 1:2, mean)
dfpr32 <- adply(vidurpr32, 2:1, identity)
dfpr32$paj13 <- paj13_tikri$V1
names(dfpr32) <- c("firmos", "veikla", "prognoze", "tikros")
dataPro32 <- dfpr32[apply(dfpr32[,3:4], 1, function(row) all(row !=0 )),]

dataProgn32 <- data.frame(numeris = ordered$numeris,
  veikla = ordered$veikla, logy = dataPro32$tikros, prognoze = dataPro32$prognoze)

xtable(accuracy(dataProgn12$logy, dataProgn12$prognoze), digits = 5) 
xtable(accuracy(dataProgn22$logy, dataProgn22$prognoze), digits = 5)
xtable(accuracy(dataProgn32$logy, dataProgn32$prognoze), digits = 5)

#------------------------------------Maþosios sritys--------------------------------------------

dataPan13 <- dataPan[grep("2013-12", dataPan$menesiai), ] 
data2013 <- dataPan13[dataPan13$numeris %in% fns,] 
ordered <- data2013[ with( data2013, order(veikla)) , ]
vepr <- unique(ordered$veikla)  
pop13 <- dataOut[dataOut$menesiai == "2013-12-01" & dataOut$veikla%in% vepr,] 

dens <- density(pop13$z)
df <- data.frame(x = dens$x, y = dens$y)
ggplot(df, aes(x, y)) + 
  geom_area(data = subset(df, x  >= 31), fill = "lightgreen") +
  geom_line() + theme_bw() + ylab("Tankis") + ggtitle("") + 
  xlab("Darbuotojø  skaiè ius")

veikl <- unique(pop13$veikla)
mz <- rep(0, 71)
for(i in 1:71){
  mz[i] <- mean(pop13$z[pop13$veikla == veikl[i]])
}
dfr <- data.frame(mz)
ggplot(dfr, aes(x = mz)) + geom_histogram(binwidth = 2.5) + theme_bw() + ylab("Daþnis") +
  xlab("Vidutinis darbuotojø  skaiè ius veikloje")

library(xtable)
xtable(table(pop13$apskritis))
 sort(c("alytus", "kaunas", "klaipeda", "marijampole", "panev", "siauliai", 
   "taurage", "telsiai", "utena", "vilnius")) 

#--------------------------------3 MODELIS----------------------------------------------
#--------------------------------Didelës firmos----------------------------------------
maziduomD1 <- pop13[pop13$z > 31,]
maziD1 <- maziduomD1[c("numeris", "veikla", "y")]
mazunrD1 <- unique(maziD1$numeris)
turimnrD1 <- dataProgn3$numeris[dataProgn3$numeris %in% mazunrD1] 
neturimnrD1 <- mazunrD1[!(mazunrD1 %in% turimnrD1)]
veiklD1 <- unique(maziD1$veikla)

thetaD1 <- rep(0, length(veiklD1))
for( i in 1:length(veiklD1)){
  thetaD1[i] <- sum(maziduomD1$y[maziduomD1$veikla == veiklD1[i]])
}

vevidD3 <- rep(0, length(veiklD1))
for(i in 1:length(veiklD1)){
  vevidD3[i] <- mean(exp(dataProgn3$prognoze[dataProgn3$veikla == veiklD1[i]]))
}

dfvidD1 <- data.frame(veiklD1, vevidD3)

dfproD1 <- data.frame(numeris = c(turimnrD1, neturimnrD1), turim = c(rep(TRUE, length(turimnrD1)), 
  c(rep(FALSE, length(neturimnrD1)))))

dfproD1$veikla <- rep(0, nrow(dfproD1))
for(i in 1:length(turimnrD1)){
  dfproD1$veikla[i] <- unique(maziD1$veikla[maziD1$numeris == turimnrD1[i]])
}
for(i in 1:length(neturimnrD1)){
  dfproD1$veikla[i + length(turimnrD1)] <- unique(maziD1$veikla[maziD1$numeris == neturimnrD1[i]])
}
dfproD1$haty3 <- rep(0, nrow(dfproD1))
for(i in 1:length(turimnrD1)){
  dfproD1$haty3[i] <- exp(dataProgn3$prognoze[dataProgn3$numeris == turimnrD1[i]])
}

dfvD1 <- dfproD1$veikla
for(i in 1:length(neturimnrD1)){
  dfproD1$haty3[i + length(turimnrD1)] <- dfvidD1$vevidD3[dfvidD1$veiklD1 == dfvD1[i + length(turimnrD1)]]
}

thetahatD3 <- rep(0, length(veiklD1))
for( i in 1:length(veiklD1)){
  thetahatD3[i] <- sum(dfproD1$haty3[dfproD1$veikla == veiklD1[i]])
}

sm1 <- data.frame(veikla = veiklD1, thetaD1, thetahatD3, Indeksas = 1:56, 
  santykis3 = thetahatD3/thetaD1)

tbl <- table(dfproD1$veikla[dfproD1$turim == "TRUE"])
tbl2 <- table(maziD1$veikla)*0
pop <-  table(maziD1$veikla)
imt <- c(tbl[intersect(names(tbl), names(tbl2))], tbl2[setdiff(names(tbl2), names(tbl))])
lentaD1 <- data.frame(pop)
names(lentaD1) <- c("veikla", "pop")
lentaD1$imtis <- imt[sort(names(pop))]

D1 <- data.frame(veikla = as.factor(c(lentaD1$veikla, lentaD1$veikla)),
  info = c(lentaD1$pop, lentaD1$imtis), kas = c(rep("populiacija", 56), rep("imtis", 56)))
ggplot(D1, aes(info, fill = kas)) + geom_bar() + theme_bw() + ylab("Daþnis") + 
  xlab("Firmø  skaiè ius") + scale_fill_grey(start = 0.3, end = .6) + 
  theme(legend.title = element_blank()) + ggtitle("Didelë s firmos")

#---------------------------Vilniaus didelës firmos-----------------------------------------

maziduomVD1 <- pop13[pop13$z > 31 & pop13$apskritis == 10,]
maziVD1 <- maziduomVD1[c("numeris", "veikla", "y")]
mazunrVD1 <- unique(maziVD1$numeris)
turimnrVD1 <- dataProgn3$numeris[dataProgn3$numeris %in% mazunrVD1] 
neturimnrVD1 <- mazunrVD1[!(mazunrVD1 %in% turimnrVD1)] 
veiklVD1 <- unique(maziVD1$veikla)

thetaVD1 <- rep(0, length(veiklVD1))
for( i in 1:length(veiklVD1)){
  thetaVD1[i] <- sum(maziVD1$y[maziVD1$veikla == veiklVD1[i]])
}

vevidVD3 <- rep(0, length(veiklVD1))
for(i in 1:length(veiklVD1)){
  vevidVD3[i] <- mean(exp(dataProgn3$prognoze[dataProgn3$veikla == veiklVD1[i]]))
}
dfvidVD1 <- data.frame(veiklVD1, vevidVD3)

dfproVD1 <- data.frame(numeris = c(turimnrVD1, neturimnrVD1), turim = c(rep(TRUE, length(turimnrVD1)), 
  c(rep(FALSE, length(neturimnrVD1)))))
dfproVD1$veikla <- rep(0, nrow(dfproVD1))
for(i in 1:length(turimnrVD1)){
  dfproVD1$veikla[i] <- maziVD1$veikla[maziVD1$numeris == turimnrVD1[i]]
}
for(i in 1:length(neturimnrVD1)){
  dfproVD1$veikla[i + length(turimnrVD1)] <- maziVD1$veikla[maziVD1$numeris == neturimnrVD1[i]]
}

dfproVD1$haty3 <- rep(0, nrow(dfproVD1))
for(i in 1:length(turimnrVD1)){
  dfproVD1$haty3[i] <- exp(dataProgn3$prognoze[dataProgn3$numeris == turimnrVD1[i]])
}

dfvVD1 <- dfproVD1$veikla

for(i in 1:length(neturimnrVD1)){
  dfproVD1$haty3[i + length(turimnrVD1)] <- dfvidVD1$vevidVD3[dfvidVD1$veiklVD1 == dfvVD1[i + length(turimnrVD1)]]
}

thetahatVD3 <- rep(0, length(veiklVD1))
for( i in 1:length(veiklVD1)){
  thetahatVD3[i] <- sum(dfproVD1$haty3[dfproVD1$veikla == veiklVD1[i]])
}

sm2 <- data.frame(veikla = veiklVD1, thetaVD1, thetahatVD3, Indeksas = 1:39, 
  santykis3 = thetahatVD3/thetaVD1)

tblV <- table(dfproVD1$veikla[dfproVD1$turim == "TRUE"])
tblV2 <- table(maziVD1$veikla)*0
popV <-  table(maziVD1$veikla)
imtV <- c(tbl[intersect(names(tblV), names(tblV2))], tblV2[setdiff(names(tblV2), names(tblV))])
lentaVD1 <- data.frame(popV)
names(lentaVD1) <- c("veikla", "pop")
lentaVD1$imtis <- imtV[sort(names(popV))]

D2 <- data.frame(veikla = as.factor(c(lentaVD1$veikla, lentaVD1$veikla)),
  info = c(lentaVD1$pop, lentaVD1$imtis), kas = c(rep("populiacija", 39), rep("imtis", 39)))
ggplot(D2, aes(info, fill = kas)) + geom_bar() + theme_bw() + ylab("Daþnis") + 
  xlab("Firmø skaièius") + scale_fill_grey(start = 0.3, end = .6) + 
  theme(legend.title = element_blank()) + ggtitle("Vilniaus didelë s firmos")

#------------------------Vilniaus ir Kauno didelës firmos------------------------------------

maziduomVKD1 <- pop13[pop13$z > 31 & pop13$apskritis%in%c(2, 10),]
maziVKD1 <- maziduomVKD1[c("numeris", "veikla", "y")]
mazunrVKD1 <- unique(maziVKD1$numeris)
turimnrVKD1 <- dataProgn3$numeris[dataProgn3$numeris %in% mazunrVKD1] 
neturimnrVKD1 <- mazunrVKD1[!(mazunrVKD1 %in% turimnrVKD1)]
veiklVKD1 <- unique(maziVKD1$veikla)

thetaVKD1 <- rep(0, length(veiklVKD1))
for( i in 1:length(veiklVKD1)){
  thetaVKD1[i] <- sum(maziVKD1$y[maziVKD1$veikla == veiklVKD1[i]])
}

vevidVKD3 <- rep(0, length(veiklVKD1))
for(i in 1:length(veiklVKD1)){
  vevidVKD3[i] <- mean(exp(dataProgn3$prognoze[dataProgn3$veikla == veiklVKD1[i]]))
}
dfvidVKD1 <- data.frame(veiklVKD1, vevidVKD3)

dfproVKD1 <- data.frame(numeris = c(turimnrVKD1, neturimnrVKD1), turim = c(rep(TRUE, length(turimnrVKD1)), 
  c(rep(FALSE, length(neturimnrVKD1)))))
dfproVKD1$veikla <- rep(0, nrow(dfproVKD1))
for(i in 1:length(turimnrVKD1)){
  dfproVKD1$veikla[i] <- maziVKD1$veikla[maziVKD1$numeris == turimnrVKD1[i]]
}
for(i in 1:length(neturimnrVKD1)){
  dfproVKD1$veikla[i + length(turimnrVKD1)] <- maziVKD1$veikla[maziVKD1$numeris == neturimnrVKD1[i]]
}

dfproVKD1$haty3 <- rep(0, nrow(dfproVKD1))
for(i in 1:length(turimnrVKD1)){
  dfproVKD1$haty3[i] <- exp(dataProgn3$prognoze[dataProgn3$numeris == turimnrVKD1[i]])
}
dfvVKD1 <- dfproVKD1$veikla
for(i in 1:length(neturimnrVKD1)){
  dfproVKD1$haty3[i + length(turimnrVKD1)] <- dfvidVKD1$vevidVKD3[dfvidVKD1$veiklVKD1 == dfvVKD1[i + length(turimnrVKD1)]]
}

thetahatVKD3 <- rep(0, length(veiklVKD1))
for( i in 1:length(veiklVKD1)){
  thetahatVKD3[i] <- sum(dfproVKD1$haty3[dfproVKD1$veikla == veiklVKD1[i]])
}

sm3 <- data.frame(veikla = veiklVKD1, thetaVKD1, thetahatVKD3, Indeksas = 1:47, 
  santykis3 = thetahatVKD3/thetaVKD1)

tblVK <- table(dfproVKD1$veikla[dfproVKD1$turim == "TRUE"])
tblVK2 <- table(maziVKD1$veikla)*0
popVK <-  table(maziVKD1$veikla)
imtVK <- c(tbl[intersect(names(tblVK), names(tblVK2))], 
  tblVK2[setdiff(names(tblVK2), names(tblVK))])
lentaVKD1 <- data.frame(popVK)
names(lentaVKD1) <- c("veikla", "pop")
lentaVKD1$imtis <- imtVK[sort(names(popVK))]
D3 <- data.frame(veikla = as.factor(c(lentaVKD1$veikla, lentaVKD1$veikla)),
  info = c(lentaVKD1$pop, lentaVKD1$imtis), kas = c(rep("populiacija", 47), rep("imtis", 47)))
ggplot(D3, aes(info, fill = kas)) + geom_bar() + theme_bw() + ylab("Daþnis") + 
  xlab("Firmø  skaiè ius") + scale_fill_grey(start = 0.3, end = .6) + 
  theme(legend.title = element_blank()) + ggtitle("Vilniaus ir Kauno didelë s firmos")

#-----------------------------Kauno didelës firmos--------------------------------------

maziduomKD1 <- pop13[pop13$z > 31 & pop13$apskritis == 2,]
maziKD1 <- maziduomKD1[c("numeris", "veikla", "y")]
mazunrKD1 <- unique(maziKD1$numeris)
turimnrKD1 <- dataProgn3$numeris[dataProgn3$numeris %in% mazunrKD1]
neturimnrKD1 <- mazunrKD1[!(mazunrKD1 %in% turimnrKD1)]
veiklKD1 <- unique(maziKD1$veikla)

thetaKD1 <- rep(0, length(veiklKD1))
for( i in 1:length(veiklKD1)){
  thetaKD1[i] <- sum(maziKD1$y[maziKD1$veikla == veiklKD1[i]])
}

vevidKD3 <- rep(0, length(veiklKD1))
for(i in 1:length(veiklKD1)){
  vevidKD3[i] <- mean(exp(dataProgn3$prognoze[dataProgn3$veikla == veiklKD1[i]]))
}
dfvidKD1 <- data.frame(veiklKD1, vevidKD3)

dfproKD1 <- data.frame(numeris = c(turimnrKD1, neturimnrKD1), 
  turim = c(rep(TRUE, length(turimnrKD1)), c(rep(FALSE, length(neturimnrKD1)))))
dfproKD1$veikla <- rep(0, nrow(dfproKD1))
for(i in 1:length(turimnrKD1)){
  dfproKD1$veikla[i] <- maziKD1$veikla[maziKD1$numeris == turimnrKD1[i]]
}
for(i in 1:length(neturimnrKD1)){
  dfproKD1$veikla[i + length(turimnrKD1)] <- maziKD1$veikla[maziKD1$numeris == neturimnrKD1[i]]
}

dfproKD1$haty3 <- rep(0, nrow(dfproKD1))
for(i in 1:length(turimnrKD1)){
  dfproKD1$haty3[i] <- exp(dataProgn3$prognoze[dataProgn3$numeris == turimnrKD1[i]])
}
dfvKD1 <- dfproKD1$veikla
for(i in 1:length(neturimnrKD1)){
  dfproKD1$haty3[i + length(turimnrKD1)] <- dfvidKD1$vevidKD3[dfvidKD1$veiklKD1 == dfvKD1[i + length(turimnrKD1)]]
}

thetahatKD3 <- rep(0, length(veiklKD1))
for( i in 1:length(veiklKD1)){
  thetahatKD3[i] <- sum(dfproKD1$haty3[dfproKD1$veikla == veiklKD1[i]])
}

sm4 <- data.frame(veikla = veiklKD1, thetaKD1, thetahatKD3, Indeksas = 1:33, 
  santykis3 = thetahatKD3/thetaKD1)

tblK <- table(dfproKD1$veikla[dfproKD1$turim == "TRUE"])
tblK2 <- table(maziKD1$veikla)*0
popK <-  table(maziKD1$veikla)
imtK <- c(tbl[intersect(names(tblK), 
  names(tblK2))], tblK2[setdiff(names(tblK2), names(tblK))])
lentaK1 <- data.frame(popK)
names(lentaK1) <- c("veikla", "pop")
lentaK1$imtis <- imtK[sort(names(popK))]
D4 <- data.frame(veikla = as.factor(c(lentaK1$veikla, lentaK1$veikla)),
  info = c(lentaK1$pop, lentaK1$imtis), kas = c(rep("populiacija", 33), rep("imtis", 33)))
ggplot(D4, aes(info, fill = kas)) + geom_bar() + theme_bw() + ylab("Daþnis") + 
  xlab("Firmø  skaiè ius") + scale_fill_grey(start = 0.3, end = .6) + 
  theme(legend.title = element_blank()) + ggtitle("Kauno didelë s firmos")

#----------------------- Klaipëdos, Panevëþio ir Ðiauliø didelës firmos--------------------------

maziduomKPS1 <- pop13[pop13$z > 31 & pop13$apskritis%in%c(3, 5, 6),]
maziKPS1 <- maziduomKPS1[c("numeris", "veikla", "y")]
mazunrKPS1 <- unique(maziKPS1$numeris)
turimnrKPS1 <- dataProgn3$numeris[dataProgn3$numeris %in% mazunrKPS1]
neturimnrKPS1 <- mazunrKPS1[!(mazunrKPS1 %in% turimnrKPS1)] 
veiklKPS1 <- unique(maziKPS1$veikla)

thetaKPS1 <- rep(0, length(veiklKPS1))
for( i in 1:length(veiklKPS1)){
  thetaKPS1[i] <- sum(maziKPS1$y[maziKPS1$veikla == veiklKPS1[i]])
}

vevidKPS3 <- rep(0, length(veiklKPS1))
for(i in 1:length(veiklKPS1)){
  vevidKPS3[i] <- mean(exp(dataProgn3$prognoze[dataProgn3$veikla == veiklKPS1[i]]))
}
dfvidKPS1 <- data.frame(veiklKPS1, vevidKPS3)

dfproKPS1 <- data.frame(numeris = c(turimnrKPS1, neturimnrKPS1), 
  turim = c(rep(TRUE, length(turimnrKPS1)), c(rep(FALSE, length(neturimnrKPS1)))))
dfproKPS1$veikla <- rep(0, nrow(dfproKPS1))
for(i in 1:length(turimnrKPS1)){
  dfproKPS1$veikla[i] <- maziKPS1$veikla[maziKPS1$numeris == turimnrKPS1[i]]
}
for(i in 1:length(neturimnrKPS1)){
  dfproKPS1$veikla[i + length(turimnrKPS1)] <- maziKPS1$veikla[maziKPS1$numeris == neturimnrKPS1[i]]
}

dfproKPS1$haty3 <- rep(0, nrow(dfproKPS1))
for(i in 1:length(turimnrKPS1)){
  dfproKPS1$haty3[i] <- exp(dataProgn3$prognoze[dataProgn3$numeris == turimnrKPS1[i]])
}
dfvKPS1 <- dfproKPS1$veikla
for(i in 1:length(neturimnrKPS1)){
  dfproKPS1$haty3[i + length(turimnrKPS1)] <- dfvidKPS1$vevidKPS3[dfvidKPS1$veiklKPS1 == dfvKPS1[i + length(turimnrKPS1)]]
}

thetahatKPS3 <- rep(0, length(veiklKPS1))
for( i in 1:length(veiklKPS1)){
  thetahatKPS3[i] <- sum(dfproKPS1$haty3[dfproKPS1$veikla == veiklKPS1[i]])
}

sm5 <- data.frame(veikla = veiklKPS1, thetaKPS1, thetahatKPS3, Indeksas = 1:35, 
  santykis3 = thetahatKPS3/thetaKPS1)

tblKPS <- table(dfproKPS1$veikla[dfproKPS1$turim == "TRUE"])
tblKPS2 <- table(maziKPS1$veikla)*0
popKPS <-  table(maziKPS1$veikla)
imtKPS <- c(tbl[intersect(names(tblKPS), 
  names(tblKPS2))], tblKPS2[setdiff(names(tblKPS2), names(tblKPS))])
lentaKPS1 <- data.frame(popKPS)
names(lentaKPS1) <- c("veikla", "pop")
lentaKPS1$imtis <- imtKPS[sort(names(popKPS))]
D5 <- data.frame(veikla = as.factor(c(lentaKPS1$veikla, lentaKPS1$veikla)),
  info = c(lentaKPS1$pop, lentaKPS1$imtis), kas = c(rep("populiacija", 35), rep("imtis", 35)))
ggplot(D5, aes(info, fill = kas)) + geom_bar() + theme_bw() + ylab("Daþnis") + 
  xlab("Firmø  skaiè ius") + scale_fill_grey(start = 0.3, end = .6) + 
  theme(legend.title = element_blank()) + ggtitle("Klaipë dos, Panevë þio, Ðiauliø  didelë s firmos")

#1 sritis
sd(thetahatD3)/mean(thetahatD3)
summary(thetaD1)/summary(thetahatD3)
mean(thetaD1/thetahatD3*100-100) #13.69422
mean(abs(thetaD1/thetahatD3*100-100)) #45.92745

#2 sritis
sd(thetahatVD3)/mean(thetahatVD3)
summary(thetaVD1)/summary(thetahatVD3)
mean(thetaVD1/thetahatVD3*100-100) #4.115098
mean(abs(thetaVD1/thetahatVD3*100-100)) #50.45139

#3 sritis
sd(thetahatKD3)/mean(thetahatKD3)
summary(thetaKD1)/summary(thetahatKD3)
mean(thetaKD1/thetahatKD3*100-100) #4.493751
mean(abs(thetaKD1/thetahatKD3*100-100)) #56.6392

#4 sritis
sd(thetahatVKD3)/mean(thetahatVKD3)
summary(thetaVKD1)/summary(thetahatVKD3)
mean(thetaVKD1/thetahatVKD3*100-100) #17.88202
mean(abs(thetaVKD1/thetahatVKD3*100-100)) #56.39683

# 5 sritis
sd(thetahatKPS3)/mean(thetahatKPS3)
summary(thetaKPS1)/summary(thetahatKPS3)
mean(thetaKPS1/thetahatKPS3*100-100) #5.230619
mean(abs(thetaKPS1/thetahatKPS3*100-100)) #42.97497

#---------------------------2 MODELIS--------------------------------------
#--------------------------Didelës firmos----------------------------------

vevidD2 <- rep(0, length(veiklD1))
for(i in 1:length(veiklD1)){
  vevidD2[i] <- mean(exp(dataProgn2$prognoze[dataProgn2$veikla == veiklD1[i]]))
}
dfvidD1$vevidD2 <- vevidD2
dfproD1$haty2 <- rep(0, nrow(dfproD1))
for(i in 1:length(turimnrD1)){
  dfproD1$haty2[i] <- exp(dataProgn2$prognoze[dataProgn2$numeris == turimnrD1[i]])
}
for(i in 1:length(neturimnrD1)){
  dfproD1$haty2[i + length(turimnrD1)] <- dfvidD1$vevidD2[dfvidD1$veiklD1 == dfvD1[i + length(turimnrD1)]]
}

thetahatD2 <- rep(0, length(veiklD1))
for( i in 1:length(veiklD1)){
  thetahatD2[i] <- sum(dfproD1$haty2[dfproD1$veikla == veiklD1[i]])
}

sm1$thetahatD2 <- thetahatD2

smdf1 <- data.frame(th = c(sm1$thetaD1, sm1$thetahatD2, sm1$thetahatD3), idx = rep(sm1$Indeksas, 3),
  mod = factor(c(rep("Tikrosios", 56), rep("2 modelis", 56), rep("3 modelis", 56))))
ggplot(data = smdf1, aes(x = idx, y = th, colour = mod)) + geom_line(size = 0.7) + 
  xlab("Veikla") + ylab("Pajamos") + theme_bw() +
  scale_x_continuous(breaks = c(1:56), labels = c(as.character(veiklD1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title = element_blank()) 

sm1$santykis2 <- thetahatD2/thetaD1
san1 <- data.frame(san = c(sm1$santykis2, sm1$santykis3), idx = rep(1:56, 2),
  mod = factor(c(rep("2 modelis", 56), rep("3 modelis", 56))))
ggplot(data = san1, aes(x = idx, y = san, colour = mod)) + geom_line(size = 0.7) + 
  xlab("Veikla") + ylab("Santykis") + theme_bw() +
  scale_x_continuous(breaks = c(1:56), labels = c(as.character(veiklD1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.title = element_blank()) 

#--------------------------------Vilniaus didelës firmos------------------------------------

vevidVD2 <- rep(0, length(veiklVD1))
for(i in 1:length(veiklVD1)){
  vevidVD2[i] <- mean(exp(dataProgn2$prognoze[dataProgn2$veikla == veiklVD1[i]]))
}
dfvidVD1$vevidVD2 <- vevidVD2
dfproVD1$haty2 <- rep(0, nrow(dfproVD1))
for(i in 1:length(turimnrVD1)){
  dfproVD1$haty2[i] <- exp(dataProgn2$prognoze[dataProgn2$numeris == turimnrVD1[i]])
}
for(i in 1:length(neturimnrVD1)){
  dfproVD1$haty2[i + length(turimnrVD1)] <- dfvidVD1$vevidVD2[dfvidVD1$veiklVD1 == dfvVD1[i + length(turimnrVD1)]]
}

thetahatVD2 <- rep(0, length(veiklVD1))
for( i in 1:length(veiklVD1)){
  thetahatVD2[i] <- sum(dfproVD1$haty2[dfproVD1$veikla == veiklVD1[i]])
}

sm2$thetahatVD2 <- thetahatVD2
sm2$santykis2 <- thetahatVD2/thetaVD1
smdf2 <- data.frame(th = c(sm2$thetaVD1, sm2$thetahatVD2, sm2$thetahatVD3), 
  idx = rep(sm2$Indeksas, 3),
  mod = factor(c(rep("Tikrosios", 39), rep("2 modelis", 39), rep("3 modelis", 39))))
ggplot(data = smdf2, aes(x = idx, y = th, colour = mod)) + geom_line(size = 0.7) + 
  xlab("Veikla") + ylab("Pajamos") + theme_bw() +
  scale_x_continuous(breaks = c(1:39), labels = c(as.character(veiklVD1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title = element_blank()) 

san2 <- data.frame(san = c(sm2$santykis2, sm2$santykis3), idx = rep(sm2$Indeksas, 2),
  mod = factor(c(rep("2 modelis", 39), rep("3 modelis", 39))))
ggplot(data = san2, aes(x = idx, y = san, colour = mod)) + geom_line(size = 0.7) + 
  xlab("Veikla") + ylab("Santykis") + theme_bw() +
  scale_x_continuous(breaks = c(1:39), labels = c(as.character(veiklVD1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.title = element_blank())  

#------------------------Vilniaus ir Kauno didelës firmos----------------------------

vevidVKD2 <- rep(0, length(veiklVKD1))
for(i in 1:length(veiklVKD1)){
  vevidVKD2[i] <- mean(exp(dataProgn2$prognoze[dataProgn2$veikla == veiklVKD1[i]]))
}
dfvidVKD1$vevidVKD2 <- vevidVKD2
dfproVKD1$haty2 <- rep(0, nrow(dfproVKD1))
for(i in 1:length(turimnrVKD1)){
  dfproVKD1$haty2[i] <- exp(dataProgn2$prognoze[dataProgn2$numeris == turimnrVKD1[i]])
}
for(i in 1:length(neturimnrVKD1)){
  dfproVKD1$haty2[i + length(turimnrVKD1)] <- dfvidVKD1$vevidVKD2[dfvidVKD1$veiklVKD1 == dfvVKD1[i + length(turimnrVKD1)]]
}

thetahatVKD2 <- rep(0, length(veiklVKD1))
for( i in 1:length(veiklVKD1)){
  thetahatVKD2[i] <- sum(dfproVKD1$haty2[dfproVKD1$veikla == veiklVKD1[i]])
}

sm3$thetahatKVD2 <- thetahatVKD2
sm3$santykis2 <- thetahatVKD2/thetaVKD1

smdf3 <- data.frame(th = c(sm3$thetaVKD1, sm3$thetahatKVD2, sm3$thetahatVKD3), 
  idx = rep(sm3$Indeksas, 3),
  mod = factor(c(rep("Tikrosios", 47), rep("2 modelis", 47), rep("3 modelis", 47))))
ggplot(data = smdf3, aes(x = idx, y = th, colour = mod)) + geom_line(size = 0.7) + 
  xlab("Veikla") + ylab("Pajamos") + theme_bw() +
  scale_x_continuous(breaks = c(1:47), labels = c(as.character(veiklVKD1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title = element_blank()) 

san3 <- data.frame(san = c(sm3$santykis2, sm3$santykis3), idx = rep(sm3$Indeksas, 2),
  mod = factor(c(rep("2 modelis", 47), rep("3 modelis", 47))))
ggplot(data = san3, aes(x = idx, y = san, colour = mod)) + geom_line(size = 0.7) + 
  xlab("Veikla") + ylab("Santykis") + theme_bw() +
  scale_x_continuous(breaks = c(1:47), labels = c(as.character(veiklVKD1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.title = element_blank()) 

#----------------------------Kauno didelës firmos----------------------------

vevidKD2 <- rep(0, length(veiklKD1))
for(i in 1:length(veiklKD1)){
  vevidKD2[i] <- mean(exp(dataProgn2$prognoze[dataProgn2$veikla == veiklKD1[i]]))
}
dfvidKD1$vevidKD2 <- vevidKD2
dfproKD1$haty2 <- rep(0, nrow(dfproKD1))
for(i in 1:length(turimnrKD1)){
  dfproKD1$haty2[i] <- exp(dataProgn2$prognoze[dataProgn2$numeris == turimnrKD1[i]])
}
for(i in 1:length(neturimnrKD1)){
  dfproKD1$haty2[i + length(turimnrKD1)] <- dfvidKD1$vevidKD2[dfvidKD1$veiklKD1 == dfvKD1[i + length(turimnrKD1)]]
}

thetahatKD2 <- rep(0, length(veiklKD1))
for( i in 1:length(veiklKD1)){
  thetahatKD2[i] <- sum(dfproKD1$haty2[dfproKD1$veikla == veiklKD1[i]])
}

sm4$thetahatKD2 <- thetahatKD2
sm4$santykis2 <- thetahatKD2/thetaKD1

smdf4 <- data.frame(th = c(sm4$thetaKD1, sm4$thetahatKD2, sm4$thetahatKD3), 
  idx = rep(sm4$Indeksas, 3),
  mod = factor(c(rep("Tikrosios", 33), rep("2 modelis", 33), rep("3 modelis", 33))))
ggplot(data = smdf4, aes(x = idx, y = th, colour = mod)) + geom_line(size = 0.7) + 
  xlab("Veikla") + ylab("Pajamos") + theme_bw() +
  scale_x_continuous(breaks = c(1:33), labels = c(as.character(veiklKD1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title = element_blank()) 

san4 <- data.frame(san = c(sm4$santykis2, sm4$santykis3), idx = rep(sm4$Indeksas, 2),
  mod = factor(c(rep("2 modelis", 33), rep("3 modelis", 33))))
ggplot(data = san4, aes(x = idx, y = san, colour = mod)) + geom_line(size = 0.7) + 
  xlab("Veikla") + ylab("Santykis") + theme_bw() +
  scale_x_continuous(breaks = c(1:33), labels = c(as.character(veiklKD1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.title = element_blank())  

#------------------------Klaipëdos, Panevëþio ir Ðiauliø didelës firmos----------------------------

vevidKPS2 <- rep(0, length(veiklKPS1))
for(i in 1:length(veiklKPS1)){
  vevidKPS2[i] <- mean(exp(dataProgn2$prognoze[dataProgn2$veikla == veiklKPS1[i]]))
}
dfvidKPS1$vevidKPS2 <- vevidKPS2
dfproKPS1$haty2 <- rep(0, nrow(dfproKPS1))
for(i in 1:length(turimnrKPS1)){
  dfproKPS1$haty2[i] <- exp(dataProgn2$prognoze[dataProgn2$numeris == turimnrKPS1[i]])
}
for(i in 1:length(neturimnrKPS1)){
  dfproKPS1$haty2[i + length(turimnrKPS1)] <- dfvidKPS1$vevidKPS2[dfvidKPS1$veiklKPS1 == dfvKPS1[i + length(turimnrKPS1)]]
}

thetahatKPS2 <- rep(0, length(veiklKPS1))
for( i in 1:length(veiklKPS1)){
  thetahatKPS2[i] <- sum(dfproKPS1$haty2[dfproKPS1$veikla == veiklKPS1[i]])
}

sm5$thetahatKPS2 <- thetahatKPS2
sm5$santykis2 <- thetahatKPS2/thetaKPS1

smdf5 <- data.frame(th = c(sm5$thetaKPS1, sm5$thetahatKPS2, sm5$thetahatKPS3), 
  idx = rep(sm5$Indeksas, 3),
  mod = factor(c(rep("Tikrosios", 35), rep("2 modelis", 35), rep("3 modelis", 35))))
ggplot(data = smdf5, aes(x = idx, y = th, colour = mod)) + geom_line(size = 0.7) + 
  xlab("Veikla") + ylab("Pajamos") + theme_bw() +
  scale_x_continuous(breaks = c(1:35), labels = c(as.character(veiklKPS1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(legend.title = element_blank()) 

san5 <- data.frame(san = c(sm5$santykis2, sm5$santykis3), idx = rep(sm5$Indeksas, 2),
  mod = factor(c(rep("2 modelis", 35), rep("3 modelis", 35))))
ggplot(data = san5, aes(x = idx, y = san, colour = mod)) + geom_line(size = 0.7) + 
  xlab("Veikla") + ylab("Santykis") + theme_bw() +
  scale_x_continuous(breaks = c(1:35), labels = c(as.character(veiklKPS1))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  theme(legend.title = element_blank())  

#1 sritis
sd(thetahatD2)/mean(thetahatD2)
summary(thetaD1)/summary(thetahatD2)
mean(thetaD1/thetahatD2*100-100) # 18.23581
mean(abs(thetaD1/thetahatD2*100-100)) #47.15487

#2 sritis
sd(thetahatVD2)/mean(thetahatVD2)
summary(thetaVD1)/summary(thetahatVD2)
mean(thetaVD1/thetahatVD2*100-100) #9.963168
mean(abs(thetaVD1/thetahatVD2*100-100)) #51.91808

#3 sritis 
sd(thetahatKD2)/mean(thetahatKD2)
summary(thetaKD1)/summary(thetahatKD2)
mean(thetaKD1/thetahatKD2*100-100)  #7.907968
mean(abs(thetaKD1/thetahatKD2*100-100)) #54.54248

#4 sritis
sd(thetahatVKD2)/mean(thetahatVKD2)
summary(thetaVKD1)/summary(thetahatVKD2)
mean(thetaVKD1/thetahatVKD2*100-100) #22.50138
mean(abs(thetaVKD1/thetahatVKD2*100-100)) # 57.3537

# 5 sritis
sd(thetahatKPS2)/mean(thetahatKPS2)
summary(thetaKPS1)/summary(thetahatKPS2)
mean(thetaKPS1/thetahatKPS2*100-100)#7.844988
mean(abs(thetaKPS1/thetahatKPS2*100-100)) #43.0625
