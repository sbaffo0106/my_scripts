#input libraries
library(dftools)
library(celestial)
library(data.table)
library(mvtnorm)
library(magicaxis)
library(plotrix)
library(xtable)
#
#define selection functions:20 BINS
analytic_codist0 <- function(x){ #all
  a=1.20483274
  b=2.41901608
  as.numeric(x>=8)*pmin(867.0,exp((x-b)/a))
}
analytic_codist1 <- function(x){ #sf
  a=1.16305515
  b=2.58691212
  as.numeric(x>=8)*pmin(867.0,exp((x-b)/a))
}
analytic_codist2 <- function(x){ #p
  a=1.0546817
  b=3.53108493
  as.numeric(x>=8)*pmin(867.0,exp((x-b)/a))
}
#
analytic_vol1 <- function(x){
  totarea1 = (54.93 + 57.44 + 56.93) 
  Omega1 = totarea1 / 3282.806
  Area1 = Omega1 * (analytic_codist0(x)^2)
  Volume1 = Area1 * (analytic_codist0(x)) / 3
  return(Volume1)
}
analytic_vol2 <- function(x){
  totarea2 = (54.93 + 57.44 + 56.93) 
  Omega2 = totarea2 / 3282.806
  Area2 = Omega2 * (analytic_codist1(x)^2)
  Volume2 = Area2 * (analytic_codist1(x)) / 3
  return(Volume2)
}
analytic_vol3 <- function(x){ 
  totarea3 = (54.93 + 57.44 + 56.93) 
  Omega3 = totarea3 / 3282.806
  Area3 = Omega3 * (analytic_codist2(x)^2)
  Volume3 = Area3 * (analytic_codist2(x)) / 3
  return(Volume3)
} 
#
actual_vol1 <- function(x){
  Area1 = Omega1 * (x^2)
  Volume1 = Area1 * (x)/3 
  return(Volume1)
}
actual_vol2 <- function(x){
  Area2 = Omega2 * (x^2)
  Volume2 = Area2 * (x)/3 
  return(Volume2)
}
actual_vol3 <- function(x){
  Area3 = Omega3 * (x^2)
  Volume3 = Area3 * (x)/3 
  return(Volume3)
}
# Start of main code
#
stub="/Users/antonio/Desktop/my_scripts/GSMF_R/data/"
#
# Read in catalogue
morphcat1=fread(paste0(stub,"allgalaxies.csv"))
morphcat2=fread(paste0(stub,"sfgalaxies.csv"))
morphcat3=fread(paste0(stub,"pgalaxies.csv"))
#
morphcat1$redshift=morphcat1$Zfinal
morphcat1$cd=morphcat1$comovingdist
morphcat1$mass=morphcat1$logmstarFINAL
morphcat1$masserr=morphcat1$dellogmstarFINAL
morphcat2$redshift=morphcat2$Zfinal
morphcat2$cd=morphcat2$comovingdist
morphcat2$mass=morphcat2$logmstarFINAL
morphcat2$masserr=morphcat2$dellogmstarFINAL
morphcat3$redshift=morphcat3$Zfinal
morphcat3$cd=morphcat3$comovingdist
morphcat3$mass=morphcat3$logmstarFINAL
morphcat3$masserr=morphcat3$dellogmstarFINAL
#
DoubleSchechter <- function(x, p) {
  mu = (10^(x-p[2]))
  log(10) * exp(-mu) * (10^p[1] * (mu)^(p[3]+1) + 10^p[4] * mu^(p[5]+1))
}
#
rho=1.0
f1 = function(x,r) as.numeric(r<=analytic_codist0(x))
f2 = function(x,r) as.numeric(r<=analytic_codist1(x))
f3 = function(x,r) as.numeric(r<=analytic_codist2(x))
totarea1 = (54.93 + 57.44 + 56.93) 
totarea2 = (54.93 + 57.44 + 56.93) 
totarea3 = (54.93 + 57.44 + 56.93) 
Omega1 = totarea1 / 3282.806  #[deg2 to sr] FOV
Omega2 = totarea2 / 3282.806  
Omega3 = totarea3 / 3282.806  
dVdr1 = function(r) r^2*Omega1*rho
dVdr2 = function(r) r^2*Omega2*rho
dVdr3 = function(r) r^2*Omega3*rho
#
morphx1=morphcat1[redshift <= 0.213]
morphx1$dist = cosdist(ref = 737, z = morphx1$redshift)$CoDist
morphx1$distlim=analytic_codist0(morphx1$mass)
morphx1$volume = actual_vol1(morphx1$dist)
morphx1$volumelim=analytic_vol1(morphx1$mass)
morph1=morphx1[morphx1$dist <= morphx1$distlim]
morphx2=morphcat2[redshift <= 0.213]
morphx2$dist = cosdist(ref = 737, z = morphx2$redshift)$CoDist
morphx2$distlim=analytic_codist1(morphx2$mass)
morphx2$volume = actual_vol2(morphx2$dist)
morphx2$volumelim=analytic_vol2(morphx2$mass)
morph2=morphx2[morphx2$dist <= morphx2$distlim]
morphx3=morphcat3[redshift <= 0.213]
morphx3$dist = cosdist(ref = 737, z = morphx3$redshift)$CoDist
morphx3$distlim=analytic_codist2(morphx3$mass)
morphx3$volume = actual_vol3(morphx3$dist)
morphx3$volumelim=analytic_vol3(morphx3$mass)
morph3=morphx3[morphx3$dist <= morphx3$distlim]
#
pinitt = c(log10(16.27e-4), 10.97, -0.53, log10(9.47e-4), -1.37)
pinitB = c(log10(9.47e-4), 10.97, -1.37)
#
survey1 = dffit(morph1$mass, x.err=morph1$masserr, r=morph1$dist, correct.lss.bias=T,
                keep.eddington.bias = FALSE, xmin = 7, xmax = 13, dx = 0.01,
                gdf = DoubleSchechter , p.initial=pinitt, #n.bootstrap = 1e2,
                selection=list(f1, dVdr1, min(morph1$dist),max(morph1$dist)))
survey2a = dffit(morph2$mass, x.err=morph2$masserr, r=morph2$dist, correct.lss.bias=T,
                 keep.eddington.bias = FALSE, xmin = 7, xmax = 13, dx = 0.01,
                 gdf = 'Schechter' , p.initial=pinitB, #n.bootstrap = 1e2,
                 selection=list(f2, dVdr2, min(morph2$dist),max(morph2$dist)))
survey2b = dffit(morph2$mass, x.err=morph2$masserr, r=morph2$dist, correct.lss.bias=T,
                 keep.eddington.bias = FALSE, xmin = 7, xmax = 13, dx = 0.01,
                 gdf = DoubleSchechter , p.initial=pinitt, #n.bootstrap = 1e2,
                 selection=list(f2, dVdr2, min(morph2$dist),max(morph2$dist)))
survey3 = dffit(morph3$mass, x.err=morph3$masserr, r=morph3$dist, correct.lss.bias=T,
                keep.eddington.bias = FALSE, xmin = 7, xmax = 13, dx = 0.01,
                gdf = DoubleSchechter , p.initial=pinitt, #n.bootstrap = 1e2,
                selection=list(f3, dVdr3, min(morph3$dist),max(morph3$dist)))
#
x_vals <- seq(8, 12, length.out = 100)
Schechter1 <- function(x, p) {
  mu = (10^(x - p[2]))
  log(10) * exp(-mu) * (10^p[1] * (mu)^(p[3] + 1))
}
Schechter2 <- function(x, p) {
  mu = (10^(x - p[2]))
  log(10) * exp(-mu) * (10^p[4] * mu^(p[5] + 1))
}
#
#DS components
pblack <- c(survey1[["fit"]][["p.best"]][1], survey1[["fit"]][["p.best"]][2],
            survey1[["fit"]][["p.best"]][3], survey1[["fit"]][["p.best"]][4],
            survey1[["fit"]][["p.best"]][5])
pblue <- c(survey2b[["fit"]][["p.best"]][1], survey2b[["fit"]][["p.best"]][2],
           survey2b[["fit"]][["p.best"]][3], survey2b[["fit"]][["p.best"]][4],
           survey2b[["fit"]][["p.best"]][5])
pred <- c(survey3[["fit"]][["p.best"]][1], survey3[["fit"]][["p.best"]][2],
          survey3[["fit"]][["p.best"]][3], survey3[["fit"]][["p.best"]][4],
          survey3[["fit"]][["p.best"]][5])
y1black <- sapply(x_vals, Schechter1, p = pblack)
y2black <- sapply(x_vals, Schechter2, p = pblack)
y1blue <- sapply(x_vals, Schechter1, p = pblue)
y2blue <- sapply(x_vals, Schechter2, p = pblue)
y1red <- sapply(x_vals, Schechter1, p = pred)
y2red <- sapply(x_vals, Schechter2, p = pred)
#
#Plot GSMF
#
png(filename=paste0(stub,"GSMFs.png"),width=22.0,height=15.0,units="cm",pointsize = 18,res=240)
#
mfplot(survey1, xlab=expression("Stellar mass [M"["\u0298"]~"h"[70]^-2~"]"),
       ylab=expression("Number density [Mpc"^-3*"dex"^-1~"h"[70]^3~"]"),
       xlim=c(1E8,1E12), ylim=c(1E-6,1E-1), p=NULL, veff=NULL,
       show.input.data = FALSE, show.data.histogram = FALSE,
       uncertainty.type=1, nbins=20, bin.xmin=8.3, bin.xmax=12,  col.fit='black',
       lwd.fit = 2, col.data.posterior = "black", margins = c(4.1, 4.1, 3.1, 2.1))
mfplot(survey2a, p=NULL, veff=NULL, show.input.data = TRUE, col.data.input='magenta',
       show.data.histogram = FALSE, uncertainty.type=1, add = TRUE, nbins=20,
       bin.xmin=8.3, bin.xmax=12,  col.fit='purple', lwd.fit = 2,
       col.data.posterior = "purple", axes=FALSE)
mfplot(survey2b, p=NULL, veff=NULL, show.input.data = TRUE, col.data.input='cyan',
       show.data.histogram = FALSE, uncertainty.type=1, add = TRUE, nbins=20,
       bin.xmin=8.3, bin.xmax=12,  col.fit='blue', lwd.fit = 2,
       col.data.posterior = "blue", axes=FALSE)
mfplot(survey3, p=NULL, veff=NULL, show.input.data = FALSE,
       show.data.histogram = FALSE, uncertainty.type=1, add = TRUE, nbins=20,
       bin.xmin=8.3, bin.xmax=12,  col.fit='red', lwd.fit = 2,
       col.data.posterior = "red", axes=FALSE)
#
lines(10^x_vals, y1black, col = "black", lwd = 2, lty = 2)
lines(10^x_vals, y2black, col = "black", lwd = 2, lty = 2)
lines(10^x_vals, y1blue, col = "blue", lwd = 2, lty = 2)
lines(10^x_vals, y2blue, col = "blue", lwd = 2, lty = 2)
lines(10^x_vals, y1red, col = "red", lwd = 2, lty = 2)
lines(10^x_vals, y2red, col = "red", lwd = 2, lty = 2)
dev.off()
#
#legend("bottomleft", legend=c("logM<=12.5","12.5<logM<=13","logM>13"), col=c("red","orange","pink"), lty=1, lwd = 2, cex = 0.75)
#legend("bottomleft", legend=c("logM<=12","12<logM<=12.5","12.5<logM<=13","13<logM<=13.5","logM>13.5"), col=c("red","orange","pink","magenta","purple"), lty=1, lwd = 2, cex = 0.75)
legend("bottomleft", 
       legend = c(expression(total),
                  expression(star-forming SS),
                  expression(star-forming DS),
                  expression(passive)),
       col = c("black", "purple", "blue", "red"), 
       lty = 1, 
       lwd = 2, 
       cex = 0.7)
#
dev.off()
###################################
#1. SS & DS
SingleSchechter <- function(x, p) {
  # p = (logPhi, logMstar, alpha)
  mu <- 10^(x - p[2])
  phi <- log(10) * 10^(p[1]) * mu^(p[3] + 1) * exp(-mu)
  return(phi)
}
DoubleSchechter <- function(x, p) {
  # p = (logPhi1, logMstar, alpha1, logPhi2, alpha2)
  mu <- 10^(x - p[2])
  phi1 <- 10^(p[1]) * mu^(p[3] + 1)
  phi2 <- 10^(p[4]) * mu^(p[5] + 1)
  phi <- log(10) * (phi1 + phi2) * exp(-mu)
  return(phi)
}

#2. Log-likelihood
logLikFunction <- function(p, x_data, y_data, model) {
  #p: SS/DS parameters; x_data&y_data: log(mass) and log(phi); model: gdf (SS/DS)
  y_pred <- sapply(x_data, model, p = p) #application to my data
  if(any(!is.finite(y_pred)) || any(y_pred <= 0)) return(1e10) #for NA, Inf, <=-0: it gives 1e10
  
  residuals <- y_data - y_pred
  sigma2 <- var(residuals) #variance
  if(sigma2 <= 0 || !is.finite(sigma2)) return(1e10)
  
  n <- length(y_data)
  logLik <- -0.5 * n * log(2 * pi * sigma2) - 0.5 * sum(residuals^2) / sigma2
  #compute the log-likelihood under the assumption of Gaussian (normal) errors 
  #with constant variance, i.e. standard log-likelihood formulafor a normal 
  #distribution, based on the residuals
  return(-logLik)  ##negative log-likelihood (can be minimized by optim())
}

#3. Data
set.seed(123) #ensures that random numbers (e.g. noise) are generated consistently 
#every time the script is run
x_data <- seq(8, 12, length.out = 50)
#true_params_double <- c(-2.76560169, 10.48320844, -0.07651236, -2.91736930, -1.44214279) #totalsf
#true_params_double <- c(-3.158797, 10.525404,  0.107700, -3.265090, -1.315379) #ggsf
#true_params_double <- c(-2.97373474, 10.37835055,  0.08056186, -3.07609870, -1.49222317) #ugsf
true_params_double <- c(-2.5909541, 10.8348893, -0.4760112, -6.0896810, -2.3031407) #totalpassive
y_true <- sapply(x_data, DoubleSchechter, p = true_params_double)

#A) additive Gaussian noise with a small sd
#noise <- rnorm(length(y_true), sd = 0.001 * max(y_true))
#y_data <- y_true + noise
#B) (multiplicative) log-normal noise model, with a small logaritmic sd (MORE REALISTIC)
sigma_noise <- 0.01  #1% logaritmic noise
z_true <- log(y_true)
z_noisy <- z_true + rnorm(length(z_true), mean = 0, sd = sigma_noise) #log space
y_data <- exp(z_noisy) #back to lineare space

#4. Fit with optim()
#initial parameters (giving p.best from dffit to help optimization!)
#init_single <- c(-2.821306, 10.744488, -1.291005) #totalsf
#init_single <- c(-3.161653, 10.789595, -1.173339) #ggsf
#init_single <- c(-3.048569, 10.683104, -1.352923) #ugsf
#init_double <- c(-2.76560169, 10.48320844, -0.07651236, -2.91736930, -1.44214279) #totalsf
#init_double <- c(-3.158797, 10.525404,  0.107700, -3.265090, -1.315379) #ggsf
#init_double <- c(-2.97373474, 10.37835055,  0.08056186, -3.07609870, -1.49222317) #ugsf
init_double <- c(-2.5909541, 10.8348893, -0.4760112, -6.0896810, -2.3031407) #totalpassive
#to minimize the negative log-likelihood for SS
fit_single <- optim(par = init_single, fn = logLikFunction,
                    x_data = x_data, y_data = y_data, model = SingleSchechter,
                    method = "L-BFGS-B",  #fast optimization algorithm that finds the 
                    #minimum of a function while specifying 
                    #lower and upper bounds on parameters
                    lower = c(-10, 8, -5), upper = c(1, 12, 1))
#to minimize the negative log-likelihood for DS
fit_double <- optim(par = init_double, fn = logLikFunction,
                    x_data = x_data, y_data = y_data, model = DoubleSchechter,
                    method = "L-BFGS-B",
                    lower = c(-10, 8, -5, -10, -5), upper = c(1, 12, 1, 1, 1))

#5. Log-likelihood calculus
logL_single <- -fit_single$value #in order to maximize Log-likelihood
logL_double <- -fit_double$value

#6. Log-likelihood ratio
LRT_stat <- -2 * (logL_single - logL_double)
p_value <- pchisq(LRT_stat, df = 2, lower.tail = FALSE)
#chi-square distribution with degrees of freedom
#lower.tail = FALSE is necessary to get the right-tail p-value, i.e the probability 
#of observing a test statistic as extreme or more extreme than the calculated LRT value

#7. Outputs
cat("Log-Likelihood (Single Schechter):", logL_single, "\n")
cat("Log-Likelihood (Double Schechter):", logL_double, "\n")
cat("Likelihood Ratio Statistic:", LRT_stat, "\n")
cat("p-value:", p_value, "\n")

#8. Plots
#A) log-log space
plot(10^x_data, y_data, pch=20,
     xlab="Mass", ylab=expression(log[10](Phi)), log="xy", col="black")
lines(10^x_data, sapply(x_data, SingleSchechter, p = fit_single$par), col="purple", lwd=2)
lines(10^x_data, sapply(x_data, DoubleSchechter, p = fit_double$par), col="blue", lwd=2)
legend("topright", legend=c("Data", "SS", "DS"),
       col=c("black", "purple", "blue"), pch=c(20, NA, NA), lty=c(NA, 1, 1))

#B) log-linear space
plot(10^x_data, y_data, pch=20,
     xlab="Mass", ylab=expression(Phi), log="x", col="black")
lines(10^x_data, sapply(x_data, SingleSchechter, p = fit_single$par), col="purple", lwd=2)
lines(10^x_data, sapply(x_data, DoubleSchechter, p = fit_double$par), col="blue", lwd=2)
legend("topright", legend=c("Data", "SS", "DS"),
       col=c("black", "purple", "blue"), pch=c(20, NA, NA), lty=c(NA, 1, 1))


#BIC/AIC (Weigel et al. 2016)
n <- length(y_data)
k_single <- length(fit_single$par)  # 3 per SS
k_double <- length(fit_double$par)  # 5 per DS

#AIC/BIC calculus & print results
AIC_single <- 2 * k_single - 2 * logL_single
AIC_double <- 2 * k_double - 2 * logL_double
BIC_single <- k_single * log(n) - 2 * logL_single
BIC_double <- k_double * log(n) - 2 * logL_double
cat("AIC Single Schechter:", AIC_single, "\n")
cat("AIC Double Schechter:", AIC_double, "\n")
cat("BIC Single Schechter:", BIC_single, "\n")
cat("BIC Double Schechter:", BIC_double, "\n")

#AIC/BIC ratios calculus & print results
RAIC <- AIC_single - AIC_double
RBIC <- BIC_single - BIC_double
cat(sprintf("RAIC = %.4f - (%.4f) = %.3f\n", AIC_single, AIC_double, RAIC))
cat(sprintf("RBIC = %.4f - (%.4f) = %.3f\n", BIC_single, BIC_double, RBIC))
#