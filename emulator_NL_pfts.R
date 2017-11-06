# -------------------------------------------------------------------
library(DiceKriging)
library(devtools)
library(fields)
library(sensitivity)

source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R")

# -------------------------------------------------------------------
# 1. Data loading
# -------------------------------------------------------------------

# 
X = read.table(file = 'paramsets8_nl.txt', header = TRUE)
# regional mean of sum of 2 year NPP in gC/m2
npp = read.table(file = 'npp.txt', header = TRUE)
#temp = read.table(file = 'temp.txt', header = TRUE)
#precip = read.table(file = 'precip.txt', header = TRUE)

#XX = cbind(X,npp,temp,precip)
XX = cbind(X,npp)
Xnl = XX[,-c(7,16)]
#Xnl = X
y = read.table(file = 'nl_pft.txt')
colnames(y) = 'nl_pft'

# Normalization
X.norm = normalize(Xnl)

ndims = ncol(X.norm)
nens = nrow(X.norm)

# -------------------------------------------------------------------
# The following section does the one-at-a-time sensitivity analysis, 
# you need to call this function.
#--------------------------------------------------------------------
oaat.design <- function(design, n, med = TRUE, hold = NULL){
  oamat <- NULL
  nd <- ncol(design)
  
  if(med){
    meandes <- apply(design,2,median)  
  }
  else{
    meandes <- apply(design,2,mean)
  }
  
  if(is.null(hold) == FALSE){
    meandes <- hold
  }
  
  mindes <- apply(design,2,min)
  maxdes <- apply(design,2,max)
  
  for (j in 1:nd){
    basemat <- matrix(meandes, nrow = n, ncol = nd , byrow = TRUE)
    vec <- seq(from = mindes[j], to = maxdes[j], length.out = n)
    basemat[ ,j] <- vec
    oamat <- rbind(oamat,basemat)
  
  }
  oamat
}

# -------------------------------------------------------------------
# 2. Fit emulator
# -------------------------------------------------------------------
# prediction points
n <- 21

X.oat <- oaat.design(X.norm, n = n)
colnames(X.oat) <- colnames(X.norm)

oat.mean.mat <- matrix(nrow = nrow(X.oat), 1)
oat.sd.mat <- matrix(nrow = nrow(X.oat), 1)

looSummary = function(loo, y){
  out =  sqrt(mean((loo$mean - y)^2))
  out
}


# linear 
fit1 = km(formula=~., design=X.norm, response=y,
         control = list(maxit=1e4))

pdf(file = 'NL_fit_linear.pdf', width = 3, height = 6)
plot(fit1)
dev.off()

loo1 = leaveOneOut.km(fit1, type = 'UK')
looSummary(loo1, y)
pred1 = predict(fit1, newdata = X.oat, type = 'UK')
hist(pred1$mean)


# Constant 
fit0 = km(formula=~1, design=X.norm, y,
         control = list(maxit=1e4))
pred0 = predict(fit0, newdata = X.oat, type = 'UK')
loo0 = leaveOneOut.km(fit0, type = 'UK')
looSummary(loo0, y)

pdf(file = 'NL_fit_const.pdf', width = 3, height = 6)
plot(fit0)
dev.off()


# constant with sqrt transform
fit2 = km(formula=~1, design=X.norm, response=sqrt(c(y, recursive = TRUE)),covtype = 'gauss',
                control = list(maxit=1e4))
pred2 = predict(fit2, newdata = X.oat, type = 'UK')
loo2 = leaveOneOut.km(fit2, type = 'UK')
sqrt(mean((loo2$mean^2 - y)^2))

pdf(file = 'NL_fit_const_sqrtTRANS.pdf', width = 3, height = 6)
plot(fit2)
dev.off()

# constant with log transform
fit3 = km(formula=~1, design=X.norm, response=log(c(y, recursive = TRUE)),
          control = list(maxit=1e4))
pred3 = predict(fit3, newdata = X.oat, type = 'UK')
loo3 = leaveOneOut.km(fit3, type = 'UK')

pdf(file = 'NL_fit_const_logTRANS.pdf', width = 3, height = 6)
plot(fit3)
dev.off()


# linear with sqrt transform
fit4 = km(formula=~., design=X.norm, response=sqrt(c(y, recursive = TRUE)),
          control = list(maxit=1e4))
pred4 = predict(fit4, newdata = X.oat, type = 'UK')
loo4 = leaveOneOut.km(fit4, type = 'UK')

pdf(file = 'NL_fit_linear_sqrtTRANS.pdf', width = 3, height = 6)
plot(fit4)
dev.off()

# linear with cube transform
fit5 = km(formula=~., design=X.norm, response=(c(y, recursive = TRUE))^3,
          control = list(maxit=1e4))
pred5 = predict(fit5, newdata = X.oat, type = 'UK')
loo5 = leaveOneOut.km(fit5, type = 'UK')

pdf(file = 'NL_fit_linear_cubeTRANS.pdf', width = 3, height = 6)
plot(fit5)
dev.off()


# RMSE
sqrt(mean((loo1$mean - y)^2))
sqrt(mean((loo0$mean - y)^2))
sqrt(mean((loo2$mean^2 - y)^2))
sqrt(mean((exp(loo3$mean) - y)^2))
sqrt(mean((loo4$mean^2 - y)^2))
sqrt(mean((loo5$mean^(1/3) - y)^2))


# constant prior
plot(c(y, recursive = TRUE), loo0$mean, col = 'black', ylim = c(0,0.5))
segments(c(y, recursive = TRUE), loo0$mean - (2*loo0$sd), c(y, recursive = TRUE), loo0$mean + (2*loo0$sd))

points(c(y, recursive = TRUE), loo2$mean^2, col = 'red')
segments(c(y, recursive = TRUE),
         (loo2$mean - 2*(loo2$sd))^2,
         c(y, recursive = TRUE),
         (loo2$mean + 2*(loo2$sd))^2, col = 'red'
)

points(c(y, recursive = TRUE), exp(loo3$mean), col = 'blue')
segments(c(y, recursive = TRUE),
         exp(loo3$mean - 2*(loo3$sd)),
         c(y, recursive = TRUE),
         exp(loo3$mean + 2*(loo3$sd)), col = 'blue'
)
abline(0,1)

# linear prior
plot(c(y, recursive = TRUE), loo1$mean, col = 'black', ylim = c(0,0.5))
segments(c(y, recursive = TRUE), loo1$mean - (2*loo1$sd), c(y, recursive = TRUE), loo1$mean + (2*loo1$sd))

points(c(y, recursive = TRUE), loo4$mean^2, col = 'red')
segments(c(y, recursive = TRUE),
         (loo4$mean - 2*(loo4$sd))^2,
         c(y, recursive = TRUE),
         (loo4$mean + 2*(loo4$sd))^2, col = 'red'
)

points(c(y, recursive = TRUE), exp(loo5$mean), col = 'blue')
segments(c(y, recursive = TRUE),
         exp(loo5$mean - 2*(loo5$sd)),
         c(y, recursive = TRUE),
         exp(loo5$mean + 2*(loo5$sd)), col = 'blue'
)
abline(0,1)

# -------------------------------------------------------------------
# choose the model
# -------------------------------------------------------------------
fit = fit4
pred = pred4
ex = 2

oat.mean.mat = pred$mean
oat.sd.mat = pred$sd

# -------------------------------------------------------------------
# plot
# -------------------------------------------------------------------
ccol = 1
ccol.fill = 1
pdf(file = 'NL_oaat_linear_sqrtTRANS.pdf',width = 12, height = 7)
par(mfrow = c(4,7), las = 1, mar = c(4,3,1.5,1), oma = c(0,5,0,0), fg = 'grey')
for(i in 1: ncol(X.norm)){
  
  ix <- seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oat[ix, i], oat.mean.mat[ix]^ex,
       pch = 1, col = (ccol)[1], type = 'l', cex = 1.2, lwd = 3,xlab = colnames(X.norm)[i],
       ylim = c(0.15,.35), 
       cex.lab = 1.5, cex.axis = 1.5)
  
  ccol.fill <- c(1)
  col.transp <- adjustcolor(ccol.fill, alpha = 0.5)
  
  polygon(x = c(X.oat[ix, i], rev(X.oat[ix, i])),
          y =c((oat.mean.mat[ix] - oat.sd.mat[ix])^ex, rev((oat.mean.mat[ix] + oat.sd.mat[ix])^ex)),
          col = col.transp, border = col.transp)
  
  if (i==1){axis (2, cex.axis = 1.5)
    mtext(side = 2, line = 3.5, text = 'needleleaf', las = 0, col = 'black')
  }
  
  if (i==8){axis (2, cex.axis = 1.5)
    mtext(side = 2, line = 3.5, text = 'needleleaf', las = 0, col = 'black')
  }
  
  if (i==15){axis (2, cex.axis = 1.5)
    mtext(side = 2, line = 3.5, text = 'needleleaf', las = 0, col = 'black')
  }
  
  if (i==22){axis (2, cex.axis = 1.5)
    mtext(side = 2, line = 3.5, text = 'needleleaf', las = 0, col = 'black')
  }
}
dev.off()

# -------------------------------------------------------------------
# 3. FAST sensitivity analysis
# -------------------------------------------------------------------
# generate the design to run the emulator at, using fast99
x = fast99(model = NULL, factors = colnames(X.norm), n = 1000,
           q = "qunif", q.arg = list(min = 0, max = 1))

fast.pred = predict(fit, newdata = x$X, type = 'UK')
fast = tell(x, fast.pred$mean^ex)
pdf(file = 'fast_nl_linear_sqrtTRANS.pdf', width = 12, height = 6)
par(las = 2, mar = c(10,4,2,1))
plot(fast)
dev.off()

