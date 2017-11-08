# exeter_meeting.R
# Code written during the Exeter FMEC meeting,
# 6th-8th November 2017


library(DiceKriging)
library(devtools)
library(fields)
library(sensitivity)
library(RColorBrewer)
library(MASS)

source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R")

source("fmecR.R")
# -------------------------------------------------------------------
# 1. Data loading
# -------------------------------------------------------------------

# 
X = read.table(file = 'paramsets_nl.txt', header = TRUE)
# regional mean of sum of 2 year NPP in gC/m2
npp = read.table(file = 'npp.txt', header = TRUE)
#temp = read.table(file = 'temp.txt', header = TRUE)
#precip = read.table(file = 'precip.txt', header = TRUE)

XX = cbind(X,npp)
Xnl = XX[,-c(7,16)]
#Xnl = X
y = read.table(file = 'nl_pft.txt')
colnames(y) = 'nl_pft'

#obs = readnl_pft_regionalmean.txt

# Normalization
X.norm = normalize(Xnl)

ndims = ncol(X.norm)
nens = nrow(X.norm)

# -------------------------------------------------------------------
# 2. Fit Emulator
# -------------------------------------------------------------------

# linear with Cubic transform
fit5 = km(formula=~., design=X.norm, response=(c(y, recursive = TRUE))^3,
          control = list(maxit=1e4))
loo5 = leaveOneOut.km(fit5, trend.reestim = TRUE, type = 'UK')

# plot in cube space
pdf(file = 'NL_fit_linear_cubeTRANS.pdf', width = 3, height = 6)
plot(fit5)
dev.off()
# plot in original space
pdf(file = 'NL_fit_linear_cubeTRANS_original.pdf', width = 6, height = 6)
plot(c(y, recursive = TRUE), loo5$mean^(1/3), main = 'original scale', xlim = c(0,0.3))
abline(0,1)
dev.off()


# linear with sqrt transform
fit4 = km(formula=~., design=X.norm, response=sqrt(c(y, recursive = TRUE)),
          control = list(maxit=1e4))
loo4 = leaveOneOut.km(fit4, trend.reestim = TRUE, type = 'UK')

#sqrt space
pdf(file = 'NL_fit_linear_sqrtTRANS.pdf', width = 3, height = 6)
plot(fit4)
dev.off()
#original space
pdf(file = 'NL_fit_linear_sqrtTRANS_original.pdf', width = 6, height = 6)
plot(c(y, recursive = TRUE), loo4$mean^2, main = 'original scale', xlim = c(0,0.3))
abline(0,1)
dev.off()


# -------------------------------------------------------------------
# 3. Choose fit
# -------------------------------------------------------------------
fit = fit5

# Plot a pairs plot with density
X.unif = samp.unif(10000, mins = rep(0, ncol(X.norm)), maxes = rep(1, ncol(X.norm)))
colnames(X.unif) <- colnames(X.norm)

pred.unif = predict(fit, newdata = X.unif, type = 'UK')

X.kept = X.unif[pred.unif$mean < 0.5, ]


rb = brewer.pal(9, "RdBu")
br = rev(rb)

pdf(file = 'pairs_dens.pdf', width = 8, height = 8)
par(oma = c(0,0,0,3))

# Emulate all input space and keep only those inputs which match a 
# criteria (such as having absolute error below a threshold)

test = pairs(X.kept,
             gap = 0, lower.panel = NULL, xlim = c(0,1), ylim = c(0,1),
             panel = dfunc.up)

image.plot(legend.only = TRUE,
           zlim = c(0,1),
           col = rb,
           legend.args = list(text = 'Density of model runs matching the criteria', side = 3, line = 1),
           horizontal = TRUE
)
dev.off()


# -------------------------------------------------------------------
# One-at-a-time sensitivity analysis, with inputs held at
# non-central values
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


n <- 21
hold = rep(0.5, ncol(X.norm))
hold[which(colnames(X.norm) =='a_wl')] = 0.1 
hold[which(colnames(X.norm) =='b_wl')] = 0.1
hold[which(colnames(X.norm) =='g_area')] = 0.2
hold[which(colnames(X.norm) =='lai_min')] = 0.2

X.oat = oaat.design(X.norm, n = n, med = FALSE, hold = hold)
colnames(X.oat) <- colnames(X.norm)
oat.pred = predict(fit4, newdata = X.oat, type = 'UK') 
  
#oat.mean.mat <- matrix(nrow = nrow(X.oat), 1)
#oat.sd.mat <- matrix(nrow = nrow(X.oat), 1)

pdf(file = 'oat.pdf',width = 12, height = 7)
par(mfrow = c(4,7), las = 1, mar = c(4,3,1.5,1), oma = c(0,5,0,0), fg = 'grey')
for(i in 1: ncol(X.norm)){
  
  ix <- seq(from = ((i*n) - (n-1)), to =  (i*n), by = 1)
  plot(X.oat[ix, i], oat.pred$mean[ix]^2,
       pch = 1, col = 'black', type = 'l', cex = 1.2, lwd = 3,
       ylim = c(0.05,0.25),
       xlab = colnames(X.norm)[i],
       cex.lab = 1.5, cex.axis = 1.5)
  
  ccol.fill <- 'grey'
  col.transp <- adjustcolor(ccol.fill, alpha = 0.5)
  
  polygon(x = c(X.oat[ix, i], rev(X.oat[ix, i])),
          y =c((oat.pred$mean[ix] - oat.pred$sd[ix])^2, rev((oat.pred$mean[ix] + oat.pred$sd[ix])^2)),
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



# -----------------------------------------------------------------------------
# Build an emulator in a two-step process - using both a stepwise linear
# model, and an initial "flat prior" Gaussian process, to initialise the
# hyperparameters and linear coefficients of a GP with a linear prior.
# -----------------------------------------------------------------------------

# transform the output for accuracy and convenience 
y_sqrt = c(sqrt(y), recursive = TRUE)

tsfit = twoStep(X.norm, y_sqrt, nugget=NULL, nuggetEstim=FALSE, noiseVar=NULL, seed=NULL, trace=FALSE, maxit=100,
                   REPORT=10, factr=1e7, pgtol=0.0, parinit=NULL, popsize=100)

# Initial leave-one-out prediction test of the two emulators
tsloo = leaveOneOut.km(tsfit$emulator, trend.reestim = TRUE, type = 'UK')
fit4loo = leaveOneOut.km(fit4, trend.reestim = TRUE, type = 'UK')

loo.rmse(fit4loo, y_sqrt)
loo.rmse(tsloo, y_sqrt)

gp.col = 'skyblue2'
ts.col = 'tomato2'

pdf(file = 'fits.pdf', width = 8, height = 5)
par(mfrow = c(1,2))
plot(y_sqrt, tsloo$mean, main = 'sqrt transform', col = ts.col)
points(y_sqrt, fit4loo$mean, col = gp.col)
abline(0,1)
legend('topleft', legend = c('GP', 'twoStep'), col = c(gp.col, ts.col), pch = 21)


plot(c(y, recursive = TRUE), tsloo$mean^2, col = ts.col, main = 'original scale')
points(c(y, recursive = TRUE), fit4loo$mean^2, col = gp.col)
abline(0,1)
legend('topleft', legend = c('GP', 'twoStep'), col = c(gp.col, ts.col), pch = 21)
dev.off()


# Where do we go in the parameter space (including to the edges)
# in order to make needleleaf fraction as large as possible.
# [the 'Optimisation code']

# attach output to inputs and make a data.frame
X.y_sqrt = data.frame(cbind(X.norm,y_sqrt))
initfit = lm(y_sqrt ~ ., data = X.y_sqrt)
stepfit = step(initfit, direction="both", k=log(length(y_sqrt)), trace=TRUE)

# This is the cost function, gets fed to optim
fn.step = function(newdata, cn){
  newdata.df  = data.frame(matrix(newdata, nrow = 1))
  colnames(newdata.df) = cn
  out = predict(stepfit, newdata = newdata.df)
  out
}

# initial values for optim in the middle of the design
startin.mat <- matrix(rep(0.5, ncol(X.norm)), nrow = 1)
startin <- data.frame(startin.mat)
colnames(startin) <- colnames(X.norm)

# Find the values which minimise needleleaf absolute error
best.X <- optim(par = startin,
              fn = fn.step,
              method = "L-BFGS-B",
              lower = rep(0,ncol(X.norm)),
              upper = rep(1,ncol(X.norm)),
              control = list(maxit = 2000),
              cn = colnames(X.norm)
)

best.X$par[best.X$par!=0.5]

pdf(file = 'bestX.pdf', width = 7, height = 6)
par(mfrow = c(2,1))
plot(best.X$par, axes = FALSE, xlab = '', ylab = 'normalised best value')
axis(1, at = 1:length(best.X$par), labels = names(best.X$par), las = 2)
axis(2)

nd = data.frame(matrix(best.X$par, nrow = 1))
colnames(nd) = colnames(X.norm)
best_y_sqrt = predict(stepfit, newdata = nd)
best_y = best_y_sqrt^2

hist(c(y, recursive = TRUE), xlim = c(0,0.4), xlab = 'Needleleaf fraction abs. error',
     main = '')
rug(best_y, col = 'red', lwd = 3)
legend('topleft', legend = 'emulated best value', pch = '|', col = 'red', bty = 'n')

dev.off()


# Quilt plots that visualise how emulated absolute model error changes
# across parameter space.

X.unif = samp.unif(100000, mins = rep(0, ncol(X.norm)), maxes = rep(1, ncol(X.norm)) )
pred.unif = predict(fit4, newdata = X.unif, type = 'UK')
pred.unif.y = pred.unif$mean^2

pdf(file = 'quilts.pdf', width = 8, height = 8)
par(mfrow = c(2,2), oma = c(1,1,1,1))

quilt.plot(x = X.unif[ ,1],
           y = X.unif[, 2],
           z = pred.unif.y,
           xlab = colnames(X.norm)[1],
           ylab = colnames(X.norm)[2],
           col = byr
)

quilt.plot(x = X.unif[ ,6],
           y = X.unif[, 7],
           z = pred.unif.y,
           xlab = colnames(X.norm)[6],
           ylab = colnames(X.norm)[7],
           col = byr
)

quilt.plot(x = X.unif[ ,9],
           y = X.unif[, 10],
           z = pred.unif.y,
           xlab = colnames(X.norm)[9],
           ylab = colnames(X.norm)[10],
           col = byr
)

quilt.plot(x = X.unif[ ,14],
           y = X.unif[, 17],
           z = pred.unif.y,
           xlab = colnames(X.norm)[14],
           ylab = colnames(X.norm)[17],
           col = byr
)

dev.off()




