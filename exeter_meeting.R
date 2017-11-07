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

# Fit emulator
# linear with sqrt transform
fit4 = km(formula=~., design=X.norm, response=sqrt(c(y, recursive = TRUE)),
          control = list(maxit=1e4))
pred4 = predict(fit4, newdata = X.oat, type = 'UK')
loo4 = leaveOneOut.km(fit4, type = 'UK')


# Plot a pairs plot with density
X.unif = samp.unif(10000, mins = rep(0, ncol(X.norm)), maxes = rep(1, ncol(X.norm)))
colnames(X.unif) <- colnames(X.norm)

pred.unif = predict(fit4, newdata = X.unif, type = 'UK')

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







