# fmecR.R
# R functions for FMEC
# Doug McNeall dougmcneall@gmail.com

library(devtools)
install_github(repo = "dougmcneall/hde")

library(hde)
library(RColorBrewer)
library(fields)
library(MASS)
library(DiceKriging)
library(stepwise)

source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/emtools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/imptools.R")
source("https://raw.githubusercontent.com/dougmcneall/packages-git/master/vistools.R")


# Example of adding  density plots to a pairs plot
dfunc.up <- function(x,y,...){
  require(MASS)
  require(RColorBrewer)
  
  rb = brewer.pal(9, "RdBu")
  br  = rev(rb)
  # function for plotting 2d kernel density estimates in pairs() plot.
  kde = kde2d(x,y)
  image(kde, col = rb, add = TRUE)
}

yg = brewer.pal(9, "YlGn")
ryb = brewer.pal(11, "RdYlBu")
byr = rev(ryb)
rb = brewer.pal(11, "RdBu")
br = rev(rb)

true.loo = function(X,y, lm = FALSE){
  # True Leave-one-out wrapper function
  out.mean = rep(NA, length(y))
  out.sd = rep(NA, length(y))
  
  for(i in 1:nrow(X)){
    X.trunc = X[-i, ]
    y.trunc = y[-i]
    
    X.target = matrix(X[i, ], nrow = 1)
    colnames(X.target) <- colnames(X)
    X.target.df = data.frame(X.target)
    
    
    if(lm){
      data.trunc = data.frame(y.trunc, X.trunc)
      colnames(data.trunc) = c('y', colnames(X.trunc))
      fit = lm(y~., data = data.trunc)
      pred = predict(fit, newdata = X.target.df)
      out.mean[i] = pred
    }
    else{
      
      fit = km(~., design = X.trunc, response = y.trunc)
      pred = predict(fit,newdata = X.target, type = 'UK')
      out.mean[i] = pred$mean
      out.sd[i = pred$sd]
    }
    
    
  }
  return(list(mean = out.mean, sd = out.sd))
}



loo.rmse = function(loo,y){
  out = sqrt(mean((loo$mean - y)^2))
  out
}


# David Sexton's two-step approach to building an emulator.
twoStep = function(X, y, nugget=NULL, nuggetEstim=FALSE, noiseVar=NULL, seed=NULL, trace=FALSE, maxit=100,
                   REPORT=10, factr=1e7, pgtol=0.0, parinit=NULL, popsize=100){
  
  control_list = list(trace=trace, maxit=maxit, REPORT=REPORT, factr=factr, pgtol=pgtol, pop.size=popsize)
  
  xvars = colnames(X)
  data = data.frame(response=y, x=X)
  colnames(data) <- c("response", xvars)
  nval = length(y)
  
  # Build the first emulator with a flat prior
  m0 = km(y ~ 1, design=X, response=y, nugget=nugget,
          nugget.estim=nuggetEstim, noise.var=noiseVar, control=control_list)
  
  coefs0 = m0@covariance@range.val 
  print('coefs0')
  print(coefs0)
  
  start.form = as.formula(paste("y ~ ", paste(xvars, collapse= "+")))
  
  # use BIC so that model is parsimonious and allow GP to pick up any other behaviour not
  # explained by the key linear terms      
  startlm = lm(start.form, data=data)
  #      print('before step')
  #      print(startlm)
  steplm = step(startlm, direction="both", k=log(nval), trace=TRUE)
  print('after step')
  print(steplm)
  form = as.formula(steplm)
  print('Formula')
  print(form)
  data$response = NULL
  labels = labels(terms(steplm))
  labels = labels[!(labels %in% c('response'))]
  if (length(labels) > 0) {
    start.form = as.formula(paste("~ ", paste(labels, collapse= "+")))
  } else {
    start.form = as.formula("~ 1")
  }    
  print("Step has found formula:")
  print(start.form)
  if (!is.null(seed)) {set.seed(seed)}
  m = km(start.form, design=X, response=y, nugget=nugget, parinit=coefs0,
         nugget.estim=nuggetEstim, noise.var=noiseVar, control=control_list)
  
  return(list(x=X, y=y, nugget=nugget, nugget.estim=nuggetEstim,
              noise.var=noiseVar, emulator=m, seed=seed, coefs=m@covariance@range.val,
              trends=m@trend.coef, meanTerms=all.vars(start.form), steplm = steplm))
  
}





