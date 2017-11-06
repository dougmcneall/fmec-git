# exeter_meeting.R
# Code written during the Exeter FMEC meeting,
# 6th-8th November 2017


Doug testing

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

rb = brewer.pal(9, "RdBu")
br = rev(rb)

pdf(file = 'pairs_dens.pdf', width = 8, height = 8)
par(oma = c(0,0,0,3))

# Emulate all input space and keep only those inputs which match a 
# criteria (such as having absolute error below a threshold)

test = pairs(X.unif[ix.kept, ],
             gap = 0, lower.panel = NULL, xlim = c(0,1), ylim = c(0,1),
             panel = dfunc.up)

image.plot(legend.only = TRUE,
           zlim = c(0,1),
           col = rb,
           legend.args = list(text = 'Density of model runs matching the criteria', side = 3, line = 1),
           horizontal = TRUE
)
dev.off()
