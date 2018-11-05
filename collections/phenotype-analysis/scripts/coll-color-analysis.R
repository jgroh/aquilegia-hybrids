# ===== COMMENTS =====
# Data are values of red and green reflectance in RGB color space from selections
# of photgraphs of Aquilegia sepals in natural populations. 
#This script calculates a metric (log R/G) to reflect color variation
# between red-flowered A. formosa, yellow-flowered A. flavescens, and their hybrids. 
# This metric is shown to correlate with crude visually-assigned integer values representing
# relative extent of red over yellow on 1-5 scale. Summary statistics are then calculated, 
# and data are plotted with means and 95% CIs by population. 


# ===== READ DATA =====
cl <- read.csv("phenotype-analysis/data/collections-phenotypes.csv", stringsAsFactors = FALSE)

# create log(r/g) variable
log.rg <- log( cl$red.mean / cl$green.mean)

# filter out NAs
ids <- which(!is.na(log.rg))
log.rg <- log.rg[ids]

# new data frame
cldat <- cbind.data.frame(population = cl$site[ids], log.rg)

# quick look
stripchart(log.rg ~ population, data = cldat, vertical = TRUE, method = "jitter", col = 1:4)

# ===== COMPARE COLOR METRICS =====

plot(log.rg ~ cl$integer.sepal.color[ids], lwd = 2, 
     xlab = "Visually-assigned score", ylab = "Log R/G metric")

abline(lm(log.rg ~ cl$integer.sepal.color[ids]))

cor(log.rg, cl$integer.sepal.color[ids], use = "complete.obs")
# 0.853178

# ===== CALCULATE SUMMARIES BY POPULATION =====

# fit model
z <- lm(log.rg ~ 1 + population, data = cldat)

#  obtain means
est <- z$coefficients
means <- est + c(0, rep(est[1], 3))

# inspect
means

# put in desired order for plot
means <- means[c(2, 3, 1, 4)]

# get confidence intervals
cnf <- confint(z)
cnf.mat <- matrix(nrow = 4, ncol = 2)

for(i in 1:4){
  if(i == 1){
    cnf.mat[i, ] <- cnf[i,]
  } else{ 
    cnf.mat[i, ] <- cnf[i, ] + est[1] }
}
rownames(cnf.mat) <- rownames(cnf)

# put in desired order for plot
cnf.mat <- cnf.mat[c(2, 3, 1, 4), ]

#inspect
cnf.mat

# ===== STRIPCHART =====

par(bty = "l", mar = c(5, 6, 4, 6))
plot.new()
plot.window(xlim = c(.1, 1.2), ylim = c(-.1, .8) )

# set up plot
pop <- table(cldat$population)
pop <- pop[c(2, 3, 1, 4)] # reorder for plot
at <- c(.25, .55, .8, .95)
bg <- c("#F5FF0A", "#ff7880", "#F50000", "#F50000")
pch <- c(24, 4, 23, 21)

# stripchart of color values by population
for(i in 1:4){
  if (i != 2) {
    points(x = jitter(rep(at[i], pop[i]), amount = .03), 
           y = cldat$log.rg[ which( cldat$population == names(pop)[i])], 
           pch = pch[i], bg = bg[i]) 
   } else { points(x = jitter(rep(at[i], pop[i]), amount = .05), 
                y = cldat$log.rg[ which( cldat$population == names(pop)[i])], 
                pch = pch[i], col = bg[i], lwd = 2)
  }
}

# means
points(x = at + 0.08, y = means, pch = 16)

# now add CIs from cnf.mat
for(i in 1:4){
  lines(x = c(at[i] + .08, at[i] + .08), y = c(cnf.mat[i,1], cnf.mat[i,2]))
  lines(x = c(at[i] + 0.07, at[i]+ 0.09), y = c(cnf.mat[i,1], cnf.mat[i,1]))
  lines(x = c(at[i] + 0.07, at[i]+ 0.09), y = c(cnf.mat[i,2], cnf.mat[i,2]))
}


# garnishments
axis(1, lwd=0,lwd.ticks=0, font=3,labels=c("A. flavescens", "A. formosa"), 
     at=c(.25, .9), line = 0, cex.axis = 1.3)
axis(1, lwd=0,lwd.ticks=0, font=1, labels=c("hybrid"), at = .55, 
     line = 0, cex.axis = 1.3)
axis(2, at = seq(-.1, .8, .1), las = 2, cex.axis = 1.3)
mtext(2, text = "Relative red vs. green reflectance (log R/G)", cex = 1.3, line = 3.5)
box()

# ===== STOP =====

