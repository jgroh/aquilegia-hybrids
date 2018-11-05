# ===== READ DATA, INDEX =====

nek <- read.csv("phenotype-analysis/data/nectar.csv", stringsAsFactors = F)

#check population sample sizes
table(nek$population)

# order populations in order of decreasing mean nectar concentration
nek$population <- factor(nek$population, 
                         levels = levels(as.factor(nek$population) )[order( tapply( 
                           nek$sucrose.concentration, nek$population, mean) ) ] )

# ===== FIG-3B =====

# fit models
z <- lm(sucrose.concentration ~ 1 + population, data = nek)

#obtain estimates
est <- z$coefficients
means <- est + c(0, rep(est[1], 3))

cnf <- confint(z)
cnf[-1, ] <- cnf[-1, ] + est[1]
cnf

# reorder objects for plotting
means <- means[c(3, 4, 2, 1)]
cnf <- cnf[c(3, 4, 2, 1), ]


# ===== RESULTS FROM Bacon, 2010 ===== 

p.fl <- 44.15 # reported estimate of mean concentration
n.fl <- 33 # reported sample size
se.fl <- 1.43 # reported standard error

# obtain critical value
cv <- abs( qt( .025, df = n.fl - 1) )

# calculate 95% confidenc interval
CI.fl.up <- p.fl + cv*se.fl 
CI.fl.low <- p.fl - cv*se.fl 


# ===== STRIPCHART =====

dev.off()
par(bty = "l", mar = c(5, 6, 4, 6))
plot.new()
plot.window(xlim = c(.1,1.2), ylim = c(0, 50) )

pop <- table(nek$population)[c(3,4,2,1)]
at <- seq(.2,1, length.out = 5)
bg <- c("#F5FF0A", "#ff7880", "#F50000", "#F50000")
pch <- c(24,4,23,21)

for(i in 1:4){
  if(i != 2){
points(x = jitter(rep(at[i+1], pop[i]), amount = .03), 
       y = nek$sucrose.concentration[which(nek$population == names(pop[i]))], 
       pch = pch[i], bg = bg[i])
  } else{points(x = jitter(rep(at[i+1], pop[i]), amount = .05), 
                y = nek$sucrose.concentration[which(nek$population == names(pop[i]))], 
                pch = pch[i], col = bg[i], lwd = 2) }
}

#add mean and CI from Bacon
points(x = .2, y = p.fl, pch = 16)
lines(c(.2, .2), c(CI.fl.low, CI.fl.up))
lines(c(0.19, 0.21), c(CI.fl.low, CI.fl.low))
lines(c(0.19, 0.21), c(CI.fl.up, CI.fl.up))

# now add CIs from cnf.mat
for(i in 2:5){
  lines(x = c(at[i] + .08, at[i] + .08), y = c(cnf[i-1,1], cnf[i-1,2]))
  lines(x = c(at[i] + 0.07, at[i]+ 0.09), y = c(cnf[i-1,1], cnf[i-1,1]))
  lines(x = c(at[i] + 0.07, at[i]+ 0.09), y = c(cnf[i-1,2], cnf[i-1,2]))
}

# means
points(x = at[-1] + 0.08, y = means, pch = 16)

# garnishments
axis(1, lwd=0,lwd.ticks=0, font=3,labels=c("A. flavescens", "A. formosa"), 
     at=c(.3,.9), line = 0, cex.axis = 1.3)
axis(1, lwd=0,lwd.ticks=0, font=1,labels=c("hybrid"), at=.6, 
     line = 0, cex.axis = 1.3)
box()

axis(2, las = 2, cex.axis = 1.3)
mtext(2, text = "% sucrose concentration of nectar", cex = 1.3, line = 3.5)

# ===== STOP =====
