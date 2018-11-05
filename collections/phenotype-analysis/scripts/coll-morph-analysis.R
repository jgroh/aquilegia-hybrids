require(MASS)
library(binom)

# ===== COMMENTS =====

# This script conducts floral morphometric analyses done on natural populations of
# Aquilegia formosa, A. flavescens, and hybrids. 
# Results of these analyses are published in 

# ===== READ AND FORMAT DATA =====

d <- read.csv("collections/phenotype-analysis/data/collections-phenotypes.csv", stringsAsFactors = FALSE)

# reformat y/n variable
d$blade.cleft[which(d$blade.cleft == "y")] <- 1
d$blade.cleft[which(d$blade.cleft == "n")] <- 0
d$blade.cleft <- as.numeric(d$blade.cleft)

# index by species
hy <- which(d$species == "hybrid") 
fo <- which(d$species == "A. formosa")
fl <- which(d$species == "A. flavescens")

# training set for LDA
prnt <- d[-hy, c(2, 4:11)] 

#remove missing data rows
na.pheno <- unique(which(is.na(prnt), arr.ind = TRUE)[, 1])
prnt <- prnt[-na.pheno, ]

fo.prnt <- which(prnt$species == "A. formosa") 
fl.prnt <- which(prnt$species == "A. flavescens") 

# CROSS VALIDATION FUNCTION -----

Vlda = function(v, formula, data, cl){
  # This will split the data randomly into v groups, 
  # for each of 1:v, fit model based on training data (-v) and predict supp data (v),
  # then determine the missclassification rate.

  grps = cut(1:nrow(data), v, labels=FALSE)[ sample(1:nrow(data)) ]
  
  pred = lapply(1:v, function(i, formula, data){
    omit = which(grps == i)
    z = lda(formula, data=data[-omit, ])
    predict(z, data[omit, ])
  }, formula, data)
  
  wh = unlist(lapply(pred, function(pp)pp$class))
  table(wh, cl[order(grps)]) 
}

# FLORAL PHENOTYPE LDA --------------------------------------------------------

# construct model
z <- lda(species ~ ., data = prnt, prior = c(.5,.5))
z$scaling #gives the coefficients of the linear combination of traits
sc <- predict(z)$x

# apply model to hybrid phenotypes
hyb <- predict(z, d[hy, c(2,4:11)], prior = c(.5,.5))
hyb.sc <- hyb$x


# FL LDA MISSCLASSIFICATION RATE -----


mis <- vector()
for(i in 1:1000){
 # calculate the average misclassification rate over 1000 iterations 
  
 xval <- Vlda(5, species ~ . , data = prnt, cl = prnt$species)
 mis[i] <- sum(xval[row(xval) != col(xval)]) / sum(xval)
}
mean(mis) #0.0002564103

# check error rate with smaller sample size (~= to that in genetic analysis)
fl.sub <- prnt[grep("flav", prnt$species), ]
fo.sub <- prnt[grep("form", prnt$species), ]

#subset phenotype data set
set.seed(123)
prnt.sub <- rbind(fl.sub[sample(1:nrow(fl.sub), 15), ], fo.sub[sample(1:nrow(fo.sub), 25), ])

# fit lda model
z.sub <- lda(species ~ ., data = prnt.sub, prior = c(.5,.5))

mis.sub <- vector()
for(i in 1:1000){
  xval <- Vlda(5, species ~ . , data = prnt.sub, cl = prnt.sub$species)
  mis.sub[i] <- sum(xval[row(xval) != col(xval)]) / sum(xval)
}
mean(mis.sub) # 0.00075



# AoBP FIGURE 4A: DENSITY PLOT OF FLORAL MORPHOLOGY LDA SCORES ---------

dev.off()
plot.new()
par(bty = "l", mar = c(5,6,4,6))
plot.window(xlim = range(sc) + c(0.1,0.3), ylim = c(0,1))

# show densities

# A. flavescens
x.fla <- seq(min( sc[fl.prnt] ), to=max( sc[fl.prnt] ), length.out = 512)
y.fla <- density( sc[fl.prnt] )$y
polygon(x.fla, y.fla, col = rgb(245/255,1,10/255,0.6), lwd = 2)

#. A. formosa
x.fo <- seq(min( sc[fo.prnt] ), to=max( sc[fo.prnt] ), length.out = 512)
y.fo <- density(  sc[fo.prnt] )$y
polygon(x.fo, y.fo, col = rgb(1,0,0,0.6), lwd = 2)

# hybrids
x.hy <- seq(min(hyb.sc), to=max(hyb.sc), length.out = 512)
y.hy <- density(hyb.sc)$y
polygon(x.hy, y.hy, col = rgb(1,120/255,128/255,0.5), lwd = 2)

# overlay rugs to show data
rug(sc[fo.prnt, 1], lwd = 3, line = -10, col = "#F50000")
rug(hyb$x, col = "#ff7880", lwd = 3, line = -11)
rug(sc[fl.prnt, 1], lwd = 3, line = -12, col = "yellow3")

axis(1, line = -.5, cex.axis = 1.3)
axis(2, at = c(0,.1,.2,.3,.4,.5), labels = TRUE, las = 2, cex.axis = 1.3)

mtext(text = "Density", side = 2, line = 3, adj =.25, cex = 1.5)
mtext(text = "Floral morphology discriminant axis", side = 1, line =2, cex = 1.5)


# ===== CORRELATION OF LOG(R/G) WITH FLORAL MORPHOLOGY LDA =====

# create log.rg variable
log.rg <- log(d$red.mean/ d$green.mean)

# subset to hybrids
log.rg.hyb <- log.rg[hy]

# calculate correlation
cor(log.rg.hyb, hyb.sc, use = "complete.obs")

# use permutation to assess significance
set.seed(123)
iter <- 9999
x <- log.rg.hyb
cnt <-0
cor.obs <- cor(log.rg.hyb, hyb.sc, use="complete.obs")

for(i in 1:iter){
  z <- x
  z <- z[sample(1:length(z))]
  cor.est <- cor(x=z, y= hyb.sc, use="complete.obs")
  
  if (abs(cor.est) >= abs(cor.obs)) { 
    cnt <- cnt + 1
    }
  
}
cnt
p.value <- round(as.numeric((cnt +1 )/(iter + 1)),digits = 4)
p.value #0.109

# ===== INDIVIDUAL TRAIT DISTRIBUTIONS =====

species <- factor(d$species, levels(factor(d$species))[c(1, 3, 2)])
bg <- c("yellow3", "#ff7880", "#F50000")

traits <- c("blade.length", "blade.width", "corolla.width", "spur.length", "sepal.width", "anther.exsertion", "sepal.length")
y.axis.labels <- c("Blade length (cm)", "Blade width (cm)", "Corolla width (cm)", "Spur length (cm)", "Sepal width (cm)", "Anther exsertion (cm)", "Sepal length (cm)")

# create multipanel plot showing trait distributions for each of 7 traits

par(mfrow=c(4,2),bty="l",mar=c(2,5,2,2),oma=c(1,0,0,0))

for(i in 1:7){
stripchart(d[[traits[i]]] ~ species, vertical = TRUE, method = "jitter",
           pch = 21, col = bg, lwd = 2, at = c(1,2,3), xlim = c(.5, 3.5), 
           xaxt = "n", yaxt = "n", ylab = y.axis.labels[i], cex.lab = 1.2)
axis(2, las = 2, cex.axis = 1.2)

z <- lm(d[[traits[i]]] ~ 1 + d$species)
est <- z$coefficients[c(1, 3, 2)]
means <- est + c(0, est[1], est[1])

points(x = 1:3 + .3, y = means, pch = 16)

cnf <- confint(z)
cnf[c(2, 3), ] <- cnf[c(2, 3), ] + est[1]
cnf <- cnf[c(1, 3, 2), ]

at <- 1:3 + .3

for(j in 1:3){
  lines(x = c(at[j], at[j]), y = c(cnf[j, 1], cnf[j, 2]))
}
}

# Proportion of individuals with cleft laminae in each species

fo.bl <- d$blade.cleft[fo]
binom.agresti.coull(x = sum(fo.bl, na.rm = TRUE), n = length(which(!is.na(fo.bl))))

hy.bl <- d$blade.cleft[hy]
binom.agresti.coull(x = sum(hy.bl), n = length(hy.bl))

# Are these proportions significantly different?
prop.test(x = c( sum(fo.bl, na.rm = TRUE), sum(hy.bl)), 
          n =  c(length(which(!is.na(fo.bl))), length(hy.bl)), 
          alternative = "two.sided" )

# ===== STOP =====


