# ===== Read & Index Data File =====

source('collections/microsat-analysis/scripts/microsat-src.R')

# in 'microsat-src.R' data is read in and indexed, and functions for cross-validation are defined.
# sample sizes are returned when sourcing this script. 
# population labels refer to sampling sites. 1=Beehive Lakes, 2=Mt. Kobau (Aquilegia flavescens), 
# 3=Porcupine Ridge (hybrid); 4=Robert's lake, 5=Clearwater (A. formosa).

# ===== Prepare for DAPC Analysis =====

# define training and prediction data sets
train <- str[-PP, ] #training set does not contain hybrids
supp <- str[PP, ] #supplementary data set, to be predicted

# before running the function function defined below, make sure level in 
# clusters correspond to  SPECIES, not POPULATION
train$pop

# define parental populations as a priori clusters
clust <- train$pop

# what is the optimal number of PC axes to retain?
DAPC <- dapc(train, clust, n.pca = 100, n.da = 1, prior = c(0.5,.5))
npca <- optim.a.score(DAPC)$best
npca

npca <- 5


# ===== Execute DAPC =====

DAPC <- dapc(train, clust, n.pca = npca, n.da = 1, prior = c(0.5,.5))

# what is the misclassification rate?
# err <- DapcError(npca, 1000) 
# mean(err) #0.04023256

#obtain scores for parental specimens along discriminant axis
sc <- DAPC$ind.coord

#predict hybrid scores
pred <- predict(DAPC, supp)
hyb.sc <- pred$ind.scores


# ===== AoBP Figure 4B: Density Plot of Genotype Discriminant Scores =====

dev.off()
par(bty = "l", mar = c(5,6,4,5))
plot.new()
plot.window(xlim = range(sc) + c(0.1,0.3), ylim = c(0,1))

# add densities for A. flavescens
x.fla <- seq(from = min( sc[ grep("BL|MK", rownames(sc)) ]), to = max( sc[grep("BL|MK", rownames(sc)), 1]), length.out = 512)
y.fla <- density( sc[ grep("BL|MK", rownames(sc)), ])$y
polygon(x.fla, y.fla, col = rgb(245/255,1,10/255,0.6), lwd = 2)

x.fo <- seq(from = min( sc[grep("RL|WG", rownames(sc)), 1]), to = max( sc[grep("RL|WG", rownames(sc)), 1]), length.out = 512)
y.fo <- density( sc[ grep("RL|WG", rownames(sc)) ] )$y
polygon(x.fo, y.fo, col = rgb(1,0,0,0.6), lwd = 2)

hy.coords <- pred$ind.scores
x.hy <- seq(min(hy.coords), to=max(hy.coords),
            length.out = 512)
y.hy <- density(hy.coords)$y
polygon(x.hy, y.hy, col = rgb(1,120/255,128/255,0.6), lwd = 2)

# overlay rugs overtop of density plot to show data

rug(sc[grep("RL|WG", rownames(sc)), 1], lwd = 3, line = -10, col = "#F50000")
rug(hyb.sc, col = "#ff7880", lwd = 3, line = -11)
rug(sc[grep("BL|MK", rownames(sc)), 1], lwd = 3, line = -12, col = "yellow3")


axis(2, at = c(0,.1,.2,.3,.4,.5), labels = TRUE, las = 2, cex = .5, cex.axis = 1.3)
axis(1, line = -.5, cex.axis = 1.3, at = -4:4)

mtext(text = "Density", side = 2, line = 3, adj = .25, cex = 1.5)
mtext(text = "Microsatellite genotype discriminant axis", side = 1, line =2, cex = 1.5)

# ===== Stop