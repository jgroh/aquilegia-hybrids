library (adegenet)
library (pegas)
library (ape)
library(poppr)
library(polysat)
library(phangorn)
library(RColorBrewer)
library(MASS)


# ===== Read & Index Data File =====

# data contains genotypes of 72 individuals for 11 microsatellite loci. 
# population labels refer to sampling sites. 1=Beehive Lakes, 2=Mt. Kobau (Aquilegia flavescens), 
# 3=Porcupine Ridge (hybrid); 4=Robert's lake, 5=Clearwater (A. formosa).

# read structure file
str <- read.structure(file = "collections/microsat-analysis/data/Aquilegia-final-micros-AoBP.str", 
                      n.ind = 72, n.loc = 11, onerowperind = T, col.lab = 1, col.pop = 2,
                      col.others = 3, row.marknames = 1, NA.char = "-9")

tab <- str@tab

BL <- grep("BL", rownames(tab))
MK <- grep("MK", rownames(tab))
PP <- grep("PP|St", rownames(tab))
RL <- grep("RL", rownames(tab))
WG <- grep("WG", rownames(tab))

pops <- vector()
pops[BL] <- 1
pops[MK] <- 2
pops[PP] <- 3
pops[RL] <- 4
pops[WG] <- 5


# ===== Define Functions =====

XvalDapc = function(v, gen, clust, n.pca){ 
  # This function is for determining the missclassification rate of a DAPC model.
  #The output is a table containing correct vs. incorrect sortings
  # v = number partitions of data set
  # gen = genind object
  # clust = factor corresponding to a priori groupings (species in this case)
  # n.pca = number of PC axes to be retained in DAPC
  
  s <- sample(1:nrow(tab(gen))) # index for randomly subsetting split data set, also to index "clust"
  
  clu <- clust[s] 
  
  grps = cut(1:nrow(tab(gen)), v, labels=FALSE)[s] #split data 5 sets, subset by random index
  
  pred = lapply(1:5,function(v, gen, n.pca){ #apply to each subset 1-5
    
    omit <- which(grps == v) #for 1:v, function will omit 1...v 
    z <- dapc(gen[-omit, ], grp=clu[-omit], n.pca = n.pca, n.da = 1) #function will run dapc on remaining data
    # (n.pca will later be indexed in a for() loop)
    
    predict(z,gen[omit, ])
  }, gen, n.pca) #function requires arguments "genind" and "n.pca"
  
  #pred is a list of length(5)
  # pred[[1]] will contain predicted values for the data omitted by the first index
  # which(grps == 1) will correspond to the IDs in pred[[1]]
  
  wh = unlist(lapply(pred, function(pp)pp$assign)) #unlists list of species assignment
  table(wh, clust[order(grps)]) #very important that this is clust, not clu
}

DapcError <- function(npca, n){
  # input number of pca axes retained, and number of iterations
  
  a <- vector() #vector for holding error rate across 50 trials
  
  for(j in 1:n){
  x <- XvalDapc(5, train, clust, n.pca = npca)
  a[j] <- sum(x[row(x) != col(x)]) / sum(x)
  }
  
return(a)
} 


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

#################################################################################

# EXECUTE DAPC -----

DAPC <- dapc(train, clust, n.pca = npca, n.da = 1, prior = c(0.5,.5))

# what is the misclassification rate?
err <- DapcError(npca, 1000) #0.04023256
mean(err)


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