library(poppr)
library(ape)
library(phangorn)

source("collections/microsat-analysis/scripts/microsat-src.R")

# ===== Principal Coordinate Analysis =====

b <- bruvo.dist(str, replen = c(3,3,3,2,3,2,3,2,2,3,5))

# principal coordinate analysis
pco <- pcoa(b, correction = "lingoes")

dev.off()
barplot(pco$values[,1]) #this shows all eigenvalues, uncorrected

#set up PCoA ordination biplot
x <- pco$vectors.cor[, 1]
y <- pco$vectors.cor[, 2]
ylim <- range(y) + c(-.1, .1)
xlim <- range(x) + c(-.1, .1)
pch <- c(22, 24, 4, 21, 23)
bg <- c("#F5FF0A", "#F5FF0A", "black", "#F50000", "#F50000")
col <- c("black", "black", "#FF3399","black", "black")
lwd <- c(1, 1, 2, 1, 1)

# plot results
dev.off()
plot.new()
par(bty = "l")
plot(x,y, ylim = ylim, xlim = xlim, 
     lwd = lwd[pops], col = col[pops], pch = pch[pops], bg = bg[pops],
     xlab = "PCo 1", ylab = "PCo 2")


# ===== Tree Visualization =====

# visualize results as tree
b.mat <- as.matrix(b)

#plot nj tree 
dev.off();plot.new()
plot(nj(b), type = "unrooted",  no.margin=F, lab4ut = "axial", show.tip.label = FALSE)
tiplabels(bg = bg[pops], col = col[pops], pch = pch[pops], lwd = lwd[pops])

#plot upgma tree
plot(upgma(b),type = "unrooted",  no.margin=F, lab4ut = "axial", show.tip.label = FALSE)
tiplabels(bg = bg[pops], col = col[pops], pch = pch[pops], lwd = lwd[pops])


#examine nj tree with bootstrap values at nodes
tree <- bruvo.boot(str, replen = c(3,3,3,2,3,2,3,2,2,3,5), tree = nj, showtree = FALSE)
tree$tip.label <- rep(" ", 72)

par(mar= c(0,3,0,3))
plot(tree,  type = "radial", no.margin=F, lab4ut = "axial", align.tip.label = TRUE, use.edge.length = T,
     show.tip.label = T, show.node.label = T, node.pos = 1)
tiplabels(bg = bg[pops], col = col[pops], pch = pch[pops], lwd = lwd[pops])

# ===== Stop