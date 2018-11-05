source("collections/microsat-analysis/scripts/microsat-src.R")

# ===== Comments =====
# Analysis and visualization of STRUCTURE results. 
# Input files are (1) str file with microsatellite genoytypes (11 microsat loci fr 72 Aquilegia individuals from 5 pops)
# and (2) output file from CLUMPP giving mean Q-matrix from 20 replicate runs. 

# ===== Read Files =====

# The .str file containing genotypes is read in by the sourced script. 

# read CLUMPP output file
q <- read.table("collections/microsat-analysis/data/CLUMPP-files/aq_outfile")
q <- as.matrix(q[, c(6:7)])


# ===== Calculate Avg. Ancestry Per Parent in Hybrids =====

mean(q[grep("PP", id), 1]) # mean A. flavescens ancestry
mean(q[grep("PP", id), 2]) # mean A. formosa ancestry


# ===== AoBP Figure 5: Structure Ancestry Plot =====

par(mfrow = c(2,1), oma = c(10,2,.5,.5), mar = c(.3,2,.3,2))
compoplot(q, col.pal = c("#F5FF0A99", "#FF0000"), las = 2, 
          legend = FALSE, border = NA, space = 0, yaxt = "n")
axis(2, las = 2, line = -.99, cex = 1.2)
lines(x = c(0,72), y = c(0,0))
lines(x = c(0,72), y = c(1,1))
lines(x = c(72,72), y = c(0,1))

#delimit populations
consec <- c(11,23,52,58)
for(i in 1:4){
  lines(x = c(consec[i], consec[i]), y = c(0,1), lty = 2)
}

mtext(text = "Q value", side = 2, cex = 1.2, outer = TRUE, adj = .85)

# ===== Stop