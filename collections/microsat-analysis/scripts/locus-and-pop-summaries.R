library(adegenet)
library(pegas)
library(hierfstat)

# ===== Read & Index Data File =====
str <- read.structure(file = "collections/microsat-analysis/data/Aquilegia-final-micros-AoBP.str", 
                      n.ind = 72, n.loc = 11, onerowperind = T, col.lab = 1, col.pop = 2,
                      col.others = 3, row.marknames = 1, NA.char = "-9")

id <- rownames(str$tab)

BL <- grep("BL", id)
MK <- grep("MK", id)
PP <- grep("PP|St", id)
RL <- grep("RL", id)
WG <- grep("WG", id)

pops <- vector()
pops[BL] <- 1
pops[MK] <- 2
pops[PP] <- 3
pops[RL] <- 4
pops[WG] <- 5

# ===== Locus Summaries =====

# general summaries, number of alleles, Hobs and Hexp
summary(str)
mean(summary(str)$loc.n.all)
Ho <- summary(str)$Hobs
Ho
mean(Ho)
He <- summary(str)$Hexp
He
mean(He)


#allele bp ranges
sapply(str@all.names, function(x)range(as.numeric(x)))

# ===== Population Summaries =====

str@pop <- as.factor(pops)
hf <- genind2hierfstat(str)
bs <- basic.stats(hf)
bs

# Fis (Inbreeding coefficient)
bs <- basic.stats(hf)
fis <- colMeans(bs$Fis)
fis
boot.ppfis(hf, 1000)

# Ho (Observed Heterozygosity)
Ho.pop <- colMeans(bs$Ho)
Ho.pop

# Hs (gene diversity)
Hs.pop <- colMeans(bs$Hs)
Hs.pop

# H (Shannon Diversity)
tbl <- mlg.table(str)
diversity_stats(tbl)
diversity_ci(str, raw = F)

# ===== STOP