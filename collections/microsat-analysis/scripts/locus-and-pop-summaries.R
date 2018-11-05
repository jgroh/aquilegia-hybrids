library(adegenet)
library(hierfstat)
library(poppr)

source('collections/microsat-analysis/scripts/microsat-src.R')

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
