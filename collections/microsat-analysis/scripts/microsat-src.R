library (adegenet)

# ===== Comments =====

# data contains genotypes of 72 individuals for 11 microsatellite loci. 
# population labels refer to sampling sites. 1=Beehive Lakes, 2=Mt. Kobau (Aquilegia flavescens), 
# 3=Porcupine Ridge (hybrid); 4=Robert's lake, 5=Clearwater (A. formosa).

# ===== Read microsat .str file and index for analyses =====
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

# ===== Check Sample Sizes =====

id <- substr(rownames(str@tab), start = 1, stop = 2)
id[which(id == "St")] <- "PP"
sample.sizes <- table(id)
print(sample.sizes)

# ===== Define Functions Used in DAPC =====

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

# ===== STOP