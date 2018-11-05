library(MASS)
library(binom)

# ===== DEFINE FUNCTIONS =====

TraitSummary <- function(x){
  # This will give summaries of individual traits for putting into table format by species
  
  n <- length(which( !is.na(x)))
  z <- lm(x ~ 1)
  s <- summary(z)$coefficients
  mn <- round(s[1], 2)
  se <- round(s[2], 2)
  ci.lower <- round(confint(z)[1], 2)
  ci.upper <- round(confint(z)[2], 2)
  summ <- c(n, mn, se, ci.lower, ci.upper)
  names(summ) <- c("n", "Mean", "SE", "CI.lower", "CI.upper")
  summ
}


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

# ===== READ DATA, REFORMAT & INDEX FOR DOWNSTREAM ANALYSES =====
hb <- read.csv("herbarium/data/herbarium-data.csv", na.strings = c("","NA"), stringsAsFactors = FALSE)

#subset columns for analysis
hb.sub <- hb[, c("species.label","corolla.width", "spur.length","blade.length","blade.width","sepal.length","sepal.width","anther.exsertion","blade.cleft")]

#change categorical "y, n" variable to "1, 0"
a <- hb.sub$blade.cleft
for(i in 1:nrow(hb.sub) ){
  if( is.na(a[i]) ){
    a[i] <- NA
  } else if (!is.na(a[i]) & (a[i] == "Y" | a[i] == "1")){
    a[i] <- 1 
  } else if (!is.na(a[i]) & (a[i] == "N" | a[i] == "0")){
    a[i] <- 0
  }
}
hb.sub$blade.cleft <- as.numeric(a)

#create row indices for species
hb.fo <- which(hb.sub$species.label == "formosa")
hb.fl <- which(hb.sub$species.label =="flavescens")

# subset whole data set (not only individuals with complete traits) by species
fo.dat <- hb.sub[hb.fo, ]
fl.dat <- hb.sub[hb.fl, ] 

#get indices of rows containing any missing values
na <- unique( which(is.na(hb.sub), arr.ind = TRUE)[ ,1]) 

#subset data to include only individuals with complete measurements
tr.compl <- hb.sub[-na, ]
ind.compl <- as.numeric(rownames(tr.compl)) #index for subsetting in downstream analysis

#by species
fo.compl <- which(tr.compl[, 1] == "formosa")
fl.compl <- which(tr.compl[, 1] == "flavescens")

  
# ===== LDA OF FLORAL MORPHOLOGY =====

# construct model
LDA.herb <- lda(species.label ~ . , data = tr.compl)

# obtain phenotype scores for downstream analysis
LDA.sc <- predict(LDA.herb)$x

# loadings of traits onto discriminant axis
LDA.loadings <- LDA.herb$scaling

# order by decreasing absolute value
ldngs.ind <- order(abs(LDA.loadings), decreasing=TRUE) #order traits in order of loadings onto RDA1
ldngs.ord <- LDA.loadings[ldngs.ind, ]
ldngs.ord

# ===== FLORAL MORPHOLOGY LDA MISSCLASSIFICATION RATE =====

mis <- vector()

#for(i in 1:1000){
  # calculate the average misclassification rate over 1000 iterations 
#  xval <- Vlda(5, species.label ~ . , data = tr.compl, cl = tr.compl$species.label)
#  mis[i] <- sum(xval[row(xval) != col(xval)]) / sum(xval)
#}

mean(mis) #0.07574172


# ===== AoBP TABLE 2: TRAIT SUMMARIES =====

# apply TraitSummary to each trait

# 1. lamina length
tapply(hb.sub$blade.length, hb.sub$species.label, TraitSummary)

# 2. lamina width
tapply(hb.sub$blade.width, hb.sub$species.label, TraitSummary)

# 3. corolla width
tapply(hb.sub$corolla.width, hb.sub$species.label, TraitSummary)

# 4. spur length
tapply(hb.sub$spur.length, hb.sub$species.label, TraitSummary)

# 5. sepal.width
tapply(hb.sub$sepal.width, hb.sub$species.label, TraitSummary)

# 6. anther exsertion
tapply(hb.sub$anther.exsertion, hb.sub$species.label, TraitSummary)

# 7. sepal length
tapply(hb.sub$sepal.length, hb.sub$species.label, TraitSummary)

# 8. Cleft laminae
hb.bl <- hb.sub$blade.cleft

# for A. formosa
x1 <- sum(hb.bl[hb.fo], na.rm = TRUE)
n1 <- length(which( !is.na(hb.bl[hb.fo])))
b <- binom.agresti.coull(x = x1, n = n1)
b
p <- b[4] 
se <- sqrt(p*(1-p)/b[3])
se

#method   n     mean     lower    upper
# agresti-coull 36 97 0.371134 0.2815011 0.470585

# for A. flavescens (p = 0)
length(which( !is.na( hb.bl[hb.fl])))


