library(geosphere)

# ===== SOURCE MORPHOLOGY ANALYSIS =====
source('herbarium/scripts/herb-morph-analysis.R')

# ===== CALCULATE RANGE CENTROIDS ===== 
# get coordinates for specimens with complete trait data
lon <- hb$Geo_LongDegree[ind.compl]
lat <- hb$Geo_LatDegree[ind.compl]

#convert from degrees to radians
lon.rad <- lon*(pi/180)
lat.rad <- lat*(pi/180)

#convert to cartesian coordinates
x <- cos(lon.rad) *cos(lat.rad)
y <- cos(lat.rad)*sin(lon.rad)
z <- sin(lat.rad)

#take average of coordinates
x.fl <- mean( x[fl.compl])
y.fl <- mean( y[fl.compl])
z.fl <- mean( z[fl.compl])

#convert to radians
hyp.fl <- sqrt(x.fl^2 + y.fl^2 + z.fl^2) 
theta.fl <- asin(z.fl/ hyp.fl) #inclination
phi.fl <- atan2(y.fl,x.fl) #rotation

#convert to degrees
latcen.fl <- theta.fl*(180/pi)
loncen.fl <- phi.fl*(180/pi)

#repeat for A. formosa
x.fo <- mean( x[fo.compl])
y.fo <- mean( y[fo.compl])
z.fo <- mean( z[fo.compl])
hyp.fo <- sqrt(x.fo^2 + y.fo^2 + z.fo^2)

theta.fo <- asin(z.fo/ hyp.fo)
phi.fo <- atan2(y.fo,x.fo)

latcen.fo <- theta.fo*(180/pi)
loncen.fo <- phi.fo*(180/pi)

###################################################################################
# CALCULATE DISTANCE TO CENTROIDS OF EACH SPECIES ----- 

# calculate distance of A. flav specimens to centroid of A. formosa, centroid of A. flav
fl.dist.from.fo <- vector()
fl.dist.from.fl <- vector()

for(i in 1:length(fl.compl) ){
  # for each A. flavescens specimen, calculate dist. from heterospecific centroid
  fl.dist.from.fo[i] <- distm( c(loncen.fo, latcen.fo), c(lon[fl.compl][i], lat[fl.compl][i]), fun=distHaversine)
}

for(i in 1:length(fl.compl)){
  # for each A. flavescens specimen, calculate dist. from conspecific centroid
  fl.dist.from.fl[i] <- distm(c(loncen.fl, latcen.fl), c(lon[fl.compl][i], lat[fl.compl][i]), fun=distHaversine)
}

# convert to km
fl.dist.from.fo <- fl.dist.from.fo/1000 
fl.dist.from.fl <- fl.dist.from.fl/1000 


fl.dist.ratio <- log(fl.dist.from.fl/fl.dist.from.fo)
fl.dist.ratio 
# should be mostly negative, meaning the ratio is less than 1 
# i.e. specimens are generally closer to the centroid of their own species

# repeat calculations for A. formosa
fo.dist.from.fl <- vector()
fo.dist.from.fo <- vector()

for(i in 1:length(fo.compl)){
  fo.dist.from.fl[i] <- distm(c(loncen.fl, latcen.fl), c(lon[fo.compl][i], lat[fo.compl][i]), fun=distHaversine)
}

for(i in 1:length(fo.compl)){
  fo.dist.from.fo[i] <- distm(c(loncen.fo, latcen.fo), c(lon[fo.compl][i], lat[fo.compl][i]), fun=distHaversine)
}

# convert to km
fo.dist.from.fl <- fo.dist.from.fl/1000 
fo.dist.from.fo <- fo.dist.from.fo/1000

fo.dist.ratio  <- log(fo.dist.from.fo/fo.dist.from.fl)

# combine dist.ratios to create variable used in regression
# the indexing should match the order of points from LDA in the morphology script.

dist.ratio <- vector()
dist.ratio[fo.compl] <- fo.dist.ratio
dist.ratio[fl.compl] <- fl.dist.ratio

###################################################################################
# AoBP FIG6: GEOSPATIAL ANALYSIS OF FLORAL MORPHOLOGY ----- 


#plot linear discriminant score against dist.ratio

col <- as.character(tr.compl[, 1])
col[which(col == "flavescens")] <- "#F5FF0A"
col[which(col == "formosa")] <- "#F50000"

# set up plot
dev.off()
par(mar=c(7,5,2,2), bty = "l")
plot.new()
plot.window( xlim = range(dist.ratio) + c(0, .6), ylim = range(LDA.sc) + c(-.4,.3) )

# add points
points(LDA.sc ~ dist.ratio,  bg = col, ylab = "", xlab = "", yaxt = "n", pch = 21)


# garnishments
axis(2, las = 2, lwd=0, lwd.ticks=1, cex.axis = 1.2)
axis(1, lwd=0,lwd.ticks=1, cex.axis = 1.2)
box()

mtext(side = 1, expression(log(frac("Distance from conspecific centroid", 
                                    "Distance from heterospecific centroid"))), line = 5, cex = 1.5)

mtext(side = 2, text = "Floral morphology discriminant axis", line = 3, cex = 1.5)

# create data frames for calculating regression models and confidence bands
fldat <- cbind.data.frame(a = LDA.sc[fl.compl], b = dist.ratio[fl.compl])
fodat <- cbind.data.frame(a = LDA.sc[fo.compl], b = dist.ratio[fo.compl])

# add regression line and confidence bands for A. flavescens
z.fl <- lm(a ~ b, data = fldat)
xpts.fl <- range(fldat$b, na.rm = T)
yhat.fl <- predict(z.fl, newdata = data.frame(b = xpts.fl))
newx <- seq(min(fldat$b, na.rm = TRUE), max(fldat$b, na.rm = TRUE), length.out=100)
preds <- predict(z.fl, newdata = data.frame(b=newx),  interval = 'confidence')

polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = rgb(190/255,190/255,190/255, .6), border = NA)
lines(yhat.fl ~ xpts.fl, lwd =2)

# repeat for A. formosa
z.fo <- lm(a ~ b, data = fodat)
xpts.fo <- range(fodat$b, na.rm = T)
yhat.fo <- predict(z.fo, newdata = data.frame(b = xpts.fo))
newx <- seq(min(fodat$b, na.rm = TRUE), max(fodat$b, na.rm = TRUE), length.out=100)
preds <- predict(z.fo, newdata = data.frame(b=newx),  interval = 'confidence')

polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = rgb(190/255,190/255,190/255, .6), border = NA)
lines(yhat.fo ~ xpts.fo, lwd =2)

# add points again 
points(LDA.sc ~ dist.ratio, pch=21, bg = col, ylab = "", xlab = "", yaxt = "n")

# ===== STOP 