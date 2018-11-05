library(mapproj)

# ===== Read Data, Index =====

d <- read.csv("herbarium/data/herbarium-data.csv", stringsAsFactors = FALSE)

lon <- d[,"Geo_LongDegree"]
lat <- d[,"Geo_LatDegree"]

fl <- which(d$species.label=="flavescens")
fo <- which(d$species.label=="formosa")


# ===== Map Making =====

dev.off()
par(oma=c(4,1,1,1))
par(fig = c(0,1,0,1), oma=c(0,0,0,0), mar=c(0,0,0,0))

map(database= "world", ylim=c(45,90), xlim=c(-160,-50),
    col=rgb(.5,.5,.8, 0.6), fill=TRUE, border = "white",lwd = .5,
    projection="gilbert", orientation= c(90,0,225))

coord.all <- mapproject(lon, lat, proj="gilbert", orientation=c(90, 0, 225))
coord.fla <- mapproject(lon[fl], lat[fl], proj="gilbert", orientation=c(90, 0, 225))
coord.fo <- mapproject(lon[fo], lat[fo], project="gilbert",orientation=c(90,0,225))

# ===== Species Distributions

#All flavescens
points(coord.fla, cex = .8, pch = 16, col = rgb(245/255,1,10/255, .5))

#All formosa
points(coord.fo, cex = .8, pch = 16, col = rgb(1,0,0, .5))

# ===== Sampling Locations =====

#WG 
points(mapproject(-120, 51.90, proj = "gilbert", orientation = c(90, 0, 225)),
       pch = 23, cex = 1.2, lwd =.5, bg = "#F50000",col="black")

#PP
points(mapproject(-121.828345, 51.112431, proj = "gilbert", orientation = c(90, 0, 225)),
       pch = 16, cex = 1.2, col = "#ff7880", lwd = 1)
points(mapproject(-121.828345,51.112431, proj = "gilbert", orientation = c(90,0,225)),
       pch = 13, cex = 1.2, lwd = 1)

#MK
points(mapproject(-119.666, 49.111, proj = "gilbert", orientation = c(90,0,225)),
       pch=24, cex=1, lwd = 1, bg = "#F5FF0A", col = "black")

#BL
points(mapproject(-116.654, 48.655, proj="gilbert", orientation = c(90,0,225)),
       pch = 22, cex = 1, lwd = 1, bg = "#F5FF0A", col = "black")

#RL
points(mapproject(-125.546, 50.215, proj = "gilbert", orientation = c(90,0,225)),
       pch=21, cex = 1, lwd = 1, bg = "#F50000", col="black")

