# ===== COMMENT =====

# This script will plot collection date (a proxy for flowering time) against altitude for 
# herbarium specimens of A. formosa and A. flavescens

# ===== READ DATA, INDEX =====

d <- read.csv("herbarium/data/herbarium-data.csv", stringsAsFactors = FALSE)

#obtain dates
dates <- as.Date(d$Date, format = '%Y %b %d')

# remove year info, not of direct interest
days <- as.Date(format(dates, "%b %d"), format = "%b %d")

# combine variables into data frame for each species model fitting
days.fla <- days[which(d$species.label == "flavescens")]
alt.fla <- d$altitude.corrected[which(d$species.label == "flavescens")]
fla <- cbind.data.frame(day = as.numeric(days.fla), alt = alt.fla)

days.fo <- days[which(d$species.label == "formosa")]
alt.fo <- d$altitude.corrected[which(d$species.label == "formosa")]
fo <- cbind.data.frame(day = as.numeric(days.fo), alt = alt.fo)


# ===== FIT MODEL =====

days <- as.numeric(days)
z <- lm(days ~ altitude.corrected + species.label, data = d, na.action = na.exclude)
summary(z)
anova(z) # difference between species is not significant
confint(z)

# fit model without species factor
z2 <- lm(days ~ altitude.corrected, data = d, na.action = na.exclude)
summary(z2)
anova(z2)
confint(z2)


# ===== MAKE PLOT =====

# set up plot
par(bty = "l", mar =c(5.1,6.1,4.1,2.1))

col <- as.character(d[, 3])
col[which(col == "formosa")] <- "#F50000"
col[which(col == "flavescens")] <- "#F5FF0A"

plot(days ~ d$altitude.corrected, 
     type = "n", xlab ="", ylab = "", 
     yaxt = "n", xaxt = "n", axes = F)

# add confidence band
xpts <- range(d$altitude.corrected)
yhat <- predict(z2, newdata = data.frame(altitude.corrected = xpts))

newx <- seq(min(d$altitude.corrected), max(d$altitude.corrected), length.out = 100)
preds <- predict(z2, newdata = data.frame(altitude.corrected = newx), interval = "confidence")
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = rgb(190/255,190/255,190/255, .6), border = NA)

# regression line
lines(yhat ~ xpts, lwd =2)

# overlay points
points(days ~ d$altitude.corrected,   pch = 21, bg = col)

# add garnishments
at <- as.numeric(as.Date(c("2018-04-01", "2018-05-01","2018-06-01",
                           "2018-07-01","2018-08-01","2018-09-01")))
labels <- c("Apr 1", "May 1", "Jun 1", "Jul 1", "Aug 1", "Sep 1")

mtext(side = 2, text = "Collection date", line = 4, cex = 1.2)
mtext(side = 1, text = "Altitude (m)", line = 3, cex = 1.2)
axis(2, at = at, labels = labels, las = 2, cex.axis = 1.1)
axis(1, las = 1, cex.axis = 1.2)


# ===== STOP