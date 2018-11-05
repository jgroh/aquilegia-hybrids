# ===== Comments =====
# This script plots a comparison of hummingbird bill lengths to Aquilegia spur lengths

# ===== Read Data =====

# Aquilegia data
a <- read.csv("herbarium/data/herbarium-data.csv", stringsAsFactors = FALSE)

fo <- which(a$species.label == "formosa")
fl <- which(a$species.label == "flavescens")

# hummingbird data
h <- read.csv("herbarium/data/hummingbird-bills.csv", stringsAsFactors = FALSE)

an <- which(h$species == "anna's")
ru <- which(h$species == "rufous")
ca <- which(h$species == "calliope")


# ===== Make Plot =====

dev.off()
par(oma = c(2,4,0,0), bty = "l")
h <- h[-an, ]

with(h, boxplot(bill.length ~ species, 
                horizontal = TRUE, pch = 16, xlim = c(0,5), at = c(3,2), 
                cex.axis = 1.2, ylim = c(5,30), yaxt = "n", lty = 1,
                col = rgb(30/255,144/255,1,.6)))

boxplot(a$spur.length[fo]*10, at = 1, add = TRUE, 
       horizontal = T, col = "#F50000", pch = 16, lty = 1, axes=F)

boxplot(a$spur.length[fl]*10, horizontal=TRUE, at = 4,
        col = "#F5FF0A", pch = 16, add = TRUE, lty=1, axes = F )

# garnishments
axis(2, at = c(2,3), labels = c( "Rufous\nhummingbird", "Calliope\nhummingbird"), 
     cex.axis = 1.2, las = 2, cex.lab = 1.1)

axis(2, at = c(1,4), labels = c("A. formosa", "A. flavescens"), 
     cex.axis = 1.2, font = 3, las = 2, cex.lab = 1.1)

mtext(side = 1, text = "Spur or bill length (mm)", line = 3, cex = 1.2 )

# ===== Stop