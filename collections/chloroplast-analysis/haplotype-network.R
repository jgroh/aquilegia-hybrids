require(haplotypes)
require(plotrix) 

# ===== Comments =====

# The purpose of this script is to read in aligned plastid sequences and output a haplotype network
# using the package 'haplotypes'. Data are concatenated rps16 and trnTL sequences from 
# Aquilegia formosa, A. flavescens, and hybrids. 

# ===== Fix read.fas() bug =====

read.fas<-function(file){
  fas<-scan(file, what="char",quote="", sep="\n", strip.white=TRUE, quiet= TRUE)
  
  l<-length(fas)
  if(!l) stop(paste(file,"file is empty."))
  heads<-grep(">",fas)
  if(!length(heads)) stop("invalid fasta format, description lines of each sequence must start with a '>' (greater-than) symbol")
  tes<- diff(heads)
  test1<-any(tes==1)
  if(test1) stop("invalid fasta format, description lines of each sequence must start with a '>' (greater-than) symbol and file must not contain empty sequences")
  seqnames<-gsub(">","",fas[heads])
  len<-length(seqnames)
  seqlist<-vector("list", len)
  names(seqlist)<-seqnames
  for(i in 1:(len-1))
  {
    seqlist[[i]]<-strsplit(paste(fas[(heads[i]+1):(heads[i+1]-1)],collapse=""),split="")[[1]]
  }
  seqlist[[i+1]]<-strsplit(paste(fas[(heads[i+1]+1):(l)],collapse=""),split="")[[1]]
  as.dna(seqlist)
  
}


# ===== Read Data =====

#this generates a warning message about invalid characters, can be ignored
fast <- read.fas("collections/chloroplast-analysis/Aquilegia-concat-rps16-trnTL-ALL-AoBP-samples.fas")

# index by species affinity
fo <- grep("WG|RL|PC|PB", names(fast))

fl <- grep("BL|MK", names(fast))
hy <- grep("PP|ST", names(fast))

# check sample sizes
length(fo)
length(fl)
length(hy)

#create vector with species affinities
pops <-vector()
pops[fo] <- "fo"
pops[fl] <- "fl"
pops[hy] <- "hy"
pops

# determine haplotypes
h <- haplotypes::haplotype(fast)

# make table of haplotype groupings
g <- grouping(h, pops)$hapmat

# ===== Plot Haplotype Network ===== 

dev.off()
par(mar = c(3,3,3,3), xpd = NA)
plot.new()

# statistical parsimony network. generates an error but OK
p <- parsimnet(fast,indels="sic",prob=.95)

#plot parsimony network scaffold
set.seed(1)
vcoor <- plot(p, interactive = F, vertex.cex = .4, 
              vertex.col =1, displaylabels = FALSE)


#total number of haplotypes
nhap <- p@nhap

vcoor <- vcoor[1:nhap,]

cols <- c("#F5FF0A99","#FF0000", "#FF788099")

#define radius to be proportional to number of individuals per haplotype
rads<-apply(g,1,sum)/15

# plot white nodes to cover up background lines
labs<-paste("","",sep="")
for(i in 1: nhap){
  floating.pie(vcoor[i,1], vcoor[i,2],g[i,],edges=200,radius=rads[i],
               col= "white",startpos=0,shadow=FALSE, border = NA)
}

# plot colored nodes
for(i in 1: nhap){
  floating.pie(vcoor[i,1], vcoor[i,2],g[i,],edges=200,radius=rads[i],
               col= cols[g[i,]>0],startpos=0,shadow=FALSE)
}