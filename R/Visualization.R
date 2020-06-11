# visualization of my MBar data... inspired by Daniel Marquina's script



### Load the libraries
### Packages not installed? -> "Error in library(pepe) : there is no package called" ...dplyr...
library(vegan)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(wesanderson)
library("ggsunburst")

library(sunburstR)
library(circlize)
library(gplots)
library(reshape)
library(reshape2)
library(gridExtra)
library(grid)
lenght = length
### Daniel often writes lenght instead of length - maybe find+replace all later? except this lenght = length

### What sets of data are u gonna visualize?
sets<-c("All Taxa","All Taxa Quality Taxonomy", "Animalia","Animalia Quality Taxonomy") #list names of sets of results that you wanna visualize
nsets<-length(sets)
#sets2<-c("All Taxa ee1.5", "All Taxa ee1.5 Quality Alignment","All Taxa ee0.5", "All Taxa ee0.5 Quality Alignment")


### Load data
setwd("C:/Users/ntb/Documents/MBar-Data/20_01_10_Tatry-Dunaj/Tatry/X_Detours/E_Cluster_otus")
s1<-(read.csv("5_OTU_table_0.01_ZERO_023c_TaxAll.csv", stringsAsFactors = FALSE, row.names = 1))
s2<-(read.csv("5_OTU_table_0.01_ZERO_023c_TaxAllHQ.csv", stringsAsFactors = FALSE, row.names = 1))

s3<- (read.csv("5_OTU_table_0.01_ZERO_023c_TaxAnimalia.csv", stringsAsFactors = FALSE, row.names = 1))
s4<- (read.csv("5_OTU_table_0.01_ZERO_023c_TaxAnimaliaHQ.csv", stringsAsFactors = FALSE, row.names = 1))

list1<- list(s1,s2,s3,s4)
nona<-list1
for (i in 1:nsets) {
  nona[[i]][is.na(list1[[i]])]<-0
}

### What samples you have? And which columns are they in?
samcols <- 10:22# adjust to fit columns of samples read abundance
samples<-colnames(nona[[1]][,samcols]) 
nsamp<-length(samples)

####################
# tables of kingdoms and how many HQ reads they had with taxonomic assignment to lower levels
taxq = data.frame(matrix(nrow = 6, ncol = 5))
colnames(taxq) = c("Animalia","Chromista","Bacteria","Plantae","Fungi")
rownames(taxq) = c("Phylum","Class","Order","Family","Genus","Species")
for (i in 1:nrow(taxq)){
  for (j in 1:ncol(taxq)){
    taxq[i,j] = sum(colSums(nona[[2]][which(nona[[2]][,1+i]!="" & nona[[2]]$Kingdom==colnames(taxq)[j]),samcols]))
  }
}
#write.csv(taxq,"reads_with_HQTaxAs_per_Kingdom_per_taxLevel.csv")

####################
# Rarefaction curve
t<-list(0)
out<-t
col <- wes_palette("Zissou1", 7, type = "continuous")
pars <- expand.grid(col = col, stringsAsFactors = FALSE)
par(mfrow=c(1,2))
for (i in 1:nsets) {
  t[[i]] <- as.data.frame(t(nona[[i]][,samcols]))
  out[[i]] <- rarecurve(t[[i]], step=20, xlab="No. reads", ylab="No. MOTUs", main=sets[i], col=col, label =TRUE, fill = rgb(red =1,green = 1, blue = 1, alpha=0.5))
}

#######################
# BARPLOTS for richness

MOTU_sample <- as.data.frame(matrix(ncol = 3, nrow = nsets*nsamp))
colnames(MOTU_sample) <- c("Sample", "MOTUs", "Dataset")
MOTU_sample$Sample <- factor(rep(samples, nsets), levels=samples)
for (h in 1:nsets){
  MOTU_sample$Dataset[(1+nsamp*(h-1)):(nsamp*h)] <- c(rep(sets[h],nsamp))
  for (i in 1:nsamp){
    MOTU_sample$MOTUs[i+nsamp*(h-1)] <- length(which(nona[[h]][,samcols[1]-1+i]!=0))
  }
}
MOTU_sample$Dataset <- factor(MOTU_sample$Dataset, levels = sets)

ggplot(MOTU_sample, aes(Sample, MOTUs, fill=Dataset)) +
  geom_bar(stat="identity") +
  facet_wrap(~Dataset, scales = "free") + 
  scale_fill_manual(values = c(wes_palette("Darjeeling1")[1:nsets]))

###########################################
### Number of samples each MOTU is found in
p<-list(0)
sample.per.MOTU<-p
for (h in 1:nsets){
  sample.per.MOTU[[h]]<- as.data.frame(matrix(ncol = 2, nrow = length(nona[[h]][,1])))
  colnames(sample.per.MOTU[[h]]) <- c("MOTU", "inSamples")
  for (i in 1:length(nona[[h]][,1])){
    sample.per.MOTU[[h]]$MOTU[i] <- rownames(nona[[h]])[i]
    sample.per.MOTU[[h]]$inSamples[i] <- length(which(nona[[h]][i,samcols]!=0))
  }
  p[[h]] <- ggplot(dplyr::count(sample.per.MOTU[[h]], var=inSamples), aes(var, n)) +
    geom_bar(stat="identity", fill=wes_palette("Darjeeling1")[h]) + theme(legend.position = 'none') + labs(title = sets[h], x="No. samples", y="No. OTUs")
}
nCol <- floor(sqrt(length(p)))
do.call("grid.arrange", c(p, ncol=nCol))


#############################
# BARPLOTS for read abundance and phyla-OTU composition

proportions<-nona
for (h in 1:nsets) {
  for (i in samcols){
    total_sample <- sum(nona[[h]][,i])
    for (j in 1:length(nona[[h]][,1])){
      proportions[[h]][j,i] <- nona[[h]][j,i]/total_sample
    }
  }
}

phl=0
phyla <- list()
for (i in 1:nsets) {
  phyla[[i]] = as.data.frame(proportions[[i]] %>% group_by(Phylum) %>% summarise_if(is.numeric, funs(sum), order=F))
  phl[i]=length(phyla[[i]][,1])
} 


Phyla_distribution <- as.data.frame(matrix(nrow=nsamp*sum(phl), ncol =5))
colnames(Phyla_distribution) <- c("Phylum", "Abundance", "MOTUs", "Sample", "Dataset")

for (i in 1:nsets) {
  tmp1<-sum(!is.na(Phyla_distribution$Dataset))
  tmp2<-tmp1+1:(phl[[i]]*nsamp)
  Phyla_distribution$Dataset[tmp2]<- rep(sets[i], phl[[i]]*nsamp)
  Phyla_distribution$Phylum[tmp2] <- rep(as.character(phyla[[i]]$Phylum), nsamp) 
  for (j in 1:nsamp) {
    tmp3<-tmp1+(1+phl[i]*(j-1)):(phl[i]*j)
    Phyla_distribution$Sample[tmp3] <- as.character(factor(rep(samples[j], phl[[i]]), levels=samples))
    Phyla_distribution$Abundance[tmp3] <- as.numeric(phyla[[i]][,j+2]) ### !!! or j+3?
    tmp4<-(subset(nona[[i]], nona[[i]][,samcols[j]]!=0) %>% group_by(Phylum) %>% dplyr::count())
    Phyla_distribution$MOTUs[tmp1+phl[[i]]*(j-1)+match(tmp4$Phylum, Phyla_distribution$Phylum[tmp3])] <- tmp4$n
  }
} #...

Phyla_distribution$Dataset <- factor(Phyla_distribution$Dataset, levels = sets)
Phyla_distribution$MOTUs[is.na(Phyla_distribution$MOTUs)] <- 0

p1 <- ggplot(Phyla_distribution, aes(factor(Sample,levels=samples), Abundance, fill=Phylum)) +
  geom_bar(position = position_fill(reverse = FALSE), stat = "identity") + #the "position = ... FALSE)" seems unnecessary
  facet_wrap(~Dataset) +
  theme(legend.position = "right") +
  xlab("Sample") +
  scale_fill_manual(values = rep(c(brewer.pal(12,"Paired"),brewer.pal(8,"Set3"), brewer.pal(8,"Accent"), brewer.pal(7,"Dark2"), brewer.pal(3,"Pastel1") ))) #adjust number of colors to number of phyla # c(brewer.pal(12,"Paired"),brewer.pal(8,"Set3"), brewer.pal(8,"Accent"), brewer.pal(7,"Dark2"), brewer.pal(3,"Pastel1") )))

p2 <- ggplot(Phyla_distribution, aes(factor(Sample,levels=samples), MOTUs, fill=Phylum)) +
  geom_bar(position = position_fill(reverse = FALSE), stat = "identity") +
  facet_wrap(~Dataset) +
  xlab("Sample") +
  scale_fill_manual(values = rep(c(brewer.pal(12,"Paired"),brewer.pal(8,"Set3"), brewer.pal(8,"Accent"), brewer.pal(7,"Dark2"), brewer.pal(3,"Pastel1") )))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}
#grid.arrange(p1, p2, nrow = 1)

# if too many phyla:
mylegend<-g_legend(p1)

p3 <- grid.arrange(arrangeGrob(p1 + theme(legend.position="none"),
                               p2 + theme(legend.position="none"),
                               nrow=1),
                   nrow=2,heights=c(10, 1)) #mylegend, 
p3
grid.arrange(mylegend)

### sunburst

agg = list()
for (i in 1:nsets){
  tmp1 = paste(nona[[i]]$Kingdom,nona[[i]]$Phylum,nona[[i]]$Class,nona[[i]]$Order,nona[[i]]$Family,nona[[i]]$Genus,nona[[i]]$Species,sep=" - ")
  agg[[i]] = data.frame(matrix(0, nrow = length(tmp1), ncol = nsamp+1))
  colnames(agg[[i]])=c("Taxonomy",samples)
  agg[[i]]$Taxonomy = tmp1
  agg[[i]][,2:(nsamp+1)] = nona[[i]][,samcols]
  agg[[i]] = aggregate(agg[[i]][,-1], by = list(Taxonomy = agg[[i]][,1]), FUN = sum)
  means = vector()
  for (j in 1:nrow(agg[[i]])){
    means[j] = sum(agg[[i]][j,1:nsamp+1])/nsamp
  }
  agg[[i]][,nsamp+2] = means
  colnames(agg[[i]])[ncol(agg[[i]])]="Means"
}
i=2 # which dataset do you want to plot?
s_2_plot = "Means" # which SAMPLE to plot?
sunburst(agg[[i]][,c(1,which(colnames(agg[[i]])==s_2_plot))])

### INSECTA:
insecta = nona[[2]][which(nona[[2]]$Class=="Insecta"),]
tmp1 = paste(insecta$Class,insecta$Order,insecta$Family,insecta$Genus,insecta$Species,sep=" - ")
agg_i = data.frame(matrix(0, nrow = length(tmp1), ncol = nsamp+1))
colnames(agg_i)=c("Taxonomy",samples)
agg_i$Taxonomy = tmp1
agg_i[,2:(nsamp+1)] = nona[[2]][which(nona[[2]]$Class=="Insecta"),samcols]
agg_i = aggregate(agg_i[,-1], by = list(Taxonomy = agg_i[,1]), FUN = sum)
means = vector()
for (j in 1:nrow(agg_i)){
  means[j] = sum(agg_i[j,1:nsamp+1])/nsamp
}
agg_i[,nsamp+2] = means
colnames(agg_i)[ncol(agg_i)]="Means"
agg_i2 = agg_i[which(agg_i$Means<12),]
s_2_plot = "Means" # which SAMPLE to plot?
sunburst(agg_i2[,c(1,which(colnames(agg_i)==s_2_plot))])


agg_phy = list()
for (i in 1:nsets){
  tmp1 = paste(nona[[i]]$Kingdom,nona[[i]]$Phylum, sep=" - ")
  agg_phy[[i]] = data.frame(matrix(0, nrow = length(tmp1), ncol = nsamp+1))
  colnames(agg_phy[[i]])=c("Taxonomy",samples)
  agg_phy[[i]]$Taxonomy = tmp1
  agg_phy[[i]][,2:(nsamp+1)] = nona[[i]][,samcols]
  agg_phy[[i]] = aggregate(agg_phy[[i]][,-1], by = list(Taxonomy = agg_phy[[i]][,1]), FUN = sum)
  means = vector()
  for (j in 1:nrow(agg_phy[[i]])){
    means[j] = sum(agg_phy[[i]][j,1:nsamp+1])/nsamp
  }
  agg_phy[[i]][,nsamp+2] = means
  colnames(agg_phy[[i]])[ncol(agg_phy[[i]])]="Means"
}
i=2 # which dataset do you want to plot?
s_2_plot = "Means" # which SAMPLE to plot?
sunburst(agg_phy[[i]][,c(1,which(colnames(agg_phy[[i]])==s_2_plot))])


## ggsunburst - not used eventually - ! REMOVE NegControl COLUMN !
#
nw <- "(((a:1, b:2, c:3):2, (d:5, e:6, f:5, g:6):1):2, (f, i, h):3);"

# extract data
sb <- sunburst_data(nw)
icicle(sb)
sunburst(sb)

# demonstrate use of data with ggplot()
p <- ggplot(sb$rects)
