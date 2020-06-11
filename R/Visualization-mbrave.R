# visualization of mBrave data ... then continue with Visualization.R script

### Load the libraries
### Packages not installed? -> "Error in library(pepe) : there is no package called" ...dplyr...
library(vegan)
library(ggplot2)
library(dplyr)
library(wesanderson)
library(sunburstR)
library(circlize)
library(gplots)
library(reshape)
library(reshape2)
library(RColorBrewer)
library(gridExtra)
library(grid)

sets = c("mBrave All", "mBrave Animalia")
nsets<-length(sets)

setwd("C:/Users/ntb/Documents/MBar-Data/20_01_10_Tatry-Dunaj/Tatry/Y_MBRAVE/tsv")
s1<-(read.csv("mbrave_c2&3Reps_+k.csv", stringsAsFactors = FALSE))
s2 = s1[which(s1$Kingdom=="Animalia"),]

list1<- list(s1,s2)
nona<-list1
for (i in 1:nsets) {
  nona[[i]][is.na(list1[[i]])]<-0
}

### What samples you have? And which columns are they in?
samcols <- 11:23# adjust to fit columns of samples read abundance
samples<-colnames(nona[[1]][,samcols]) 
nsamp<-length(samples)
