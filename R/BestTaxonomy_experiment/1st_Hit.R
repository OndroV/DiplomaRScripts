# take only 1st row where an element occurs (e.g. OTU_1) and save new table



### many 1st_hit csvs?
setwd("OTUsLostForLQTax/ht")
csvs = list.files()
csvs = csvs[order(as.numeric(substr(csvs,1,nchar(csvs)-4)))]
merged = list()
ht = data.frame()
for (i in 1:length(csvs)){
  merged[[i]] = read.csv(csvs[i],header = FALSE)
  tmp1 = nrow(merged[[i]])
  if (i == 1) {ht = merged[[i]]} else {ht[(nrow(ht)+1):(nrow(ht)+tmp1),] = merged[[i]]}
}


#setwd("losses")#("C:/Users/ntb/Documents/MBar-Data/20_01_10_Tatry-Dunaj/Tatry/X_Test/cluster_otus0.01")
#tabname<-"23c_myRef.csv"
#ht <- read.csv(tabname, header = FALSE)

ht = ht[which(ht$V12>300),]
ht$sort = 1:nrow(ht)

lev<- levels(as.factor(as.character(ht[,1]))) # changes order!! ned to use sort
ht1 <- ht[1,]
ht2<-ht[1,]
n=1
for (i in 1:length(lev)){
  ht1[i,]<-ht[match(lev[i], ht[,1]),]
  if (nchar(as.character(ht[match(lev[i], ht[,1]),2]))!=10) {
    ht2[n,]<-ht[match(lev[i], ht[,1]),]
    n=n+1
  }
}
ht1=ht1[order(ht1$sort),] #reorder
ht2=ht2[order(ht2$sort),] #reorder

write.table(ht1, paste(substr(tabname,1,nchar(tabname)-4),"_1st.csv",sep = ""), quote=FALSE, sep = ",", qmethod = "double", row.names = FALSE, col.names = FALSE)
write.table(ht2, paste(substr(tabname,1,nchar(tabname)-4),"_1st_suspect.csv",sep = ""), quote=FALSE, sep = ",", qmethod = "double", row.names = FALSE, col.names = FALSE)

#wanna save by x rows?
for (i in 1:ceiling(nrow(ht1)/10)){
  write.table(ht1[((i-1)*10+1):(i*10),], paste("LostTax_LQ_1st_",i,".csv",sep = ""), quote=FALSE, sep = ",", qmethod = "double", row.names = FALSE, col.names = FALSE)
}
# after annotate.py load them back
annot = list.files(pattern = "*annotated.csv")
num_start = rep(nchar("LostTax_LQ_1st_x"),length(annot))
num_end = nchar(annot)-14
annot = annot[order(as.numeric(substr(annot,num_start,num_end)))]
merged = list()
ht1tax = data.frame(matrix(nrow = nrow(ht1), ncol = 60))
tmp2 = 0
for (i in 1:length(annot)){
  merged[[i]] = read.csv(annot[i],header = FALSE, col.names = 1:60)
  tmp1 = nrow(merged[[i]])
  ht1tax[(tmp2+1):(tmp2+tmp1),] = merged[[i]]
  tmp2 = length(which(!is.na(ht1tax$X1)))
}
write.table(ht1tax, "LostTax_LQ_1sANNOT.csv", quote=FALSE, sep = ",", qmethod = "double", row.names = FALSE, col.names = FALSE)

# which myref were top / good?
myref_top_otus=ht1$V1[which(ht1$V3>93)]
myref_good_otus=ht1$V1[which(ht1$V3<93)]
# were they better than in BOLDigger ?
taxotuTOPmyRef = taxotu[which(taxotu$OTU %in% myref_top_otus),] # need to run TaxonomyTab_Combine1.R first
taxotuGOODmyRef = taxotu[which(taxotu$OTU %in% myref_good_otus),]
