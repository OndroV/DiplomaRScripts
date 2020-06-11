# this script takes (Gurdhhu/bioinf_scripts)-annotated HitTable from blast 1st hit per query and store only a few relevant columns into a new table

#load table and throw away empty columns
setwd("C:/Users/ntb/Documents/MBar-Data/20_01_10_Tatry-Dunaj/Tatry/O_U_cluster_otus")

anhits <- read.table("0.01_miss_HT_1st_annotated.csv", header = FALSE, sep = ",", col.names = paste("v",1:100), fill = TRUE)
anhits <- anhits[,which(!is.na(anhits[1,]))]
taxhits <-anhits[,c(1,3,12,14,2,17:(length(names(anhits))))]

#take first hit per OTU into new hit1 table
r = 2
hit1 <- taxhits[1,]
otu <- as.character(taxhits[1,1])
for (i in 1:nrow(taxhits)){
  if (otu != taxhits[i,1]){
    otu <- as.character(taxhits[i,1])
    hit1[r,] <-taxhits[i,]
    r = r+1
  }
}

hit1f<-hit1
#again throw away empty columns
#hit1f <- hit1[,1:5]
#r=6
#for (i in 6:ncol(hit1)){
#  if (T %in% (hit1[,i] != "")){
#    hit1f[,r] <- hit1[,i]
#    r=r+1
#  }
#}

colnames(hit1f) <- c("OTU", "Identity%", "Score","Name","BLAST_Access", paste("Tax", 1:(ncol(hit1f)-5), sep = ""))

# align blast taxonomy (some miss the "Eukaryota" level)
#WARNING!!! Few rows can lead to lack of comparability WORKS ONLY IF THE MISALIGNED TAXON IS BOTH IN Tax1 and Tax2!!!
#it overwrites last Tax column with penultimate #but it doesnt matter i for Genus and species u use the Name column
tmp1=7:ncol(hit1f)
a=vector()
for (i in 1:length(levels(hit1f$Tax2))){
  a<-append(a,which(hit1f$Tax1==levels(hit1f$Tax2)[i])) #save rows whose Tax1 matches Tax2 in another row
}
for (i in a){ #go to each misaligned line
  for (j in tmp1) {#go to each Tax column
    tmp2 = hit1f[i,j-1] #take that element
    tmp3 = levels(hit1f[,j]) #take levels at that column (summarised values in the column)
    if (!tmp2 %in% tmp3) { #if the element value is not yet among the levels, add it
      levels(hit1f[,j]) = append(tmp3, levels(tmp2)[as.numeric(tmp2)])
    }
  }
  hit1f[i,tmp1] <- hit1f[i,tmp1-1] #move elements of the columns in that row to right 
  hit1f[i,6] <- hit1f[match(hit1f[i,6],hit1f[-i,7]), 6] #fill Tax1 element with value of another row matching in Tax2 column 
}
levels(hit1f[,7]) = append(levels(hit1f[,7]),"Animalia")
hit1f[,7]<- replace(hit1f[,7], hit1f[,7]=="Metazoa", "Animalia")

#again throw away empty columns
hit1f2 <- hit1f[,1:5]
r=6
for (i in 6:ncol(hit1f)){
  if (T %in% (hit1f[,i] != "")){
    hit1f2[,r] <- hit1f[,i]
    r=r+1
  }
}
colnames(hit1f2) <- c("OTU", "Identity%", "Score","Name","BLAST_Access", paste("Tax", 1:(ncol(hit1f2)-5), sep = ""))

write.csv(hit1f2,"blast_1st_hits.csv", row.names = F)
