###### this script merges sample-replicates in OTU tables, including subtraction of reads found in negative control

### SET constant variables
setwd("C:/Users/ntb/Documents/MBar-Data/20_01_10_Tatry-Dunaj/Tatry/X_Detours/E_Cluster_otus")

tabname <-"5_OTU_table_0.01_ZERO.csv"
tab <- read.csv(tabname) #otu table withOUT taxonomy
tab = tab[-nrow(tab),c(3:(ncol(tab)-1),2,1,ncol(tab))] # if your data has NAs, also remove them

nreps <- 3 # SET number of replicates
nsamp <- (length(tab[1,])-3)/nreps # -3 because of 3 more columns : sort, OTU ID (beginning), sequence (end)
samnames<-substr(colnames(tab)[(1:nsamp)*nreps], 1,3) #prepare new colnames without the suffix of replicates

nc="BNK" # SET negative control (nc) - good if it's named the same in old and new names
ncn=which(samnames==nc)
no_nc = (1:nsamp)[-ncn]
ncreps = which(substr(colnames(tab),1,3) == nc)
no_nc_cols = (1:(nsamp*nreps))[-ncreps]


thr = 1 # SET threshold what number of ZERO-reads is ACCEPTED to say that OTU was found in the sample

# premenovanie je len formalita...:
### RENAME samples to new system: 
newnames=c("BA_2", "BNK", "VA_3", "MY_9", "MY_13", "ZL_1", "MS_1", "SL_1", "VS_5", "MY_6", "ME_4", "MS_4", "TE_1")
if (length(samnames)==length(newnames) & ncn==which(newnames==nc)) {
  samnames = newnames
  message("Pray that you made no mistake in order of newnames :D")
} else {message("Something's wrong with the names. Old names were retained. You might need to:
                a) include / exclude somehing from the newnames or samnames
                b) change this script
                c) or set newname of negative control the same as in samnames")}

### SET variables evolving during the process:
decontF <- data.frame(matrix(vector(),0,nsamp,dimnames = list(c(),samnames)), stringsAsFactors = F) #empty data.frame for output of replicate merging and 0-rep-filter
ncount = tab[0,c(ncreps[1],length(tab)-2:0)] # to ncount we will save max of OTUs' reads in nc, so that later we can find which OTUs got lost due to decontamination (later)
decontA = tab # all replicates will be saved here
decontA[,1:(ncol(tab)-3)] = 0

### "DECONTAMINATION" and MERGING of replicates with FILTERING
# ** This happens in the order of samples in the input csv. If you change order of samnames in advance, you mess up the data!!
for (i in 1:(nrow(tab))){ # for each OTU
  for (j in c(ncn ,no_nc)){ # each column, but start with nc (assuming it has same number of replicates as the rest), so that "contamination" reads can be subtracted from true samples
    if (j==ncn){
      tmp1 = max(tab[i,ncreps], na.rm = TRUE) # save max(no. reads of ith OTU) among nc replicates
      if (tmp1>0){ ncount[nrow(ncount)+1,] = c(tmp1,tab[i,ncol(tab)-2:0]) } # save OTUs of contamination  to ncount
      tmp2=(1:(nsamp*nreps))
      decontA[i,tmp2] = tab[i,tmp2] - tmp1 # subtract it from ith OTU's reads in all replicates of all samples except nc
      decontA[i, which(decontA[i,tmp2]<0 )] = 0 # replace negative numbers for 0
    } else {
      tmp2 = (j*nreps)-nreps
      z=0 # zero counter is first set to 0
      for (k in 1:nreps){
        if (decontA[i,(tmp2+k)]==0){
          z=z+1 # count zeros in ith OTU read abundances of replicates from jth sample
        }
      }
      if (z>thr) {tmp1 = 0} else {
        tmp1 = sum(decontA[i,tmp2+1:nreps], na.rm = TRUE) # sum reads from jth sample only if it passes the threshold: not too many zeros in read abundances among replicates of jth sample
      } 
    }
    decontF[i,j] = tmp1 # save the sum (or 0 or in case of nc:max) to new table - thus, replicates are merged
  }
}
decontF[,nsamp+1:3]<-tab[,ncol(tab)-2:0] ### #join with last 3 columns
decontA = decontA2
decontF = decontF2

### REMOVE OTUs that didn't "survive" in any sample
decontF<- decontF[-which(rowSums(decontF[,no_nc], na.rm = TRUE)==0),]
decontA<- decontA[-which(rowSums(decontA[,(1:(nsamp*nreps))[-ncreps]], na.rm = TRUE)==0),]
#reads_lost_by_decont = sum(colSums(tab[,no_nc_cols])) - sum(colSums(decontA[,no_nc_cols]))
#reads_lost_by_merging = sum(colSums(decontA[,no_nc_cols])) - sum(colSums(decontF[,-c(2,ncol(decontF)-0:2)]))
otus_lost_in_decont = as.character(ncount$ID[which(!(ncount$ID %in% decontF$ID))]) #OTUs present in nc, but not in final decont, because they had 0s all samples except nc

### REORDER columns as you wish BUT DON'T MESS UP :D
ord=c(3,5,4,10,13,11,6,1,9,8,12,7,2) # i put nc into last column
decontF[,1:nsamp] <- decontF[ord]
colnames(decontF)[1:nsamp] = newnames[ord]
for (i in 1:nreps){
  decontA[,(1:nsamp-1)*nreps+i] = decontA[,(ord-1)*nreps+i]
  colnames(decontA)[(1:nsamp-1)*nreps+i] = paste(newnames[ord], LETTERS[i], sep = "_")
}

### export
write.csv(decontF, paste(substr(tabname,1,nchar(tabname)-4), "_023c.csv", sep=""))
write.csv(decontA, paste(substr(tabname,1,nchar(tabname)-4), "_0123c.csv", sep=""))
if ( length(lost_in_decont) != 0 ) { write.csv(decontA, paste(substr(tabname,1,nchar(tabname)-4), "_023LIDC.csv", sep="")) }


### Produce .fasta file with sequences of OTUs which passed
fasta <- ""
for (i in 1:nrow(decontF)) {
  fasta[i*2-1]<- paste(">", as.character(decontF$ID[i]),sep = "")
  fasta[i*2]<- as.character(decontF[i, ncol(decontF)])
}
writeLines(fasta, "5_OTU_sub_0.01_23c.fasta")

#############
### OTUs Lost in merging replicates:
lost_1 = decontA[which(!(decontA[,ncol(decontA)-2] %in% decontF[,ncol(decontF)-2])),]
lost_1$sums = rowSums(lost_1[,1:(nsamp*nreps)])
lost_1<-lost_1[order(lost_1$sums, decreasing = TRUE),]
# Fasta
fasta <- ""
for (i in 1:nrow(lost_1)) {
  fasta[i*2-1]<- paste(">", as.character(lost_1$ID[i]),sep = "")
  fasta[i*2]<- as.character(lost_1$sequence[i])
}
writeLines(fasta, "OTUsLostInMergRep.fasta")

### LOSSES of READS in surviving OTUs
decontF$before = rep(0, nrow(decontF))
decontF$after = rep(0, nrow(decontF))
decontF$loss = rep(0, nrow(decontF))
decontA$before = rowSums(decontA[,1:((nsamp-1)*nreps)])
decontA$after = rep(0, nrow(decontA))
decontA$loss = rep(0, nrow(decontA))

for (i in 1:nrow(decontF)){
  decontF$before[i] = sum(decontA[which(decontA$ID == decontF$ID[i]),1:((nsamp-1)*nreps)]) ### no_nc_cols NOT CORRECT ANYMORE!!!
  
  decontF$after[i] = sum(decontF[i,1:(nsamp-1)])
  decontA$after[which(decontA$ID == decontF$ID[i])] = decontF$after[i]
  
  decontF$loss[i] = decontF$before[i]-decontF$after[i]
}
decontA$loss = decontA$before - decontA$after


losses_reads_surviving = decontA[which(decontA$ID %in% head(decontF$ID[order(decontF$loss,decreasing = TRUE)],20)),]
