# similar purpose as MergeReplicates.R, but with mBrave data
###### this script merges sample-replicates in BIN tables, including subtraction of reads found in negative control

### LOAD DATA
setwd("C:/Users/ntb/Documents/MBar-Data/scripts/SkriptyOndroCheck/skript1")
tsvs <- list.files()
tsvs = tsvs[sub('.*(?=.{3}$)', '', tsvs, perl=T)=="tsv"] # only tsv files
a <- list()
for (i in 1:length(tsvs)) {
  a[[i]] <- read.csv(tsvs[i], sep = "\t", stringsAsFactors = F)
  a[[i]] <- a[[i]][a[[i]][16]>200,] # remove BINs with low overlap
}

# NASLEDUJUCI ODSEK NETREBA KONTROLOVAT
## I need to rename files from Starolesnianske pleso to get 3-letter code followed by replicate code (Stl-A instead of Stles-A)
x = which(substr(tsvs, 19,23)=="Stles") # find Stles tsvs
y = tsvs[x]
substring(y,first = 22)=substring(y,first = 24) # move filename by 2 letters to the left
tsvs[x] = substr(y,1,nchar(y)-2) # remove the extra "sv" from the filename and replace

### basic processing
##  establish variables used later
names(a)=substr(tsvs, 19,23) # name each tsv with its replicate name - names of mymbrave runs are 19th-23rd letter in the file name
rep=levels(as.factor((names(a)))) # save replicate names
names3=substr(names(a),1,3) # substring each tsv's name to sample name
samp<-levels(as.factor(names3)) # sample names
nc = "BNK" # SET NAME OF NC
ncn = which(samp==nc) # position of negative control in samp 
sampnew<-c("BA_2", "BNK", "VA_3", "MY_9", "MY_13", "ZL_1", "MS_1", "SL_1", "VS_5", "MY_6", "ME_4", "MS_4", "TE_1") # new names in alphabetic order of sample names
ord=c(3,5,4,10,13,11,6,1,9,8,12,7,2) # order of samples for final table

### JOIN TABLES
tab1=a[[1]]
for (i in 2:length(a)){
  tmp1=a[[i]][!(a[[i]][,1] %in% tab1[,1]),] # take those BINs from ith sample, which are not yet in tab2
  tab1=rbind(tab1,tmp1)
}
tab1=tab1[,c(1:8,15)] # Remove columns from Sequences to the end, except similarity
tab1[,rep]=NA # add empty columns with names of replicates
for (i in 1:length(a)){ # each replicate has its own column
  tmp2=which(colnames(tab1)==names(a)[i]) # remember column number for current replicate
  for (j in 1:nrow(a[[i]])){ # Take each BIN from the sample, 
    tab1[ match(a[[i]][j,1],tab1[,1]), tmp2 ] = a[[i]][j,9] # find which row includes the BIN and add its number of Sequences to that sample's column.
  }
}

### REMOVE CONTAMINATION, KEEP only BINs in 2/3 REPLICATES, SUM their READS
b = list()
for (i in c(ncn ,(1:length(samp))[-ncn])) { # go sample by sample, but take nc first
  tmp1=which(names3==samp[i]) # take replicates of the sample
  tmp2=data.frame()
  for (j in 1:length(tmp1)){ # put their important columns in common data.frame tmp2
    tmp4=(nrow(tmp2)+1) # save number of 1st row to fill
    tmp5=nrow(tmp2)+nrow(a[[tmp1[j]]]) # save number of last row to fill
    tmp2[tmp4:tmp5,1:10] = a[[tmp1[j]]][,c(1:8,15,9)] # fill the range with relevant columns
  }
  
  if (i==ncn) { # NEGATIVE CONTROL (nc)
    func=max # reads in nc are not summed, but maximum will be recorded per BIN
    ## (would averaging nc be better? some kind of custom mean function(my_mean){sum(x)/nreps} where nreps=length(rep)/length(samp))
    tmp5=tmp2[!duplicated(tmp2[,1]),] # for BNK I need all BINs, but dont want them to repeat
    ncount = tmp5 # save ncount to be used in later rounds of the current "for" cycle
    } else { # OTHER SAMPLES
    func=sum # will be summed
    for (j in 1:nrow(tmp2)) { # for each BIN
      tmp3=ncount$Seq[which(ncount[,1]==tmp2[j,1])] # save no. reads of current BIN from nc
      if (length(tmp3)==0) {tmp3=0} # if the BIN didnt appear in nc, set 0 (otherwise error)
      tmp2$Sequences[j] = tmp2$Sequences[j]-tmp3 # subtract no. seq in nc from BIN's no. sequences
      if (tmp2$Sequences[j]<0){tmp2$Sequences[j]=0} # if result is <0 , set it to 0 (otherwise negatives spoil the sum)
    }
    tmp2 = tmp2[which(tmp2$Sequences>0),] # subset rows with at least 1 read
    tmp4 = tmp2[which(duplicated(tmp2[,1])),] # subset rows with at least 1 read in at least 2 replicates
    tmp5 = tmp4[duplicated(tmp4[,1]),] # subset to BINs found in 3 replicates
    tmp4 = tmp4[!(tmp4[,1] %in% tmp5[,1]),] # subset to BINs found in exactly 2 replicates
    if (nrow(tmp4)>0) { #this condition is to prevent error
      tmp5[(nrow(tmp5)+1):(nrow(tmp5)+nrow(tmp4)),]=tmp4 # merge 3rep and 2rep back to one table
    }  # no BINs repeat in table tmp5 anymore
  } 
  
  for (k in 1:length(tmp5[,1])){
      tmp5$Sequences[k]=func(tmp2$Sequences[which(tmp2[,1]==tmp5[k,1])], na.rm=TRUE) # sum reads of the BINs, but from the table with all non-0 replicates still present // in case of nc it is not sum, but max
  }
  b[[i]]=tmp5 # save the current table of read numbers of the current sample for later
}

### remove contamination also from tab1
for (i in 1:nrow(tab1)){
  tmp1 = which(ncount[,1] == tab1[i,1]) # save rownumber in table nc where the current BIN is
  tmp2 = 10:ncol(tab1) # read abundance columns
  if (length(tmp1)==1) { # but hey, does it even appear in nc? (this together with previous line is to prevent errors)
    tab1[i,tmp2] = tab1[i,tmp2]-ncount$Sequences[tmp1] # if it appears in nc, then subtract nc reads also from tab1
    tab1[i,which(tab1[i,]<0)] = 0
  }
}


### JOIN the tmp5 TABLES
tab2=b[[1]]
for (i in 2:length(b)){
  tmp1=b[[i]][!(b[[i]][,1] %in% tab2[,1]),] # add those BINs from ith sample, which are not yet in tab2
  tab2=rbind(tab2,tmp1)
}
tab2=tab2[,-ncol(tab2)] # Remove Sequences column.
for (i in 1:length(b)){ # For each sample there will be special column.
  tmp2=ncol(tab2)+1
  for (j in 1:nrow(b[[i]])){ # Take each BIN from the sample, 
  tab2[ match(b[[i]][j,1],tab2[,1]), tmp2 ] = b[[i]][j,10] # find which row includes the BIN and add its number of Sequences to that sample's column.
  }
}

### renaming & reordering of tab2
tmp1=tab2 # save the current version of tab2 to prevent mistakes
for (i in 1:length(samp)){ # each sample
  tab2[,i+9]=tmp1[ord[i]+9] #reorder
  colnames(tab2)[i+9]=sampnew[ord[i]] #rename
}

# removing 0-read-BINs
d=list(tab1,tab2)
for (i in 1:2) {
  if (i==1) { 
    tmp1 = which(substr(colnames(tab1),1,3) == samp[ncn]) # which are nc reps? they will not be counted
  } else {tmp1 = which(ord==ncn)+9} # which is the new nc?
  tmp2 = c(10:ncol(d[[i]])) # non-nc columns
  tmp2 = tmp2[-(which(tmp2==tmp1))] # non-nc columns...
  tmp3 = rowSums(d[[i]][,tmp2], na.rm = TRUE) # sum non-nc columns
  d[[length(d)+1]] = d[[i]][which(tmp3==0),] # save BINs that finally remain only in nc
  d[[i]] = d[[i]][-which(tmp3==0),] # remove them from final table
}

### EXPORT tables
write.csv(d[[1]],"mbrave_cAllReps.csv")
write.csv(d[[2]],"mbrave_c2&3Reps.csv")
write.csv(d[[3]],"mbrave_AllLost_in_decontamination.csv")
write.csv(d[[4]],"mbrave_2&3Lost_in_decontamination.csv")