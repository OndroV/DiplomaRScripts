# split a file into several smaller files
setwd("C:/Users/ntb/Documents/MBar-Data/20_01_10_Tatry-Dunaj/Tatry/X_Detours/E_Cluster_otus")

# fasta 100 sequences per file
orig <- readLines("OTUsLostForLQTax.fasta")
len <- length(orig)
len2 <- ceiling(len/200)
  
for (i in 1:len2){
  temp <- orig[(1+((i-1)*200)):(i*200)]
  temp <- temp[!is.na(temp)]
  filename <- paste("OTUs", i, "100.fasta", sep = "_")
  writeLines(temp, filename)
}

# csv HitTable 1000 rows per file (this is 100 OTUs if you decreased to 10 hits per OTU before downloading the HitTable)
orig <- readLines("0.1HitTable.csv")
len <- length(orig)
len2 <- ceiling(len/200)

for (i in 1:len2){
  temp <- orig[(1+((i-1)*200)):(i*200)]
  temp <- temp[which(temp[!is.na(temp)] != "")]
  temp2<- temp[1]
  temp5<- temp[1] ###
  for (j in 2:length(temp)){
    temp3 <- strsplit(temp[j-1], ",")
    temp4 <- strsplit(temp[j], ",")
    if (temp3[[1]][1] != temp4[[1]][1]) {
      temp2[j]<-temp[j]
    } else {
        if (temp3[[1]][3] != temp4[[1]][3]) { 
          temp2[j]<-temp[j]
        }  else { 
            temp5[j]<-temp[j] 
            }
    }
  }
  temp2 <- temp2[!is.na(temp2)]
  temp5 <- temp5[!is.na(temp5)]
  
  filename <- paste("AHT", i, "2.csv", sep = "_")
  writeLines(temp2, filename)
}

# reduce size of HitTable by removing sequential rows with same values

