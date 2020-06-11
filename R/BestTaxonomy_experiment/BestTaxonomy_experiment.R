# 1. Combine results of blasting against my own reference tables, against GBIF, NCBI (blast) and BOLD databases.
# 2. Then choose the taxonomic assignment with highest quality (slovami taxonomic assignment myslim alignment, len sa mi to nechcelo prepisovat)


setwd("C:/Users/ntb/Documents/MBar-Data/20_01_10_Tatry-Dunaj/Tatry/X_Test/cluster_otus0.01/skript3")
tabname<-"5_OTU_table_0.01_ZERO_PCRmerged023.csv"
otus <- read.csv(tabname)[,c(2:6,8,7)] # [,c(2:15,17,16)] #

myref<-read.csv("0.01_miss_MyRef_HT_1st.csv", header = FALSE, stringsAsFactors = FALSE)[,c(1:3,12)] # to get 1 row per otu use 1st_Hit.R
colnames(myref)<-c("OTU","BARCODE","Similarity","Score")

gbif <- read.csv("0.01_miss_gbif_split.csv", stringsAsFactors = F)[,c(1,3,4,7:14)] # poznamka pre mna: separate "Ta_Xo_No_No my" first: e.g. in excel select 6 columns -> insert new columns; select taxonomy -> Data -> text to columns ... (sep="_") and remove genus names from species column (ctrl+H -> search"* ", replace"",maybe need to tick a checkbox allowing "zastupne znaky", in english sth like "substitute symbols")
colnames(gbif) <- c(colnames(gbif)[1:4],"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
for (i in 2:3){gbif[,i] <- replace(gbif[,i],which(is.na(gbif[,i])),0)} #replace NAs 

blast <- read.csv("blast_1st_hits.csv", stringsAsFactors = F) # use output of Blast_1hit.R even if u used 1st_Hit.R to ensure correct taxonomy
#to obtain blast_1st_hits.csv: blast your OTU file, download HitTable, use annotate_blast_hits.py from https://github.com/Gurdhhu/bioinf_scripts (maybe need to subsample the table, otherwise the script fails), and use Blast_1hit.R  
blast[which("environmental" == sapply(strsplit(blast$Name, split = " "), "[[", 2)),2:3] <- 0 #environmental samples are not useful for us :)

bold <- read.csv("BOLDwh_taxonomy.csv", stringsAsFactors = F)

### merge the taxonomy tables
colnames(bold)[1]="OTU"
colnames(gbif)[1]="OTU"
colnames(otus)[5]="OTU" # or change according to your table :) e.g. colnames(otus)[14]="OTU" # 
all1 <- merge(bold, gbif, by="OTU")
all2<-merge(all1,blast, by="OTU")

### MyRef merging
tmp3<- ncol(all2)
for (i in 1:nrow(myref)){
  tmp<- strsplit(myref$BARCODE[i], split = "_")
  tmp2<- c(myref[i,3:4], tmp[[1]][2])
  tmp2[4:10]<-c(rep(NA,5),tmp[[1]][4:5])
  if (!is.na(match(tmp[[1]][4],gbif$Genus))){
    tmp2[4:8]<- c(gbif[match(tmp[[1]][4],gbif$Genus), c(5:9)])
  } else { if(!is.na(match(tmp[[1]][3],gbif$Order))) {
      tmp2[4:7]<- c(gbif[match(tmp[[1]][3],gbif$Order), c(5:8)])
    } else { if(tmp[[1]][3]!="Haplotaxida"){tmp2[4:7]<-c("Animalia", "Arthropoda", "Insecta",tmp[[1]][3])
      } else {tmp2[4:7]<-c("Animalia", "Annelida", "Clitellata", "Haplotaxida")} ## only for MyRef from march 2020
    }
  }
  all2[which(all2$OTU==myref[i,1]), (1:10)+tmp3]<- tmp2
}
#check if some OTUs lost by merging with BLAST are in MyRef! If yes, there will be NAs in: for (i in 1:length(myref[,1])){tmp3<-append(tmp3, match(myref$OTU[i],all2$OTU))}
  
taxbest<-data.frame(matrix(all2[,1]))

#choose "the most reliable taxonomy" - dalo by sa vyhrat s poziadavkami na "most reliable" ale na teraz mam taketo:
#priorita je BOLD, potom zarovnanie s najvyssm Score. BOLD totiz neuvadza Score.
for (i in 1:length(all2[,1])){
  if (!is.na(all2$Phylum.x[i])){ # BOLD has priority, but it has many NA rows, because only species-level data is used
    tmp <- match(all2$Phylum.x[i],all2$Tax4) # filtering NA rows
    if(!is.na(tmp)) {
      taxbest[i,2]<-all2$Tax2[match(all2$Phylum.x[i],all2$Tax4)] #assign blast kingdom whose phylum matches
    } else {taxbest[i,2]<-"not_assigned"}
    taxbest[i,3:12]<-c(all2[i,2:8],777,all2$Status[i], "bold") #assign bold taxonomy + made-up score "777" - so that it has top priority
  } else {
    if (all2$bitScore[i]<all2$Score[i]){
      tmp<-strsplit(all2$Name[i], split = " ")
      tmp2<-c(tmp[[1]][1], paste(tmp[[1]][2:length(tmp[[1]])], collapse = "_"))
      taxbest[i,2:12]<-c(all2[i,c(25,27,29,33,36)],tmp2,all2[i,c(20,21,23)],"blast")
      #assign blast taxonomy + %ID + accession number + score
      #the selected columns might not be ideal, but for arthropoda it worked fine
    } else {
      taxbest[i,2:12]<-c(all2[i,c(13:19,10:12)],"gbif") #assign gbif taxonomy + %ID + accession number + score
    }
  if ((!is.na(all2$Phylum.x[i]) & all2$Score.1[i]>=taxbest[i,10])) {
    taxbest[i,2:12]<-c(all2[i,c(44:51, 41:43)],"MyRef")
    }
  } 
}
colnames(taxbest)<-c("OTU", "Kingdom", "Phylum", "Class", "Order","Hopefully_Family","Genus","Species","Similarity","SCORE","db_Number","which_db")


### add abundance in samples + sequences + sort number
taxotu<-merge(taxbest,otus, by="OTU")
taxotu<-taxotu[order(taxotu$sort),]

### remove low quality taxonomy assignments
taxotuHQ<-taxotu[which(taxotu$SCORE>400 & taxotu$Similarity>80),]
taxotuAnimalia<-taxotu[which(taxotu$Kingdom=="Animalia"),]
taxotuAnimaliaHQ<-taxotuAnimalia[which(taxotuAnimalia$SCORE>400 & taxotuAnimalia$Similarity>80),]

### export
write.csv(taxotuAnimaliaHQ,paste(substr(tabname,1,nchar(tabname)-4),"_TaxAnimaliaHQ.csv",sep = ""), row.names = F)
write.csv(taxotuAnimalia,paste(substr(tabname,1,nchar(tabname)-4),"_TaxAnimalia.csv",sep = ""), row.names = F)
write.csv(taxotuHQ,paste(substr(tabname,1,nchar(tabname)-4),"_TaxAllHQ.csv",sep = ""), row.names = F)
write.csv(taxotu,paste(substr(tabname,1,nchar(tabname)-4),"_TaxAll.csv",sep = ""), row.names = F)