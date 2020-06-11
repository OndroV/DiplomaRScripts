# Combine OTU read abundance tables with tables of taxonomy 

setwd("C:/Users/ntb/Documents/MBar-Data/20_01_10_Tatry-Dunaj/Tatry/X_Detours/E_Cluster_otus")
tabname<-"5_OTU_table_0.01_ZERO_023c.csv"

### Read abundance table
otus <- read.csv(tabname)[,c(2:15,17,16)] #[,c(2:6,8,7)]

### BOLDigger taxonomy table
bold <- read.csv("BOLDigger_023c.csv", stringsAsFactors = F) # if u use BOLDigger, first save the hit tab as *.csv
bold <- bold[,-10]
x = gsub("," , "." , bold$Similarity)   #need to switch to numeric 
x = gsub("No Match" , NA, x)
bold$Similarity = as.numeric(x)
bold$ID=gsub(">", "", bold$ID)          #remove the ">" from ID

colnames(bold)[1]="OTU"
colnames(otus)[14]="OTU"
taxbest<-bold # zdanlivo zbytocne, ale je to kvoli kompatibilite s nasledujucou castou, ktora je z ineho skriptu

### merge tables + sort 
taxotu<-merge(taxbest,otus, by="OTU")
taxotu<-taxotu[order(taxotu$sort),]

### add Kingdom column
my_phy_lev = levels(as.factor(taxotu$Phylum)) # phyla present in current dataset
# with each new dataset check & UPDATE lists of "phyla"
Animalia = c("Annelida","Arthropoda","Bryozoa","Cercozoa","Cnidaria","Echinodermata","Gastrotricha","Chordata","Mollusca","Nematoda","Nemertea","Onychophora","Placozoa","Porifera","Rotifera","Tardigrada")
Fungi = c("Ascomycota","Basidiomycota","Myxomycota")
Bacteria = c("Actinobacteria","Firmicutes","Proteobacteria")
Plantae = c("Bryophyta","Charophyta","Chlorophyta","Magnoliophyta","Rhodophyta")
Protozoa = c("Amoebozoa")
Chromista = c("Bacillariophyta","Cryptophyta","Haptophyta","Heterokontophyta","Ochrophyta")
Archaea = c("Crenarchaeota","Euryarchaeota")
Other = c("No Match","Eukarya_unassigned")

phy = c(Animalia,Fungi,Bacteria,Plantae,Protozoa,Chromista,Archaea,Other)
k_list = list(Animalia,Fungi,Bacteria,Plantae,Protozoa,Chromista,Archaea,Other)
k_lev = c("Animalia", "Fungi", "Bacteria", "Plantae", "Protozoa", "Chromista", "Archaea", "Other")
king=vector()
for (i in 1:length(k_list)){
  king[(length(king)+1):(length(king)+length(k_list[[i]]))] = rep(k_lev[i], length(k_list[[i]]))
}
#phy_king = data.frame(phy,king) # just for easier checking

Kingdom=factor(rep("",nrow(taxotu)),levels = k_lev)
taxotu = data.frame(taxotu[,1], Kingdom, taxotu[,2:ncol(taxotu)]) #mbrave = data.frame(mbrave[,1], Kingdom, mbrave[,2:ncol(mbrave)])
colnames(taxotu)[1]="OTU"
for (i in 1:nrow(taxotu)){
  taxotu[i,2] = king[which(phy==taxotu[i,3])]  
}

### same but for mbrave table - it is already merged 
mbrtab = "mbrave_c2&3Reps.csv"
mbrave = read.csv(mbrtab, stringsAsFactors = F)[,-1]

my_phy_lev = levels(as.factor(mbrave$Phylum)) # phyla present in current dataset
Kingdom=factor(rep("",nrow(mbrave)),levels = k_lev)
mbrave = data.frame(mbrave[,1], Kingdom, mbrave[,2:ncol(mbrave)])
colnames(mbrave)[1]="OTU"
for (i in 1:nrow(mbrave)){
  mbrave[i,2] = king[which(phy==mbrave[i,3])]  
}

### subsettinghigh quality assignments and certain taxonomic groups 
taxotuHQ<-taxotu[which(taxotu$Similarity>80),]
taxotuAnimalia <-taxotu[which(!is.na(match(taxotu$Phylum, Animalia))),]
taxotuAnimaliaHQ <-taxotuAnimalia[which(taxotuAnimalia$Similarity>80),]

### write out
write.csv(taxotuAnimaliaHQ,paste(substr(tabname,1,nchar(tabname)-4),"_TaxAnimaliaHQ.csv",sep = ""), row.names = F)
write.csv(taxotuAnimalia,paste(substr(tabname,1,nchar(tabname)-4),"_TaxAnimalia.csv",sep = ""), row.names = F)
write.csv(taxotuHQ,paste(substr(tabname,1,nchar(tabname)-4),"_TaxAllHQ.csv",sep = ""), row.names = F)
write.csv(taxotu,paste(substr(tabname,1,nchar(tabname)-4),"_TaxAll.csv",sep = ""), row.names = F)
write.csv(mbrave,paste(substr(mbrtab,1,nchar(mbrtab)-4),"k.csv",sep = "_"), row.names = F) # k ako kingdom - ze sme tam pridali kingdom column

#############
### OTUs Lost for LOW tax QUALITY:
lost_LQ_Tax = taxotu[which(!(taxotu$OTU %in% taxotuHQ$OTU)),]
lost_LQ_Tax$sums = rowSums(lost_LQ_Tax[,11:22])
lost_LQ_Tax<-lost_LQ_Tax[order(lost_LQ_Tax$sums, decreasing = TRUE),]
# Fasta
fasta <- ""
for (i in 101:nrow(lost_LQ_Tax)) { # alebo 1:nrow(lost_LQ_Tax)
  fasta[i*2-1]<- paste(">", as.character(lost_LQ_Tax$OTU[i]),sep = "")
  fasta[i*2]<- as.character(lost_LQ_Tax$sequ[i])
}
writeLines(fasta, "OTUsLostForLQTaxAfter100.fasta")


################################ do the following ONLY AFTER BLASTing lost_LQ_Tax fasta, annotating it and merging to 1 table, 1 hit per OTU!!
library(stringr)
### compare taxotu with NEW blast

# load
ht1tax = read.csv("losses/OTUsLostForLQTax/ht/LostTax_LQ_1sANNOT.csv", header = FALSE)[,-c(2,4:14,17,18)]
ht1tax = ht1tax[which(ht1tax$V3 >= 80),]
# remove useless columns
lastcol = vector()
for (i in 1:nrow(ht1tax)) {
  lastcol[i] = min(c(which(is.na(ht1tax[i,]))[1]-1,which(ht1tax[i,]=="")[1]-1), na.rm = TRUE)
}
ht1tax = ht1tax[,1:max(lastcol)]
#reach species
slipped = which(str_detect(ht1tax[,4],regex(".-...-.."))) # some are moved 1 column left
ht1tax[slipped,4] = ht1tax[slipped,3]
cand = which(str_detect(ht1tax[,4],regex("candidat", ignore_case = TRUE))) # some are "Candidatus Planktophila dulcis"...
for (i in 1:length(cand)){ # remove word Candidatus
    tmp1 = strsplit(ht1tax[cand[i],4]," ")[[1]]
    ht1tax[cand[i],4] = paste(tmp1[2], tmp1[3], sep=" ")
}
ht1tax[,4] = sub(".*? ","",ht1tax[,4]) # remove first word from here (most probably genus)
for (i in 5:ncol(ht1tax)){
  ht1tax[which(ht1tax[,i]=="Metazoa"),i] = "Animalia"
  ht1tax[which(ht1tax[,i]=="Viridiplantae"),i] = "Plantae"
}
ht1tax$stop = 8
ht1tax$stop[which(ht1tax$V3 <85)] = 4 # to know how many tax levels to fill :)
ht1tax$stop[which(ht1tax$V3 >85 & ht1tax$V3<90)] = 5
ht1tax$stop[which(ht1tax$V3 >90 & ht1tax$V3<95)] = 6
ht1tax$stop[which(ht1tax$V3 >95 & ht1tax$V3<97)] = 7

in_blast = which(taxotu$OTU %in% ht1tax$V1)
# taxotu2 = taxotu # SAVE
# taxotu = taxotu2 # LOAD
taxotu$Status[in_blast] = "NCBI"
for (i in 2:8){ # each taxonomic level # wanna retry? LOAD taxotu
  taxotu[,i] = factor(taxotu[,i],levels = c(levels(as.factor(taxotu[,i])),"to_be_continued")) # allow failure of matching
  for (j in 1:nrow(ht1tax)){
    if (i <= ht1tax$stop[j]) {
      tmp2 = which(taxotu$OTU == ht1tax[j,1]) #  position of current OTU  in taxotu table
      if (i == 2) { 
        taxotu$Similarity[tmp2] = ht1tax[j,2]} # in 1st cycle replace Similarity
        taxotu[tmp2, which(taxotu[tmp2,] == "No Match")] = "" # remove "No Match" values
      if (i == 8) { 
        if (!(ht1tax[j,4] %in% levels(taxotu[,i]))) { # need to allow new species level
          taxotu[,i] = factor(taxotu[,i],levels = c(levels(as.factor(taxotu[,i])), ht1tax[j,4])) 
        }
        taxotu[tmp2,i] = ht1tax[j,4]  
      } else {tmp3 = which(ht1tax[j,5:lastcol[j]] %in% levels(taxotu[,i])) # find which column !!(-tcols[1]+1)!! has the taxonomic info that we want
        if (length(tmp3) == 0) { taxotu[tmp2,i] = "to_be_continued"  # if there's no such, write t_b_c
        } else {  taxotu[tmp2,i] = ht1tax[j,tcols[1]-1+tmp3] # else assign it
        }
      }
    }
  }
}
### subsettinghigh quality assignments and certain taxonomic groups 
taxotuHQ<-taxotu[which(taxotu$Similarity>80),]
taxotuAnimalia <-taxotu[which(!is.na(match(taxotu$Phylum, Animalia))),]
taxotuAnimaliaHQ <-taxotuAnimalia[which(taxotuAnimalia$Similarity>80),]

### write out
write.csv(taxotuAnimaliaHQ,paste(substr(tabname,1,nchar(tabname)-4),"_TaxAnimaliaHQ+blast.csv",sep = ""), row.names = F)
write.csv(taxotuAnimalia,paste(substr(tabname,1,nchar(tabname)-4),"_TaxAnimalia+blast.csv",sep = ""), row.names = F)
write.csv(taxotuHQ,paste(substr(tabname,1,nchar(tabname)-4),"_TaxAllHQ+blast.csv",sep = ""), row.names = F)
write.csv(taxotu,paste(substr(tabname,1,nchar(tabname)-4),"_TaxAll+blast.csv",sep = ""), row.names = F)


# better_blast = which(taxotu$Similarity[in_blast]<ht1$V3[ht1$V1 == taxotu$OTU[in_blast]]) # I guess all in ht1 are better, because they have >80% and they were discarded from taxotus because of LQ :D
# copy taxonomic information by finding which column in BLAST fits the particular tax level in taxotu