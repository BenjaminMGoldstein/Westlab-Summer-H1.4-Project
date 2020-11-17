library(stringr)

#insert sequence here, can repeat for different sequences
sequence="ATGTCTGAAACAGTGCCTCCCGCCCCCGCCGCTTCTGCTGCTCCTGAGAAACCTTTAGCTGGCAAGAAGGCAAAGAAACCTGCTAAGGCTGCAGCAGCCTCCAAGAAAAAACCCGCTGGCCCTTCCGTGTCAGAGCTGATCGTGCAGGCTGCTTCCTCCTCTAAGGAGCGTGGTGGTGTGTCGTTGGCAGCTCTTAAAAAGGCGCTGGCGGCCGCAGGCTACGACGTGGAGAAGAACAACAGCCGCATTAAGCTGGGCATTAAGAGCCTGGTAAGCAAGGGAACGTTGGTGCAGACAAAGGGTACCGGAGCCTCGGGTTCCTTCAAGCTCAACAAGAAGGCGTCCTCCGTGGAAACCAAGCCCGGCGCCTCAAAGGTGGCTACAAAAACTAAGGCAACGGGTGCATCTAAAAAGCTCAAAAAGGCCACGGGGGCTAGCAAAAAGAGCGTCAAGACTCCGAAAAAGGCTAAAAAGCCTGCGGCAACAAGGAAATCCTCCAAGAATCCAAAAAAACCCAAAACTGTAAAGCCCAAGAAAGTAGCTAAAAGCCCTGCTAAAGCTAAGGCTGTAAAACCCAAGGCGGCCAAGGCTAGGGTGACGAAGCCAAAGACTGCCAAACCCAAGAAAGCGGCACCCAAGAAAAAGTAA"

#clean up
seq_split=unlist(strsplit(sequence,""))
seq_split=seq_split[seq_split !="\n"]
seq_joined=paste(seq_split,collapse="")

######### INSERTIONS

#Generate dataframe with all insertions
seq_ins=as.data.frame(matrix(,nrow=length(seq_split)+1,ncol=length(seq_split)))

for (i in 1:length(seq_split)){
  seq_ins[,i]=append(seq_split,seq_split[i],after=i)
}

#If want all insertions joined

Ins_Joined=as.data.frame(matrix(,nrow=10,ncol=length(seq_split)))

for (i in 1:ncol(Ins_Joined)){
  Ins_Joined[1,i]=paste(seq_ins[,i],collapse="")
  
}

#Find what position a stop codon is at, might need to numbers in here depending on length of amino acid and bp sequences
Ins_codons=as.data.frame(matrix(,nrow=length(seq_split)/3,ncol=length(seq_split)))
for (i in 1:ncol(Ins_Joined)){
  Ins_codons[1:216,i]=sapply(seq(1,648,by=3),function(j) substr(Ins_Joined[1,i],j,j+2))
}

#This sets na's just to 1000, can change if desired
for (i in 1:ncol(Ins_Joined)){
  Ins_codons[217,i]=grep("TAG",Ins_codons[,i])[1]
  Ins_codons[218,i]=grep("TAA",Ins_codons[,i])[1]
  Ins_codons[219,i]=grep("TGA",Ins_codons[,i])[1]
  stringg=c(Ins_codons[217,i],Ins_codons[218,i],Ins_codons[219,i])
  stringg[is.na(stringg)]=1000
  Ins_codons[220,i]=min(as.numeric(stringg))
}

d=c(1:648)
plot(d,Ins_codons[220,],xlab="bp",ylab="location of first stop codon")  
title("H1E insertions")

Ins_Joined=Ins_Joined[1,]
write.csv(Ins_Joined, "/path.csv", row.names=FALSE)

#Convert all dna sequences to amino acids, for insertions

Ins_AAs=Ins_codons
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TTT",replacement="F",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TTC",replacement="F",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TTA",replacement="L",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TTG",replacement="L",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CTT",replacement="L",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CTC",replacement="L",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CTA",replacement="L",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CTG",replacement="L",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="ATT",replacement="I",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="ATC",replacement="I",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="ATA",replacement="I",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="ATG",replacement="M",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GTT",replacement="V",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GTC",replacement="V",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GTA",replacement="V",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GTG",replacement="V",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TCT",replacement="S",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TCC",replacement="S",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TCA",replacement="S",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TCG",replacement="S",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CCT",replacement="P",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CCC",replacement="P",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CCA",replacement="P",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CCG",replacement="P",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="ACT",replacement="T",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="ACC",replacement="T",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="ACA",replacement="T",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="ACG",replacement="T",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GCT",replacement="A",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GCC",replacement="A",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GCA",replacement="A",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GCG",replacement="A",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TAT",replacement="Y",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TAC",replacement="Y",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TAA",replacement="*",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TAG",replacement="*",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CAT",replacement="H",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CAC",replacement="H",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CAA",replacement="Q",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CAG",replacement="Q",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="AAT",replacement="N",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="AAC",replacement="N",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="AAA",replacement="K",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="AAG",replacement="K",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GAT",replacement="D",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GAC",replacement="D",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GAA",replacement="E",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GAG",replacement="E",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TGT",replacement="C",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TGC",replacement="C",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TGA",replacement="*",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="TGG",replacement="W",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CGT",replacement="R",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CGC",replacement="R",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CGA",replacement="R",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="CGG",replacement="R",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="AGT",replacement="S",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="AGC",replacement="S",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="AGA",replacement="R",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="AGG",replacement="R",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GGT",replacement="G",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GGC",replacement="G",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GGA",replacement="G",fixed=TRUE)
Ins_AAs=lapply(Ins_AAs,gsub,pattern="GGG",replacement="G",fixed=TRUE)

Ins_AAs=as.data.frame(matrix(unlist(Ins_AAs),nrow=length(unlist(Ins_AAs[1]))))
Ins_AAs=Ins_AAs[-(217:220),]

#Join all amino acids
AA_Joined=as.data.frame(matrix(,nrow=1,ncol=length(seq_split)))
for (i in 1:ncol(AA_Joined)){
  AA_Joined[,i]=paste(Ins_AAs[,i],collapse="")
}

#delete everything after first stop codon
for (i in 1:ncol(AA_Joined)){
  AA_Joined[,i]=gsub("[*].*$","",AA_Joined[,i])
}

######To calculate net charge for insertions
AAs_forcharge=as.data.frame(matrix(,nrow=(nrow(Ins_AAs[1])+1),ncol=length(seq_split)))
for (i in 1:ncol(AA_Joined)){
  AAs_forcharge[1:nchar(AA_Joined[,i]),i]=unlist(strsplit(AA_Joined[,i],""))
}
AAs_forcharge[is.na(AAs_forcharge)]=0
for (i in 1:ncol(AAs_forcharge)){
  Ni=c(0.997)
  Nj=c(0.99997)
  for (y in 1:nrow(AAs_forcharge[1])){
    if (AAs_forcharge[y,i]=="R"){
      Ni=append(Ni,0.99999)
    } else if (AAs_forcharge[y,i]=="K"){
      Ni=append(Ni,0.9997)
    } else if (AAs_forcharge[y,i]=="H"){
      Ni=append(Ni,0.0909)
    } else if (AAs_forcharge[y,i]=="D"){
      Nj=append(Nj,0.9992)
    } else if (AAs_forcharge[y,i]=="E"){
      Nj=append(Nj,0.9982)
    } else if (AAs_forcharge[y,i]=="C"){
      Nj=append(Nj,0.0446)
    } else if (AAs_forcharge[y,i]=="Y"){
      Nj=append(Nj,0.00085)
    }
  }
  Z=sum(Ni)-sum(Nj)
  AAs_forcharge[217,i]=Z
}

d=c(1:648)
plot(d,AAs_forcharge[217,],xlab="bp",ylab="Net charge protein")
title("insertions- Net charge")

#######END INSERTIONS


####BEGIN DELETIONS
#Generate dataframe with all deletions
seq_del=as.data.frame(matrix(,nrow=length(seq_split)-1,ncol=length(seq_split)))

for (i in 1:length(seq_split)){
  seq_del[,i]=seq_split[-i]
}

#If want all deletions joined
Del_Joined=as.data.frame(matrix(,nrow=10,ncol=length(seq_split)))

for (i in 1:ncol(Del_Joined)){
  Del_Joined[,i]=paste(seq_del[,i],collapse="")
}

#Locate stop codons for deletions
Del_codons=as.data.frame(matrix(,nrow=216,ncol=length(seq_split)))
for (i in 1:ncol(Del_codons)){
  Del_codons[1:216,i]=sapply(seq(1,648,by=3),function(j) substr(Del_Joined[1,i],j,j+2))
}

for (i in 1:ncol(Del_Joined)){
  Del_codons[217,i]=grep("TAG",Del_codons[,i])[1]
  Del_codons[218,i]=grep("TAA",Del_codons[,i])[1]
  Del_codons[219,i]=grep("TGA",Del_codons[,i])[1]
  stringg=c(Del_codons[217,i],Del_codons[218,i],Del_codons[219,i])
  stringg[is.na(stringg)]=1000
  Del_codons[220,i]=min(as.numeric(stringg))
}
d=c(1:648)
plot(d,Del_codons[220,],xlab="bp",ylab="location of first stop codon")
title("Deletions")

Del_Joined=Del_Joined[1,]
write.csv(Del_Joined, "/path.csv", row.names=FALSE)

#Convert all dna sequences to amino acids, for deletions
Del_AAs=Del_codons
Del_AAs=lapply(Del_AAs,gsub,pattern="TTT",replacement="F",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TTC",replacement="F",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TTA",replacement="L",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TTG",replacement="L",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CTT",replacement="L",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CTC",replacement="L",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CTA",replacement="L",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CTG",replacement="L",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="ATT",replacement="I",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="ATC",replacement="I",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="ATA",replacement="I",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="ATG",replacement="M",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GTT",replacement="V",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GTC",replacement="V",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GTA",replacement="V",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GTG",replacement="V",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TCT",replacement="S",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TCC",replacement="S",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TCA",replacement="S",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TCG",replacement="S",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CCT",replacement="P",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CCC",replacement="P",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CCA",replacement="P",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CCG",replacement="P",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="ACT",replacement="T",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="ACC",replacement="T",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="ACA",replacement="T",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="ACG",replacement="T",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GCT",replacement="A",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GCC",replacement="A",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GCA",replacement="A",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GCG",replacement="A",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TAT",replacement="Y",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TAC",replacement="Y",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TAA",replacement="*",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TAG",replacement="*",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CAT",replacement="H",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CAC",replacement="H",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CAA",replacement="Q",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CAG",replacement="Q",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="AAT",replacement="N",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="AAC",replacement="N",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="AAA",replacement="K",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="AAG",replacement="K",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GAT",replacement="D",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GAC",replacement="D",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GAA",replacement="E",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GAG",replacement="E",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TGT",replacement="C",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TGC",replacement="C",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TGA",replacement="*",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="TGG",replacement="W",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CGT",replacement="R",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CGC",replacement="R",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CGA",replacement="R",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="CGG",replacement="R",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="AGT",replacement="S",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="AGC",replacement="S",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="AGA",replacement="R",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="AGG",replacement="R",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GGT",replacement="G",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GGC",replacement="G",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GGA",replacement="G",fixed=TRUE)
Del_AAs=lapply(Del_AAs,gsub,pattern="GGG",replacement="G",fixed=TRUE)

Del_AAs=as.data.frame(matrix(unlist(Del_AAs),nrow=length(unlist(Del_AAs[1]))))
Del_AAs=Del_AAs[-(217:220),]

#Join all amino acids
AA_Del_Joined=as.data.frame(matrix(,nrow=1,ncol=length(seq_split)))
for (i in 1:ncol(AA_Del_Joined)){
  AA_Del_Joined[,i]=paste(Del_AAs[,i],collapse="")
}

#delete everything after first stop codon
for (i in 1:ncol(AA_Del_Joined)){
  AA_Del_Joined[,i]=gsub("[*].*$","",AA_Del_Joined[,i])
}

write.csv(AA_Joined, "/path.csv", row.names=FALSE)

#To calculate net charge for insertions
AAs_Del_forcharge=as.data.frame(matrix(,nrow=(nrow(Del_AAs[1])+1),ncol=length(seq_split)))
for (i in 1:ncol(AA_Del_Joined)){
  AAs_Del_forcharge[1:nchar(AA_Del_Joined[,i]),i]=unlist(strsplit(AA_Del_Joined[,i],""))
}
AAs_Del_forcharge[is.na(AAs_Del_forcharge)]=0
for (i in 1:ncol(AAs_Del_forcharge)){
  Ni=c(0.997)
  Nj=c(0.99997)
  for (y in 1:nrow(AAs_Del_forcharge[1])){
    if (AAs_Del_forcharge[y,i]=="R"){
      Ni=append(Ni,0.99999)
    } else if (AAs_Del_forcharge[y,i]=="K"){
      Ni=append(Ni,0.9997)
    } else if (AAs_Del_forcharge[y,i]=="H"){
      Ni=append(Ni,0.0909)
    } else if (AAs_Del_forcharge[y,i]=="D"){
      Nj=append(Nj,0.9992)
    } else if (AAs_Del_forcharge[y,i]=="E"){
      Nj=append(Nj,0.9982)
    } else if (AAs_Del_forcharge[y,i]=="C"){
      Nj=append(Nj,0.0446)
    } else if (AAs_Del_forcharge[y,i]=="Y"){
      Nj=append(Nj,0.00085)
    }
  }
  Z=sum(Ni)-sum(Nj)
  AAs_Del_forcharge[217,i]=Z
}

d=c(1:648)
plot(d,AAs_Del_forcharge[217,],xlab="bp",ylab="Net charge protein")
title("H1E deletions- Net charge")


######END DELETIONS##########

###BEGIN WT

#Compare to WT to correct for length
#Generate joined dataframe for WT and then codons
seq_WT=as.data.frame(matrix(,nrow=1,ncol=length(seq_split)))
for (i in 1:ncol(seq_WT)){
  seq_WT[1,i]=paste(seq_split,collapse="")
  
}

WT_codons=as.data.frame(matrix(,nrow=length(seq_split)/3,ncol=length(seq_split)))
for (i in 1:ncol(WT_codons)){
  WT_codons[,i]=sapply(seq(1,648,by=3),function(j) substr(seq_WT[1,i],j,j+2))
}

WT_AAs=WT_codons
WT_AAs=lapply(WT_AAs,gsub,pattern="TTT",replacement="F",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TTC",replacement="F",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TTA",replacement="L",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TTG",replacement="L",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CTT",replacement="L",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CTC",replacement="L",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CTA",replacement="L",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CTG",replacement="L",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="ATT",replacement="I",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="ATC",replacement="I",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="ATA",replacement="I",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="ATG",replacement="M",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GTT",replacement="V",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GTC",replacement="V",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GTA",replacement="V",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GTG",replacement="V",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TCT",replacement="S",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TCC",replacement="S",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TCA",replacement="S",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TCG",replacement="S",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CCT",replacement="P",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CCC",replacement="P",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CCA",replacement="P",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CCG",replacement="P",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="ACT",replacement="T",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="ACC",replacement="T",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="ACA",replacement="T",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="ACG",replacement="T",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GCT",replacement="A",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GCC",replacement="A",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GCA",replacement="A",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GCG",replacement="A",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TAT",replacement="Y",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TAC",replacement="Y",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TAA",replacement="*",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TAG",replacement="*",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CAT",replacement="H",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CAC",replacement="H",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CAA",replacement="Q",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CAG",replacement="Q",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="AAT",replacement="N",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="AAC",replacement="N",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="AAA",replacement="K",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="AAG",replacement="K",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GAT",replacement="D",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GAC",replacement="D",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GAA",replacement="E",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GAG",replacement="E",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TGT",replacement="C",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TGC",replacement="C",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TGA",replacement="*",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="TGG",replacement="W",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CGT",replacement="R",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CGC",replacement="R",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CGA",replacement="R",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="CGG",replacement="R",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="AGT",replacement="S",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="AGC",replacement="S",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="AGA",replacement="R",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="AGG",replacement="R",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GGT",replacement="G",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GGC",replacement="G",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GGA",replacement="G",fixed=TRUE)
WT_AAs=lapply(WT_AAs,gsub,pattern="GGG",replacement="G",fixed=TRUE)

WT_AAs=as.data.frame(matrix(unlist(WT_AAs),nrow=length(unlist(WT_AAs[1]))))
WT_AAs=WT_AAs[-(217:220),]

###WT length for insertions

WT=as.data.frame(matrix(,nrow=(nrow(WT_AAs[1])+1),ncol=length(seq_split)))
for (i in 1:ncol(WT)){
  WT[1:nchar(AA_Joined[,i]),i]=as.character(WT_AAs[1:nchar(AA_Joined[,i]),i])
}

WT[is.na(WT)]=0
for (i in 1:ncol(WT)){
  Ni=c(0.997)
  Nj=c(0.99997)
  for (y in 1:nrow(WT[1])){
    if (WT[y,i]=="R"){
      Ni=append(Ni,0.99999)
    } else if (WT[y,i]=="K"){
      Ni=append(Ni,0.9997)
    } else if (WT[y,i]=="H"){
      Ni=append(Ni,0.0909)
    } else if (WT[y,i]=="D"){
      Nj=append(Nj,0.9992)
    } else if (WT[y,i]=="E"){
      Nj=append(Nj,0.9982)
    } else if (WT[y,i]=="C"){
      Nj=append(Nj,0.0446)
    } else if (WT[y,i]=="Y"){
      Nj=append(Nj,0.00085)
    }
  }
  Z=sum(Ni)-sum(Nj)
  WT[217,i]=Z
}

points(d,WT[217,],col="red")

#Subtract WT from Mut
Net_ins=unlist(as.numeric(AAs_forcharge[217,]))-unlist(as.numeric(WT[217,]))
write.csv(Net_ins, "/path.csv", row.names=FALSE)
plot(d,Net_ins,ylim=c(-40,20),ylab="Net charge (Mut-WT)",xlab="BP")
title("Insertions")

#WT length for deletions

WT_Del=as.data.frame(matrix(,nrow=(nrow(WT_AAs[1])+1),ncol=length(seq_split)))
for (i in 1:ncol(WT_Del)){
  WT_Del[1:nchar(AA_Del_Joined[,i]),i]=as.character(WT_AAs[1:nchar(AA_Del_Joined[,i]),i])
}

WT_Del[is.na(WT_Del)]=0
for (i in 1:ncol(WT_Del)){
  Ni=c(0.997)
  Nj=c(0.99997)
  for (y in 1:nrow(WT_Del[1])){
    if (WT[y,i]=="R"){
      Ni=append(Ni,0.99999)
    } else if (WT_Del[y,i]=="K"){
      Ni=append(Ni,0.9997)
    } else if (WT_Del[y,i]=="H"){
      Ni=append(Ni,0.0909)
    } else if (WT_Del[y,i]=="D"){
      Nj=append(Nj,0.9992)
    } else if (WT_Del[y,i]=="E"){
      Nj=append(Nj,0.9982)
    } else if (WT_Del[y,i]=="C"){
      Nj=append(Nj,0.0446)
    } else if (WT_Del[y,i]=="Y"){
      Nj=append(Nj,0.00085)
    }
  }
  Z=sum(Ni)-sum(Nj)
  WT_Del[217,i]=Z
}

plot(d,AAs_Del_forcharge[217,],xlab="bp",ylab="Net charge protein")
points(d,WT_Del[217,],col="red")
title("deletions - net charge")

#Subtract WT from Mut
Net_dels=unlist(as.numeric(AAs_Del_forcharge[217,]))-unlist(as.numeric(WT_Del[217,]))
write.csv(Net_dels, "/path.csv", row.names=FALSE)
plot(d,Net_dels,ylim=c(-40,20),ylab="Net charge (Mut-WT)",xlab="BP")
title("Deletions")

