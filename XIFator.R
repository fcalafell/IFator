#Change directory and filenames as necessary
inputfile<-"xtestdata.csv"
filenameforallelesnotfound<-"xallelesnotfound.csv"
filenameforseqgenotypes<-"xtestseqgenotypes.csv"


setwd("c:/docum180612/Francesc/Illumina Forense")

indata<-read.csv(inputfile, header=FALSE)
matrixdata<-as.matrix(indata)
nlines<-nrow(indata)
seqaldef<-read.csv("xseqaldefinitions_rev.csv", header=FALSE)
nl<-nrow(seqaldef)
maxal<-200

#find individuals in sample
firstcolumn<-array(indata$V1,dim=nlines)
ind<-unique(firstcolumn)
nind<-dim(ind)

#find different sequence alleles in external definition
fc<-array(seqaldef$V1,dim=nl)
loci<-unique(fc)
nloci<-dim(loci)
refalseq<-array(rep(""),dim=c(maxal,nloci))
refalseq[1,1]<-as.character(seqaldef$V3[1])
kseq<-array(c(0:0), dim=nloci)
ct<-1
cl<-1
for (i in c(2:nl)){
  if (seqaldef$V1[i]==seqaldef$V1[i-1]){
    ct<-ct+1
    refalseq[ct,cl]<-as.character(seqaldef$V3[i])
  }
  else{
    kseq[cl]<-ct
    ct<-1
    cl<-cl+1
    refalseq[ct,cl]<-as.character(seqaldef$V3[i])
  }
}
kseq[cl]<-ct

#find different length alleles in external definition
allen<-array(c(0:0),dim=c(maxal,nloci))
klen<-array(c(0:0), dim=nloci)
ct0<-1
ct1<-0
for (i in c(1:nloci)){
  ct0<-ct1+1
  ct1<-ct0+kseq[i]-1
  tall<-array(c(seqaldef$V2[ct0:ct1]),dim=ct1-ct0+1)
  alll<-unique(tall)
  klen[i]<-dim(alll)
  
  for (j in c(1:klen[i])){
    allen[j,i]<-alll[j]
  }
}

#sequence allele calling
notfound<-0
nf<-array(c(0:0),dim=c(nlines,4))
alseq<-array(c(0:0),dim=nlines)

for (i in c(1:nlines)){
  found<-0
  ct<-0
  while ((found==0) & (ct<nl)){
    ct<-ct+1
    if((as.matrix(indata[i,2])==as.matrix(seqaldef[ct,1])) & (as.matrix(indata[i,4])==as.matrix(seqaldef[ct,4]))){
      alseq[i]<-as.character(seqaldef$V3[ct])
      found<-1                       
    }
  }
  if (found==0){notfound<-notfound+1
                nf[notfound,]<-as.matrix(indata[i,])}
}
if (notfound>0){
  cadout="sample;locus;length;sequence"
  write(cadout,file=filenameforallelesnotfound)
  for (z in (1:notfound)){
    cadout=""
    for (zz in (1:4)){
      cadout=paste0(cadout,nf[z,zz],";")
    }
    write(cadout,file=filenameforallelesnotfound,append=TRUE)  
  }

  stop(paste0("new allele found. Check ",filenameforallelesnotfound))
}

#construct haplotypes: skip missings and aggregate duplicates
haplotypeseq<-array(rep(" "),dim=c(nind,nloci))

i<-1
j<-1
k<-1
while (i<nlines){
  #skip missing values
  while (matrixdata[i,2]!=loci[k]){
    k<-k+1
    if (k %% nloci == 1){
      j<-j+1
      k<-k %% nloci}
  }
  haplotypeseq[j,k]<-alseq[i]
  
  
  
  
  sk<-1
  while (matrixdata[i+sk,2]==matrixdata[i,2]){
    haplotypeseq[j,k]<-paste0(haplotypeseq[j,k],",",alseq[i+sk])
   
    
    sk<-sk+1}
  i<-i+sk
  
  k<-k+1
  if (k %% nloci == 1){
    j<-j+1
    k<-k %% nloci}
}

#write haplotype with sequence alleles
cadout="Individual;"
for (k in c(1:nloci)){
  cadout<-paste0(cadout, loci[k],";")
}
write(cadout,file=filenameforseqgenotypes)
for (j in c(1:nind)){
  cadout=ind[j]
  for (k in c(1:nloci)){
    
    cadout=paste0(cadout,"; ",haplotypeseq[j,k])
    
  }
  write(cadout,file=filenameforseqgenotypes, append=TRUE)  
}
