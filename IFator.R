#Change directory and filenames as necessary
inputfile="testdata.csv"
filenameforallelesnotfound="allelesnotfound.csv"
filenameforseqgenotypes="testseqgenotypes.csv"
filenameforseqallelefrequencies="testseqfreq.csv"
filenameforlenallelefrequencies="testlenfreq.csv"
filenameforstats="teststats.csv"

setwd("c:/docum180612/Francesc/Illumina Forense")

maxal<-100

indata<-read.csv(inputfile, header=FALSE)
matrixdata<-as.matrix(indata)
nlines<-nrow(indata)
seqaldef<-read.csv("seqaldefinitions_rev.csv", header=FALSE)
nl<-nrow(seqaldef)

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
      if(i<nlines-1){
      if((as.matrix(indata[i+2,1])==as.matrix(indata[i,1])) & (as.matrix(indata[i+2,2])==as.matrix(indata[i,2]))){
      cadout=paste0("Three alleles found in individual ",indata[i,1]," at locus ",indata[i,2])
      stop(cadout)
      } 
      }
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
  
  stop("new allele found. Check the filenameforallelesnotfound")
}

genotype<-array(rep(" "),dim=c(nind,nloci,2))
nchr<-array(c(0:0), dim=c(nind,nloci))
freqallen<-array(c(0:0), dim=c(maxal,nloci))
freqalseq<-array(c(0:0), dim=c(maxal,nloci))
freqallenrel<-array(c(0:0), dim=c(maxal,nloci))
freqalseqrel<-array(c(0:0), dim=c(maxal,nloci))
samplesize<-array(c(0:0),dim=nloci)
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
  genotype[j,k,1]<-alseq[i]
  nchr[j,k]<-2
  
  #adds to seq and len allele frequencies
  m<-1
  while (alseq[i]!=refalseq[m,k]){
    m<-m+1
  }
  freqalseq[m,k]<-freqalseq[m,k]+1
  
  mm<-1
  while (as.numeric(matrixdata[i,3])!=allen[mm,k]){
    mm<-mm+1
  }
  freqallen[mm,k]<-freqallen[mm,k]+1
  
  if (matrixdata[i+1,2]!=matrixdata[i,2]){
    genotype[j,k,2]<-" "
    freqalseq[m,k]<-freqalseq[m,k]+1
    freqallen[mm,k]<-freqallen[mm,k]+1
  i<-i+1} 
  else {
    genotype[j,k,2]<-alseq[i+1]
    
    #adds to seq and len allele frequencies
    m<-1
    while (alseq[i+1]!=refalseq[m,k]){
      m<-m+1
    }
    freqalseq[m,k]<-freqalseq[m,k]+1
    
    mm<-1
    while (as.numeric(matrixdata[i+1,3])!=allen[mm,k]){
      mm<-mm+1
    }
    freqallen[mm,k]<-freqallen[mm,k]+1
    
  i<-i+2
  }
  k<-k+1
  if (k %% nloci == 1){
    j<-j+1
   k<-k %% nloci}
}

#write sequence genotypes
cadout=" ;"
for (k in c(1:nloci)){
  cadout<-paste0(cadout, loci[k],"; ;")
}
write(cadout,file=filenameforseqgenotypes)
for (j in c(1:nind)){
  cadout=ind[j]
  for (k in c(1:nloci)){
    for (l in c(1:2)){
      cadout=paste0(cadout,"; ",genotype[j,k,l])
    }
  }
  write(cadout,file=filenameforseqgenotypes, append=TRUE)  
}


for (k in c(1:nloci)){
  samplesize[k]<-sum(nchr[,k])
  freqalseqrel[,k]<-freqalseq[,k]/samplesize[k]
  freqallenrel[,k]<-freqallen[,k]/samplesize[k]
}

#writes allele frequencies
write(loci[1],file=filenameforseqallelefrequencies)
write(loci[1],file=filenameforlenallelefrequencies)
for (m in c(1:kseq[1])){
  cadout=paste0(refalseq[m,1],"; ",freqalseq[m,1],"; ", round(freqalseqrel[m,1],digits=4))
  write(cadout,file=filenameforseqallelefrequencies,append=TRUE)
 
}

for (m in c(1:klen[1])){
  cadout=paste0(allen[m,1],"; ",freqallen[m,1],"; ", round(freqallenrel[m,1],digits=4))
  write(cadout,file=filenameforlenallelefrequencies,append=TRUE)
}

  for (k in c(2:nloci)){
  write(loci[k],file=filenameforseqallelefrequencies, append=TRUE)
  write(loci[k],file=filenameforlenallelefrequencies, append=TRUE)
  for (m in c(1:kseq[k])){
    cadout=paste0(refalseq[m,k],"; ",freqalseq[m,k],"; ", round(freqalseqrel[m,k],digits=4))
    write(cadout,file=filenameforseqallelefrequencies,append=TRUE)
  }
  for (m in c(1:klen[k])){
    cadout=paste0(allen[m,k],"; ",freqallen[m,k],"; ", round(freqallenrel[m,k],digits=4))
    write(cadout,file=filenameforlenallelefrequencies,append=TRUE)
    }
  }

#computing informativity statistics
acklen<-array(c(0:0),dim=nloci)
ackseq<-array(c(0:0),dim=nloci)
hetlen<-array(c(0:0),dim=nloci)
hetseq<-array(c(0:0),dim=nloci)
podlen<-array(c(0:0),dim=nloci)
podseq<-array(c(0:0),dim=nloci)
celen<-array(c(0:0),dim=nloci)
ceseq<-array(c(0:0),dim=nloci)

for (k in c(1:nloci)){
  sumpod<-0.0
  sumce<-0
  sumce2<-0
  ssumpod<-0
  ssumce<-0
  ssumce2<-0
  sp<-0
  ssp<-0
for(m in c(1:klen[k])){
  if (freqallen[m,k]>0) {
    acklen[k]<-acklen[k]+1
    hetlen[k]<-hetlen[k]+freqallenrel[m,k]^2
    p=freqallenrel[m,k]
    sumce<-sumce+p*(1-p+p^2)*(1-p)^2
    sp<-sp+p^4
  }
}

for (m1 in 1:klen[k]){
  
  for (m2 in 1:klen[k]){
    pi=freqallenrel[m1,k]
    pj=freqallenrel[m2,k]
    if (m1!=m2){
    sumpod<-sumpod+4*pi^2*pj^2
    sumce2<-sumce2+pi*pj*(pi+pj)*(1-pi-pj)^2
    }
    }
}
sumpod<-sumpod/2
sumce2<-sumce2/2

for(m in c(1:kseq[k])){
  if (freqalseq[m,k]>0) {
    ackseq[k]<-ackseq[k]+1
    hetseq[k]<-hetseq[k]+freqalseqrel[m,k]^2
    ssp<-ssp+freqalseqrel[m,k]^4
    ssumce<-ssumce+freqalseqrel[m,k]*(1-freqalseqrel[m,k]+freqalseqrel[m,k]^2)*(1-freqalseqrel[m,k])^2}
}

for (m1 in c(1:kseq[k])){
  for (m2 in c(1: kseq[k])){
    if (m1!=m2){
    ssumpod<-ssumpod+4*freqalseqrel[m1,k]^2*freqalseqrel[m2,k]^2
    ssumce2<-ssumce2+freqalseqrel[m1,k]*freqalseqrel[m2,k]*(freqalseqrel[m1,k]+freqalseqrel[m2,k])*(1-freqalseqrel[m1,k]-freqalseqrel[m2,k])^2
  }
  }
}
ssumpod<-ssumpod/2
ssumce2<-ssumce2/2

podlen[k]=1-(sp+sumpod)
celen[k]=sumce+sumce2
podseq[k]=1-(ssp+ssumpod)
ceseq[k]=ssumce+ssumce2

hetlen[k]<-1-hetlen[k]
hetseq[k]<-1-hetseq[k]
}

cadout="Locus; 2n;K (len); K(seq); het(len); het(seq); POD(len); POD(seq); CE(len); CE(seq)"
write(cadout,file=filenameforstats)
for (k in 1:nloci){
  cadout=paste0(loci[k],"; ", samplesize[k],"; ",acklen[k],"; ",ackseq[k],"; ",round(hetlen[k],digits=4),"; ",round(hetseq[k],digits=4),"; ",round(podlen[k],digits=4),"; ",round(podseq[k],digits=4),"; ",round(celen[k],digits=4),"; ",round(ceseq[k],digits=4))
  write(cadout,file=filenameforstats,append=TRUE)
}

