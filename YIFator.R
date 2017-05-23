#Change directory and filenames as necessary
inputfile<-"ytestdata.csv"
filenameforallelesnotfound<-"yallelesnotfound.csv"
filenameforseqhaplotypes<-"ytestseqhaplotypes.csv"
filenameforseqallelefrequencies<-"ytestseqfreq.csv"
filenameforlenallelefrequencies<-"ytestlenfreq.csv"
filenameforseqhapfreq<-"ytesthapseqfreq.csv"
filenameforlenhapfreq<-"ytesthaplenfreq.csv"
filenameforstats<-"yteststats.csv"

setwd("c:/docum180612/Francesc/Illumina Forense")

indata<-read.csv(inputfile, header=FALSE)
matrixdata<-as.matrix(indata)
nlines<-nrow(indata)
seqaldef<-read.csv("yseqaldefinitions_rev.csv", header=FALSE)
nl<-nrow(seqaldef)
maxal<-100

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
  
  stop("new allele found. Check the filenameforallelesnotfound")
}

#construct haplotypes: skip missings and aggregate duplicates
haplotypeseq<-array(rep(" "),dim=c(nind,nloci))
haplotypelen<-array(rep(" "),dim=c(nind,nloci))
hapseqnm<-array(rep(" "),dim=c(nind,nloci))
haplennm<-array(rep(" "),dim=c(nind,nloci))
hapsequ<-array(rep(" "),dim=c(nind,nloci))
haplenu<-array(rep(" "),dim=c(nind,nloci))
hapseqns<-array(rep(" "),dim=c(nind,nloci))
haplenns<-array(rep(" "),dim=c(nind,nloci))
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
  haplotypelen[j,k]<-as.character(matrixdata[i,3])
  
  
  
  sk<-1
  while (matrixdata[i+sk,2]==matrixdata[i,2]){
    haplotypeseq[j,k]<-paste0(haplotypeseq[j,k],",",alseq[i+sk])
    haplotypelen[j,k]<-paste0(haplotypelen[j,k],",",as.character(matrixdata[i+sk,3]))
    
    sk<-sk+1}
  i<-i+sk
  
  k<-k+1
  if (k %% nloci == 1){
    j<-j+1
    k<-k %% nloci}
}

#allele frequencies, considering duplicates as different alleles
alseqdp<-array(rep(" "),dim=c(nloci,maxal))
allendp<-array(rep(" "),dim=c(nloci,maxal))
alss<-array(rep(" "),dim=c(nloci,maxal))
alsss<-array(rep(" "),dim=c(nloci,maxal))
alls<-array(rep(" "),dim=c(nloci,maxal))
allss<-array(rep(" "),dim=c(nloci,maxal))
alseqnomis<-array(rep(" "),dim=c(nloci,maxal))
allennomis<-array(rep(" "),dim=c(nloci,maxal))
freqseqdp<-array(c(0:0),dim=c(nloci,maxal))
freqlendp<-array(c(0:0),dim=c(nloci,maxal))
kseqdp<-array(c(0:0),dim=nloci)
klendp<-array(c(0:0),dim=nloci)
samplesize<-array(c(0:0),dim=nloci)
freqseqrel<-array(c(0:0),dim=c(nloci,maxal))
freqlenrel<-array(c(0:0),dim=c(nloci,maxal))
freqhapseq<-array(c(0:0),dim=nind)
freqhaplen<-array(c(0:0),dim=nind)
freqhapseqrel<-array(c(0:0),dim=nind)
freqhaplenrel<-array(c(0:0),dim=nind)

for (k in (1:nloci)){
alseqdp[k,1]<-haplotypeseq[1,k]
allendp[k,1]<-haplotypelen[1,k]
}

for (k in (1:nloci)){
  
  ck<-1
  ckl<-1
  for (j in (1:nind)){
    f<-0
    c<-1
    while ((f==0) & (c<=ck)){
      if (haplotypeseq[j,k]==alseqdp[k,c]){
        f<-1
        
      }
      else {c<-c+1}
    }
    if (f==0){
      ck<-ck+1
      alseqdp[k,ck]<-haplotypeseq[j,k]
      
    }
    kseqdp[k]<-ck
    
    fl<-0
    cl<-1
    while ((fl==0) & (cl<=ckl)){
      if (haplotypelen[j,k]==allendp[k,cl]){
        fl<-1
        
      }
      else {cl<-cl+1}
    }
    if (fl==0){
      ckl<-ckl+1
      allendp[k,ckl]<-haplotypelen[j,k]
      
    }
    klendp[k]<-ckl
  
  }
  
}
#sort alleles
for (k in (1:nloci)){
  alss[k,]<-sort(alseqdp[k,])
  alls[k,]<-sort(allendp[k,])
}

#move sorted alleles to the top of the array
for (k in (1:nloci)){
  for (l in (1:kseqdp[k])){
  alsss[k,l]<-alss[k,maxal-kseqdp[k]+l]
  }
  for (l in (1:klendp[k])){
    allss[k,l]<-alls[k,maxal-klendp[k]+l]
  }
  }

#remove missing allele
for (k in (1:nloci)){
  if (alsss[k,1]==" "){
    for (l in (1:kseqdp[k]-1)){
      alseqnomis[k,l]<-alsss[k,l+1]
      allennomis[k,l]<-allss[k,l+1]
    }
    kseqdp[k]<-kseqdp[k]-1
    klendp[k]<-klendp[k]-1
  }
  else {
    alseqnomis[k,]<-alsss[k,]
    allennomis[k,]<-allss[k,]
  }
}

#absolute allele frequencies
for (k in (1:nloci)){
  for (j in (1:nind)){
    for (l in (1:kseqdp[k])){
    if (haplotypeseq[j,k]==alseqnomis[k,l]){freqseqdp[k,l]<-freqseqdp[k,l]+1}
  }
  for (l in (1:klendp[k])){
    if (haplotypelen[j,k]==allennomis[k,l]){freqlendp[k,l]<-freqlendp[k,l]+1}
  }
}
}

#relative allele frequencies
for (k in (1:nloci)){
  samplesize[k]<-sum(freqseqdp[k,])
}
for (k in (1:nloci)){
    freqseqrel[k,]<-freqseqdp[k,]/samplesize[k]
    freqlenrel[k,]<-freqlendp[k,]/samplesize[k]
}

#forgo haplotypes with missing values

cnm<-0
for (j in (1:nind)){
  nm<-0
  for (k in (1:nloci)){
    if (haplotypeseq[j,k]==" "){nm<-nm+1}
  }
  if(nm==0){
    cnm<-cnm+1
    hapseqnm[cnm,]<-haplotypeseq[j,]
    haplennm[cnm,]<-haplotypelen[j,]
    
  }
}

#find unique haplotypes and their frequencies
hapsequ[1,]<-hapseqnm[1,]
haplenu[1,]<-haplennm[1,]
freqhapseq[1]<-1
freqhaplen[1]<-1
  
khaps<-1
khapl<-1
for (j in (2:cnm)){
  fh<-0
  ch<-1
  while ((fh==0) & (ch<=khaps)){
    mt<-0
    for (k in (1:nloci)){
    if (hapseqnm[j,k]==hapsequ[ch,k]){mt<-mt+1}
    }
     if(mt==nloci) {
      freqhapseq[ch]<-freqhapseq[ch]+1
      fh<-1}
    else{ch<-ch+1}
  }
  if (fh==0){
    khaps<-khaps+1
    hapsequ[khaps,]<-hapseqnm[j,]
    freqhapseq[khaps]<-1
      }
  
  fhl<-0
  chl<-1
  while ((fhl==0) & (chl<=khapl)){
    mtl<-0
    for (k in (1:nloci)){
      if (haplennm[j,k]==haplenu[chl,k]){mtl<-mtl+1}
    }
    if(mtl==nloci) {
      fhl<-1
      freqhaplen[chl]<-freqhaplen[chl]+1}
    else{chl<-chl+1}
  }
  if (fhl==0){
    khapl<-khapl+1
    haplenu[khapl,]<-haplennm[j,]
    freqhaplen[khapl]<-1
    
  }
}


# relative haplotype frequencies
sshap<-sum(freqhapseq)
freqhapseqrel<-freqhapseq/sshap
freqhaplenrel<-freqhaplen/sshap

# haplotype diversity

sps2<-0
sps3<-0
spl2<-0
spl3<-0

for (m in (1:khaps)){
  sps2<-sps2+freqhapseqrel[m]^2
  sps3<-sps3+freqhapseqrel[m]^3
 }
for (m in (1:khapl)){
  spl2<-spl2+freqhaplenrel[m]^2
  spl3<-spl3+freqhaplenrel[m]^3
}
divhs<-(cnm/(cnm-1))*(1-sps2)
divhl<-(cnm/(cnm-1))*(1-spl2)
sdhs<-sqrt((2/(cnm*(cnm-1)))*(2*(cnm-2)*(sps3-sps2^2)+sps2-sps2^2))
sdhl<-sqrt((2/(cnm*(cnm-1)))*(2*(cnm-2)*(spl3-spl2^2)+spl2-spl2^2))

#write haplotype with sequence alleles
cadout="Individual;"
for (k in c(1:nloci)){
  cadout<-paste0(cadout, loci[k],";")
}
write(cadout,file=filenameforseqhaplotypes)
for (j in c(1:nind)){
  cadout=ind[j]
  for (k in c(1:nloci)){
    
      cadout=paste0(cadout,"; ",haplotypeseq[j,k])
    
  }
  write(cadout,file=filenameforseqhaplotypes, append=TRUE)  
}

#write lenghth and sequence allele frequencies

write(loci[1],file=filenameforseqallelefrequencies)
write(loci[1],file=filenameforlenallelefrequencies)
for (m in c(1:kseqdp[1])){
  cadout=paste0(alseqnomis[1,m],"; ",freqseqdp[1,m],"; ", round(freqseqrel[1,m],digits=4))
  write(cadout,file=filenameforseqallelefrequencies,append=TRUE)
  
}

for (m in c(1:klendp[1])){
  cadout=paste0(allennomis[1,m],"; ",freqlendp[1,m],"; ", round(freqlenrel[1,m],digits=4))
  write(cadout,file=filenameforlenallelefrequencies,append=TRUE)
}

for (k in c(2:nloci)){
  write(loci[k],file=filenameforseqallelefrequencies, append=TRUE)
  write(loci[k],file=filenameforlenallelefrequencies, append=TRUE)
  for (m in c(1:kseqdp[k])){
    cadout=paste0(alseqnomis[k,m],"; ",freqseqdp[k,m],"; ", round(freqseqrel[k,m],digits=4))
    write(cadout,file=filenameforseqallelefrequencies,append=TRUE)
  }
  for (m in c(1:klendp[k])){
    cadout=paste0(allennomis[k,m],"; ",freqlendp[k,m],"; ", round(freqlenrel[k,m],digits=4))
    write(cadout,file=filenameforlenallelefrequencies,append=TRUE)
  }
}

#write sequence and lenght haplotype frequencies
cadout<-"Haplotye;abs. freq.;rel.freq;"
for (k in c(1:nloci)){
  cadout<-paste0(cadout, loci[k],";")
}
write(cadout,file=filenameforseqhapfreq)
for (j in c(1:khaps)){
  cadout<-paste0("H",j,";",freqhapseq[j],";",round(freqhapseqrel[j],digits=4))
  for (k in c(1:nloci)){
    
    cadout<-paste0(cadout,"; ",hapseqnm[j,k])
    
  }
  write(cadout,file=filenameforseqhapfreq, append=TRUE)  
}
cadout<-"Haplotye;abs. freq.;rel.freq;"
for (k in c(1:nloci)){
  cadout<-paste0(cadout, loci[k],";")
}
write(cadout,file=filenameforlenhapfreq)
for (j in c(1:khapl)){
  cadout<-paste0("H",j,";",freqhaplen[j],";",round(freqhaplenrel[j],digits=4))
  for (k in c(1:nloci)){
    
    cadout=paste0(cadout,"; ",haplennm[j,k])
    
  }
  write(cadout,file=filenameforlenhapfreq, append=TRUE)  
}

#write seq and len statistics

write(paste0("Number of haplotypes (without missing alleles)=",sshap),file=filenameforstats)
write(paste0("Number of different sequence haplotypes=",khaps),file=filenameforstats,append=TRUE)
write(paste0("Sequence haplotype diversity=",round(divhs,digits=4),", sd=", round(sdhs,digits=4)),file=filenameforstats,append=TRUE)
write(paste0("Number of different length haplotypes=",khapl),file=filenameforstats,append=TRUE)
write(paste0("Length haplotype diversity=",round(divhl,digits=4),", sd=", round(sdhl,digits=4)),file=filenameforstats,append=TRUE)
