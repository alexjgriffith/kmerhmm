library(mulcal)
     library(kmerhmm)
     fastaFile<-"~/Dropbox/UTX-Alex/jan/combined.fasta"
     fastaFile<-"~/Dropbox/UTX-Alex/jan/combined.fasta"
     heightFile<-"~/Dropbox/UTX-Alex/jan/combined_heights.bed"
     temp<-read.table(heightFile)
     data<-as.matrix(temp[,4:dim(temp)[2]])
     dataLocs<-list(list("eryt", 1, "top"),
                   list("tall", 1,"bottom"),
                   list("ecfc",3,"top"),
                   list("other",c(3,5),c("top","top")),
                   list("meka",c(3,7),c("top","top")),
                   list("hspc",c(3,7),c("top","bottom")),
                   list("diff", 7, "top"))
     normData<-qn(data)
     prc<-prcomp(t(normData))$rotation
     mers<-countNMers(fastaFile,prc[,1],8)

wei2<-mers[,2]
names(wei2)<-as.character(mers[,1])
medWei<-median(wei2)

MAD<-median(abs(wei2-medWei))

sigma=4*MAD/0.6745

erythroid<-sort(wei2[wei2>sigma])

leukemic<-sort(wei2[wei2<(-sigma)])

library(msa)
library(abind)

msaResults<-msa(names(leukemic),"Muscle",cluster="neighborjoining",gapOpening = 420, gapExtension =10,type="dna")
consensusString(msaResults)

levels<-c("A","C","G","T","-")
nodes<-30
iterations<-5
#starts<-10
sequences<-as.matrix(msa(names(leukemic),type="dna"))

hmm<-buildHMM(sequences,levels,nodes)


C2DD<-function(n,dim,varname){
    paste("const double",varname, "[][]= {{",paste(apply(matrix(runif(n*dim),ncol=dim),2,function(x) paste(x/sum(x),collapse=", " )),collapse ="},{"),"}};\n")
}
