#' Count N Mers
#' Find the weighted value of each nmer
#' @useDynLib kmerhmm , .registration =TRUE
#' @param fastaFile The fasta file to be analyzed, each region must be represented in a single line
#' @param weights a vector of weight the same length as the number of members in the fasta file
#' @param n the number of ks
#' @examples
#' library(mulcal)
#' library(kmerhmm)
#' fastaFile<-"~/Dropbox/UTX-Alex/jan/combined.fasta"
#' heightFile<-"~/Dropbox/UTX-Alex/jan/combined_heights.bed"
#' temp<-read.table(heightFile)
#' data<-as.matrix(temp[,4:dim(temp)[2]])
#' dataLocs<-list(list("eryt", 1, "top"),
#'               list("tall", 1,"bottom"),
#'               list("ecfc",3,"top"),
#'               list("other",c(3,5),c("top","top")),
#'               list("meka",c(3,7),c("top","top")),
#'               list("hspc",c(3,7),c("top","bottom")),
#'               list("diff", 7, "top"))
#' normData<-qn(data)
#' prc<-prcomp(t(normData))$rotation
#' projection<-t(normData)%*%prc
#' weights<-prc[,1]
#' countNMers(fastaFile,weights,4)
#' @export
countNMers<-function(fastaFile,weights,n=6){    
    if((n-floor(n/2)*2)==0){
        #even
        c=4^(n/2)
        b=(4^n-c)/2
        length=c+b
    }
    else{
        #odd
        c=0
        b=(4^n+c)/2
        length=b
    }
    file<-normalizePath(fastaFile)
    results<-.C("countNMers",                
                file,as.integer(n),weights,
                numeric(length),character(length))
    data.frame(kmer=results[[5]],score=results[[4]])
  
}
