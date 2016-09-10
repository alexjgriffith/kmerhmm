#' @export
nsum1<-function(x,n=1)
    apply(x,n,function(y) y/sum(y))

#' @export
forward<-function(t,p,a,b,y,n=1,alpha=c()){    
    N<-length(p)
    if(n==1){
        alphaP<-b[(y[n]+1),]*p
        forward(t,p,a,b,y,n+1,cbind(alpha,alphaP))
    }
    else if(n>1 & n<=t){

        alphaP<-sapply(seq(N),function(i) b[(y[n]+1),i]*sum( alpha[,n-1]*a[i,]))
        forward(t,p,a,b,y,n+1,cbind(alpha,alphaP))
    }
    else if(n==t+1)
        alpha
}

#' @export
backward<-function(t,a,b,y,n=length(y),count=0,beta=c()){
    N=dim(b)[2]
    if(count==0)
        betaP<-rep(1,N)
    else if (count==1)
        betaP<-sapply(seq(N),function(i) sum(beta[i]*a[,i]*b[y[n-count]+1,i]))
    else
        betaP<-sapply(seq(N),function(i) beta[i,count]*sum(a[,i]*b[y[n-count]+1,i]))
    if((n-count)==t)
        t(apply(beta,1,rev))
    else
        backward(t,a,b,y,n,count+1, cbind(beta,betaP))
    
}

#' @export
weightedBW<-function(y,a,b,p){
    M<-dim(b)[1]
    N<-dim(a)[1]
    T<-length(y)
    beta<-backward(0,a,b,y)    
    alpha<-forward(T,p,a,b,y)    
    denom<-sum(alpha[,T])
    gamaIT<-do.call(rbind,lapply(seq(N), function(i) alpha[i,]*beta[i,]/denom))
    zetaIJT<-do.call(abind,append(lapply(seq(length(T)-1) ,function(l) outer(seq(N),seq(N),Vectorize(function(i,j,t=l){alpha[i,t]*a[i,j]*beta[j,t+1]*b[(y[t+1]+1),j]/denom}))),list(along=3)))
    ap<-nsum1(outer(seq(N),seq(N),Vectorize(function(i,j) sum(zetaIJT[i,j,])/sum(gamaIT[i,]))))
    bp<-nsum1(do.call(rbind,lapply(seq(M),function(l) apply(gamaIT,1,function(x)sum(x[y+1==l])/sum(x)))),n=2)
    pp<-nsum1(matrix(gamaIT[,1],ncol= N))
    rbind(ap,bp,t(pp))
}

#' @export
manyWB<-function(sequence,a,b,p){
    Num<-dim(sequence)[1]
    M<-dim(b)[1]
    N<-dim(a)[1]
    s<-do.call(abind,append(lapply(seq(Num),function(i) weightedBW(sequence[i,],a,b,p)),list(along=3)))

    ret<-outer(seq(N+M+1),seq(N), Vectorize(function(i,j)mean(s[i,j,])))
    a<-ret[1:N,]
    b<-ret[(N+1):(N+M),]
    p<-matrix(ret[(N+M+1):dim(ret)[1],],nrow=N)
    list(a,b,p)
}

#' @export
inferenceWB<-function(sequence,a,b,p,n){
    if(n>0){
        new<-manyWB(sequence,a,b,p)
        print(n)
        do.call(inferenceWB,append(list(sequence),append(new,list(n-1))))
    }
    else
        list(a,b,p)        
}

#' @export
buildHMM<-function(sequences,levels,nodes,iterations=5){

    sequencep<-apply(sequences,2,function(x) sapply(x,function(y) which(y==levels)-1))

    a=nsum1(matrix(runif(nodes^2),nrow=nodes))

    b=nsum1(matrix(runif(nodes*length(levels)),nrow=nodes))

    p=nsum1(matrix(runif(nodes),ncol=nodes))

    inferenceWB(sequencep,a,b,p,iterations)
}
