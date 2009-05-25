`duembgen.shape` <-
function(X,init=NULL,steps=Inf,eps=1e-6,maxiter=100,in.R=FALSE,na.action=na.fail,...)
    {
    X<-na.action(X)
    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)
    
    
    p <- dim(X)[2]
    if (p<2) stop("'X' must be at least bivariate")  
    
    data2<-pair.diff(X)
    colnames(data2) <- colnames(X)
    duembgen<-tyler.shape(data2,0,init,steps,eps,maxiter,in.R,...)
    return(duembgen)
    }
