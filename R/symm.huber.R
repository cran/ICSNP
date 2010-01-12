symm.huber <- function(X, qg=0.9, init=NULL, eps = 1e-06, maxiter = 100, na.action=na.fail)
    {
    X <- na.action(X)
    if (!all(sapply(X, is.numeric))) 
        stop("'X' must be numeric")
    X <- as.matrix(X)
    CNAMES <- colnames(X)
    p <- ncol(X)
    
    c.square <- 2 * qchisq(qg, p)
    sigma.square <- 2 * pchisq(c.square / 2, p + 2) + (c.square / p) * (1 - qg)
    
    if (p<2) stop("'X' must be at least bivariate")  
    
    X<-pair.diff(X)
    n<-dim(X)[1]

    if (is.numeric(init)) V.0 <- init else V.0<- crossprod(X) / n
 
    iter<-0
    differ<-Inf
    
    while(differ>eps)
        {
        if (iter==maxiter)
            {
             stop("maxiter reached without convergence")
            }
        iter<-iter+1
        d2<- mahalanobis(X,0,V.0)
        w2<-ifelse(d2<=c.square, 1/sigma.square, (c.square/d2)/sigma.square)
        w.X <- sqrt(w2) * X
        V.new<- (1/n) * crossprod(w.X)
        differ <- ICS:::frobenius.norm(V.new-V.0)    
        V.0<-V.new
        }
    colnames(V.new) <- CNAMES
    rownames(V.new) <- CNAMES
    V.new
    }
