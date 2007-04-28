`spatial.median.step` <-
function(y,datas)
    { 
    eta.y <- 1-as.numeric(apply(datas,1,setequal,y))
    if (sum(eta.y)==0) stop("all observations are equal")
    datas.y <- sweep(datas,2,y) 
    d.y <- norm(datas.y)
    T.tilde.y <- (sum(1/d.y[eta.y==1]))^(-1)*colSums((sweep(datas[eta.y==1,],1,d.y[eta.y==1],"/")))
    R.tilde.y <- colSums((sweep(datas.y[eta.y==1,],1,d.y[eta.y==1],"/")))
    r.y <- sqrt(sum(R.tilde.y^2))
    if (length(eta.y)==length(eta.y[eta.y==1])) y.new <- T.tilde.y
    else y.new <- max(c(0,(1-sum(eta.y)/r.y)))*T.tilde.y+ min(c(0,(sum(eta.y)/r.y)))*y
    return(y.new)
    }

