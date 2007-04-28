`rank.ctest`<-function(X,...)
    {
    UseMethod("rank.ctest")
    }

`rank.ctest.default` <-
function(X,Y=NULL,mu=NULL,scores="rank",na.action=na.fail,...)
    {
    if (is.null(Y)) 
       {
       DNAME<-deparse(substitute(X))
       }
    else
       {
       DNAME<-paste(deparse(substitute(X)),"and",deparse(substitute(Y)))
       }
       
    p<-dim(X)[2]
   
    if (is.null(mu)) mu<-rep(0,p)
    else if (length(mu)!=p) stop("length of 'mu' must equal the number of columns of 'X'")

    X<-na.action(X)

    if(!all(sapply(X, is.numeric))) stop("'X' must be numeric")
    X<-as.matrix(X)

    if (!is.null(Y)) 
        {
        Y<-na.action(Y)
        if(!all(sapply(Y, is.numeric))) stop("'Y' must be numeric")
        Y<-as.matrix(Y)
        if (p!=dim(Y)[2]) stop("'X' and 'Y' must have the same number of columns")
        if (dim(X)[1]<2 | dim(Y)[1]<2) stop("both 'X' and 'Y' must have at least two observations")   
        }
    else
        {
        if (dim(X)[1]<2) stop("'X' must have at least two observations")  
        }
     

    X<-sweep(X,2,mu)
    
    scores<-match.arg(scores,c("sign","rank","normal"))
    
    if (is.null(Y) & scores=="sign") version<-"one.sample.sign"
    if (is.null(Y) & scores=="rank") version<-"one.sample.rank"
    if (is.null(Y) & scores=="normal") version<-"one.sample.normal"
    if (!is.null(Y) & scores=="sign") version<-"two.sample.sign"
    if (!is.null(Y) & scores=="rank") version<-"two.sample.rank"
    if (!is.null(Y) & scores=="normal") version<-"two.sample.normal"
    
    res1<-switch(version,
            "one.sample.sign"={
                result<-MaRaTe.internal(X,test=scores)
                STATISTIC<-result$test.statistic
                names(STATISTIC)<-"T"
                PVAL<-result$p.value
                METHOD<-"Marginal One Sample Sign Test"
                PARAMETER<-p
                names(PARAMETER)<-"df"
                RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
               
               RVAL}
            ,
            "one.sample.rank"={
                result<-MaRaTe.internal(X,test=scores)
                STATISTIC<-result$test.statistic
                names(STATISTIC)<-"T"
                PVAL<-result$p.value
                METHOD<-"Marginal One Sample Signed Rank Test"
                PARAMETER<-p
                names(PARAMETER)<-"df"
                RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
               
               RVAL}
            ,
            "one.sample.normal"={
                result<-MaRaTe.internal(X,test=scores)
                STATISTIC<-result$test.statistic
                names(STATISTIC)<-"T"
                PVAL<-result$p.value
                METHOD<-"Marginal One Sample Normal Scores Test"
                PARAMETER<-p
                names(PARAMETER)<-"df"
                RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
               
               RVAL}
            ,
            "two.sample.sign"={
                result<-MaRaTe.internal(X,Y,scores)
                STATISTIC<-result$test.statistic
                names(STATISTIC)<-"T"
                PVAL<-result$p.value
                METHOD<-"Marginal Two Sample Median Test"
                PARAMETER<-p
                names(PARAMETER)<-"df"
                RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
               
               RVAL}
            ,
            "two.sample.rank"={
                result<-MaRaTe.internal(X,Y,scores)
                STATISTIC<-result$test.statistic
                names(STATISTIC)<-"T"
                PVAL<-result$p.value
                METHOD<-"Marginal Two Sample Rank Sum Test"
                PARAMETER<-p
                names(PARAMETER)<-c("df")
                RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
               
               RVAL}
            ,
            "two.sample.normal"={
                result<-MaRaTe.internal(X,Y,scores)
                STATISTIC<-result$test.statistic
                names(STATISTIC)<-"T"
                PVAL<-result$p.value
                METHOD<-"Marginal Two Sample Normal Scores Test"
                PARAMETER<-p
                names(PARAMETER)<-c("df")
                RVAL<-list(statistic=STATISTIC,p.value=PVAL,method=METHOD,parameter=PARAMETER)
               
               RVAL}
            
            )
    ALTERNATIVE="two.sided"
    NVAL<-paste("c(",paste(mu,collapse=","),")",sep="")
    if (is.null(Y)==TRUE) names(NVAL)<-"location" else names(NVAL)<-"location difference"
    res<-c(res1,list(data.name=DNAME,alternative=ALTERNATIVE,null.value=NVAL))
    class(res)<-"htest"
    return(res)
    }

`rank.ctest.formula` <-
function(formula, na.action = na.fail, ...)
   {
    if(missing(formula)
       || (length(formula) != 3)
       || (length(attr(terms(formula[-2]), "term.labels")) != 1))
        stop("'formula' missing or incorrect")
    m <- match.call(expand.dots = FALSE)
    if(is.matrix(eval(m$data, parent.frame())))
        m$data <- as.data.frame(data)
    m[[1]] <- as.name("model.frame")
    m$... <- NULL
    mf <- eval(m, parent.frame())
    DNAME <- paste(names(mf), collapse = " by ")
    names(mf) <- NULL
    response <- attr(attr(mf, "terms"), "response")
    g <- factor(mf[[-response]])
    if(nlevels(g) != 2)
        stop("grouping factor must have exactly 2 levels")
    DATA <- split(as.data.frame(mf[[response]]), g)
    names(DATA) <- c("X", "Y")
    y <- do.call("rank.ctest", c(DATA, list(...)))
    y$data.name <- DNAME
    return(y)
  }

rank.ctest.ics <-
function(X, g = NULL, index = NULL, na.action = na.fail, ...)
    {
    
    if (!is.null(g)) DNAME<-paste(deparse(substitute(X)),"and",deparse(substitute(g)))
    else DNAME<-deparse(substitute(X))
    if (class(X) != "ics") stop("'ics' must be of class 'ics'") 
    if (!is.null(g) & !is.factor(g)) stop("'g' must be a factor with two levels")
    if (!is.null(g) & nlevels(g) != 2) stop("'g' must be a factor with two levels")
    
    Z <- ics.components(X)
    p <- dim(Z)[2]
    if (is.null(index)) index <-  1:p
    if(is.null(g))
        {
        DATA<-list(X=Z[,index])
        y <- do.call("rank.ctest", c(DATA, list(...)))
        y$data.name <- DNAME
        }
    if(!is.null(g))
        {
        DATA <- split(Z[,index], g)
        names(DATA) <- c("X", "Y")
        y <- do.call("rank.ctest", c(DATA, list(...)))
        y$data.name <- DNAME
        }
    return(y)    
    }
