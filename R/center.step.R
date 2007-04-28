`center.step` <-
function(datas,A,maxiter,eps.center)                                            
    {
     solve(A)%*%spatial.median(datas%*%t(A),init=NULL,maxiter=maxiter,eps=eps.center,print.it=FALSE) 
    }

