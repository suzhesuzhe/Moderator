
pvalue_permute <- function(trt, response, designMatrix, type)
{
    
    dat1 <- cbind(trt, response, designMatrix)
    colnames(dat1)[1] <- "Treatment"
    dat1[,1] <- factor(dat1[,1],labels=c(0,1))
    N <- ddply(dat1, .(Treatment), nrow)[,2]
    mm <- gemFunction(trt, response, designMatrix, type)
    p <- c(mm[[4]])
  
    for(i in 1:5000)
    {
        
        
        #permutation part
        sel <- sample(1:sum(N),N[1], replace=F)
        new_trt <- rep(1,sum(N))
        new_trt[sel] <- 0
        
        mmm <- gemFunction(new_trt, response, designMatrix, type)
        p <- c(p, mmm[[4]])
        
    }
    p_value <- sum(p[2:5001] <= p[1])/5000
    return(p_value)
    
}

