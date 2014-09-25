Sign <- function(x)
{
    if(x>=0) sign <- 1
    if(x < 0) sign <- -1
    return(sign)
}


Eigen <- function(dat)
{
    m1 <- eigen(dat)$values
    m2 <- eigen(dat)$vectors
    m2 <- apply(m2, 2, function(x){x * Sign(Re(x[1]))})
    results <- list("values"=m1,
                    "vectors"=m2)
    return(results)
}

#this calculates the effect size of a moderator 

effectSize <- function(response, treatment, moderator)
{
    mod_scale <- scale(moderator, center=TRUE, scale=TRUE)
    treat_recode <- as.numeric(as.character(factor(treatment, labels=c(-0.5, 0.5))))
    mm <- lm(response~treat_recode + mod_scale + treat_recode * mod_scale)
    
    eff_size <- as.numeric(mm[[1]][4]/2/sqrt(summary(mm)[[6]]^2+mm[[1]][3]^2+mm[[1]][4]^2/4))
    
    return(eff_size)
    
}

g <- function (eta, x)
{
    if(t(eta) %*% x >=0) assignment <- 1
    else assignment <- 0
    return(assignment)
}
