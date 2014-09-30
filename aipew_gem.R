
#this is a general function for the gem

#it should has the same argument with the gemFuncion.R

require("RcppArmadillo")
g <- function (eta, x)
{
    if(t(eta) %*% x >=0) assignment <- 1
    else assignment <- 0
    return(assignment)
}

gem_aipew <- function(trt, response, designMatrix,type, both)
{
     
    gem <- gemFunction(trt, response, designMatrix,type)
    dat <- gem[[16]]   #This selected uncentered designMatrix determine that this eta criteria is uncentered
    dat[,1] <- factor(dat[,1])
    mod <- lm(uncenteredYVector~trt+uncenteredZ+trt:uncenteredZ,data=dat)
    
    alpha_gem <- c(mod[[1]][2],mod[[1]][4]*gem[[2]])
	signZ <- Sign(mod[[1]][4])
    alpha_gem <- alpha_gem/sqrt(sum(alpha_gem^2))
    if (both=="F")
    {
    	alpha_aipew <-NA
    	aipew <-NA
    	constraintAIPEW <- NA
    }
    if (both=="T")
	{
    	eta_matrix <- cbind(seq(-1,1,0.01),0,0)
    	eta_matrix[,2] <- signZ * sqrt(1-eta_matrix[,1]^2)
    	for(i in 1:nrow(eta_matrix))
    	{
    		
    		datt <- dat
    		eta <- eta_matrix[i,1:2]
    		datt$A_opt <- apply(datt,1,function(x){g(eta, c(1, as.numeric(x["uncenteredZ"])) )})
    		datt$C_eta <- (as.numeric(datt[,1])-1)*datt$A_opt + (2-as.numeric(datt[,1]))*(1-datt$A_opt)
    		
    		sm <- as.matrix(cbind(1,datt$uncenteredZ,as.numeric(datt[,1])-1,datt$uncenteredZ*(as.numeric(datt[,1])-1)))
    		yy <- datt[,2]
    		posit_misspecified <- fastLm(sm,yy)
    		
    		mu_1_misspecified <- cbind(1,datt$uncenteredZ,1,datt$uncenteredZ) %*% posit_misspecified[[1]]
    		mu_0_misspecified <- cbind(1,datt$uncenteredZ,0,0) %*% posit_misspecified[[1]]
    		
    		m_misspecified <- mu_1_misspecified *datt$A_opt + mu_0_misspecified *(1-datt$A_opt)
    		
    		eta_matrix[i,3] <- sum(2*datt$C_eta*datt[,2]-(2*datt$C_eta-1)*m_misspecified)/nrow(datt)
    	}
    	
    	s <- which(eta_matrix[,3]==max(eta_matrix[,3]))[1]
    	aipew <- eta_matrix[s,3]
    	
    	alpha_aipew <- c(eta_matrix[s,1], eta_matrix[s,2]*gem[[2]])
    	alpha_aipew <- alpha_aipew/sqrt(sum(alpha_aipew^2))    
    	constraintAIPEW= -eta_matrix[s,1]/eta_matrix[s,2]    ###warning, this constraint is indeed for the uncentered Z
    	dat$opAIPEWTrt <- apply(cbind(1,dat[,c(-1,-2,-ncol(dat),-ncol(dat)+1)]),1,function(x){g(alpha_aipew,x)})

	}   
    
    #The following part give the population averge under the gem funciton criteria


    result <- list("eta_gem" = alpha_gem,  # the transformed gem criteria(this is from the uncentered designMatrix )
     			   "alpha_gem" = gem[[2]],   # the original weight by the gem
                   "eta_aipew" = alpha_aipew, # the weight by aipew
                   "intersectGem" = gem[[14]],  #constraint for the gem with the uncenteredZ
                   "intersectAIPEW"= constraintAIPEW,
    			   "E_aipew" = aipew,
    			   "optPopuAverage"=gem[[13]],
    			   "extendData"=dat,
    			   "gemObject" =gem)
    return(result)
}




####THE FOLLOWING FUNCTION IS NOT GENERAL ONE, BASICALLY FOR THE CALCULATION OF THE SIMULATED DATA OF EXPONENTIAL TYPE

gem_aipew_exponetial <- function(type)
{
    datt <- data
    gem <- gem_aipew(datt[,4],datt[,9],datt[,1:3],type,"T")
  
    eta_gem <- gem[[1]]
    gem_result <- eta_opt(datt,eta_gem)
    p_gem <- sum(gem_result[[1]][,10]==true_opt)/500
    
    eta_aipew <- gem[[3]]
    aipew_result <- eta_opt(datt,eta_aipew)
    p_aipew <- sum(aipew_result[[1]][,10]==true_opt)/500
    
    result <- list("E_gem" = gem_result[[2]],
                   "eta_gem" = eta_gem,
                   "E_aipew" = aipew_result[[2]],
                   "eta_aipew" = eta_aipew,
                   "p_gem" = p_gem,
                   "p_aipew" = p_aipew)
    return(result)
}
