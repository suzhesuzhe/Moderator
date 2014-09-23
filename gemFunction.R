

source("~/Documents/New_York_University/moderator/code/gem4_09:13:2014/dependency.R")
#keep in mind that each argument may not in order, say the trt can be 010010101...
gemFunction<- function(trt, response, designMatrix, type)
{
    #the preprocessing of the data
    K <- length(unique(trt))           #group categories
    p <- ncol(designMatrix)            #p: the number of the predictors  
    
    #create the dat matrix and the variable Treatment
    dat <- as.data.frame(cbind(trt, response, designMatrix))    
    colnames(dat)[1] <- "Treatment"    

    #centered and uncentered list for the baseline predictors
    centeredXList <- dlply(dat, .(Treatment),function(x)
    	{
    		apply(x[,2:ncol(dat)], 2, function(y){y-mean(y)})[,2:(p+1)]
    	})  
    uncenteredXList <-dlply(dat, .(Treatment),function(x){x[,c(-1,-2)]})
    centeredXFrame <- ldply(centeredXList, function(x){x})[,-1]
    uncenteredXFrame <- ldply(uncenteredXList, function(x){x})[,-1]

    #centered list , centered and uncentered vector for the response
    centeredYList <- dlply(dat, .(Treatment), function(x)
    	{
    		apply(x[,2:ncol(dat)], 2, function(y){y-mean(y)})[,1]
    	}) 
    uncenteredYVector <- unlist(dlply(dat, .(Treatment),function(x){x[,2]}))
    centeredYVector <- unlist(centeredYList)    

    #recording the outcome variable for each group
    orderedTrt <- ldply(centeredXList, function(x){x})[,1] 
    
    
    co <- cov(centeredXFrame)
    e <- Eigen(co)
    sqrtco <- e$vectors%*%diag(sqrt(e$values))%*%t(e$vectors)
    N <- ddply(dat, .(Treatment), nrow)[,2]    #a vector recording the number of observation for each group
    pro <- N/sum(N)
    pro1<-(N-1)/sum(N)
    Beta <- vector("list",K)
    for(i in 1:K)
    {
        Beta[[i]] <- as.matrix(lm(centeredYList[[i]]~centeredXList[[i]])$coeff[-1])
    }
    
    Beta.bar<-Reduce('+',lapply(1:K, function(j){pro[j]*Beta[[j]]}))
    B<-Reduce('+',lapply(1:K, function(j)
    	{
    		pro[j] * (Beta[[j]]- Beta.bar) %*% t(Beta[[j]] - Beta.bar)
    	}))
    # B is the matrix for the numerator method   
    D<-Reduce('+',lapply(1:K, function(j){pro1[j] * (Beta[[j]]) %*% t(Beta[[j]])}))
    # D is the matrix for the denominator method   
    A<-Reduce('+',lapply(1:K, function(j)
    	{
    		Reduce('*',(pro1[j]/(N[j]-1))*t(centeredYList[[j]])%*%centeredYList[[j]],diag(p))
    	}))-sqrtco%*%D%*%sqrtco
    # A is the matrix in the denominator for the F method  
    
    ######################  the three different methods  #########
    # maximizing the interaction term (nu)
    if(type=="nu")
    {                
        astar <- Eigen(sqrtco%*%B%*%sqrtco)$vector[,1]
        alpha_nu <- solve(sqrtco)%*%astar
        alpha_nu_2 <- Sign(Re(alpha_nu[1,1]))*Re(alpha_nu)
        alpha<-alpha_nu_2
    }
    
    # minimizing the error sum of squares (de)
    if(type=="de")
    {        
        astar <- Eigen(sqrtco%*%D%*%sqrtco)$vector[,1]
        alpha_de <- solve(sqrtco) %*% astar
        alpha_de_2 <- Sign(Re(alpha_de[1,1]))*Re(alpha_de)
        alpha<-alpha_de_2
    }
    
    # maximizing the F ratio (F)    
    if(type=="F")
    {
        astar<- Eigen(solve(A)%*%sqrtco%*%B%*%sqrtco)$vector[,1]
        alpha_F <- solve(sqrtco)%*%astar
        c <- sqrt(as.numeric(t(alpha_F)%*%co%*%alpha_F))
        alpha_F <- Re(alpha_F)/c
        
        alpha_F_2 <- Sign(Re(alpha_F[1,1])) * Re(alpha_F)
        alpha<-alpha_F_2
    }
    
    # compute the angle in case of K=2
    if (p==2) 
    {
        if(alpha[2,1]>=0) angle<- atan(alpha[2,1]/alpha[1,1])
        if(alpha[2,1]<0) angle <- atan(alpha[2,1]/alpha[1,1])+pi
    }
    if(p!=2) angle <- NA
    
    # The centered model for the combined data       
    ZZinter <- matrix(0, sum(N), K)
    ZZinter[1:N[1],1]<-1
    for(j in 2:K)
    {
        ZZinter[(sum(N[1:(j-1)])+1):sum(N[1:j]), j] <- 1
    }   
    
    centeredZZ <- matrix(0, sum(N), K)
    centeredZZ[1:N[1],1]<- as.matrix(centeredXList[[1]]) %*% alpha
    for(j in 2:K)
    {
        centeredZZ[(sum(N[1:(j-1)])+1):sum(N[1:j]), j] <- as.matrix(centeredXList[[j]]) %*% alpha
    }   
    centeredZZ_red <- as.matrix(apply(centeredZZ,1,sum))
    
    centeredGemReduced <- lm(uncenteredYVector ~ -1+cbind(ZZinter,centeredZZ_red))
    centeredGemFull <- lm(uncenteredYVector~-1+cbind(ZZinter,centeredZZ))
    centeredAnova <- anova(centeredGemReduced, centeredGemFull)
    
    # The uncentered model for the combined data       
    uncenteredZZ <- matrix(0, sum(N), K)
    uncenteredZZ[1:N[1],1]<- as.matrix(uncenteredXList[[1]]) %*% alpha
    for(j in 2:K)
    {
        uncenteredZZ[(sum(N[1:(j-1)])+1):sum(N[1:j]), j] <- as.matrix(uncenteredXList[[j]]) %*% alpha
    }   
    uncenteredZZ_red <- as.matrix(apply(uncenteredZZ,1,sum))
    
    uncenteredGemReduced <- lm(uncenteredYVector ~ -1+cbind(ZZinter,uncenteredZZ_red))
    uncenteredGemFull <- lm(uncenteredYVector~-1+cbind(ZZinter,uncenteredZZ))
    uncenteredAnova <- anova(uncenteredGemReduced, uncenteredGemFull)
    
    
    
    
    #we <- lm(uncenteredY~orderedTrt*ZZ_red)
    #cbind(1,0,ZZ_red,0)%*%we$coef
    #cbind(1,1,ZZ_red,ZZ_red)%*%we$coef
	if(K==2)
	{
		centeredValue0 <- cbind(1,centeredZZ_red) %*% centeredGemFull$coef[c(1,K+1)]
    	centeredValue1 <- cbind(1,centeredZZ_red) %*% centeredGemFull$coef[c(2,K+2)]
		
		uncenteredValue0 <- cbind(1,uncenteredZZ_red) %*% uncenteredGemFull$coef[c(1,K+1)]
    	uncenteredValue1 <- cbind(1,uncenteredZZ_red) %*% uncenteredGemFull$coef[c(2,K+2)]
		
		
		centeredOptGemAss <- rep(0,nrow(designMatrix))
    	centeredOptGemAss[centeredValue1>=centeredValue0] <- 1
		
		uncenteredOptGemAss <- rep(0,nrow(designMatrix))
    	uncenteredOptGemAss[uncenteredValue1>=uncenteredValue0] <- 1
			
    	centeredOptProTrt <- sum(orderedTrt==centeredOptGemAss)/length(trt)
    	uncenteredOptProTrt <- sum(orderedTrt==uncenteredOptGemAss)/length(trt)

		centeredOptPopuAverage <- sum(uncenteredYVector[orderedTrt==centeredOptGemAss])/sum(orderedTrt==centeredOptGemAss)
		uncenteredOptPopuAverage <- sum(uncenteredYVector[orderedTrt==uncenteredOptGemAss])/sum(orderedTrt==uncenteredOptGemAss)

		centeredIntersectGem <- (centeredGemFull$coef[1]-centeredGemFull$coef[2])/(centeredGemFull$coef[4]-centeredGemFull$coef[3])
		uncenteredIntersectGem <- (uncenteredGemFull$coef[1]-uncenteredGemFull$coef[2])/(uncenteredGemFull$coef[4]-uncenteredGemFull$coef[3])

		effect_size <- abs(effectSize(centeredYVector,orderedTrt,centeredZZ_red))
	}
    if(K!=2)
    {
    	optProAss <- NA
    	optPopuAverage <- NA
    	constraintGem <- NA
    	 effect_size <- NA
    }
    
    
    centeredPvalue=centeredAnova[2,6]
    uncenteredPvalue=uncenteredAnova[2,6]


    #XX <- matrix(0, sum(N), K*p)
    #XX[1:N[1],1:p]<- X[[1]]
    #for(j in 2:K)
    #{
    #    XX[(sum(N[1:(j-1)])+1):sum(N[1:j]), ((j-1)*p+1):(j*p)] <- X[[j]]
    #}
    
    #XX_red <-as.matrix(ldply(X, function(x){x})[,-1])
    #fitfull <- lm(uncenteredY ~ cbind(ZZinter,XX))
    #fitreduced <- lm(uncenteredY ~ XX_red)
    #anovafit <- anova(fitfull,fitreduced)
    
    results <- list("type"= type,      #1 The method for the gem implementation
                    "alpha"=alpha,     #2 The weight
                    "angle"=angle,     #3 if two baseline, this is the corresponding angle between them
                    "centeredPvalue"=centeredPvalue,           #4 the p value for whether the centered combined baseline is significant
                    "uncenteredPvalue"=uncenteredPvalue,       #5 the p value for whether the uncentered combined baseline is significant
    				"effect.size"=effect_size,                 #6 the effect size of the combined baseline, this model's response is centered
                    #"p.X" = anovafit[2,6],#6 the p value for the interaction with multiple baseline predictor
                    "centeredGemFull" = centeredGemFull,       #7 the model of centered combined baseline with uncentered response
    				"uncenteredGemFull"=uncenteredGemFull,     #8 the model of uncentered combined baseline with uncentered response
    				
    				"optProTrt_centered"=centeredOptProTrt,                #9  the coincidence proportion for observed treatment with the gem criteria(centered designMatrix)
    				"optPopuAverage_centered"= centeredOptPopuAverage,     #10 the population average for those treatment assignment coincide with gem criteria(centered designMatrix)
    				"intersectGem_centered"=centeredIntersectGem,          #11 the threshold for the combined baseline in terms of the uncentered y(centered designMatrix)
    				
    				"optProTrt_uncentered"=uncenteredOptProTrt,            #12  the coincidence proportion for observed treatment with the gem criteria(uncentered designMatrix)
    				"optPopuAverage_uncentered"= uncenteredOptPopuAverage, #13 the population average for those treatment assignment coincide with gem criteria(uncentered designMatrix)
    				"intersectGem_uncentered"=uncenteredIntersectGem,      #14 the threshold for the combined baseline in terms of the uncentered y(uncentered designMatrix)
    				
    				"centeredExtendData" =cbind("trt"=orderedTrt,uncenteredYVector,centeredXFrame,
    									"centeredZ"=centeredZZ_red,
    									"centeredOpGemTrt"=centeredOptGemAss), #15 dataset with uncentered y, centered x, centered z, ordered treatment, optimal treatment 
    				"uncenteredExtendData" =cbind("trt"=orderedTrt,uncenteredYVector,uncenteredXFrame,
    									"uncenteredZ"=uncenteredZZ_red,
    									"uncenteredOpGemTrt"=uncenteredOptGemAss)) #16 dataset with uncentered y, uncentered x, uncentered z, ordered treatment, optimal treatment                  
    
    return(results)
}




