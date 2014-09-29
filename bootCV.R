# this is a piece of code for the bootstrap and cross validation to 
#get ab interval for the estimated statistics





bootCv <- function(trt, response, designMatrix,type)
{
	dat <- cbind(trt,response,designMatrix)
	datList <- dlply(dat,.(trt))
	N <- lapply(datList,nrow)
	
	result <- as.data.frame(matrix(0,10000,6))
	colnames(result) <- c("bootstrap","CV","nGem","nAIPEW","apbGem","apbAIPEW")
	
	
	for (i in 1:1000)
	{
		cat(i)
		set.seed(i)
		datReorder <- llply(datList,function(x){x[sample(1:nrow(x)),]})
		split0 <- split(1:N[[1]],rep(1:10,length.out=N[[1]]))
		split1 <- split(1:N[[2]],rep(1:10,length.out=N[[2]]))
		dat0List <- llply(split0,function(x){datReorder[[1]][x,]})
		dat1List <- llply(split1,function(x){datReorder[[2]][x,]})
		for (j in 1:10)
		{
			datTrain <- rbind(ldply(dat0List[-j]),ldply(dat1List[-j]))
			datTest <-  rbind(ldply(dat0List[j]),ldply(dat1List[j]))
			datTest$.id <- NULL
			datTrain$.id <- NULL
			trtTrain <- datTrain[,1]
			resTrain <- datTrain[,2]
			designMatrixTrain <- datTrain[,c(-1,-2)]
			model <- gem_aipew(trtTrain,resTrain,designMatrixTrain,type,"T")
			trtTestGem <- apply(datTest[,c(-1,-2)],1,function(x){g(c(1,x),model[[1]])})
			trtTestAIPEW <- apply(datTest[,c(-1,-2)],1,function(x){g(c(1,x),model[[3]])})

			result[(i-1)*10+j,1] <- i
			result[(i-1)*10+j,2] <- j
			result[(i-1)*10+j,3] <- sum(datTest$trt==trtTestGem)
			result[(i-1)*10+j,4] <- sum(datTest$trt==trtTestAIPEW)
			result[(i-1)*10+j,5] <- sum(datTest[datTest$trt==trtTestGem,2])/sum(datTest$trt==trtTestGem)
			result[(i-1)*10+j,6] <- sum(datTest[datTest$trt==trtTestAIPEW,2])/sum(datTest$trt==trtTestAIPEW)
		}		
	}
	
	return(result)	
}