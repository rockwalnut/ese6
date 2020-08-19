### Project : Modified Simulated Annealing Algorithm 
### Filename : MSA.R
### Created : 19.08.2009
### Last edited : 10.01.2010

randomMatrix <- function(run, var)														### random LHD
{
    LHD <- matrix(0, run, var); 
    x <- seq(0, 1, (1 - 0)/(run - 1));  
    for(i in 1 : var) 
    {
        LHD[ ,i] <- sample(x);
    }	
    return(LHD);                
}

euclideanMatrix <- function(LHD)														### generate euclidean matrix 
{
    D <- matrix(NA, dim(LHD)[1], dim(LHD)[1])
    for(i in 1 : dim(LHD)[1]) 
    {
		for(j in i : dim(LHD)[1]) 
		{
			if(i == j) 
			{
				sum <- NA;  
			} 
			else 
			{
				sum <- 0;
				for(k in 1 : dim(LHD)[2]) 
				{
					sum <- sum + ((LHD[i,k] - LHD[j,k])^2)
				}  
			}
            D[i,j] <- sqrt(sum);    
        }
    } 
    return(D);      
}

calculateMaximin <- function(D)															### find maximin 
{
    TEMP <- D;			
    TEMP[is.na(TEMP)] <- 9;	
    return(min(TEMP));	
}

calculatePhiP <- function(LHD)															### calculate phi-p 
{
    p <- 5;
    D <- euclideanMatrix(LHD);
    TEMP <- D;
    TEMP[is.na(TEMP)] <- 0;
    TEMP <- 1/TEMP;	
    TEMP[is.infinite(TEMP)] <- 0;
    TEMP <- TEMP^p;
    PhiP <- sum(TEMP);
    PhiP <- PhiP^(1/p);
    list(PhiP = PhiP, D = D);   	
}

randomRows <- function(LHD)																### exechange-element
{
    rcol <- (runif(1)*ncol(LHD))+1; 		 
    randRow <- sample(dim(LHD)[1],2);		
    RAND <- c(randRow[1], randRow[2], rcol); 
    RAND <- floor(RAND);
    return(RAND);   
}

elementExchenges <- function(LHD, col, row1, row2)										### exechange-element to create new LHD 
{
    m_col <- col; 		 
    m_row1 <- row1;		 
    m_row2 <- row2; 		
    NLHD <- LHD;

    while(NLHD[m_row1,m_col] == NLHD[m_row2,m_col]) 
    {
        m_row2 <- (runif(1)*nrow(LHD))+1;       
    }

    temp <- NLHD[m_row1,m_col];
    NLHD[m_row1,m_col] <- NLHD[m_row2,m_col];   
    NLHD[m_row2,m_col] <- temp;
    return(NLHD);   
}

calculatePhiPNew <- function(LHD, D, RAND, PhiP)										### calculate new phi-p by new method
{
    p <- 5;
    D[is.na(D)] <- 0;
    row1 <- RAND[1];
    row2 <- RAND[2];
    rcol <- RAND[3];
    sum1 <- 0;
    sum2 <- 0;

    for(loop in 1 : dim(LHD)[1]) 
    {
        if(loop != row1 && loop != row2) 
		{
			d <- abs(LHD[row2,rcol] - LHD[loop,rcol])^2 - abs(LHD[row1,rcol] - LHD[loop,rcol])^2; 
			if(loop < row1) 
			{
				dDat11 <- sqrt(D[loop, row1]^2 + d);
				sum1 <- sum1 + (dDat11^(-p) - D[loop, row1]^(-p)); 		
				D[loop, row1] <- dDat11; 
            } 
			else 
			{
				dDat12 <- sqrt(D[row1, loop]^2 + d);
				sum1 <- sum1 + (dDat12^(-p) - D[row1, loop]^(-p));
				D[row1, loop] <- dDat12;
			}    
			if(loop < row2) 
			{
				dDat21 <- sqrt(D[loop, row2]^2 - d);	
				sum2 <- sum2 + (dDat21^(-p) - D[loop, row2]^(-p));
				D[loop, row2] <- dDat21;
            } 
			else 
			{
				dDat22 <- sqrt(D[row2, loop]^2 - d);	
				sum2 <- sum2 + (dDat22^(-p) - D[row2, loop]^(-p));
				D[row2, loop] <- dDat22;
			}  
		}
    }
    DTry <- D;
    PhiPTry <- (PhiP^p + sum1 + sum2)^(1/p);     
    list(PhiPTry = PhiPTry, DTry = DTry);   
}

MSAOutter <- function(LHD)																### MSA Main
{
    m_iteration <- 500;   					
    c_rate <- 1;		
	tol <- 0.0001;
	
    D <- euclideanMatrix(LHD);				
    TEMP <- D;
    TEMP[is.na(TEMP)] <- 0;
    temperature <- mean(TEMP);																		

    lhd_best <- LHD;
    TEMP <- calculatePhiP(LHD);
    phi_p <- TEMP$PhiP;						
    phip_best <- phi_p;

    flag <- TRUE;
    label <- FALSE; 

	Time <- 0;
	LoopCount <- 0;

    while(flag == TRUE || label == TRUE) 
    { 
        temperature <- temperature * c_rate;	
 		c_rate <- 0.95;
        label <- FALSE; 
        flag <- FALSE; 
		
		i <- 1;
        while(i < m_iteration) 
		{ 	
			LoopCount <- LoopCount + 1;
			RAND <- randomRows(LHD);
			ptm <- proc.time();
			lhd_try  <- elementExchenges(LHD, RAND[3], RAND[1], RAND[2]);				
			temp <- calculatePhiPNew(LHD, D, RAND, phi_p);								### calculate new phi-p with new method	
			phip_try <- temp$PhiPTry; 
			d_try <- temp$DTry;

            if(phi_p - phip_try >= tol || rbinom(1,1,exp(-1*((abs(phip_try - phi_p))/temperature))) == 1) 
			{
				phi_p <- phip_try;
				D <- d_try;
				LHD <- lhd_try;															### replace old LHD if phi-p new is better than
  				label <- TRUE;
				flag <- TRUE;
            } 

			if(phip_try < phip_best)													### track best desig so far
			{
				i <- 1;
				lhd_best <- lhd_try;
				phip_best <- phip_try;
			} 
			else 
			{ 
                i <- i + 1;  
 			}
			gap <- proc.time() - ptm;
			Time <- Time + gap[1];
        }
    }
    list(PhiP = phip_best, Time = Time, Loop = LoopCount);      
}

LHDBestMSA <- function(run, var)														### write output file
{
    PhiPBest <- matrix(0, 10, 1);
    TimeElapsed <- matrix(0, 10, 1);
    LoopExchange <- matrix(0, 10, 1);

    for(i in 1 : 10)
    {
        LHD <- randomMatrix(run, var);
        temp <-  MSAOutter(LHD);
        PhiPBest[i, 1] <- temp$PhiP;
        TimeElapsed[i, 1] <- temp$Time;
		LoopExchange[i, 1] <- temp$Loop;

        resultMSA <- data.frame(PhiPBest, TimeElapsed, LoopExchange);
    }

	fileName <- paste("e:\\MSAdata", substr(date(), 5, 7),substr(date(), 9, 10), ".csv", sep="") 
	write.table(substr(date(), 5, 19), file = fileName, append = FALSE, row.names = FALSE, col.names = FALSE);
    write.table(resultMSA, file = fileName, append = TRUE, sep = ",", row.names = FALSE, col.names = c("Phi-P Best", "Time Elapsed", "Loop Exchange"));
    return(resultMSA);
}