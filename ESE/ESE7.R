### Project : Modified Enhanced Stochastic Evolutionary Algorithm 
### Filename : ESE7.R
### Created :
### Last edited : 

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
	
euclideanMatrix <- function(LHD)														### euclidean matrix 
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

elementExchenge <- function(LHD, col, row1, row2)										### element-exchange to create new LHD 
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
    d_try <- D;
    PhiPTry <- (PhiP^p + sum1 + sum2)^(1/p);     
    list(PhiPTry = PhiPTry, d_try = d_try);   
}

checkContain <- function(array, item)
{
    for(i in 1 : dim(array)[1])
    {
        if((item %in% array[i,])[1] == TRUE && (item %in% array[i,])[2] == TRUE)
		{
			return(TRUE); 
		}
    }
    return(FALSE);
}

ESEInner <- function(LHD, phiP, D, th)													
{										
    n_acpt <- 0;
    n_imp <- 0;

   	j_distinct <- floor(((dim(LHD)[1]*(dim(LHD)[1]-1))/2)/5);
	if(j_distinct > 50)
    {
        j_distinct <- 50;
    }

	inn_iteration <- (2*((dim(LHD)[1]*(dim(LHD)[1]-1))/2)*dim(LHD)[2])/j_distinct;					##(2*(dim(LHD)[1]*dim(LHD)[2]))/j;
    if(inn_iteration > 50)
    {
        inn_iteration <- 50;
    }

    lhd_best <- LHD;																	### best LHD
    phip_best <- phiP;
    d_best <- D;
	timeInner <- 0;

    for(i in 1 : inn_iteration)															### start M iteration 
    {
        col <- (i %% (dim(LHD)[2]) + 1);									
		ptm <- proc.time();
		phip_try_best <- Inf;
        for(k in 1 : j_distinct)														### update elements and choose best LHD try
		{  	
			randRow <- sample(dim(LHD)[1],2);
			RAND <- c(randRow[1], randRow[2], col); 	
			lhd_try <- elementExchenge(LHD, col, randRow[1], randRow[2]);
			temp_phip_try <- calculatePhiPNew(LHD, D, RAND, phiP);						### use new method to calculate new phi-p
			phip_try <- temp_phip_try$PhiPTry;
			d_try <- temp_phip_try$d_try;

			##tol_inn <- 0.001 ;														

			if(phip_try < phip_try_best)										
			{
				phip_try_best <- phip_try;
				lhd_try_best <- lhd_try;
				d_try_best <- d_try;
			}
		}

		if((phip_try_best - phiP) <= (th*runif(1)))										### replace current LHD
		{
			LHD <- lhd_try_best ;
			D <- d_try_best;
			phiP <- phip_try_best;
			n_acpt <- n_acpt + 1;
	    
			if(phiP < phip_best)														### replace best LHD
			{						
				lhd_best <- LHD;								
                phip_best <- phiP;
                d_best <- D;
				n_imp <- n_imp + 1;
			}	
		}
		gap <- proc.time() - ptm;
		timeInner <- timeInner + gap[1]; 
	}
    list(LHD = LHD, D = D, phiP = phiP, lhd_best = lhd_best, phip_best = phip_best, d_best = d_best, n_acpt = n_acpt, n_imp = n_imp, loopCount = inn_iteration * j_distinct, TimeInner = timeInner, inn_iteration = inn_iteration); 
}

ESEOuter <- function(lhd_ini)															
{
    temp_phip <- calculatePhiP(lhd_ini);												### calculate phi-p with initial LHD
    phip_ini <- temp_phip$PhiP;
    d_ini <- temp_phip$D;
	
    LHD <- lhd_ini;																		### current LHD
    phiP <- phip_ini;
    D <- d_ini;

    lhd_best <- LHD;																	### best LHD
    phip_best <- phiP;
    d_best <- D;

    phip_global_best <- phip_best; 

    th_ini <- phip_ini*0.005;
    th <- th_ini;

    n_acpt <- 0;
    n_imp <- 0;

    tol_out <- 0.0001;
    out_iteration <- 20;

	cf <- 0;																			

    alpha_1 <- 0.8;
    alpha_2 <- 0.9 + cf;
    alpha_3 <- 0.7 - cf;

	beta_1 <- 0.1;																		
	beta_2 <- 0.8;
	
    flag_imp <- 1;
    flag <- 1;
    LoopExchange <- 0;
	Time <- 0;																			
	i <- 0;
	flag_loop <- 0;
	notImprove <- 10;
	j <- 0;

    while(i < out_iteration && j < notImprove)															### outter loop iteration 
    {
		ptm <- proc.time();
        if(flag == 1)
		{	
            phip_old_best <- phip_best;													### old best LHD   
			flag <- 0;	
		}
		else
		{
			phip_old_best <- temp_inner$phip_best;
		}

		gap <- proc.time() - ptm;
		temp_inner <- ESEInner(LHD, phiP, D, th);										### call Inner loop
		inn_iteration_out <- temp_inner$inn_iteration;
		ptm <- proc.time();
		TimeInner <- temp_inner$TimeInner;
		LoopExchange <- LoopExchange + temp_inner$loopCount;
		LHD <- temp_inner$LHD;
		D <- temp_inner$D;
		phiP <- temp_inner$phiP;

		n_acpt <- temp_inner$n_acpt;
		n_imp <- temp_inner$n_imp;
		
		phip_best <- temp_inner$phip_best;

        if(phip_old_best - phip_best > tol_out)											### improving process/exploration process
        {
            flag_imp <- 1; 
        }
        else
        {
			flag_imp <- 0; 
        }
		if( flag_imp == 1)																### update th on improving process
		{
			if((n_acpt/inn_iteration_out) > 0.1)
			{
				if((n_imp/inn_iteration_out) < (n_acpt/inn_iteration_out))
				{
					th = th*alpha_1;   
				}
				else if((n_imp/inn_iteration_out) == (n_acpt/inn_iteration_out))
				{
					th = th*1;   
				}
				else
				{
					th = th/alpha_1; 
				}                
			}
			else
			{
					th = th/alpha_1; 
			}                
		}
		else																			### update th on exploration process
		{
			if((n_acpt/inn_iteration_out) >= beta_1 && (n_acpt/inn_iteration_out) <= beta_2 && flag_loop == 0)
			{
				th = th/alpha_3;
			}
			else if((n_acpt/inn_iteration_out) >= beta_1 && (n_acpt/inn_iteration_out) <= beta_2 && flag_loop == 1)
			{
				th = th/alpha_3;
			}
			else if((n_acpt/inn_iteration_out) >= beta_1 && (n_acpt/inn_iteration_out) <= beta_2 && flag_loop == 2)
			{
				th = th*alpha_2;
			}	
			else if((n_acpt/inn_iteration_out) < beta_1)
			{
				th = th/alpha_3;
				flag_loop <- 1;
			}	
			else if((n_acpt/inn_iteration_out) > beta_2)
			{
				th = th*alpha_2;
				flag_loop <- 2;
			}		
		}
		
		if(phip_best < phip_global_best)												### keep track of best design
		{
			phip_global_best <- phip_best;												### modified 2 -> stopping rule add 
			j <- 0;
		}
		else
		{	
			j <- j + 1;
		}

		i <- i + 1;
		gap <- gap + (proc.time() - ptm);
		Time <- Time + gap[1] + TimeInner;
    } 
    list(PhiP = phip_global_best, Time = Time, Loop = LoopExchange);
}

LHDBestESE <- function(run, var)														### write output file
{
    PhiPBest <- matrix(0, 10, 1);
    TimeElapsed <- matrix(0, 10, 1);
    LoopExchange <- matrix(0, 10, 1);

    for(i in 1 : 10)
    {
        LHD <- randomMatrix(run, var);
        temp <-  ESEOuter(LHD);
        PhiPBest[i, 1] <- temp$PhiP;
        TimeElapsed[i, 1] <- temp$Time;
		LoopExchange[i, 1] <- temp$Loop;

        resultESE <- data.frame(PhiPBest, TimeElapsed, LoopExchange);
    }

	fileName <- paste("c:\\E-ESEdata", substr(date(), 5, 7),substr(date(), 9, 10), ".csv", sep="") 
	write.table(substr(date(), 5, 19), file = fileName, append = FALSE, row.names = FALSE, col.names = FALSE);
    write.table(resultESE, file = fileName, append = TRUE, sep = ",", row.names = FALSE, col.names = c("Phi-P Best", "Time Elapsed", "Loop Exchange"));
    return(resultESE);
}