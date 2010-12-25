


phenoBoots <- function(LM, B, Contrasts, refpluun)
{
	#The phenoBoots function creates 
	#the empirical distribution F hat, 
	#samples from it, re-computes 
	#the estimator c^t beta hat some 
	#B number of times.
	
	#LM is the LM object produced
	#by running LMFit without 
	#specifying a contrasts.
	#LMFit is a list which contains
	#The quantities of interest 
	#produced by a standard linear regression.
	
	#B is the number of bootstraps	
	
	#Contrasts is contrast vector
	
	#Design is design matrix
	
	#refpluun is the index 
	#(in the contrasts vector) 
	#of source in the
	#same group as reference.
	
	
	f <- dim(LM$coef)[1]
	k <- dim(Contrasts)[2]

	#f is the number of features per array
	#k is the number of contrasts
	
	#The output of this function is a list of 	#k matrices, each of dimension f x B
	
	
	zero <- rep(0, f) 

	#A vector of zeros

	beta <- cbind(zero, LM$coef)

	#This holds the original coefficients
	#The column of zeros at the beginning 
	#corresponds ot the reference source.
	
	re.beta <- beta

	#this will hold each 
	#resampling of 
	#the original coefficients. 
	#The column of zeros is 
	#added so the reference
	#source is resampled in 
	#the bootstrap as well.
	
	distributions <- list(rep(0, k))

	#distributions is a list
	#which holds the empirical 
	#distributions for 
	#each given contrast.
	#(it is possible to give 
	#more than one contrasts 
	#vector at once)

	for(z in 1:k) 
	{
		#Moving through each contrast

	Contrasts0 <- c(Contrasts[refpluun,z], Contrasts)

	#Normally, the reference 
	#source would be given weight
	#zero in the contrasts 
	#vector. However, this function
	# identifies groups by their 
	#value in the contrast vector,
	#so here it is given the same contrast 
	#value as some other member of 
	#its same group, so that 
	#it is counted properly, and
	#concatenated to the beginning
	#of the contrasts vector.
	
	s <- list() 

	#s holds the indices of 
	#the coefficients chosen
	# in the resampling.

	u <- unique(Contrasts0)

	#u picks out each of the
	#different values of 
	#Contrasts0 exactly once.

	dist <- matrix(NA, nrow = f, ncol = B)

	#dist holds the empirical 
	#distribution for this contrast.
	#Each row is a empirical distribution for 
	#the feature of corresponding index.
	
	for(i in 1:B)
	
	#The bootstrap itself
	
	{ 

		for(j in 1:length(u) )
		{ 
		
			#This loop picks out 
			#groups. We resample from each group 
			#seperately, so any between-group
			#differences are preserved.
			
			who <- Contrasts0 == u[j]
					
			#who is in the jth group?
					
			s[[ j ]] <- as.vector(sample( 1:sum(who),
			sum(who), replace = TRUE) )
				
			#rolling the die to 
			#pick out the indices for 
			#the jth group
			
			re.beta[ , who] <- as.matrix( beta[ , 
			who])[ , s[[j]] ] 

			#resamples from the 
			#coefficient vector beta, 
			#according to which 
			#entries of beta belongs 
			#to which phenotype.
		}
	dist[,i] <- re.beta %*% as.matrix(Contrasts0)
	#The i-th bootstrap replicate for all f features.
	}	
	distributions[[z]] <- dist
	}	
	
	#The rows of dist correspond to features. 
	#The columns are bootstrapped contrast 
	#coefficient.
	
	names(distributions) <- colnames(Contrasts)
	distributions
}	



bootEQ <- function(epsilon, B, LM, Contrasts, refpluun)

#The BootEQ function carries out the Bootstrapped test of
#equivalence, calculating a significance level for the 
#equivalence of each feature based on the empirical 
#bootstap distribution. Significance level is defined as the probability of rejecting a true null hypothesis (type one error). Hence, the function returns the probability that the true value really is outside epsilon.

{
	f <- dim(LM$coef)[1]
	k <- dim(Contrasts)[2]

	#f is the number of features per array
	#k is the number of contrasts
	
	pboots <- phenoBoots(LM, B, Contrasts, refpluun)
		#The resampling step.
	ctbeta <- LM$coef %*% Contrasts
		#Calculating the
		#contrast coefficients.
		#ctbeta is m x k, where m
		#is the number of sources
		#and k is the number of contrasts	
	p <- rep(NA, f)
	
		# p will hold to p-values
	
	sig <- matrix(NA, f, k )
		
		#sig will hold a p-value for each
		#feature and each contrasts

	for(z in 1:k)
	{	
		boots <- pboots[[z]]
		cbeta <- as.matrix(ctbeta[,z])
		
		#boots and cbeta correspond to the z-th contrast
		
		Each.pos <- matrix( boots > epsilon, nrow = f, ncol = B )
		Each.neg <- matrix( -boots < -epsilon, nrow = f, ncol = B )
			#Each has the same dimension as boots
			#This returns TRUE if a 
			#point of the empirical
			#distribution is greater than
			#epsilon, FALSE
			#if it is less than
			#epsilon, and NA if there 
			#was an NA in that feature
			#for cbeta.
			
		
		b.pos <- apply(Each.pos, 1, sum, na.rm = TRUE)
		b.neg <- apply(Each.neg, 1, sum, na.rm = TRUE)

			#How many points in the
			#empirical distribution 
			#are outside cbeta and zero?
			
		p <- (b.neg + b.pos) / B
			#What fraction is this of total?
			
		rej <- as.vector( (abs(cbeta) >= 
			epsilon) & !is.na(cbeta) )
			
			#The rej are the contrast 
			#coefficients outside (-epsilon, epsilon)
			#and which are not NA
			
		p[rej] <- 1
			#If the contrast coefficient
			# is not inside the interval,
			# it is not equivalent
		p[is.na(cbeta)] <- NA	
				
		sig[,z] <- p
			
		}
	colnames(sig) <- paste(colnames(Contrasts),
					".tEQ", sep = "")
	#Making sure the names are consisitent					
	sig
}		