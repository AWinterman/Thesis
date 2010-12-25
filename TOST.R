
tEQ <- function(Fit, eps )
{
	
#tEQ is a function of Fit and eps. Fit 
#is the output of LMFit or EBayes, two 
#modified LIMMA  functions, (code 
#available in the  appendix) which 
#first run linear  regression on 
#Microarray data and then  carry out 
#"empirical Bayes shrinkage  of the 
#standard errors towards a common 
#value" (LIMMA).

#eps is the epsilon from the discussion  
#above: the distance from zero to the  
#endpoint of the equivalence intervals. 

#eps can be specified once for every  
#feature on the array, or once for each  
#feature on the array. i.e. it must  
#either be a single number or a vector  
#with length equal to the number of  
#features on the array.

#If more than one contrast is specified  
#for LMFit (which supports computing  
#several different contrasts at once),  
#this funtion will calculate  p-values  
#for each of them.

	eps <- abs(eps)
	d <- as.matrix(Fit$coeff) 
	
#d holds the contrast 
#coefficients, c^t Beta 

if(is.null(Fit$s2.post))
#Checks to see if Empirical Bayes
#analysis was done
	{SE <- as.matrix(Fit$sigma*Fit$stdev.unscaled)}
else
	{SE <- as.matrix(sqrt(Fit$s2.post)*Fit$stdev.unscaled)}
	
#SE holds the standard errors. If #empirical bayes analysis was done
#uses a 'better'  
#estimator for standard deviation  
#than the one described  
#in the appendix

	dgf <- Fit$df.resid
	
#dgf is the degrees of freedom of  
#the residuals

	if( length(eps) != 1 & length(eps) != nrow(d) )
	{
		
#Stops the computation if eps  
#is not specified per-feature   
#or universally.

		stop("eps must either have length one or nrow(d)")
	}			
	p <- matrix(NA, nrow = nrow(d), ncol = ncol(d))
	
#p will hold significance levels   
#for the TOST test of equivalence 

	for ( i in 1:ncol(d) ) 
	{
		
#Moving through each
#contrast (colums of d) 

		SEc <- SE[,i]
		dc <- d[,i]
		
#dc and SEc are the ith column 
#of d and SE respectively

		j <-  (abs(dc) < eps) & !is.na(SEc) & (dgf != 0)
		
#which contrast coefficients - 
#i.e. measurements of c^t beta 
#hat - are inside the 
#equivalence region, not 
#missing due to data deletion, 
#and with nonzero degrees of freedom

		h <-  (abs(dc) >= eps)  & !is.na(SEc) & (dgf != 0)
		
#which contrast coefficients	
#are outside equivalence 	
#region, not missing due to 	
#data deletion, and with 	
#nonzero degrees of freedom
		
		####
		
		p[h,i] <- 1
		
#Contrast coefficients outside 
#the equivalence region are 
#not rejected. Hence the 
#probability of a false 
#rejection is zero, and the 
#significance level is 1.

if (length(eps)>1)

#This if() else statement 
#checks the length of eps, and 
#then computes the distance 
#(in standard error units) 
#between the contrast 
#coefficient, dc (in the i-th 
#contrast of course), and the 
#epsilon for that feature.

			{
			k <- ( eps[j] - abs( dc[j] )) / SEc[j]
			} 
else		
			{
			k <- ( eps - abs( dc[j] )) / SEc[j]

			}
			
#k is the distance between the edge of 
#the equivalence region (in the first 
#case, the region is specified per 
#feature, in the second, for every 
#feature), and the contrast coefficient 
#in units of standard error of the 
#contrast coefficient.
			p[j,i] <- pt(k, dgf[j], lower.tail = FALSE) 
			
#Finding  the probability that k < t, 
#where t is student's t distributed 
#with the appropriate number of degrees 
#of freedom.

		}
colnames(p) <- paste(colnames(Fit$contrasts),".tEQ", sep 
= "")
#Names of contrasts, so we know 
#which is which if more than 
#one was specified.
p

#This is the significance level with 
#which each feature is equivalent. 
#Specify a p-value, that is a cutoff 
#value above which we say "this feature 
#is not significantly equivalent" and 
#below which we say "this feature is 
#significantly equivalent."
}