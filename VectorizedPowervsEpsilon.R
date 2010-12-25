
PowerB<- function(D, n, B, alpha, epsilon)
{ 
x <- rnorm(n) 
y <- rnorm(n, D) 
#The data

X <- matrix(rep(x,B), n, B)
Y <-  matrix(rep(y,B), n, B) 
 #Repeat the data B times and put it in a matrix


Xprime <- apply(X, 2, sample, n, replace = TRUE)
Yprime <- apply(Y, 2, sample, n, replace = TRUE)
#Resample from the columns of the matrix
se <- apply( Xprime - Yprime, 1, sd)/sqrt(n)

Bool <- apply(abs(Xprime - Yprime), 1, '<', epsilon )
#It doesn't matter if '<' is applied to rows or columns so long as it is applied to every element of the matrix.

p <- apply(Bool, 2, sum, na.rm = TRUE)	/B
powerB <- sum(p > 1 - alpha)/n
powerB
}


D <- 0
n <- 500
alpha <- .05
epsilon <- 1
B <- 1000
k <- 100
q <- 200
	
TB <- rep(NA, k)
E <- rep(NA,k)
for(i in 0:k)
{ 
	E[i] <- i/q
	TB[i] <- PowerB(D, n, B, alpha, i/q)
	
}	
	
plot(E, TB, cex = .5, pch = 16, main = "Normal Data: Power vs. Epsilon", xlab = "Epsilon", ylab = "Power")