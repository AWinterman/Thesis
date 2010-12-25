n <- 100
x <- rnorm(n, sd = 1)
y <- rnorm(n, sd = 1)
t.test(x,y, var.equal = TRUE)

abs(mx-my) + c(-1,1)*Qt * se

toeTOST <- function(x,y,eps=1,alpha=.05)
{
mx <- mean(x)
my <- mean(y)
dfx <- length(x)-1
dfy <- length(y)-1
dft <- dfx+dfy
se <- sqrt( ( SSQ(x)+SSQ(y) ) / dft  ) * sqrt( (1/ (dfx+1)) + (1/(dfy+1)) )
Qt <- qt(1-alpha/2,dft)
if(eps < abs(mx-my)) FALSE
else ( eps > abs(mx-my) + Qt * se  )
}

SSQ <- function(x) sum( (x-mean(x))^2)

toeBoot <- function(x,y,eps=1,alpha=.05,B=1000)
{

if( abs(mean(x)-mean(y)) > eps ) FALSE
else
{
nx <- length(x)
ny <- length(y)	
bx <- sample(x, size = B*nx, replace = TRUE)
by <- sample(y, size = B*ny, replace = TRUE)
Mx <- matrix(bx, nrow = nx, ncol = B)
My <- matrix(by, nrow = ny, ncol = B)	
onex <- matrix(1, nrow = 1, ncol = nx)
oney <- matrix(1, nrow = 1, ncol = ny) 
Ex <- onex %*% Mx / nx
Ey <- oney %*% My / ny
Delta <- Ex - Ey
Delta <- sort(Delta)
n1 <- ceiling(B*(alpha/2))
n2 <- ceiling(B*(1-alpha/2))
(Delta[n1] > -eps) && (Delta[n2] < eps)
}
}


toePower <- function(n, d=0, Reps=100, eps=1, alpha=.05, B=10000, data = "normal")
{
	
Result <- cbind( rep(FALSE, Reps), rep(FALSE, Reps))

if(data == "normal") {
x <- matrix(rnorm(n * Reps), ncol = Reps, nrow = n)
y <- matrix(rnorm(n * Reps) + d, ncol = Reps, nrow = n) 
}

if(data == "student's_t") {
x <- matrix(rt(n,6), ncol = Reps, nrow = n)
y <- matrix(rt(n,6) + d, ncol = Reps, nrow = n) 
}


if(data == "gamma") {
x <- matrix( rgamma(n, shape = 9 , rate = 2), ncol = Reps, nrow = n)
y <- matrix( rgamma(n, shape = 9 , rate = 2) + d, ncol = Reps, nrow = n) 
}


for(i in 1:Reps)
 	{
	Result[i,1] = toeBoot(x[,i], y[,i], eps, alpha, B)
	Result[i,2] = toeTOST(x[,i], y[,i], eps, alpha)
	}	
pB <- sum(Result[,1], na.rm = TRUE)/Reps
pT <- sum(Result[,2], na.rm = TRUE)/Reps
BTpower <- cbind(pB, pT)
colnames(BTpower) <- c("Bootstrap", "TOST")
BTpower
}

date()
p0 <- cbind(numeric(40), numeric(40) )
for(n in 1:40)  p0[n,] <- toePower(n+5, d = 0, Reps = 1000, eps = 1, alpha = 0.05, B = 1000)
date()


p.25 <- cbind(numeric(40), numeric(40) )
for(n in 1:40)  p.25[n,] <- toePower(n+5, d = 0.25, Reps = 1000, eps = 1, alpha = 0.05, B = 10000)


date()
p.5 <- cbind(numeric(40), numeric(40) )
for(n in 1:40)  p.5[n,] <- toePower(n+5, d = 0.5, Reps = 1000, eps = 1, alpha = 0.05, B = 10000)


date()
p.75 <- cbind(numeric(40), numeric(40) )
for(n in 1:40)  p.75[n,] <- toePower(n+5, d = 0.75, Reps = 1000, eps = 1, alpha = 0.05, B = 10000)

date()


plot(6:45, p0[,1], type = "o",col =  "firebrick", cex = .7, pch = 1, main = "Normal Data: Power vs. Sample Size", xlab = "Sample Size", ylab = "Power", ylim = c(0,1))

lines(x = 6:45, y = p0[,2], col= "red", type = "o" , cex = .7, pch = 3)

lines(x = 6:45, y = p.25[,1], col = "purple", type ="o", cex = .7, pch = 1)

lines(x = 6:45, y = p.25[,2], col = "violet", type ="o", cex = .7, pch = 3)

lines(x = 6:45, y = p.5[,1], col = "darkgreen", type ="o", cex = .7, pch = 1)

lines(x = 6:45, y = p.5[,2], col = "green", type ="o", cex = .7, pch = 3)

lines(x = 6:45, y = p.75[,1], col = "black", type ="o", cex = .7, pch = 1)

lines(x = 6:45, y = p.75[,2], col = "dimgrey", type ="o", cex = .65, pch = 3)


legend( x = 6, y = 1, legend = c( c(expression( paste("Bootstrap, ", Delta, " = 0") )), c(expression( paste("t test, ", Delta, " = 0") )), c(expression( paste("Bootstrap, ", Delta, " = 0.25") )), c(expression( paste("t test, ", Delta, " = 0.25") )), 
c(expression( paste("Bootstrap, ", Delta, " = 0.5") )), c(expression( paste("t test, ", Delta, " = 0.5") )), c(expression( paste("Bootstrap, ", Delta, " = 0.75") )), c(expression( paste("t test, ", Delta, " = 0.75") ))  ), col = c("firebrick", "red", "purple", "violet", "darkgreen", "green","black","dimgrey"), pch = c(1, 3), cex = .65  )

getwd()
 dev.copy(pdf, "NormalPowerEpsilon.pdf")
 dev.off()


save.image(file = "PowerEpsilon.RData")


q0 <- cbind(numeric(40),numeric(40) )
for(n in 1:40)  q0[n,] <- toePower(n+5, d = 0, Reps = 1000, eps = 1, alpha = 0.05, B = 1000, data = "student's_t")


q.25 <- cbind(numeric(40), numeric(40) )
for(n in 1:40)  q.25[n,] <- toePowerT(n+5, d = 0.25, Reps = 1000, eps = 1, alpha = 0.05, B = 1000, data = "student's_t")

q.5 <- cbind(numeric(40), numeric(40) )
for(n in 1:40)  q.5[n,] <- toePowerT(n+5, d = 0.5, Reps = 1000, eps = 1, alpha = 0.05, B = 1000, data = "student's_t")

q.75 <- cbind(numeric(40), numeric(40) )
for(n in 1:40)  q.75[n,] <- toePowerT(n+5, d = 0.75, Reps = 1000, eps = 1, alpha = 0.05, B = 1000, data = "student's_t")


plot(6:45, q0[,1], type = "o",col =  "firebrick", cex = .7, pch = 1, main = "Student's t data: Power vs. Sample Size", xlab = "Sample Size", ylab = "Power", ylim = c(0,1))

lines(x = 6:45, y = q0[,2], col= "red", type = "o" , cex = .7, pch = 3)

lines(x = 6:45, y = q.25[,1], col = "purple", type ="o", cex = .7, pch = 1)

lines(x = 6:45, y = q.25[,2], col = "violet", type ="o", cex = .7, pch = 3)

lines(x = 6:45, y = q.5[,1], col = "darkgreen", type ="o", cex = .7, pch = 1)

lines(x = 6:45, y = q.5[,2], col = "green", type ="o", cex = .7, pch = 3)

lines(x = 6:45, y = q.75[,1], col = "black", type ="o", cex = .7, pch = 1)

lines(x = 6:45, y = q.75[,2], col = "dimgrey", type ="o", cex = .65, pch = 3)


legend( x = 6, y = 1, legend = c( c(expression( paste("Bootstrap, ", Delta, " = 0") )), c(expression( paste("t test, ", Delta, " = 0") )), c(expression( paste("Bootstrap, ", Delta, " = 0.25") )), c(expression( paste("t test, ", Delta, " = 0.25") )), 
c(expression( paste("Bootstrap, ", Delta, " = 0.5") )), c(expression( paste("t test, ", Delta, " = 0.5") )), c(expression( paste("Bootstrap, ", Delta, " = 0.75") )), c(expression( paste("t test, ", Delta, " = 0.75") ))  ), col = c("firebrick", "red", "purple", "violet", "darkgreen", "green","black","dimgrey"), pch = c(1, 3), cex = .65  )


getwd()
 dev.copy(pdf, "Student'sTPowerEpsilon.pdf")
 dev.off()
 
 
save.image(file = "PowerEpsilon.RData")




toePowerG <- function(n, d=0, Reps=10, eps=1, alpha=.05, B=10000)
{
Result <- cbind( rep(FALSE, Reps), rep(FALSE, Reps))

for(i in 1:Reps)
 {
 	x <- rgamma(n, shape = 9 , rate = 2)
	y <- rgamma(n, shape = 9 , rate = 2) + d
	Result[i,1] = toeBoot(x, y, eps, alpha, B)
	Result[i,2] = toeTOST(x,y,eps,alpha)
}	
	pB <- sum(Result[,1])/Reps
	pT <- sum(Result[,2])/Reps
	BTpower <- cbind(pB, pT)
	colnames(BTpower) <- c("Bootstrap", "TOST")
	BTpower
}

date()
g0 <- cbind(numeric(40),numeric(40) )
for(n in 1:40)  g0[n,] <- toePowerG(n+5, d = 0, Reps = 1000, eps = 1, alpha = 0.05, B = 10000)
date()
g.25 <- cbind(numeric(40),numeric(40) )
for(n in 1:40)  g.25[n,] <- toePowerG(n+5, d = 0.25, Reps = 1000, eps = 1, alpha = 0.05, B = 10000)
date()
g.5 <- cbind(numeric(40),numeric(40) )
for(n in 1:40)  g.5[n,] <- toePowerG(n+5, d = 0.5, Reps = 1000, eps = 1, alpha = 0.05, B = 10000)
date()
g.75 <- cbind(numeric(40),numeric(40) )
for(n in 1:40)  g.75[n,] <- toePowerG(n+5, d = 0.75, Reps = 1000, eps = 1, alpha = 0.05, B = 10000)
date()

plot(6:45, g0[,1], type = "o",col =  "firebrick", cex = .7, pch = 1, main = "Gamma data: Power vs. Sample Size", xlab = "Sample Size", ylab = "Power")

lines(x = 6:45, y = g0[,2], col= "red", type = "o" , cex = .7, pch = 3)

lines(x = 6:45, y = g.25[,1], col = "purple", type ="o", cex = .7, pch = 1)

lines(x = 6:45, y = g.25[,2], col = "violet", type ="o", cex = .7, pch = 3)

lines(x = 6:45, y = g.5[,1], col = "darkgreen", type ="o", cex = .7, pch = 1)

lines(x = 6:45, y = g.5[,2], col = "green", type ="o", cex = .7, pch = 3)

lines(x = 6:45, y = q.75[,1], col = "black", type ="o", cex = .7, pch = 1)

lines(x = 6:45, y = g.75[,2], col = "dimgrey", type ="o", cex = .65, pch = 3)

getwd()
 dev.copy(pdf, "GammaPowerEpsilon.pdf")
 dev.off()
save.image(file = "PowerEpsilon.RData")
date()


Ep0 = p0[,1] - p0[,2]
sum(Ep0[1:20]/20)

Ep.25 = p.25[,1] - p.25[,2]
sum(Ep.25[1:20]/20)

Ep.5 = p.5[,1] - p.5[,2]
sum(Ep.5[1:20]/20)


Epb0 = p0[,1] - p0[,2]
sum(Ep0[21:40]/20)

Epb.25 = p.25[,1] - p.25[,2]
sum(Ep.25[21:40]/20)

Epb.5 = p.5[,1] - p.5[,2]
sum(Ep.5[21:40]/20)






Eq0 = q0[,1] - q0[,2]
sum(Eq0[1:20]/20)

Eq.25 = q.25[,1] - q.25[,2]
sum(Eq.25[1:20]/20)

Eq.5 = q.5[,1] - q.5[,2]
sum(Eq.5[1:20]/20)


Eqb0 = q0[,1] - q0[,2]
sum(Eqb0[21:40]/20)

Eqb.25 = q.25[,1] - q.25[,2]
sum(Eqb.25[21:40]/20)

Eqb.5 = q.5[,1] - q.5[,2]
sum(Eqb.5[21:40]/20)