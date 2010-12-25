#normal diagram

plot(function(x) dnorm(x,0,1), -4,4 , bty = "l")
abline(h = 0, lty = 1, col = "red")
segments(x0 = -1, y0 = -0.005, y1 = 0.005, lty = 2)
segments(x0 = 1, y0 = -0.005, y1 = 0.005, lty = 2)
segments(x0 = 0, y0 = -0.005, y1 = 0.05, lty = 2)
segments(x0 = 0, y0 = 0.075, y1 = 0.6, lty = 2)
segments(x0 = -1, y0 = 0.025, y1 = 0.6, lty = 2)
segments(x0 = 1, y0 = 0.025, y1 = 0.6, lty = 2)
text(x = -1.08, y = .015, c(expression( paste("-", epsilon) )) , cex = 2)
text(x = 1.0, y = .015, expression(epsilon), cex = 2)
text(x = 0, y = 0.062, expression(Delta), cex = 2 )

#This plot gives an example of a possible measurement of D and the probability maximizing choice of Delta allowable under the null of the two one sided t test.

plot(function(x) dnorm(x,1,1), -4,4, bty = "l", main = "Distribution of D and Equivalence region", ylab = "Probability of x" )
abline(h = 0, lty = 1, col = "red")
segments(x0 = -1, y0 = -1, y1 = 0.005, lty = 2)
segments(x0 = 1, y0 = -1, y1 = 0.005, lty = 2)
text(x = -1.08, y = 0.015, c(expression( paste("-", epsilon) )) , cex = 1.5)
text(x = 1, y = 0.015, c(expression( paste(Delta, " = ", epsilon) )), cex = 1.5 )
text(x = 0.3, y = 0.03, "D", cex = 1.5)
segments(x0 = 0.3, y0 = -0.014, y1 = 0.015, lty = 1)

segments(x0 = 1, y0 = 0.025, y1 = .6, lty = 2)
segments(x0 = -1, y0 = 0.025, y1 = 0.6, lty = 2)
abline(v = 0, lty = 2)

normal <- function(x) dnorm(x,1,2)
curve(normal, add = TRUE, lty = 3, col = "darkblue")
legend(x = 2.5, y = 0.35, col = c("black", "darkblue"),lty = c(1,3), c( "N(1,1)", "N(1,2)" ))


setwd("/Users/andrewwinterman/documents/Fall 2010/Test of Equivalence/Paper/Thesis")

dev.copy(pdf, "DiagramNormal.pdf")
dev.off()

