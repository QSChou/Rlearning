
## binominal test by tossing the coin for 5 times and get (H,H,H,T,T)
## the likelihood function L(θ) = θ^3*θ^2
curve(x^3*(1-x)^2,from = 0, to= 1,xlab = "Parameter",ylab = "Likelihood")
library("ggplot2")
eq1 = function(x){x^3*(1-x)^2}
eq2 = function(x){6*x^2*(1-x)}
eq3 = function(x){3*log(x)+2*log(1-x)}

ggplot(data.frame(x=c(0,1)), aes(x=x)) +
  stat_function(fun=eq1, geom="line") +
  stat_function(fun=eq2, geom="line") +
  stat_function(fun=eq3, geom="line") +
  xlab("x") + ylab("y")
