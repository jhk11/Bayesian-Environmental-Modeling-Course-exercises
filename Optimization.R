#Initial draft by John H. Kim INCOMPLETE!!!
#Optimization using hydromad packag using the shuffled complex evolution
#Using Les Vidourlades flood example
#May 4 2017 UFZ Leipzig Environmental modeling course
library(hydromad)

Data <- read.table("C:/Users/John/Desktop/2017_UFZ_Springschool/ratingcurve_data/Ardeche_gaugings.txt",header=T)
h <- Data$H.m.

#Define the function we want to optimize
#In our case, the likelihood function of the probability of observed Y given the parameter set
likelihood <- function(x){
x1 <- x[1]
x2 <- x[2]
x3 <- x[3]

h = Data$H.m.
Q_obs = Data$Q.m3.s.

1/sqrt(2*200^2*pi)*exp(-(Q_obs - x1*(h-x2)^x3 -0)^2/(2*200^2))
}

SCEoptim(FUN, par, ..., lower = -Inf, upper = Inf, fnscale = -1, control = list(fnscale=1))

## Rosenbrock Banana function
Rosenbrock <- function(x){
  x1 <- x[1]
  x2 <- x[2]
  100 * (x2 - x1 * x1)^2 + (1 - x1)^2
}
#lower <- c(-10,-10)
#upper <- -lower
ans <- SCEoptim(Rosenbrock, c(-1.2,1), control = list(fnscale = 1))
str(ans)
