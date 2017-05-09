#Initial draft by John H. Kim
#Bayesian inference for Floods of les Vidourlades example
#Create a posterior distribution for return period given 5 flood dates
#May 3 2017 UFZ Leipzig Environmental modeling course
rm(list = ls(all = TRUE))
rm(list=ls())
library(matrixStats)

par(mfrow=c(3,1))

#Dates <- as.Date(c("10-16-1907", 09-27-xxx))
#four inter-event duration are given from the five flood dates recorded on the side of a building
Duration <- c(25.95, 25.02, 43.93, 0.26)

#Prior based on the dates is that if we assume the building to be 200 years old, 5 flood dates mean that the return period on average should be 50 years
#We will use a log-normal as a prior
mu <- log(50) #mean
sig <- 0.6 #standard deviation
x <- c(1:200) #inter-event duration (return periods) that we want to test

prior <- dlnorm(x, meanlog = mu, sdlog = sig)
plot(x,prior)
#We can also try normal distribution (below)
#prior <- 1/(sig*sqrt(2*pi))*exp(-(log(x)-mu)^2/(2*sig^2))
#plot(x,prior)

#likelihood fuction is given by our model. Our model is that we expect inter-event duration to be exponentially distributed
#pdf of an exponential distribution is given as
#likelihood <- 1/theta*exp(-d/theta) theta is the mean, standard deviation and scale parameter of the distribution (wikipedia). In our case, it is the return period vector, x that we are interested in finding the distribution of. d is the data of duration (one of the four data points that we have).

temp <- c();
for (i in 1: length(Duration)) {
  Exp <- 1/x*exp(-Duration[i]/x) #likelihood probability distribution function
  temp <- cbind(temp, Exp)
}
ProbDtheta <- rowProds(temp) #probability of the observed duration given our guessed mean return period (theta) of 50 years, likelihood PDF
plot(x,ProbDtheta)

#Posterior is the product of prior and likelihood
posterior <- prior*ProbDtheta
plot(x,posterior)

x[which.min(abs(cumsum(posterior/sum(posterior))-0.5))] #mean of the return period according to our posterior
