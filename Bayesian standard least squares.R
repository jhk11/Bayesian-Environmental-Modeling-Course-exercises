#Initial draft by John H. Kim
#Bayesian inference on model error
#implements Bayesian least squares to find uncertainties associated with parameters in hydrologic rating curve
#example from France (gage height vs. discharge)
#May 3 2017 UFZ Leipzig Environmental modeling course

rm(list = ls(all = TRUE))
rm(list=ls())

#The file path is for my PC laptop but for portability, uncomment the lines below to read input files.
Data <- read.table("C:/Users/John/Desktop/2017_UFZ_Springschool/ratingcurve_data/Ardeche_gaugings.txt",header=T)
#Data <- read.table(file.choose())
h <- Data$H.m.

#discharge is a function of measured water height and parameters
theta1 <- 130; theta2 <- 0.5; theta3 <- 1.7
Q <- theta1*(h-theta2)^theta3
plot(h,Q,type="l")
points(h,Data$Q.m3.s.)

#create parameter values we will test
x1 <- seq(0,500,by=1); x3 <- seq(0,3,by=0.01)

#Prior based on what we know--not much so assume uniform distribution
prior1 <- dunif(x1, min = 0, max = 500, log = FALSE)
prior3 <- dunif(x3, min = 1, max = 3, log = FALSE)

#Residual is the difference of observed - estimated (error)
Residual <- Data$Q.m3.s. - Q

likelihood <- c()
for (j in 1:length(x3)){
  temp <- c()
  for (i in 1:length(h)) {
    Q1 <- x1*(h[i]-theta2)^x3[j] #using all possible values of theta1 (x1[]), simulate Q at one height value and one theta3 value (x3[])
    Residual1 <- Data$Q.m3.s.[i] - Q1
    #assume residuals are Gaussian distributed with zero mean and s.d. of 200 m3/sec
    #Normal PDF = 1/sqrt(2*sig^2*pi)*exp(-(x-mu)^2/(2*sig^2))
    #temp <- cbind(temp,1/sqrt(2*200^2*pi)*exp(-(Residual1-0)^2/(2*200^2)))
    temp <- cbind(temp,dnorm(Residual1,mean=0,sd=200,log=F))
  }
  #multiply to get likelihood PDF
  #temp2 <- rowProds(temp) #multiplying small probabilities become 0 quickly. log solves this.
  temp2 <- rowSums(log(temp))
  temp3 <- cbind.data.frame(x3=x3[j],x1=x1,PDF=temp2)
  likelihood <- rbind.data.frame(likelihood,temp3)
}
HeatCol <- rev(heat.colors(20, alpha = 1))
#arrange our probability density values in length(x1) by length(x3) matrix for graphing
z <- matrix(likelihood$PDF,nrow=length(x1),ncol=length(x3))
z <- -log(-z) #rescale for better visualization
dev.new()
image(x=x1, y=x3, z=z, col=HeatCol)

#graph Q using parameter values estimated with Bayesian least squares and compare it to the data and the Q calculated from parameter values given at the beginning of the exercise
theta1_est <- likelihood[which.max(z),]$x1
theta3_est <- likelihood[which.max(z),]$x3
Q <- theta1*(h-theta2)^theta3
dev.new()
plot(h,Q,type="l")
points(h,Data$Q.m3.s.)
Q_est <- theta1_est*(h-theta2)^theta3_est
lines(h,Q_est,lty=2)

#correlation between data and Q calculated from parameter values given at the beginning of the exercise is higher than that of data and Q using parameter values estimated with Bayesian least squares, because correlation is most likely not a good measure of fit given what we are interested in (uncertainty in parameter estimation). But standard least squares should give us the identical parameter estimate to our Bayesian method.
cor(Data$Q.m3.s.,Q_est)
cor(Data$Q.m3.s.,Q)

posterior1 <- prior1 * exp(likelihood$PDF[which(round(likelihood$x3,3) == round(theta3_est,3))])
posterior3 <- prior3 * exp(likelihood$PDF[which(round(likelihood$x1,3) == round(theta1_est,3))])

dev.new()
par(mfrow=c(2,1))
plot(x1,posterior1)
plot(x3,posterior3)

