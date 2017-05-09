#Initial draft by John H. Kim
#Monte Carlo Makov Chain sampling
#Sampling random parameter values based on its distribution (posterior) calcualted from the data
#PDF of prior and posterior are produced to be compared to the frequency distribution of sampled and accepted values
#Example for Floods of les Vidourlades
#May 4 2017 UFZ Leipzig Environmental modeling course

rm(list = ls(all = TRUE))
rm(list=ls())
par(mfrow=c(2,1))

Duration <- c(25.95, 25.02, 43.93, 0.26)

#Prior based on the dates is that if we assume the building to be 200 years old, 5 flood dates mean that the return period on average should be 50 years
#We will use a log-normal as a prior
mu <- log(50) #mean
sig <- 0.6 #standard deviation

#select a starting value for the parameter (return period) and specify others based on a gaussian distribution
StartingMean <- 50; #the mean has to change with every jump!
SD <- 8; #size of the jump from from the mean based on sd
n <- 10000 #number of samples to be generated


temp <- c();
for (i in 1: length(Duration)) {
  likelihood <- 1/Norm*exp(-Duration[i]/Norm) #likelihood probability distribution function
  temp <- cbind(temp, likelihood)
}
ProbDtheta <- rowProds(temp) #probability of the observed duration given our guessed mean return period (theta) of 50 years, likelihood PDF

Sampled_Final <- c(); prior_all <- c(); posterior_all <- c();
for (i in 1:(n+1)) {
  if (i==1) {
    Norm <- StartingMean #We start our simulation at a pre-defined place, our starting mean
    Sampled_Final[i] <- StartingMean #We will get rid of this pre-defined place later but will accept it for now for downstream calculations for MCMC
  } else if (Ratio >= Rand) {  #Is the probability of posterior[n+1] higher than posterior[n]? If so, accept Norm[n+1]. If not, accept of reject based on its probability. If rejected, use Norm[n]
    Sampled_Final[i] <- Norm[i]
  } else {
    Sampled_Final[i] <- Sampled_Final[i-1]
  }

  Norm[i+1] <- abs(rnorm(1,mean=Sampled_Final[i],sd=SD)) #inter-event duration selected from a gaussian distribution around mean and sd of our choice. When elements of Norm is negative, it creates NA so use absolute value.

  prior <- dlnorm(Norm[i:(i+1)], meanlog = mu, sdlog = sig) #PDF of our prior, assuming log-normal distribution
  prior_all <- c(prior_all,prior[2])

  temp <- c();
  #likelihood fuction is given by our model. Our model is that we expect inter-event duration to be exponentially distributed
  #pdf of an exponential distribution is given as
  #likelihood <- 1/theta*exp(-d/theta) theta is the mean, standard deviation and scale parameter of the distribution. In our case, it is the randomly generated, normally distributed return period vector, Norm. d is the four data points of duration that we have.

  for (j in 1: length(Duration)) {
    Exp <- 1/Norm[i:(i+1)]*exp(-Duration[j]/Norm[i:(i+1)]) #likelihood probability distribution function
    temp <- cbind(temp, Exp)
  }
  Likelihood <- rowProds(temp) #probability of the observed duration given our guessed mean return period (theta), likelihood PDF
  
  #Posterior is the product of prior and likelihood
  posterior <- prior*Likelihood
  posterior_all <- c(posterior_all,posterior[2])
  
  Ratio <- posterior[2]/posterior[1]
  Rand <- runif(1, min = 0, max = 1) #random uniform values between 0 and 1
}

if(0){
#tried to avoid using loop with ifelse function but haven't been able to make it work yet
Ratio <- posterior[2:(n+1)]/posterior[1:n]
Rand <- runif(n, min = 0, max = 1) #random uniform values between 0 and 1
Accept.Reject <- ifelse(Ratio > runif(n, min = 0, max = 1),1,0)
#Accept.Reject[which(is.na(posterior))-1] <- 0 
Sampled <- ifelse(Accept.Reject == 1,Norm[2:(n+1)],Norm[1:n]) #ifelse function cannot handle when there are consecutive 0's.
Sampled_Final <- Sampled

#  Norm <- rnorm(1,mean=Sampled_Final[i],sd=SD)
#  if(Ratio[i] >= Rand[i]) {
#    Sampled_Final[i] <- Norm[i+1]
#  } else if(i==1 & Ratio[i] < Rand[i]) {
#    Sampled_Final[i] <- Norm[1]
#  } else if (Ratio[i] < Rand[i]) {
#    Sampled_Final[i] <- Sampled_Final[i-1]
#  }

#  if(Accept.Reject[i] == 1) {
#    Sampled_Final[i] <- Norm[i+1]
#  } else if(i==1) {
#    Sampled_Final[i] <- Norm[1]
#  } else {
#    Sampled_Final[i] <- Sampled_Final[i-1]
#  }
#}
}

Norm <- Norm[2:length(Norm)]
par(mfrow=c(2,1))
plot(Norm,prior_all,xlim=c(0,80))
plot(Norm,posterior_all,xlim=c(0,80))
dev.new()
par(mfrow=c(2,1))
hist(Norm, breaks = 100,xlim=c(0,80))
hist(Sampled_Final, breaks = 100,xlim=c(0,80))
