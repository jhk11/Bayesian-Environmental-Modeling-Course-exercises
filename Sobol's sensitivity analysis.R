#Initial draft by John H. Kim
#Sobol's senstivity analysis of screened parameters. Main effect is calculated as a ratio of variance of the model results if one parameter is variable and variance of model results if all parameters are variable. The total effect is calculated as the ratio of variance of model output using variations in one parameter and any one other interacting parameters and variance of model results if all parameters are variable.
#May 2 2017 UFZ Leipzig Environmental modeling course

#Load two files, both with parameter sets for our model
Afiles <- c('Sobol_sequences_group_2_A_100sets.dat','Sobol_sequences_group_2_A_1000sets.dat','Sobol_sequences_group_2_A_2000sets.dat')
Bfiles <- c('Sobol_sequences_group_2_B_100sets.dat','Sobol_sequences_group_2_B_1000sets.dat','Sobol_sequences_group_2_B_2000sets.dat')
N <- c(100,1000,2000)

for (j in 1:3) {
  #A and B are parameter sets and are independent of each other. The file path is for my PC laptop but for portability, uncomment the lines below to read input files.
  A <- read.table(paste("C:/Users/John/Desktop/2017_UFZ_Springschool/Exercise_Sensitivity/samples_SI/",Afiles[j],sep=""))
  B <- read.table(paste("C:/Users/John/Desktop/2017_UFZ_Springschool/Exercise_Sensitivity/samples_SI/",Bfiles[j],sep=""))
  #A <- read.table(file.choose())
  #B <- read.table(file.choose())
  
  Data_A <- A*2*pi-pi; Data_B <- B*2*pi-pi #parameter values are normalized values and need to be converted to radians because our model has trig. functions.
  a1 <- 0.5; a2 <- 2.0
  
  P1 <- Data_A[,1]; P2 <- Data_A[,2]; P3 <- Data_A[,3]
  
  #Model output
  M_A <- sin(P1)+a1*sin(P2)*sin(P2)+a2*P3^4*sin(P1)
  
  P1_B <- Data_B[,1]; P2_B <- Data_B[,2]; P3_B <- Data_B[,3]
  M_B <- sin(P1_B)+a1*sin(P2_B)*sin(P2_B)+a2*P3_B^4*sin(P1_B)
  
  S_all <- var(c(M_A,M_B))
  Si <- c(); St <- c()
  for (i in 1:3) { #for each of the three parameters
    #change parameter of A with that of B, one parameter at a time
    Data <- A
    Data[,i] <- B[,i]
    #rescale the parameter variations to -pi to pi
    Data_rad <- Data*2*pi-pi
    P1 <- Data_rad[,1]; P2 <- Data_rad[,2]; P3 <- Data_rad[,3]
    
    #Model output
    M <- sin(P1)+a1*sin(P2)*sin(P2)+a2*P3^4*sin(P1)
    #As given in supplementary information of Cuntz and Mai et al 2015 WRR
    Si <- c(Si,sum(M_B*(M-M_A))/length(M_A))
    St <- c(St,sum((M_A-M)^2)/(2*length(M_A)))
  }
  #Main effect
  Si <- Si/S_all
  #Total effect
  St <- St/var(M_A)
  print(c("Si",N[j],signif(Si,3)))
  print(c("St",N[j],signif(St,3)))
}
#colnames for the printed outputs
print(c("Effect type","N","Parameter A","Parameter B","Parameter C"))
