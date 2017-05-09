#Initial draft by John H. Kim
#Sensitivity analysis using Morris's Elementary Effects (for acreening parameters)
#May 2 2017 UFZ Leipzig Environmental modeling course
rm(list = ls(all = TRUE))
rm(list=ls())

#Datafiles are names of files containing the Morris trajectories for calculating the elementary effects. Trajfiles contain matrix of same dimension that specifies which of the three parameters are being changed in the Datafiles
Datafiles <- c('Morris_trajectories_group_2_10traj.dat','Morris_trajectories_group_2_50traj.dat','Morris_trajectories_group_2_100traj.dat')
Trajfiles <- c('Morris_parachanged_group_2_10traj.dat','Morris_parachanged_group_2_50traj.dat','Morris_parachanged_group_2_100traj.dat')
N <- c(10,50,100); EE_summary <- c()
for (j in 1:3) {
  #The file path is for my PC laptop but for portability, uncomment the lines below to read input files.
  Data <- read.table(paste("C:/Users/John/Desktop/2017_UFZ_Springschool/Exercise_Sensitivity/samples_EE/",Datafiles[j],sep=""))
  #Data <- read.table(file.choose())

              
  #rescale the parameter variations to -pi to pi. Parameter values are normalized values and need to be converted to radians because our model has trig. functions.
  Data_rad <- Data*2*pi-pi
  P1 <- Data_rad[,1]; P2 <- Data_rad[,2]; P3 <- Data_rad[,3]
  a1 <- 0.5; a2 <- 2.0
                                                                                                                          
  #Model output
  M <- sin(P1)+a1*sin(P2)*sin(P2)+a2*P3^4*sin(P1)
  
  #calculate the change in the model output along the trajectory (Elementary Effect)
  M_Diff <- abs(diff(M))
  
  #file of which parameters were changed in the Morris trajectory changes
  #The file path is for my PC laptop but for portability, uncomment the lines below to read input files.
  Traj_change <- read.table(paste("C:/Users/John/Desktop/2017_UFZ_Springschool/Exercise_Sensitivity/samples_EE/",Trajfiles[j],sep=""))
  #Traj_change <- read.table(file.choose())
  Traj_change <- Traj_change[,1]
  
  temp <- c()
  #Can we do this without loops? Create two vectors with indices [1:39] and [2:40] and get the diff()
  
  for (i in 1:length(Traj_change)) {
    if (i%%4 == 0) {
    temp <- c(temp,NA)
    } else {
    temp <- c(temp,abs(Data[i,Traj_change[i]]-Data[i+1,Traj_change[i]])) #how much did the particular parameter change by?
    }
  }
  
  EE <- M_Diff/temp[1:(length(temp)-1)] #the difference in model output based on the Morris trajectory is divided by the normalized change in the parameter to give elementary effect.
  EE <- cbind(N=N[j],EE,Traj_change=Traj_change[1:(length(Traj_change)-1)])
  EE_summary <- rbind(EE_summary,aggregate(EE~Traj_change+N,EE,mean))
}
