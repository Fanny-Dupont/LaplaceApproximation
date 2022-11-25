# ================================================
# Implementation Laplace Approximation
# 
# STAT 548
# 
# Fanny Dupont
# Department of Statistics
# University of British Columbia
# ================================================
library("madness")
library("numDeriv")
library("maxLik")
library("mvtnorm")
library("stats")
library("geostats")
library("CircStats")
library("cubature")

#numDeriv::hessian(f,x,method="Richardson"), #Uses Richardson's extrapolation as it is said in Koyama
#x must be an array
source("DataBear.R")
Get_Predictive_function<- function(mu,Sigma,dim,state_distribution,state_space){
  if (dim == 1){
    hat_p <- function(x){
      dnorm(x,mean=mu,sd=sqrt(Sigma))
    }
    f <- function(x,y){
      state_distribution(y,x)*hat_p(x)
    }
    
    hat_p_predict <- function(y){
      g <- function(x){
        f(x,y)
      }
      return(integrate(g,state_space[1],state_space[2]))
    }
  }
  else{
    if(dim == 2){
      hat_p_predict <- function(y){
        
        gg<- function(x){
          if(abs(det(Sigma))>10000){
            scal <- 0
          }else{
            scal <- 1/(2*pi*det(Sigma)^{1/2})
          }
          d <- matrix(c(x[1],x[2])-mu, ncol=2,nrow=1)
          term <-exp(-0.5*d%*%solve(Sigma)%*%t(d))
          scal*term
          #d <- matrix(c(x1,x2), ncol=2,nrow=1)
          return(scal*term*state_distribution(y[1],y[2],x[1],x[2]))
        }
        return(hcubature(gg,lowerLimit = c(state_space[,1]),upperLimit=c(state_space[,2]))$integral)
        
      }
    }else{
      if(dim == 3){
        hat_p_predict <- function(y){
          gg<- function(x){
            if(abs(det(Sigma))>10000){
              scal <- 0
            }else{
              scal <- 1/(2*pi*det(Sigma)^{1/2})
            }
            d <- matrix(c(x[1],x[2],x[3])-mu, ncol=2,nrow=1)
            term <-exp(-0.5*d%*%solve(Sigma)%*%t(d))
            scal*term
            #d <- matrix(c(x1,x2), ncol=2,nrow=1)
            return(scal*term*state_distribution(y[1],y[2],y[3],x[1],x[2],x[3]))
          }
          return(hcubature(gg,lowerLimit = c(state_space[,1]),upperLimit=c(state_space[,2]))$integral)
          
        }
      }
      
    }
  }
  
  
  
  return(hat_p_predict)
  
}


Laplace_Approximation <- function(obs,log_like1,dim,state_distribution,state_space,LikeY,A,B){
  #Initialize
  log_like <- log_like1
  T <- length(obs[,1])
  grad_log_like<-function(x){
    numDeriv::grad(log_like,x,method="Richardson")
  }
  x <- obs
  x[1,]<-obs[2,]
  hessian_log_like <- function(x){
    numDeriv::hessian(log_like,x,method="Richardson")}
  HatPPredict <- list()
  for(t in 1:T){
    
    #Compute
    MaxParam <- maxLik(log_like,grad=grad_log_like,start=x[t,], constraints=list(ineqA=A, ineqB=B))
    result_max <-list(x_estimate =c(MaxParam$estimate),Log_Like_Value=c(MaxParam$maximum))
    cov <- solve(-hessian_log_like(result_max$x_estimate))
    if(det(cov)<0){
      cat("Iteration ", t, "\n", sep = "", red("The covariance matrix is not positive-definite"))
      break()
    }
    approx <- list(mean =c(result_max$x_estimate), Sigma = cov)
    f <- Get_Predictive_function(approx$mean,approx$Sigma,dim,state_distribution,state_space)
    HatPPredict<- append(HatPPredict,f)
    
    #Incrementation 
    Predict<-HatPPredict[[t]]
    log_like <- function(xt) log(LikeY(obs[t,],xt)*Predict(xt))
    grad_log_like<-function(x){
      numDeriv::grad(log_like,x,method="Richardson")
    }
    hessian_log_like <- function(x){
      numDeriv::hessian(log_like,x,method="Richardson")}
    cat("Iteration ", t, "\n", sep = "")
    
  }
  
  return(HatPPredict)
}
#obs is a dataframe of observations. It can have multiple columns is the states are multidimensional.
#rp generates random deviantes from the initial distribution p with parameters "parameters". 
#log_like takes 3 arguments, xt,yt and y1:t-1
#State-Space is a matrix with 2 columns, where each line is the boundary for the first element of the state.
#State_distribution is a function of x_{t+1} and x_t. First argument is x_{t+1} and second x_{t}
#If dim(x_t) = 2 state_distribution has the following form:
#State_distribution(x_1{t+1},x_2{t+1},x_1{t},x_2{t})
#A and B are for if the maximization requires constraints.



#Initialization
state_distribution <- function(y1,y2,x1,x2){#We fix here kappa=1. First coordinate of the state is the angle and the
  #y1 is the current angle, x1 is the previous one.
  #y2 is the current step-length, x2 is the previous one.
  if ((y1 > pi) & (y1 <= -pi)){
    O
  }else{
    vonMises(y1, mu = x1, kappa = KMLE, degrees = FALSE)*dgamma(y2,shape=GammaMomentshape,scale=GammaMomentscale)#We used Bayes rule to obtain this formula.
  }
}
LikeY <- function(yt,xt){
  (2*pi)^{-1}*det(SigmaY)^{-1/2}*exp(-0.5*t(yt-xt)%*%solve(SigmaY)%*%(yt-xt))
}
dim<-2
SigmaY <- 2*diag(2)
state_space<- matrix(c(-2*pi,0,2*pi,20),nrow=2,ncol=2)

set.seed(0)
obs <- cbind(BearARGOS$Tangle,as.numeric(BearARGOS$distancekm))
T<-length(obs[,1])
obs <- obs[2:20,]
#Adapt the dataset to the model we chose. Some step-length in the observation are to high for the gamma distribution we picked.
#Adapting dataset

#The angles modulo 2pi
idx <-which(obs[,2]>=40)
for (m in (1:length(idx))){
  loc <- idx[m]
  a<-rgamma(1,shape=GammaMomentshape,scale=GammaMomentscale)
  obs[loc,2]<-a
}

#Initial log-likelihood, with normal density for state. Gamma did not work.
log_like1 <- function(xt){
  return(log(LikeY(obs[1,],xt)*(2*pi)^{-1}*exp(-0.5*(xt[1]-1)^2)*(2*pi)^{-1}*exp(-0.5*(xt[2]-1)^2)))
}

#Build the constraints matrices
#Here, step length >=0.001 and angles between -2*pi and 2pi, we let the angles
#with a wider space than usual since [0, 2pi] is to hard of a constraint and the MaxLik 
#function does not work well with it.
A <-  matrix(c(-1, 1,0,0,0,1), 3,2)
B <- c(2*pi,2*pi,0.001)

#We perform Laplace approximation
f <- Laplace_Approximation(obs,log_like1,2,state_distribution,state_space,LikeY,A,B)


#Get the turning angles######
N <- 20
l<- abs(state_space[1,2]-state_space[1,1])
deltax <- l/N

Ex <-array(dim=c(T))
density <- array(dim=c(T))

for(i in (1:T)){
  densityx <- function(x){
    g <- function(y){
      f[[i]](c(x,y))
    }
    return(hcubature(g, lowerLimit = state_space[2,1], upperLimit=state_space[2,2])$integral)
  }
  F <- array(dim=c(N))
  for (j in 0:(N-1)){
    F[j+1] <- densityx(0.5*(state_space[1,1]+j*deltax+state_space[1,1]+(j+1)*deltax))
    cat("",j, "\n", sep = "")
  }
  
  density[i] <- abs((state_space[1,1]-state_space[1,2]))*sum(F)
  F <- array(dim=c(N))
  
  for (j in 0:(N-1)){
    F[j+1] <- c(densityx(0.5*(state_space[1,1]+j*deltax+state_space[1,1]+(j+1)*deltax))*(0.5*(state_space[1,1]+j*deltax+state_space[1,1]+(j+1)*deltax)))/density[i]
    cat("",j, "\n", sep = "")
  }
  
  Ex[i] <- abs((state_space[1,1]-state_space[1,2]))*sum(F)
  cat("Iteration ", i, "\n", sep = "")
  cat("Expectation ", Ex[i], "\n", sep = "")
  
}


Ey <-array(dim=c(T))
ydensity <- array(dim=c(T))
for(i in (1:T)){
  densityy <- function(y){
    g <- function(x){
      f[[i]](c(x,y))
    }
    for (j in 0:(N-1)){
      F[j+1] <- c(g(0.5*(state_space[1,1]+j*deltax+state_space[1,1]+(j+1)*deltax))*(0.5*(state_space[1,1]+j*deltax+state_space[1,1]+(j+1)*deltax)))
    }
    return(abs((state_space[1,1]-state_space[1,2]))*sum(F))
  }
  
  F <- array(dim=c(N))
  for (j in 0:(N-1)){
    F[j+1] <- densityy(0.5*(state_space[2,1]+j*deltax+state_space[2,1]+(j+1)*deltax))
    cat("",j, "\n", sep = "")
  }
  
  ydensity[i] <- abs((state_space[2,1]-state_space[2,2]))*sum(F)
  
  for (j in 0:(N-1)){
    F[j+1] <- c(densityy(0.5*(state_space[2,1]+j*deltax+state_space[2,1]+(j+1)*deltax))*(0.5*(state_space[2,1]+j*deltax+state_space[2,1]+(j+1)*deltax)))/ydensity[i]
    cat("",j, "\n", sep = "")
  }
  
  Ey[i] <- abs((state_space[2,1]-state_space[2,2]))*sum(F)
  cat("Iteration ", i, "\n", sep = "")
  cat("Expectation ", Ey[i], "\n", sep = "")
  
  
}
#Analysis of the Result####



Ex <- array(c(-1.422077,-1.782781,-1.532892,-1.124888,-0.7784948,-0.6810644,-0.5389538,-0.642001,-0.6263413,-0.6708166,-0.7960095,-0.8017846,-0.7489581,-0.6543089,-0.3683804,-0.3441443,-0.3442728,-0.3924621,-0.3876097))
Ey <- array(rep(6.731358,19))

Result <- cbind(data.frame(Tangle=c(Ex),Distance=c(Ey)))
Result$MotionAngle <- array() 
Result$MotionAngle[1]<-Result$Tangle[1]
for (i in (2:19)){
  Result$MotionAngle[i] <- Result$Tangle[i]+Result$MotionAngle[i-1]%%2*pi
}
Result$Lon <- Result$Lat <- array()
Result$LonR[1] <- BearARGOS$Lon[1]* pi/180
Result$LatR[1] <- BearARGOS$Lat[1]*pi/180
Result$AD <- Result$Distance/6371
for (i in (2:19)){
  Result$LatR[i] <- asin(sin(Result$LatR[i-1])*cos(Result$AD[i-1])+cos(Result$LatR[i-1])*sin(Result$AD[i-1])*cos(Result$MotionAngle[i-1]))
  Result$LonR[i]<- Result$LonR[i-1]+ atan2(sin(Result$MotionAngle[i-1])*sin(Result$AD[i-1])*cos(Result$LatR[i-1]),cos(Result$AD[i-1])-sin(Result$LatR[i-1])*sin(Result$LatR[i]))
  }
Result$Lat<-(180/pi)*Result$LatR
Result$Lon<-(180/pi)*Result$LonR

# creating a data.frame with your lat/lon points

png("resultdddcomp.png")
df <- as.data.frame(cbind(Lon = Result$Lon,Lat =Result$Lat))
df1 <- as.data.frame(cbind(Lon = BearGPS$Lon,Lat =BearGPS$Lat))
#df2 <- as.data.frame(cbind(Lon = BearARGOS$Lon, Lat =BearARGOS$Lat))
bbox <-c(min(df$Lon,df2$Lon)+7,max(df$Lon,df2$Lon)-30,min(df$Lat,df2$Lat)+5,max(df$Lat,df2$Lat)-3)


basemap(limits = c(bbox)) + 
  #geom_spatial_path(data = df[1:20,], aes(x = Lon, y = Lat),crs = 4326, color = "blue")+
  #geom_spatial_path(data = df1[1:20,], aes(x = Lon, y = Lat),crs = 4326, color = "green")
  #geom_spatial_path(data = df2[1:20,], aes(x = Lon, y = Lat),crs = 4326, color = "orange")
  geom_spatial_point(data = df[1,], aes(x = Lon+2, y = Lat),crs=4326, color = "red")
geom_spatial_path(data = df[1,], aes(x = Lon, y = Lat),crs=4326, color = "blue")

dev.off()


delta <- (Result$Lat-BearGPS[1:19,]$Lat)^2+(Result$Lon-BearGPS[1:19,]$Lon)^2
delta <- sqrt(sum((Result$Lon-BearGPS[1:19,]$Lon)^2))

MSE <- sqrt(sum(delta))
acf(delta)
arima(delta, order = c(1, 0, 0))
 


