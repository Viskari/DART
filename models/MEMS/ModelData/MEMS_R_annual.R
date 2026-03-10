run_MEMS_init <- function(inPara,env_init,soil,LQ){
  
  Nyear <- 700
  time <- 1

  C <- rep(0.,11)
  dC <- rep(0.,11)  
  litter <- rep(0,4)
  
  C[1] <- 0.35
  C[2] <- 0.35
  C[3] <- 0.15
  C[6] <- 0.15
  C[9] <- 3000.
  
  hemis <- -1.5
  
  mean_temp <- sum((env_init[,1]+env_init[,2])/2)/365.
  
  LitN <- LQ[1]
  fsi <- LQ[2]
  flig <- LQ[3]
  fdoc <- LQ[4]
  
  jj <- 0
  LCI <- 0.

  tmod_opt <- 45.
  tmod_lag <- 4.
  tmod_shp <- 15.
  tmod_Q10 <- 2.
  tmod_Tref <- 13.5

  lstep <- Nyear*365
    
  for(ii in 1:lstep){
    #for(ii in 1:Nyear){
    jj <- jj+1
    
    if(jj > 365){
      jj <- 1
    }

    C[7] <- 0.
    
    NPP <- env_init[jj,4]
    litter[1] = NPP*fsi - NPP*fsi*fdoc
    litter[2] = NPP - NPP*(fsi+flig)
    litter[3] = NPP*flig
    litter[4] = NPP*fsi*fdoc
    
    tmp_range <- (env_init[jj,1] - env_init[jj,2])/2.
    tmp <- tmp_range*sin((2*jj/365.*pi)+hemis) + mean_temp
        
    tmp_mod = exp(-(tmp/(tmod_opt+tmod_lag))**tmod_shp)*(tmod_Q10**((tmp-tmod_Tref)/tmod_Tref))
        
    dC <- MEMS_func(inPara,C,litter,LitN,tmp_mod,soil)
    C = C + dC*time

    C[7] <- 0.
    C[11] <- 0.
        
    C[C < 0.] <- 0.
  }
  
  #  print(dC)

  return(C)
    
}


run_MEMS_ann <- function(N,inPara,init,env_in,soil,LQ){
  
  Nyear <- N
  time <- 1
  
  C <- init
  dC <- rep(0.,11)  
  litter <- rep(0,4)

  Cann <- rep(0,N)

  Cdaily <- matrix(0,nrow=5,ncol=N*365)
      
  hemis <- -1.5
  
  LitN <- LQ[1]
  fsi <- LQ[2]
  flig <- LQ[3]
  fdoc <- LQ[4]
  
  jj <- 0
  LCI <- 0.
  
  tmod_opt <- 45.
  tmod_lag <- 4.
  tmod_shp <- 15.
  tmod_Q10 <- 2.
  tmod_Tref <- 13.5
  
  lstep <- 365
  ll <- 0.
    
  for(kk in 1:N){

  mean_temp <- 0.
  for(jj in ((kk-1)*365+1):(kk*365)){
    mean_temp <- mean_temp + (env_in[jj,1]+env_in[jj,2])/2
  }

  mean_temp <- mean_temp/365.
    
  for(ii in 1:lstep){
    
    ll <- ll + 1
    
    NPP <- env_in[(kk-1)*365+ii,4]
    litter[1] = NPP*fsi - NPP*fsi*fdoc
    litter[2] = NPP - NPP*(fsi+flig)
    litter[3] = NPP*flig
    litter[4] = NPP*fsi*fdoc
        
    tmp_range = (env_in[(kk-1)*365+ii,1] - env_in[(kk-1)*365+ii,2])/2.
    tmp = tmp_range*sin((2*ii/365.*pi)+hemis) + mean_temp
        
    tmp_mod = exp(-(tmp/(tmod_opt+tmod_lag))**tmod_shp)*(tmod_Q10**((tmp-tmod_Tref)/tmod_Tref))
    
    dC <- MEMS_func(inPara,C,litter,LitN,tmp_mod,soil)
    C = C + dC*time

    Cdaily[1,ll] <- NPP
    Cdaily[2,ll] <- C[5] + C[8] + C[9] + C[10]
    Cdaily[3,ll] <- C[5] + C[8] + C[10]
    Cdaily[4,ll] <- C[9]
    Cdaily[5,ll] <- C[7]
        
    C[7] <- 0.
    C[11] <- 0.
            
    C[C < 0.] <- 0.
  }
  
  Cann[kk] <- sum(C)
  }
  
  #  print(dC)

  Cout <- list("Cann"=Cann,"Cdaily"=Cdaily)
    
  return(Cout)
  
}


MEMS_func <- function(inPara,C,litter,LitN,tmp_mod,soil){
  
  Para <- rep(0,25)
  
  Para[1] <- inPara[9]
  Para[2] <- inPara[10] 
  Para[3] <- 0.         
  Para[4] <- inPara[11]
  Para[5] <- inPara[31]
  Para[6] <- inPara[23] 
  Para[7] <- inPara[24]
  Para[8] <- inPara[25]
  Para[9] <- inPara[6]
  Para[10] <- inPara[7] 
  Para[11] <- inPara[8] 
  Para[12] <- inPara[29] 
  Para[13] <- inPara[28] 
  Para[14] <- inPara[12]
  Para[15] <- inPara[13] 
  Para[16] <- inPara[14] 
  Para[17] <- inPara[15]
  Para[18] <- inPara[18]
  Para[19] <- inPara[19]
  Para[20] <- inPara[20]
  Para[21] <- inPara[16]
  Para[22] <- inPara[17]
  Para[23] <- inPara[1]
  Para[24] <- inPara[2]
  Para[25] <- inPara[22]
  
  POMsplit <- 0.3
  LCI <- 0.
  
  if(C[3] > 0.){
    LCI <- C[3]/(C[2]+C[3])
  }
  
  gamma <- 1/(1+exp(-Para[9]*(LitN-Para[10])))
  
  if(LCI < 0.7){
    Para[3] <- Para[2]*(0.2/(1+(200/exp(8.15*LCI))))
    ek <- exp(-3*LCI)  
    eb <- 1-exp(-0.7*(0.7-LCI)*10)
  } else {
    Para[3] <- Para[2]*exp(-2.1)
    ek <- exp(-3*0.7)
    eb <- 0.
  }

  uk <- min(gamma,ek)
  ub <- min(gamma,eb)
  
  la1 = min((Para[15]-((Para[15]-Para[14])/Para[11])*LCI),(Para[15]-((Para[15]-Para[14])/Para[9])*LitN))
  la4 = min((Para[17]-((Para[17]-Para[16])/Para[11])*LCI),(Para[17]-((Para[17]-Para[16])/Para[9])*LitN))
  
  if(la1 < 0){
    la1 <- 0.
  } 
  
  if(la4 < 0){
    la4 <- 0.
  } 
  
  Klm <- exp(-0.186*soil[4]-0.216) 
  ScConc <- Para[24]*(100-soil[1]) + Para[23] 
  Qmax <- soil[7]*100*100*soil[6]/1000*ScConc*(1-soil[3]) 
  sorption <- C[8]*((Klm*Qmax*C[8])/(1+Klm*C[8])-C[9])/Qmax
    
  dC <- rep(0.,11)
  
  dC[1] <- litter[1] - Para[1]*uk*tmp_mod*C[1]
  dC[2] <- litter[2] - Para[2]*uk*tmp_mod*C[2] - Para[12]*C[2]
  dC[3] <- litter[3] - Para[3]*tmp_mod*C[3] - Para[12]*C[3]
  dC[4] <- ub*Para[19]*(1-la1)*Para[2]*uk*tmp_mod*C[2] + ub*Para[18]*(1-la4)*Para[1]*uk*tmp_mod*C[1] - Para[4]*tmp_mod*C[4]
  dC[5] <- (Para[12]*C[2] + Para[12]*C[3])*POMsplit + Para[20]*(1-Para[21])*Para[4]*tmp_mod*C[4] - Para[5]*tmp_mod*C[5]
  dC[6] <- litter[4] + la1*Para[2]*uk*tmp_mod*C[2]+ Para[22]*Para[3]*tmp_mod*C[3]+ Para[21]*Para[4]*tmp_mod*C[4]+ la4*Para[1]*uk*tmp_mod*C[1]- Para[13]*C[6]
  dC[7] <- ((1-ub*Para[18])*(1-la4))*Para[1]*uk*tmp_mod*C[1]+ ((1-ub*Para[19])*(1-la1))*Para[2]*uk*tmp_mod*C[2]+ (1-Para[22])*Para[3]*tmp_mod*C[3]+ ((1-Para[20])*(1-Para[21]))*Para[4]*tmp_mod*C[4]+ (1-Para[22])*Para[5]*tmp_mod*C[5]+ Para[6]*tmp_mod*C[8] + Para[7]*tmp_mod*C[9]+ (1-Para[22])*Para[8]*tmp_mod*C[10]
  dC[8] <- C[6]*Para[13] + Para[22]*Para[5]*C[5]*tmp_mod + Para[22]*Para[8]*C[10]*tmp_mod - sorption - Para[25]*C[8] - Para[6]*tmp_mod*C[8]
  dC[9] <- sorption - Para[7]*tmp_mod*C[9]
  dC[10] <- (1-POMsplit)*Para[12]*C[2] + (1-POMsplit)*Para[12]*C[3]- Para[8]*tmp_mod*C[10]
  dC[11] <- Para[25]*C[8]

  #print(c(Para[22]*Para[8]*C[10]*tmp_mod, sorption))
    
  return(dC)
    
}