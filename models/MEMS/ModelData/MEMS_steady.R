MEMS_steady <- function(vPar,slcs){

#args <- commandArgs(trailingOnly = TRUE)
#tile<-as.numeric(args[1])+1
#sim<-as.numeric(args[2])
#fold<- tile

###set dir reading raster
dirSL<-("/eos/jeodpp/data/projects/SOIL-NACA/MODEL/SptLYR/")
dirSC<-("/eos/jeodpp/data/projects/SOIL-NACA/MODEL/statCRP/")

###read LUCAS all LU 2009
LUsub<-c("C10", "C20", "C30", "E10", "E20", "E30", "CR")


####grid meteo############################################################################
  
  EOBS01<-read.dbf(paste0(dirSL, "EOBS_0.1reg_land.dbf"))
  EOBS11<-read.dbf(paste0(dirSL, "EOBS_11rot_land.dbf"))
 
 lat_long<-as.array(stack(paste0(dirSL,"lat_long.tif")))/100	##*100
  grid2D.rcm<- cbind(EOBS01$x, EOBS01$y) 		####0.1d regular
  grid2D.rcm_<- cbind(EOBS11$x, EOBS11$y) 		####11' rotated

  ###conversion LUCAS LU to MEMS plant type
  PLT<-read.table("/eos/jeodpp/data/projects/SOIL-NACA/MEMS/mems_lucas/DB/LU_SCH.txt", header = TRUE, sep="\t")
  
  ####tiles 
  #lo<- round(seq(1,nrow(lcs), nrow(lcs)/sim),0)
  #up<- c(lo[2:(sim)]-1, nrow(lcs))
  
  out <- matrix(0.,nrow=nrow(slcs),ncol=11)


###########################################################################################

for (i in 1:nrow(slcs)) {
   #lcs<-subset(lcs, LU %in% LUsub)		###subset LU
   if (!(slcs$LU[i]  %in% c("C10", "C20", "C30", "E10", "E20", "E30", "CR"))) next	
   if (is.na(slcs$NPP[i])) next
  
  C <- rep(0,11)
  soil <- rep(1,7)
  LQ <- rep(1,4)
  env_init <- matrix(0.,nrow=365,ncol=4)

  LU <- slcs$LU[i]

  if(LU == "E10" | LU == "E20" | LU == "E30" | LU == "CR"){
    LQ[1] <- 1.1   # LitN
    LQ[2] <- 0.35   # fsi
    LQ[3] <- 0.15  # flig
    LQ[4] <- 0.35  # fdoc
  }

  if(LU == "C10"){
    LQ[1] <- 1.32   # LitN
    LQ[2] <- 0.4   # fsi
    LQ[3] <- 0.27  # flig
    LQ[4] <- 0.15  # fdoc
  }

  if(LU == "C20"){
    LQ[1] <- 0.41   # LitN
    LQ[2] <- 0.35  # fsi
    LQ[3] <- 0.32  # flig
    LQ[4] <- 0.15  # fdoc
  }
  
  if(LU == "C30"){
    LQ[1] <- 0.87   # LitN
    LQ[2] <- 0.375   # fsi
    LQ[3] <- 0.295  # flig
    LQ[4] <- 0.15  # fdoc
  }
  
   
### in_Soil ##############################################################################
 BDm<- 0.80806 + (0.823844 * exp(-0.27993 * slcs$OC09[i] * 0.1))+(0.0014065 * slcs$sand[i])-(0.0010299 * slcs$clay[i])   #BD Fine Earth 
 BDm_c<- (100- slcs$coarse[i])/((slcs$coarse[i]/2.65)+(100-slcs$coarse[i])/BDm)							 			#BD FE without gravel
  PORm<- (1-BDm/2.65)
  sk_p<- (BDm - BDm_c)/BDm																 							#volume gravel

 WC33   <-   (0.2576 +(-0.002 * slcs$sand[i])+(0.0036 * slcs$clay[i])+ 0.0299 * min(slcs$OC[i]* 0.1 * 1.72, 6))
 WC1500 <-   (0.026 + (0.005 * slcs$clay[i])+ 0.0158 * min(slcs$OC[i]* 0.1 * 1.72, 6)) 
 
 Ks<- (1930*(PORm-WC33)^(3- log(WC1500/WC33)/log(33/1500)) )/10 													##Ks mm/h
 
 WC33s<- WC33 * (1-sk_p)																							##correction for gravel
 WC1500s <- WC1500 * (1-sk_p)
 #Kss <- Ks *(1-sk_p)
 
 soil[1] <- slcs$sand[i]
 soil[2] <- slcs$clay[i]
 soil[3] <- slcs$coarse[i]
 soil[4] <- slcs$pH_in_H2O[i]
 soil[5] <- slcs$OC09[i]
 soil[6] <- BDm
 soil[7] <- 20.
  
 ### in_Soil ##############################################################################
 
 
 ### in_site ###########################################################################
 lat<- slcs$lat[i]
 long<-slcs$long[i]
 
 ### in_site ###########################################################################
 
 
 
 ### in_WTH ###########################################################################
 f<-spDistsN1(pts=grid2D.rcm, pt=as.matrix(c(long,lat)), longlat=TRUE) 
 EOB01r<-EOBS01[which.min(f),1]
 
 f_<-spDistsN1(pts=grid2D.rcm_, pt=as.matrix(c(long,lat)), longlat=TRUE) 
 EOB11r<-EOBS11[which.min(f_),3]
 
 wth<- paste("/eos/jeodpp/data/projects/SOIL-NACA/MODEL/meteo/_clima/EOBS0.1reg_day_v25/", EOB01r, ".wth", sep="") 
 if (!file.exists(wth)) next 
 
 
 obs<-read.table(wth)[10594:14245, 3:7]                 ###read 2009-2018
 
 NPPin<-PLT[PLT$LU==as.character(slcs$LU[i]),4]					###NPPin after harvest ratio

 NPPd<-dnorm(seq(1:365), 200, 50) * (slcs$NPP[i]*100)*NPPin		###NPP (g/m2)
 NPPdd<-dnorm(seq(1:366), 200, 50) * (slcs$NPP[i]*100)*NPPin
 obs[,6]<-c(NPPd, NPPd,NPPd, NPPdd, NPPd, NPPd,NPPd, NPPdd, NPPd, NPPd)
 
 env_init[,1] <- obs[1:365,3]
 env_init[,2] <- obs[1:365,4]
 env_init[,3] <- obs[1:365,5]
 env_init[,4] <- obs[1:365,6]

 init <- run_MEMS_init(vPar,env_init,soil,LQ)

 out[i,] <- init

}
 
return(out)

}