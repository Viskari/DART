library(raster)
library(foreign)
library(stringr)

source("MEMS_steady.R")
source("MEMS_R_annual.R")
source("NCinput.R")

default_par <- read.csv(
  file = "defaults_and_priors_MEMS.csv",
  colClasses = c("character", rep("numeric", 3))
)

# Line numbers indicating which time series is used from the LU filtered table
lp <- c(1)

# The starting year
year <- 2009

# Number of ensemble members
N_M <- 50

# MEMS SOC state vector length
# Here 5: SOC_tot, C5, C8, C9, C10
N_V <- 5


par_chosen <- c(1:2,23:25,31)
l_chosen <- length(par_chosen)

bPar <- default_par$value
tPar <- c(7.16,0.54,0.00053,0.000055,0.00021,0.00043) # 4DEnVar 0.35

Par <- bPar

for (jj in 1:l_chosen){
  Par[par_chosen[jj]] <- tPar[jj]
}

LUf <- c("C10", "C20", "C30", "E10", "E20", "E30", "CR")

lc <- read.csv("LUCAS.csv")

flc <- lc[lc$LU %in% LUf,]

slc <-flc[lp,]

ic <- MEMS_steady(Par,slc)

V_SOC <- c(sum(ic[5]+ic[8]+ic[9]+ic[10])/100.,ic[5],ic[8],ic[9],ic[10])

NCinput(lp,year,N_M,N_V,V_SOC)
