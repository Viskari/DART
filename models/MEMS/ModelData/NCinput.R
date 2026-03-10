NCinput <- function(sites,year,N_M,N_V,VSOC){

# Location vector
  loc <- rep(0.,sites*N_V)
  
  for (jj in 1:length(sites)) {
    loc[((jj-1)*N_V+1):N_V] <- rep(sites[jj],N_V)    
  }

# Number of sites
N_S <- length(sites)

# Inflation vectors
prior_inf_mean <- rep(1.,(N_S*N_V))
prior_inf_sd <- rep(1.,(N_S*N_V))
post_inf_mean <- rep(1.,(N_S*N_V))
post_inf_sd <- rep(1.,(N_S*N_V))

# CDL file name
cdlname <- "filter_input.cdl"

write("netcdf filter_input {",cdlname)
write("dimensions:",cdlname,append=TRUE)
write(paste("      member = ", N_M," ;",sep=""),cdlname,append=TRUE)
write("      metadatalength = 32 ;",cdlname,append=TRUE)
write(paste("      location = ", N_V," ;",sep=""),cdlname,append=TRUE)
write("      time = UNLIMITED ; // (1 currently)",cdlname,append=TRUE)

write("variables:",cdlname,append=TRUE)
write("",cdlname,append=TRUE)
write("      char MemberMetadata(member, metadatalength) ;",cdlname,append=TRUE)
write('            MemberMetadata:long_name = "description of each member" ;',cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write("      double location(location) ;",cdlname,append=TRUE)
write('            location:short_name = "loc1d" ;',cdlname,append=TRUE)
write('            location:long_name = "The line in the database" ;',cdlname,append=TRUE)
write('            location:dimension = 1 ;',cdlname,append=TRUE)
write(paste("            location:valid_range = 0., ", N_S,". ;",sep=""),cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write("      double state(time, member, location) ;",cdlname,append=TRUE)
write('            state:long_name = "the ensemble of model states" ;',cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write("      double state_priorinf_mean(time, location) ;",cdlname,append=TRUE)
write('            state_priorinf_mean:long_name = "prior inflation value" ;',cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write("      double state_priorinf_sd(time, location) ;",cdlname,append=TRUE)
write('            state_priorinf_sd:long_name = "prior inflation standard deviation" ;',cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write("      double state_postinf_mean(time, location) ;",cdlname,append=TRUE)
write('            state_postinf_mean:long_name = "posterior inflation value" ;',cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write("      double state_postinf_sd(time, location) ;",cdlname,append=TRUE)
write('            state_postinf_sd:long_name = "posterior inflation standard deviation" ;',cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write("      double time(time) ;",cdlname,append=TRUE)
write('            time:long_name = "valid time of the model state" ;',cdlname,append=TRUE)
write('            time:axis = "T" ;',cdlname,append=TRUE)
write('            time:cartesian_axis = "T" ;',cdlname,append=TRUE)
write('            time:calendar = "none" ;',cdlname,append=TRUE)
write('            time:units = "year" ;',cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write("      double advance_to_time ;",cdlname,append=TRUE)
write('            advance_to_time:long_name = "desired time at end of the next model advance" ;',cdlname,append=TRUE)
write('            time:axis = "T" ;',cdlname,append=TRUE)
write('            time:cartesian_axis = "T" ;',cdlname,append=TRUE)
write('            time:calendar = "none" ;',cdlname,append=TRUE)
write('            time:units = "year" ;',cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write("      // global attributes:",cdlname,append=TRUE)
write('            :title = "an ensemble of spun-up model states" ;',cdlname,append=TRUE)
write('            :version = "$Id: $" ;',cdlname,append=TRUE)
write('            :model = "MEMS v1" ;',cdlname,append=TRUE)
write('            :history = "None" ;',cdlname,append=TRUE)

write("data:",cdlname,append=TRUE)
write("",cdlname,append=TRUE)
write(" MemberMetadata =",cdlname,append=TRUE)

for (ii in 1:(N_M-1)) {
  write(paste('  "ensemble member      ',ii,'",',sep=""),cdlname,append=TRUE)
}
write(paste('  "ensemble member      ',N_M,'" ;',sep=""),cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write(paste('location =  ',paste(loc,collapse=", "),' ;',collapse="",sep=""),cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write(" state =",cdlname,append=TRUE)

P_SOC <- rep(1., (N_M*N_V))

for (ii in 1:N_M) {
  C5 <- rnorm(1, VSOC[2], 0.1*VSOC[2])
  C8 <- rnorm(1, VSOC[3], 0.1*VSOC[3])
  C9 <- rnorm(1, VSOC[4], 0.1*VSOC[4])
  C10 <- rnorm(1, VSOC[5], 0.1*VSOC[5])
  CTot <- C5+C8+C9+C10

  P_SOC[((ii-1)*N_V+1):(ii*N_V)] <- c(CTot,C5,C8,C9,C10)  
}

write(paste(paste(P_SOC,collapse=", "),' ;',collapse="",sep=""),cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write(" state_priorinf_mean =",cdlname,append=TRUE)
write(paste("  ",paste(prior_inf_mean,collapse=", "),' ;',sep=""),cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write(" state_priorinf_sd =",cdlname,append=TRUE)
write(paste("  ",paste(prior_inf_sd,collapse=", "),' ;',sep=""),cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write(" state_postinf_mean =",cdlname,append=TRUE)
write(paste("  ",paste(post_inf_mean,collapse=", "),' ;',sep=""),cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write(" state_postinf_sd =",cdlname,append=TRUE)
write(paste("  ",paste(post_inf_sd,collapse=", "),' ;',sep=""),cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write(paste(" time = ", year," ;",sep=""),cdlname,append=TRUE)
write("",cdlname,append=TRUE)

write(paste(" advance_to_time = ", year," ;",sep=""),cdlname,append=TRUE)
write("}",cdlname,append=TRUE)

}