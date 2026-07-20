source("functions.R")
library(reticulate)

# settings
res=40 # resolution 40 by 40
rangex=c(0,0.2)
rangey=c(0,0.2)
h=0.05
n=50
##load data 
pickle <- import("pickle")
py <- import_builtins()


f <- py$open('../data/orbit_data/small_PD.pkl', "rb")
data <- pickle$load(f)

small_PD_r1 <- data[["small_PD_r1"]]
small_PD_r2 <- data[["small_PD_r2"]]
small_PD_r3 <- data[["small_PD_r3"]]

second_small_PD_r1 <- data[["second_small_PD_r1"]]
second_small_PD_r2 <- data[["second_small_PD_r2"]]
second_small_PD_r3 <- data[["second_small_PD_r3"]]


compute_and_save_PI <- function(pd_list, group_name, res, h, rangex, rangey) {
  n <- length(pd_list)
  
  pi_constant <- vector("list", n)
  pi_linear   <- vector("list", n)
  pi_arctan   <- vector("list", n)
  
  for (jj in 1:n) {
    pd <- pd_list[[jj]]
    pdmat <- data.frame(t(matrix(pd, nrow=2, byrow=TRUE)))
    colnames(pdmat) <- c("birth","death")
    pdmat$death <- pdmat$death - pdmat$birth
    
    # weights
    arctanweight <- arctanw(pdmat, 0.5, 0.5)
    linearweight <- linearw(pdmat)
    
    # persistence image
    pi_constant[[jj]] <- pers.image(pd=pdmat, rangex=rangex, rangey=rangey,
                                    wgt=rep(1, nrow(pdmat)), nbins=res, h=h)
    pi_arctan[[jj]]   <- pers.image(pd=pdmat, rangex=rangex, rangey=rangey,
                                    wgt=arctanweight, nbins=res, h=h)
    pi_linear[[jj]]   <- pers.image(pd=pdmat, rangex=rangex, rangey=rangey,
                                    wgt=linearweight, nbins=res, h=h)
  }
  
  
  pi_constant_array <- array(0, dim=c(res, res, n))
  pi_linear_array   <- array(0, dim=c(res, res, n))
  pi_arctan_array   <- array(0, dim=c(res, res, n))
  
  for (k in 1:n) {
    pi_constant_array[,,k] <- pi_constant[[k]]
    pi_arctan_array[,,k]   <- pi_arctan[[k]]
    pi_linear_array[,,k]   <- pi_linear[[k]]
  }
  

  
  return(list(constant=pi_constant_array, arctan=pi_arctan_array, linear=pi_linear_array))
}
orbit_r1_arrays <- compute_and_save_PI(small_PD_r1, "r1", res, h, rangex, rangey)
orbit_r2_arrays <- compute_and_save_PI(small_PD_r2, "r2", res, h, rangex, rangey)
orbit_r3_arrays <- compute_and_save_PI(small_PD_r3, "r3", res, h, rangex, rangey)

orbit_second_r1_arrays <- compute_and_save_PI(second_small_PD_r1, "second_r1", res, h, rangex, rangey)
orbit_second_r2_arrays <- compute_and_save_PI(second_small_PD_r2, "second_r2", res, h, rangex, rangey)
orbit_second_r3_arrays <- compute_and_save_PI(second_small_PD_r3, "second_r3", res, h, rangex, rangey)

save(orbit_r1_arrays, file="PI_orbit5k_data/PI_orbit_r1_arrays.Rdata")
save(orbit_r2_arrays, file="PI_orbit5k_data/PI_orbit_r2_arrays.Rdata")
save(orbit_r3_arrays, file="PI_orbit5k_data/PI_orbit_r3_arrays.Rdata")

save(orbit_second_r1_arrays, file="PI_orbit5k_data/PI_orbit_second_r1_arrays.Rdata")
save(orbit_second_r2_arrays, file="PI_orbit5k_data/PI_orbit_second_r2_arrays.Rdata")
save(orbit_second_r3_arrays, file="PI_orbit5k_data/PI_orbit_second_r3_arrays.Rdata")

