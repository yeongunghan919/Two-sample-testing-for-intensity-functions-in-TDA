source("functions.R")
res=40 # resolution 40 by 40

load("PI_orbit5k_data/PI_orbit_r1_arrays.Rdata")
load("PI_orbit5k_data/PI_orbit_r2_arrays.Rdata")
load("PI_orbit5k_data/PI_orbit_r3_arrays.Rdata")
load("PI_orbit5k_data/PI_orbit_second_r1_arrays.Rdata")
load("PI_orbit5k_data/PI_orbit_second_r2_arrays.Rdata")
load("PI_orbit5k_data/PI_orbit_second_r3_arrays.Rdata")

orbit_r1_arrays_constant=orbit_r1_arrays[['constant']]
orbit_r1_arrays_arctan=orbit_r1_arrays[['arctan']]
orbit_r1_arrays_linear=orbit_r1_arrays[['linear']]

orbit_r2_arrays_constant=orbit_r2_arrays[['constant']]
orbit_r2_arrays_arctan=orbit_r2_arrays[['arctan']]
orbit_r2_arrays_linear=orbit_r2_arrays[['linear']]

orbit_r3_arrays_constant=orbit_r3_arrays[['constant']]
orbit_r3_arrays_arctan=orbit_r3_arrays[['arctan']]
orbit_r3_arrays_linear=orbit_r3_arrays[['linear']]

orbit_second_r1_arrays_constant=orbit_second_r1_arrays[['constant']]
orbit_second_r1_arrays_arctan=orbit_second_r1_arrays[['arctan']]
orbit_second_r1_arrays_linear=orbit_second_r1_arrays[['linear']]

orbit_second_r2_arrays_constant=orbit_second_r2_arrays[['constant']]
orbit_second_r2_arrays_arctan=orbit_second_r2_arrays[['arctan']]
orbit_second_r2_arrays_linear=orbit_second_r2_arrays[['linear']]

orbit_second_r3_arrays_constant=orbit_second_r3_arrays[['constant']]
orbit_second_r3_arrays_arctan=orbit_second_r3_arrays[['arctan']]
orbit_second_r3_arrays_linear=orbit_second_r3_arrays[['linear']]


groups <- c("r1","r2","r3")
weights <- c("constant","linear","arctan")

ttest_results <- list()

for (w in weights) {
  for (i in 1:length(groups)) {
    for (j in i:length(groups)) {  
      group_i <- groups[i]
      group_j <- groups[j]
      
      arr_i <- get(paste0("orbit_", group_i, "_arrays_", w))
      arr_j <- get(paste0("orbit_", group_j, "_arrays_", w))
      
      if (i == j) {
        arr_j <- get(paste0("orbit_second_", group_j, "_arrays_", w))
      }
      test_name <- paste0(group_i, "_vs_", group_j, "_", w)
      ttest_results[[test_name]] <- ttestpi(arr_i, arr_j, res, cc = 0.2)
    }
  }
}

ttest_results
write.csv(ttest_results, file = "../Rdata/orbit5k_PI_result.csv")
