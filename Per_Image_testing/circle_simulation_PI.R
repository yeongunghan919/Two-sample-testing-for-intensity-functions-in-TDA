source("functions.R")

# constant
load("PI_circle_data/data_pi_onecircle_constant.Rdata")
load("PI_circle_data/data_pi_twocircle_constant.Rdata")


# settings
sig=seq(0.05,0.20,length.out = 4)
npc=10
nset=100
range=c(0,2)
res=40
alpha=0.05

## Power.  
pvals_power <- matrix(NA, nrow = length(sig), ncol =nset)

test_results=ts_main_power(onepi_constant,twopi_constant,sig,npc,nset,range,res,alpha)

rowMeans(test_results) #0.968 0.694 0.236 0.098
saveRDS(rowMeans(test_results), "../results/simulation_results/Circle_results/Circles_PI_power_result_diffnoise.rds")
