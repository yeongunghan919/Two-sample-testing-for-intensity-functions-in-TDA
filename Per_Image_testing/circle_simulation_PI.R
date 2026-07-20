source("functions.R")

# settings
sig=seq(0.05,0.20,length.out = 4)
npc=50
nset=100
range=c(0,2)
res=40
alpha=0.05


# constant
load("PI_circle_data/data_pi_onecircle_constant.Rdata")
load("PI_circle_data/data_pi_twocircle_constant.Rdata")


# settings with different sample sizes
sig =c(0.20)
list_npc=c(20,40,70,100)
nset=50
range=c(0,4)
res=40
alpha=0.05

## Power.  
pvals_power <- matrix(NA, nrow = length(list_npc), ncol =nset)
for (n in 1:length(list_npc)){
  npc=list_npc[n]
  test_results=ts_main_power(onepi_constant,twopi_constant,sig,npc,nset,range,res,alpha)
  pvals_power[n,]=test_results[,,5]  
}

rowMeans(pvals_root_power)
saveRDS(rowMeans(pvals_root_power), "../results/simulation_results/Circle_results/Circles_PI_power_result_diffsamples.rds")
