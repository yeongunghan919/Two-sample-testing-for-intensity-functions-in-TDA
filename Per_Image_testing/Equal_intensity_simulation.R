setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("functions.R")

# settings
sig=c(0.0)
npc=50
nset=100
range=c(0,4)
res=40
alpha=0.05

load("PI_Equal_intensity_data/P_pi_torus_linear.Rdata")
load("PI_Equal_intensity_data/Q_pi_torus_linear.Rdata")
length(P_pi_linear)
length(Q_pi_linear)
P_pi_linear = list(P_pi_linear)
Q_pi_linear = list(Q_pi_linear)
pair_blocks <- as.matrix(
  read.csv("../simulations/Equal_intensity_simulation/selected_pairs_Equal_intensity_simulation.csv", header = FALSE)
)
## Power.  
pvals_linear_power=ts_main_power_using_pair_block(P_pi_linear,Q_pi_linear,sig,npc,nset,range,res,alpha,pair_blocks)
mean(pvals_linear_power[,,2]) #C=0.2 <-> ,,2 
saveRDS(mean(pvals_linear_power[,,2]), "../results/simulation_results/Equal_intensity_results/Equal_intensity_PI_power_result.rds")
