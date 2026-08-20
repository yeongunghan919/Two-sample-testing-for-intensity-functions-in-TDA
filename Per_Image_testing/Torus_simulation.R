setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("functions.R")

# settings
sig=c(0.01,0.02,0.03,0.04,0.05)
npc=50
nset=100
range=c(0,4)
res=40
alpha=0.05

load("PI_Torus_data/data_pi_torus_root.Rdata")
pair_blocks <- as.matrix(
  read.csv("../simulations/Torus/selected_pairs_torus_diffnoise.csv", header = FALSE)
)
## Power.  
torus_pi_root_noised = torus_pi_root[1:5]
torus_pi_root_not_noised = rep(torus_pi_root[6],5)
pvals_root_power=ts_main_power_using_pair_block(torus_pi_root_noised,torus_pi_root_not_noised,sig,npc,nset,range,res,alpha,pair_blocks)
rowMeans(pvals_root_power[,,2]) #C=0.2 <-> ,,2 
saveRDS(rowMeans(pvals_root_power[,,2]), "../results/simulation_results/Torus_results/Torus_PI_power_result.rds")

###############################################################################

# settings with different sample sizes
sig =c(0.02)
list_npc=c(20,40,70,100)
nset=50
range=c(0,4)
res=40
alpha=0.05

load("PI_Torus_data/data_pi_torus_root.Rdata")
pair_blocks <- as.matrix(
  read.csv("../simulations/Torus/selected_pairs_torus_diffsamplesize.csv", header = FALSE)
)
## Power  
torus_pi_root_noised     = torus_pi_root[2]
torus_pi_root_not_noised = torus_pi_root[6]

pvals_root_power <- matrix(
  NA,
  nrow = length(list_npc),
  ncol = nset
)

for (n in 1:length(list_npc)) {
  
  npc = list_npc[n]
  
  test_results = ts_main_power_using_pair_block(
    torus_pi_root_not_noised,   # X = no-noise
    torus_pi_root_noised,       # Y = noise
    sig,
    npc,
    nset,
    range,
    res,
    alpha,
    pair_blocks
  )
  
  pvals_root_power[n,] = test_results[,,2]
}
rowMeans(pvals_root_power)
saveRDS(rowMeans(pvals_root_power), "../results/simulation_results/Torus_results/Torus_PI_power_result_diffsamples.rds")

###############################################################################
# validity 
sig = c(0.0)
npc = 25
nset = 200
npair = 1000
range = c(0,4)
res = 40
alpha = 0.05

load("PI_Torus_data/data_pi_torus_root.Rdata")
pair_blocks <- as.matrix(
  read.csv("../simulations/Torus/selected_pairs_torus_fpr.csv", header = FALSE)
)
target_data <- torus_pi_root[6]

test_results = ts_main_fpr(
  target_data,
  sig,
  npc,
  nset,
  npair,
  range,
  res,
  alpha,
  pair_blocks
)

mean(test_results[,,2]) # 0.021 

type1_error <- mean(test_results[,,2])
saveRDS(type1_error, "../results/simulation_results/Torus_results/Torus_PI_validity_result.rds")

