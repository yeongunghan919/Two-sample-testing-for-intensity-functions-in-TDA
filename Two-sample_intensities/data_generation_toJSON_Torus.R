library(TDA)
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)
# generate data
sig=c(0.01,0.02,0.03,0.04,0.05,0.0) # noise level, the last noise is not noise 
npc=50 # number of point clouds per set
nset=100 # number of sets


set.seed(47) # seed

#function of homo Poisson point process on Torus 
Poisson_homo_torus = function(lamb,R,r){
  total_intensity = 4*(pi**2)*R*r
  n = rpois(1,lamb*total_intensity)
  return(torusUnif(n,a=r,c=R))
}
torus_lamb=5
torus_R=2
torus_r=1
# generate data
torus_pc=torus_pd=list()
for (ii in 1:length(sig)){
  torus_pc[[ii]]=torus_pd[[ii]]=list()
  for (jj in 1:(npc*nset)) {
    # generate point cloud
    not_noise_data=Poisson_homo_torus(lamb=torus_lamb,R=torus_R,r=torus_r)
    torus_pc[[ii]][[jj]]=not_noise_data+matrix(rnorm(prod(dim(not_noise_data)),0,sig[ii]),dim(not_noise_data)[1],dim(not_noise_data)[2])
    # compute persistence diagram
    torus_diag=ripsDiag(torus_pc[[ii]][[jj]],maxdimension=1,maxscale=3)
    df= as.data.frame(unclass(torus_diag$diagram))
    colnames(df) = NULL 
    torus_pd[[ii]][[jj]]=unlist(df,use.names=FALSE) #we transform each diagram to data.frame for transforming json 
  }
}



saveRDS(torus_pd, "Rdata/npc=50_nset=100torus_pd.Rdata")
library(jsonlite)

# JSON
writeLines(toJSON(torus_pd),"json/npc=50_nset=100torus_pd.json")


