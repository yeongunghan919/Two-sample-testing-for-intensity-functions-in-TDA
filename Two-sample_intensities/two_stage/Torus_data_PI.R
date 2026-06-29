script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)
source("functions.R")

# settings
res=40 # resolution 40 by 40
rangex=c(0,2)
rangey=c(0,2)
h=0.15



## load data
torus_pd = readRDS("../Rdata/npc=50_nset=100torus_pd.Rdata")

torus_pi_root = list()

for (ii in 1:6){
  torus_pi_root[[ii]] = list()
  for (jj in 1:5000){
    pd=torus_pd[[ii]][[jj]]
    pdmat=data.frame(t(matrix(pd, nrow = 3, byrow = TRUE))) #revised version 
    colnames(pdmat)=c("dimension","birth","death") 
    pd1 = pdmat %>%  
      dplyr::filter(dimension==1) %>%
      dplyr::select(-dimension)
    pd1$death = pd1$death - pd1$birth
    
    # weights
    arctanweight = arctanw(pd1,0.5,0.5)
    linearweight = linearw(pd1)
    rootweight = polyweight(pd1,n=1/2)
    
    # compute persistence image
    pi1_root=pers.image(pd=pd1, rangex=rangex, rangey=rangey, wgt=rep(1,nrow(pd1)), nbins=res, h=h)
    torus_pi_root[[ii]][[jj]] = pi1_root
  }
}

## save persistence images
save(torus_pi_root, file=paste0("PI_Torus_data/data_pi_torus_root.Rdata"))
