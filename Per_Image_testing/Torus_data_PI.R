script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)
source("functions.R")

# settings
res=40 # resolution 40 by 40
rangex=c(0,2)
rangey=c(0,2)
h=0.15



## load data
torus_pd = readRDS("../data/Rdata/npc=50_nset=100torus_pd.Rdata")

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
    rootweight = polyweight(pd1,n=1/2)
    
    # compute persistence image
    pi1_root=pers.image(pd=pd1, rangex=rangex, rangey=rangey, wgt=rootweight, nbins=res, h=h)
    torus_pi_root[[ii]][[jj]] = pi1_root
  }
}

## save persistence images
save(torus_pi_root, file=paste0("PI_Torus_data/data_pi_torus_root.Rdata"))


######new
## load data

torus_pd = readRDS("../data/Rdata/npc=50_nset=100torus_pd.Rdata")

## weight settings
weight_list = list(
  constant = 0,
  root     = 1/2,
  linear   = 1
)

## compute and save persistence images for each weight

for (weight_name in names(weight_list)) {
  
  weight_power = weight_list[[weight_name]]
  
  torus_pi = list()
  
  for (ii in 1:6) {
    
    torus_pi[[ii]] = list()
    
    for (jj in 1:5000) {
      
      pd = torus_pd[[ii]][[jj]]
      
      pdmat = data.frame(
        t(matrix(pd, nrow = 3, byrow = TRUE))
      ) # revised version
      
      colnames(pdmat) = c("dimension", "birth", "death")
      
      pd1 = pdmat %>%
        dplyr::filter(dimension == 1) %>%
        dplyr::select(-dimension)
      
      pd1$death = pd1$death - pd1$birth
      
      ## weights
      weight = polyweight(pd1, n = weight_power)
      
      ## compute persistence image
      pi1 = pers.image(
        pd = pd1,
        rangex = rangex,
        rangey = rangey,
        wgt = weight,
        nbins = res,
        h = h
      )
      
      torus_pi[[ii]][[jj]] = pi1
    }
  }
  
  ## save persistence images
  save(
    torus_pi,
    file = paste0(
      "PI_Torus_data/data_pi_torus_",
      weight_name,
      ".Rdata"
    )
  )
}





