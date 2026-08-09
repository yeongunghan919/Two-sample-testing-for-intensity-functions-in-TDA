script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)
source("functions.R")

# settings
res=40 # resolution 40 by 40
rangex=c(0,2)
rangey=c(0,2)
h=0.15



## load data
load("../data/Rdata/diagrams_for_equal_intensity_simul.RData")
length(P_diagrams)
length(Q_diagrams)

P_pi_linear = list()
Q_pi_linear = list()

for (jj in 1:length(P_diagrams)){
    pdP=P_diagrams[[jj]]
    pdQ=Q_diagrams[[jj]]
    
    pd1P=data.frame(matrix(pdP, ncol = 2))
    pd1Q=data.frame(matrix(pdQ, ncol = 2))
    colnames(pd1P)=c("birth","death") 
    colnames(pd1Q)=c("birth","death") 
    pd1P$death = pd1P$death - pd1P$birth
    pd1Q$death = pd1Q$death - pd1Q$birth
    
    
    # weights
    linearweightP = polyweight(pd1P,n=1)
    linearweightQ = polyweight(pd1Q,n=1)
    
    # compute persistence image
    pi1P_linear=pers.image(pd=pd1P, rangex=rangex, rangey=rangey, wgt=linearweightP, nbins=res, h=h)
    pi1Q_linear=pers.image(pd=pd1Q, rangex=rangex, rangey=rangey, wgt=linearweightQ, nbins=res, h=h)
    P_pi_linear[[jj]] = pi1P_linear
    Q_pi_linear[[jj]] = pi1Q_linear
}


## save persistence images
save(P_pi_linear, file=paste0("PI_Equal_intensity_data/P_pi_torus_linear.Rdata"))
save(Q_pi_linear, file=paste0("PI_Equal_intensity_data/Q_pi_torus_linear.Rdata"))
