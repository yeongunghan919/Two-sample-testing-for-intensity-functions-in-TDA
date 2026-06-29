library(TDA)
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(script_dir)
# generate data
sig=c(0.05,0.10,0.15,0.20) # noise level
npc=50 # number of point clouds per set
nset=100 # number of sets

set.seed(47) # seed

# function to generate unit circle
onecircle = function(){
  theta=runif(50)*2*pi
  x=cos(theta)
  y=sin(theta)
  circledat=cbind(x,y)
  return(circledat)
}

# generate data
onepc=onepd=list()
for (ii in 1:length(sig)){
  onepc[[ii]]=onepd[[ii]]=list()
  for (jj in 1:(npc*nset)) {
    # generate point cloud
    onepc[[ii]][[jj]]=onecircle()+matrix(rnorm(50,0,sig[ii]),50,2)
    # compute persistence diagram
    onediag=ripsDiag(onepc[[ii]][[jj]],maxdimension=1,maxscale=3)
    df= as.data.frame(unclass(onediag$diagram))
    colnames(df) = NULL 
    onepd[[ii]][[jj]]=unlist(df,use.names=FALSE) #we transform each diagram to data.frame for transforming json 
  }
}




set.seed(45) # seed

# function to generate two circles
twocircle = function(){
  num_small=sum(runif(50)>0.5)
  theta=runif(50)*2*pi
  x1=cos(theta[1:num_small])*0.9
  y1=sin(theta[1:num_small])*0.9
  
  x2=cos(theta[(num_small+1):50])*1.1
  y2=sin(theta[(num_small+1):50])*1.1
  circledat=cbind(c(x1,x2),c(y1,y2))
  return(circledat)
}

# generate data
twopc=twopd=list()
for (ii in 1:length(sig)){
  twopc[[ii]]=twopd[[ii]]=list()
  for (jj in 1:(npc*nset)) {
    # generate point cloud
    twopc[[ii]][[jj]]=twocircle()+matrix(rnorm(50,0,sig[ii]),50,2)
    # compute persistence diagram
    twodiag=ripsDiag(twopc[[ii]][[jj]],maxdimension=1,maxscale=3)
    df =as.data.frame(unclass(twodiag$diagram)) 
    colnames(df)=NULL
    twopd[[ii]][[jj]]=unlist(df,use.names=FALSE) #we transform each diagram to data.frame for transforming json
  }
}

saveRDS(onepd, "Rdata/npc=50_nset=100onepd.Rdata")
saveRDS(twopd, "Rdata/npc=50_nset=100twopd.Rdata")
library(jsonlite)

# JSON으로 변환
writeLines(toJSON(onepd),"json/npc=50_nset=100onepd.json")
writeLines(toJSON(twopd),"json/npc=50_nset=100twopd.json")


