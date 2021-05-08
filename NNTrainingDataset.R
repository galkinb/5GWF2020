##UAV is aware of 10 candidate BSs, and 20 interfering BSs
#Here we generate the dataset for training the NN to determine which candidate BS will have the best directional antenna SINR
rm(list=ls())
library(spatstat)
library(VGAM)
library(matlab)
#library(keras)
library(parallel)
library(hypergeo)
source("MiscFunctions.R")
#source("ClosedFormFunctions.R")

MCtrials = 500000

#building params
alpha = 0.5
beta = 300
gamma = 20

#building parameters
buildDens = 300/(1000^2) #density
buildWidth = 40 
buildR = buildWidth/(2*sin(pi/4)) #radius of the circular footprint of the building

rdist = 300
etilt = 10

#width of the simulation window, density of BSs
windowWidth = 5000
BHdensity = 5/(1000^2)


#number of antenna elements in the BS antennas
Nt= 8

#coverage probability threshold
Tu = 0

#LOS/NLOS pathloss exponents, LOS/NLOS fading parameter, BS transmit power (Watts), 
al = 2.1
an = 4
mal = 1
man = 1
BStx = 40

##Transmit frequency (for near-field pathloss calculation)
Freq = 2*10^9

#Noise power
N = -174+10*log10(20*10^6)+10
N = 10^(N/10)/1000

#MC cores for multi-core calculation
cores =8

#number BS candidates
BScand = 10
BSi = 20

##near-field pathloss parameter
K = (((3*10^8)/Freq)/(4*pi))^2

#Raytracing function for checking LOS
isLOS = function(buildings,buildR,buildH,x0,y0,x,y,h,BSh){
  angle = atan2((y-y0),(x-x0))
  dist = sqrt((x-x0)^2+(y-y0)^2)
  
  build = buildings
  build = shift(build,c(-x0,-y0))
  build = rotate(build,angle=-angle)
  
  buildX = build$x
  buildY = build$y
  foo = which(buildX<dist)
  buildX = buildX[foo]
  buildY = buildY[foo]
  buildH = buildH[foo]
  
  foo = which(buildX>0)
  buildX = buildX[foo]
  buildY = buildY[foo]
  buildH = buildH[foo]
  
  foo = which(abs(buildY)<=buildR)
  buildX = buildX[foo]
  buildY = buildY[foo]
  buildH = buildH[foo]
  
  foo = buildH>((abs(h-BSh)*(buildX/dist)+min(BSh,h)))
  if(length(which(foo==TRUE))>0){
    return(FALSE)
  }
  else{
    return(TRUE)
  }
}

##dataset generation
iteration = function(m){
 if(m%%1000==0){
    print(m)
 }
  
  UAVBHBW = pi/4
  UAVBHgain = 4*pi/((UAVBHBW/2)^2)
  
 
  gamma= 20
  h = runif(n=1,min=0,max=300)
  BSh=30

  ##window in which we simulate BS distribution
  sWindow = owin(xrange=c(-windowWidth/2,windowWidth/2),yrange=c(-windowWidth/2,windowWidth/2))
  gen=FALSE
  while(gen==FALSE){
    BHppp = rpoispp(lambda=BHdensity,win=sWindow)
    if(BHppp$n>5){
      gen=TRUE
    }
  }
  
  buildings = gridcenters(window=sWindow,nx=floor(sqrt(buildDens*(5000^2))),ny=ceil(sqrt(buildDens*(5000^2))))#rpoispp(lambda=buildDens,win=windowBS)
  build = rpoispp(lambda=buildDens,win=sWindow)
  build$x = buildings$x
  build$y = buildings$y
  build$n = length(buildings$x)
  buildings = build
  buildH = rrayleigh(n=buildings$n,scale=gamma)
  
  ##which network a given BS belongs to
  ##Note, I wrote this code for multi-network scenarios where the UAV can choose from several different operator networks
  ##In the paper it's just one network
  ##Nonetheless the dataset and the resulting neural network is capable of distinguishing between different networks, we simply don't use that function in the simulations
  #all the networks set to ID 1
  whichNetwork=ones(nrow=BHppp$n,ncol=1) #floor(runif(n=BHppp$n,min=1,max=4))
  
  LOS = vector(length=BHppp$n)

  measuredpower = vector(length=BHppp$n)
  load = round(runif(n=BScand,min=0,max=100)) ##the number of RBs that are in use for a BS (we don't use this in the current version)
  connectedto=floor(runif(n=1,min=1,max=(BScand+1))) #the UAV is connected to one of the BScand closest BSs (so it then decides whether to stay connected or change)
  softhandoverpenalty=runif(n=1,min=0.1,max=1) #soft handover penalty (intra-operator)
  hardhandoverpenalty = runif(n=1,min=0,max=softhandoverpenalty) #penalty for inter-operator handover
  numInterferers = zeros(nrow=BScand,ncol=BSi)
  distances = vector(length=BScand)
  
  
  #for each BS get the distance, channel type (LOS/NLOS) and the signal power received by the omnidirectional antenna
  rdist = vector(length=BHppp$n)
  for(i in 1:BHppp$n){
    foo = runif(n=1,min=0,max=1)
    if(isLOS(buildings=buildings,buildR=buildR,buildH=buildH,x0=0,y0=0,x=BHppp$x[i],y=BHppp$y[i],h=h,BSh=BSh)){
      LOS[i]=TRUE
    }
    else{LOS[i]=FALSE
    }
    
    #get vertical antenna gain
    angle = (atan2(BSh-h,sqrt((BHppp$x[i])^2+(BHppp$y[i])^2)))
    rdist[i] = (sqrt((BHppp$x[i])^2+(BHppp$y[i])^2))
    g=(1/Nt)*((sin(Nt*pi*(sin(angle))/2)^2)/(sin(pi*(sin(angle))/2)^2))
    
    #calculate measured power
    g = g*BStx*K
    if(LOS[i]==TRUE){
      measuredpower[i]=g*(sqrt((BHppp$x[i])^2+(BHppp$y[i])^2+(BSh-h)^2))^(-al)  
    }else{measuredpower[i]=g*(sqrt((BHppp$x[i])^2+(BHppp$y[i])^2+(BSh-h)^2))^(-an)}
  }
  
  #now get the omnidirectional SINR for each BS
  oSINR = vector(length=BHppp$n)
  for(i in 1:BHppp$n){
  foo = 1:BHppp$n
  foo = foo[foo!=i]
  #foo = (whichNetwork[foo]==whichNetwork[i])
  oSINR[i] = measuredpower[i]/(sum(measuredpower[foo])+N)
  }
  
  
  iBSBHHeight = vector(length=BHppp$n)
  iBSBHHeight[LOS==TRUE]=0
  iBSBHHeight[LOS==FALSE]=Inf
  
  order =order(rdist,decreasing=FALSE)
  
  
  achieveableRate = vector(length=BScand)
  SINR = vector(length=BScand)
  RBsneeded = vector(length=BScand)
  
  ##get the rate that would be achieved from each BS through the directional antenna
  for(i in 1:BScand){
    BHdist = sqrt((BHppp$x)^2+(BHppp$y)^2)
    BHBS = c(BHppp$x[order[i]],BHppp$y[order[i]])

    ind = order[i]
    BHint = 1:BHppp$n
    
    
    BHint = BHint[BHint!=ind] #indices of interfering BSs
    BSdist = BHdist[ind] #distance to serving BS
    BHdist = BHdist[BHint] #distances to interfering BSs
    
    distances[i]=BSdist/1000
    
    hopt = h
    angle = atan2(hopt-BSh,BSdist)
    
    #exclude the BSs that are outside the antenna radiation lobe
    if(UAVBHBW<pi/2){
      if((angle < (pi/2-(UAVBHBW/2))) && (angle > (UAVBHBW/2))){
        uthreshold = (hopt-BSh)/tan(angle-(UAVBHBW/2))  
      }
      else if(angle>(pi/2-(UAVBHBW/2))){
        uthreshold = (hopt-BSh)/tan(pi/2 - UAVBHBW)  
      }
      else{
        uthreshold = windowWidth
      }
    }
    else{uthreshold = windowWidth}
    
    if(angle<(pi/2-UAVBHBW/2)){
      lthreshold = (hopt-BSh)/tan(angle+(UAVBHBW/2))  
    }else{lthreshold=0}
    
    BHint = BHint[find(BHdist<=uthreshold)]
    BHdist = BHdist[BHdist<=uthreshold]
    BHint = BHint[find(BHdist>=lthreshold)]
    BHdist = BHdist[BHdist>=lthreshold]
    
    BHint = getInt2(x=c(0,0),int=BHint,BHBS=BHBS,grid=BHppp,UAVBHBW=UAVBHBW)
    BHH=iBSBHHeight[BHint]
    BHint = cbind(BHppp$x[BHint],BHppp$y[BHint])
    
    #store the distances of the BSi closest interfering BSs in the numInterfers matrix
    if(length(BHint)>0){
    iDist = sqrt((BHint[,1]/1000)^2+(BHint[,2]/1000)^2)
    foo = order(iDist,decreasing=FALSE)
    for(j in 1:min(BSi,length(foo))){
      numInterferers[i,j] = 1/(iDist[foo[j]])
    }
    }
    
    ##get spectral efficiency
    specEff = getDualRateRamyAntenna(x=c(0,0,hopt),BHBS=BHBS,BSh=BSh,withFading=FALSE,LOS=LOS[ind],iBSBH=cbind(BHint[,1],BHint[,2]),BHtilt=BHtilt,iBSBHHeight=BHH,Nt=Nt,al=al,an=an,mal=mal,man=man,PWRgain=BStx*UAVBHgain*K,N=N,alpha=alpha,beta=beta,gamma=gamma)
    RBsneeded[i] = 
    SINR[i]=2^(specEff)-1 #observed SINR

  }
  foo = 1:BScand
  ##return the observations for the 5 BSs, and label of which BS has the highest directional SINR
  ##note that we're storing the resource block load metrics as well and the network the BS belongs to, even though we don't use those parameters
  return(c(oSINR[order[foo]],measuredpower[order[foo]],load[foo],as.vector(t(numInterferers)),whichNetwork[order[foo]],distances[foo],softhandoverpenalty,hardhandoverpenalty,h,gamma,BSh,UAVBHBW,SINR[connectedto],which.max(SINR))) #return observed data
}

X=1:MCtrials
opt = mclapply(X=X,FUN=iteration,mc.cores=cores)

  results=zeros(nrow=MCtrials,(BScand*(5+BSi)+8))
#  
  for(k in 1:MCtrials){
    results[k,] = opt[[k]]

  }
save(results,file="trainingDataset.RData",version=2)

