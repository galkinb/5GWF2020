##Evaluate the neural network. Calculate the nearest-BS association coverage probability using the stochastic geometry derivations using our prior papers
##Then use simulations to evaluate the trained neural network, and the strongest omnidirectional SINR association
rm(list=ls())
library(spatstat)
library(VGAM)
library(hypergeo)
library(keras)
library(parallel)
library(matlab)
source("MiscFunctions.R")
source("ClosedFormFunctions.R")

#number of MC trials and cores to use for the calculation
MCtrials = 5000
cores=16

alpha = 0.5
beta = 300
gamma = 20


#building parameters
buildDens = 300/(1000^2)
buildWidth = 40
buildR = buildWidth/(2*sin(pi/4))
heightParam = 20

#UAV velocity
velocity = 10

#number of antenna elements
Nt= 8

BSh = 30
windowWidth = 5000
BHdensity = 5/(1000^2)

#UAVdensity = 25/(1000^2)#5*BSdensity
skipN = 1

UAVBHBW = pi*1/4

#coverage threshold
Tu = 0

#channel parameters, noise power, BS transmit power
al = 2.1
an = 4
mal = 1
man = 1
BStx = 40


#carrier frequency
Freq = 2*10^9

#noise power
N = -174+10*log10(20*10^6)+10
N = 10^(N/10)/1000

#BS antenna downtilt
BHtilt=-10


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

antennaGain = function(r,BSh,h,Nt){
  angle = (atan2(BSh-h,r))
  g=(1/Nt)*((sin(Nt*pi*(sin(angle))/2)^2)/(sin(pi*(sin(angle))/2)^2))
  g=g*BStx*K
  return(g)
}


maxSINR = 10.45007
#normalise the neural network input data (same was we did for the training)
normaliseData = function(P,Int,Dist,h){
  
  #normalise observed power
  foo=min(P)
  P= log10(P/foo)  
  
  mnfoo=min(P)
  mxfoo=max(P)
  if(mxfoo>0){
    P=(P-mnfoo)/(mxfoo-mnfoo)
  }
  
  #normalise interference score
  mnfoo=min(Int)
  mxfoo=max(Int)
  if(mxfoo>0){
    Int=(Int-mnfoo)/(mxfoo-mnfoo)
  }
  
  mnfoo=min(Dist)
  mxfoo=max(Dist)
  if(mxfoo>0){
    Dist=(Dist-mnfoo)/(mxfoo-mnfoo)
  }
  
  h = h/300

  
  return(c(P,Int,Dist,h))
}


#generate the model
model <- keras_model_sequential()
model %>%
  layer_dense(units=221,input_shape = c(221)) %>%
  layer_dense(units = 500, activation = 'relu') %>%
  layer_dense(units = 500, activation = 'relu') %>%
  layer_dense(units = 10, activation = 'softmax')




model %>% compile(
  optimizer = 'adamax', 
  loss = 'categorical_crossentropy',
  metrics = c('accuracy')
)


#load the weights from the training script
model %>% load_model_weights_hdf5("NNweights.h5")



#testing the model performance
h=seq(from=20,to=200,by=10)
covProb= zeros(nrow=length(h),ncol=length(Tu))
covProbnearest = zeros(nrow=length(h),ncol=length(Tu))
covProbmulti= zeros(nrow=length(h),ncol=length(Tu))
covProbML= zeros(nrow=length(h),ncol=length(Tu))
covProbMLmulti= zeros(nrow=length(h),ncol=length(Tu))


h=seq(from=20,to=200,by=6)
CFcoverageProb = zeros(nrow=length(h),ncol=length(Tu))
CFcoverageProbMulti = zeros(nrow=length(h),ncol=length(Tu))

nthBS = zeros(nrow=length(h),ncol=10)


h=seq(from=20,to=200,by=6)
#calculate nearest distance association using the analytical results
#we calculate it over a range of heights
for(j in 1:length(h)){
  for(l in 1:length(Tu)){ #this code will let us also calculate over a range of coverage thresholds (not used in these results)
    print(h[j])
    UAVBHgain = 4*pi/((UAVBHBW/2)^2)
    K = (((3*10^8)/Freq)/(4*pi))^2
  
    CFcoverageProb[j,l]=MCCovProbRamyNakApprox(rxh=h[j],txh=BSh,angle=UAVBHBW,rmax=windowWidth/2,density=BHdensity,uptilt=BHtilt,T=Tu[l],gain=BStx*UAVBHgain*K,Nt=Nt,al=al,an=an,mal=mal,man=man,alpha=alpha,beta=beta,gamma=gamma,dx=1000,cores=cores)
    save(CFcoverageProb,covProb,covProbML,h,file="HeightEvaluation.RData")
  }
}


h=seq(from=20,to=200,by=10)

##UAV antenna gain
UAVBHgain = 4*pi/((UAVBHBW/2)^2)
K = (((3*10^8)/Freq)/(4*pi))^2

for(j in 1:(length(h))){
  for(l in 1:length(Tu)){#this code will let us also calculate over a range of coverage thresholds (not used in these results)
    print(h[j])
    if(h[j]==30){h[j]=30.1} #if the UAV height is exactly equal to the BS height this causes an error in the calculation. I just offset the UAV height by a tiny amount of get around this

    mcS =zeros(nrow=MCtrials,ncol=10)
    mcInt =zeros(nrow=MCtrials,ncol=200)
    mcDist =zeros(nrow=MCtrials,ncol=10) 
    mcCov =zeros(nrow=MCtrials,ncol=10)
    mcbestSINR =zeros(nrow=MCtrials,ncol=1)
    mcnearestSINR =zeros(nrow=MCtrials,ncol=1)
    
    cPm= vector(length=MCtrials)
    cPmn = vector(length=MCtrials)
    cPmML= vector(length=MCtrials)
    
    whichNthBS = vector(length=MCtrials)
    
    uavx= 0
    uavy=0
    
    num = vector(length=MCtrials)
    
    BScand = 10
    BScandi = 20
    
    #generate MC trials for the NN
    iteration = function(m){
      if(m%%10000==0){
        print(m)
      }
      
      ##window in which we simulate BS distribution
      sWindow = owin(xrange=c(-windowWidth/2,windowWidth/2),yrange=c(-windowWidth/2,windowWidth/2))
      gen=FALSE
      while(gen==FALSE){
        BHppp = rpoispp(lambda=BHdensity,win=sWindow)
        if(BHppp$n>BScand){
          gen=TRUE
        }
      }
      
      buildings = gridcenters(window=sWindow,nx=floor(sqrt(buildDens*(5000^2))),ny=ceil(sqrt(buildDens*(5000^2))))#rpoispp(lambda=buildDens,win=windowBS)
      build = rpoispp(lambda=buildDens,win=sWindow)
      build$x = buildings$x
      build$y = buildings$y
      build$n = length(buildings$x)
      buildings = build
      buildH = rrayleigh(n=buildings$n,scale=heightParam)
      
      ##which network a given BS belongs to
      ##Note, I wrote this code for multi-network scenarios where the UAV can choose from several different operator networks
      ##I've decided to drop that in the paper we're writing as it gives us nothing of value. In the paper it's just one network
      ##Nonetheless the dataset and the resulting neural network is capable of distinguishing between different networks, we simply don't use that function in the simulations
      whichNetwork=ones(nrow=BHppp$n,ncol=1)#floor(runif(n=BHppp$n,min=1,max=4))
      
      LOS = vector(length=BHppp$n)
      
      measuredpower = vector(length=BHppp$n)
      load = round(runif(n=BScand,min=0,max=100)) ##the number of RBs that are in use for a BS (we don't use this in the current version)
      connectedto=floor(runif(n=1,min=1,max=(BScand+1))) #the UAV is connected to one of the BScand closest BSs (so it then decides whether to stay connected or change)
      numInterferers = zeros(nrow=BScand,ncol=BScandi)
      distances = vector(length=BScand)
      cov = vector(length=BScand)
      
      
      #for each BS get the distance, channel type (LOS/NLOS) and the signal power received by the omnidirectional antenna
      rdist = vector(length=BHppp$n)
      for(i in 1:BHppp$n){
        foo = runif(n=1,min=0,max=1)
        if(isLOS(buildings=buildings,buildR=buildR,buildH=buildH,x0=0,y0=0,x=BHppp$x[i],y=BHppp$y[i],h=h[j],BSh=BSh)){
          LOS[i]=TRUE
        }
        else{LOS[i]=FALSE
        }
        
        angle = (atan2(BSh-h[j],sqrt((BHppp$x[i])^2+(BHppp$y[i])^2)))
        rdist[i] = (sqrt((BHppp$x[i])^2+(BHppp$y[i])^2))
        g=(1/Nt)*((sin(Nt*pi*(sin(angle))/2)^2)/(sin(pi*(sin(angle))/2)^2))
        g = g*BStx*K
        if(LOS[i]==TRUE){
          measuredpower[i]=g*(sqrt((BHppp$x[i])^2+(BHppp$y[i])^2+(BSh-h[j])^2))^(-al)  
        }else{measuredpower[i]=g*(sqrt((BHppp$x[i])^2+(BHppp$y[i])^2+(BSh-h[j])^2))^(-an)}
      }
      
      #now get the SINR for each BS
      oSINR = vector(length=BHppp$n)
      for(i in 1:BHppp$n){
        foo = 1:BHppp$n
        foo = foo[foo!=i]
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
        
        hopt = h[j]
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
          iD = sqrt((BHint[,1]/1000)^2+(BHint[,2]/1000)^2)
          f = order(iD,decreasing=FALSE)
        for(k in 1:min(BScandi,length(f))){
            numInterferers[i,k] = 1/(iD[f[k]])
          }
        }
      #  numInterferers[i] = sum(1/sqrt((BHint[,1]/1000)^2+(BHint[,2]/1000)^2))
        
        ##get get the observed omnidirectional SINR with averaged out multipath effects (as the UAV is assumed to measure the values over a period of time and cancel out the small scale fading effects)
        specEff = getDualRateRamyAntenna(x=c(0,0,hopt),BHBS=BHBS,BSh=BSh,withFading=FALSE,LOS=LOS[ind],iBSBH=cbind(BHint[,1],BHint[,2]),BHtilt=BHtilt,iBSBHHeight=BHH,Nt=Nt,al=al,an=an,mal=mal,man=man,PWRgain=BStx*UAVBHgain*K,N=N,alpha=alpha,beta=beta,gamma=gamma)
        RBsneeded[i] = 0#for resource blocks if we were considering throughput with RBs, can ignore
        SINR[i]=2^(specEff)-1
        
        #get the instantaneous SINR (which IS affected by multipath effects) and which determines our coverage probability
        specEff = getDualRateRamyAntenna(x=c(0,0,hopt),BHBS=BHBS,BSh=BSh,withFading=TRUE,LOS=LOS[ind],iBSBH=cbind(BHint[,1],BHint[,2]),BHtilt=BHtilt,iBSBHHeight=BHH,Nt=Nt,al=al,an=an,mal=mal,man=man,PWRgain=BStx*UAVBHgain*K,N=N,alpha=alpha,beta=beta,gamma=gamma)
        RBsneeded[i] = 0#
        cov[i]=2^(specEff)-1
        
      }
        order =order(rdist,decreasing=FALSE)
      
      foo = 1:BScand
  
      ##return the observations for the 5 BSs, and label of which BS has the highest directional SINR
      ##note that we're storing the resource block load metrics as well and the network the BS belongs to, even though we don't use those parameters
      return(c(measuredpower[order[foo]],as.vector(t(numInterferers)),distances[foo],cov[foo],cov[1]))
    }
      
    #run the above function in multi-core over many iterations
      X=1:MCtrials
      opt = mclapply(X=X,FUN=iteration,mc.cores=cores)
      
    #as above, but this time we're generating MC results for the strongest BS association case
      iteration = function(m){
        if(m%%10000==0){
          print(m)
        }
        
        ##window in which we simulate BS distribution
        sWindow = owin(xrange=c(-windowWidth/2,windowWidth/2),yrange=c(-windowWidth/2,windowWidth/2))
        gen=FALSE
        while(gen==FALSE){
          BHppp = rpoispp(lambda=BHdensity,win=sWindow)
          if(BHppp$n>BScand){
            gen=TRUE
          }
        }
        
        buildings = gridcenters(window=sWindow,nx=floor(sqrt(buildDens*(5000^2))),ny=ceil(sqrt(buildDens*(5000^2))))#rpoispp(lambda=buildDens,win=windowBS)
        build = rpoispp(lambda=buildDens,win=sWindow)
        build$x = buildings$x
        build$y = buildings$y
        build$n = length(buildings$x)
        buildings = build
        buildH = rrayleigh(n=buildings$n,scale=heightParam)
        
        ##which network a given BS belongs to
        ##Note, I wrote this code for multi-network scenarios where the UAV can choose from several different operator networks
        ##I've decided to drop that in the paper we're writing as it gives us nothing of value. In the paper it's just one network
        ##Nonetheless the dataset and the resulting neural network is capable of distinguishing between different networks, we simply don't use that function in the simulations
        whichNetwork=ones(nrow=BHppp$n,ncol=1)
        
        LOS = vector(length=BHppp$n)
        
        measuredpower = vector(length=BHppp$n)
        load = round(runif(n=BScand,min=0,max=100)) ##the number of RBs that are in use for a BS (we don't use this in the current version)
        connectedto=floor(runif(n=1,min=1,max=(BScand+1))) #the UAV is connected to one of the BScand closest BSs (so it then decides whether to stay connected or change)
        softhandoverpenalty=runif(n=1,min=0.1,max=1) #soft handover penalty (intra-operator)
        hardhandoverpenalty = runif(n=1,min=0,max=softhandoverpenalty) #penalty for inter-operator handover
        numInterferers = zeros(nrow=BScand,ncol=BScand)
        distances = vector(length=BScand)
        cov = vector(length=BScand)
        
        
        #for each BS get the distance, channel type (LOS/NLOS) and the signal power received by the omnidirectional antenna
        rdist = vector(length=BHppp$n)
        for(i in 1:BHppp$n){
          foo = runif(n=1,min=0,max=1)
          if(isLOS(buildings=buildings,buildR=buildR,buildH=buildH,x0=0,y0=0,x=BHppp$x[i],y=BHppp$y[i],h=h[j],BSh=BSh)){
            LOS[i]=TRUE
          }
          else{LOS[i]=FALSE
          }
          
          angle = (atan2(BSh-h[j],sqrt((BHppp$x[i])^2+(BHppp$y[i])^2)))
          rdist[i] = (sqrt((BHppp$x[i])^2+(BHppp$y[i])^2))
          g=(1/Nt)*((sin(Nt*pi*(sin(angle))/2)^2)/(sin(pi*(sin(angle))/2)^2))
          g = g*BStx*K
          if(LOS[i]==TRUE){
            measuredpower[i]=g*(sqrt((BHppp$x[i])^2+(BHppp$y[i])^2+(BSh-h[j])^2))^(-al)  
          }else{measuredpower[i]=g*(sqrt((BHppp$x[i])^2+(BHppp$y[i])^2+(BSh-h[j])^2))^(-an)}
        }
        
        #now get the SINR for each BS
        oSINR = vector(length=BHppp$n)
        for(i in 1:BHppp$n){
          foo = 1:BHppp$n
          foo = foo[foo!=i]
          oSINR[i] = measuredpower[i]/(sum(measuredpower[foo])+N)
        }
        
        
        iBSBHHeight = vector(length=BHppp$n)
        iBSBHHeight[LOS==TRUE]=0
        iBSBHHeight[LOS==FALSE]=Inf
      #  order =order(rdist,decreasing=FALSE)
        
        
        achieveableRate = vector(length=BScand)
        SINR = vector(length=BScand)
        RBsneeded = vector(length=BScand)
        
        #determine which BS has the strongest omnidirectional SINR as connect UAV to that
        order =which.max(oSINR)
        
        
        BHdist = sqrt((BHppp$x)^2+(BHppp$y)^2)
        BHBS = c(BHppp$x[order],BHppp$y[order])
        
        ind = order
        BHint = 1:BHppp$n
        
        
        BHint = BHint[BHint!=ind]
        BSdist = BHdist[ind]
        BHdist = BHdist[BHint]
          
        hopt = h[j]
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
        
  
        ##get the directional SINR of the BS with the strongest omnidirectional SINR
        specEff = getDualRateRamyAntenna(x=c(0,0,hopt),BHBS=BHBS,BSh=BSh,withFading=TRUE,LOS=LOS[ind],iBSBH=cbind(BHint[,1],BHint[,2]),BHtilt=BHtilt,iBSBHHeight=BHH,Nt=Nt,al=al,an=an,mal=mal,man=man,PWRgain=BStx*UAVBHgain*K,N=N,alpha=alpha,beta=beta,gamma=gamma)
        strongestSINR=2^(specEff)-1
        
      #  order =order(rdist,decreasing=FALSE)
        
      #  foo = 1:BScand
        
        #return the strongest SINR value
        return(c(strongestSINR))
      }
      
      X=1:MCtrials
      opt2 = mclapply(X=X,FUN=iteration,mc.cores=cores)
      
      mP = 1:(BScand)
      Int = (mP[BScand]+1):(mP[BScand]+BScand*BScandi)
      Distances = (Int[BScand*BScandi]+1):(Int[BScand*BScandi]+BScand)
      C = (Distances[BScand]+1):(Distances[BScand]+BScand)
      sSINR = C[BScand]+1
     # nSINR = C[BScand]+2
    
      
      for(k in 1:MCtrials){
      mcS[k,]=opt[[k]][mP] #observed RSRP
      mcInt[k,]=opt[[k]][Int] #observed int distances
      mcDist[k,]=opt[[k]][Distances] #candidate BS distances
      mcCov[k,] =opt[[k]][C] #the corresponding directional SINR values
      mcbestSINR[k] =opt2[[k]] #the vector of the SINR values for the strongest oSINR association
      mcnearestSINR[k] =opt[[k]][sSINR] #nearest SINR values (not used here, because we use the )
      }
    
      #Now we go through the iterations and find the coverage probability, for the dumb strongest SINR association, and the neural network result
      associated=0
      associatedML=0
      pastSINRML = 0
      
      for(m in 1:MCtrials){
        if(mcbestSINR[m]>10^(Tu[l]/10)){cPm[m]=TRUE}
        else{cPm[m]=FALSE} 
        
        if(mcnearestSINR[m]>10^(Tu[l]/10)){cPmn[m]=TRUE}
        else{cPmn[m]=FALSE} 
        
        #get the normalised observed state for the given MC trial
        test = normaliseData(P=mcS[m,],Int=mcInt[m,],Dist=mcDist[m,],h=h[j])
        foo = rbind(test[1:221],vector(length=221))
        
        #predict the best BS via the NN
        newAss = model %>% predict_classes(foo)
        newAss=newAss[1]+1
        
        if(newAss>10 | newAss < 1){newAss=1}
       
        #get the directional SINR of that chosen BS
        cov = mcCov[m,(newAss)]
        if(cov>10^(Tu[l]/10)){cPmML[m]=TRUE}
        else{cPmML[m]=FALSE}
        
        #which BS did UAV associate with, closest, 2nd closest...?
        order = order(mcDist[k,],decreasing=FALSE)
        whichNthBS[m] = which(order==newAss)
      }
    
    #get coverage over the MC trials
    covProb[j,l] = sum(cPm)/(MCtrials)
    covProbnearest[j,l]= sum(cPmn)/(MCtrials)
    covProbML[j,l] = sum(cPmML)/(MCtrials)
    
    for(i in 1:10){
    nthBS[j,i] = sum(whichNthBS==i)/MCtrials
    }
    
    
    h=seq(from=20,to=200,by=6)
    plot(x=h,y=CFcoverageProb[,1],type='l',lwd=3,lty=1,cex.lab=1.2,cex.axis=1.2,xlab="UAV height (m)",ylab="Coverage Probability",ylim=c(0,1))
    h=seq(from=20,to=200,by=10)
    lines(x=h,y=covProb[,1],lwd=3,col='red')
    lines(x=h,y=covProbML[,1],lwd=3,col='blue')
    
  save(nthBS,covProbnearest,CFcoverageProb,covProb,covProbML,covProbMLmulti,h,file="HeightEvaluation.RData")
    
    }
    }  

foo = rbind(nthBS[4,],nthBS[9,],nthBS[14,])

#plot the barchart showing which BS the NN chose to associate with
#barplot(height=foo,beside=TRUE,col=c("black","red","blue"),legend=c("50m","100m","150m"),xlab=("n-th Closest BS"),ylab=("Probability of Association"),names.arg = 1:10,cex.names=1.2,cex.axis = 1.2,cex.lab=1.2)

grid(nx = NULL, ny = NULL, col = "darkgray", lty = "dotted",
     lwd = par("lwd"), equilogs = TRUE)


legend("bottomleft",                       # x-y coordinates for location of the legend  
       legend=c("Nearest Association","SINR Association (omni)","Neural Network Association"),      # Legend labels  
       col=c("black","red","blue"),   # Color of points or lines  
       lty=c(1,1,1,1,1,1),                    # Line type  
       lwd=c(8,8,8,8,8),                    # Line width  
       #       pch=c(15,16,17),
       cex=1.2
) 

