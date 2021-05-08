
##Coverage probability using Ramy Amer's BS antenna model, and the Nakagami CDF approximation.
MCCovProbRamyNakApprox = function(rxh,txh,angle,rmax,density,uptilt,T,gain,Nt,al,an,mal,man,alpha,beta,gamma,dx,cores){
  Tlin = 10^(T/10)
  M = (rmax^2+(rxh-txh)^2)
  
  rm = sqrt(-log(0.001)/(pi*density))
  
  r = seq(from=rm/dx,to=rm,by=rm/dx)
  
  iteration=function(q){
    Cov = zeros(nrow=1,ncol=2)
    if(angle<pi/2){
      ang = atan2(abs(rxh-txh),r[q])
      if((ang < (pi/2-(angle/2))) && (ang > (angle/2))){
        threshold = abs(rxh-txh)/tan(ang-(angle/2))  
      }
      else if(ang>(pi/2-(angle/2))){
        threshold = abs(rxh-txh)/tan(pi/2 - angle)  
      }
      else{
        threshold=rmax
      }
    }else{threshold=rmax}
    ds = (r[q]^2+(rxh-txh)^2)
    
    
    for (p in 1:2){
      ang = atan2((rxh-txh),r[q])
      gn=(1/Nt)*((sin(Nt*pi*(sin(ang))/2)^2)/(sin(pi*(sin(ang))/2)^2))
      if(p==1){
        s = (1/gn)*Tlin*ds^(al/2)
        m=mal
      }
      else{
        s = (1/gn)*Tlin*ds^(an/2)
        m=man
      }
      Cov[1,p] = 0
      
      bx= 0
      if(m==1){
        bx=1  
      }else if(m==2){
        bx=1.487
      }else if(m==3){
        bx = 1.81
      }else if(m==10){
        bx = 2.872
      }
      
      for(j in 1:m){
      dr = seq(from=r[q],to=threshold,by=(threshold-r[q])/200) 
      ##Laplace of LOS and NLOS interferers 
      Ll = vector(length=length(dr))
      Ln = vector(length=length(dr)) 
      for(i in 1:length(dr)){
        P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=dr[i],alpha=alpha,beta=beta,gamma=gamma)
        ang = atan2((rxh-txh),dr[i])
        gn=(1/Nt)*((sin(Nt*pi*(sin(ang))/2)^2)/(sin(pi*(sin(ang))/2)^2))
        g = bx*j*s*gn*(dr[i]^2+(rxh-txh)^2)^(-al/2)
        Ll[i] = (1-1/(g+1))*P*dr[i]*(dr[2]-dr[1])
        g = bx*j*s*gn*(dr[i]^2+(rxh-txh)^2)^(-an/2)
        Ln[i] = (1-1/(g+1))*(1-P)*dr[i]*(dr[2]-dr[1])
      }
      
      Cov[1,p] = Cov[1,p]+choose(m,j)*((-1)^(j+1))*exp(-density*angle*sum(Ll))*exp(-density*angle*sum(Ln))*exp(-bx*j*s*N/gain)
      }
      
    }
    
    P = Plos(htx=max(rxh,txh),hrx=min(rxh,txh),r=r[q],alpha=alpha,beta=beta,gamma=gamma)
    
    Cov[1,1] = Cov[1,1]*P
    Cov[1,2] = Cov[1,2]*(1-P)
    Cov[1,1] = (Cov[1,1]+Cov[1,2])*2*pi*density*r[q]*exp(-pi*density*r[q]^2)*(r[2]-r[1])
    
    return(Cov[1,1])
  }
  
  X=1:length(r)
  c = mclapply(X=X,FUN=iteration,mc.cores=cores)
  
  Cov = 0
  for(k in 1:length(r)){
    Cov = Cov+c[[k]]
  }
  
  return(Cov)
}



