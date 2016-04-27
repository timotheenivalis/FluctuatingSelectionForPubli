F_initPl<-function(MaxLitter=4,mums=mums,mother=mother){
  mothernb<-length(mums)
  pld<-matrix(data=0,nrow=nind,ncol=MaxLitter*mothernb)
  pli<-matrix(data=NA,nrow=nind,ncol=MaxLitter*mothernb)
  
  for(m in 1:mothernb)#this assumes that the first individual must come from the first litter!! (as it should be indeed now that data are re-ordered)
  {
    ind<-which(mother==m)
    for (i in 1:length(ind))
    {
      if (i==1)
      {
        pld[ind[i],((mother[ind[i]]-1)*MaxLitter+1):(mother[ind[i]]*MaxLitter)]<-c(1,rep(x=0,times=MaxLitter-1))
        pli[ind[i],((mother[ind[i]]-1)*MaxLitter+1):(mother[ind[i]]*MaxLitter)]<-NA
      }else{
        pld[ind[i],((mother[ind[i]]-1)*MaxLitter+1):(mother[ind[i]]*MaxLitter)]<-NA
        pli[ind[i],((mother[ind[i]]-1)*MaxLitter+1):(mother[ind[i]]*MaxLitter)]<-runif(MaxLitter,0,1)
      }
    }
  }
  
  initmeanLB<-vector(length=mothernb*MaxLitter)
  littermum<-vector()
  for(i in 1:length(mums))
  {
    littermum[(1+(i-1)*MaxLitter):(1+i*MaxLitter)] <-  i
  }
  
  for (i in 1:length(initmeanLB))
  {
    if(i==1)
    {
      initmeanLB[i]<-runif(1,0,40)
    }
    else
    {
      if(littermum[i]!=littermum[i-1])
      {
        initmeanLB[i]<-runif(1,0,40)
      }
      else
      {
        initmeanLB[i]<-runif(1,initmeanLB[i-1]+20,initmeanLB[i-1]+60)  
      }
    }  
  }
  toreturn<-list(pld=pld,pli=pli,initmeanLB=initmeanLB,littermum=littermum)
  return(toreturn)
}

##############################################%%
F_initPlsparse<-function(MaxLitter=4,mums=mums,mother=mother){
  mothernb<-length(mums)
  pld<-matrix(data=0,nrow=nind,ncol=MaxLitter)
  pli<-matrix(data=NA,nrow=nind,ncol=MaxLitter)
  
  for(m in 1:mothernb)#this assumes that the first individual must come from the first litter!! (as it should be indeed now that data are re-ordered)
  {
    ind<-which(mother==m)
    for (i in 1:length(ind))
    {
      if (i==1)
      {
        pld[ind[i],1:MaxLitter]<-c(1,rep(x=0,times=MaxLitter-1))
        pli[ind[i],1:MaxLitter]<-NA
      }else{
        pld[ind[i],1:MaxLitter]<-NA
        pli[ind[i],1:MaxLitter]<-runif(MaxLitter,0,1)
      }
    }
  }
  
  initmeanLB<-matrix(data = 0,nrow = mothernb,ncol = MaxLitter)
  
  for (i in 1:nrow(initmeanLB))
  {
    initmeanLB[i,1]<-runif(1,0,40)
    for (j in 2: ncol(initmeanLB))
      {
        initmeanLB[i,j]<-runif(1,initmeanLB[i,j-1]+20,initmeanLB[i,j-1]+60)
      }
  }
  toreturn<-list(pld=pld,pli=pli,initmeanLB=initmeanLB)
  return(toreturn)
}