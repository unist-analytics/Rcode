#########################
### functions
#########################

hrha_from_hrha=function(hrha,n_L,n_H){
  n_H<-nrow(hrha)
  tmp=hrha  
  for (i in 1:n_L){
    ii=(n_H/n_L)*(i-1)+1
    ii1= (n_H/n_L)*(i)
    for (j in 1:n_L){
      jj=(n_H/n_L)*(j-1)+1
      jj1= (n_H/n_L)*(j)
      tmp[ii:ii1,jj:jj1]<-hrha[ii:ii1,jj:jj1]
      p1<-sample(c(ii:ii1))
      p2<-sample(c(jj:jj1))
      for(k in 1:(length(p1)/4*(n_H/n_L))){ 
        tmp[p1[k],p2[k]]<-0
      }
    }
  }
  return(tmp)
}

lrla_from_hrla=function(hrla){
  lrla=matrix(rep(0,n_L^2),n_L,n_L)  
  for (i in 1:n_L){
    ii=(n_H/n_L)*(i-1)+1
    ii1= (n_H/n_L)*(i)
    for (j in 1:n_L){
      jj=(n_H/n_L)*(j-1)+1
      jj1= (n_H/n_L)*(j)
      lrla[i,j]=sum(hrla[ii:ii1,jj:jj1])
    }
  }
  return(lrla)
}

delta_low2high<-function(delta){
  w2_hat=c()
  tmp=rep(delta/((n_H/n_L)^2),each=(n_H/n_L))
  for(i in 1:n_L){
    tmp2=tmp[((i-1)*n_H+1):(i*n_H)]
    tmp3<-c()
    for(j in 1:(n_H/n_L)){
      tmp3<-rbind(tmp3,tmp2)
    }
    w2_hat=rbind(tmp3,w2_hat)
  }
  return(as.vector(qX_mat(w2_hat)))
}

qX_mat=function(dd){
  t(dd[nrow(dd):1,])
}

data_generation_L<-function(mvariance,Ipp_real){
  
  xy_lowA=cbind(x=Ipp_real$x+rnorm(Ipp_real$n,0,mvariance),y=Ipp_real$y+rnorm(Ipp_real$n,0,mvariance))
  xy_lowA1=xy_lowA[(xy_lowA[,1]<=10)&(xy_lowA[,1]>=-5)&(xy_lowA[,2]<=15)&(xy_lowA[,2]>=0),]
  ppp.hrla=as.ppp(as.points(xy_lowA1[,1],xy_lowA1[,2]),c(-5,10,0,15))
  hrla=quadratcount(ppp.hrla,n_H,n_H)  #high resolution, low accuracy
  lrla=lrla_from_hrla(hrla)
  
  return(list(lrla=lrla,hrla=hrla))
}

mglmm_method<-function(lrha,lrla,hrla){
  mse=c() 
  X1<-seq(-5,10,length.out=(n_L+1))[1:n_L]+15/n_L/2
  X2<-seq(0,15,length.out=(n_L+1))[1:n_L]+15/n_L/2
  coords <-cbind(X1=rep(X1,n_L),X2=rep(X2,each=n_L))

  y=as.vector(qX_mat(lrha))  ## y_H
  X=as.vector(qX_mat(lrla))  ## x_H 
  
  
  df<-data.frame(y=y,X,coords)
  
  fit<-fitme(y~-1+X+Matern(1|X1+X2),data=df,family=Poisson(link="identity"), method="ML")
  random.effect<-ranef(fit)$`Matern(1 | X1 + X2)`
  
  w2.hat<-delta_low2high(random.effect)
  fixed.effect.hr<-cbind(as.vector(qX_mat(hrla)))%*%fixef(fit)
  
  return(fixed.effect.hr+w2.hat)
  
}

kdemethod<-function(Ipp_high_s){
  xx<-Ipp_high_s$x
  yy<-Ipp_high_s$y
  
  tmp1<-seq(-5,10,length.out=(n_H+1))[1:n_H]+15/n_H/2
  tmp2<-seq(0,15,length.out=(n_H+1))[1:n_H]+15/n_H/2
  
  f1 <- kde2d(xx,yy , n = c(n_H,n_H),lims=c(range(tmp1),range(tmp2)))
  
  result<-as.vector(f1$z*sum(true_hrha)/sum(f1$z))
  return(list(result=result,coords=expand.grid(f1$x,f1$y)))
}

pkrigemethod<-function(grids,hrha){
  gd<-as.geodata(data.frame(grids,as.vector(qX_mat(hrha))))
  
  pk<-pois.krige(gd,locations =gd$coords, krige = list(cov.pars = c(1,1), beta = 1),
                 mcmc.input = mcmc.control(S.scale = 0.2, thin = 1))
  return(pk$predict)
}

draw_contour_plot<-function(dd){
  dn<-density(dd)
  v<-ggplot(data.frame(dn))+ aes(x, y,z =value) 
  v+  geom_contour()+geom_text_contour(size = 5)+theme(axis.title.x=element_blank(),
                                                       axis.text.x=element_blank(),
                                                       axis.ticks.x=element_blank(),
                                                       axis.title.y=element_blank(),
                                                       axis.text.y=element_blank(),
                                                       axis.ticks.y=element_blank())
}

draw_density_plot<-function(dd,val,brks2,M=220){
  ggplot(dd)+
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),text = element_text(size = 18))+
    aes(x, y, fill = Freq) +geom_tile() +
    scale_fill_viridis_c(values=val,breaks=brks2, labels=brks2,limits = c(0,M))  
}
