#################################
## Packages
#################################

library(MASS)
library(spatstat)
library(splancs)
library(plyr)
library(geoRglm) #geoRglm_0.9-16.tar.gz #https://cran.r-project.org/src/contrib/Archive/geoRglm/
library(spaMM)
library(metR)
library(ggplot2)

source("Functions.R")

#########################
### Parameters
#########################

n_H<<-40
n_L<<-20
r<-0.125 
svar<-0.1
brks2<-c(0,25,50,75,100,125)

#########################
### Data Loading
#########################

load("simulated_data_example_2.RData")

true_hrha<- quadratcount(d0$Ipp_high, n_H, n_H)  
hrha<-quadratcount(d0$Ipp_high_s, n_H, n_H)
lrha<-quadratcount(d0$Ipp_high_s, n_L, n_L)
AA=as.vector(qX_mat(true_hrha))

MM<-125

#########################
### Simulated Data Plotting
#########################
draw_density_plot(data.frame(qX_mat(true_hrha)), c(0,0.5,1),brks2,M=MM) 
draw_contour_plot(d0$Ipp_high)  



#########################
### Model: KDE
#########################


fit_kde<-kdemethod(d0$Ipp_high_s,FF=bw.CvL)
pred_kde<-fit_kde$result
mean((AA-pred_kde)^2)


F2<-data.frame(qX_mat(true_hrha))
F2$Freq<-pred_kde
draw_density_plot(F2, c(0,0.5,1),brks2,M=MM)

#########################
### Model: poisson krige
#########################

grids<-fit_kde$coords

fit_pkrige<-pkrigemethod(grids,hrha)
pred_pkrige<-fit_pkrige*sum(true_hrha)/sum(fit_pkrige)
mean((AA-pred_pkrige)^2)

F3<-data.frame(qX_mat(true_hrha))
F3$Freq<-pred_pkrige
draw_density_plot(F3,c(0,0.5,1),brks2,M=MM)


#########################
### Model: Proposed Method
#########################

d1<-data_generation_L(svar,d0$Ipp_real)

ptm <- proc.time()
lambda<-mglmm_method(lrha,d1$lrla,d1$hrla)
t1<-proc.time() - ptm
t1

lambda[lambda<0]<-0
pred_mglmm<-lambda*(sum(true_hrha)/sum(lambda))
mean((AA-pred_mglmm)^2)

F1<-data.frame(qX_mat(true_hrha))
F1$Freq<-pred_mglmm
draw_density_plot(F1,c(0,0.5,1),brks2,M=MM) 


