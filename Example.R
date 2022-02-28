###### Setting 1

library(foreach)
library(doParallel)
library(doRNG)

cl=makeCluster(10, outfile="")
registerDoParallel(cl)

library(MASS)
library(Matrix)
library(glmnet)
library(compositions)
library(lpSolve)
library(scout)
library(hdi)

source('imputeglm.predict.R')
source('cv.MBI_p1_c.R')
source('cv.MBI_p2_osqp_2_c.R')
source('myFDR.R')

n=150
ratio_x=3
nx=ratio_x*n
n1=ratio_x*30
n2=ratio_x*70
n3=ratio_x*25
n4=ratio_x*25

p1=40
p2=45
p3=115
p=p1+p2+p3
q1=2
q2=2
q3=5
q=q1+q2+q3  

a=c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)
beta=c(0.2, 0.5, 0.8, 1, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7, 8, 9, 10) 

TT=250
nmethod=12

alpha=0.05
debiased=c(p1+p2+1,p)
pd=length(debiased)

set.seed(1000)
seed_each=sample.int(5000000, 100000)
seed_count=1

s1=2
b=1

sigma1=diag(1,p,p)
sigma1[1:(p1+p2),1:(p1+p2)]=diag(1-a[s1],p1+p2,p1+p2)+matrix(a[s1],p1+p2,p1+p2)

true=c(beta[b]*rep(1,q1), rep(0,p1-q1), 
       beta[b]*rep(1,q2), rep(0,p2-q2),
       beta[b]*rep(1,q3), rep(0,p3-q3))


set.seed(seed_each[seed_count])
x00=mvrnorm(nx,rep(0,p),sigma1)          ## original data
x=x00
noise=rnorm(nx,0,1)
y0 = noise+as.numeric(x%*%true)

prob0=rep(1,n)
prob0=exp(-10*y0)
if (sum(prob0==Inf)>0){
  prob0[prob0==Inf]=max(prob0[!prob0==Inf])
}
group1=sample.int(nx,n1, prob = prob0)
group2=sample.int(nx-n1,n2)
group3=sample.int(nx-n1-n2,n3)

x1=x
index_g1=NULL
index_g2=NULL
index_g3=NULL
index_g4=NULL
for (k1 in 1:nx) {
  if (sum((1:nx)[-group1][group2]==k1)>0) {
    x1[k1,(p1+p2+1):(p1+p2+p3)]=NA
    index_g2=c(index_g2, k1)
  } else if (sum((1:nx)[-group1][-group2][group3]==k1)>0) {
    x1[k1,(p1+1):(p1+p2)]=NA
    index_g3=c(index_g3, k1)
  } else if (sum((1:nx)[-group1][-group2][-group3]==k1)>0) {
    x1[k1,1:p1]=NA
    index_g4=c(index_g4, k1)
  } else {
    index_g1=c(index_g1, k1)
  }
}

if ((length(index_g1)!=n1)|(length(index_g2)!=n2)|(length(index_g3)!=n3)|(length(index_g4)!=n4)) {
  print('Error index_g')
}

y=rep(NA, nx)
index_y=c(index_g1[1:(n1/ratio_x)], index_g2[1:(n2/ratio_x)],
          index_g3[1:(n3/ratio_x)], index_g4[1:(n4/ratio_x)])
y[index_y] = y0[index_y]


model1_1=cv.MBI_p1(x=x1, y=y, x.ori=x, start_source=c(1,p1+1, p1+p2+1), standardize = TRUE,
                   lambda=exp(seq(log(0.1),log(1),length.out = 10)), 
                   lambda2=sqrt(log(p)/n)/10,  
                   debiased=debiased, parallel = TRUE, parallel_LP=TRUE, alpha=alpha)
model1=cv.MBI_p2(x=x1, y=y, start_source=c(1,p1+1, p1+p2+1), standardize = TRUE,
                 x4=model1_1$x4, y4=model1_1$y4, x6=model1_1$x6, model=model1_1$fit,
                 best=model1_1$best, ug = model1_1$ug, pat = model1_1$pat, 
                 pat_ug = model1_1$pat_ug, N_g = model1_1$N_g,
                 lambda=model1_1$lambda, 
                 lambda2=sqrt(log(p)/n)/10,  
                 debiased=debiased, parallel = TRUE, alpha=alpha)

beta_d=model1$beta_final[debiased]
upb=model1$upperbound[debiased]
lwb=model1$lowerbound[debiased]


stopCluster(cl)


