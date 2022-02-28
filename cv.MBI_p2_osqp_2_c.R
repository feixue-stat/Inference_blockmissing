library(compositions)
library(osqp)
library(foreach)
library(doParallel)
library(Matrix)

cv.MBI_p2 <- function (x, y, x.ori=NULL, start_source, standardize=TRUE, 
                       lambda, lambda2, r_va=0.2, nfolds=10, alpha=0.05, 
                       debiased=1:p, parallel=FALSE, b_count_max=10, 
                       use_true=FALSE, bound0=0.1,
                       x4, y4, x6, model, best, ug, pat, pat_ug, N_g) {
  
  # x: Block-wise missing covariates/design matrix
  # y: Vector of a partially observed response
  
  
  n=dim(x)[1]
  p=dim(x)[2]
  
  centerx=rep(0,p)
  scalex=rep(1,p)
  centery=0
  x2=x
  y2=y
  if (!is.null(x.ori)) {
    x.ori1=x.ori
  }
  if (standardize) {
    centery=mean(y, na.rm=TRUE)
    y2=y-centery
    for (j in 1:p) {
      centerx[j]=mean(x[,j], na.rm=TRUE)
      scalex[j]=sd(x[,j], na.rm = TRUE)
      x2[,j]=(x[,j]-centerx[j])/scalex[j]
    }
    if (!is.null(x.ori)) {
      x.ori1=x.ori
      for (j in 1:p) {
        x.ori1[,j]=(x.ori[,j]-centerx[j])/scalex[j]
      }
    }
  }
  
  # miss=is.na(x2)
  # miss_num=apply(miss,1,function(x) unbinary(paste(as.numeric(x), collapse="")))
  # 
  # x5=x2[order(miss_num),]
  # y5=y2[order(miss_num)]
  # if (!is.null(x.ori)) {
  #   x.ori2=x.ori1[order(miss_num),]
  # }
  # 
  # bi_miss=is.na(x5)[,start_source]
  # id=apply(bi_miss,1,function(x) unbinary(paste(as.numeric(x), collapse="")))
  # n_pat=length(unique(id))
  # pat_start=rep(1,n_pat)
  # table_pat=table(id)
  # count_tmp=1
  # for (i in 2:n_pat){
  #   count_tmp=count_tmp+table_pat[i-1]
  #   pat_start[i]=count_tmp
  # }
  # 
  # 
  # 
  # miss_pat=is.na(x5[pat_start, start_source])
  # Im=0
  # ug=rep(list(NULL),n_pat)
  # pat=NULL
  # pat_ug=NULL
  # for (i in 1:n_pat) {
  #   if (sum(miss_pat[i,])>0) {
  #     which_miss=which(miss_pat[i,])
  #     for (j in (1:n_pat)[-i]) {
  #       if (sum(miss_pat[j,which_miss])==0 & 
  #           sum(miss_pat[j,-which_miss])<length(miss_pat[j,-which_miss])) {
  #         ug[[i]]=c(ug[[i]], j)
  #       }
  #     }
  #   } else {
  #     ug[[i]]=i
  #   }
  #   if (is.null(ug[[i]])) {
  #     print('Error: each missing pattern should be able to be imputed')
  #   }
  #   Im=Im+length(ug[[i]])
  #   pat=c(pat, rep(i, length(ug[[i]])))
  #   pat_ug=c(pat_ug, ug[[i]])
  # } 
  # 
  # x3=array(0,dim=c(n,p,Im))
  # y3=matrix(0,n,Im)
  # index1=matrix(FALSE,n,Im)            # use which sample
  # for (j in 1:Im) {
  #   if (pat[j]<n_pat) {
  #     index_temp=pat_start[pat[j]]:(pat_start[pat[j]+1]-1)
  #   } else {
  #     index_temp=pat_start[pat[j]]:n
  #   }
  #   y3[index_temp,j]=y5[index_temp]
  #   x3[index_temp,,j] = x5[index_temp,]
  #   index1[index_temp,j]=TRUE
  # } 
  # 
  # 
  # ########### Multiple Blockwise Imputation ############
  # miss2=is.na(x5)
  # N_g=0              ######## Number of all estimating equations
  # for (i in 1:n_pat) {
  #   if (sum(miss_pat[i,])>0) {
  #     ind_y=which(is.na(x5[pat_start[i],]))
  #     for (j in ug[[i]]) {
  #       N_g=N_g+sum(!is.na(x5[pat_start[j],]))
  #       k=which((pat==i) & (pat_ug==j))
  #       #temp1 = which(!is.na(x5[pat_start[j],]))
  #       ind_x = which((!is.na(x5[pat_start[i],])) & (!is.na(x5[pat_start[j],])))
  #       #temp1[!temp1 %in% ind_y]
  #       if (length(ind_x)==0 | length(ind_y)==0 | k==0) {
  #         print('Error')
  #       }
  #       x3[index1[,k]==TRUE,ind_y,k] = imputeglm.predict(X=x5, ind_y=ind_y, ind_x=ind_x, miss=miss2, 
  #                                                        newdata = x3[index1[,k]==TRUE,ind_x,k])$PRED
  #     }
  #   } else {
  #     for (j in ug[[i]]) {
  #       N_g=N_g+sum(!is.na(x5[pat_start[j],]))
  #     }
  #   }
  # }
  # # > sum(is.na(x3))
  # # [1] 0
  # 
  # impest=NULL
  # if (!is.null(x.ori)) {
  #   for (j in 1:Im) {
  #     if (pat[j]<n_pat) {
  #       index_temp=pat_start[pat[j]]:(pat_start[pat[j]+1]-1)
  #     } else {
  #       index_temp=pat_start[pat[j]]:n
  #     }
  #     if (sum(miss2[pat_start[pat[j]],])>0) {
  #       im_temp=mean((x3[index_temp, miss2[pat_start[pat[j]],], j]
  #                     -x.ori2[index_temp, miss2[pat_start[pat[j]],]])^2)
  #       impest=c(impest, im_temp)
  #       
  #       if (use_true) {
  #         x3[index_temp, miss2[pat_start[pat[j]],], j] = x.ori2[index_temp, miss2[pat_start[pat[j]],]]
  #       }
  #     }
  #   }
  # }
  # 
  # 
  # ############## remove NAs in y ###############
  # index_y=which(!is.na(y5))
  # y6=y5[index_y]
  # x6=x5[index_y,]
  # # dim(x6)
  # y4=y3[index_y,]
  # x4=x3[index_y,,]
  # # dim(x4)
  # 
  # ################## Update n, n_pat, pat_start, ug, pat_ug, N_g ############
  # n=length(y6)
  # 
  # bi_miss=is.na(x6)[,start_source]
  # id=apply(bi_miss,1,function(x) unbinary(paste(as.numeric(x), collapse="")))
  # n_pat=length(unique(id))
  # pat_start=rep(1,n_pat)
  # table_pat=table(id)
  # count_tmp=1
  # for (i in 2:n_pat){
  #   count_tmp=count_tmp+table_pat[i-1]
  #   pat_start[i]=count_tmp
  # }
  # 
  # miss_pat=is.na(x6[pat_start, start_source])
  # Im=0
  # ug=rep(list(NULL),n_pat)
  # pat=NULL
  # pat_ug=NULL
  # N_g=0              ######## Number of all estimating equations
  # for (i in 1:n_pat) {
  #   if (sum(miss_pat[i,])>0) {
  #     which_miss=which(miss_pat[i,])
  #     for (j in (1:n_pat)[-i]) {
  #       if (sum(miss_pat[j,which_miss])==0 & 
  #           sum(miss_pat[j,-which_miss])<length(miss_pat[j,-which_miss])) {
  #         ug[[i]]=c(ug[[i]], j)
  #         N_g=N_g+sum(!is.na(x6[pat_start[j],]))
  #       }
  #     }
  #   } else {
  #     ug[[i]]=i
  #     N_g=N_g+sum(!is.na(x6[pat_start[i],]))
  #   }
  #   if (is.null(ug[[i]])) {
  #     print('Error: each missing pattern should be able to be imputed')
  #   }
  #   Im=Im+length(ug[[i]])
  #   pat=c(pat, rep(i, length(ug[[i]])))
  #   pat_ug=c(pat_ug, ug[[i]])
  # } 
  # 
  # 
  # foldid=NULL
  # for (i in 1:n_pat) {
  #   if (i<n_pat) {
  #     templ=length(pat_start[i]:(pat_start[i+1]-1))
  #   } else {
  #     templ=length(pat_start[i]:n)
  #   }
  #   foldid = c(foldid, sample(rep(seq(nfolds), length = templ)))
  # }
  # if (length(foldid)!=n) {
  #   print('Error')
  # }
  # 
  # # test_index=NULL
  # # for (i in 1:n_pat) {
  # #   if (i<n_pat) {
  # #     index_temp=pat_start[i]:(pat_start[i+1]-1)
  # #   } else {
  # #     index_temp=pat_start[i]:n
  # #   }
  # #   test_index = c(test_index, sample(index_temp, round(r_va*length(index_temp))))
  # # }
  # 
  # nlam=length(lambda)
  # PMSE=rep(0,nlam)
  # PMSE_fold=matrix(0,nfolds,nlam)
  # status_fold=matrix(0,nfolds,nlam)
  nlam2=length(lambda2)
  Pvar=matrix(0,nlam2,p)
  
  # for (t in 1:nfolds){
  #   test_index=which(foldid==t)
  #   
  #   x3_train=x4[-test_index,,]
  #   y_train=y4[-test_index,]
  #   x5_train=x6[-test_index,]
  #   x5_test=x6[test_index,]
  #   x3_test=x4[test_index,,]
  #   y_test=y4[test_index,]
  #   
  #   # sum(is.na(x3_train))
  #   # [1] 0
  #   # sum(is.na(x3_test))
  #   # [1] 0
  #   
  #   if (nlam>1) {
  #     model_fold = MBI_ini(x_train = x3_train, y_train = y_train, x5_train = x5_train, 
  #                          x_test = x3_test, y_test = y_test,
  #                          start_source = start_source, lambda = lambda,
  #                          ug = ug, pat = pat, pat_ug = pat_ug)
  #     PMSE_fold[t,]=model_fold$PMSE
  #     status_fold[t,]=model_fold$status
  #     #PMSE = PMSE + model_fold$PMSE
  #   }
  #   
  #   if (nlam2>1) {
  #     weight_fold = weight.v(x_train = x3_train, y_train = y_train, x5_train = x5_train, 
  #                            x_test=x3_test, y_test=y_test, x5_test=x5_test, 
  #                            start_source = start_source, lambda2 = lambda2,
  #                            ug = ug, pat = pat, pat_ug = pat_ug, N_g = N_g, 
  #                            debiased=debiased, parallel=parallel, b_count_max=b_count_max,
  #                            bound0=bound0)
  #     Pvar = Pvar + weight_fold$Pvar
  #   }
  #   print(paste('t=',t,sep = ''))
  # }
  # 
  # model = MBI_ini(x_train = x4, y_train = y4, x5_train = x6,
  #                 start_source = start_source, lambda = lambda,
  #                 ug = ug, pat = pat, pat_ug = pat_ug, res=TRUE,
  #                 scalex = scalex, centerx = centerx, centery = centery)
  
  # PMSE=apply(PMSE_fold, 2, sum)
  # best=which(PMSE==min(PMSE[model$status==0]))[1]
  lambda.min=lambda[best]
  beta.min = model$beta[best,]
  beta_out.min = model$beta_out[best,]
  intercept.min = model$intercept[best]
  sigma_squared = model$TMSE[best]
  
  n=dim(x6)[1]
  
  #print('MBI_ini finished')
  
  
  weight_model = weight.v(x_train = x4, y_train = y4, x5_train = x6,
                          start_source = start_source, lambda2 = lambda2,
                          ug = ug, pat = pat, pat_ug = pat_ug, N_g = N_g,
                          minimum=TRUE, debiased=debiased, parallel=parallel,
                          b_count_max=b_count_max, bound0=bound0)
  
  best2=rep(0,p)
  lambda2.min=rep(0,p)
  v.min=matrix(0,weight_model$N_g2,p)
  beta_final=beta.min
  upperbound=rep(0,p)
  lowerbound=rep(0,p)
  s_square=rep(0,p)
  pvalues=rep(0,p)
  test_stats=rep(0,p)
  for (j in debiased) {
    if (sum(weight_model$v_status[,j]==TRUE)==0) {
      best2[j]=nlam2
    } else {
      best2[j]=which(Pvar[,j]==min(Pvar[weight_model$v_status[,j]==TRUE,j]))[1]
    }
    lambda2.min[j]=lambda2[best2[j]]
    v.min[,j]=weight_model$v[best2[j],j,]
    beta_final[j]=t(v.min[,j])%*%(weight_model$b1-weight_model$A1[,-j]%*%beta.min[-j])/(t(v.min[,j])%*%weight_model$A1[,j])
    s_square[j]=sigma_squared*weight_model$Tvar[best2[j],j]
    upperbound[j]=beta_final[j]+sqrt(s_square[j])*qnorm(alpha/2, lower.tail=FALSE)/n
    lowerbound[j]=beta_final[j]-sqrt(s_square[j])*qnorm(alpha/2, lower.tail=FALSE)/n
    test_stats[j]=n*beta_final[j]/sqrt(s_square[j])
    pvalues[j]=2*pnorm(abs(test_stats[j]), lower.tail = FALSE)
  }
  
  beta_final_out=beta_final/scalex
  intercept_final = centery-crossprod(beta_final_out,as.numeric(centerx))
  upperbound_out=upperbound/scalex
  lowerbound_out=lowerbound/scalex
  
  returnlist=list("beta.min"=beta_out.min, "beta_final"=beta_final_out,
                  "beta_final0"=beta_final,
                  "fit"=model, "weight_model"=weight_model,
                  #"x4"=x4, "y4"=y4, "x6"=x6, "ug"=ug, "pat"=pat, "pat_ug"=pat_ug, "N_g"=N_g,
                  "best"=best, "best2"=best2,
                  "upperbound"=upperbound_out, "lowerbound"=lowerbound_out,
                  "centerx"=centerx, "centery"=centery, "scalex"=scalex,
                  "intercept.min"=intercept.min, "intercept_final"=intercept_final,
                  "lambda.min"=lambda.min, "lambda2.min"=lambda2.min,
                  #"PMSE"=PMSE, 
                  "Pvar"=Pvar, "v.min"=v.min,
                  "lambda"=lambda, "lambda2"=lambda2, "s_square"=s_square,
                  #"impest"=impest, 
                  "debiased"=debiased, "sigma_squared"=sigma_squared,
                  "test_stats"=test_stats, "pvalues"=pvalues)
  #"PMSE_fold"=PMSE_fold, "status_fold"=status_fold
  return(returnlist)
} 




MBI_ini <- function (x_train, y_train, x5_train, x_test=NULL, y_test=NULL, 
                     start_source, lambda, ug, pat, pat_ug, 
                     scalex=NULL, centerx=NULL, centery=NULL, res=FALSE, parallel_LP) {
  n_train=dim(x5_train)[1]
  p=dim(x5_train)[2]
  
  bi_miss=is.na(x5_train)[,start_source]
  id=apply(bi_miss,1,function(x) unbinary(paste(as.numeric(x), collapse="")))
  n_pat=length(unique(id))
  pat_start=rep(1,n_pat)
  table_pat=table(id)
  count_tmp=1
  for (i in 2:n_pat){
    count_tmp=count_tmp+table_pat[i-1]
    pat_start[i]=count_tmp
  }
  
  nlam=length(lambda)
  beta=matrix(0, nlam, p)
  PMSE=rep(0,nlam)
  TMSE=rep(0,nlam)
  status=rep(0,nlam)
  beta_out=matrix(0, nlam, p)
  intercept=rep(0,p)
  
  ############# Vectors and matrices for constraints ############
  b1=NULL
  A1=NULL
  Im=0
  for (i in 1:n_pat) {
    if (i<n_pat) {
      index_temp=pat_start[i]:(pat_start[i+1]-1)
    } else {
      index_temp=pat_start[i]:n_train
    }
    for (j in ug[[i]]) {
      ak = which(!is.na(x5_train[pat_start[j],]))
      k=which((pat==i) & (pat_ug==j))
      b1=c(b1, t(x_train[index_temp,ak,k])%*%y_train[index_temp,k]/length(index_temp))
      A1=rbind(A1, t(x_train[index_temp,ak,k])%*%x_train[index_temp,,k]/length(index_temp))
      # if (sum(x3[index_temp,,k]==0)>1) {
      #   print(j)
      # }
    }
    Im=Im+length(ug[[i]])
  }
  # length(b1)
  # [1] 328
  # dim(A1)
  # [1] 328  40
  
  # f.con=rbind(cbind(diag(1,p,p), diag(-1,p,p)), 
  #             cbind(diag(-1,p,p), diag(-1,p,p)),
  #             cbind(-A1, matrix(0,dim(A1)[1],p)),
  #             cbind(A1, matrix(0,dim(A1)[1],p)))
  # dim(f.con)
  # [1] 736  80
  # f.dir = rep('<=', dim(f.con)[1])
  # f.obj = c(rep(0,p), rep(1,p))
  
  f.mat= rbind(cbind(-A1, A1),
               cbind(A1, -A1))
  f.obj = rep(1, 2*p)
  f.dir = rep('<=', dim(f.mat)[1])
  
  if (parallel_LP) {
    results_l = foreach (l=1:nlam, .packages='Rglpk') %dopar% {
      # f.rhs = c(rep(0,2*p), lambda[l]-b1, lambda[l]+b1)
      f.rhs = c(lambda[l]-b1, lambda[l]+b1)
      
      results=Rglpk_solve_LP(f.obj, f.mat, f.dir, f.rhs, max = FALSE)
      #lp ("min", f.obj, f.con, f.dir, f.rhs)
      beta_l=results$solution[1:p] - results$solution[(p+1):(2*p)]
      status_l=results$status
      
      TMSE_l=0
      PMSE_l=0
      if (res) {
        count=0
        for (i in 1:n_pat) {
          if (i<n_pat) {
            index_temp=pat_start[i]:(pat_start[i+1]-1)
          } else {
            index_temp=pat_start[i]:n_train
          }
          for (j in ug[[i]]) {
            k=which((pat==i) & (pat_ug==j))
            TMSE_l=TMSE_l+sum((y_train[index_temp,k]-x_train[index_temp,,k]%*%beta_l)^2)
            count=count+length(index_temp)
          }
        }
        TMSE_l=TMSE_l/count
      }
      if (!is.null(x_test) & !is.null(y_test)) {
        for (i in 1:Im) {
          PMSE_l=PMSE_l+sum((y_test[,i]-x_test[,,i]%*%beta_l)^2)
        }
      }
      
      beta_out_l=rep(0,p)
      intercept_l=0
      if (!is.null(scalex) & !is.null(centerx) & !is.null(centery)) {
        beta_out_l=beta_l/scalex
        intercept_l = centery-crossprod(beta_out_l,as.numeric(centerx))
      }
      return(list('beta_l'=beta_l, 'status_l'=status_l, 'TMSE_l'=TMSE_l, 'PMSE_l'=PMSE_l,
                  'beta_out_l'=beta_out_l, 'intercept_l'=intercept_l))
    }
    for (l in 1:nlam) {
      beta[l,]=results_l[[l]]$beta_l
      status[l]=results_l[[l]]$status_l
      TMSE[l]=results_l[[l]]$TMSE_l
      PMSE[l]=results_l[[l]]$PMSE_l
      beta_out[l,]=results_l[[l]]$beta_out_l
      intercept[l]=results_l[[l]]$intercept_l
    }
  } else {
    for (l in 1:nlam) {
      # f.rhs = c(rep(0,2*p), lambda[l]-b1, lambda[l]+b1)
      f.rhs = c(lambda[l]-b1, lambda[l]+b1)
      
      results=Rglpk_solve_LP(f.obj, f.mat, f.dir, f.rhs, max = FALSE)
      #lp ("min", f.obj, f.con, f.dir, f.rhs)
      beta[l,]=results$solution[1:p] - results$solution[(p+1):(2*p)]
      status[l]=results$status
      
      if (res) {
        count=0
        for (i in 1:n_pat) {
          if (i<n_pat) {
            index_temp=pat_start[i]:(pat_start[i+1]-1)
          } else {
            index_temp=pat_start[i]:n_train
          }
          for (j in ug[[i]]) {
            k=which((pat==i) & (pat_ug==j))
            TMSE[l]=TMSE[l]+sum((y_train[index_temp,k]-x_train[index_temp,,k]%*%beta[l,])^2)
            count=count+length(index_temp)
          }
        }
        TMSE[l]=TMSE[l]/count
      }
      if (!is.null(x_test) & !is.null(y_test)) {
        for (i in 1:Im) {
          PMSE[l]=PMSE[l]+sum((y_test[,i]-x_test[,,i]%*%beta[l,])^2)
        }
      }
      
      if (!is.null(scalex) & !is.null(centerx) & !is.null(centery)) {
        beta_out[l,]=beta[l,]/scalex
        intercept[l] = centery-crossprod(beta_out[l,],as.numeric(centerx))
      }
    }
  }
  
  
  
  returnlist=list("beta"=beta, "beta_out"=beta_out, "intercept"=intercept,
                  "PMSE"=PMSE, "TMSE"=TMSE, 
                  "lambda"=lambda, "status"=status,
                  "A1"=A1, "b1"=b1)
  return(returnlist)
}



weight.v <- function (x_train, y_train, x5_train, x_test=NULL, y_test=NULL, x5_test=NULL, 
                      start_source, lambda2, ug, pat, pat_ug, N_g, minimum=FALSE, 
                      debiased, parallel, b_count_max, bound0) {
  n_train=dim(x5_train)[1]
  p=dim(x5_train)[2]
  
  bi_miss=is.na(x5_train)[,start_source]
  id=apply(bi_miss,1,function(x) unbinary(paste(as.numeric(x), collapse="")))
  n_pat=length(unique(id))
  pat_start=rep(1,n_pat)
  table_pat=table(id)
  count_tmp=1
  for (i in 2:n_pat){
    count_tmp=count_tmp+table_pat[i-1]
    pat_start[i]=count_tmp
  }
  
  N_g2=0
  for (i in 1:n_pat) {
    for (j in ug[[i]]) {
      N_g2_temp=sum(!is.na(x5_train[pat_start[j],]) & !is.na(x5_train[pat_start[i],]))
      N_g2=N_g2+N_g2_temp
    }
  }
  
  ############# Vectors and matrices for constraints ############
  ############# This does not depend on lambda
  A1=NULL
  b1=NULL
  W=matrix(0,N_g2,N_g2)
  if (!is.null(x_test) & !is.null(y_test) & !is.null(x5_test)) {
    W_test=matrix(0,N_g2,N_g2)
    n_test=dim(x_test)[1]
    
    bi_miss_test=is.na(x5_test)[,start_source]
    id_test=apply(bi_miss_test,1,function(x) unbinary(paste(as.numeric(x), collapse="")))
    pat_start_test=rep(1,n_pat)
    table_pat_test=table(id_test)
    count_tmp_test=1
    for (i in 2:n_pat){
      count_tmp_test=count_tmp_test+table_pat_test[i-1]
      pat_start_test[i]=count_tmp_test
    }
  }
  count=0
  for (i in 1:n_pat) {
    if (i<n_pat) {
      index_temp=pat_start[i]:(pat_start[i+1]-1)
    } else {
      index_temp=pat_start[i]:n_train
    }
    if (!is.null(x_test) & !is.null(y_test)) {
      if (i<n_pat) {
        index_temp_test=pat_start_test[i]:(pat_start_test[i+1]-1)
      } else {
        index_temp_test=pat_start_test[i]:n_test
      }
    }
    for (j in ug[[i]]) {
      ak = which(!is.na(x5_train[pat_start[j],]) & !is.na(x5_train[pat_start[i],]))
      k=which((pat==i) & (pat_ug==j))
      b1=c(b1, t(x_train[index_temp,ak,k])%*%y_train[index_temp,k]/length(index_temp))
      A1=rbind(A1, t(x_train[index_temp,ak,k])%*%x_train[index_temp,,k]/length(index_temp))
      count_temp=sum(!is.na(x5_train[pat_start[j],]) & !is.na(x5_train[pat_start[i],]))
      W[(count+1):(count+count_temp),(count+1):(count+count_temp)] = 
        (n_train/length(index_temp))^2*t(x_train[index_temp,ak,k])%*%x_train[index_temp,ak,k]
      if (!is.null(x_test) & !is.null(y_test)) {
        W_test[(count+1):(count+count_temp),(count+1):(count+count_temp)] = 
          (n_test/length(index_temp_test))^2*t(x_test[index_temp_test,ak,k])%*%x_test[index_temp_test,ak,k]
      }
      count=count+count_temp
    }
  }
  # dim(A1)
  # [1] 328  40
  
  nlam2=length(lambda2)
  v=array(0, dim=c(nlam2, p, N_g2))
  Tvar=matrix(0,nlam2,p)
  Pvar=matrix(0,nlam2,p)
  v_status=matrix(TRUE,nlam2,p)
  bound_status=matrix(TRUE,nlam2,p)
  
  W1 <- as(W, "dgCMatrix")
  #Dmat <- W
  #Dmat <- nearPD(W)$mat
  dvec <- rep(0,N_g2)
  #Amat <- rbind(-t(A1), t(A1))
  #Emat <- rbind(diag(1,p,p), -diag(1,p,p))
  Amat <- t(A1)
  Amat1 <- as(Amat, "dgCMatrix")
  Emat <- diag(1,p,p)
  
  # start.time = proc.time()
  # for (j in 1:p) {
  #   for (l in 1:nlam2) {
  #     #results=solve.QP(Dmat,dvec,Amat,bvec=Emat[,j]-lambda2[l])
  #     #v[l,j,]=results$solution
  #     mycop = cop(f=quadfun(Dmat), 
  #                 lc=lincon(A=Amat, d=Emat[,j]+lambda2[l], 
  #                           dir=rep(">=",nrow(Amat)), name=1:nrow(Amat)))
  #     #oldw <- getOption("warn")      ### do not print warning
  #     #options(warn = -1)
  #     results <- solvecop(mycop, quiet=TRUE, trace=FALSE)
  #     #options(warn = oldw)
  #     Evaluation <- validate(mycop, results, quiet=TRUE)
  #     v[l,j,]=results$x
  #     if (sum(Evaluation$summary[-1,5])!=nrow(Amat)) {
  #       v_status[l,j]=FALSE
  #     }
  #     if (!is.null(x_test) & !is.null(y_test)) {
  #       Pvar[l,j]=v[l,j,]%*%W_test%*%v[l,j,]
  #     }
  #   }
  # }
  # time1 = proc.time() - start.time
  ########## Parallel
  #start.time = proc.time()
  if (parallel) {
    results_j=foreach (j = debiased, .packages="osqp") %dopar% {
      v_j=matrix(0,nlam2,N_g2)
      Tvar_j=rep(0,nlam2)
      Pvar_j=rep(0,nlam2)
      v_status_j=rep(TRUE,nlam2)
      bound_status_j=rep(TRUE,nlam2)
      for (l in 1:nlam2) {
        #results=solve.QP(Dmat,dvec,t(Amat),bvec=-Emat[,j]-lambda2[l])
        #v[l,j,]=results$solution
        # mycop = cop(f=quadfun(Dmat), 
        #             lc=lincon(A=Amat, d=Emat[,j]+lambda2[l], 
        #                       dir=rep(">=",nrow(Amat)), name=1:nrow(Amat)))
        # #oldw <- getOption("warn")      ### do not print warning
        # #options(warn = -1)
        # results <- solvecop(mycop, quiet=TRUE, trace=FALSE)
        # #options(warn = oldw)
        # Evaluation <- validate(mycop, results, quiet=TRUE)
        # v_j[l,]=results$x
        # bound=bound0
        # b_step=0.1
        # b_count=1
        # while (b_count<=b_count_max) {
        #start.time = proc.time()
        # sv <- ipop(c=dvec, H=Dmat, A=Amat, b=Emat[,j]-rep(lambda2[l],p), 
        #            l=rep(-bound,N_g), u=rep(bound,N_g), r=2*rep(lambda2[l],p))
        #time1 = proc.time() - start.time
        settings <- osqpSettings(eps_abs =1e-10, eps_rel=1e-5, eps_prim_inf = 1e-10, max_iter=5000, verbose=FALSE)
        sv <- solve_osqp(P=W1, q=dvec, A=Amat1, l=Emat[,j]-rep(lambda2[l],p), 
                         u=Emat[,j] + rep(lambda2[l],p), settings)
        #v_j[l,]=sv@primal
        v_j[l,]=sv$x
        # if (max(abs(v_j[l,]))<bound) {
        #   break
        # } 
        #   opt_value=v_j[l,]%*%W%*%v_j[l,]
        #   if (b_count>1) {
        #     if (opt_value>=opt_value_pre) {
        #       break
        #     }
        #   }
        #   opt_value_pre=opt_value
        #   bound=bound+b_step
        #   b_count=b_count+1
        # }
        # if (b_count==(b_count_max+1)) {
        #   bound_status_j[l]=FALSE
        #   #print("Error: Bound is not large enough")
        # }
        #if (sum(Evaluation$summary[-1,5])!=nrow(Amat)) {
        # if (sv@how!="converged") {
        if (sv$info$status!='solved') {
          v_status_j[l]=FALSE
        }
        if (minimum) {
          Tvar_j[l]=v_j[l,]%*%W%*%v_j[l,]
        }
        if (!is.null(x_test) & !is.null(y_test)) {
          Pvar_j[l]=v_j[l,]%*%W_test%*%v_j[l,]
        }
      }
      return(list('v_j'=v_j, 'Pvar_j'=Pvar_j, 'Tvar_j'=Tvar_j, 
                  'v_status_j'=v_status_j, 'bound_status_j'=bound_status_j))
    }
  } else {
    results_j=foreach (j = debiased, .packages="osqp") %do% {
      v_j=matrix(0,nlam2,N_g2)
      Tvar_j=rep(0,nlam2)
      Pvar_j=rep(0,nlam2)
      v_status_j=rep(TRUE,nlam2)
      bound_status_j=rep(TRUE,nlam2)
      for (l in 1:nlam2) {
        #results=solve.QP(Dmat,dvec,t(Amat),bvec=-Emat[,j]-lambda2[l])
        #v[l,j,]=results$solution
        # mycop = cop(f=quadfun(Dmat), 
        #             lc=lincon(A=Amat, d=Emat[,j]+lambda2[l], 
        #                       dir=rep(">=",nrow(Amat)), name=1:nrow(Amat)))
        # #oldw <- getOption("warn")      ### do not print warning
        # #options(warn = -1)
        # results <- solvecop(mycop, quiet=TRUE, trace=FALSE)
        # #options(warn = oldw)
        # Evaluation <- validate(mycop, results, quiet=TRUE)
        # v_j[l,]=results$x
        # bound=bound0
        # b_step=0.1
        # b_count=1
        # while (b_count<=b_count_max) {
        #start.time = proc.time()
        # sv <- ipop(c=dvec, H=Dmat, A=Amat, b=Emat[,j]-rep(lambda2[l],p), 
        #            l=rep(-bound,N_g), u=rep(bound,N_g), r=2*rep(lambda2[l],p))
        #time1 = proc.time() - start.time
        settings <- osqpSettings(eps_abs =1e-10, eps_rel=1e-5, eps_prim_inf = 1e-10, max_iter=5000, verbose=FALSE)
        sv <- solve_osqp(P=W1, q=dvec, A=Amat1, l=Emat[,j]-rep(lambda2[l],p), 
                         u=Emat[,j] + rep(lambda2[l],p), settings)
        #v_j[l,]=sv@primal
        v_j[l,]=sv$x
        # if (max(abs(v_j[l,]))<bound) {
        #   break
        # } 
        #   opt_value=v_j[l,]%*%W%*%v_j[l,]
        #   if (b_count>1) {
        #     if (opt_value>=opt_value_pre) {
        #       break
        #     }
        #   }
        #   opt_value_pre=opt_value
        #   bound=bound+b_step
        #   b_count=b_count+1
        # }
        # if (b_count==(b_count_max+1)) {
        #   bound_status_j[l]=FALSE
        #   #print("Error: Bound is not large enough")
        # }
        #if (sum(Evaluation$summary[-1,5])!=nrow(Amat)) {
        # if (sv@how!="converged") {
        if (sv$info$status!='solved') {
          v_status_j[l]=FALSE
        }
        if (minimum) {
          Tvar_j[l]=v_j[l,]%*%W%*%v_j[l,]
        }
        if (!is.null(x_test) & !is.null(y_test)) {
          Pvar_j[l]=v_j[l,]%*%W_test%*%v_j[l,]
        }
      }
      return(list('v_j'=v_j, 'Pvar_j'=Pvar_j, 'Tvar_j'=Tvar_j, 
                  'v_status_j'=v_status_j, 'bound_status_j'=bound_status_j))
    }
  }
  for (j in 1:length(debiased)) {
    v[,debiased[j],]=results_j[[j]]$v_j
    Tvar[,debiased[j]]=results_j[[j]]$Tvar_j
    Pvar[,debiased[j]]=results_j[[j]]$Pvar_j
    v_status[,debiased[j]]=results_j[[j]]$v_status_j
    bound_status[,debiased[j]]=results_j[[j]]$bound_status_j
  }
  #time2 = proc.time() - start.time
  
  returnlist=list("v"=v, "Pvar"=Pvar, "Tvar"=Tvar, "lambda2"=lambda2, 
                  "v_status"=v_status, "bound_status"=bound_status, 
                  "A1"=A1, "b1"=b1, "N_g2"=N_g2)
  return(returnlist)
} 