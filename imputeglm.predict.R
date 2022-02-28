library(MASS)
library(glmnetcr)
imputeglm.predict <- function (X, ind_y, ind_x=-ind_y, miss, newdata, family="gaussian") {
  ny=length(ind_y)
  nx=length(X[1,ind_x])+1
  B=matrix(0,nx,ny)
  PRED=matrix(0,dim(newdata)[1],ny)
  for (l in 1:ny) {
    ind_obs=!miss[,ind_y[l]]
    #ind_tmp=ind_y
    #ind_tmp2=11:30
    
    x.t=X[ind_obs,ind_x]
    y.t=X[ind_obs,ind_y[l]]
    x.train=x.t[apply(miss[ind_obs,ind_x],1,sum)==0,]
    y.train=y.t[apply(miss[ind_obs,ind_x],1,sum)==0]
    #cvfit <- cv.glmnet(x.train,y.train, nfolds = 3)
    #fit <- cvfit$glmnet.fit
    #coeff=coef(fit, s=cvfit$lambda.1se)
    #y.test=coeff[1]+x.test%*%coeff[2:length(coeff)]
    data0=data.frame(y.train, x.train)
    newx=as.data.frame(newdata)
    colnames(newx)=colnames(data0)[-1]
    if (family=="ordinal") {
      data0[,1]=as.factor(data0[,1])
      if (dim(x.train)[1]>dim(x.train)[2]) {
        fit=polr(y.train~., data = data0)
        coeff=c(0,fit$coef)                       ##### The intercept may be wrong
        pred=as.numeric(as.character(predict(fit, type = 'class', newdata=newx)))
      } else {
        fit <- glmnetcr(data0[,-1], data0[,1], maxit=500)
        select = select.glmnetcr(fit)
        coeff=fit$beta[1:nx,select]               ##### The intercept may be wrong
        pred=as.numeric(fitted(fit, newx=newx, s=select)$class)
      }
    } else if (family=="binomial") {
      data0[,1]=as.factor(data0[,1])
      if (dim(x.train)[1]>dim(x.train)[2]) {
        fit=glm(y.train~., family=family, data=data0)
        coeff=fit$coef
        prob=predict(fit,newdata=newx, type = 'response')
        pred=rep(as.numeric(levels(data0$y.train)[1]),dim(newdata)[1])
        pred[prob>0.5]=as.numeric(levels(data0$y.train)[2])
      } else {
        fit=glmnet(x.train, y.train, family=family)
        k <- fit$df
        n <- fit$nobs
        select=which.min(log(n)*k+deviance(fit))
        # cvfit <- cv.glmnet(x.train, y.train, nfolds=3, family=family)
        # fit <- cvfit$glmnet.fit
        coeff=t(as.numeric(coef(fit, s=select)))
        pred=as.numeric(predict(fit, newx=newdata, s=select, type = 'class'))
      }
    } else {
      if (dim(x.train)[1]>dim(x.train)[2]) {
        fit=glm(y.train~., family=family, data=data0)
        coeff=fit$coef
        pred=predict(fit,newdata=newx)
      } else {
        cvfit <- cv.glmnet(x.train, y.train, nfolds=3, family=family)
        fit <- cvfit$glmnet.fit
        coeff=t(as.numeric(coef(fit, s=cvfit$lambda.min)))
        pred=predict(fit, newx=newdata, s=cvfit$lambda.min)
      }
    }
    B[,l]=coeff
    PRED[,l]=pred
  }
  returnlist=list("B"=B, "PRED"=PRED)
  return(returnlist)
}