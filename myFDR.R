myFDR <- function (test_stats, alpha) {
  p=length(test_stats)
  b_p=sqrt(2*log(p)-2*log(log(p)))
  test_stats_2=sort(abs(test_stats))
  test_stats_2[p]=Inf
  t_hat=NULL
  for (i in 1:p) {
    t=test_stats_2[i]
    tmp=p*(2-2*pnorm(t))/(max(c(1,sum(test_stats_2>=t))))
    if (tmp<=alpha) {
      tmp2=1-alpha*(max(c(1,sum(test_stats_2>=t))))/(2*p)
      t_hat=qnorm(tmp2)
      if (t_hat>t | t_hat<test_stats_2[i-1] | sum(test_stats_2>=t)!=sum(test_stats_2>=t_hat)) {
        print('FDR error')
      }
      break
    }
  }
  if (is.null(t_hat)) {
    print('FDR error')
  } else if (t_hat>b_p | t_hat<0) {
    t_hat=sqrt(2*log(p))
  }
  return(t_hat)
}