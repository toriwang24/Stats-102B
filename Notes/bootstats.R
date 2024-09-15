# Calculating statistics of interest


boot.stats =  function(boot.sample,myfun,...){
  if(is.matrix(boot.sample)){
    theta = apply(boot.sample,2,myfun,...)
  } else {
    theta = apply(boot.sample,3,myfun,...)
  }
  if(is.matrix(theta)){
    return(list(theta=theta,cov=cov(t(theta))))
  } else{
    return(list(theta=theta,se=sd(theta)))
  } }
