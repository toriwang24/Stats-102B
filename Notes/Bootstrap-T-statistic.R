# Bootstrap Sampling Distribution of T-Statistic
source("bootsample-1.R")
source("bootstats.R")

set.seed(100)
n = 100 # Size of a single sample  
B = c(200,500,1000,2000,5000,10000) #,50000,100000,500000)
# Number of bootstrap replicates

#xseq = seq(-0.7+2,0.7+2,length=200)
quartz(width=12,height=8)
par(mfrow=c(2,3))


x = rnorm(n,3,1) # Data from N(3,1)
y = rnorm(n,1,1) # Data from N(1,1)
z=cbind(x,y)

for(j in 1:6){
  boot.samples = bootsampling(z,B[j])
  
  t.statistic = function(w){(mean(w[,1])-mean(w[,2]))/
            (sqrt((var(w[,1])/length(w[,1]))+(var(w[,2])/length(w[,2]))))}
  
  tstat = boot.stats(boot.samples,t.statistic)
  
  hist(tstat$theta, breaks=40,freq=F, xlim=c(-3.4,3.4), ylim=c(0,0.55),
       main=paste("Bootstrap Sampling Distribution: B =",B[j]))
  legend("topright",expression(tstat*" PDF"),lty=1,bty="n")
}

