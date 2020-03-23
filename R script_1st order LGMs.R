library(lavaan)

#mean and covariance matrix of composite scores computed as summed item scores
sigma_sum=matrix(c(21.60,4.704,5.488,6.272,4.704,25.72,7.84,9.408,5.488,7.84,26.52,12.544,6.272,9.408,12.544,30.44),4,4)
mu_sum=c(5.2,8.0,10.80,13.6)

names(mu_sum)=c("x1","x2","x3","x4")
rownames(sigma_sum)=colnames(sigma_sum)=c("x1","x2","x3","x4")

#mean and covariance matrix of composite scores computed as mean item scores
sigma_mean=matrix(c(1.35,.294,.343,.392,.294,1.6075,.490,.588,.343,.490,1.6575,.784,.392,.588,.784,1.9025),4,4)
mu_mean=c(1.3,2.0,2.7,3.4)

names(mu_mean)=c("x1","x2","x3","x4")
rownames(sigma_mean)=colnames(sigma_mean)=c("x1","x2","x3","x4")


#specify the 1st-order LGM

LGM="
level =~ 1*x1 + 1*x2 + 1*x3 + 1*x4
slope =~ 0*x1 + 1*x2 + 2*x3 + 3*x4
level ~~ NA*level
slope ~~ NA*slope
level ~~ NA*slope

x1 ~~ NA*x1
x2 ~~ NA*x2
x3 ~~ NA*x3
x4 ~~ NA*x4

x1 ~0*1
x2 ~0*1
x3 ~0*1
x4 ~0*1

level ~ NA*1
slope ~ NA*1"

#This is to fit the 1st-order LGM using composite scores computed as summed item scores
sum.fit=sem(LGM,sample.mean=mu_sum,sample.cov=sigma_sum,sample.nobs=1000000)
summary(sum.fit,fit.measures=TRUE,rsquare=TRUE,modindices=FALSE,standardized=TRUE)

#This is to fit the 1st-order LGM using composite scores computed as mean item scores
mean.fit=sem(LGM,sample.mean=mu_mean,sample.cov=sigma_mean,sample.nobs=1000000)
summary(mean.fit,fit.measures=TRUE,rsquare=TRUE,modindices=FALSE,standardized=TRUE)