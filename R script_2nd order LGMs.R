library(MASS)
library(lavaan)

#This is to compute means and covariance matrix among items based on the selected parameter values

#intercepts of items
tau=c(.5,.5,.7,.7,.5,.5,.7,.7,.5,.5,.7,.7,.5,.5,.7,.7)
#intercepts of latent factor over time
a=c(0,0,0,0)
#mean of growth factors
k=c(1,1)

#first-order factor loadings matrix
lambda=matrix(0,16,4)
lambda[c(1:4),c(1)]=c(.6,.6,.8,.8)
lambda[c(5:8),c(2)]=c(.6,.6,.8,.8)
lambda[c(9:12),c(3)]=c(.6,.6,.8,.8)
lambda[c(13:16),c(4)]=c(.6,.6,.8,.8)

#item residual covariance matrix
res=c(.6,.6,.4,.4,.5,.5,.6,.6,.7,.7,.8,.8,.8,.8,.7,.7)
psi=diag(16)*res

#second-order factor loadings (for growth trajectory)
gamma=matrix(c(1,1,1,1,0,1,2,3),4,2)

#covariance matrix of growth factors
phi=matrix(c(.5,.1,.1,.1),2,2)

#covariance matrix of factor uniqueness
fecov=(c(2.0,2.2,1.7,1.5))
theta=diag(4)*fecov


#compute the means and covariance matrix among items
mu=tau+lambda%*%(a+gamma%*%k)
mu=t(mu)
fcov=gamma%*%phi%*%t(gamma)+theta
sigma=lambda%*%fcov%*%t(lambda)+psi


colnames(mu)=c("x11","x12","x13","x14","x21","x22","x23","x24","x31","x32","x33","x34","x41","x42","x43","x44")
rownames(sigma)=colnames(sigma)=c("x11","x12","x13","x14","x21","x22","x23","x24","x31","x32","x33","x34","x41","x42","x43","x44")



#This is to analyze 11 second-order LGMs.

#The Marker-1 approach: x.1 was chosen as the marker variable.
#The loading and the intercept of this item was fixed at 1 and 0, respectively.

Marker1="
eta1 =~1*x11+L2*x12+L3*x13+L4*x14
eta2 =~1*x21+L2*x22+L3*x23+L4*x24
eta3 =~1*x31+L2*x32+L3*x33+L4*x34
eta4 =~1*x41+L2*x42+L3*x43+L4*x44
level =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4
slope =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4

x11 ~~ NA*x11
x12 ~~ NA*x12
x13 ~~ NA*x13
x14 ~~ NA*x14
x21 ~~ NA*x21
x22 ~~ NA*x22
x23 ~~ NA*x23
x24 ~~ NA*x24
x31 ~~ NA*x31
x32 ~~ NA*x32
x33 ~~ NA*x33
x34 ~~ NA*x34
x41 ~~ NA*x41
x42 ~~ NA*x42
x43 ~~ NA*x43
x44 ~~ NA*x44


level ~~ NA*level
slope ~~ NA*slope
level ~~ NA*slope

eta1 ~~ NA*eta1
eta2 ~~ NA*eta2
eta3 ~~ NA*eta3
eta4 ~~ NA*eta4


x11 ~ 0*1
x12 ~ tau2*1
x13 ~ tau3*1
x14 ~ tau4*1
x21 ~ 0*1
x22 ~ tau2*1
x23 ~ tau3*1
x24 ~ tau4*1
x31 ~ 0*1
x32 ~ tau2*1
x33 ~ tau3*1
x34 ~ tau4*1
x41 ~ 0*1
x42 ~ tau2*1
x43 ~ tau3*1
x44 ~ tau4*1

eta1 ~0*1
eta2 ~0*1
eta3 ~0*1
eta4 ~0*1

level ~ NA*1
slope ~ NA*1"

Marker1.fit=sem(Marker1,sample.mean=mu,sample.cov=sigma,sample.nobs=1000000)
summary(Marker1.fit,fit.measures=TRUE,rsquare=TRUE,modindices=FALSE,standardized=TRUE)


#The Marker-2 approach: x.1 was chosen as the marker variable.
#The loading and the intercept of this item was fixed at 0.6 and 0, respectively.

Marker2="
eta1 =~0.6*x11+L2*x12+L3*x13+L4*x14
eta2 =~0.6*x21+L2*x22+L3*x23+L4*x24
eta3 =~0.6*x31+L2*x32+L3*x33+L4*x34
eta4 =~0.6*x41+L2*x42+L3*x43+L4*x44
level =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4
slope =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4

x11 ~~ NA*x11
x12 ~~ NA*x12
x13 ~~ NA*x13
x14 ~~ NA*x14
x21 ~~ NA*x21
x22 ~~ NA*x22
x23 ~~ NA*x23
x24 ~~ NA*x24
x31 ~~ NA*x31
x32 ~~ NA*x32
x33 ~~ NA*x33
x34 ~~ NA*x34
x41 ~~ NA*x41
x42 ~~ NA*x42
x43 ~~ NA*x43
x44 ~~ NA*x44


level ~~ NA*level
slope ~~ NA*slope
level ~~ NA*slope

eta1 ~~ NA*eta1
eta2 ~~ NA*eta2
eta3 ~~ NA*eta3
eta4 ~~ NA*eta4


x11 ~ 0*1
x12 ~ tau2*1
x13 ~ tau3*1
x14 ~ tau4*1
x21 ~ 0*1
x22 ~ tau2*1
x23 ~ tau3*1
x24 ~ tau4*1
x31 ~ 0*1
x32 ~ tau2*1
x33 ~ tau3*1
x34 ~ tau4*1
x41 ~ 0*1
x42 ~ tau2*1
x43 ~ tau3*1
x44 ~ tau4*1

eta1 ~0*1
eta2 ~0*1
eta3 ~0*1
eta4 ~0*1

level ~ NA*1
slope ~ NA*1"

Marker2.fit=sem(Marker2,sample.mean=mu,sample.cov=sigma,sample.nobs=1000000)
summary(Marker2.fit,fit.measures=TRUE,rsquare=TRUE,modindices=FALSE,standardized=TRUE)



#The Marker-3 approach: x.1 was chosen as the marker variable.
#The loading and the intercept of this item was fixed at 0.6 and 0.5, respectively.

Marker3="
eta1 =~0.6*x11+L2*x12+L3*x13+L4*x14
eta2 =~0.6*x21+L2*x22+L3*x23+L4*x24
eta3 =~0.6*x31+L2*x32+L3*x33+L4*x34
eta4 =~0.6*x41+L2*x42+L3*x43+L4*x44
level =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4
slope =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4

x11 ~~ NA*x11
x12 ~~ NA*x12
x13 ~~ NA*x13
x14 ~~ NA*x14
x21 ~~ NA*x21
x22 ~~ NA*x22
x23 ~~ NA*x23
x24 ~~ NA*x24
x31 ~~ NA*x31
x32 ~~ NA*x32
x33 ~~ NA*x33
x34 ~~ NA*x34
x41 ~~ NA*x41
x42 ~~ NA*x42
x43 ~~ NA*x43
x44 ~~ NA*x44


level ~~ NA*level
slope ~~ NA*slope
level ~~ NA*slope

eta1 ~~ NA*eta1
eta2 ~~ NA*eta2
eta3 ~~ NA*eta3
eta4 ~~ NA*eta4


x11 ~ 0.5*1
x12 ~ tau2*1
x13 ~ tau3*1
x14 ~ tau4*1
x21 ~ 0.5*1
x22 ~ tau2*1
x23 ~ tau3*1
x24 ~ tau4*1
x31 ~ 0.5*1
x32 ~ tau2*1
x33 ~ tau3*1
x34 ~ tau4*1
x41 ~ 0.5*1
x42 ~ tau2*1
x43 ~ tau3*1
x44 ~ tau4*1

eta1 ~0*1
eta2 ~0*1
eta3 ~0*1
eta4 ~0*1

level ~ NA*1
slope ~ NA*1"

Marker3.fit=sem(Marker3,sample.mean=mu,sample.cov=sigma,sample.nobs=1000000)
summary(Marker3.fit,fit.measures=TRUE,rsquare=TRUE,modindices=FALSE,standardized=TRUE)


#The Marker-4 approach: x.4 was chosen as the marker variable.
#The loading and the intercept of this item was fixed at 1 and 0, respectively.

Marker4="
eta1 =~1*x14+L1*x11+L2*x12+L3*x13
eta2 =~1*x24+L1*x21+L2*x22+L3*x23
eta3 =~1*x34+L1*x31+L2*x32+L3*x33
eta4 =~1*x44+L1*x41+L2*x42+L3*x43
level =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4
slope =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4

x11 ~~ NA*x11
x12 ~~ NA*x12
x13 ~~ NA*x13
x14 ~~ NA*x14
x21 ~~ NA*x21
x22 ~~ NA*x22
x23 ~~ NA*x23
x24 ~~ NA*x24
x31 ~~ NA*x31
x32 ~~ NA*x32
x33 ~~ NA*x33
x34 ~~ NA*x34
x41 ~~ NA*x41
x42 ~~ NA*x42
x43 ~~ NA*x43
x44 ~~ NA*x44


level ~~ NA*level
slope ~~ NA*slope
level ~~ NA*slope

eta1 ~~ NA*eta1
eta2 ~~ NA*eta2
eta3 ~~ NA*eta3
eta4 ~~ NA*eta4

x11 ~ tau1*1
x12 ~ tau2*1
x13 ~ tau3*1
x14 ~ 0*1
x21 ~ tau1*1
x22 ~ tau2*1
x23 ~ tau3*1
x24 ~ 0*1
x31 ~ tau1*1
x32 ~ tau2*1
x33 ~ tau3*1
x34 ~ 0*1
x41 ~ tau1*1
x42 ~ tau2*1
x43 ~ tau3*1
x44 ~ 0*1

eta1 ~0*1
eta2 ~0*1
eta3 ~0*1
eta4 ~0*1

level ~ NA*1
slope ~ NA*1"

Marker4.fit=sem(Marker4,sample.mean=mu,sample.cov=sigma,sample.nobs=1000000)
summary(Marker4.fit,fit.measures=TRUE,rsquare=TRUE,modindices=FALSE,standardized=TRUE)



#The Marker-5 approach: x.4 was chosen as the marker variable.
#The loading and the intercept of this item was fixed at 0.8 and 0, respectively.

Marker5="
eta1 =~0.8*x14+L1*x11+L2*x12+L3*x13
eta2 =~0.8*x24+L1*x21+L2*x22+L3*x23
eta3 =~0.8*x34+L1*x31+L2*x32+L3*x33
eta4 =~0.8*x44+L1*x41+L2*x42+L3*x43
level =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4
slope =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4

x11 ~~ NA*x11
x12 ~~ NA*x12
x13 ~~ NA*x13
x14 ~~ NA*x14
x21 ~~ NA*x21
x22 ~~ NA*x22
x23 ~~ NA*x23
x24 ~~ NA*x24
x31 ~~ NA*x31
x32 ~~ NA*x32
x33 ~~ NA*x33
x34 ~~ NA*x34
x41 ~~ NA*x41
x42 ~~ NA*x42
x43 ~~ NA*x43
x44 ~~ NA*x44


level ~~ NA*level
slope ~~ NA*slope
level ~~ NA*slope

eta1 ~~ NA*eta1
eta2 ~~ NA*eta2
eta3 ~~ NA*eta3
eta4 ~~ NA*eta4


x11 ~ tau1*1
x12 ~ tau2*1
x13 ~ tau3*1
x14 ~ 0*1
x21 ~ tau1*1
x22 ~ tau2*1
x23 ~ tau3*1
x24 ~ 0*1
x31 ~ tau1*1
x32 ~ tau2*1
x33 ~ tau3*1
x34 ~ 0*1
x41 ~ tau1*1
x42 ~ tau2*1
x43 ~ tau3*1
x44 ~ 0*1

eta1 ~0*1
eta2 ~0*1
eta3 ~0*1
eta4 ~0*1

level ~ NA*1
slope ~ NA*1"

Marker5.fit=sem(Marker5,sample.mean=mu,sample.cov=sigma,sample.nobs=1000000)
summary(Marker5.fit,fit.measures=TRUE,rsquare=TRUE,modindices=FALSE,standardized=TRUE)


#The Marker-6 approach: x.4 was chosen as the marker variable.
#The loading and the intercept of this item was fixed at 0.8 and 0.7, respectively.

Marker6="
eta1 =~0.8*x14+L1*x11+L2*x12+L3*x13
eta2 =~0.8*x24+L1*x21+L2*x22+L3*x23
eta3 =~0.8*x34+L1*x31+L2*x32+L3*x33
eta4 =~0.8*x44+L1*x41+L2*x42+L3*x43
level =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4
slope =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4

x11 ~~ NA*x11
x12 ~~ NA*x12
x13 ~~ NA*x13
x14 ~~ NA*x14
x21 ~~ NA*x21
x22 ~~ NA*x22
x23 ~~ NA*x23
x24 ~~ NA*x24
x31 ~~ NA*x31
x32 ~~ NA*x32
x33 ~~ NA*x33
x34 ~~ NA*x34
x41 ~~ NA*x41
x42 ~~ NA*x42
x43 ~~ NA*x43
x44 ~~ NA*x44

level ~~ NA*level
slope ~~ NA*slope
level ~~ NA*slope

eta1 ~~ NA*eta1
eta2 ~~ NA*eta2
eta3 ~~ NA*eta3
eta4 ~~ NA*eta4


x11 ~ tau1*1
x12 ~ tau2*1
x13 ~ tau3*1
x14 ~ 0.7*1
x21 ~ tau1*1
x22 ~ tau2*1
x23 ~ tau3*1
x24 ~ 0.7*1
x31 ~ tau1*1
x32 ~ tau2*1
x33 ~ tau3*1
x34 ~ 0.7*1
x41 ~ tau1*1
x42 ~ tau2*1
x43 ~ tau3*1
x44 ~ 0.7*1

eta1 ~0*1
eta2 ~0*1
eta3 ~0*1
eta4 ~0*1

level ~ NA*1
slope ~ NA*1"

Marker6.fit=sem(Marker6,sample.mean=mu,sample.cov=sigma,sample.nobs=1000000)
summary(Marker6.fit,fit.measures=TRUE,rsquare=TRUE,modindices=FALSE,standardized=TRUE)



#The Marker-7 approach: x.1 was chosen as the marker variable.
#The loading and the intercept of this item was fixed at 0.9487 and 1.1, respectively.
#Consequently, the latent factor at time 1 was standardized.

Marker7="
eta1 =~0.9487*x11+L2*x12+L3*x13+L4*x14
eta2 =~0.9487*x21+L2*x22+L3*x23+L4*x24
eta3 =~0.9487*x31+L2*x32+L3*x33+L4*x34
eta4 =~0.9487*x41+L2*x42+L3*x43+L4*x44
level =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4
slope =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4

x11 ~~ NA*x11
x12 ~~ NA*x12
x13 ~~ NA*x13
x14 ~~ NA*x14
x21 ~~ NA*x21
x22 ~~ NA*x22
x23 ~~ NA*x23
x24 ~~ NA*x24
x31 ~~ NA*x31
x32 ~~ NA*x32
x33 ~~ NA*x33
x34 ~~ NA*x34
x41 ~~ NA*x41
x42 ~~ NA*x42
x43 ~~ NA*x43
x44 ~~ NA*x44

level ~~ NA*level
slope ~~ NA*slope
level ~~ NA*slope

eta1 ~~ NA*eta1
eta2 ~~ NA*eta2
eta3 ~~ NA*eta3
eta4 ~~ NA*eta4

x11 ~ 1.1*1
x12 ~ tau2*1
x13 ~ tau3*1
x14 ~ tau4*1
x21 ~ 1.1*1
x22 ~ tau2*1
x23 ~ tau3*1
x24 ~ tau4*1
x31 ~ 1.1*1
x32 ~ tau2*1
x33 ~ tau3*1
x34 ~ tau4*1
x41 ~ 1.1*1
x42 ~ tau2*1
x43 ~ tau3*1
x44 ~ tau4*1

eta1 ~0*1
eta2 ~0*1
eta3 ~0*1
eta4 ~0*1

level ~ NA*1
slope ~ NA*1"

Marker7.fit=sem(Marker7,sample.mean=mu,sample.cov=sigma,sample.nobs=1000000)
summary(Marker7.fit,fit.measures=TRUE,rsquare=TRUE,modindices=FALSE,standardized=TRUE)


#The Marker-8 approach: x.4 was chosen as the marker variable.
#The loading and the intercept of this item was fixed at 1.265 and 1.5, respectively.
#Consequently, the latent factor at time 1 was standardized.

Marker8="
eta1 =~1.265*x14+L1*x11+L2*x12+L3*x13
eta2 =~1.265*x24+L1*x21+L2*x22+L3*x23
eta3 =~1.265*x34+L1*x31+L2*x32+L3*x33
eta4 =~1.265*x44+L1*x41+L2*x42+L3*x43
level =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4
slope =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4

x11 ~~ NA*x11
x12 ~~ NA*x12
x13 ~~ NA*x13
x14 ~~ NA*x14
x21 ~~ NA*x21
x22 ~~ NA*x22
x23 ~~ NA*x23
x24 ~~ NA*x24
x31 ~~ NA*x31
x32 ~~ NA*x32
x33 ~~ NA*x33
x34 ~~ NA*x34
x41 ~~ NA*x41
x42 ~~ NA*x42
x43 ~~ NA*x43
x44 ~~ NA*x44

level ~~ NA*level
slope ~~ NA*slope
level ~~ NA*slope

eta1 ~~ NA*eta1
eta2 ~~ NA*eta2
eta3 ~~ NA*eta3
eta4 ~~ NA*eta4


x11 ~ tau1*1
x12 ~ tau2*1
x13 ~ tau3*1
x14 ~ 1.5*1
x21 ~ tau1*1
x22 ~ tau2*1
x23 ~ tau3*1
x24 ~ 1.5*1
x31 ~ tau1*1
x32 ~ tau2*1
x33 ~ tau3*1
x34 ~ 1.5*1
x41 ~ tau1*1
x42 ~ tau2*1
x43 ~ tau3*1
x44 ~ 1.5*1

eta1 ~0*1
eta2 ~0*1
eta3 ~0*1
eta4 ~0*1

level ~ NA*1
slope ~ NA*1"

Marker8.fit=sem(Marker8,sample.mean=mu,sample.cov=sigma,sample.nobs=1000000)
summary(Marker8.fit,fit.measures=TRUE,rsquare=TRUE,modindices=FALSE,standardized=TRUE)



#The ECI-1 approach: the sum of loadings=4; the sum of item intercept=0

ECI1="
eta1 =~L1*x11+L2*x12+L3*x13+L4*x14
eta2 =~L1*x21+L2*x22+L3*x23+L4*x24
eta3 =~L1*x31+L2*x32+L3*x33+L4*x34
eta4 =~L1*x41+L2*x42+L3*x43+L4*x44
level =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4
slope =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4

L1 == 4-L2-L3-L4

x11 ~~ NA*x11
x12 ~~ NA*x12
x13 ~~ NA*x13
x14 ~~ NA*x14
x21 ~~ NA*x21
x22 ~~ NA*x22
x23 ~~ NA*x23
x24 ~~ NA*x24
x31 ~~ NA*x31
x32 ~~ NA*x32
x33 ~~ NA*x33
x34 ~~ NA*x34
x41 ~~ NA*x41
x42 ~~ NA*x42
x43 ~~ NA*x43
x44 ~~ NA*x44

eta1 ~~ NA*eta1
eta2 ~~ NA*eta2
eta3 ~~ NA*eta3
eta4 ~~ NA*eta4

level ~~ NA*level
slope ~~ NA*slope
level ~~ NA*slope

x11 ~ tau1*1
x12 ~ tau2*1
x13 ~ tau3*1
x14 ~ tau4*1
x21 ~ tau1*1
x22 ~ tau2*1
x23 ~ tau3*1
x24 ~ tau4*1
x31 ~ tau1*1
x32 ~ tau2*1
x33 ~ tau3*1
x34 ~ tau4*1
x41 ~ tau1*1
x42 ~ tau2*1
x43 ~ tau3*1
x44 ~ tau4*1

tau1 == 0-tau2-tau3-tau4


eta1 ~0*1
eta2 ~0*1
eta3 ~0*1
eta4 ~0*1

level ~ NA*1
slope ~ NA*1"
ECI1.fit=sem(ECI1,sample.mean=mu,sample.cov=sigma,sample.nobs=1000000,std.lv=TRUE)
summary(ECI1.fit,fit.measures=TRUE,rsquare=TRUE,modindices=FALSE,standardized=TRUE)


#The ECI-2 approach: the sum of loadings=1; the sum of item intercept=0
ECI2="
eta1 =~L1*x11+L2*x12+L3*x13+L4*x14
eta2 =~L1*x21+L2*x22+L3*x23+L4*x24
eta3 =~L1*x31+L2*x32+L3*x33+L4*x34
eta4 =~L1*x41+L2*x42+L3*x43+L4*x44
level =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4
slope =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4

L1 == 1-L2-L3-L4

x11 ~~ NA*x11
x12 ~~ NA*x12
x13 ~~ NA*x13
x14 ~~ NA*x14
x21 ~~ NA*x21
x22 ~~ NA*x22
x23 ~~ NA*x23
x24 ~~ NA*x24
x31 ~~ NA*x31
x32 ~~ NA*x32
x33 ~~ NA*x33
x34 ~~ NA*x34
x41 ~~ NA*x41
x42 ~~ NA*x42
x43 ~~ NA*x43
x44 ~~ NA*x44

eta1 ~~ NA*eta1
eta2 ~~ NA*eta2
eta3 ~~ NA*eta3
eta4 ~~ NA*eta4

level ~~ NA*level
slope ~~ NA*slope
level ~~ NA*slope

x11 ~ tau1*1
x12 ~ tau2*1
x13 ~ tau3*1
x14 ~ tau4*1
x21 ~ tau1*1
x22 ~ tau2*1
x23 ~ tau3*1
x24 ~ tau4*1
x31 ~ tau1*1
x32 ~ tau2*1
x33 ~ tau3*1
x34 ~ tau4*1
x41 ~ tau1*1
x42 ~ tau2*1
x43 ~ tau3*1
x44 ~ tau4*1

tau1 == 0-tau2-tau3-tau4


eta1 ~0*1
eta2 ~0*1
eta3 ~0*1
eta4 ~0*1

level ~ NA*1
slope ~ NA*1"

ECI2.fit=sem(ECI2,sample.mean=mu,sample.cov=sigma,sample.nobs=1000000,std.lv=TRUE)
summary(ECI2.fit,fit.measures=TRUE,rsquare=TRUE,modindices=FALSE,standardized=TRUE)



#The ECI-3 approach: the sum of loadings=4.427; the sum of item intercept=5.2

ECI3="
eta1 =~L1*x11+L2*x12+L3*x13+L4*x14
eta2 =~L1*x21+L2*x22+L3*x23+L4*x24
eta3 =~L1*x31+L2*x32+L3*x33+L4*x34
eta4 =~L1*x41+L2*x42+L3*x43+L4*x44
level =~ 1*eta1 + 1*eta2 + 1*eta3 + 1*eta4
slope =~ 0*eta1 + 1*eta2 + 2*eta3 + 3*eta4

L1 == 4.427-L2-L3-L4

x11 ~~ NA*x11
x12 ~~ NA*x12
x13 ~~ NA*x13
x14 ~~ NA*x14
x21 ~~ NA*x21
x22 ~~ NA*x22
x23 ~~ NA*x23
x24 ~~ NA*x24
x31 ~~ NA*x31
x32 ~~ NA*x32
x33 ~~ NA*x33
x34 ~~ NA*x34
x41 ~~ NA*x41
x42 ~~ NA*x42
x43 ~~ NA*x43
x44 ~~ NA*x44

eta1 ~~ NA*eta1
eta2 ~~ NA*eta2
eta3 ~~ NA*eta3
eta4 ~~ NA*eta4

level ~~ NA*level
slope ~~ NA*slope
level ~~ NA*slope

x11 ~ tau1*1
x12 ~ tau2*1
x13 ~ tau3*1
x14 ~ tau4*1
x21 ~ tau1*1
x22 ~ tau2*1
x23 ~ tau3*1
x24 ~ tau4*1
x31 ~ tau1*1
x32 ~ tau2*1
x33 ~ tau3*1
x34 ~ tau4*1
x41 ~ tau1*1
x42 ~ tau2*1
x43 ~ tau3*1
x44 ~ tau4*1

tau1 == 5.2-tau2-tau3-tau4


eta1 ~0*1
eta2 ~0*1
eta3 ~0*1
eta4 ~0*1

level ~ NA*1
slope ~ NA*1"
ECI3.fit=sem(ECI3,sample.mean=mu,sample.cov=sigma,sample.nobs=1000000,std.lv=TRUE)
summary(ECI3.fit,fit.measures=TRUE,rsquare=TRUE,modindices=FALSE,standardized=TRUE)



fit_Marker1 <- fitMeasures(Marker1.fit, c("chisq", "df", "pvalue", "cfi","tli", "rmsea"))
fit_Marker2 <- fitMeasures(Marker2.fit, c("chisq", "df", "pvalue", "cfi","tli", "rmsea"))
fit_Marker3 <- fitMeasures(Marker3.fit, c("chisq", "df", "pvalue", "cfi","tli", "rmsea"))
fit_Marker4 <- fitMeasures(Marker4.fit, c("chisq", "df", "pvalue", "cfi","tli", "rmsea"))
fit_Marker5 <- fitMeasures(Marker5.fit, c("chisq", "df", "pvalue", "cfi","tli", "rmsea"))
fit_Marker6 <- fitMeasures(Marker6.fit, c("chisq", "df", "pvalue", "cfi","tli", "rmsea"))
fit_Marker7 <- fitMeasures(Marker7.fit, c("chisq", "df", "pvalue", "cfi","tli", "rmsea"))
fit_Marker8 <- fitMeasures(Marker8.fit, c("chisq", "df", "pvalue", "cfi","tli", "rmsea"))

fit_ECI1 <- fitMeasures(ECI1.fit, c("chisq", "df", "pvalue", "cfi","tli", "rmsea"))
fit_ECI2 <- fitMeasures(ECI2.fit, c("chisq", "df", "pvalue", "cfi","tli", "rmsea"))
fit_ECI3 <- fitMeasures(ECI3.fit, c("chisq", "df", "pvalue", "cfi","tli", "rmsea"))

fit_Marker1
fit_Marker2
fit_Marker3
fit_Marker4
fit_Marker5
fit_Marker6
fit_Marker7
fit_Marker8
fit_ECI1
fit_ECI2
fit_ECI3
