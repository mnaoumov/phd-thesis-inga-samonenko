library(gtools)
library(maxLik)

#    xd is a standardized vector of length N
generateVector = function(randomGenerator, N) {
 xd = randomGenerator(N)
 N = length(xd)
 xd = (xd - mean(xd)) / sqrt(var(xd) * (N - 1) / N)
 xd
}
#    The following functions are defined to make the code more simple to read
getMean = function(b, f) {
 k = length(f) + 1
 N = length(b)
 f1 = c(f, 1 - sum(f))
 n = f1 * N
 means = numeric(k - 1)
 lastIndex = 0
 for (i in 1 : (k - 1)) {
  nextIndex = lastIndex + n[i]
  means[i] = sum(b[(lastIndex + 1) : nextIndex]) / N;
  lastIndex = nextIndex
 }
  means
}


kcgf= function(tau, xd, f){
#    This function computes an average cumulant generating function defined in the equation (4.2) and its derivatives
#    Inputs
#       tau    2*(k-1)-dimensional vector where tau = (t_0,t_1)
#       xd     the population a=1,a=2, ..., a=N
#       f      vector of ratios of group sizes (n_1,..  .,n_{k-1})
#    Outputs
#       kcgf$cgf      k(t), as defined in (4.2), where t is (k-1)-dimensional vector
#       kcgf$dcgf     (k-1)-dimensional vector k'(t)
#       kcgf$ddcgf    (k-1) by (k-1) matrix of k''(t)
 taum=matrix(tau,length(f),2)
 nn = length(xd)
 xsd = sqrt((var(xd) * (nn - 1))/nn)
 xdb = mean(xd)
 x = (xd - xdb)/xsd
 f1=1-sum(f)
 nn1=rep(1,nn)
 x2=rbind(nn1,x)
#   k(t), as defined in (4.2)
 kappa = mean(log(f1 + f%*%exp(taum%*%x2)))
 h=t(t(f * exp(taum%*%x2))/c(f1 + f %*% exp(taum%*%x2)))
 hx=t(x*t(f * exp(taum%*%x2))/c(f1 + f %*% exp(taum%*%x2)))
 hhx=rbind(h,hx)
 kappad = c(h%*%t(x2)/nn)
 kappad2= (cbind(rbind(diag(c(h%*%nn1)),diag(c(h%*%x))),rbind(diag(c(h%*%x)),diag(c(hx%*%x))))-hhx%*%t(hhx))/nn
 list(cgf = kappa, dcgf = kappad, ddcgf = kappad2)
}

kLam=function(x1, xd, f){
#    This function computes the test statistic defined in the equation (4.3)
#    This function solves  K'(t)=x for a given x=(f,x1) using Newton-Raphson method
#    Input
#        x1     (k-1)-dimension vector
#      xd     the population a=1,a=2, ..., a=N
#        f      vector of ratios of group sizes (n_1,..  .,n_{k-1})
#    Outputs
#        kLam$lam   value of the test statistic definned in (4.3) at x1
#        kLam$t     t(x), (k-1)-dimension vector, t(x)%*%x-K(t(x))=sup{t.x-K(t): t}
#        kLam$kdd   2*(k-1) by 2*(k-1) matrix K"(t(x))
#        kLam$err   calculate the accuracy of approximation using Newton-Raphson menthod
#        kLam$checkNA    check of convergance
 fx=c(f,x1)
 checkNA=0
 lf=length(f)
 tt=rep(0,2*lf)
 for (i in 1:10)
  {
   ktt=kcgf(tt,xd,f)
   if (!is.na(det(ktt$ddcgf)) && det(ktt$ddcgf)>.00000001)
    {tt=tt+solve(ktt$ddcgf)%*%(fx-ktt$dcgf)}
   else {{checkNA=1}&{break}}#t0=t0+solve(ktt$ddcgf)%*%(fx-ktt$dcgf)
  }
   ktt=kcgf(tt,xd,f)
   list(lam=-ktt$cgf+fx%*%tt,t=tt,kdd=ktt$ddcgf,err=fx-ktt$dcgf,checkNA=checkNA)
}



kHu=function(xd, f, v,M){
#     this function calculates SP tail probs for k-sample perm test based on Lam(x)
#    Input
#      xd     the population a=1,a=2, ..., a=N
#        f      vector of ratios of group sizes (n_1,..  .,n_{k-1})
#        v      parameter in P(Lam(X)>v)
#        M      number of MC replicates to approximate theintegral on the sphere
#    Outputs
#        kHu$tailpchi2    chi squared approximation
#        kHu$tailp        saddlepoint approximations as in Theorem 2 of Kolassa and Robinson (2011)
 nn = length(xd)
 n = f * nn
 xsd = sqrt((var(xd) * (nn - 1))/nn)
 xdb = mean(xd)
 x = (xd - xdb)/xsd
 f1=1-sum(f)
 lf=length(f)
 tailp1=1-pchisq(nn*v^2,lf)#first order approximation
 VV=kcgf(rep(0,2*lf),x,f)$ddcgf[(lf+1):(2*lf),(lf+1):(2*lf)]#xd
 svdV=svd(VV)
 Vh=svdV$u%*%diag(svdV$d^(1/2))%*%t(svdV$v)
 Gu=0
 sdelus=matrix(0,M,lf+4)
 Mm=M
 for (i in 1:M){
  s=rnorm(lf)
  s=s/sqrt(sum(s^2))
  s=c(Vh%*%s)
  r=v
   for (j in 1:6){
     kl=kLam(r*s,x,f)
     if (kl$checkNA == 1)
     break;
     r=r+(v-sqrt(2*kl$lam))*sqrt(2*kl$lam)/sum(s*kl$t[(lf+1):(2*lf)])
    }
 kl=kLam(r*s,x,f)#xd
 if (kl$checkNA==1){delus=NA}else {delus= f1*prod(f)*r^2/(v*sum(s*kl$t[(lf+1):(2*lf)])*(det(kl$kdd))^.5)}
 #Gu=ifelse(abs(v-sqrt(2*kl$lam))>.00001,Gu,Gu+delus)
 if (kl$checkNA==0) Gu=Gu+delus
 if (kl$checkNA==1) Mm=Mm-1
 sdelus[i,]=c(delus,r*s,v-sqrt(2*kl$lam),sum(s*kl$t[(lf+1):(2*lf)]),det(kcgf(kl$t,xd,f)$ddcgf))
 }
 cn=(nn^(lf/2))/((2)^(lf/2-1)*gamma(lf/2))
 tailpLR=ifelse(Mm==0,tailp1,tailp1+nn^(-1)*cn*v^(lf-2)*exp(-nn*v^2/2)*(Gu/Mm-1))
 tailpBN=ifelse(Mm==0,tailp1,1-pchisq(nn*(v-log(Gu/Mm)/(nn*v))^2,lf))
 list(tailpchi2=tailp1,tailp=c(tailpLR,tailpBN))
}


pvalcal=function(xd,f,M=30,MC,vv=0){
#   this function computes MC number of replicates and is used in PowerForMu function
 nn=length(xd)
 n=f*nn
 xsd = sqrt((var(xd) * (nn - 1))/nn)
 xdb = mean(xd)
 x = (xd - xdb)/xsd
 f1=1-sum(f)
 lf=length(f)
 N=c(0,cumsum(n))
 xb=orderedMeans=rep(0,lf)
 for(i in 1:lf){xb[i]=sum(x[(N[i]+1):N[i+1]])/nn}
 vl=kLam(xb,xd,f)
 v=ifelse(vv==0,sqrt(2*vl$lam),vv)
 out=ifelse(M==0,1-pchisq(nn*v^2,lf),kHu(xd,f,v,M)$tailp)
 MCpv=rep(0,MC)
 MCsqpv=rep(0,MC)
 MCm=MC
 for (k in 1:MC){
  xs=sample(x,nn)
  xsb=rep(0,lf)
  for(i in 1:lf){xsb[i]=sum(xs[(N[i]+1):N[i+1]])/nn }
  MCpv[k]=kLam(xsb,xd,f)$lam
  xsbk=c(xsb,-sum(xsb))
  MCsqpv[k]=nn*xsbk%*%diag(1/c(f,f1))%*%xsbk
 }
 list(SPpv=out,MCpv=MCpv,MCsqpv=MCsqpv,xb=xb)
}

PowerForMu = function(mu, distr, f, Z, M) {
#    PowerForMu gives the power of the F test and Lambda test, using approximations as in Theorem 1 #of K&R
#        mu     k-dimensional vector parameter
#        distr  distribution function, e.g. rnorm, rexp
#        Z      integer number of random permutations to obtain one p-value
#        M      integer number of p-values for each power
 k = length(f) + 1
 N = length(mu)
 n = f*N

 pFValues = numeric(Z)
 pSPLRValues = numeric(Z)
 pSPBNValues = numeric(Z)

 for (i in 1 : Z) {
  a = generateVector(distr, N)
  a = a + mu
  a = (a - mean(a)) / sqrt(var(a) * (N - 1) / N)
  Xbar = getMean(a, f)
  Fbar = c(Xbar, -sum(Xbar))
  LVal = kLam(Xbar, a, f)
  if (LVal$checkNA == 1) {
   pFValues[i] = NA
   pSPLRValues[i] = NA
   pSPBNValues[i] = NA
   next
  }
  u = sqrt(2*LVal$lam)
  FVal = c(N*Fbar%*%diag(1/c(n,N-sum(n)))%*%Fbar)
  FValues = pvalcal(a,f,MC=M)$MCsqpv
  pFValues[i] = mean(FValues > FVal)
  tt2 = kHu(a, f, u, 40)$tailp
  pSPLRValues[i] = tt2[1]
  pSPBNValues[i] = tt2[2]
  }
 PowerF = mean(pFValues[!is.na(pFValues)] < .05)
 PowerSPLR = mean(pSPLRValues[!is.na(pSPLRValues)] < .05)
 PowerSPBN = mean(pSPBNValues[!is.na(pSPBNValues)] < .05)
 list (PowerF = PowerF, PowerSPLR = PowerSPLR,PowerSPBN = PowerSPBN)
}

#    the following are functions which generate different distributions and can be input in #PowerForMu function
rexp2 = function(N) rexp(N) ^ 2
rgamma5 = function(N) rgamma(N, 5)
rgamma05 = function(N) rgamma(N, 0.5)

