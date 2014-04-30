#    The following function gives a standardized l by N matrix 
generateMatrix = function(distr, l, N) {
 A = matrix(distr(l * N), l, N)
 A = A - apply(A, 1, mean)
 B = A %*% t(A)
 s = svd(B)
 sqrtB = s$u %*% diag(s$d ^ (1 / 2)) %*% t(s$v)
 A = sqrt(N) * solve(sqrtB) %*% A
 A
}
#    The following functions are defined to make the code more simple to read
getMean = function(xd, f) {
 k = length(f) + 1
 l = dim(xd)[1]
 N = dim(xd)[2]
 f1 = c(f, 1 - sum(f))
 n = f1 * N
 means = matrix(0, l, k - 1)
 lastIndex = 0
 for (i in 1 : (k - 1)) {
  nextIndex = lastIndex + n[i]
  means[,i] = apply(xd[,((lastIndex + 1) : nextIndex), drop = FALSE], 1, sum) / N;
  lastIndex = nextIndex
 }
 means
}
getRandomMean = function(xd, f) {
 N = dim(xd)[2]
 getMean(xd[,sample.int(N)], f)
}



kcgf = function (tau, xd, f) {
#    This function computes an average cumulant generating function defined in the equation (5.2) and its derivatives
#    Inputs
#       tau    (l+1)*(k-1)-dimensional vector where tau = (t_0,t_1)
#       xd     the population of l-dimensional vectors a=1,a=2, ..., a=N
#       f      vector of ratios of group sizes (n_1,..  .,n_{k-1})
#    Outputs
#       kcgf$cgf      k(t_0,t_1), as defined in (5.2)
#       kcgf$dcgf     k'(t_0,t_1)
#       kcgf$ddcgf    k''(t_0,t_1)
 l = dim(xd)[1]
 N = dim(xd)[2]
 lf = length(f)
 f1 = 1 - sum(f)
 taum=matrix(tau,lf, l+1)
 nn1=rep(1,N)
 x2 = rbind(nn1, xd) #x2 = rbind(N1, xm)
 kt = mean(log(f1 + f %*% exp(taum %*% x2)))
 h = t(t(f * exp(taum %*% x2)) / c(f1 + f %*% exp(taum %*% x2)))
 dk = c(h%*%t(x2) / N)

 hx = list()
 hhx = NULL
 for (i in 1 : (l + 1)) {
  hx[[i]] = t(x2[i,] * t(h))
  hhx = rbind(hhx, hx[[i]])
 }

 ddk = NULL
 for (i in 1 : (l + 1)) {
  columnBlock = NULL
   for (j in 1 : (l + 1)) {
      columnBlock = rbind(columnBlock, diag(c(hx[[i]] %*% x2[j,])))
    }
  ddk = cbind(ddk, columnBlock)
 }

 ddk = (ddk - hhx %*% t(hhx)) / N
 list(cgf = kt, dcgf = dk, ddcgf = ddk)
}




kLam = function(x1,xd,f) {
#    This function computes the test statistic defined in the equation (4.3)
#    This function solves  K'(t)=x for a given x=(f,x1) using Newton-Raphson method
#    Input
#        x1     (k-1)-dimension vector
#        xd     the population of l-dimensional vectors a=1,a=2, ..., a=N
#        f      vector of ratios of group sizes (n_1,..  .,n_{k-1})
#    Outputs
#        kLam$lam   value of the test statistic definned in (4.3) at x1
#        kLam$t     t(x), (k-1)-dimension vector, t(x)%*%x-K(t(x))=sup{t.x-K(t): t}
#        kLam$kdd   (l+1)*(k-1) by (l+1)*(k-1) matrix K"(t(x))
#        kLam$err   calculate the accuracy of approximation using Newton-Raphson menthod
#        kLam$checkNA    check of convergance
  k = length(f) + 1
  fx = c(f, x1)
  checkNA = 0
  tt = rep(0, length(fx))
  for (i in 1 : 10) {
    ktt = kcgf(tt, xd, f)
    if (!is.na(det(ktt$ddcgf)) && det(ktt$ddcgf) > .00000001) {
      tt = tt + solve(ktt$ddcgf) %*% (fx - ktt$dcgf)
    } else {
      checkNA = 1
      break
    }
  }
  ktt = kcgf(tt, xd, f)
  list(lam = -ktt$cgf + fx %*% tt, t = tt, kdd = ktt$ddcgf, err = fx - ktt$dcgf, checkNA = checkNA)
}


kHu = function(xd, f, v, M) {
#     this function calculates SP tail probs for k-sample perm test based on Lam(x)
#    Input
#        xd     the population of l-dimensional vectors a=1,a=2, ..., a=N
#        f      vector of ratios of group sizes (n_1,..  .,n_{k-1})
#        v      parameter in P(Lam(X)>v)
#        M      number of MC replicates to approximate theintegral on the sphere
#    Outputs
#        kHu$tailpchi2    chi squared approximation
#        kHu$tailpLR      saddlepoint approximations in form of Lugananni and Rice
#        kHu$tailpBN      saddlepoint approximations in form of Barndoff and Nielson
 l = dim(xd)[1]
 N = dim(xd)[2]
 lf = length(f)
 k = lf + 1
 tailp1 = 1 - pchisq(N * v ^ 2, lf*l) #first order approximation
#   Vlr = \kappa_{00}'' = lower right corner of \kappa''(0)
 V = kcgf(rep(0, lf * (l + 1)), xd, f)$ddcgf
 Vlr = V[(lf + 1) : (lf * (l + 1)), (lf + 1) : (lf * (l + 1))]
 Vul = V[1 : lf, 1 : lf]
#   Vh = [\kappa_{00}'']^{1/2}
 svdV = svd(Vlr)
 Vh = svdV$u %*% diag(svdV$d ^ (1 / 2)) %*% t(svdV$v)
 cn = (N ^ (l * lf / 2)) / (2 ^ (l * lf / 2 - 1) * gamma(l * lf / 2))
#   Integral on sphere Gu = \int_{S_{l * lf}\delta(u,s)ds}
 Gu = 0
 Mm = M
 for (i in 1 : M) {
  s = rnorm(l * lf)
  s = s / sqrt(sum(s ^ 2))
  st = c(Vh %*% s)
#   r - solution of \Lambda(r [\kappa_{00}'']^{1/2} s) = v
  r = v
  for (j in 1 : 50) {
   kl = kLam(r * st, xd, f)
   if (kl$checkNA == 1)
         break;
   r = r - (kl$lam - v^2/2) / c(kl$t[(lf + 1) : (lf * (l + 1))] %*% st)
   }
  x = r * st
#   that = \hat t = (\hat t_0, \hat t_1): \kappa(\hat t_0, \hat t_1) = (p, x)
  that = kLam(x, xd, f)
  if (that$checkNA == 1) {
   Mm = Mm - 1
     next;
   }
#   t0hat = \hat t_0
  t0hat = that$t[1 : lf]
#   t1hat = \hat t_1
  t1hat = that$t[(lf + 1) : (lf * (l + 1))]
  delta = (det(Vul) ^ 0.5*
  det(kcgf(that$t, xd, f)$ddcgf) ^ (-0.5) *det(Vlr) ^ .5 *r ^ (l * lf - 1)) /
     (v ^ (l * lf - 2) * abs(t(s) %*% Vh %*% t1hat))
  Gu = Gu + delta
 }
 Gu = Gu / Mm
 if (Mm == 0) {
  tailpLR = tailp1
  tailpBN = tailp1
  } 
 else {
  tailpLR = tailp1 + N ^ (-1) * cn * v ^ (l * lf - 2) * exp(-N * v ^ 2 /2) * (Gu - 1)
  tailpBN = 1 - pchisq(N * (v - log(Gu) / (N * v)) ^ 2, lf*l)
  }
 list(tailp1 = tailp1, tailpLR = tailpLR, tailpBN = tailpBN)
}

getV=function(xx,f){
#   getV computes the test statistic V
 k = length(f) + 1
 l = dim(xx)[1]
 N = dim(xx)[2]
 f1 = c(f, 1 - sum(f))
 n = f1 * N

 S = matrix(0,l,k)
 S[,-k] = getMean(xx,f) 
 S[,k] = -apply(S[,-k],1,sum)

 Xbar = apply(xx,1,sum)/N
 T=xx%*%t(xx)-N*outer(Xbar,Xbar,"*")
 B=S%*%diag(1/n)%*%t(S)-N*outer(Xbar,Xbar,"*")
 W=T-B
 V=(N-k)*sum(diag(solve(W)%*%B))
 list(B=B,W=W,V=V)
}


PValuesForMu = function(mu, distr, l, N, f, Z, V) {
#   this function computes p-values for each statistic
#   Input
#     mu      l*k matrix  
#     distr   distribution function, e.g. rnorm, rexp
#     Z       integer number of p-values
#     V       integer number of replicates for each p-values of V test
 k = length(f) + 1
 f1 = c(f, 1 - sum(f))
 n = f1 * N

   addMu=function(xd,mu){
     M=matrix(0,l,N)
     lastIndex = 0

     for (i in 1:k){
       nextIndex = lastIndex + n[i]
       M[,((lastIndex + 1) : nextIndex)] = xd[,((lastIndex + 1) : nextIndex)]+mu[,i]
       lastIndex = nextIndex
     }
    M
   }

 pVValues = NULL
 pLRValues = NULL
 pBNValues = NULL
 for (i in 1 : Z) {
  xd = generateMatrix(distr, l, N)
  xd = addMu(xd, mu)
  xd = xd - apply(xd, 1, mean)
  B = xd %*% t(xd)
  s = svd(B)
  sqrtB = s$u %*% diag(s$d ^ (1 / 2)) %*% t(s$v)
  xd = sqrt(N) * solve(sqrtB) %*%xd

  Mean = getMean(xd, f)/N
  LVal = kLam(t(Mean), xd, f)$lam
  u = sqrt(2*LVal)
####get observed V
  VVal=getV(xd,f)$V

  VValues=NULL
  for(j in 1:V) VValues[j]=getV(xd[,sample(1:N,N)],f)$V#This gets rid of randomMeans
  pVValues[i] = mean(VValues>VVal)
    tt = kHu(xd, f, u, M=10)
    pLRValues[i] = c(tt$tailpLR)
    pBNValues[i] = c(tt$tailpBN)
 }
list(pVValues=pVValues,pLRValues=pLRValues,pBNValues=pBNValues)
}

outpow0=PValuesForMu(mu, distr, l, N, f, Z, V)
PowerV=mean(outpow0$pVValues<.05)
PowerLR=mean(outpow0$pLRValues<.05)
PowerBN=mean(outpow0$pBNValues<.05)

rbind(PowerV=PowerV,PowerLR=PowerLR,PowerBN=PowerBN)

#    the following are functions which generate different distributions and can be input in #PowerForMu function
rexp2 = function(N) rexp(N) ^ 2
rgamma5 = function(N) rgamma(N, 5)
rgamma05 = function(N) rgamma(N, 0.5)

