# used for permutations
library(gtools)

# used for maxNR
library(maxLik)

#    A is a standardized matrix from an assignment of k different treatments to each block, b is total number of blocks
generateRandomMatrix = function(randomGenerator, b, k) {
#    This function generates the matrix A
  A = matrix(randomGenerator(b * k), b, k)
  A = A - apply(A, 1, mean)
  A = A / sqrt(mean(A ^ 2))
}

#    The following functions are defined to make the code more simple to read
getb = function(A) dim(A)[1]
getk = function(A) dim(A)[2]
perm = function(A, i) {
  k = getk(A)
  permutations(k, k - 1, A[i,])
 }
expb = function(A, i, beta) exp(perm(A, i) %*% beta)
orderedMeans = function(A) {
  b = getb(A)
  k = getk(A)
  ASorted = matrix(numeric(), b, k)
  for(i in 1 : b) ASorted[i,] = sort(A[i,])
  colMeans(ASorted)
 }
getRandomMean = function(A) {
  b = getb(A)
  k = getk(A)
  B = array(dim = c(b, k))
  for(i in 1 : b) B[i,] = sample(A[i,], k)
  mean = apply(B, 2, mean)[1 : (k - 1)]
  mean
 }
addMu = function(A, mu) {
  B = A
  b = dim(B)[1]
  for (i in 1 : b) B[i,] = B[i,] + mu
  B/sqrt(mean(B^2))
 }
Sqr = function(M) {
 s = svd(M)
 s$u %*% diag(sqrt(s$d)) %*% t(s$v)
 }
UniformOnSphere = function(n) {
  x = rnorm(n)
  r = sqrt(sum(x ^ 2))
  x / r
 }
IntOnSphereMonteCarlo = function(h, d, M) {
  sum = 0
  pointsSkipped = 0
  for(i in 1 : M) {
    funcValue = h(UniformOnSphere(d))
    if (is.na(funcValue)) pointsSkipped = pointsSkipped + 1
    else sum = sum + funcValue
   }
  if (pointsSkipped == M) return (0)
  sum / (M - pointsSkipped)
 }
isPointInsideDomain = function(A, x, tolerance = 1) {
#   this function perform a check whether a given point x belonds to the admissible domain
 k = getk(A)
 y = x
 if (length(y) == k - 1)
 y[length(y) + 1] = -sum(y)
 stopifnot(length(y) == k)
  for (usedIndices in 1 : (k - 1)) {
   subVectors = combinations(k, usedIndices, y, set = FALSE)
   subVectorSums = apply(subVectors, 1, sum)
   sumMins = sum(orderedMeans(A)[1 : usedIndices])
   sumMaxs = sum(orderedMeans(A)[(k - usedIndices + 1) : k])
   if (!all(subVectorSums > sumMins * tolerance))
     return(FALSE);
   if (!all(subVectorSums < sumMaxs * tolerance))
     return(FALSE);
  }
 return(TRUE)
}



cgf = function(A, beta){
#    This function computes an average cumulant generating function defined in the equation (3.1) and its derivatives
#    kappa      k(beta), where beta is (k-1)-dimensional vector
#    kappad     (k-1)-dimensional vector k'(beta)
#    kappad2    (k-1) by (k-1) matrix of k''(beta)
  b = getb(A)
  k = getk(A)
  kappa = 0
  kappad = 0
  kappad2 = 0
  for (i in 1 : b) {
    perm = perm(A, i)
    expb = expb(A, i, beta)
    kappa = kappa + log(sum(expb) / prod(1 : k))
    kappad = kappad + t(perm) %*%  expb / sum(expb)
    kappad2 = kappad2 - (t(perm) %*% expb %*% t(expb) %*% perm) / sum(expb) ^ 2 +
      (t(perm) %*% diag(c(expb)) %*% perm) / sum(expb)
  }
  kappa = kappa / b
  kappad = kappad / b
  kappad2 = kappad2 / b

  list(kappa = kappa, kappad = kappad, kappad2 = kappad2)
}

getBoundaryEdge=function(A, x, tolerance = 0.0001) {
#    This function determines parameters s_1,...,s_{k-1} and l used in Theorem 3.3
 k = getk(A)
 y = x
 if (length(y) == k - 1) y[length(y) + 1] = -sum(y)
 stopifnot(length(y) == k)
 for (l in 1 : (k - 1)) {
  om = orderedMeans(A)
  sumMins = sum(om[1 : l])
  sumMaxs = sum(om[(k - l + 1) : k])
  subIndices = combinations(k, l)
  for (j in 1 : dim(subIndices)[1]) {
   s = subIndices[j, ]
   subVectorSum = sum(y[s])
   if (subVectorSum < sumMins-10^-12 || subVectorSum > sumMaxs+10^-12)
       stop("this point is not inside the closure of the domain")
   if ((subVectorSum < sumMins + tolerance) || (subVectorSum > sumMaxs - tolerance)){
        min=(subVectorSum < sumMins + tolerance)
        return (list(l = l, s = s, min=min))
    }
   }
  }
 stop("point is not near the boundaries")
}

Lambda = function(A, x, tolerance = 0.9999) {
#    Lambda functions compute the test statistic defined in the equation (3.2)
  k = getk(A)
  #x = putInsideDomain(A, x, tolerance)

  ll = function(beta) { beta %*% x - cgf(A, beta)$kappa }
  lll = function(beta) { x - cgf(A, beta)$kappad[,1] }
  llll = function(beta) { -cgf(A, beta)$kappad2 }
  lam = try(maxNR(ll, grad = lll, hess = llll, start = rep(0, k - 1), tol = 1 - tolerance), silent = TRUE)

  if ((length(lam) == 1) && (class(lam) == "try-error")) {
    return (list(checkNA =  1))
  }

  lam$checkNA = 0

  lam
}

V0h = function(A) {
#   this function computes a constant used in saddlepoint approximations
  k = getk(A)
  Sqr(cgf(A, rep(0, k - 1))$kappad2)
 }
transformIntoSphere = function(A, r, s) (r * V0h(A) %*% s)[,1]
#   this function gives point on the surface of the unit sphere corresponding to a given point r


delta = function(A, u, s, deltaIterations = 5) {
#    delta computes the equation (3.14) which used in the saddlepoint approximations given in Theorem 3.4
 k = getk(A)
 r = u
 V0h = V0h(A)
 for (i in 1 : deltaIterations) {
  tr = transformIntoSphere(A, r, s)
#   if point tr ends up on the boundary L$est becomes infinite so we choose point inside the domain close to the boundary
  while (isPointInsideDomain(A,tr)==FALSE) {
    r=.9*r
    tr = transformIntoSphere(A, r, s)
  }
  L = Lambda1(A, tr)
  r = r - ((L$max - u ^ 2 / 2) / (L$est %*% V0h %*% s))[1,1]
 }
 betah = L$est
 (det(cgf(A, betah)$kappad2) ^ (- 1 / 2) * r ^ (k - 2) * det(V0h)) / (u ^ (k - 3) * abs(t(s) %*% V0h %*% betah))
}

SpApprox = function(A, u, M, deltaIterations = 5) {
#    This function computes saddlepoint approximations to the tail probabilities of (Lambda(x) > u^2/2) given in Theorem 3.4
#    A     observed experimental values in a matrix form
#    u     parameter in the tail probability of (Lambda(x) > u^2/2)
#    M     integer number of points taken to approximate the integral over the surface of the unit sphere using MC method
#if choosing Gb = IntOnSphereMonteCarlo(delta2, k - 1, M) or M is integer number of degree parameter
#if choosing Genz method Gb = sphrul(k-1,M,delta2)

 k = getk(A)
 b = getb(A)
 delta2 = function(s) delta(A, u, s, deltaIterations)
 cb = (b ^ ((k - 1) / 2)) / (2 ^ ((k - 1) / 2 - 1) * gamma((k - 1) / 2))
# Gb = IntOnSphereMonteCarlo(delta2, k - 1, M)
 Gb = sphrul(k-1,M,delta2)$intval/( 2*pi^((k-1)/2)/(gamma((k-1)/2)))
 LugRice = 1 - pchisq(b * u ^ 2, k - 1) + cb / b * u ^ (k - 3) * exp(-b * u ^ 2 / 2) * (Gb - 1)
 ustar = u - log(Gb) / (b * u)
 NielCox = 1 - pchisq(b * ustar ^ 2, k - 1)
 list(LR=c(LugRice), BN=c(NielCox))
}



PowerForMu = function(mu, distr, Z, V) {
#    PowerForMu gives the power of the test
#    mu     k-dimensional vector parameter
#    distr  distribution function, e.g. rnorm, rexp
#    Z      integer number of random permutations of each block to obtain one p-value
#    V      integer number of p-values for each power
 pFValues = numeric(Z)
 pLRValues = numeric(Z)
 pNCValues = numeric(Z)
 for (i in 1 : Z) {
#    in this loop we calculate p-values
  A = generateRandomMatrix(distr, b, k)
  A = addMu(A, mu)
  A = A - apply(A, 1, mean)
  A = A / sqrt(mean(A ^ 2))

  Xbar = apply(A,2,mean)
  LVal = Lambdagen(A,Xbar[-k])
  u = sqrt(2*LVal)
  FVal = b*sum(Xbar ^ 2)/(mean(A^2))

  xxx = NULL; for (j in 1 : V) xxx = c(xxx, apply(t(apply(A, 1, sample, size = k)), 2, mean))
  xxx = matrix(xxx, k , V)

  FValues = b * apply(xxx ^ 2, 2, sum)

  pFValues[i] = mean(FValues > FVal)
  tt = SpApprox(A,u,M=20,deltaIterations = 5)
  pLRValues[i] = tt$LR
  pNCValues[i] = tt$NC
  }
#     calculation of power ofthe test
  PowerF = mean(pFValues < .05)
  PowerLR = mean(pLRValues < .05)
  PowerNC = mean(pNCValues < .05)

  list (PowerF = PowerF, PowerLR = PowerLR, PowerBN = PowerNC)
}

#    the following are functions which generate different distributions and can be input in #PowerForMu function
rexp2 = function(N) rexp(N) ^ 2
rgamma5 = function(N) rgamma(N, 5)
rgamma05 = function(N) rgamma(N, 0.5)

