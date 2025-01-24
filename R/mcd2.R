#' Multinomial Canonical Decomposition Model for a multinomial outcome
#'
#' The function mcd2 fits the multinomial canonical decomposition model to a multinomial outcome
#' i.e. a double constrained reduced rank multinomial logistic model
#'
#' @param G An N times C class indicator matrix
#' @param X An N by P matrix with predictor variables
#' @param Z design matrix for response
#' @param S Positive number indicating the dimensionality of teh solution
#' @param trace whether progress information should be printed on the screen
#' @param maxiter maximum number of iterations
#' @param dcrit convergence criterion
#'
#' @return This function returns an object of the class \code{mcd} with components:
#' \item{call}{function call}
#' \item{Xoriginal}{Matrix X from input}
#' \item{G}{Class indicator matrix G}
#' \item{X}{Scaled X matrix}
#' \item{mx}{Mean values of X}
#' \item{sdx}{Standard deviations of X}
#' \item{pnames}{Variable names of profiles}
#' \item{xnames}{Variable names of predictors}
#' \item{znames}{Variable names of responses}
#' \item{Z}{Design matrix Z}
#' \item{m}{main effects}
#' \item{Bx}{regression weights for X}
#' \item{Bz}{regression weights for Z}
#' \item{A}{regression weights (Bx Bz')}
#' \item{U}{matrix with coordinates for row-objects}
#' \item{V}{matrix with coordinates for column-objects}
#' \item{Ghat}{Estimated values of G}
#' \item{deviance}{value of the deviance at convergence}
#' \item{df}{number of paramters}
#' \item{AIC}{Akaike's informatoin criterion}
#' \item{iter}{number of main iterations from the MM algorithm}
#' \item{svd}{Singular value decomposition in last iteration}
#'
#' @examples
#' \dontrun{
#' data(dataExample_lpca)
#' Y = as.matrix(dataExample_lpca[ , 1:5])
#' X = as.matrix(dataExample_lpca[ , 9:13])
#' #unsupervised
#' output = mcd1(X, Y, S = 2, ord.z = 2)
#' }
#'
#' @importFrom stats plogis model.matrix formula
#' @importFrom nnet class.ind
#'
#'
#' @export
#' 
#' 
mcd2 = function(X, G, Z, S = 2, trace = TRUE, maxiter = 65536, dcrit = 1e-6){
  # multinomial canonical decomposition for multinomial outcome
  # i.e. double constrained reduced rank multinomial logistic model
  # preparation
  cal = match.call()
  
  n = nrow(X)
  ones.n = matrix(1,n,1)
  P = ncol(X)
  C = ncol(G)
  Q = ncol(Z)
  # G = Y # class indicator matrix
  
  # out = make.profile(Y = Y, ord = ord.z)
  # if(is.null(Z)){ Z = out$Z }
  # G = out$G
  # A = out$A
  # Q = ncol(Z)  
  
  # if(is.null(W)){
  #   if(ord.m ==1){
  #     W = model.matrix(~., data = A)[, -1]
  #     qrz = qr(model.matrix(~., data = A)[, -1])
  #   }
  #   else{
  #     W = model.matrix(formula(paste("~.^", ord.m, sep = "")), data = A)[, -1]
  #     qrz = qr(model.matrix(formula(paste("~.^", ord.m, sep = "")), data = A)[ , -1])
  #   }
  # }
  # else{
  #   qrz = qr(W)
  # }
  # iWWW = solve( t(W) %*% W + 1e-6 * diag(ncol(W)) ) %*% t(W)
  # bm = solve.qr(qrz, m)
  # m = qr.fitted(qrz, m)
  
  # scaling of predictor matrices
  Xoriginal = X
  outx = procx(X)
  X = outx$X
  mx = outx$mx
  sdx = outx$sdx
  
  if(P == 1){
    iRx = 1/sqrt(t(X) %*% X)
  }
  else{
    eig.out = eigen(t(X) %*% X)
    iRx = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
  }
  
  if(ncol(Z) == 1){
    iRz = 1/sqrt(t(Z) %*% Z)
  }
  else{
    eig.out = eigen(t(Z) %*% Z)
    iRz = eig.out$vectors %*% diag(1/sqrt(eig.out$values)) %*% t(eig.out$vectors)
  }
  
  # initialization
  m = colMeans(G)
  udv = svd(1/sqrt(n) * iRx %*% t(X) %*% G %*% Z %*% iRz)
  Bx = iRx %*% matrix(udv$u[, 1:S], P, S) %*% diag(udv$d[1:S], nrow = S, ncol = S)
  Bz = iRz %*% matrix(udv$v[, 1:S], Q, S) 
  theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
  Ghat = exp(theta) / rowSums(exp(theta))
  dev.old = -2 * sum(log(Ghat[G == 1]))
  
  # iteration
  iter = 0; dif = 1
  while(dif > dcrit){
    iter = iter + 1
    # update m
    ZZ = theta + 4 * (G - Ghat)
    m = colMeans(ZZ - X %*% Bx %*% t(Z %*% Bz))

    # update Bx (U) and Bz (V)
    theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
    Ghat = exp(theta) / rowSums(exp(theta))
    ZZ = theta + 4 * (G - Ghat) - ones.n %*% t(m)
    udv = svd(iRx %*% t(X) %*% ZZ %*% Z %*% iRz)
    Bx = iRx %*% matrix(udv$u[, 1:S], P, S) %*% diag(udv$d[1:S], nrow = S, ncol = S) 
    Bz = iRz %*% matrix(udv$v[, 1:S], Q, S) 
    
    # deviance
    theta = ones.n %*% t(m) + X %*% Bx %*% t(Z %*% Bz)
    Ghat = exp(theta) / rowSums(exp(theta))
    dev.new = -2 * sum(log(Ghat[G == 1]))
    
    # convergence
    dif = 2 * (dev.old - dev.new)/ ((dev.old + dev.new))
    if ( trace ) cat(iter, dev.old, dev.new, dif, "\n")
    if ( dif < dcrit ) break
    if ( iter == maxiter ) warning("Maximum number of iterations reached - not converged (yet)")
    dev.old = dev.new
  } # end iteration
  
  BB = Bx %*% t(Bz)
  
  # name giving
  rownames(Bz) = colnames(Z)
  rownames(Bx) = colnames(X)
  rownames(BB) = colnames(X)
  colnames(BB) = colnames(Z)
  # rownames(bm) = colnames(W)
  
  npar = C + (P + Q - S) * S
  # create output object
  results = list(
    call = cal,
    Xoriginal = Xoriginal,
    X = X,
    mx = mx,
    sdx = sdx,
    G = G,
    xnames = colnames(X),
    znames = colnames(Z),
    Z = Z,
    m = m,
    Bx = Bx,
    Bz = Bz,
    A = BB,
    U = X %*% Bx,
    V = Z %*% Bz,
    Ghat = Ghat,
    deviance = dev.new,
    df = npar,
    AIC = dev.new + 2 * npar,
    iter = iter,
    svd = udv,
    type = "mcd2"
  )
  class(results) <- "mcd"
  return(results)
}
