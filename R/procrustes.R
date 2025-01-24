#' Two procedures for procrustes analysis
#'
#' procrustes1 does orthogonal procrustes analysis
#' procrustes2 does similarity procrustes analysis
#' 
#' @param V1 first coordinate matrix
#' @param V2 second coordinate matrix
#' @return either a rotation matrix (procrustes1) or a list with a rotation matrix (R), a stretching factor (s) and a translation (t) (procrustes2)
#'
#' @export


procrustes1 = function(V1, V2){
  # orthogonal procrustes analysis
  # returns a rotation matrix
  pq = svd(t(V1) %*% V2)
  R = pq$v %*% t(pq$u)
  return(R)
}

procrustes2 = function(V1, V2){
  # similarity procrustes analysis
  # returns a rotation matrix (R), a stretching factor (s) and a translation (t)
  # V = s * V %*% R + outer(rep(1,nrow(V)), t)
  
  C = nrow(V1)
  J = diag(C) - 1/C
  pq = svd(t(V1) %*% J %*% V2)
  R = pq$v %*% t(pq$u)
  s = sum(diag(t(V1) %*% J %*% V2 %*% R))/sum(diag(t(V2) %*% J %*% V2))
  t = colMeans( V1 - s * V2 %*% R)
  output = list(s = s, t = t, R = R)
  return(output)
}