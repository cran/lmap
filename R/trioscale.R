#'  Function for TRIOSCALE
#'
#' @param y A response formula with 3 classes
#' @param X A predictor matrix
#'
#' @return This function returns an object of class trioscale
#' \item{data}{data}
#' \item{mlr}{Output object from ts.mlr}
#' \item{Q}{result from P2Q}
#' \item{X}{X matrix with coordinates}
#' \item{Xdf}{X as a data frame}
#'
#' @examples
#' \dontrun{
#' data(diabetes)
#' output = trioscale(y = diabetes$y, X = diabetes$X)
#' plot(output)
#' }
#'
#'
#' @export

trioscale = function(y, X)
  {
  # ----------------------------------------------------------------
  # ----------------------------------------------------------------
  # help functions
  # ----------------------------------------------------------------
  # ----------------------------------------------------------------
  make.dfs.for.predictors = function(mlr.output){
    
    # get some info from object
    B = mlr.output$A
    Xo = mlr.output$Xoriginal
    mx = mlr.output$mx
    sdx = mlr.output$sdx
    P = nrow(B)
    xnames = mlr.output$xnames
    
    # for solid line
    MCx1 <- data.frame(labs=character(),
                       varx = integer(),
                       dim1 = double(),
                       dim2 = double(), stringsAsFactors=FALSE)
    
    # for markers continuous variables
    MCx2 <- data.frame(labs=character(),
                       varx = integer(),
                       dim1 = double(),
                       dim2 = double(), stringsAsFactors=FALSE)
    
    # for markers dichotomous variables
    MCx3 <- data.frame(labs=character(),
                       varx = integer(),
                       dim1 = double(),
                       dim2 = double(), stringsAsFactors=FALSE)
    
    ll = 0
    lll = 0
    llll = 0
    
    id.pdichotomous = rep(0, P) # indicator for dichotomous predictors
    
    for(p in 1:P){
      markers = outer(rep(1,2), mx)
      # --------------------------------------------------------------------------
      # solid line
      # --------------------------------------------------------------------------
      minx = min(Xo[, p])
      maxx = max(Xo[, p])
      
      # no solid lines for dichotomous predictors
      if((minx == 0) & (maxx == 1)){id.pdichotomous[p] = 1}
      else{
        m.x1 = c(minx,maxx)
        markers[, p] = m.x1 # (m.x1 - mx[p])/sdx[p]
        PP = predict.mlr(mlr.output, newX = markers)
        markerscoord1 = ts.p2q2x(PP)
        MCx1[(ll + 1): (ll + 2), 1] = paste0(c("min", "max"), p)
        MCx1[(ll + 1): (ll + 2), 2] = p
        MCx1[(ll + 1): (ll + 2), 3:4] = markerscoord1
        ll = ll + 2
      }
      
      # --------------------------------------------------------------------------
      # markers
      # --------------------------------------------------------------------------
      
      if((minx == 0) & (maxx ==1)){
        # for dichotomous variables only a point indicating the 1 category
        # assumes the name of the variable indicates the 1 category
        # that is: not gender but female
        markers3 = outer(rep(1,2), mx)
        markers3[ , p] = c(0, 1)
        PPP = predict.mlr(mlr.output, newX = markers3)
        if(llll == 0){
          # ook een punt in de oorsprong
          markerscoord3 = ts.p2q2x(PPP)
          MCx3[1:2, 1] = xnames[p]
          MCx3[1:2, 2] = p
          MCx3[1:2, 3:4] = markerscoord3
          MCx3[1, 1] = NA # hoe gaat ggplot om met NA namen??
          MCx3[1, 2] = NA # 
          # llll = llll + 1
        }
        else{
          markerscoord3 = ts.p2q2x(PPP[2, , drop = FALSE])
          MCx3[(llll + 1), 1] = xnames[p]
          MCx3[(llll + 1), 2] = p
          MCx3[(llll + 1), 3:4] = markerscoord3
        }
        llll = llll + 1
      } # dichotomous
      else{
        m.x2 = pretty(Xo[, p], n = 7)
        m.x2 = m.x2[which(m.x2 > minx & m.x2 < maxx)]
        m.x2 = m.x2[ which((m.x2 < (mx[p] - sdx[p])) | (m.x2 > (mx[p] + sdx[p])))]
        l.m = length(m.x2)
        markers2 = outer(rep(1,l.m), mx)
        markers2[ , p] = m.x2 # (m.x2 - mx[p])/sdx[p]
        PPP = predict.mlr(mlr.output, newX = markers2)
        markerscoord2 = ts.p2q2x(PPP)
        MCx2[(lll + 1): (lll + l.m), 1] = paste(m.x2)
        MCx2[(lll + 1): (lll + l.m), 2] = p
        MCx2[(lll + 1): (lll + l.m), 3:4] = markerscoord2
        lll = lll + l.m
      } # ! dichotomous
    } # loop p
    
    id.pdichotomous = (id.pdichotomous == 1)
    
    output = list(MCx1 = MCx1,
                  MCx2 = MCx2,
                  MCx3 = MCx3,
                  dichotomous = id.pdichotomous)
    
    return(output)
  }
  # 
  ts.p2q2x = function(P){
    Q = cbind(log(P[, 2] / P[, 1]), log(P[, 3] / P[, 1]), log(P[, 3] / P[, 2]))
    X = cbind(Q[, 1], (2 * Q[, 2] - Q[, 1]) / sqrt(3))
    return(X)
  } 
  # ----------------------------------------------------------------
  # ----------------------------------------------------------------
  # COMPUTATIONS
  # ----------------------------------------------------------------
  # ----------------------------------------------------------------

  # fit a multinomial logistic regression
  mlr.output = mlr(y = y, X = X)
  
  # transform fitted probabilities to trioscale coordinates
  U = ts.p2q2x(mlr.output$Ghat)
  Udf = as.data.frame(U)
  colnames(Udf) = c("dim1", "dim2")
  
  # make data frames for predictor variable axes
  dfxs = make.dfs.for.predictors(mlr.output)
  
  # ----------------------------------------------------------------
  # ----------------------------------------------------------------
  # OUTPUT OBJECT
  # ----------------------------------------------------------------
  # ----------------------------------------------------------------
  output = list(
    y = y,
    X = X,
    mlr = mlr.output,
    U = U,
    Udf = Udf,
    dfxs = dfxs
  )
  # output = list(
  #   mlr = mlr.output,
  #   Q = Q,
  #   X = X,
  #   Xdf = Xdf,
  #   dens = dens, 
  #   output.lraxes = output.lraxes,
  #   output.predaxes = output.predaxes,
  #   smooth = smooth
  # )
  class(output) = "trioscale"
  return(output) 
}