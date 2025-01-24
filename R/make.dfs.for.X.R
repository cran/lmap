#' Helper function for the plot functions
#'
#' Helper function for the plot functions to add variable markers and labels 
#' to the predictor variable axes
#'
#' @param Xo Original predictor matrix 
#' @param P an integer indicating the number of predictor variables
#' @param B matrix (P x S) with weights
#' @param xnames a vector with variable names
#' @param mx averages of the original predictor matrix
#' @param sdx standard deviations of the original predictor matrix
#' @return output with information for the placement of variable markers and labels
#'
#' @export
make.dfs.for.X = function(Xo, P, B, xnames, mx, sdx){
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
    b = matrix(B[p , ], 2, 1)

    # --------------------------------------------------------------------------
    # solid line
    # --------------------------------------------------------------------------
    minx = min(Xo[, p])
    maxx = max(Xo[, p])

    # no solid lines for dichotomous predictors
    if((minx == 0) & (maxx ==1)){id.pdichotomous[p] = 1}
    else{
      m.x1 = c(minx,maxx)
      markers1 = matrix((m.x1 - mx[p])/sdx[p], 2, 1)
      markerscoord1 = outer(markers1, b) # markers1 %*% t(b %*% solve(t(b) %*% b))
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
      if(llll == 0){
        # ook een punt in de oorsprong
        MCx3[(llll + 1), 1] = NA # hoe gaat ggplot om met NA namen??
        MCx3[(llll + 1), 2] = NA #
        MCx3[(llll + 1), 3:4] = c(0,0) # referentie punt in de oorsprong
        llll = llll + 1
      }
      m.x2 = c(0, 1)
      MCx3[(llll + 1), 1] = xnames[p]
      MCx3[(llll + 1), 2] = p
      MCx3[(llll + 1), 3:4] = b
      llll = llll + 1
    } # dichotomous
    else{
      m.x2 = pretty(Xo[, p])
      m.x2 = m.x2[which(m.x2 > minx & m.x2 < maxx)]
      # to avoid clutter around origin
      m.x2 = m.x2[ which((m.x2 < (mx[p] - sdx[p])) | (m.x2 > (mx[p] + sdx[p])))]
      l.m = length(m.x2)
      markers2 = matrix((m.x2 - mx[p])/sdx[p], l.m, 1)
      markerscoord2 = outer(markers2, b) # markers2 %*% t(b %*% solve(t(b) %*% b))
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
