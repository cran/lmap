#' This function compares the predictive performance of several models fitted on the same data
#'
#' The number of bootstraps should be the same for each model
#' Ideally, the seed used in bootstrapping should also be the same
#'
#' @param objectlist An list with output objects from the bootstrap functions in lmap
#' @param xlabel A character object, specifying the label on the horizontal axis. Default is "Model"
#' @return plot A boxplot with prediction errors for each model
#' @return pe A data frame with average prediction error for each bootstrap
#' @return fit A matrix with prediction error statistics for each model
#'
#' @examples
#' \dontrun{
#' data(dataExample_mru)
#' y = as.matrix(dataExample_mru[ , 1])
#' X = as.matrix(dataExample_mru[ , 2:6])
#' output2 = mrrr(y = y, X = X, S = 2)
#' b2 = bootstrap.mrrr(output2)
#' output3 = mrrr(y = y, X = X, S = 3)
#' b3 = bootstrap.mrrr(output3)
#' myobjects = list(b2, b3)
#' comparison = oos.comparison(myobjects)
#' comparison$plot
#' comparison$fit
#' }
#' @import ggplot2
#' @importFrom grDevices rgb
#' @importFrom stats median sd
#'
#' @export

oos.comparison = function(objectlist, xlabel = "Model"){

  M = length(objectlist)
  Bsamples = length(objectlist[[1]]$sdev.oos)

  df = data.frame(numeric(),numeric())
  colnames(df) = c("model", "oos")


  for(m in 1:M){
    df[((m-1)*Bsamples + 1) : (m * Bsamples) , 1] = m
    df[((m-1)*Bsamples + 1) : (m * Bsamples) , 2] = objectlist[[m]]$sdev.oos
  }

  # make boxplot
  plt = ggplot(df, aes(x = factor(.data$model), y = .data$oos, fill=factor(.data$model))) +
    geom_boxplot(aes(group = factor(.data$model))) +
    geom_jitter(width = 0.05, height = 0, colour = rgb(0,0,0,.3)) +
    ylab("Prediction error") + xlab(xlabel) +
    theme_bw() +
    theme(legend.position = "none")

  # make table
  fit = matrix(NA, M, 6)
  colnames(fit) = c("model", "mean", "median", "min", "max", "sd")
  fit[ , 1] = 1:M
  fit[ , 2] = sapply(objectlist, function(l){mean(l$sdev.oos)})
  fit[ , 3] = sapply(objectlist, function(l){median(l$sdev.oos)})
  fit[ , 4] = sapply(objectlist, function(l){min(l$sdev.oos)})
  fit[ , 5] = sapply(objectlist, function(l){max(l$sdev.oos)})
  fit[ , 6] = sapply(objectlist, function(l){sd(l$sdev.oos)})

  return = list(plot = plt, pe = df, fit = fit)
}
