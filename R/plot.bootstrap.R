#' Plot an object obtained using one of the bootstrap functions
#'
#' @param x an object of type bootstrap
#' @param level level of confidence regios
#' @param type choose between "Bag" (default) or "norm"
#' @param \dots additional arguments to be passed.
#' @return Plot of the results obtained from bootstrap
#'
#'
#' @import ggplot2
#' @import ggpubr
#' @importFrom grid polygonGrob gpar grobName
#' @importFrom stats setNames
#' @importFrom grDevices chull
#'
#' @export
plot.bootstrap = function(x, level = 0.95, type = "Bag", ...){

  StatBag <- ggproto("Statbag", Stat,
                     compute_group = function(data, scales, prop = 0.5) {

                       #################################
                       #################################
                       # originally from aplpack package, plotting functions removed
                       plothulls_ <- function(x, y, fraction, n.hull = 1,
                                              col.hull, lty.hull, lwd.hull, density=0, ...){
                         # function for data peeling:
                         # x,y : data
                         # fraction.in.inner.hull : max percentage of points within the hull to be drawn
                         # n.hull : number of hulls to be plotted (if there is no fractiion argument)
                         # col.hull, lty.hull, lwd.hull : style of hull line
                         # plotting bits have been removed, BM 160321
                         # pw 130524
                         if(ncol(x) == 2){ y <- x[,2]; x <- x[,1] }
                         n <- length(x)
                         if(!missing(fraction)) { # find special hull
                           n.hull <- 1
                           if(missing(col.hull)) col.hull <- 1
                           if(missing(lty.hull)) lty.hull <- 1
                           if(missing(lwd.hull)) lwd.hull <- 1
                           x.old <- x; y.old <- y
                           idx <- chull(x,y); x.hull <- x[idx]; y.hull <- y[idx]
                           for( i in 1:(length(x)/3)){
                             x <- x[-idx]; y <- y[-idx]
                             if( (length(x)/n) < fraction ){
                               return(cbind(x.hull,y.hull))
                             }
                             idx <- chull(x,y); x.hull <- x[idx]; y.hull <- y[idx];
                           }
                         }
                         if(missing(col.hull)) col.hull <- 1:n.hull
                         if(length(col.hull)) col.hull <- rep(col.hull,n.hull)
                         if(missing(lty.hull)) lty.hull <- 1:n.hull
                         if(length(lty.hull)) lty.hull <- rep(lty.hull,n.hull)
                         if(missing(lwd.hull)) lwd.hull <- 1
                         if(length(lwd.hull)) lwd.hull <- rep(lwd.hull,n.hull)
                         result <- NULL
                         for( i in 1:n.hull){
                           idx <- chull(x,y); x.hull <- x[idx]; y.hull <- y[idx]
                           result <- c(result, list( cbind(x.hull,y.hull) ))
                           x <- x[-idx]; y <- y[-idx]
                           if(0 == length(x)) return(result)
                         }
                         result
                       } # end of definition of plothulls
                       #################################


                       # prepare data to go into function below
                       the_matrix <- matrix(data = c(data$x, data$y), ncol = 2)

                       # get data out of function as df with names
                       setNames(data.frame(plothulls_(the_matrix, fraction = prop)), nm = c("x", "y"))
                       # how can we get the hull and loop vertices passed on also?
                     },
                     required_aes = c("x", "y")
  )

  #' @param prop Proportion of all the points to be included in the bag (default is 0.5)
  stat_bag <- function(mapping = NULL, data = NULL, geom = "polygon",
                       position = "identity", na.rm = FALSE, show.legend = NA,
                       inherit.aes = TRUE, prop = 0.5, alpha = 0.3, ...) {
    layer(
      stat = StatBag, data = data, mapping = mapping, geom = geom,
      position = position, show.legend = show.legend, inherit.aes = inherit.aes,
      params = list(na.rm = na.rm, prop = prop, alpha = alpha, ...)
    )
  }


  geom_bag <- function(mapping = NULL, data = NULL,
                       stat = "identity", position = "identity",
                       prop = 0.5,
                       alpha = 0.3,
                       ...,
                       na.rm = FALSE,
                       show.legend = NA,
                       inherit.aes = TRUE) {
    layer(
      data = data,
      mapping = mapping,
      stat = StatBag,
      geom = GeomBag,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        na.rm = na.rm,
        alpha = alpha,
        prop = prop,
        ...
      )
    )
  }

  #' @rdname ggplot2-ggproto
  #' @format NULL
  #' @usage NULL
  #' @export
  GeomBag <- ggproto("GeomBag", Geom,
                     draw_group = function(data, panel_scales, coord) {
                       n <- nrow(data)
                       if (n == 1) return(zeroGrob())

                       munched <- coord_munch(coord, data, panel_scales)
                       # Sort by group to make sure that colors, fill, etc. come in same order
                       munched <- munched[order(munched$group), ]

                       # For gpar(), there is one entry per polygon (not one entry per point).
                       # We'll pull the first value from each group, and assume all these values
                       # are the same within each group.
                       first_idx <- !duplicated(munched$group)
                       first_rows <- munched[first_idx, ]

                       ggname("geom_bag",
                                        polygonGrob(munched$x, munched$y, default.units = "native",
                                                           id = munched$group,
                                                           gp = gpar(
                                                             col = first_rows$colour,
                                                             fill = alpha(first_rows$fill, first_rows$alpha),
                                                             lwd = first_rows$size * .pt,
                                                             lty = first_rows$linetype
                                                           )
                                        )
                       )


                     },
                     default_aes = aes(colour = "NA", fill = "grey20", size = 0.5, linetype = 1,
                                       alpha = NA, prop = 0.5),
                     handle_na = function(data, params) {
                       data
                     },
                     required_aes = c("x", "y"),
                     draw_key = draw_key_polygon
  )

  object <- x

  B = as.data.frame(object$obj$B)
  colnames(B) = paste0("dim", 1:ncol(B))
  rownames(B) = object$obj$xnames
  P = nrow(B)
  B$Predictor = c(1:P)
  B$Predictor = factor(B$Predictor, levels = 1:P, labels = rownames(B))

  V = as.data.frame(object$obj$V)
  colnames(V) = paste0("dim", 1:ncol(V))
  rownames(V) = object$obj$ynames
  R = nrow(V)
  V$Response = c(1:R)
  V$Response = factor(V$Response, levels = 1:R, labels = rownames(V))

  Bplot = ggplot(object$BBdf, aes(.data$dim1, .data$dim2, color = .data$Predictor), show.legend = F) +
    geom_point(shape = ".", show.legend = F) +
    stat_ellipse(type = "norm", show.legend = F) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    coord_fixed() + theme_bw() +
    labs(title = "Weights plot")# , x = "Dimension 1", y = "Dimension 2")
  
  if(type == "Bag"){
    Bplot = Bplot + geom_bag(prop = level, aes(fill = .data$Predictor, colour = .data$Predictor), alpha = 0.1, show.legend = FALSE) 
  }
  else if(type == "norm"){
    Bplot = Bplot + stat_ellipse(type = "norm", level = level, show.legend = F)
  }
    geom_text(data = B, aes(x = .data$dim1, y = .data$dim2, label = .data$Predictor), 
							color = "black", family = 'mono', fontface = 'bold', show.legend = F)

  Vplot = ggplot(object$BVdf, aes(.data$dim1, .data$dim2, color = .data$Response), show.legend = F) +
    geom_point(shape = ".", show.legend = F) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    coord_fixed() + theme_bw() +
    labs(title = "Loadings plot") # , x = "Dimension 1", y = "Dimension 2")

  Vplot = Vplot +
    geom_text(data = V, aes(x = .data$dim1, y = .data$dim2, label = .data$Response), 
              color = "black", family = 'mono', fontface = 'bold', show.legend = F)
  
  if(type == "Bag"){
    Vplot = Vplot + geom_bag(prop = level, aes(fill = .data$Response, colour = .data$Response), alpha = 0.1, show.legend = FALSE)
  }
  else if(type == "norm"){
    Vplot = Vplot + stat_ellipse(type = "norm", level = level, show.legend = F)
  }
  
  Bplot = Bplot + theme_lmda()
  Vplot = Vplot + theme_lmda()
  
  # combine the two plots
  BVplot = ggarrange(Bplot, Vplot, ncol = 2, align = "h")
  print(BVplot)

  output = list(
    BVplot = BVplot,
    Bplt = Bplot,
    Vplt = Vplot
  )
  return(output)
}



ggname<-function(prefix, grob){
  grob$name <- grobName(grob, prefix)
  grob
}
