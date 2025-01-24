#' Theme_lmda
#'
#' produces the logistic MDA theme for ggplots
#'
#' @export

theme_lmda = function(){
  theme(axis.line = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x=element_blank(), # removes dimension 1 label
        axis.title.y=element_blank(), # removes dimension 2 label
        panel.border = element_blank(),
        panel.background = element_blank())
}