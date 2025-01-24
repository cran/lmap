#' Function for reading the drug consumption data from the UCI repository
#'
#' @return X first coordinate matrix
#' @return Y matrix with response variables (Alcohol,Am,Amyl,Be,Caff,Ca,Choc,Co,Crack,Ex,Heroin,Ke,Le,LSD,Me,Mu,Ni,Semer,VSA)
#' @return idx indicator which response variables have a probability between 0.1 and 0.9
#'
#' @references Fehrman, E., Muhammad, A. K., Mirkes, E. M., Egan, V., & Gorban, A. N. (2017).
#' The five factor model of personality and evaluation of drug consumption risk. In Data science:
#' innovative developments in data analysis and clustering (pp. 231-242). Springer International Publishing.
#'
#' @importFrom utils read.table
#' @export

read_drugdata = function(){
  drugdat <- read.table('https://archive.ics.uci.edu/ml/machine-learning-databases/00373/drug_consumption.data', sep = ",")
  for (v in 14:32){
    drugdat[,v] = ifelse(drugdat[, v] == "CL3", 1,
                         ifelse(drugdat[, v] == "CL4", 1,
                                ifelse(drugdat[, v] == "CL5",1,
                                       ifelse(drugdat[, v] == "CL6", 1, 0))))
  }


  colnames(drugdat) = c("id", "age", "gender", "educ", "country", "ethnic",
                        "N","E","O","A","C","impulse","SS",
                        "Alcohol","Am","Amyl","Be","Caff","Ca","Choc","Co","Crack","Ex","Heroin",
                        "Ke","Le","LSD","Me","Mu","Ni","Semer","VSA")

  X = as.matrix(drugdat[,c(2,3,7:13)])
  Y = as.matrix(drugdat[,14:32])
  idx = which(colMeans(Y) > 0.1 & colMeans(Y) < 0.9)

  output = list(X = X, Y = Y, idx = idx)
  return(output)

}

# add variable names


