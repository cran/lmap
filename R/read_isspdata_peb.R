#' Function to read in the ISSP data
#' It requires the file ZA7650_v1-0-0.sav to be on your computer
#' this file can be obtained from
#' /www.gesis.org/en/issp/modules/issp-modules-by-topic/environment/2020
#' ZA7650 Data file Version 1.0.0, https://doi.org/10.4232/1.13921.
#'
#' @param path path where the ZA7650_v1-0-0.sav file is saved
#' @return X matrix containing the predictor variables
#' @return Y matrix with response variables
#'
#' @references ISSP Research Group (2022). International social survey programme:
#' Environment IV - issp 2020. GESIS, Cologne.
#'
#' @importFrom stats na.omit
#' @importFrom haven read_sav
#' @export

read_isspdata_peb = function(path){

  originalpath = getwd() # gets original working directory
  setwd(path) # set path to working directory of the data
  ZA7650data <- read_sav("ZA7650_v1-0-0.sav")
  setwd(originalpath) # set working directory back to original

  # selection of variables from data base
  mydata = ZA7650data[,c(4:68, 70:89)]
  mydata = mydata[,-c(68:81)]
  mydata = mydata[, c(50, 54, 56, 57, 1, 18, 33:39, 65, 66, 67, 70)] # with education years
  mydata = na.omit(mydata)
  colnames(mydata) = c("PEB1", "PEB2", "PEB3", "PEB4", "COUNTRY", "EC", "EE1", "EE2", "EE3", "EE4", "EE5", "EE6", "EE7", "SEX", "AGE", "EDU", "WORK")
  mydata$EE = apply(mydata[,7:13], 1, mean)

  # get response variables for pro-environmental behaviour
  # recode - more pro-environmental behavior higher score
  Y = mydata[, 1:4]
  Y = as.matrix(Y)
  Y[ , 1] = 6 - Y[, 1]
  Y[ , 2] = 8 - Y[, 2]
  Y[ , 3] = 5 - Y[, 3]
  Y[ , 4] = 5 - Y[, 4]

  # get predictor variables country, gender, education, age, education, EE and EC
  X = mydata[, -c(1:4, 7:13)]
  X = as.matrix(X)

  # create dummies for X
  Xgirl = X[,3] - 1 # sexe = female
  Xcountry = class.ind(X[, 1])[, -13] # Thailand = baseline
  colnames(Xcountry) = c("AT", "TW", "FI", "DE", "HU", "IS", "JP", "NZ", "PH", "RU", "SL", "CH")

  Xedu = scale(X[, 5], center = TRUE, scale = TRUE)
  colnames(Xedu) = "Eduyrs"

  Xwork = class.ind(X[, 6])[, -3] # never paid work = baseline
  colnames(Xwork) = c("W1", "W2")

  # standardized
  # X = as.matrix(cbind(Xcountry, Xgirl, Xedu, Xwork, scale(X[, 4]), scale(X[, 2]), scale(X[, 7])))
  # not standardized
  X = as.matrix(cbind(Xcountry, Xgirl, Xedu, Xwork,  X[, c(4, 2, 7)]))
  colnames(X) = c(colnames(Xcountry), "Female", colnames(Xedu), colnames(Xwork), "age", "EC", "EE")

  output = list(X = X, Y = Y)
  return(output)
}
