#' Liver
#'
#' Data description
#'
#' @name  liver
#' @docType data
#'
#' @usage data(liver)
#'
#' @keywords dataset
#'
#' @format A list of 3 matrices
#' \itemize{
#'   \item X: A 218 x 3 matrix containing observed values on three liver function test (predictor variables).
#'   ASpartate aminotransferase (AS), ALanine aminotransferase (AL), and Glutamate Dehydrogenase (GD)
#'   \item G: A 218 x 4 indicator matrix containing the responses on
#'   Acute Viral Hepatitis (AVH), Persistent Chronic Hepatitis (PCH), Aggressive Chronic Hepatitis (ACH), Post-Necrotic Cirrhosis (PNC).
#'   \item y: A vector of length 218 containing the responses (AVH, PCH, )
#' }
#'
#' @references Plomteux, G. (1980). Multivariate analysis of an enzymic profile for the differential diagnosis of viral hepatitis.
#' Clinical Chemistry, 26(13), 1897–1899.
"liver"
