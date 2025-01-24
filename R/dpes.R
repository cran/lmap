#' Dutch Parliamentary Election Study
#'
#' Data description
#'
#' @name  dpes
#' @docType data
#'
#' @usage data(dpes)
#'
#' @keywords dataset
#'
#' @format A list of 3 matrices
#' \itemize{
#'   \item X: A 275 x 5 matrix containing observed values on five predictor variables.
#'   E = Euthanasia; ID = Income differences; AS = Asylum Seekers; C = Crime; LR = Left-right scaling
#'   \item G: A 275 x 8 indicator matrix containing the responses on eight response classes
#'   PvdA, CDA, VVD, D66, GL, CU, LPF, SP.
#'   \item y: A vector of length 275 containing the responses (PvdA, CDA, VVD, D66, GL, CU, LPF, SP)
#' }
#'
#' @references Irwin, G., van Holsteyn, J., and den Ridder, J. (2003). Nationaal Kiezersonderzoek, NKO 2002 2003. DANS.
"dpes"
