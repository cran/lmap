% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kieskompas.R
\docType{data}
\name{kieskompas}
\alias{kieskompas}
\title{Kieskompas data}
\format{
A list of 7 matrices
\itemize{
  \item G: an indicator matrix of dimension 25001 by 21 indicating the vote intention of the participants
  \item Xs: responses to the 30 propositions of the participants
  \item Xs2: responses to the 30 propositions of the participants, with NA's coded as neutral (3)
  \item Xb: background variables Age, Gender (male = 0, female = 1), and Education (theoretical = 0, practical = 0)
  \item Z: position of the 18 political parties on the 30 propositions
  \item stellingen: the exact Dutch wording of the propositions
  \item grouping: the 30 propositions can be grouped in themes (mk, bb, ef, sc, ai, bo, et, oz).
  mk: Environment & Climate; bb: Foreign Policy; ef: Economy & Finance; sc: Social Affairs & Culture
  ai: Asylum & Immigration; cc: Construction & Environment; et: Ethics; oz: Education & Healthcare
}
}
\usage{
kieskompas
}
\description{
Data description
}
\references{
van Lindert, J., Meijer, S., Etienne, T., van der Steen, S., Kutiyski, Y., Moreda Laguna, O., Brousianou, A.,
 Blanken, W., & Krouwel, A. (2023). Het Kieskompas voor de Nederlandse Tweede Kamerverkiezingen van 2023 [dataset].
 Kieskompas, Amsterdam, the Netherlands.
}
\keyword{dataset}
