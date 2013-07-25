\name{di.fn}
\alias{di.fn}
\title{
\bold{di} function
}
\description{
Calculates a drug's dosing interval.
}
\usage{
di.fn(msc, mec, ke)
}
\arguments{
  \item{msc}{
Drug's maximum safe concentration or Cmax (peak) concentration.
}
  \item{mec}{
Drug's minimum effective concentration or Cmin (trough) concentration.
}
  \item{ke}{
Drug's total elimination rate constant.
}
}

\value{
Returns the dosing interval (h).
}
\references{
See \bold{cpk-package} help.
}

\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}

\examples{

msc <- 50; mec <- 20; ke <- 0.2078;
di  <- di.fn(msc, mec, ke)

}
