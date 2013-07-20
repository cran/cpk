\name{ar.fn}
\alias{ar.fn}
\title{ Accumulation ratio function }

\description{Calculates drug's accumulation ratio}

\usage{ ar.fn(ke, di) }
\arguments{
  \item{ke}{
       Drug's total elimination rate constant ke (h^-1).
}
  \item{di}{
       Dosing interval (h).
}
}

\value{
		Returns the calculated accumulation ratio
}

\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}

\examples{
ke <- 0.2350; di <- 4;
ar <- ar.fn(ke, di)
}

