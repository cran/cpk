\name{dm.fn}
\alias{dm.fn}
\title{
Maintenance dose function.
}
\description{
Calculates drug's maintenance dose.
}
\usage{
dm.fn(dr, di)
}
\arguments{
  \item{dr}{
is the dose rate.
}
  \item{di}{
is the dosing interval, which can be either computed or assigned (e.g., di <- 4)
}
}

\value{
Returns the maintenance dose.
}
\references{
See \bold{cpk-package} help.
}

\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}


\examples{
dr <- 1.42; di <- 4; 
dm <- dm.fn(dr, di)
}
