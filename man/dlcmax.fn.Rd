\name{dlcmax.fn}
\alias{dlcmax.fn}
\title{ Loading dose function }

\description{Calculates a drug's loading dose based on its cmax}

\usage{ dlcmax.fn(cmax, vd) }
\arguments{
  \item{cmax}{
       Drug's maximum concentration (cmax).
}
  \item{vd}{
       Drug's apparent volume of distribution (vd).
}
}

\value{
		Returns the calculated loading dose
}

\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}

\examples{
cmax <- 18.93; vd <- 35;
dlcmax <- dlcmax.fn(cmax, vd)
}

