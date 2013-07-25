\name{dlar.fn}
\alias{dlar.fn}
\title{ Loading dose function }

\description{Calculates a drug's loading dose based on its accumulation ratio}

\usage{ dlar.fn(dm, ar) }
\arguments{
  \item{dm}{
       Drug's maintenance dose.
}
  \item{ar}{
       Drug's accumulation ratio.
}
}

\value{
		Returns the calculated loading dose
}

\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}

\examples{
dm <- 276; ar <- 2.4;
dlar <- dlar.fn(dm, ar)
}

