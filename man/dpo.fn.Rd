\name{dpo.fn}
\alias{dpo.fn}
\title{
Oral dose of drug.
}
\description{
Calculates the oral dose of drug administered.
}
\usage{
dpo.fn(dr, di)
}
\arguments{
  \item{dr}{
is the dose rate in mg/h.
}
  \item{di}{
is the dosing interval in h.
}
}

\value{
Returns oral dose administered in ug.
}

\references{
See \bold{cpk-package} help.
}

\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}

\examples{
dr <- 0.86; di <- 4;
dpo <- dpo.fn(dr, di)
}
