\name{dr.fn}
\alias{dr.fn}
\title{
Dose rate function
}
\description{
Calculates dose rate.
}
\usage{
dr.fn(ttc, cl, wtkg, f)
}
\arguments{
  \item{ttc}{
target therapeutic concentration.
}
  \item{cl}{
clearance rate.
}
  \item{wtkg}{
patient's weight in kilograms.
}
  \item{f}{
bioavailability.
}
}

\value{
Returns the dose rate.
}
\references{
See \bold{cpk-package} help.
}
\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}

\examples{

ttc <- 32.74; cl <- 0.59; wtkg <- 86; f <- 0.74; 
dr <- dr.fn(ttc, cl, wtkg, f)

}
