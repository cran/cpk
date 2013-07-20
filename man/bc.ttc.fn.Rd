\name{bc.ttc.fn}
\alias{bc.ttc.fn}

\title{
Back-calculation function
}

\description{
Performs a back-calculation to get the target therapeutic concentration (TTC) based on the administered dose.
}

\usage{
bc.ttc.fn(dr, f, cl, wtkg = 86)
}
\arguments{
  \item{dr}{
dr is }
  \item{f}{
bioavailability 
}
  \item{cl}{
clearance rate
}
  \item{wtkg}{
weight in kg average man in US (NHANES III)
}
}

\value{
Returns the back-calculated TTC (ug/L).
}
\references{
See \bold{cpk-package} help.
}

\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}

\examples{

ttc <- 20; dr <- 2.26; cl <- 0.57; f <- 0.74; wtkg <- 86;
dr  <- dr.fn(ttc,cl,wtkg,f)

}
