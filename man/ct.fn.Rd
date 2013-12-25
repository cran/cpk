\name{ct.fn}
\alias{ct.fn}
\title{
\bold{ct} function
}
\description{
Predicts drug concentration at specified time (default time=0 h) [ng/mL].
}
\usage{
ct.fn(cmax, ke, time=0)
}
\arguments{
  \item{cmax}{
is the maximum concentration (cmax).
}
  \item{ke}{
is the first-order elimination rate constant (ke).
}
  \item{time}{
is the time on the clearance curve for which a drug concentration is desired.
}
}

\value{
Returns the concentration at specified time.
}
\references{
See \bold{cpk-package} help.
}

\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}


\examples{
cmax <- 15; ke <- 0.1333; time <- 4; 
ct <- ct.fn(cmax, ke, time=0)
}
