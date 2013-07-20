\name{cmin.fn}
\alias{cmin.fn}
\title{
Calculates \bold{cmin}
}

\description{
Calculates cmin (trough) drug concentration in ug/L, which is the same as ng/ml.
}
\usage{
cmin.fn(cmax, ke, di)
}

\arguments{
  \item{cmax}{
is the maximum concentration of drug (peak level).
}
  \item{ke}{
Drug's total elimination rate constant.
}
  \item{di}{
is the dosing interval in h.
}
}

\value{
Returns cmin ug/L.
}

\references{
See \bold{cpk-package} help.
}

\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}

\examples{
cmax <- 25.37; ke <- 0.1333; di <- 4;
cmin <- cmin.fn(cmax, ke, di)
}
