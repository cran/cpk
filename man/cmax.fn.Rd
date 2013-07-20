\name{cmax.fn}
\alias{cmax.fn}
\title{
Calculates \bold{cmax}
}
\description{
Calculates cmax (peak) drug concentration in ug/L, which is the same as ng/ml.
}
\usage{
cmax.fn(f, dpo, vd, ar, wtkg = 86)
}
\arguments{
  \item{f}{
is bioavailability.
}
  \item{dpo}{
is oral dose administered in ug.
}
  \item{vd}{
is apparent volume of distribution in L.
}
  \item{ar}{
is accumulation ratio.
}
  \item{wtkg}{
is patient weight in kg.
}
}

\value{
Returns cmax ug/L.
}

\references{
See \bold{cpk-package} help.
}

\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}

\examples{
f = 0.74; dpo <- 3440; vd = 2.8; ar = 2.4;
cmax <- cmax.fn(f, dpo, vd, ar)
}

