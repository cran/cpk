\name{css.fn}
\alias{css.fn}
\title{
\bold{css} function
}
\description{
Calculates css (steady-state) drug concentration in ng/mL, which is the same as ug/mL.
}
\usage{
css.fn(f, dpo, di, cl, ar, wtkg)
}
\arguments{
  \item{f}{
is bioavailability.
}
  \item{dpo}{
is oral dose administered in ug.
}
  \item{di}{
is the dosing interval in h.
}
  \item{cl}{
is the total clearance rate in L/kg/h.
}
  \item{ar}{
is stepwise accumulation ratio.
}
  \item{wtkg}{
is patient weight in kg.
}
}

\value{
Returns css in mg/L and ug/L
}

\references{
See \bold{cpk-package} help.
}

\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}

\examples{
f <- 0.74; dpo <- 3440; di <- 4; cl <- 0.37; ar <-1.4; wtkg <- 86;
css <-css.fn(f, dpo, di, cl, ar, wtkg)
}

