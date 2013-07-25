\name{R0.fn}
\alias{R0.fn}
\title{
\bold{R0} function
}
\description{
Predicts drug infusion rate [ug/h].
}
\usage{
R0.fn(css, cl)
}
\arguments{
  \item{css}{
is the steady-state concentration.
}
  \item{cl}{
is the clearance rate.
}
}

\value{
Returns the infusion rate.
}
\references{
See \bold{cpk-package} help.
}

\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}


\examples{
css <- 14.43; cl <- 3.2; 
R0  <- R0.fn(css, cl)
}
