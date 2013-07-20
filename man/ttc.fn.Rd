\name{ttc.fn}
\alias{ttc.fn}
\title{
target therapeutic concentration function.
}
\description{
Calculates a drug's target therapeutic concentration (ttc) based on it's msc/Cmax and mec/Cmin.
}
\usage{
ttc.fn(msc, mec)
}
\arguments{
  \item{msc}{
maximum safe concentration or cmax
}
  \item{mec}{
minimum effective concentration of cmin
}
}

\value{
Returns target therapeutic concentration (ttc).
}
\references{
See \bold{cpk-package} help.
}

\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}
\examples{
msc = 50; mec = 20;
ttc <- ttc.fn(msc, mec)
}
