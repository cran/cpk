\name{cpk-package}
\alias{cpk-package}
\alias{cpk}
\docType{package}
\title{
Clinical Pharmacokinetics.
}

\description{
This package was written to teach concepts and techniques for individualizing drug therapy to clinicians and pharmacists. While the package may be used by clinicians for actual dosage regimen analysis, design, and modification, caution should be exercised when applying the techniques described in \R package cpk.
}
\details{
\tabular{ll}{
Package: \tab cpk\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2013-07-19\cr
License: \tab GPL-2\cr
}
}
\author{Oscar A. Linares MD and David T. Daly JD/MBA.\cr
Maintainer: Oscar A. Linares MD <OALinaresMD@gmail.com>
}

\references{
Holford, N and Sheiner, LB. Understanding the dose-effect relationship: clinical applications of pharmacokinetic-pharmacodynamic models. Clin Pharmacokinet 1981;6:429-453. 

Greenblatt D, Koch-Weser J. Clinical pharmacokinetics (First of two parts). New Engl J Med 1975;293:702-705.

Greenblatt D, Koch-Weser J. Clinical pharmacokinetics (Second of two parts). New Engl J Med 1975;293:964-970

Perlin E, Taylor BK, Peck CC. Clinical pharmacokinetics: a simplified approach, part 1. J Natl Med Assoc 1985;77:475.

Perlin E, Taylor R, Peck CC. Clinical pharmacokinetics: a simplified approach, part 2. J Natl Med Assoc 1986;78:835.

Shafer SL, Varvel JR. Pharmacokinetics, pharmacodynamics, and rational opioid selection. Anesthesiology 1991;74:53-63.

Mehvar R. Pharmacokinetic-based design and modification of dosage regimens. Am J Pharmaceut Educ 1998;62:189-195.

Brocks DR, Mehvar R. Rate and extent of drug accumulation after multiple dosing revisited. Clin Pharmacokinet 2010;49:421-438.

Jannetto PJ, Bratanow NC. Utilization of pharmacogenomics and therapeutic drug monitoring for opioid pain management. Pharmacogenomics 2009;10(7):1157-67.

Jannetto PJ, Bratanow NC. Pain management in the 21st century: Utilization of pharmacogenomics and therapeutic drug monitoring. Expet Opin Drug Metabol Toxicol 2012;7(6):745-52.

Linares OA. Pharmacokinetics-based opioid prescribing: Case report. The Pain Practitioner 2010;20:19-21.

Linares O, Linares A. Computational opioid prescribing: A novel application of clinical pharmacokinetics. J Pain Palliat Care Pharmacother 2011;25:125-135.

Linares OA, Daly D, Stefanovski D, Boston RC. A new model for using quantitative urine testing as a diagnostic tool for oxycodone treatment and compliance. J Pain Palliat Care Pharmacother 2013;27:(in press).
}

\keyword{ cpk package }

\examples{
ke <- 0.2350; di <- 4;
ar <- ar.fn(ke, di)

ttc <- 20; dr <- 2.26; cl <- 0.57; f <- 0.74; wtkg <- 86;
dr  <- dr.fn(ttc,cl,wtkg,f)

f = 0.74; dpo <- 3440; vd = 2.8; ar = 2.4;
cmax <- cmax.fn(f, dpo, vd, ar)

f <- 0.74; dpo <- 3440; di <- 4; cl <- 0.37;
css <-css.fn(f, dpo, di, cl)

cmax <- 25.37; ke <- 0.1333; di <- 4;
cmin <- cmin.fn(cmax, ke, di)

msc <- 50; mec <- 20; ke <- 0.2078;
di  <- di.fn(msc, mec, ke)

dr <- 1.42; di <- 4; 
dm <- dm.fn(dr, di)

dr <- 0.86; di <- 4;
dpo <- dpo.fn(dr, di)

ttc <- 32.74; cl <- 0.59; wtkg <- 86; f <- 0.74; 
dr <- dr.fn(ttc, cl, wtkg, f)

msc = 50; mec = 20;
ttc <- ttc.fn(msc, mec)

}
