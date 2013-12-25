\name{cpk-package}
\alias{cpk-package}
\alias{cpk}
\docType{package}
\title{
Clinical Pharmacokinetics.
}

\description{
This package was written to teach concepts and techniques to clinicians and pharmacists for individualizing drug therapy. While the package may be used by clinicians for actual dosage regimen analysis, design, and modification, caution should be exercised when applying the techniques described in \R package cpk.
}
\details{
\tabular{ll}{
Package: \tab cpk\cr
Type: \tab Package\cr
Version: \tab 1.3\cr
Date: \tab 2013-07-25\cr
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

wtkg = 181;       # kg

# Drug Disposition Parameters
thalf = 4;        # h
ke    = 0.3180;   # h^-1
vd    = 4.5;      # L/kg
cl    = 1.43;     # L/h/kg
f     = 1.00;     # dpo range: 15-64% (avg 38%), IV f = 1.

#### Set TTC
ttc = 25;

#### Calculate dose rate (mg/h)
dr  <- dr.fn (ttc, cl, wtkg, f)

#### Set dosing interval based on thalf
di = 4;   # h

#### Calculate dose (if f = 1, dose IV)
dpo <- dpo.fn (dr, di)
dpo/1000 # convert to mg

#### Calculate loading dose based on cmax
cmax <- 18.93; vd <- 35;
dlcmax <- dlcmax.fn(cmax, vd)

#### Calculate loading dose based on ar
dm <- 276; ar <- 2.4;
dlar <- dlar.fn(dm, ar)

#------------------------------------
# ANALYSIS
#------------------------------------
bc.ttc <- bc.ttc.fn (dr, f, cl, wtkg)

#### predict avg steady-state concentration
# accumulation ratio
ar  <- ar.fn (ke, di)

css <- css.fn (f, dpo, di, cl, ar, wtkg)
css # mg/L

# predict fluctuation about steady-state (ug/L)
cmax <- cmax.fn (f, dpo, vd, ar, wtkg)
cmin <- cmin.fn (cmax, ke, di)

ct <- ct.fn(cmax, ke, time=0)
ct <- ct.fn(cmax, ke, time=4)

# infusion rate
css <- 14.43; cl <- 3.2; 
R0  <- R0.fn(css, cl)

}
