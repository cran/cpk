# Oscar A. Linares MD and David T. Daly JD/MBA
# Plymouth Pharmacokinetic Modeling Study Group
#
# cpk Version 1.2
# 7/20/2013

#####################################
# R FUNCTIONS
#####################################

# target therapeutic concentration
ttc.fn <- function(msc, mec) {
  numerator   <- msc - mec
  denominator <- log(msc/mec)
  result      <- numerator/denominator
  print(paste("The value of ttc (ug/L) is", round(result,2), sep=" "))
  return(round(result,0))
}

# dose rate
dr.fn <- function(ttc, cl, wtkg, f) {
  
  numerator   <- ttc * cl * wtkg
  denominator <- f
  result      <- numerator/denominator * 0.001
  print(paste("The value of dr (mg/h) is", round(result,2), sep=" "))
  return(round(result,2))
  
}

# dosing interval
di.fn <- function (msc, mec, ke) {
  numerator   <- log(msc/mec)
  denominator <- ke
  result      <- numerator/denominator
  print(paste("The value of di (h) is", round(result,1), sep=" "))
  return(round(result,1))
  
}

# maintenance dose
dm.fn <- function (dr, di) {
  result <- dr * di
  print(paste("The value of dm (mg every di h) is", round(result,2), sep=" "))
  return(round(result,2))
}

# back calculation TTC
bc.ttc.fn <- function (dr, f, cl, wtkg = 86) {
  numerator   <- dr * 1000 * f
  denominator <- cl * wtkg
  result      <- numerator/denominator
  print(paste("The value of bc.ttc (ug/L) is", round(result,1), sep=" "))
  return(round(result,1))
}

# accumulation rate
ar.fn <- function(ke, di) {
  numer  <- 1
  denom  <- 1 - exp(-ke * di)
  result <- numer/denom
  print(paste("The value of ar is", round(result,1), sep=" "))
  return(round(result,1))
}

# Oral dose
dpo.fn <- function(dr, di) {
  result      <- dr * di * 1000
  print(paste("The value of dpo (ug) is", round(result,2), sep=" "))
  return(round(result,2))
}

# Cmax
cmax.fn <- function(f, dpo, vd, ar, wtkg = 86) {
  numer <- f * dpo
  denom <- vd * wtkg
  result      <- (numer/denom) * ar
  print(paste("The value of cmax (ug/L) is", round(result,2), sep=" "))
  return(round(result,2))
}

# Cmin
cmin.fn <- function(cmax, ke, di) {
  result      <- cmax * exp(-ke * di) 
  print(paste("The value of cmin (ug/L) is", round(result,2), sep=" "))
  return(round(result,2))
}

# Css
css.fn <- function(f, dpo, di, cl, wtkg = 86) {
  numer  <- f * (dpo/di)
  denom  <- cl * wtkg
  result <- (numer/denom)  
  print(paste("The value of css (ug/L) is", round(result,2), sep=" "))
  return(round(result,2))
}
####################### END ####################### 
