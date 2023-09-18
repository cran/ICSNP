\name{LASERI}
\Rdversion{1.1}
\alias{LASERI}
\docType{data}
\title{
Cardiovascular Responses to Head-up Tilt
}
\description{
This data set contains the cardiovascular responses to a passive head-up tilt for 223 subjects.
}
\usage{data(LASERI)}
\format{
  A data frame with 223 observations on the following 32 variables.
  \describe{
    \item{\code{Sex}}{a factor with levels \code{Female} and \code{Male}.}
    \item{\code{Age}}{Age in years.}
    \item{\code{Height}}{Height in cm.}
    \item{\code{Weight}}{Weight in kg.}
    \item{\code{Waist}}{Waist circumference in cm.}
    \item{\code{Hip}}{Hip circumference in cm.}
    \item{\code{BMI}}{Body mass index.}
    \item{\code{WHR}}{Waist hip ratio.}
    \item{\code{HRT1}}{Average heart rate in the tenth minute of rest.}
    \item{\code{HRT2}}{Average heart rate in the second minute during the tilt.}
    \item{\code{HRT3}}{Average heart rate in the fifth minute during the tilt.}
    \item{\code{HRT4}}{Average heart rate in the fifth minute after the tilt.}
    \item{\code{COT1}}{Average cardiac output in the tenth minute of rest.}
    \item{\code{COT2}}{Average cardiac output in the second minute during the tilt.}
    \item{\code{COT3}}{Average cardiac output in the fifth minute during the tilt.}
    \item{\code{COT4}}{Average cardiac output in the fifth minute after the tilt.}
    \item{\code{SVRIT1}}{Average systemic vascular resistance index in the tenth minute of rest.}
    \item{\code{SVRIT2}}{Average systemic vascular resistance index in the second minute during the tilt.}
    \item{\code{SVRIT3}}{Average systemic vascular resistance index in the fifth minute during the tilt.}
    \item{\code{SVRIT4}}{Average systemic vascular resistance index in the fifth minute after the tilt.}
    \item{\code{PWVT1}}{Average pulse wave velocity in the tenth minute of rest.}
    \item{\code{PWVT2}}{Average pulse wave velocity in the second minute during the tilt.}
    \item{\code{PWVT3}}{Average pulse wave velocity in the fifth minute during the tilt.}
    \item{\code{PWVT4}}{Average pulse wave velocity in the fifth minute after the tilt.}
    \item{\code{HRT1T2}}{Difference \code{HRT1} - \code{HRT2}.}
    \item{\code{COT1T2}}{Difference \code{COT1} - \code{COT2}.}
    \item{\code{SVRIT1T2}}{Difference \code{SVRIT1} - \code{SVRIT2}.}
    \item{\code{PWVT1T2}}{Difference \code{PWVT1} - \code{PWVT2}.}
    \item{\code{HRT1T4}}{Difference \code{HRT1} - \code{HRT4}.}
    \item{\code{COT1T4}}{Difference \code{COT1} - \code{COT4}.}
    \item{\code{SVRIT1T4}}{Difference \code{SVRIT1} - \code{SVRIT4}.}
    \item{\code{PWVT1T4}}{Difference \code{PWVT1} - \code{PWVT4}.}
  }
}
\details{
This data is a subset of hemodynamic data collected as a part of the LASERI
study (English title: \dQuote{Cardivascular risk in young Finns study}) using
whole-body impedance cardiography and plethysmographic blood pressure
recordings from fingers. The data given here comes from 223 healthy subjects
between 26 and 42 years of age, who participated in the recording of the
hemodynamic variables both in a supine position and during a passive head-up
tilt on a motorized table. During that experiment the subject spent the
first ten minutes in a supine position, then the motorized table was tilted
to a head-up position (60 degrees) for five minutes, and for the last five
minutes the table was again returned to the supine position.

Of interest in this data is for example if the values 5 minutes after the
tilt are already returned to their pre-tilt levels.
}
\source{
Data courtesy of the LASERI study \cr (\url{https://youngfinnsstudy.utu.fi/}).
}

\examples{
# for example testing if the location before the tilt is the same as 
# 5 minutes after the tilt:
data(LASERI)
DIFFS.T1T4 <- subset(LASERI,select=c(HRT1T4,COT1T4,SVRIT1T4))
rank.ctest(DIFFS.T1T4)
rank.ctest(DIFFS.T1T4, score="s")
}
\keyword{datasets}
