% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perm-class.R
\docType{class}
\name{perm-class}
\alias{perm-class}
\alias{perm}
\alias{show,perm-method}
\alias{summary,perm-method}
\alias{plot,perm,ANY-method}
\title{The perm Class.}
\usage{
\S4method{show}{perm}(object)

\S4method{summary}{perm}(object)

\S4method{plot}{perm,ANY}(x, y, ...)
}
\arguments{
\item{object}{A perm class object}

\item{x}{A perm class object}

\item{y}{missing}

\item{...}{The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}}
}
\description{
Class of object returned by function \code{\link[MetabolicSurv]{DistHR}}.
}
\section{Slots}{

\describe{
\item{\code{HRobs}}{Estimated HR for low risk group on the original data.}

\item{\code{HRperm}}{Estimated HR for low risk group on the permuted data}

\item{\code{nperm}}{Number of permutations carried out.}

\item{\code{Validation}}{The validation scheme that was used.}
}}

\note{
The first, third and last vertical line on the plot are the lower,
median  and upper CI of the permuted data estimated HR while
the red line is the estimated HR of the original data
}
\examples{
\donttest{
## GENERATE SOME METABOLIC SURVIVAL DATA WITH PROGNOSTIC FACTORS
Data<-MSData(nPatients=100,nMet=150,Prop=0.5)

## USING THE FUNCTION
Example <- DistHR(Survival = Data$Survival,Mdata = t(Data$Mdata),
Censor = Data$Censor,Reduce=FALSE,Select=15,Prognostic=Data$Prognostic,
Quantile = 0.5, nperm=10, case=2, Validation=c("L1based"))

## GET THE CLASS OF THE OBJECT
class(Example)     # A "perm" Class

##  METHOD THAT CAN BE USED FOR THIS CLASS
show(Example)
summary(Example)
plot(Example)
}
}
\seealso{
\code{\link[MetabolicSurv]{DistHR}}, \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[MetabolicSurv]{SurvPcaClass}}, \code{\link[MetabolicSurv]{SurvPlsClass}}, \code{\link[MetabolicSurv]{Majorityvotes}}, \code{\link[MetabolicSurv]{Lasoelacox}}, \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[MetabolicSurv]{Lasoelacox}}
}
\author{
Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@uhasselt.be}

Ziv Shkedy
}
