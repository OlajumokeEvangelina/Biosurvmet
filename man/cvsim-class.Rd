% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cvsim-class.R
\docType{class}
\name{cvsim-class}
\alias{cvsim-class}
\alias{cvsim}
\alias{show,cvsim-method}
\alias{summary,cvsim-method}
\alias{plot,cvsim,missing-method}
\title{The cvsim Class.}
\usage{
\S4method{show}{cvsim}(object)

\S4method{summary}{cvsim}(object)

\S4method{plot}{cvsim,missing}(x, y, type = 1, ...)
}
\arguments{
\item{object}{A cvsim class object}

\item{x}{A cvsim class object}

\item{y}{missing}

\item{type}{Plot type. 1 distribution of the HR under test For the Top K metabolites using PCA. 2 distribution of the HR under test For the Top K metabolites using PLS.}

\item{...}{The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}}
}
\description{
Class of object returned by function \code{\link[MetabolicSurv]{CVSim}}.
}
\section{Slots}{

\describe{
\item{\code{HRpca}}{A 3-way array in which first, second, and third dimensions correspond to number of metabolites, Hazard ratio infromation(Estimated HR, LowerCI and UpperCI), and number of cross validation respectively. This contains the estimated HR on test data and dimension reduction method is PCA.}

\item{\code{HRpls}}{A 3-way array in which first, second, and third dimensions correspond to number of metabolites, Hazard ratio infromation(Estimated HR, LowerCI and UpperCI), and number of cross validation respectively. This contains the estimated HR on test data and dimension reduction method is PLS.}

\item{\code{Nmets}}{The number of metabolites in the reduced matrix}

\item{\code{Ncv}}{The number of cross validation done}

\item{\code{Top}}{A sequence of top k metabolites considered. Default is Top=seq(5,100,by=5)}
}}

\examples{
\donttest{
## GENERATE SOME METABOLIC SURVIVAL DATA WITH PROGNOSTIC FACTORS
Data<-MSData(nPatients=100,nMet=150,Prop=0.5)

## FIRST IS THE NETABOLITE BY METABOLITE ANALYSIS
w = CVMetSpecificCoxPh(Fold=3,Survival=Data$Survival,
Mdata=t(Data$Mdata),Censor= Data$Censor,Reduce=TRUE,
Select=150,Prognostic=Data$Prognostic,Quantile = 0.5,Ncv=3)

## USING THE FUNCTION
Result = CVSim(w, Top = seq(5, 100, by = 5), Survival=Data$Survival,
 Censor=Data$Censor, Prognostic = Data$Prognostic)

## GET THE CLASS OF THE OBJECT
class(Result)     # A "cvsim" Class

##  METHOD THAT CAN BE USED FOR THE RESULT
show(Result)
summary(Result)
plot(Result, type =2)
}
}
\seealso{
\code{\link[MetabolicSurv]{CVPcaPls}}, \code{\link[MetabolicSurv]{SurvPcaClass}}, \code{\link[MetabolicSurv]{SurvPlsClass}}
}
\author{
Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@uhasselt.be}

Ziv Shkedy
}
