#' The MSpecific Class.
#'
#' Class of object returned by function \code{\link[Biosurvmet]{MSpecificCoxPh}}.
#'
#' plot {signature(x = "MSpecific"): Plots for MSpecific class analysis results}
#' signature(x = "MSpecific"): Plots for MSpecific class analysis results.
#'
#' Any parameters of \code{\link[graphics]{plot.default}} may be passed on to this particular plot method.
#'
#' @usage ## S4 method for signature 'MSpecific'
#' plot(x, y, ...)
#' ## S4 method for signature 'MSpecific'
#' summary(MSpecific-object)
#' ## S4 method for signature 'MSpecific'
#' show(MSpecific-object)
#' @name MSpecific-class
#' @rdname MSpecific-class
#' @exportClass MSpecific
#' @param x	 A MSpecific class object
#' @param y	 missing
#' @param ...	The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}
#' @slot Result A list of dataframes of each output object of coxph for the metabolites.
#' @slot HRRG A dataframe with estimated metabolite-specific HR for low risk group and 95 percent CI.
#' @slot Group A matrix of the classification group a subject belongs to for each of the metabolite analysis. The metabolites are on the rows and the subjects are the columns
#' @slot Metnames The names of the metabolites for the analysis
#'
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[Biosurvmet]{MSpecificCoxPh}}
#' @examples
#' ## GENERATE SOME METABOLIC SURVIVAL DATA WITH PROGNOSTIC FACTORS
#' Data<-MSData(nPatients=100,nMet=150,Prop=0.5)
#'
#' ## DO THE METABOLITE BY METABOLITE ANALYSIS
#' Eg = MSpecificCoxPh(Survival=Data$Survival, Mdata=t(Data$Mdata),
#' Censor=Data$Censor, Reduce = FALSE, Select = 15,
#' Prognostic=Data$Prognostic, Quantile = 0.5)
#'
#' ## GET THE CLASS OF THE OBJECT
#' class(Eg)     # An "MSpecific" Class
#'
#' ##  METHOD THAT CAN BE USED FOR THIS CLASS
#' show(Eg)
#' summary(Eg)
#' plot(Eg)

setClass("MSpecific",slots = list(Result="list",HRRG="matrix",Group="matrix",Metnames="vector"),
         prototype=list(Result=list(1),HRRG=matrix(0,0,0),Group=matrix(0,0,0), Metnames = vector()))


#' Method show.
#' @name MSpecific
#' @rdname MSpecific-class
#' @exportMethod show
#setGeneric("show", function(object) standardGeneric("show"))

#' @rdname MSpecific-class
#' @aliases show,MSpecific-method
setMethod("show",signature="MSpecific"
          , function(object){
            cat("Metabolite by Metabolite CoxPh Model\n")
            cat("Number of Metabolite used: ", length(object@Metnames), "\n")
            cat("Number of Univariate coxph: ", length(object@Metnames), "\n")
          })

#' Method summary.
#' @name MSpecific-class
#' @rdname MSpecific-class
#' @exportMethod summary
setGeneric("summary", function(object,...) standardGeneric("summary"))

#' @rdname MSpecific-class
#' @aliases summary,MSpecific-method
setMethod("summary",signature="MSpecific"
          , function(object){
  cat("Summary of Metabolite by Metabolite CoxPh Models\n")
  cat("Number of Metabolites used: ", length(object@Metnames), "\n")
  cat("Top 15 Metabolites out of ", length(object@Metnames), "\n")
  cat("Estimated HR for the low risk group\n")
  # top Metabolites based on upper CI HR GS+
  Names.KMetabolites<-object@Metnames
  index.Top.KMetabolites <- order(object@HRRG[,1],decreasing =FALSE)
  index.Top.KMetabolites<-index.Top.KMetabolites[1:15]
  Top.KMetabolites.GSplus<-data.frame(Metnames=Names.KMetabolites[index.Top.KMetabolites],object@HRRG[index.Top.KMetabolites,])

  colnames(Top.KMetabolites.GSplus)<-c("Metnames","HR GS+","LowerCI","UpperCI")
  # FDR corrected CI for top k  metabolites
  cilevel <- 1-0.05*15/nrow(object@HRRG)

  HRpadj <- exp(log(object@HRRG[index.Top.KMetabolites,1]) + log(object@HRRG[index.Top.KMetabolites,c(2,3)])-log(object@HRRG[index.Top.KMetabolites,1])*qnorm(cilevel)/1.96)  #
  res.topkMetabolites<-data.frame(Top.KMetabolites.GSplus,HRpadj)

  colnames(res.topkMetabolites)<-c("Metnames","HR","LCI","UCI","FDRLCI","FDRUCI")
  print(res.topkMetabolites)
})



#' Method plot.
#' @name MSpecific-class
#' @rdname MSpecific-class
#' @exportMethod plot

#' @rdname MSpecific-class
#' @aliases plot,MSpecific-method
setMethod(f="plot", signature = "MSpecific",
          definition = function(x,y,...){
            object <-  x
            Names.KMetabolites<-object@Metnames
            index.Top.KMetabolites <- order(object@HRRG[,1],decreasing =FALSE)
            Top.KMetabolites.GSplus<-data.frame(Metnames=Names.KMetabolites,object@HRRG)

            colnames(Top.KMetabolites.GSplus)<-c("Metnames","HR GS+","LowerCI","UpperCI")
            # FDR corrected CI for top k  metabolites
            cilevel <- 1-0.05*15/nrow(object@HRRG)

            HRpadj <- exp(log(object@HRRG[index.Top.KMetabolites,1]) + log(object@HRRG[index.Top.KMetabolites,c(2,3)])-log(object@HRRG[index.Top.KMetabolites,1])*qnorm(cilevel)/1.96)  #
            res.topkMetabolites<-data.frame(Top.KMetabolites.GSplus,HRpadj)

            colnames(res.topkMetabolites)<-c("Metnames","HR","LCI","UCI","FDRLCI","FDRUCI")

            x= 1:length(object@Metnames)
            par(mfrow=c(3,1))
            plot(x=x, y = res.topkMetabolites[,2],
                 ylim= c(0,max(res.topkMetabolites[,4])),
                 pch=19, xlab="Metabolites", ylab="Hazard ratio",
                 main="Hazard ratio plot with unadjusted confidence interval")
            arrows(x, res.topkMetabolites[,3], x, res.topkMetabolites[,4], length=0.05, angle=90, code=3)
            abline(h=1,col="red2",lwd=2.0)

            plot(x=x, y = res.topkMetabolites[,2],
                 ylim= c(0,max(res.topkMetabolites[,6])),
                 pch=19, xlab="Metabolites", ylab="Hazard ratio",
                 main="Hazard ratio plot with adjusted confidence interval")
            arrows(x, res.topkMetabolites[,5], x, res.topkMetabolites[,6], length=0.05, angle=90, code=3)
            abline(h=1,col="red2",lwd=2.0)

            Group2 = object@Group
            prop = sapply(1:ncol(Group2),function(k) sum(Group2[,k]=="Low risk")/nrow(Group2))
            barplot(prop, col=rainbow(length(Group2[1,])), xlab = "Patient index", ylab = "Proportion",main = "Proportion of being classified as low risk group",ylim = c(0,max(prop)))
            return(invisible())
          })

