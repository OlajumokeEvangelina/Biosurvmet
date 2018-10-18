#'Cross validation for majority votes
#'
#' This function does cross validation for the Majority votes based classification.
#'
#' This function does cross validation for the Majority votes based classification which is a cross validated approach to \code{\link[Biosurvmet]{MajorityVotes}}.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Censor A vector of censoring indicator
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Mdata A large or small metabolic profile matrix. A matrix with metabolic profiles where the number of rows should be equal to the number of metabolites and number of columns should be equal to number of patients.
#' @param Reduce A boolean parameter indicating if the metabolic profile matrix should be reduced, default is TRUE and larger metabolic profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of metabolites (default is 15) to be selected from supervised PCA. This is valid only if th argument Reduce=TRUE
#' @param Fold Number of times in which the dataset is divided. Default is 3 which implies dataset will be divided into three groups and 2/3 of the dataset will be the train datset and 1/3 will be to train the results.
#' @param Ncv The Number of cross validation loop. Default is 50 but it is recommended to have at least 100.
#' @return A object of class \code{\link[Biosurvmet]{cvmv}} is returned with the following values
#'   \item{HRTrain}{A matrix of survival information for the training dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.}
#'   \item{HRTest}{A matrix of survival information for the test dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.}
#'   \item{Ncv}{The number of cross validation used}
#'  \item{Mdata}{The Metabolite data matrix that was used for the analysis either same as Mdata or a reduced version.}
#'   \item{Progfact}{The names of prognostic factors used}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[Biosurvmet]{Majorityvotes}}
#' @examples
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#' Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## USING THE FUNCTION
#' Result = CVMajorityVotes(Survival=Data$Survival,Censor=Data$Censor,
#' Prognostic=Data$Prognostic, Mdata=t(Data$Mdata), Reduce=FALSE, S
#' elect=15, Fold=3, Ncv=10)
#'
#' ## GET THE CLASS OF THE OBJECT
#' class(Result)     # An "cvpp" Class
#'
#' ##  METHOD THAT CAN BE USED FOR THE RESULT
#' show(Result)
#' summary(Result)
#' plot(Result)
#' @export CVMajorityVotes

CVMajorityVotes<-function(Survival,Censor, Prognostic=NULL, Mdata, Reduce=TRUE, Select=150, Fold=3, Ncv=100){

  options( warn = -1)


  if (Reduce) {
    DataForReduction<-list(x=Mdata,y=Survival, censoring.status=Censor, featurenames=rownames(Mdata))
    TentativeList<-names(sort(abs(superpc.train(DataForReduction, type="survival")$feature.scores),decreasing =TRUE))[1:Select]
    TentativeList

    ReduMdata<-Mdata[TentativeList,]
  } else {
    ReduMdata<-Mdata
  }

  n.mets<-nrow(ReduMdata)
  n.patients<-ncol(ReduMdata)
  n.train<-(n.patients-floor(n.patients/Fold))
  n.test<-floor(n.patients/Fold)
  ind.train <-matrix(0,Ncv,n.train)
  ind.test  <-matrix(0,Ncv,n.test)
  res <-  vector("list", n.mets)


  HRp.train <- matrix(0,Ncv,3)  # Training
  HRp.test <-  matrix(0,Ncv,3)     # TTsting

  set.seed(123)
  pIndex <- c(1:n.patients)
  res <-res1<-res2<-res3<- vector("list", Ncv)

  for (j in 1:Ncv){
    message('Cross validation loop ',j)
    p1<-NA
    p2<-NA
    gr.train <- matrix(0, n.mets,ncol(ind.train))
    gr.test  <- matrix(0,n.mets,ncol(ind.test))

    ind.train[j,] <-sort(sample(pIndex,n.train,replace=F) )
    ind.test[j,] <-c(1:n.patients)[-c(intersect(ind.train[j,] ,c(1:n.patients)))]


    for (i in 1:n.mets){
      genei <- ReduMdata[i,ind.train[j,]]

      if (is.null(Prognostic)) {

        cdata <- data.frame(Survival=Survival[ind.train[j,]],Censor=Censor[ind.train[j,]],genei)
        m0 <- coxph(Surv(Survival, Censor==1) ~ genei,data=cdata)
      }
      if (!is.null(Prognostic)) {
        if (is.data.frame(Prognostic)) {
          nPrgFac<-ncol(Prognostic)
          cdata <- data.frame(Survival=Survival[ind.train[j,]],Censor=Censor[ind.train[j,]],genei,Prognostic[ind.train[j,],])
          NameProg<-colnames(Prognostic)
          eval(parse(text=paste( "m0 <-coxph(Surv(Survival, Censor==1) ~ genei",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
        } else {

          stop(" Argument 'Prognostic' is NOT a data frame ")
        }

      }
      #risk Score
      TrtandPC1<-summary(m0)[[7]][c("genei"),1]
      p1 <- TrtandPC1*genei

      TempGenei <-EstimateHR(p1,Data.Survival=cdata, Prognostic=Prognostic[ind.train[j,],],Plots = FALSE, Quantile = 0.5 )
      gr.train[i,]<-TempGenei$Riskgroup


      #---------------------------------------  Testing  Set -------------------------------------------
      geneit <- ReduMdata[i,ind.test[j,]]
      if (!is.null(Prognostic)) {
        if (is.data.frame(Prognostic)) {
          nPrgFac<-ncol(Prognostic)
          cdata <- data.frame(Survival=Survival[ind.test[j,]],Censor=Censor[ind.test[j,]],geneit,Prognostic[ind.test[j,],])
          NameProg<-colnames(Prognostic)
          eval(parse(text=paste( "m0 <-coxph(Surv(Survival, Censor==1) ~ geneit",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
        } else {

          stop(" Argument 'Prognostic' is NOT a data frame ")
        }

      }
      TrtandPC1<-summary(m0)[[7]][c("geneit"),1]
      p2 <- TrtandPC1*geneit

      TempGeneit <-EstimateHR(p2,Data.Survival=cdata, Prognostic=Prognostic[ind.test[j,],],Plots = FALSE, Quantile = 0.5 )
      gr.test[i,]<-TempGeneit$Riskgroup

      #mp1 <- median(p1)
      #p2 <- TrtandPC1*geneit



      #gr.test[i,] <- gs
    } # END OF LOOP over Metabolites---------------------------------------------------------------------

    # ------------ count majority votes for jth Cross validation and estimate HR --------------

    ggr.train <- per.R<-per.NR <- NULL
    for (k in 1:ncol(ind.train)){
      per.R[k]<-sum(gr.train[,k]=="Low risk")
      ggr.train[k]<-ifelse((n.mets-per.R[k])>per.R[k],"High risk","Low risk")
    }



    #-------------------- HR estimation for Training  ----------------------
    GS<-as.factor(ggr.train)
    if (is.null(Prognostic)) {

      cdata <- data.frame(Survival=Survival[ind.train[j,]],Censor=Censor[ind.train[j,]],GS)
      mTrain <- coxph(Surv(Survival, Censor==1) ~ GS,data=cdata)
    }
    if (!is.null(Prognostic)) {
      if (is.data.frame(Prognostic)) {
        nPrgFac<-ncol(Prognostic)
        cdata <- data.frame(Survival=Survival[ind.train[j,]],Censor=Censor[ind.train[j,]],GS,Prognostic[ind.train[j,],])
        NameProg<-colnames(Prognostic)
        eval(parse(text=paste( "mTrain <-coxph(Surv(Survival, Censor==1) ~ GS",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
      } else {
        stop(" Argument 'Prognostic' is NOT a data frame ")
      }

    }
    HRp.train[j,]<-(summary(mTrain)[[8]][1,])[-2]

    #-------------------- HR estimation for Testing  ----------------------
    ggr.test <- per.R<-per.NR <- NULL
    for (k in 1:ncol(ind.test)){
      per.R[k]<-sum(gr.test[,k]=="Low risk")
      ggr.test[k]<-ifelse((n.mets-per.R[k])>per.R[k],"High risk","Low risk")
    }

    GS<-as.factor(ggr.test)
    if (is.null(Prognostic)) {

      cdata <- data.frame(Survival=Survival[ind.test[j,]],Censor=Censor[ind.test[j,]],GS)
      mTest <- coxph(Surv(Survival, Censor==1) ~ GS,data=cdata)
    }
    if (!is.null(Prognostic)) {
      if (is.data.frame(Prognostic)) {
        nPrgFac<-ncol(Prognostic)
        cdata <- data.frame(Survival=Survival[ind.test[j,]],Censor=Censor[ind.test[j,]],GS,Prognostic[ind.test[j,],])
        NameProg<-colnames(Prognostic)
        eval(parse(text=paste( "mTest <-coxph(Surv(Survival, Censor==1) ~ GS",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
      } else {
        stop(" Argument 'Prognostic' is NOT a data frame ")
      }

    }
    HRp.test[j,]<-(summary(mTest)[[8]][1,])[-2]

  }#---------------------------  END OF  FOR LOOP over Cross Validations ------------------------

  pFactors<-NA
  if (!is.null(Prognostic)) pFactors <-colnames(Prognostic)

  return(new("cvmv",HRTrain=HRp.train,HRTest=HRp.test,Ncv=Ncv,Mdata=ReduMdata, Progfact=pFactors))
  }
