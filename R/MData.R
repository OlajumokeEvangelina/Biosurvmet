#' Metabolic Data
#'
#' This is a dataset containing the metabolite information of 149 subjects. The data contains 110 different metabolic profiles
#'
#' @format A data frame with 149 rows and 110 metabolites, the rows are the subject while the columns are the metabolites represented as;
#' \describe{
#' \item{var1}{The first metabolite }
#'   \item{var2}{The second metabolite}
#'   \item{var3}{The third metabolite}
#'   \item{var4}{The fourth metabolite}
#'   ...
#'   ...
#'   ...
#'   \item{var110}{The one hundred and tenth metabolite}
#' }
#' @keywords datasets
#' @importFrom Rdpack reprompt
#' @usage data(MData)
#' @source \url{https://bmccancer.biomedcentral.com/articles/10.1186/s12885-018-4755-1}
#' @references
#' \insertRef{yeye}{Biosurvmet}
#'
#' \insertRef{yep}{Biosurvmet}
#' @examples
#' data(MData)
# # summary statistics of the first 10 metabolites
#' summary(MData[,1:10])
#' # summary statistics of all the metabolites
#' summary(MData[,])
"MData"
