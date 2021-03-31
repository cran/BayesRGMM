#' The German socioeconomic panel study data
#'
#' The German socioeconomic panel study data was taken from the first twelve annual waves (1984 through 1995) 
#' of the German Socioeconomic Panel (GSOEP) which surveys a representative sample of East and West German households. 
#' The data provide detailed information on the utilization of health care facilities, characteristics of current 
#' employment, and the insurance schemes under which individuals are covered. We consider the sample of individuals 
#' aged 25 through 65 from the West German subsample and of German nationality. 
#' The sample contained 3691 male and 3689 female individuals which make up a sample of 14,243 male and 13,794 female 
#' person-year observations. 
#'
#' @docType data
#'
#' @usage data(GSPS)
#'
#' @format A data frame with 27326 rows and 25 variables
#' \describe{
#' \item{id}{	     person - identification number}
#' \item{female}{	     female = 1; male = 0}
#' \item{year}{	     calendar year of the observation}
#' \item{age}{	     age in years}
#' \item{hsat}{	     health satisfaction, coded 0 (low) - 10 (high)}
#' \item{handdum}{	     handicapped = 1; otherwise = 0}
#' \item{handper}{	     degree of handicap in percent (0 - 100)}
#' \item{hhninc}{	     household nominal monthly net income in German marks / 1000}
#' \item{hhkids}{	     children under age 16 in the household = 1; otherwise = 0}
#' \item{educ}{	     years of schooling}
#' \item{married}{	     married = 1; otherwise = 0}
#' \item{haupts}{	     highest schooling degree is Hauptschul degree = 1; otherwise = 0}
#' \item{reals}{	     highest schooling degree is Realschul degree = 1; otherwise = 0}
#' \item{fachhs}{	     highest schooling degree is Polytechnical degree = 1; otherwise = 0}
#' \item{abitur}{	     highest schooling degree is Abitur = 1; otherwise = 0}
#' \item{univ}{	     highest schooling degree is university degree = 1; otherwise = 0}
#' \item{working}{	     employed = 1; otherwise = 0}
#' \item{bluec}{	     blue collar employee = 1; otherwise = 0}
#' \item{whitec}{	     white collar employee = 1; otherwise = 0}
#' \item{self}{	     self employed = 1; otherwise = 0}
#' \item{beamt}{	     civil servant = 1; otherwise = 0}
#' \item{docvis}{	     number of doctor visits in last three months}
#' \item{hospvis}{	     number of hospital visits in last calendar year}
#' \item{public}{	     insured in public health insurance = 1; otherwise = 0}
#' \item{addon}{	     insured by add-on insurance = 1; otherswise = 0}
#' }
#'
#' @keywords datasets
#'
#' @references{
#'   \insertRef{Riphahn:etal:2003}{BayesRGMM} 
#'}
#' @source \href{http://qed.econ.queensu.ca/jae/2003-v18.4/riphahn-wambach-million/}{JAE Archive}
"GSPS"