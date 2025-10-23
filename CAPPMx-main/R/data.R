#' Glioblastoma Data Collected by MD Anderson Cancer Center
#' 
#'  Records from 339 highly clinically and molecularly annotated GBM patients
#'  treated at the MD Anderson Cancer Center over more than 10 years comprising  11 clinically important categorical covariates
#'  and overall survival (some are right-censored) endpoints.
#'  
#' To maintain patients' privacy we are providing a resampled version of the original dataset
#' 
#' @name MDACC_reproduced
#' @docType data
#' @format A \code{data.table} with 339 rows and 13 columns:
#' \describe{
#' \item{Surgery Reason}{1=``therapeutic'' or 0=``other'' (relapse)}
#' \item{Histologic Grade}{the grade of astrocytoma: 1=``IV (GBM)'' or 0=``I-III (low-grade or anaplastic)''  }
#' \item{EOR}{the \emph{extent of resection}: 1=``gross-total''; 0=``subtotal or laser interstitial thermal therapy''}
#' \item{Gender}{1=Male or 0=Female}
#' \item{ATRX}{indicates loss of the ATRX chromatin remodeler gene: 1=Yes;  0=No}
#' \item{MGMT}{status of MGMT (\eqn{O^6}-methylguanine-DNA methyltransferase) gene: 1=Methylated;  0= Methylated or Uniterpretable}
#' \item{CT}{indicates whether the patient participated in a therapeutic trial: 1=Yes; 0=No}
#' \item{SOC}{indicates whether the patient received standard-of-care (concurrent radiation therapy and temozolomide): 1=Yes; 0=No}
#' \item{RT Dose}{radiation therapy dose 1=\eqn{\geq 50} Gray; 0=\eqn{< 50} Gray}
#' \item{KPS}{Karnofsky performance score 1=\eqn{> 60}; 0=\eqn{\leq 60}}
#' \item{Age}{ 1=\eqn{> 55} years; 0=\eqn{\leq 55} years}
#' \item{endpts}{Recorded overall survival times (OS) in weeks}
#' \item{surv_inds}{A logical vector indicates whether the corresponding recorded OS is a ovserved failue or not}
#' }
#' 
#' @source {MD Anderson Cancer Center}
#' @keywords Glioblastoma
#' 
#' @examples
#' data(MDACC_reproduced)
#' head(MDACC_reproduced)
"MDACC_reproduced"

#' Covariate names in the \code{MDACC_reproduced} dataset
#' 
#' @format A string vector with the covariate names in the \code{\link{MDACC_reproduced}} dataset
#' 
#' @examples
#' data(imp.variables)
#' imp.variables
"imp.variables"