##### entropart Methods to measure the diversity of spatialized communities ####


#' Tsallis (HCDT) Entropy of a spatialized community
#'
#' Calculates the HCDT, also known as Tsallis entropy of order \eqn{q} of a probability vector.
#'
#' Tsallis (Havrda and Charvat, 1967; Daroczy, 1970; Tsallis, 1988) generalized entropy is a generalized measure of diversity (Jost, 2006).
#' See \code{\link{Tsallis}} for more details.
#'
#' @param NorP An object of class "wmppp" (\code{\link{wmppp.object}}), with \code{PointType} values as species names.
#' @param q A number: the order of entropy. Some corrections allow only a positive number. Default is 1 for Shannon entropy.
#' @param Correction A string containing one of the possible corrections: \code{"None"} (no correction), \code{"ChaoShen"}, \code{"GenCov"}, \code{"Grassberger"}, \code{"Holste"}, \code{"Bonachela"}, \code{"ZhangGrabchak"}, or \code{"ChaoWangJost"}, \code{"Marcon"}, \code{"UnveilC"}, \code{"UnveiliC"}, \code{"UnveilJ"} or \code{"Best"}, the default value.  Currently, \code{"Best"} is \code{"ChaoWangJost"}.
#' @param Ps Ignored
#' @param Ns Ignored
#' @inheritParams as.ProbaVector.wmppp
#'
#' @return A named number equal to the calculated entropy. The name is that of the bias correction used.
#'
#' @export
#' @importFrom entropart Tsallis
#' @method Tsallis wmppp
#' @importFrom entropart bcTsallis
#'
#' @examples #TODO
Tsallis.wmppp <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  return(entropart::bcTsallis(Ns=as.AbdVector(NorP), q=q, Correction=Correction, CheckArguments=CheckArguments))
}



#' Tsallis (HCDT) Entropy of a vector of individuals
#'
#' Calculates the HCDT, also known as Tsallis entropy of order \eqn{q} of a probability vector.
#'
#' Tsallis (Havrda and Charvat, 1967; Daroczy, 1970; Tsallis, 1988) generalized entropy is a generalized measure of diversity (Jost, 2006).
#' See \code{\link{Tsallis}} for more details.
#'
#' @param NorP A vector of factors.
#' @inheritParams Tsallis.wmppp
#'
#' @return A named number equal to the calculated entropy. The name is that of the bias correction used.
#'
#' @export
#' @importFrom entropart Tsallis
#' @method Tsallis factor
#' @importFrom entropart bcTsallis
#'
#' @examples #TODO
Tsallis.factor <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  return(entropart::bcTsallis(Ns=table(NorP), q=q, Correction=Correction, CheckArguments=CheckArguments))
}



#' Tsallis (HCDT) Entropy of a vector of individuals
#'
#' Calculates the HCDT, also known as Tsallis entropy of order \eqn{q} of a probability vector.
#'
#' Tsallis (Havrda and Charvat, 1967; Daroczy, 1970; Tsallis, 1988) generalized entropy is a generalized measure of diversity (Jost, 2006).
#' See \code{\link{Tsallis}} for more details.
#'
#' @param NorP A vector of characters.
#' @inheritParams Tsallis.wmppp
#'
#' @return A named number equal to the calculated entropy. The name is that of the bias correction used.
#'
#' @export
#' @importFrom entropart Tsallis
#' @method Tsallis character
#' @importFrom entropart bcTsallis
#'
#' @examples #TODO
Tsallis.character <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  return(entropart::bcTsallis(Ns=table(NorP), q=q, Correction=Correction, CheckArguments=CheckArguments))
}


#' HCDT diversity of a spatialized community
#'
#' \code{Diversity} calls \code{\link{Tsallis}} to calculate entropy and transforms it into diversity by calculating its deformed exponential.
#'
#' @inheritParams Tsallis.wmppp
#'
#' @return A named number equal to the calculated diversity. The name is that of the bias correction used.
#'
#' @export
#' @importFrom entropart Diversity
#' @method Diversity wmppp
#' @importFrom entropart bcDiversity
#'
#' @examples #TODO
Diversity.wmppp <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  return (entropart::bcDiversity(Ns=as.AbdVector(NorP), q=q, Correction=Correction, CheckArguments=CheckArguments))
}



#' HCDT diversity of a vector of individuals
#'
#' \code{Diversity} calls \code{\link{Tsallis}} to calculate entropy and transforms it into diversity by calculating its deformed exponential.
#'
#' @inheritParams Tsallis.factor
#'
#' @return A named number equal to the calculated diversity. The name is that of the bias correction used.
#'
#' @export
#' @importFrom entropart Diversity
#' @method Diversity factor
#' @importFrom entropart bcDiversity
#'
#' @examples #TODO
Diversity.factor <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  return (entropart::bcDiversity(Ns=table(NorP), q=q, Correction=Correction, CheckArguments=CheckArguments))
}



#' HCDT diversity of a vector of individuals
#'
#' \code{Diversity} calls \code{\link{Tsallis}} to calculate entropy and transforms it into diversity by calculating its deformed exponential.
#'
#' @inheritParams Tsallis.factor
#'
#' @return A named number equal to the calculated diversity. The name is that of the bias correction used.
#'
#' @export
#' @importFrom entropart Diversity
#' @method Diversity character
#' @importFrom entropart bcDiversity
#'
#' @examples #TODO
Diversity.character <-
function(NorP, q = 1, Correction = "Best", CheckArguments = TRUE, Ps = NULL, Ns = NULL)
{
  return (entropart::bcDiversity(Ns=table(NorP), q=q, Correction=Correction, CheckArguments=CheckArguments))
}
