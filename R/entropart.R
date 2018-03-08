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
#' @param ... Further arguments. Unsused.
#' @inheritParams EntAccum
#'
#' @return A named number equal to the calculated entropy. The name is that of the bias correction used.
#'
#' @importFrom entropart Tsallis
#' @method Tsallis wmppp
#' @export
#'
#' @examples #TODO
Tsallis.wmppp <-
function(NorP, q = 1, Correction = "Best", ..., CheckArguments = TRUE)
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
#' @importFrom entropart Tsallis
#' @method Tsallis factor
#' @export
#'
#' @examples #TODO
Tsallis.factor <-
function(NorP, q = 1, Correction = "Best", ..., CheckArguments = TRUE)
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
#' @importFrom entropart Tsallis
#' @method Tsallis character
#' @export
#'
#' @examples #TODO
Tsallis.character <-
function(NorP, q = 1, Correction = "Best", ..., CheckArguments = TRUE)
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
#' @importFrom entropart Diversity
#' @method Diversity wmppp
#' @export
#'
#' @examples #TODO
Diversity.wmppp <-
function(NorP, q = 1, Correction = "Best", ..., CheckArguments = TRUE)
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
#' @importFrom entropart Diversity
#' @method Diversity factor
#' @export
#'
#' @examples #TODO
Diversity.factor <-
function(NorP, q = 1, Correction = "Best", ..., CheckArguments = TRUE)
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
#' @importFrom entropart Diversity
#' @method Diversity character
#' @export
#'
#' @examples #TODO
Diversity.character <-
function(NorP, q = 1, Correction = "Best", ..., CheckArguments = TRUE)
{
  return (entropart::bcDiversity(Ns=table(NorP), q=q, Correction=Correction, CheckArguments=CheckArguments))
}



#' Richness of a spatialized community
#'
#' \code{Richness} is the number of species.
#'
#' @inheritParams Tsallis.wmppp
#' @param Correction A string containing one of the possible corrections: \code{"None"} (no correction), \code{"Jackknife"}, \code{"iChao1"}, or \code{"Chao1"}, the default value.
#' @param Alpha The risk level, 5\% by default, used to optimize the jackknife order.
#' @param JackOver If \code{TRUE}, retain the jackknife order immediately superior to the optimal one, usually resulting in the overestimation of the number of species. Default is \code{FALSE}.
#' 
#' @return A named number equal to the calculated diversity. The name is that of the bias correction used.
#'
#' @importFrom entropart Richness
#' @method Richness wmppp
#' @export
#'
#' @examples #TODO
Richness.wmppp <-
function(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE,  ..., CheckArguments = TRUE)
{
  return (entropart::bcRichness(Ns=as.AbdVector(NorP), Correction=Correction, Alpha=Alpha, JackOver=JackOver, CheckArguments=CheckArguments))
}



#' Richness of a vector of individuals
#'
#' \code{Richness} is the number of species.
#' 
#' @inheritParams Richness.wmppp
#'
#' @return A named number equal to the calculated diversity. The name is that of the bias correction used.
#'
#' @importFrom entropart Richness
#' @method Richness factor
#' @export
#'
#' @examples #TODO
Richness.factor <-
function(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE,  ..., CheckArguments = TRUE)
{
  return (entropart::bcRichness(Ns=table(NorP), Correction=Correction, Alpha=Alpha, JackOver=JackOver, CheckArguments=CheckArguments))
}



#' Richness of a vector of individuals
#'
#' \code{Richness} is the number of species.
#'
#' @inheritParams Tsallis.factor
#'
#' @return A named number equal to the calculated diversity. The name is that of the bias correction used.
#'
#' @importFrom entropart Richness
#' @method Richness character
#' @export
#'
#' @examples #TODO
Richness.character <-
function(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE,  ..., CheckArguments = TRUE)
{
  return (entropart::bcRichness(Ns=table(NorP), Correction=Correction, Alpha=Alpha, JackOver=JackOver, CheckArguments=CheckArguments))
}


