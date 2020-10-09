##### entropart Methods to measure the diversity of spatialized communities ####


#' Tsallis (HCDT) Entropy of a spatialized community
#'
#' Calculates the HCDT, also known as Tsallis entropy of order $q$ of a spatialized community.
#'
#' Tsallis \insertCite{Havrda1967,Daroczy1970,Tsallis1988}{SpatDiv} generalized entropy is a generalized measure of diversity \insertCite{Jost2006}{SpatDiv}.
#' See [Tsallis] for more details.
#'
#' @param NorP An object of class [dbmss::wmppp], with `PointType` values as species names.
#' @param q A number: the order of entropy. Some corrections allow only a positive number. Default is 1 for Shannon entropy.
#' @param Correction A string containing one of the possible corrections: `"None"` (no correction), `"ChaoShen"`, `"GenCov"`, `"Grassberger"`, `"Holste"`, `"Bonachela"`, `"ZhangGrabchak"`, or `"ChaoWangJost"`, `"Marcon"`, `"UnveilC"`, `"UnveiliC"`, `"UnveilJ"` or `"Best"`, the default value.  Currently, `"Best"` is `"ChaoWangJost"`.
#' @param ... Further arguments. Unsused.
#' @inheritParams EntAccum
#'
#' @return A named number equal to the calculated entropy. The name is that of the bias correction used.
#'
#' @importFrom entropart Tsallis
#' @method Tsallis wmppp
#' @export
#'
#' @references
#' \insertAllCited{}
#' @examples
#' Tsallis(Paracou6)
#' 
Tsallis.wmppp <-
function(NorP, q = 1, Correction = "Best", ..., CheckArguments = TRUE)
{
  return(entropart::bcTsallis(Ns=as.AbdVector(NorP), q=q, Correction=Correction, CheckArguments=CheckArguments))
}



#' Tsallis (HCDT) Entropy of a vector of individuals
#'
#' Calculates the HCDT, also known as Tsallis entropy of order $q$ of a vector of individuals
#'
#' Tsallis \insertCite{Havrda1967,Daroczy1970,Tsallis1988}{SpatDiv} generalized entropy is a generalized measure of diversity \insertCite{Jost2006}{SpatDiv}.
#' See [Tsallis] for more details.
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
#' @references
#' \insertAllCited{}
#' @examples
#' Tsallis(Paracou6$marks$PointType)
#' 
Tsallis.factor <-
function(NorP, q = 1, Correction = "Best", ..., CheckArguments = TRUE)
{
  return(entropart::bcTsallis(Ns=as.numeric(table(NorP)), q=q, Correction=Correction, CheckArguments=CheckArguments))
}



#' Tsallis (HCDT) Entropy of a vector of individuals
#'
#' Calculates the HCDT, also known as Tsallis entropy of order $q$ of a vector of individuals.
#'
#' Tsallis \insertCite{Havrda1967,Daroczy1970,Tsallis1988}{SpatDiv} generalized entropy is a generalized measure of diversity \insertCite{Jost2006}{SpatDiv}.
#' See [Tsallis] for more details.
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
#' @references
#' \insertAllCited{}
#' @examples
#' SpeciesList <- as.character(Paracou6$marks$PointType)
#' Tsallis(SpeciesList)
#' 
Tsallis.character <-
function(NorP, q = 1, Correction = "Best", ..., CheckArguments = TRUE)
{
  return(entropart::bcTsallis(Ns=as.numeric(table(NorP)), q=q, Correction=Correction, CheckArguments=CheckArguments))
}



#' HCDT diversity of a spatialized community
#'
#' `Diversity` calls [Tsallis] to calculate entropy and transforms it into diversity by calculating its deformed exponential.
#'
#' @inheritParams Tsallis.wmppp
#'
#' @return A named number equal to the calculated diversity. The name is that of the bias correction used.
#'
#' @importFrom entropart Diversity
#' @method Diversity wmppp
#' @export
#'
#' @examples
#' Diversity(Paracou6)
#' 
Diversity.wmppp <-
function(NorP, q = 1, Correction = "Best", ..., CheckArguments = TRUE)
{
  return (entropart::bcDiversity(Ns=as.AbdVector(NorP), q=q, Correction=Correction, CheckArguments=CheckArguments))
}



#' HCDT diversity of a vector of individuals
#'
#' `Diversity` calls [Tsallis]  to calculate entropy and transforms it into diversity by calculating its deformed exponential.
#'
#' @inheritParams Tsallis.factor
#'
#' @return A named number equal to the calculated diversity. The name is that of the bias correction used.
#'
#' @importFrom entropart Diversity
#' @method Diversity factor
#' @export
#'
#' @examples
#' Diversity(Paracou6$marks$PointType)
#' 
Diversity.factor <-
function(NorP, q = 1, Correction = "Best", ..., CheckArguments = TRUE)
{
  return (entropart::bcDiversity(Ns=as.numeric(table(NorP)), q=q, Correction=Correction, CheckArguments=CheckArguments))
}



#' HCDT diversity of a vector of individuals
#'
#' `Diversity` calls [Tsallis]  to calculate entropy and transforms it into diversity by calculating its deformed exponential.
#'
#' @inheritParams Tsallis.factor
#'
#' @return A named number equal to the calculated diversity. The name is that of the bias correction used.
#'
#' @importFrom entropart Diversity
#' @method Diversity character
#' @export
#'
#' @examples
#' SpeciesList <- as.character(Paracou6$marks$PointType)
#' Diversity(SpeciesList)
#' 
Diversity.character <-
function(NorP, q = 1, Correction = "Best", ..., CheckArguments = TRUE)
{
  return (entropart::bcDiversity(Ns=as.numeric(table(NorP)), q=q, Correction=Correction, CheckArguments=CheckArguments))
}



#' Richness of a spatialized community
#'
#' `Richness` is the number of species.
#'
#' @inheritParams Tsallis.wmppp
#' @param Correction A string containing one of the possible corrections: `"None"` (no correction), `"Jackknife"`, `"iChao1"`, or `"Chao1"`, the default value.
#' @param Alpha The risk level, 5\% by default, used to optimize the jackknife order.
#' @param JackOver If `TRUE`, retain the jackknife order immediately superior to the optimal one, usually resulting in the overestimation of the number of species. Default is `FALSE`.
#' 
#' @return A named number equal to the calculated diversity. The name is that of the bias correction used.
#'
#' @importFrom entropart Richness
#' @method Richness wmppp
#' @export
#'
#' @examples
#' Richness(Paracou6, Correction="Jackknife")
#' 
Richness.wmppp <-
function(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE,  ..., CheckArguments = TRUE)
{
  return (entropart::bcRichness(Ns=as.AbdVector(NorP), Correction=Correction, Alpha=Alpha, JackOver=JackOver, CheckArguments=CheckArguments))
}



#' Richness of a vector of individuals
#'
#' `Richness` is the number of species.
#' 
#' @inheritParams Richness.wmppp
#'
#' @return A named number equal to the calculated diversity. The name is that of the bias correction used.
#'
#' @importFrom entropart Richness
#' @method Richness factor
#' @export
#'
#' @examples
#' Richness(Paracou6$marks$PointType, Correction="Jackknife")
#' 
Richness.factor <-
function(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE,  ..., CheckArguments = TRUE)
{
  return (entropart::bcRichness(Ns=as.numeric(table(NorP)), Correction=Correction, Alpha=Alpha, JackOver=JackOver, CheckArguments=CheckArguments))
}



#' Richness of a vector of individuals
#'
#' `Richness` is the number of species.
#'
#' @inheritParams Richness.wmppp
#'
#' @return A named number equal to the calculated diversity. The name is that of the bias correction used.
#'
#' @importFrom entropart Richness
#' @method Richness character
#' @export
#'
#' @examples
#' SpeciesList <- as.character(Paracou6$marks$PointType)
#' Richness(SpeciesList, Correction="Jackknife")
#' 
Richness.character <-
function(NorP, Correction = "Chao1", Alpha = 0.05, JackOver = FALSE,  ..., CheckArguments = TRUE)
{
  return (entropart::bcRichness(Ns=as.numeric(table(NorP)), Correction=Correction, Alpha=Alpha, JackOver=JackOver, CheckArguments=CheckArguments))
}



#' Shannon Entropy of a spatialized community
#'
#' Calculates the Shannon entropy of a probability vector.
#'
#'
#' @param NorP An object of class [dbmss::wmppp], with `PointType` values as species names.
#' @param Correction A string containing one of the possible corrections: see [Shannon].
#' 
#' @inheritParams Tsallis.wmppp
#'
#' @return A named number equal to the calculated entropy. The name is that of the bias correction used.
#'
#' @importFrom entropart Shannon
#' @method Shannon wmppp
#' @export
#'
#' @examples
#' Shannon(Paracou6)
#'  
Shannon.wmppp <-
function(NorP, Correction = "Best", ..., CheckArguments = TRUE)
{
  return(entropart::bcShannon(Ns=as.AbdVector(NorP), Correction=Correction, CheckArguments=CheckArguments))
}



#' Shannon Entropy of a vector of individuals
#'
#' Calculates the Shannon entropy of a vector of individuals.
#'
#'
#' @param NorP A vector of factors.
#' @inheritParams Shannon.wmppp
#'
#' @return A named number equal to the calculated entropy. The name is that of the bias correction used.
#'
#' @importFrom entropart Shannon
#' @method Shannon factor
#' @export
#'
#' @examples
#' Shannon(Paracou6$marks$PointType)
#' 
Shannon.factor <-
function(NorP, Correction = "Best", ..., CheckArguments = TRUE)
{
  return(entropart::bcShannon(Ns=as.numeric(table(NorP)), Correction=Correction, CheckArguments=CheckArguments))
}



#' Shannon Entropy of a vector of individuals
#'
#' Calculates the Shannon entropy of a vector of individuals.
#'
#'
#' @param NorP A vector of characters.
#' @inheritParams Shannon.wmppp
#'
#' @return A named number equal to the calculated entropy. The name is that of the bias correction used.
#'
#' @importFrom entropart Shannon
#' @method Shannon character
#' @export
#'
#' @examples
#' SpeciesList <- as.character(Paracou6$marks$PointType)
#' Shannon(SpeciesList)
#' 
Shannon.character <-
function(NorP, Correction = "Best", ..., CheckArguments = TRUE)
{
  return(entropart::bcShannon(Ns=as.numeric(table(NorP)), Correction=Correction, CheckArguments=CheckArguments))
}



#' Simpson Entropy of a spatialized community
#'
#' Calculates the Simpson entropy of a probability vector.
#'
#'
#' @param NorP An object of class [dbmss::wmppp], with `PointType` values as species names.
#' @param Correction A string containing one of the possible corrections: see [Simpson].
#' 
#' @inheritParams Tsallis.wmppp
#'
#' @return A named number equal to the calculated entropy. The name is that of the bias correction used.
#'
#' @importFrom entropart Simpson
#' @method Simpson wmppp
#' @export
#'
#' @examples
#' Simpson(Paracou6)
#' 
Simpson.wmppp <-
function(NorP, Correction = "Lande", ..., CheckArguments = TRUE)
{
  return(entropart::bcSimpson(Ns=as.AbdVector(NorP), Correction=Correction, CheckArguments=CheckArguments))
}



#' Simpson Entropy of a vector of individuals
#'
#' Calculates the Simpson entropy of a vector of individuals.
#'
#'
#' @param NorP A vector of factors.
#' @inheritParams Simpson.wmppp
#'
#' @return A named number equal to the calculated entropy. The name is that of the bias correction used.
#'
#' @importFrom entropart Simpson
#' @method Simpson factor
#' @export
#'
#' @examples
#' Simpson(Paracou6$marks$PointType)
#' 
Simpson.factor <-
function(NorP, Correction = "Lande", ..., CheckArguments = TRUE)
{
  return(entropart::bcSimpson(Ns=as.numeric(table(NorP)), Correction=Correction, CheckArguments=CheckArguments))
}



#' Simpson Entropy of a vector of individuals
#'
#' Calculates the Simpson entropy of a vector of individuals.
#'
#'
#' @param NorP A vector of characters.
#' @inheritParams Simpson.wmppp
#'
#' @return A named number equal to the calculated entropy. The name is that of the bias correction used.
#'
#' @importFrom entropart Simpson
#' @method Simpson character
#' @export
#'
#' @examples
#' SpeciesList <- as.character(Paracou6$marks$PointType)
#' Simpson(SpeciesList)
#' 
Simpson.character <-
function(NorP, Correction = "Lande", ..., CheckArguments = TRUE)
{
  return(entropart::bcSimpson(Ns=as.numeric(table(NorP)), Correction=Correction, CheckArguments=CheckArguments))
}
