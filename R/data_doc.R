#' An extract of the Montreal road network
#'
#' A SpatialLinesDataFrame representing a portion of the Montreal road network
#'
#' @format A SpatialLinesDataFrame with 2946 rows and 4 variables:
#' \describe{
#'   \item{IdRte}{a unique id for each road}
#'   \item{ClsRte}{the category of the road}
#'   \item{ClsFnction}{the function of the road}
#'   \item{nbAccident}{the number of
#'   bike accidents recorded on the road} ... }
#' @source \url{http://donnees.ville.montreal.qc.ca/}
"mtl_network"

#' A smaller extract of the Montreal road network
#'
#' A SpatialLinesDataFrame representing a smaller portion of the Montreal road
#'  network.
#'
#' @format A SpatialLinesDataFrame with 1244 rows and 4 variables:
#' \describe{
#'   \item{IdRte}{a unique id for each road}
#'   \item{ClsRte}{the category of the road}
#'   \item{ClsFnction}{the function of the road}
#'   \item{nbAccident}{the number of
#'   bike accidents recorded on the road} ... }
#' @source \url{http://donnees.ville.montreal.qc.ca/}
"small_mtl_network"



#' An extract of the bike accidents on the Montreal road network
#'
#' A SpatialPointsDataFrame representing bike accidents that occurred on the
#' Montreal road network in 2016
#'
#' @format A SpatialPointsDataFrame with 2946 rows and 4 variables: \describe{
#'   \item{LOC_X}{The X coordinates of the accident (EPSG : 32188)}
#'   \item{LOC_Y}{The Y coordinates of the accident (EPSG : 32188)}
#'   \item{NB_VICTIME}{the number of victims} \item{AN}{the year of the
#'   accident} ... }
#' @source \url{http://donnees.ville.montreal.qc.ca/}
"bike_accidents"


#' An extract of the Montreal road network with only the main roads
#'
#' A SpatialLinesDataFrame representing the main roads of the Montreal road
#' network.
#'
#' @format A SpatialLinesDataFrame with 17304 rows and 3 variables: \describe{
#'   \item{CLASS}{The category of the road} \item{NAME}{The name of the road}
#'   \item{TYPE}{the type of the road} }
#' @source \url{http://donnees.ville.montreal.qc.ca/}
"main_network_mtl"


#' Libraries
#'
#' A SpatialPointsDataFrame representing libraries in Montreal
#'
#' @format A SpatialPointsDataFrame with 55 rows and 2 variables
#' @source \url{http://donnees.ville.montreal.qc.ca/}
"libraries_mtl"

#' Theatres
#'
#' A SpatialPointsDataFrame representing theatres in Montreal
#'
#' @format A SpatialPointsDataFrame with 54 rows and 2 variables
#' @source \url{http://donnees.ville.montreal.qc.ca/}
"theatres_mtl"

