#' An extract of the Montreal road network
#'
#' A SpatialLinesDataFrame representing a portion of the road network of
#' Montreal.
#'
#' @format A SpatialLinesDataFrame with 2946 rows and 4 variables:
#' \describe{
#'   \item{IdRte}{a unique id for each road}
#'   \item{ClsRte}{the type of the road}
#'   \item{ClsFnction}{the function of the road}
#'   \item{nbAccident}{the number of
#'   bike accidents recorded on the road} ... }
#' @source \url{http://donnees.ville.montreal.qc.ca/}
"mtl_network"

#' A smaller extract of the Montreal road network
#'
#' A SpatialLinesDataFrame representing a smaller portion of the road
#' network of Montreal.
#'
#' @format A SpatialLinesDataFrame with 1244 rows and 4 variables:
#' \describe{
#'   \item{IdRte}{a unique id for each road}
#'   \item{ClsRte}{the type of the road}
#'   \item{ClsFnction}{the function of the road}
#'   \item{nbAccident}{the number of
#'   bike accidents recorded on the road} ... }
#' @source \url{http://donnees.ville.montreal.qc.ca/}
"small_mtl_network"



#' An extract of the bike accidents on Montreal road network
#'
#' A SpatialPointsDataFrame representing bike accidents that occured on the
#' montreal road network in 2016
#'
#' @format A SpatialLinesDataFrame with 2946 rows and 4 variables: \describe{
#'   \item{LOC_X}{The X coordinates of the accident (EPSG : 32188)}
#'   \item{LOC_Y}{The Y coordinates of the accident (EPSG : 32188)}
#'   \item{NB_VICTIME}{the number of victims} \item{AN}{the year of the
#'   accident} ... }
#' @source \url{http://donnees.ville.montreal.qc.ca/}
"bike_accidents"
