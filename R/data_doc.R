#' Road accidents including a bicyle in Montreal in 2016
#'
#' A feature collection (sf object) representing road accidents including a
#' cyclist in Montreal in 2016. The EPSG is 3797, and the data comes from the
#' Montreal OpenData website. It is only a small subset in central districts used to demonstrate the
#' main functions of spNetwork.
#'
#' @format A sf object with 347 rows and 4 variables
#' \describe{
#'   \item{NB_VICTIME}{the number of victims}
#'   \item{AN}{the year of the accident}
#'   \item{Date}{the date of the accident (yyyy/mm/dd)}
#'   \item{geom}{the geometry (points)}
#' }
#' @source <https://donnees.montreal.ca/ville-de-montreal/collisions-routieres>
"bike_accidents"

#' Road network of Montreal
#'
#' A feature collection (sf object) representing the road network of Montreal. The EPSG is 3797, and the data comes from the
#' Montreal OpenData website. It is only a small subset in central districts used to demonstrate the
#' main functions of spNetwork.
#'
#' @format A sf object with 2945 rows and 2 variables
#' \describe{
#'   \item{ClsRte}{the category of the road}
#'   \item{geom}{the geometry (linestrings)}
#' }
#' @source <https://donnees.montreal.ca/ville-de-montreal/geobase>
"mtl_network"

#' Primary road network of Montreal
#'
#' A feature collection (sf object) representing the primary road network of Montreal. The EPSG is 3797, and the data comes from the
#' Montreal OpenData website.
#'
#' @format A sf object with 2945 rows and 2 variables
#' \describe{
#'   \item{TYPE}{the type of road}
#'   \item{geom}{the geometry (linestrings)}
#' }
#' @source <https://donnees.montreal.ca/ville-de-montreal/geobase>
"main_network_mtl"

#' Smaller subset road network of Montreal
#'
#' A feature collection (sf object) representing the road network of Montreal. The EPSG is 3797, and the data comes from the
#' Montreal OpenData website. It is only a small extract in central districts used to demonstrate the
#' main functions of spNetwork. It is mainly used internally for tests.
#'
#' @format A sf object with 1244 rows and 2 variables
#' \describe{
#'   \item{TYPE}{the type of road}
#'   \item{geom}{the geometry (linestrings)}
#' }
#' @source <https://donnees.montreal.ca/ville-de-montreal/geobase>
"small_mtl_network"


#' Libraries of Montreal
#'
#' A feature collection (sf object) representing the libraries of Montreal. The EPSG is 3797 and the data comes from the
#' Montreal OpenData website.
#'
#' @format A sf object with 55 rows and 3 variables.
#' \describe{
#'   \item{CP}{the postal code}
#'   \item{NAME}{the name of the library}
#'   \item{geom}{the geometry (points)}
#' }
#' @source <https://donnees.montreal.ca/ville-de-montreal/lieux-culturels>
"mtl_libraries"

#' Theatres of Montreal
#'
#' A feature collection (sf object) representing the theatres of Montreal. The EPSG is 3797 and the data comes from the
#' Montreal OpenData website.
#'
#' @format A sf object with 54 rows and 3 variables.
#' \describe{
#'   \item{CP}{the postal code}
#'   \item{NAME}{the name of the theatre}
#'   \item{geom}{the geometry (points)}
#' }
#' @source <https://donnees.montreal.ca/ville-de-montreal/lieux-culturels>
"mtl_theatres"

