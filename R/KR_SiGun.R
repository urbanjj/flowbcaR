#' Administrative Boundaries of Si-Gun in South Korea
#'
#' An \code{sf} object containing the administrative boundaries for Si (cities)
#' and Gun (counties) in the Republic of Korea. This dataset is intended for
#' use with the spatial analysis and visualization functions within the package.
#'
#' @name KR_SiGun
#' @docType data
#' @format A simple feature collection (\code{sf}) with 159 features and 2 fields.
#' The coordinate reference system is KGD2002 / Unified CS (EPSG: 5179).
#' \describe{
#' \item{SiGun_CD}{Unique administrative code for the Si/Gun.}
#' \item{SiGun_NM}{Romanized name of the Si/Gun.}
#' \item{geom}{The \code{sfc} geometry column (POLYGON/MULTIPOLYGON).}
#' }
#'
#' @source This data is derived from public administrative boundary data of 2023 provided by
#' the Korean government. The data has been processed for use within this package.
#'
#' @keywords datasets spatial korea administrative
#'
#' @examples
#' # To load the data from the package:
#' data(KR_SiGun)
#'
"KR_SiGun"
