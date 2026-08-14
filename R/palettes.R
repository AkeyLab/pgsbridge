#' Colour palettes shared by every figure
#'
#' One palette per role, so that a colour means the same thing wherever it
#' appears.  All are drawn from Okabe and Ito's eight-colour set, which stays
#' distinguishable under the three common forms of colour vision deficiency.
#'
#' `okabe_ito` is the full set, named by colour.  The others select from it:
#'
#' \describe{
#'   \item{`pair_pal`}{the three unordered population pairs, used wherever
#'     transferability is shown in both directions.}
#'   \item{`ratio_pal`}{the rare-to-common effect-size ratio, ordered from the
#'     smallest ratio in the coolest colour to the largest in the warmest.}
#'   \item{`decomp_pal`}{the two variant classes in the heritability
#'     decomposition.}
#'   \item{`col_discovery`, `col_target`}{the discovery and target populations,
#'     used when a determinant is swept in one or the other.}
#'   \item{`col_point`}{scatter points, where no grouping variable applies.}
#' }
#'
#' @format Named character vectors of hexadecimal colours; `col_discovery` and
#'   `col_target` are single strings.
#'
#' @examples
#' okabe_ito
#' pair_pal
#'
#' @name palettes
NULL

#' @rdname palettes
#' @export
okabe_ito <- c(black      = "#000000",
               orange     = "#E69F00",
               skyblue    = "#56B4E9",
               green      = "#009E73",
               yellow     = "#F0E442",
               blue       = "#0072B2",
               vermillion = "#D55E00",
               purple     = "#CC79A7",
               grey       = "#999999")

#' @rdname palettes
#' @export
col_discovery <- unname(okabe_ito["vermillion"])

#' @rdname palettes
#' @export
col_target <- unname(okabe_ito["blue"])

#' @rdname palettes
#' @export
col_point <- unname(okabe_ito["blue"])

#' @rdname palettes
#' @export
pair_pal <- c("ASN-AFR" = unname(okabe_ito["vermillion"]),
              "EUR-AFR" = unname(okabe_ito["blue"]),
              "ASN-EUR" = unname(okabe_ito["green"]))

#' @rdname palettes
#' @export
ratio_pal <- c("1"  = unname(okabe_ito["skyblue"]),
               "2"  = unname(okabe_ito["blue"]),
               "5"  = unname(okabe_ito["orange"]),
               "10" = unname(okabe_ito["vermillion"]))

#' @rdname palettes
#' @export
decomp_pal <- c("common" = unname(okabe_ito["blue"]),
                "rare"   = unname(okabe_ito["orange"]))

#' Base font size shared by every figure
#'
#' The point size passed to the `ggplot2` theme in each figure, so that text is
#' the same size across the whole set.
#'
#' @format A single number.
#'
#' @examples
#' PRS_BASE_SIZE
#'
#' @export
PRS_BASE_SIZE <- 11
