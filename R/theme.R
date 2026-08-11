# R/theme.R
#
# One palette and one font size for every figure in the paper, so that a colour
# means the same thing wherever it appears.
#
# The palette is Okabe and Ito's eight-colour set, which stays distinguishable
# under the three common forms of colour vision deficiency.

suppressMessages(library(ggplot2))

okabe_ito <- c(black      = "#000000",
               orange     = "#E69F00",
               skyblue    = "#56B4E9",
               green      = "#009E73",
               yellow     = "#F0E442",
               blue       = "#0072B2",
               vermillion = "#D55E00",
               purple     = "#CC79A7",
               grey       = "#999999")

# Base font size, in points, shared by every figure.
PRS_BASE_SIZE <- 11

# --- fixed colour roles ----------------------------------------------------

# Figure 1, panels b to d: which population the swept determinant belongs to.
col_discovery <- unname(okabe_ito["vermillion"])   # discovery population A
col_target    <- unname(okabe_ito["blue"])         # target population B

# Figure 1, panel a: the validation scatter.
col_point <- unname(okabe_ito["blue"])

# Figures 2, S2 and S3: which pair of populations a line belongs to.
pair_pal <- c("ASN-AFR" = unname(okabe_ito["vermillion"]),
              "EUR-AFR" = unname(okabe_ito["blue"]),
              "ASN-EUR" = unname(okabe_ito["green"]))

# Figures 3 and S4 to S6: the rare-to-common effect-size ratio, ordered from the
# smallest ratio in the coolest colour to the largest in the warmest.
ratio_pal <- c("1"  = unname(okabe_ito["skyblue"]),
               "2"  = unname(okabe_ito["blue"]),
               "5"  = unname(okabe_ito["orange"]),
               "10" = unname(okabe_ito["vermillion"]))

# Figure 4: the two variant classes in the heritability decomposition.
decomp_pal <- c("common" = unname(okabe_ito["blue"]),
                "rare"   = unname(okabe_ito["orange"]))
