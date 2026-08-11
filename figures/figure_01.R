# figures/figure_01.R
#
# Figure 1. Determinants of best-case transferability and their validation.
#
#   a  Empirical against theoretical transferability, across simulated pairs of
#      distinct populations. Reads the result of
#      simulation/03_theory_validation.R.
#   b  Response of transferability to the additive genetic variance.
#   c  Response of transferability to the average self-kinship.
#   d  Response of transferability to the non-genetic variance.
#
# Panels b, c and d carry no simulation. Each curve is the closed-form
# transferability ratio evaluated on a grid, so the panels show the exact theoretical
# response with no Monte Carlo error.
#
# Usage:  Rscript figures/figure_01.R

source("config.R")
source(file.path(PRS_ROOT, "R", "theory.R"))
source(file.path(PRS_ROOT, "R", "theme.R"))

suppressMessages({
  library(ggplot2)
  library(patchwork)
})

# --- panel a: validation against forward simulation -------------------------

validation <- readRDS(file.path(PRS_DATA, "theory_validation",
                                "distinct_populations.rds"))
scenarios <- validation$scenarios
r_squared <- validation$r_squared
axis_range <- range(c(scenarios$T_theoretical, scenarios$T_empirical))

panel_a <- ggplot(scenarios, aes(T_theoretical, T_empirical)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "red") +
  geom_point(colour = col_point, alpha = 0.55, size = 1.5, stroke = 0) +
  annotate("text", x = axis_range[1], y = axis_range[2], hjust = 0, vjust = 1,
           label = sprintf("r^2 == %.3f", r_squared), parse = TRUE, size = 4) +
  labs(x = expression(paste("Theoretical transferability  ", T[theo])),
       y = expression(paste("Empirical transferability  ", T[emp]))) +
  coord_fixed(xlim = axis_range, ylim = axis_range, clip = "off") +
  theme_classic(base_size = 12) +
  theme(axis.title = element_text(size = 11),
        plot.margin = margin(6, 8, 4, 6))

# --- panels b, c, d: the closed form, one determinant at a time --------------

# The symmetric baseline. With both populations at these values the two
# accuracies coincide and transferability is exactly one, so any departure seen
# in a panel is the effect of the single determinant being swept.
BASELINE_SIGMA2 <- 0.3
BASELINE_PHI <- 0.6
BASELINE_TAU2 <- 1

#' One panel of the parameter sweep
#'
#' Overlays two curves. The discovery curve varies the determinant in the
#' discovery population A while B stays at baseline; the target curve varies it
#' in the target population B while A stays at baseline. A change in A enters
#' the transferability ratio through its denominator and a change in B through
#' its numerator, so the two curves run in opposite directions and cross at the
#' baseline.
#'
#' @param grid Values of the determinant to sweep over.
#' @param determinant One of "sigma2", "phi" or "tau2".
#' @param x_label Axis label for the swept determinant.
#' @param label_A,label_B Legend labels for the two curves.
sweep_panel <- function(grid, determinant, x_label, label_A, label_B) {
  evaluate <- function(value, population) {
    args <- list(sigma2_A = BASELINE_SIGMA2, phi_bar_A = BASELINE_PHI,
                 tau2_A = BASELINE_TAU2,
                 sigma2_B = BASELINE_SIGMA2, phi_bar_B = BASELINE_PHI,
                 tau2_B = BASELINE_TAU2)
    prefix <- switch(determinant,
                     sigma2 = "sigma2_", phi = "phi_bar_", tau2 = "tau2_",
                     stop("determinant must be sigma2, phi or tau2"))
    args[[paste0(prefix, population)]] <- value
    do.call(transferability, args)
  }

  curves <- rbind(
    data.frame(x = grid, population = "A",
               value = vapply(grid, evaluate, numeric(1), population = "A")),
    data.frame(x = grid, population = "B",
               value = vapply(grid, evaluate, numeric(1), population = "B")))

  ggplot(curves, aes(x, value, colour = population, linetype = population)) +
    geom_hline(yintercept = 1, linetype = 3, colour = "grey55") +
    geom_line(linewidth = 1.1) +
    scale_colour_manual(values = c(A = col_discovery, B = col_target),
                        breaks = c("A", "B"),
                        labels = c(label_A, label_B), name = NULL) +
    scale_linetype_manual(values = c(A = "longdash", B = "solid"),
                          breaks = c("A", "B"),
                          labels = c(label_A, label_B), name = NULL) +
    labs(x = x_label,
         y = expression(paste("Transferability  ", T[A %->% B]))) +
    theme_classic(base_size = 12) +
    theme(axis.title = element_text(size = 11),
          legend.position = c(0.04, 0.96),
          legend.justification = c(0, 1),
          legend.background = element_rect(
            fill = scales::alpha("white", 0.75), colour = NA),
          legend.key = element_blank(),
          legend.key.width = unit(1.6, "lines"),
          legend.key.height = unit(1.0, "lines"),
          legend.text = element_text(size = 10),
          legend.spacing.y = unit(0, "pt"),
          plot.margin = margin(6, 8, 4, 6))
}

panel_b <- sweep_panel(
  seq(0.1, 1.0, length.out = 300), "sigma2",
  expression(paste("Additive genetic variance  ", sigma^2)),
  expression(sigma[A]^2), expression(sigma[B]^2))

panel_c <- sweep_panel(
  seq(0.5, 1.0, length.out = 300), "phi",
  expression(paste("Mean self-kinship  ", bar(Phi))),
  expression(bar(Phi)[A]), expression(bar(Phi)[B]))

panel_d <- sweep_panel(
  seq(0.3, 3.0, length.out = 300), "tau2",
  expression(paste("Non-genetic variance  ", tau^2)),
  expression(tau[A]^2), expression(tau[B]^2))

figure <- (panel_a | panel_b) / (panel_c | panel_d) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 14))

ggsave(file.path(PRS_OUTPUT, "Fig1.pdf"), figure, width = 9.2, height = 7.4)
ggsave(file.path(PRS_OUTPUT, "Fig1.png"), figure, width = 9.2, height = 7.4,
       dpi = 300)
cat("Figure 1 written to", PRS_OUTPUT, "\n")
