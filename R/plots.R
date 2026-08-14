#' Transferability in both directions, for three population pairs
#'
#' The shape used by Figures 2, S2 and S3. Each unordered population pair is one
#' coloured line joining its two directions, so the slope of the line shows how
#' much transferability depends on which population is the discovery one and
#' which is the target. Panels are the three heritability levels.
#'
#' Each point is labelled with the ordered pair it represents and its value, the
#' left column labelled to its left and the right column to its right. Reading
#' the number off the axis is not accurate enough to compare pairs that differ
#' in the third decimal place, which several of them do.
#'
#' Where three lines lie close together the labels would print on top of one
#' another, so within each panel and each side they are pushed apart vertically
#' to a minimum separation before being drawn. The leader is the point, not the
#' label, so a nudged label still reads against its own colour.
#'
#' @param summary_df Output of `read_transferability_by_heritability()`.
#' @param label_points Whether to draw the value labels.
#' @param label_gap Horizontal offset of the labels, in units of the axis, which
#'   runs from 1 to 2. Widen it if the labels collide with the error bars.
#' @param label_margin Horizontal room reserved for the labels on each side, in
#'   the same units. Widen it if the labels are clipped by the panel edge.
#' @param min_label_separation Smallest vertical gap between two labels in the
#'   same panel and on the same side, as a fraction of the vertical range.
#' @return A ggplot object.
#' @export
plot_direction_pairs <- function(summary_df, label_points = TRUE,
                                 label_gap = 0.06, label_margin = 1.5,
                                 min_label_separation = 0.055) {
  # The horizontal axis has the discovery-to-target direction on the left. The
  # ordered pair a point represents is the unordered label for the left column
  # and its reverse for the right one.
  summary_df$on_left <- summary_df$direction == "Discovery-Target"
  summary_df$x <- ifelse(summary_df$on_left, 1, 2)

  reverse_pair <- function(label) {
    parts <- strsplit(as.character(label), "-", fixed = TRUE)
    vapply(parts, function(p) paste(rev(p), collapse = "-"), character(1))
  }
  summary_df$ordered_pair <- ifelse(summary_df$on_left,
                                    as.character(summary_df$pop_pair),
                                    reverse_pair(summary_df$pop_pair))
  summary_df$label <- sprintf(
    "%s: %.3f",
    sub("-", " to ", summary_df$ordered_pair, fixed = TRUE),
    summary_df$emp_tau)

  # Push labels apart where the lines they belong to lie close together. Working
  # upwards from the lowest label in a group, each is raised to at least
  # `separation` above the one below it.
  separation <- min_label_separation *
    diff(range(c(summary_df$emp_tau - summary_df$ci,
                 summary_df$emp_tau + summary_df$ci)))
  summary_df$label_y <- summary_df$emp_tau
  for (panel in unique(summary_df$h2)) {
    for (side in c(TRUE, FALSE)) {
      rows <- which(summary_df$h2 == panel & summary_df$on_left == side)
      rows <- rows[order(summary_df$emp_tau[rows])]
      if (length(rows) < 2) next
      for (k in seq_along(rows)[-1]) {
        lowest_allowed <- summary_df$label_y[rows[k - 1]] + separation
        if (summary_df$label_y[rows[k]] < lowest_allowed) {
          summary_df$label_y[rows[k]] <- lowest_allowed
        }
      }
    }
  }

  # A constant column giving the spanning label above the panel strips.
  summary_df$heritability_title <- factor("Heritability")

  plot <- ggplot(summary_df,
                 aes(x = x, y = emp_tau, colour = pop_pair, group = pop_pair)) +
    geom_errorbar(aes(ymin = emp_tau - ci, ymax = emp_tau + ci), width = 0.05) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.5) +
    facet_nested(. ~ heritability_title + h2) +
    scale_x_continuous(limits = c(1 - label_margin, 2 + label_margin),
                       breaks = c(1, 2)) +
    ylab("Transferability") +
    scale_color_manual("Pop Pair", values = pair_pal) +
    theme_bw(base_size = PRS_BASE_SIZE) +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.background = element_rect(fill = "grey92", colour = NA),
          legend.position = "bottom")

  if (label_points) {
    plot <- plot +
      geom_text(aes(x = x + ifelse(on_left, -label_gap, label_gap),
                    y = label_y, label = label,
                    hjust = ifelse(on_left, 1, 0)),
                colour = "black", size = 2.6, show.legend = FALSE)
  }

  plot
}

#' Transferability against the fraction of rare causal variants
#'
#' The shape used by Figures 3 and S4 to S6. One line per rare-to-common
#' effect-size ratio, with the fraction of causal variants that are rare on the
#' horizontal axis.
#'
#' @param summary_df Output of `read_rare_variant_sweep()`.
#' @param title Panel title, naming the direction of transfer.
#' @param y_max Upper limit of the vertical axis. The panels of one figure share
#'   a limit so that they can be read against each other.
#' @return A ggplot object.
#' @export
plot_rare_variant_sweep <- function(summary_df, title, y_max) {
  dodge <- position_dodge(0.01)
  ggplot(summary_df, aes(x = fraction, y = emp_tau, colour = ratio)) +
    geom_errorbar(aes(ymin = emp_tau - ci, ymax = emp_tau + ci),
                  width = 0.1, position = dodge) +
    geom_line(position = dodge, linewidth = 0.8) +
    geom_point(position = dodge) +
    labs(colour = "Rare/Common Effect Size Ratio") +
    scale_color_manual(values = ratio_pal) +
    xlab("Fraction of Rare Variants") +
    ylab("Transferability") +
    ggtitle(title) +
    ylim(0, y_max) +
    scale_x_continuous(breaks = seq(0, 1, by = 0.1)) +
    theme_minimal(base_size = PRS_BASE_SIZE) +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = PRS_BASE_SIZE))
}

#' Both directions of one population pair, side by side
#'
#' Figures 3 and S4 to S6 all show a pair of populations as two panels, the
#' first transferring from the first population to the second and the second
#' transferring back. This reads both directions and arranges them under a
#' single shared legend.
#'
#' @param data_dir Directory holding the rare-variant sweep files.
#' @param h2_tag File-name tag for the European heritability.
#' @param first,second The two population labels, for example "AFR" and "EUR".
#'   The left panel is `first` to `second`.
#' @param y_max Upper limit of the vertical axis, shared by both panels.
#' @return A patchwork arrangement of the two panels under one legend.
#' @export
arrange_rare_variant_pair <- function(data_dir, h2_tag, first, second, y_max) {
  forward <- plot_rare_variant_sweep(
    read_rare_variant_sweep(data_dir, h2_tag, paste0(first, "-", second)),
    title = paste(first, "to", second), y_max = y_max)
  backward <- plot_rare_variant_sweep(
    read_rare_variant_sweep(data_dir, h2_tag, paste0(second, "-", first)),
    title = paste(second, "to", first), y_max = y_max)
  (forward | backward) +
    plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
}
