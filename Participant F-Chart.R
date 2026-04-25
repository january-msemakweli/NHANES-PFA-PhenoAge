library(ggplot2)
library(grid)

# ---- Main vertical boxes (simplified to key analytic steps) ----
main_labels <- c(
  "Downloaded NHANES Participants\n(Selected Cycles)\n(N = 80,312)",
  "Age >= 20\n(N = 44,790)",
  "Non-missing PhenoAge, CKD,\nand Required Modeling Components\n(N = 23,387)",
  "Non-missing Selected Weight,\nPSU, and Strata (Analytic Sample)\n(N = 15,239)",
  "Complete-Case Modeling Sample\n(qgcomp-Compatible PFAS Subset)\n(N = 3,058)"
)

exclusion_labels <- c(
  "Excluded: Age < 20\n(N = 35,522)",
  "Excluded: Missing PhenoAge, CKD,\nor Required Modeling Components\n(N = 21,403)",
  "Excluded: Missing Weight/PSU/Strata\n(N = 8,148)",
  "Excluded:\nNot in Complete-Case\nModeling Sample\n(N = 12,181)"
)

top_y <- 20
step_gap <- 3
box_height <- 2

boxes <- data.frame(
  xmin = 0,
  xmax = 6,
  ymax = top_y - (seq_along(main_labels) - 1) * step_gap,
  label = main_labels
)
boxes$ymin <- boxes$ymax - box_height

# ---- Exclusion boxes ----
exclusions <- data.frame(
  xmin = 9,
  xmax = 15,
  ymax = boxes$ymax[seq_along(exclusion_labels)],
  label = exclusion_labels
)
exclusions$ymin <- exclusions$ymax - box_height

p <- ggplot() +
  
  # Main boxes
  geom_rect(data = boxes,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax),
            fill = "white", color = "black", linewidth = 1) +
  
  # Exclusion boxes
  geom_rect(data = exclusions,
            aes(xmin = xmin, xmax = xmax,
                ymin = ymin, ymax = ymax),
            fill = "white", color = "black", linewidth = 1) +
  
  # Vertical arrows
  geom_segment(
    data = data.frame(
      x = 3, xend = 3,
      y = boxes$ymin[-nrow(boxes)],
      yend = boxes$ymax[-1]
    ),
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  
  # Horizontal exclusion arrows
  geom_segment(
    data = data.frame(
      x = 6, xend = 9,
      y = (boxes$ymin[seq_len(nrow(exclusions))] + boxes$ymax[seq_len(nrow(exclusions))]) / 2,
      yend = (boxes$ymin[seq_len(nrow(exclusions))] + boxes$ymax[seq_len(nrow(exclusions))]) / 2
    ),
    aes(x = x, xend = xend, y = y, yend = yend),
    arrow = arrow(length = unit(0.25, "cm"))
  ) +
  
  # Text labels
  geom_text(data = boxes,
            aes(x = 3, y = (ymin + ymax)/2, label = label),
            size = 3.2) +
  
  geom_text(data = exclusions,
            aes(x = 12, y = (ymin + ymax)/2, label = label),
            size = 3.0, lineheight = 0.95) +
  
  theme_void() +
  theme(
    text = element_text(family = "Times New Roman"),
    plot.title = element_text(
      hjust = 0.5,
      size = 16,
      face = "bold"
    )
  ) +
  
  coord_fixed() +
  xlim(-1, 16) +
  ylim(min(boxes$ymin) - 1, top_y + 1)

# Save as publishable-quality JPEG (300 dpi, Times New Roman)
ggsave(
  "Participant_F-Chart.jpg",
  plot = p,
  width = 10,
  height = 8,
  dpi = 300,
  device = "jpeg",
  quality = 95,
  bg = "white"
)
print(p)