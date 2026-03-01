library(ggplot2)
library(grid)

# ---- Main vertical boxes ----
boxes <- data.frame(
  xmin = c(0, 0, 0, 0, 0),
  xmax = c(6, 6, 6, 6, 6),
  ymin = c(10, 7, 4, 1, -2),
  ymax = c(12, 9, 6, 3, 0),
  label = c(
    "Downloaded NHANES Participants\n(N = 80,312)",
    "Adults Age ≥ 20\n(N = 44,790)",
    "Non-missing PhenoAge, CKD,\nSurvey Weights\n(N = 23,387)",
    "Primary PFAS Complete-Case\n(N = 4,696)",
    "Final Analytic Sample\nModels A & B\n(N = 1,994)"
  )
)

# ---- Exclusion boxes ----
exclusions <- data.frame(
  xmin = c(9, 9, 9, 9),
  xmax = c(15, 15, 15, 15),
  ymin = c(10, 7, 4, 1),
  ymax = c(12, 9, 6, 3),
  label = c(
    "Excluded: Age < 20\n(N = 35,522)",
    "Excluded: Missing PhenoAge or CKD\n(N = 21,403)",
    "Excluded: Missing Primary PFAS\n(N = 18,691)",
    "Excluded: Missing Covariates\n(N = 2,702)"
  )
)

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
  geom_segment(aes(x = 3, xend = 3, y = 10, yend = 9),
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(aes(x = 3, xend = 3, y = 7, yend = 6),
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(aes(x = 3, xend = 3, y = 4, yend = 3),
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(aes(x = 3, xend = 3, y = 1, yend = 0),
               arrow = arrow(length = unit(0.25, "cm"))) +
  
  # Horizontal exclusion arrows
  geom_segment(aes(x = 6, xend = 9, y = 11, yend = 11),
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(aes(x = 6, xend = 9, y = 8, yend = 8),
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(aes(x = 6, xend = 9, y = 5, yend = 5),
               arrow = arrow(length = unit(0.25, "cm"))) +
  geom_segment(aes(x = 6, xend = 9, y = 2, yend = 2),
               arrow = arrow(length = unit(0.25, "cm"))) +
  
  # Text labels
  geom_text(data = boxes,
            aes(x = 3, y = (ymin + ymax)/2, label = label),
            size = 3.8) +
  
  geom_text(data = exclusions,
            aes(x = 12, y = (ymin + ymax)/2, label = label),
            size = 3.6) +
  
  # Title
  labs(title = "Participant Flow Chart") +
  
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
  ylim(-3, 13)

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