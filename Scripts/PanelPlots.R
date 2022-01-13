###############################
### Panel plot for figure 6 ###
###############################

panel6 = ggdraw() +
  draw_plot(p, x = 0, y = .5, width = .5, height = .5) +
  draw_plot(pRA, x = .5, y = .5, width = .5, height = .5) +
  draw_plot(Fig6C, x = 0, y = 0, width = 1, height = 0.5) +
  draw_plot_label(label = c("A", "B", "C"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))
panel6 

ggsave("Fig6.svg",
       plot = panel6,
       units="cm",
       width=18,
       height=18,
       dpi = 300)
