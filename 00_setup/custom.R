#########################
# Set up my environment #
#########################

# packages ----------------------------------------------------------------

# packages I need
library(tidyverse)
library(coolorrr)

# some customizations
set_palette() # set color palette with defaults
theme_set(ggthemes::theme_fivethirtyeight())
theme_update(
  axis.title = element_text(hjust = 0.9),
  plot.title.position = "plot",
  panel.background = element_rect(fill = "white"),
  plot.background = element_rect(fill = "white")
)
