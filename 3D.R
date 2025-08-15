# Libraries

library(fda)
library(dplyr)
library(purrr)
library(plotly)

all_smooth_x <- readRDS("smoothed_x_all_subjects.rds")
all_smooth_y <- readRDS("smoothed_y_all_subjects.rds")
all_smooth_z <- readRDS("smoothed_z_all_subjects.rds")


#  time grid from an fd object's basis range
make_time_grid <- function(fd_obj, n = 200) {
  rg <- fd_obj$basis$rangeval
  seq(rg[1], rg[2], length.out = n)
}

# Create a common time grid using X's first subject (assumes same range across axes)
unitime <- make_time_grid(all_smooth_x[[1]], n = 200)

# convert a list of fd objects (subjects) into a long data.frame
fd_list_to_long <- function(fd_list, unitime, axis_name = "axis") {
  
  subjects <- sprintf("S%03d", seq_along(fd_list))
  
  dfs <- lapply(seq_along(fd_list), function(i) {
    fd_i <- fd_list[[i]]
    mat  <- eval.fd(unitime, fd_i)  # returns length(unitime) x nrep
    if (is.null(dim(mat))) mat <- matrix(mat, ncol = 1)
    data.frame(
      time      = rep(unitime, ncol(mat)),
      amplitude = as.vector(mat),
      curve     = rep(seq_len(ncol(mat)), each = length(unitime)),
      Player_ID = subjects[i],
      axis      = axis_name
    )
  })
  bind_rows(dfs)
}

# Build long data for X, Y, Z

df_x <- fd_list_to_long(all_smooth_x, unitime, axis_name = "X")
df_y <- fd_list_to_long(all_smooth_y, unitime, axis_name = "Y")
df_z <- fd_list_to_long(all_smooth_z, unitime, axis_name = "Z")

# Subject-level mean at each time for each axis

subject_means_x <- df_x %>%
  group_by(Player_ID, time) %>%
  summarise(mean_x = mean(amplitude, na.rm = TRUE), .groups = "drop")

subject_means_y <- df_y %>%
  group_by(Player_ID, time) %>%
  summarise(mean_y = mean(amplitude, na.rm = TRUE), .groups = "drop")

subject_means_z <- df_z %>%
  group_by(Player_ID, time) %>%
  summarise(mean_z = mean(amplitude, na.rm = TRUE), .groups = "drop")

# Merge XYZ by subject and time
subject_means_xyz <- subject_means_x %>%
  left_join(subject_means_y, by = c("Player_ID", "time")) %>%
  left_join(subject_means_z, by = c("Player_ID", "time")) %>%
  arrange(Player_ID, time)


# Interactive 3D plot: one line per subject (mean X,Y,Z over time)

p <- plot_ly(
  subject_means_xyz,
  x = ~mean_x, y = ~mean_y, z = ~mean_z,
  split = ~Player_ID,  # groups lines by subject
  type = "scatter3d", mode = "lines"
) %>%
  layout(
    title = "Mean Swing Trajectories in 3D (All Subjects)",
    scene = list(
      xaxis = list(title = "X"),
      yaxis = list(title = "Y"),
      zaxis = list(title = "Z")
    ),
    showlegend = FALSE
  )

p
# Pick a nice camera angle
camera_view <- list(
  eye = list(x = 1.5, y = 1.5, z = 1)  # change values until you like the view
)

p <- p %>% layout(scene = list(camera = camera_view))

# Save as HTML 
htmlwidgets::saveWidget(p, "mean_swing_3d.html", selfcontained = TRUE)


