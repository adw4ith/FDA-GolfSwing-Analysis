library(fda)
library(MFPCA)
library(R.matlab)
library(here)

subjects <- c(4,5,7,9,10:27,29:49,52:57,59:63)
subjects <- paste0("00", subjects)
subjects[1:4] <- paste0("0", subjects[1:4])
subjects <- paste0("AMID", subjects)

# Data and Basis Setup 
npts <- 1498
unitime <- seq(0, 1, length.out = npts)
nbasis <- 30
norder <- 6
basis_fd <- create.bspline.basis(c(0, 1), nbasis = nbasis, norder = norder)
fdpar_opt <- fdPar(basis_fd, 2, 1e-6)

# Process All Subjects 
all_x <- list(); all_y <- list(); all_z <- list()
n_shots_per_subject <- integer(length(subjects))  # Number of curves per subject

for (s in seq_along(subjects)) {
  file_path <- here("Data", "Mat files", paste0(subjects[s], ".mat"))
  raw <- readMat(file_path)
  shot_nums <- as.vector(raw$Data[[2]])
  n_shots <- length(shot_nums)
  n_shots_per_subject[s] <- n_shots
}

#  Build Per-Curve Subject ID Vector
subject_ids <- rep(subjects, times = n_shots_per_subject)
drop_subject <- 3; drop_shot <- 11
drop_idx <- sum(n_shots_per_subject[1:(drop_subject-1)]) + drop_shot
subject_ids <- subject_ids[-drop_idx]

# Load Registered Curves
reg_fd <- readRDS("registered_xyz_fd.rds")  # X, Y, Z registered fd objects

# MFPCA 
time_points <- seq(0, 1, length.out = 100)
x_data <- eval.fd(time_points, reg_fd$x)
y_data <- eval.fd(time_points, reg_fd$y)
z_data <- eval.fd(time_points, reg_fd$z)

funData_list <- list(
  funData(argvals = time_points, X = t(x_data)),
  funData(argvals = time_points, X = t(y_data)),
  funData(argvals = time_points, X = t(z_data))
)
mfd <- multiFunData(funData_list)

mfpca_result <- MFPCA(mfd, M = 6, uniExpansions = list(
  list(type = "splines1D", k = 10),
  list(type = "splines1D", k = 10),
  list(type = "splines1D", k = 10)
))

# Extract eigenvalues from the MFPCA result
eigenvalues <- mfpca_result$values

# Calculate variance explained by each FPC
variance_explained <- eigenvalues / sum(eigenvalues)

# Print as percentages, rounded to two decimals
for (i in seq_along(variance_explained)) {
  cat(sprintf("FPC%d explains %.2f%% of the total variance\n",
              i, 100 * variance_explained[i]))
}

# Assume mfpca_result is your MFPCA output
n_pc <- length(mfpca_result$values)
axis_names <- c("X", "Y", "Z")

par(mfrow = c(1, 3))
for (axis in 1:3) {
  matplot(
    t(mfpca_result$functions[[axis]]@X), type = "l", lty = 1, 
    main = paste(axis_names[axis], "axis eigenfunctions"), 
    xlab = "Basis / Time Index", ylab = "Eigenfunction Value"
  )
  legend("topright", legend = paste0("PC", 1:n_pc), col = 1:n_pc, lty = 1)
}
par(mfrow = c(1,1))

# Variance Explained
var_explained <- mfpca_result$values / sum(mfpca_result$values)
cum_var <- cumsum(var_explained)
n_pc <- length(var_explained)
cum_var_percent <- 100 * cum_var
par(mfrow = c(1,2))
barplot(100 * var_explained, names.arg = paste0("PC", 1:n_pc),
        ylab = "Variance Explained (%)", main = "Scree Plot")
plot(1:n_pc, 100 * cum_var, type = "b", pch=16,
     xlab = "Number of PCs", ylab = "Cumulative Variance (%)",
     main = " Variance Explained")
abline(h=90, lty=2, col="red")
text(x = 1:n_pc, y = cum_var_percent, 
     labels = paste0(sprintf("%.1f", cum_var_percent), "%"), 
     pos = 3, cex = 0.8)
par(mfrow = c(1,1))

# Axis Contributions Table 
calculate_axis_variance <- function(mfpca_res, mfd_data) {
  n_components <- length(mfpca_res$values)
  n_axes <- length(mfd_data)
  axis_names <- c("X", "Y", "Z")
  eigenfunctions <- mfpca_res$functions
  axis_variance <- array(0, dim = c(n_components, n_axes))
  colnames(axis_variance) <- axis_names
  rownames(axis_variance) <- paste0("PC", 1:n_components)
  marginal_vars <- numeric(n_axes)
  for(j in 1:n_axes) {
    data_matrix <- t(mfd_data[[j]]@X)
    marginal_vars[j] <- var(as.vector(data_matrix))
  }
  for(i in 1:n_components) {
    axis_norms <- numeric(n_axes)
    for(j in 1:n_axes) {
      eigenfunction_j <- eigenfunctions[[j]]@X[i, ]
      axis_norms[j] <- sum(eigenfunction_j^2) * marginal_vars[j]
    }
    total_norm <- sum(axis_norms)
    if (total_norm > 0) axis_variance[i, ] <- axis_norms / total_norm
  }
  return(axis_variance)
}
library(ggpubr)
axis_contributions <- calculate_axis_variance(mfpca_result, mfd)
print(round(axis_contributions, 3))
axis_table <- round(axis_contributions, 3)
axis_table_df <- as.data.frame(axis_table)
axis_table_df <- cbind(PC = rownames(axis_table), axis_table_df)

p <- ggtexttable(axis_table_df, rows = NULL, theme = ttheme("classic"))
print(p)


# Plot Eigenfunctions/Harmonics 
n_pcs_to_plot <- min(4, length(mfpca_result$values))
axis_labels <- c("X", "Y", "Z")
axis_colors <- c("steelblue", "tomato", "seagreen")
par(mfrow = c(n_pcs_to_plot, 3), mar = c(3, 3, 2, 1))
for (i in 1:n_pcs_to_plot) {
  for (j in 1:3) {
    time <- mfpca_result$functions[[j]]@argvals[[1]]
    values <- mfpca_result$functions[[j]]@X[i, ]
    plot(time, values, type = "l", col = axis_colors[j],
         main = paste0("PC", i, " - ", axis_labels[j], " axis"),
         xlab = if(j==2) "Time" else "", ylab = "Eigenfunction", lwd = 2)
    abline(h = 0, lty = 2, col = "gray")
  }
}
par(mfrow = c(1,1))

# Mode of Variation Plot 
plot_modes_of_variation <- function(reg_fd, mfpca_result, pc = 2, k = 2) {
  time <- mfpca_result$functions[[1]]@argvals[[1]]
  mean_x <- eval.fd(time, mean.fd(reg_fd$x))
  mean_y <- eval.fd(time, mean.fd(reg_fd$y))
  mean_z <- eval.fd(time, mean.fd(reg_fd$z))
  sdev <- sqrt(mfpca_result$values[pc])
  axis_labels <- c("X", "Y", "Z")
  pos_colors <- c("steelblue", "tomato", "seagreen")
  neg_colors <- c("dodgerblue4", "firebrick", "forestgreen")
  par(mfrow=c(3,1), mar=c(4,4,2,1))
  for (j in 1:3) {
    ef <- mfpca_result$functions[[j]]@X[pc, ]
    mean_axis <- switch(j, mean_x, mean_y, mean_z)
    lims <- range(
      mean_axis + k * sdev * ef,
      mean_axis - k * sdev * ef,
      mean_axis
    )
    plot(time, mean_axis, type = "l", col = "black", lwd = 2, ylim = rev(lims),
         main = paste0("PC", pc, " - ", axis_labels[j], " axis"),
         ylab = "Value", xlab = "Time")
    lines(time, mean_axis + k * sdev * ef, col = pos_colors[j], lty = 2, lwd = 2)
    lines(time, mean_axis - k * sdev * ef, col = neg_colors[j], lty = 2, lwd = 2)
    legend("topright",
           legend = c("Mean", "+2 SD", "-2 SD"),
           col = c("black", pos_colors[j], neg_colors[j]),
           lty = c(1, 2, 2), lwd = 2)
  }
  par(mfrow = c(1,1))
}
plot_modes_of_variation(reg_fd, mfpca_result, pc = 1, k = 2)


# Mean FPC Scores Per Subject 
scores_df <- data.frame(
  Subject_ID = subject_ids,
  mfpca_result$scores
)
colnames(scores_df)[-1] <- paste0("FPC", seq_len(ncol(mfpca_result$scores)))

subject_fpc_means <- aggregate(. ~ Subject_ID, data = scores_df, FUN = mean)
print(subject_fpc_means)


barplot(
  t(axis_contributions),                           # transpose for stacking
  beside = FALSE,                                 # stacked, not grouped
  col = c("steelblue", "tomato", "seagreen"),     # colors for X, Y, Z
  names.arg = paste0("PC", 1:nrow(axis_contributions)),
  ylab = "Axis Contribution (Proportion)",
  main = "Axis Contribution for Each FPC",
  legend = c("X", "Y", "Z"),
  args.legend = list(x = "topright", bty = "n", inset = 0.01)
)
# Set margins so legend doesn't overlap bars
par(mar = c(5, 4, 6, 2))  # extra space at top for legend

bar_heights <- t(axis_contributions)
bp <- barplot(
  bar_heights,                           # transpose for stacking
  beside = FALSE,                        # stacked, not grouped
  col = c("steelblue", "tomato", "seagreen"), # colors for X, Y, Z
  names.arg = paste0("PC", 1:nrow(axis_contributions)),
  ylab = "Axis Contribution (Proportion)",
  main = "Axis Contribution for Each FPC"
)

# Add legend above plot
legend("top", legend = c("X", "Y", "Z"),
       fill = c("steelblue", "tomato", "seagreen"),
       bty = "n", horiz = TRUE, inset = c(0, -0.1), xpd = TRUE)

# Add percentage labels in the middle of each stacked segment
cum_heights <- apply(bar_heights, 2, cumsum)
mid_heights <- rbind(
  bar_heights[1, ] / 2,
  cum_heights[1, ] + bar_heights[2, ] / 2,
  cum_heights[2, ] + bar_heights[3, ] / 2
)

for (i in 1:3) { # X, Y, Z
  text(
    x = bp, 
    y = mid_heights[i, ],
    labels = sprintf("%.2f", bar_heights[i, ]),
    cex = 0.8, col = "white"
  )
}



library(rgl)

plot_mfpca_mode3D <- function(reg_fd, mfpca_result, pc = 1, k = 2) {
  time <- mfpca_result$functions[[1]]@argvals[[1]]
  mean_x <- eval.fd(time, mean.fd(reg_fd$x))
  mean_y <- eval.fd(time, mean.fd(reg_fd$y))
  mean_z <- eval.fd(time, mean.fd(reg_fd$z))
  sdev <- sqrt(mfpca_result$values[pc])
  
  ef_x <- mfpca_result$functions[[1]]@X[pc, ]
  ef_y <- mfpca_result$functions[[2]]@X[pc, ]
  ef_z <- mfpca_result$functions[[3]]@X[pc, ]
  
  # Mean
  curve_mean <- cbind(mean_x, mean_y, mean_z)
  # Mean + k SD * eigenfunction
  curve_plus <- cbind(mean_x + k*sdev*ef_x, mean_y + k*sdev*ef_y, mean_z + k*sdev*ef_z)
  # Mean - k SD * eigenfunction
  curve_minus <- cbind(mean_x - k*sdev*ef_x, mean_y - k*sdev*ef_y, mean_z - k*sdev*ef_z)
  
  # Open 3D window and plot
  open3d()
  plot3d(curve_mean, type = "l", lwd = 3, col = "black", xlab = "X", ylab = "Y", zlab = "Z",
         main = paste0("3D Mode of Variation: PC", pc))
  lines3d(curve_plus, col = "blue", lwd = 2)
  lines3d(curve_minus, col = "red", lwd = 2)
  legend3d("topright", legend = c("Mean", "+2SD", "-2SD"), col = c("black", "blue", "red"), lty = 1, lwd = 2)
}
plot_mfpca_mode3D(reg_fd, mfpca_result, pc = 1, k = 2)



