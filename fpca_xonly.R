library(fda)
library(R.matlab)
library(here)

# Setup
subjects <- c(4,5,7,9,10:27,29:49,52:57,59:63)
subjects <- paste0("00", subjects)
subjects[1:4] <- paste0("0", subjects[1:4])
subjects <- paste0("AMID", subjects)

npts <- 1498
unitime <- seq(0, 1, length.out = npts)
nbasis <- 30
norder <- 6
basis_fd <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis, norder = norder)
fdpar_opt <- fdPar(fdobj = basis_fd, Lfdobj = 2, lambda = 1e-6)  # use same optimal lambda if fixed

# Store all X smooth curves here
all_x_fd_list <- list()

for (s in seq_along(subjects)) {
  subject_file <- paste0(subjects[s], ".mat")
  mat_file_path <- here("Data", "Mat files", subject_file)
  
  raw <- readMat(mat_file_path)
  shot_nums <- as.vector(raw$Data[[2]])
  trajectories <- raw$Data[[4]]
  
  n_shots <- length(shot_nums)
  traj_length <- numeric(n_shots)
  
  for (i in seq_len(n_shots)) {
    traj_i <- unlist(trajectories[[shot_nums[i]]])
    traj_length[i] <- length(traj_i)
  }
  
  max_length <- max(traj_length)
  max_points <- max_length / 3
  
  subject_trajectories <- array(NA_real_, dim = c(max_points, n_shots, 3))
  
  for (i in seq_len(n_shots)) {
    traj_i <- unlist(trajectories[[shot_nums[i]]])
    len_i <- length(traj_i)
    dimlen <- len_i / 3
    
    if (abs(dimlen - as.integer(dimlen)) > 1e-10) next
    
    dim1 <- traj_i[seq(1, len_i, by = 3)]
    dim2 <- traj_i[seq(2, len_i, by = 3)]
    dim3 <- traj_i[seq(3, len_i, by = 3)]
    
    pad <- max_points - length(dim1)
    if (pad > 0) {
      dim1 <- c(dim1, rep(NA, pad))
      dim2 <- c(dim2, rep(NA, pad))
      dim3 <- c(dim3, rep(NA, pad))
    }
    
    subject_trajectories[, i, 1] <- dim1
    subject_trajectories[, i, 2] <- dim2
    subject_trajectories[, i, 3] <- dim3
  }
  
  # Interpolation
  uniswing <- array(NA_real_, dim = c(npts, n_shots, 3))
  for (j in seq_len(n_shots)) {
    t0 <- seq(0, 1, length.out = max_points)
    for (k in 1:3) {
      if (all(is.na(subject_trajectories[, j, k]))) next
      uniswing[, j, k] <- approx(x = t0, y = subject_trajectories[, j, k], xout = unitime, rule = 2)$y
    }
  }
  
  # Smooth and collect only valid X
  try({
    smooth_opt <- smooth.basis(argvals = unitime, y = uniswing[,,1], fdParobj = fdpar_opt)
    x_fd <- smooth_opt$fd
    
    # Remove shot 11 for subject 3 only
    if (s == 3) {
      x_fd <- fd(coef = x_fd$coefs[, -11], basisobj = x_fd$basis)
    }
    
    all_x_fd_list[[length(all_x_fd_list) + 1]] <- x_fd
  }, silent = TRUE)
}

# Combine all X curves
x_fd_combined <- do.call(fd, list(
  coef = do.call(cbind, lapply(all_x_fd_list, function(fdobj) fdobj$coefs)),
  basisobj = basis_fd
))

# Perform FPCA on all X
pca_x <- pca.fd(x_fd_combined, nharm = 4, centerfns = TRUE)
# Extract variance proportions
var_explained <- pca_x$varprop * 100  # convert to percentage
harmonics <- seq_along(var_explained)

# Make scree plot
plot(harmonics, var_explained,
     type = "b", pch = 19, lwd = 2,
     xlab = "Principal Component",
     ylab = "Variance Explained (%)",
     main = "Scree Plot - FPCA",
     ylim = c(0, max(var_explained) * 1.1))

# Add percentage labels above points
text(harmonics, var_explained,
     labels = sprintf("%.1f%%", var_explained),
     pos = 3, cex = 0.8)

# Plot eigenfunctions
plot.pca.fd(pca_x, harm = 1:3)
plot.pca.fd(pca_x, harm = 1, expand = 2)  # makes the differences more obvious

obj <- readRDS("registered_xyz_fd.rds")
class(obj)
# Load registered XYZ curves
registered_xyz_fd <- readRDS("registered_xyz_fd.rds")

# Extract X-axis registered curves
x_fd <- registered_xyz_fd$x


# Compute mean curve
mean_x_fd <- mean.fd(x_fd)

# Evaluate curves at equally spaced points
unitime <- seq(x_fd$basis$rangeval[1], x_fd$basis$rangeval[2], length.out = 200)
x_mat <- eval.fd(unitime, x_fd)  # rows: time points, cols: curves

# Identify outliers: mean displacement in second half < 0
second_half_idx <- which(unitime > 0.5)
curve_means_second_half <- colMeans(x_mat[second_half_idx, ], na.rm = TRUE)

# Logical index of curves to keep
keep_idx <- curve_means_second_half > 0

# Filter fd object
x_fd_clean <- x_fd[keep_idx]

# Recompute mean on cleaned data
mean_x_fd_clean <- mean.fd(x_fd_clean)

# Plot cleaned curves
plot(x_fd_clean,
     col = rgb(0.5, 0.5, 0.5, 0.15),
     lty = 1,
     xlab = "Normalized Time",
     ylab = "X-axis displacement",
     main = "Registered X-axis Curves with Mean Function")

# Overlay mean curve
lines(mean_x_fd_clean, col = "black", lwd = 3)

# Check count
removed_count <- sum(!keep_idx)
cat("Removed", removed_count, "outlier curves\n")

# Plot all curves 
plot(x_fd,
     col = rgb(0.5, 0.5, 0.5, 0.15),  # low transparency
     lty = 1,
     xlab = "Normalized Time",
     ylab = "X-axis displacement",
     main = "Registered X-axis Curves with Mean Function")

# Overlay mean curve in bold red
lines(mean_x_fd, col = "black", lwd = 3)

# Variance explained
round(pca_x$varprop * 100, 2)

# 1. Compute mean function across all X-axis curves
mean_x_fd <- mean.fd(x_fd_combined)

# 2. Perform continuous registration on all X curves
registered_result <- register.fd(y0fd = mean_x_fd, yfd = x_fd_combined, conv = 1e-10, crit = 2)

# 3. Extract the registered X curves
x_fd_registered <- registered_result$regfd  # functional data object

# 4. Run FPCA on registered data
pca_x_reg <- pca.fd(x_fd_registered, nharm = 4, centerfns = TRUE)

# 5. Plot results
plot.pca.fd(pca_x_reg, harm = 1:3)
plot.pca.fd(pca_x_reg, harm = 1, expand = 2)

# 6. Show variance explained
round(pca_x_reg$varprop * 100, 2)
# Evaluate the mean function over time
mean_x_reg <- eval.fd(unitime, pca_x_reg$meanfd)

# Evaluate the first harmonic (FPC1) over time
fpc1_reg <- eval.fd(unitime, pca_x_reg$harmonics[1])



# Evaluate the mean function over time

plus_curve_reg <- mean_x_reg + expand * fpc1_reg
minus_curve_reg <- mean_x_reg - expand * fpc1_reg
plot(unitime, mean_x_reg, type = "l", lwd = 2, col = "black",
     ylim = range(c(plus_curve_reg, minus_curve_reg)),
     ylab = "X Position", xlab = "Normalized Time",
     main = paste("Mean ±", expand, "× FPC1"))
lines(unitime, plus_curve_reg, col = "blue", lwd = 2, lty = 2)
lines(unitime, minus_curve_reg, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Mean", "Mean + FPC1", "Mean - FPC1"),
       col = c("black", "blue", "red"), lwd = 2, lty = c(1, 2, 2))

pca_x_reg$scores

# Step 3: Get FPC1 scores
fpc1_scores_reg <- pca_x_reg$scores[, 1]

# Step 4: Associate scores with subjects
library(dplyr)

# Combine into a data frame
fpca_df_reg <- data.frame(
  Subject_reg = swing_subject_map,
  FPC1_Score_reg = fpc1_scores_reg
)

# Step 5: Summarize by subject
summary_by_subject_reg <- fpca_df_reg %>%
  group_by(Subject_reg) %>%
  summarise(
    Mean_FPC1 = mean(FPC1_Score_reg),
    SD_FPC1 = sd(FPC1_Score_reg),
    N_Swings = n()
  )

print(summary_by_subject_reg,n=54)

mean_x <- eval.fd(unitime, pca_x$meanfd)
# Evaluate the first harmonic (FPC1) over time
fpc1 <- eval.fd(unitime, pca_x$harmonics[1])

# Choose a multiplier to exaggerate the variation
expand <- 2
# Calculate mean ± FPC1
plus_curve <- mean_x + expand * fpc1
minus_curve <- mean_x - expand * fpc1

# Plot
plot(unitime, mean_x, type = "l", lwd = 2, col = "black",
     ylim = range(c(plus_curve, minus_curve)),
     ylab = "X Position", xlab = "Normalized Time",
     main = paste("Mean ±", expand, "× FPC1"))
lines(unitime, plus_curve, col = "blue", lwd = 2, lty = 2)
lines(unitime, minus_curve, col = "red", lwd = 2, lty = 2)
legend("topright", legend = c("Mean", "Mean + FPC1", "Mean - FPC1"),
       col = c("black", "blue", "red"), lwd = 2, lty = c(1, 2, 2))

pca_x$scores


# Step 1: Count how many swings per subject 
swings_per_subject <- sapply(all_x_fd_list, function(fdobj) {
  if (is.null(fdobj)) return(0)
  ncol(fdobj$coefs)
})

# Step 2: Create subject mapping for each swing
swing_subject_map <- rep(subjects, times = swings_per_subject)

# Confirm the structure
table(swing_subject_map)

# Step 3: Get FPC1 scores
fpc1_scores <- pca_x$scores[, 1]

# Step 4: Associate scores with subjects
library(dplyr)

# Combine into a data frame
fpca_df <- data.frame(
  Subject = swing_subject_map,
  FPC1_Score = fpc1_scores
)

# Step 5: Summarize by subject
summary_by_subject <- fpca_df %>%
  group_by(Subject) %>%
  summarise(
    Mean_FPC1 = mean(FPC1_Score),
    SD_FPC1 = sd(FPC1_Score),
    N_Swings = n()
  )

print(summary_by_subject)
boxplot(FPC1_Score ~ Subject, data = fpca_df,
        las = 2,  # Rotate x-axis labels
        col = "lightblue",
        main = "FPC1 Score Distribution by Subject",
        ylab = "FPC1 Score", xlab = "Subject ID")

# Compare FPC2 scores across subjects
fpca_df$FPC2_Score <- pca_x$scores[, 2]

# Summarize by subject
sum <- fpca_df %>%
  group_by(Subject) %>%
  summarise(
    Mean_FPC2 = mean(FPC2_Score),
    SD_FPC2 = sd(FPC2_Score)
  )
print(sum)
print(sum, n=54)
# Evaluate mean ± FPC2
mean_fd <- mean.fd(x_fd_combined)
fpc2_fd <- pca_x$harmonics[2]

mean_eval <- eval.fd(unitime, mean_fd)
fpc2_eval <- eval.fd(unitime, fpc2_fd)

plot(unitime, mean_eval, type = "l", lwd = 2, col = "black",
     ylim = range(mean_eval + 2*fpc2_eval, mean_eval - 2*fpc2_eval),
     main = "Mean ± 2 × FPC2", xlab = "Normalized Time", ylab = "X")
lines(unitime, mean_eval + 2*fpc2_eval, col = "blue", lty = 2, lwd = 2)
lines(unitime, mean_eval - 2*fpc2_eval, col = "red", lty = 2, lwd = 2)
legend("topright", legend = c("Mean", "+2×FPC2", "−2×FPC2"),
       col = c("black", "blue", "red"), lty = c(1,2,2), lwd = 2)

