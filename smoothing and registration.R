# Load packages
library(R.matlab)
library(fda)
library(here)

# ----------------------
# Subject Setup
# ----------------------
subjects <- c(4,5,7,9,10:27,29:49,52:57,59:63)
subjects <- paste0("00", subjects)
subjects[1:4] <- paste0("0", subjects[1:4])
subjects <- paste0("AMID", subjects)

subject_index <- 3 # Choose subject here
subject_file <- paste0(subjects[subject_index], ".mat")
mat_file_path <- here("Data", "Mat files", subject_file)


# Read Trajectories

raw <- readMat(mat_file_path)
shot_nums <- as.vector(raw$Data[[2]])
trajectories <- raw$Data[[4]]


# Process All Trajectories

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
  
  stopifnot(abs(dimlen - as.integer(dimlen)) < 1e-10)
  
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

npts <- 1498
unitime <- seq(0, 1, length.out = npts)
uniswing <- array(NA_real_, dim = c(npts, n_shots, 3), dimnames = list(NULL, NULL, c("X", "Y", "Z")))

for (j in seq_len(n_shots)) {
  t0 <- seq(0, 1, length.out = max_points)
  for (k in 1:3) {
    uniswing[, j, k] <- approx(x = t0, y = subject_trajectories[, j, k], xout = unitime, rule = 2)$y
  }
}



# Smoothing with GCV

nbasis <- 30
norder <- 6
basis_fd <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis, norder = norder)
golfLoglam <- seq(-25, 5, by = 1)
gcvStats <- data.frame(log10.lambda = golfLoglam, df = NA_real_, gcv = NA_real_)

for (i in seq_along(golfLoglam)) {
  lambda <- 10^golfLoglam[i]
  fdpar <- fdPar(fdobj = basis_fd, Lfdobj = 2, lambda = lambda)
  tryCatch({
    sr <- smooth.basis(argvals = unitime, y = uniswing, fdParobj = fdpar)
    gcvStats$df[i] <- sr$df
    gcvStats$gcv[i] <- sum(sr$gcv)
  }, error = function(e) {
    gcvStats$df[i] <- NA
    gcvStats$gcv[i] <- NA
  })
}

realistic_gcv <- gcvStats[gcvStats$df < 80, ]
if (nrow(realistic_gcv) > 0) {
  opt <- realistic_gcv[which.min(realistic_gcv$gcv), ]
} else {
  realistic_gcv <- gcvStats[gcvStats$df < 90, ]
  if (nrow(realistic_gcv) == 0) stop("No good lambda found")
  opt <- realistic_gcv[which.min(realistic_gcv$gcv), ]
}
lambda_opt <- 10^(opt$log10.lambda)

fdpar_opt <- fdPar(fdobj = basis_fd, Lfdobj = 2, lambda = lambda_opt)
smooth_opt <- smooth.basis(argvals = unitime, y = uniswing, fdParobj = fdpar_opt)
golf_fd <- smooth_opt$fd
tmp_fit <- eval.fd(unitime, golf_fd)
mse <- mean((uniswing - tmp_fit)^2, na.rm = TRUE)
message(sprintf("Mean Squared Error: %.6f", mse))
#plot basis fubctions
plot(basis_fd, xlab = "Time (normalized)", ylab = "Basis function value", 
     main = "B-spline Basis Functions (nbasis = 30, order = 6)")


# Plotting Raw vs. Smoothed Curve


# Choose a swing and an axis 
swing_id <- 1
axis <- 1  

# Extract raw and smoothed data
raw_curve <- uniswing[, swing_id, axis]
smoothed_curve <- eval.fd(unitime, golf_fd[axis])[, swing_id]

# Use a deliberately larger lambda for visual demo
#lambda_demo <- 1e-2  # much higher than optimal, to exaggerate smoothing

fdpar_demo <- fdPar(basis_fd, Lfdobj = 2, lambda = lambda_demo)
smooth_demo <- smooth.basis(argvals = unitime, y = uniswing, fdParobj = fdpar_demo)
golf_fd_demo <- smooth_demo$fd

# Extract curves for plotting
smoothed_demo1 <- eval.fd(unitime, golf_fd_demo[axis])[, swing_id]
raw_curve <- uniswing[, swing_id, axis]

# Plot
plot(unitime, raw_curve, type = "l", col = "grey50", lty = 3, lwd = 5,
     ylab = "Position (X)", xlab = "Normalized Time",
   )
lines(unitime, smoothed_demo, col = "red", lwd = 1)
lines(unitime, smoothed_demo1, col = "blue", lwd = 1)

legend("bottomright", legend = c("Raw Data", "Smoothed Curve ","Higher lambda"),
       col = c("grey50", "red","blue"), lty = c(2, 1), lwd = c(1, 1))
# ----------------------
par(mfrow = c(3, 1), mar = c(4, 4, 2, 1))
axes <- c("X", "Y", "Z")
for (k in 1:3) {
  yr <- range(uniswing[,,k], tmp_fit[,,k], na.rm = TRUE)
  plot(unitime, uniswing[, 1, k], type = "n", ylim = yr,
       xlab = "Normalized Time", ylab = paste0(axes[k], "-coordinate"),
       main = paste(axes[k], ": Raw (red) vs Fit (blue)"))
  matlines(unitime, uniswing[,,k], col = adjustcolor("red", alpha.f = 0.5))
  matlines(unitime, tmp_fit[,,k], col = "blue", lwd = 1.2)
}
par(mfrow = c(1, 1))
x <- smooth_opt$fd[,1]
plot(x)

mean_target <- mean.fd(x)

plot(x)
lines(mean_target, lwd = 5)
x_eval <- eval.fd(evalarg = unitime, fdobj = x)
registered_contin <- register.fd(y0fd = mean.fd(smooth_opt$fd), yfd = smooth_opt$fd, conv = 1e-10, crit = 2)

par(mfrow = c(2, 1))

y <- smooth_opt$fd[, 2]
z <- smooth_opt$fd[, 3]
plot(x)
plot(registered_contin$regfd[, 1])
plot(registered_contin$warpfd[, 1])

plot(z)
plot(registered_contin$regfd[, 3])
plot(registered_contin$warpfd[, 3])

plot(y)
plot(registered_contin$regfd[, 2])