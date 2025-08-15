library(fda)
library(R.matlab)
library(here)
reg_fd <- readRDS("registered_xyz_fd.rds")
str(reg_fd)

# Setup

npts <- 1498
unitime <- seq(0, 1, length.out = npts)
nbasis <- 30
norder <- 6
basis_fd <- create.bspline.basis(c(0, 1), nbasis = nbasis, norder = norder)
fdpar_opt <- fdPar(basis_fd, 2, 1e-6)

all_x <- list(); all_y <- list(); all_z <- list()

for (s in seq_along(subjects)) {
  file_path <- here("Data", "Mat files", paste0(subjects[s], ".mat"))
  raw <- readMat(file_path)
  shot_nums <- as.vector(raw$Data[[2]])
  trajectories <- raw$Data[[4]]
  n_shots <- length(shot_nums)
  
  traj_len <- sapply(shot_nums, function(i) length(unlist(trajectories[[i]])))
  max_points <- max(traj_len) / 3
  data_arr <- array(NA, dim = c(npts, n_shots, 3))
  
  for (j in seq_along(shot_nums)) {
    coords <- unlist(trajectories[[shot_nums[j]]])
    len <- length(coords)
    if (len %% 3 != 0) next
    t0 <- seq(0, 1, length.out = len / 3)
    for (k in 1:3) {
      raw_k <- coords[seq(k, len, by = 3)]
      data_arr[, j, k] <- approx(t0, raw_k, xout = unitime, rule = 2)$y
    }
  }
  
  # Smooth
  smooth_x <- smooth.basis(unitime, data_arr[,,1], fdpar_opt)$fd
  smooth_y <- smooth.basis(unitime, data_arr[,,2], fdpar_opt)$fd
  smooth_z <- smooth.basis(unitime, data_arr[,,3], fdpar_opt)$fd
  
  # Drop shot 11 from subject 3
  if (s == 3) {
    smooth_x <- fd(coef = smooth_x$coefs[,-11], basisobj = basis_fd)
    smooth_y <- fd(coef = smooth_y$coefs[,-11], basisobj = basis_fd)
    smooth_z <- fd(coef = smooth_z$coefs[,-11], basisobj = basis_fd)
  }
  
  all_x[[s]] <- smooth_x
  all_y[[s]] <- smooth_y
  all_z[[s]] <- smooth_z
}

# Combine all subjects' data

x_fd_all <- fd(coef = do.call(cbind, lapply(all_x, function(fd) fd$coefs)), basisobj = basis_fd)
y_fd_all <- fd(coef = do.call(cbind, lapply(all_y, function(fd) fd$coefs)), basisobj = basis_fd)
z_fd_all <- fd(coef = do.call(cbind, lapply(all_z, function(fd) fd$coefs)), basisobj = basis_fd)


# Continuous registration

reg_x <- register.fd(mean.fd(x_fd_all), x_fd_all, conv = 1e-10, crit = 2)$regfd
reg_y <- register.fd(mean.fd(y_fd_all), y_fd_all, conv = 1e-10, crit = 2)$regfd
reg_z <- register.fd(mean.fd(z_fd_all), z_fd_all, conv = 1e-10, crit = 2)$regfd


mfd <- list(reg_x, reg_y, reg_z)
# Combine registered X, Y, Z into a single fd object
# reg_xyz_fd <- fd(
#   coef = array(c(reg_x$coefs, reg_y$coefs, reg_z$coefs), 
#                dim = c(dim(reg_x$coefs)[1], dim(reg_x$coefs)[2], 3)),
#   basisobj = basis_fd
# )
reg_xyz_fd <- fd(
  coef = array(
    c(reg_fd$x$coefs, reg_fd$y$coefs, reg_fd$z$coefs),
    dim = c(
      dim(reg_fd$x$coefs)[1],         # n_basis
      dim(reg_fd$x$coefs)[2],         # n_curves
      3                               # xyz
    )
  ),
  basisobj = reg_fd$x$basis  # same basis for all
)



# Save Results

saveRDS(list(x = reg_x, y = reg_y, z = reg_z), file = "registered_xyz_fd.rds")

