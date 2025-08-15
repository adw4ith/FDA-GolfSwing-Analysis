library(R.matlab)
library(fda)
library(here)

# -------- Setup --------
subjects <- c(4,5,7,9,10:27,29:49,52:57,59:63)
subjects <- paste0("00", subjects)
subjects[1:4] <- paste0("0", subjects[1:4])
subjects <- paste0("AMID", subjects)
npts <- 1498
unitime <- seq(0, 1, length.out = npts)
nbasis <- 30
norder <- 6
basis_fd <- create.bspline.basis(rangeval = c(0, 1), nbasis = nbasis, norder = norder)

# Lists to save smoothed curves
all_smooth_x <- list()
all_smooth_y <- list()
all_smooth_z <- list()

for (s in seq_along(subjects)) {
  subject_file <- paste0(subjects[s], ".mat")
  mat_file_path <- here("Data", "Mat files", subject_file)
  raw <- readMat(mat_file_path)
  shot_nums <- as.vector(raw$Data[[2]])
  trajectories <- raw$Data[[4]]
  n_shots <- length(shot_nums)
  
  # --- Build array of interpolated swings ---
  traj_len <- sapply(shot_nums, function(i) length(unlist(trajectories[[i]])))
  max_points <- max(traj_len) / 3
  uniswing <- array(NA_real_, dim = c(npts, n_shots, 3))
  
  for (j in seq_along(shot_nums)) {
    coords <- unlist(trajectories[[shot_nums[j]]])
    len <- length(coords)
    if (len %% 3 != 0) next
    t0 <- seq(0, 1, length.out = len / 3)
    for (k in 1:3) {
      raw_k <- coords[seq(k, len, by = 3)]
      uniswing[, j, k] <- approx(t0, raw_k, xout = unitime, rule = 2)$y
    }
  }
  
  # --- Smoothing with GCV (per subject) ---
  golfLoglam <- seq(-25, 5, by = 1)
  gcvStats <- data.frame(log10.lambda = golfLoglam, df = NA_real_, gcv = NA_real_)
  for (i in seq_along(golfLoglam)) {
    lambda <- 10^golfLoglam[i]
    fdpar <- fdPar(basis_fd, 2, lambda)
    tryCatch({
      sr <- smooth.basis(unitime, uniswing, fdParobj = fdpar)
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
  fdpar_opt <- fdPar(basis_fd, 2, lambda_opt)
  smooth_opt <- smooth.basis(unitime, uniswing, fdParobj = fdpar_opt)
  
  # Save smoothed X, Y, Z fd objects for this subject
  all_smooth_x[[s]] <- smooth_opt$fd[, 1]
  all_smooth_y[[s]] <- smooth_opt$fd[, 2]
  all_smooth_z[[s]] <- smooth_opt$fd[, 3]
}


# Save to RDS for later use if you want
saveRDS(all_smooth_x, "smoothed_x_all_subjects.rds")
saveRDS(all_smooth_y, "smoothed_y_all_subjects.rds")
saveRDS(all_smooth_z, "smoothed_z_all_subjects.rds")


all_smooth_x <- readRDS("smoothed_x_all_subjects.rds")
all_smooth_y <- readRDS("smoothed_y_all_subjects.rds")
all_smooth_z <- readRDS("smoothed_z_all_subjects.rds")
n_subjects <- length(all_smooth_x)
unitime <- seq(0, 1, length.out = 1498)

library(dplyr)
library(tidyr)
library(ggplot2)

# Stack all curves into one big matrix for each axis
get_all_curves <- function(fd_list) {
  do.call(cbind, lapply(fd_list, function(fdobj) eval.fd(unitime, fdobj)))
}

all_x_mat <- get_all_curves(all_smooth_x)
all_y_mat <- get_all_curves(all_smooth_y)
all_z_mat <- get_all_curves(all_smooth_z)

df_x <- as.data.frame(all_x_mat)
df_x$time <- unitime
df_x_long <- pivot_longer(df_x, -time, names_to = "curve", values_to = "amplitude")

ggplot(df_x_long, aes(x = time, y = amplitude, group = curve)) +
  geom_line(alpha = 0.2) +
  labs(title = "All Smoothed X Curves", x = "Normalized Time", y = "X Amplitude")

mean_x <- rowMeans(all_x_mat)
sd_x <- apply(all_x_mat, 1, sd)

plot(unitime, mean_x, type = "l", lwd = 2, col = "blue", ylim = range(mean_x + 2*sd_x, mean_x - 2*sd_x),
     main = "Mean X with ±1 SD Band", xlab = "Normalized Time", ylab = "X Amplitude")
lines(unitime, mean_x + sd_x, col = "black", lty = 2)
lines(unitime, mean_x - sd_x, col = "black", lty = 2)


subject_means_x <- sapply(all_smooth_x, function(fdobj) rowMeans(eval.fd(unitime, fdobj)))

matplot(unitime, subject_means_x, type = "l", lty = 1, col = rainbow(ncol(subject_means_x)),
        xlab = "Normalized Time", ylab = "Mean Y Amplitude", main = "Subject Mean Y Curves")


min_time <- apply(all_x_mat, 2, function(curve) unitime[which.min(curve)])
max_time <- apply(all_x_mat, 2, function(curve) unitime[which.max(curve)])

hist(min_time, breaks = 50, main = "Histogram of Min Timing (X)", xlab = "Normalized Time")
hist(max_time, breaks = 50, main = "Histogram of Max Timing (X)", xlab = "Normalized Time")

mean_curve <- mean_x
curve_dist <- apply(all_x_mat, 2, function(curve) sqrt(mean((curve - mean_curve)^2)))
outlier_idx <- order(curve_dist, decreasing = TRUE)[1:5]

matplot(unitime, all_x_mat[, outlier_idx], type = "l", lty = 1, col = "red",
        main = "Most Distant X Curves (Potential Outliers)", xlab = "Normalized Time", ylab = "X Amplitude")
lines(unitime, mean_curve, col = "blue", lwd = 2)


library(readxl)

# Read the full sheet
full_data <- readxl::read_excel("Data/Player data.xlsx", sheet = "Study 2 data", col_names = FALSE)

# Group 1 columns: 2-6 (B-F); Group 2 columns: 8-12 (H-L)
group1 <- full_data[3:nrow(full_data), c(2,3,4,5,6)]
group2 <- full_data[3:nrow(full_data), c(8,9,10,11,12)]

colnames(group1) <- c("Player_ID", "Hcp", "Age", "Height", "Mass")
colnames(group2) <- c("Player_ID", "Hcp", "Age", "Height", "Mass")

group1$skill_group <- "High"
group2$skill_group <- "Low"

# Remove rows where Player_ID is NA
group1 <- group1[!is.na(group1$Player_ID), ]
group2 <- group2[!is.na(group2$Player_ID), ]

player_data <- bind_rows(group1, group2)

subjects <- c(4,5,7,9,10:27,29:49,52:57,59:63)
subjects <- paste0("00", subjects)
subjects[1:4] <- paste0("0", subjects[1:4])
subjects <- paste0("AMID", subjects)

subject_df <- data.frame(subject_index = seq_along(subjects), Player_ID = subjects)
players <- merge(subject_df, player_data, by = "Player_ID", all.x = TRUE)

all_smooth_x <- readRDS("smoothed_x_all_subjects.rds")
unitime <- seq(0, 1, length.out = 1498)

basisfd <- all_smooth_x[[1]]$basis
coef_list <- lapply(all_smooth_x, function(fd) fd$coefs)  # each: nbasis x ncurves_subject
coef_all  <- do.call(cbind, coef_list)                    # nbasis x total_curves
x_fd_combined <- fd(coef_all, basisfd)
col_meta <- do.call(
  rbind,
  lapply(seq_along(all_smooth_x), function(i) {
    nci <- ncol(all_smooth_x[[i]]$coefs)
    data.frame(
      Player_ID = rep(subjects[i], nci),
      curve     = seq_len(nci),
      stringsAsFactors = FALSE
    )
  })
) %>%
  mutate(
    col_idx = row_number(),
    skill_group = players$skill_group[match(Player_ID, players$Player_ID)]
  )

# Mean template
mean_x_fd <- mean.fd(x_fd_combined)

# Warping model: smooth monotone time warps; tweak nbasis/lambda as needed
warp_basis <- create.bspline.basis(rangeval = c(0, 1), nbasis = 7)
WfdPar     <- fdPar(warp_basis, Lfdobj = 2, lambda = 1e-2)

registered_result <- register.fd(
  y0fd     = mean_x_fd,
  yfd      = x_fd_combined,
  conv     = 1e-10,
  crit     = 2   # L2 criterion
)

fd_reg <- registered_result$regfd

# Evaluate and plot by skill 
x_reg_mat <- eval.fd(unitime, fd_reg)  # matrix

df_x_reg <- data.frame(
  time      = rep(unitime, ncol(x_reg_mat)),
  amplitude = as.vector(x_reg_mat),
  col_idx   = rep(seq_len(ncol(x_reg_mat)), each = length(unitime))
) %>%
  left_join(col_meta, by = "col_idx")

ggplot(df_x_reg,
       aes(x = time, y = amplitude,
           group = interaction(Player_ID, curve),
           color = skill_group)) +
  geom_line(alpha = 0.15) +
  labs(title = "All Registered X Curves by Skill Group",
       x = "Normalized Time", y = "X Amplitude") +
  scale_color_manual(values = c("High" = "blue", "Low" = "red"), na.value = "grey50") +
  theme_minimal()



# Evaluate registered curves and their first derivatives
X  <- eval.fd(unitime, fd_reg)                     # time x curves
Xd <- eval.fd(unitime, deriv.fd(fd_reg, 1))        # derivative

# --- focus on the second half (t in [0.6, 0.95]; tweak if needed) ---
t1 <- 0.60; t2 <- 0.95
idx <- which(unitime >= t1 & unitime <= t2)

# Robust "majority" reference in that window
center_seg <- apply(X[idx, , drop = FALSE], 1, median, na.rm = TRUE)

# Per-curve correlation to the majority pattern in the window
cors <- sapply(seq_len(ncol(X)), function(j) {
  y <- X[idx, j]
  if (sd(y) == 0 || sd(center_seg) == 0) return(NA_real_)
  cor(y, center_seg, use = "complete.obs")
})

# Per-curve mean slope in the window
mean_slope <- colMeans(Xd[idx, , drop = FALSE], na.rm = TRUE)

# Majority slope sign (usually positive here)
maj_sign <- sign(median(mean_slope, na.rm = TRUE))
if (maj_sign == 0) maj_sign <- sign(mean(mean_slope, na.rm = TRUE))

# Flag “opposite” curves:
#  - correlation negative (pattern is opposite),
#  - and slope opposite to the majority.
opp_idx <- which((cors < 0) & (mean_slope * maj_sign < 0))

# Mark and drop only those
col_meta$opp <- col_meta$col_idx %in% opp_idx

X_keep <- X[, !col_meta$opp, drop = FALSE]
keep_meta <- col_meta[!col_meta$opp, ]

# Tidy for plotting
df_clean <- data.frame(
  time      = rep(unitime, ncol(X_keep)),
  amplitude = as.vector(X_keep),
  col_idx   = rep(keep_meta$col_idx, each = length(unitime))
) %>%
  left_join(keep_meta, by = "col_idx")

ggplot(df_clean, aes(time, amplitude,
                     group = interaction(Player_ID, curve),
                     color = skill_group)) +
  geom_line(alpha = 0.15) +
  labs(title = "All X Curves by Skill Group ",
       x = "Normalized Time", y = "X Amplitude") +
  scale_color_manual(values = c(High = "blue", Low = "red"), na.value = "grey50") +
  theme_minimal()

#  All X curves by skill group (separate plots) 

# Choose the registered source data frame you already built
df_src <- if (exists("df_clean")) {
  df_clean
} else if (exists("df_x_reg_clean")) {
  df_x_reg_clean
} else {
  df_x_reg
}

# Build palettes per group (one color per Player_ID)
high_ids <- df_src %>% filter(skill_group == "High") %>% pull(Player_ID) %>% unique() %>% sort()
low_ids  <- df_src %>% filter(skill_group == "Low")  %>% pull(Player_ID) %>% unique() %>% sort()

pal_high <- setNames(colorRampPalette(c("#1b6ef3", "#a8c5ff"))(length(high_ids)), high_ids)
pal_low  <- setNames(colorRampPalette(c("#e31a1c", "#f7b6b2"))(length(low_ids)),  low_ids)

# High group: all curves
ggplot(
  df_src %>%
    filter(skill_group == "High") %>%
    mutate(Player_ID = factor(Player_ID, levels = high_ids)),
  aes(x = time, y = amplitude,
      group = interaction(Player_ID, curve),
      color = Player_ID)
) +
  geom_line(alpha = 0.15) +
  scale_color_manual(values = pal_high) +
  labs(title = "All Registered X Curves – High Skill",
       x = "Normalized Time", y = "X Amplitude") +
  theme_minimal() +
  theme(legend.position = "none")

# Low group: all curves
ggplot(
  df_src %>%
    filter(skill_group == "Low") %>%
    mutate(Player_ID = factor(Player_ID, levels = low_ids)),
  aes(x = time, y = amplitude,
      group = interaction(Player_ID, curve),
      color = Player_ID)
) +
  geom_line(alpha = 0.15) +
  scale_color_manual(values = pal_low) +
  labs(title = "All Registered X Curves – Low Skill",
       x = "Normalized Time", y = "X Amplitude") +
  theme_minimal() +
  theme(legend.position = "none")



all_x <- lapply(seq_along(all_smooth_x), function(i) {
  x_mat <- eval.fd(unitime, all_smooth_x[[i]])
  if (is.null(dim(x_mat))) x_mat <- matrix(x_mat, ncol = 1)
  data.frame(
    time = rep(unitime, ncol(x_mat)),
    amplitude = as.vector(x_mat),
    curve = rep(seq_len(ncol(x_mat)), each = length(unitime)),
    Player_ID = subjects[i],
    skill_group = players$skill_group[i]
  )
})
df_x <- bind_rows(all_x)

# Overlay plot
ggplot(df_x, aes(x = time, y = amplitude, group = interaction(Player_ID, curve), color = skill_group)) +
  geom_line(alpha = 0.15) +
  labs(title = "All Smoothed X Curves by Skill Group", x = "Normalized Time", y = "X Amplitude") +
  scale_color_manual(values = c("High" = "blue", "Low" = "red")) +
  theme_minimal()

# Mean curves by group
df_group_mean <- df_x %>%
  group_by(time, skill_group) %>%
  summarise(mean_amp = mean(amplitude, na.rm = TRUE), .groups = "drop")

ggplot(df_group_mean, aes(x = time, y = mean_amp, color = skill_group)) +
  geom_line(linewidth = 1.5) +
  labs(title = "Mean X Curve by Skill Group", x = "Normalized Time", y = "Mean X Amplitude") +
  scale_color_manual(values = c("High" = "blue", "Low" = "red")) +
  theme_minimal()

# For Y
all_smooth_y <- readRDS("smoothed_y_all_subjects.rds")

all_y <- lapply(seq_along(all_smooth_y), function(i) {
  y_mat <- eval.fd(unitime, all_smooth_y[[i]])
  if (is.null(dim(y_mat))) y_mat <- matrix(y_mat, ncol = 1)
  data.frame(
    time = rep(unitime, ncol(y_mat)),
    amplitude = as.vector(y_mat),
    curve = rep(seq_len(ncol(y_mat)), each = length(unitime)),
    Player_ID = subjects[i],
    skill_group = players$skill_group[i]
  )
})
df_y <- bind_rows(all_y)

df_group_mean_y <- df_y %>%
  group_by(time, skill_group) %>%
  summarise(mean_amp = mean(amplitude, na.rm = TRUE), .groups = "drop")

ggplot(df_group_mean_y, aes(x = time, y = mean_amp, color = skill_group)) +
  geom_line(size = 1.5) +
  labs(title = "Mean Y Curve by Skill Group", x = "Normalized Time", y = "Mean Y Amplitude") +
  scale_color_manual(values = c("High" = "blue", "Low" = "red")) +
  theme_minimal()

# For Z
all_smooth_z <- readRDS("smoothed_z_all_subjects.rds")

all_z <- lapply(seq_along(all_smooth_z), function(i) {
  z_mat <- eval.fd(unitime, all_smooth_z[[i]])
  if (is.null(dim(z_mat))) z_mat <- matrix(z_mat, ncol = 1)
  data.frame(
    time = rep(unitime, ncol(z_mat)),
    amplitude = as.vector(z_mat),
    curve = rep(seq_len(ncol(z_mat)), each = length(unitime)),
    Player_ID = subjects[i],
    skill_group = players$skill_group[i]
  )
})
df_z <- bind_rows(all_z)

df_group_mean_z <- df_z %>%
  group_by(time, skill_group) %>%
  summarise(mean_amp = mean(amplitude, na.rm = TRUE), .groups = "drop")

ggplot(df_group_mean_z, aes(x = time, y = mean_amp, color = skill_group)) +
  geom_line(size = 1.5) +
  labs(title = "Mean Z Curve by Skill Group", x = "Normalized Time", y = "Mean Z Amplitude") +
  scale_color_manual(values = c("High" = "blue", "Low" = "red")) +
  theme_minimal()

# Compute mean X for each subject
subject_means_x <- lapply(seq_along(all_smooth_x), function(i) {
  x_mat <- eval.fd(unitime, all_smooth_x[[i]])
  if (is.null(dim(x_mat))) x_mat <- matrix(x_mat, ncol = 1)
  data.frame(
    time = unitime,
    mean_x = rowMeans(x_mat),
    Player_ID = subjects[i],
    skill_group = players$skill_group[i]
  )
})
df_subject_mean_x <- bind_rows(subject_means_x)

# Plot High group
ggplot(
  filter(df_subject_mean_x, skill_group == "High"),
  aes(x = time, y = mean_x, group = Player_ID, color = Player_ID)
) +
  geom_line(alpha = 0.8) +
  labs(
    title = "Subject Mean X Curves (High Skill Group)",
    x = "Normalized Time", y = "Mean X Amplitude"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(
  filter(df_subject_mean_x, skill_group == "Low"),
  aes(x = time, y = mean_x, group = Player_ID, color = Player_ID)
) +
  geom_line(alpha = 0.8) +
  labs(
    title = "Subject Mean X Curves (Low Skill Group)",
    x = "Normalized Time", y = "Mean X Amplitude"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

subject_means_y <- lapply(seq_along(all_smooth_y), function(i) {
  y_mat <- eval.fd(unitime, all_smooth_y[[i]])
  if (is.null(dim(y_mat))) y_mat <- matrix(y_mat, ncol = 1)
  data.frame(
    time = unitime,
    mean_y = rowMeans(y_mat),
    Player_ID = subjects[i],
    skill_group = players$skill_group[i]
  )
})
df_subject_mean_y <- bind_rows(subject_means_y)

# High
ggplot(
  filter(df_subject_mean_y, skill_group == "High"),
  aes(x = time, y = mean_y, group = Player_ID, color = Player_ID)
) +
  geom_line(alpha = 0.8) +
  labs(
    title = "Subject Mean Y Curves (High Skill Group)",
    x = "Normalized Time", y = "Mean Y Amplitude"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


#low
ggplot(
  filter(df_subject_mean_y, skill_group == "Low"),
  aes(x = time, y = mean_y, group = Player_ID, color = Player_ID)
) +
  geom_line(alpha = 0.8) +
  labs(
    title = "Subject Mean Y Curves (Low Skill Group)",
    x = "Normalized Time", y = "Mean Y Amplitude"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


subject_means_z <- lapply(seq_along(all_smooth_z), function(i) {
  z_mat <- eval.fd(unitime, all_smooth_z[[i]])
  if (is.null(dim(z_mat))) z_mat <- matrix(z_mat, ncol = 1)
  data.frame(
    time = unitime,
    mean_z = rowMeans(z_mat),
    Player_ID = subjects[i],
    skill_group = players$skill_group[i]
  )
})
df_subject_mean_z <- bind_rows(subject_means_z)

# High
ggplot(
  filter(df_subject_mean_z, skill_group == "High"),
  aes(x = time, y = mean_z, group = Player_ID, color = Player_ID)
) +
  geom_line(alpha = 0.8) +
  labs(
    title = "Subject Mean Z Curves (High Skill Group)",
    x = "Normalized Time", y = "Mean Z Amplitude"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggplot(
  filter(df_subject_mean_z, skill_group == "Low"),
  aes(x = time, y = mean_z, group = Player_ID, color = Player_ID)
) +
  geom_line(alpha = 0.8) +
  labs(
    title = "Subject Mean Z Curves (Low Skill Group)",
    x = "Normalized Time", y = "Mean Z Amplitude"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


# Assume df_x was built as in earlier code:

ggplot(
  filter(df_x, skill_group == "High"),
  aes(x = time, y = amplitude, group = interaction(Player_ID, curve),color = "blue")
) +
  geom_line(alpha = 0.15) +
  labs(
    title = "All X Curves (High Skill Group)",
    x = "Normalized Time", y = "X Amplitude"
  ) +
  theme_minimal() +
  theme(legend.position = "none")


ggplot(
  filter(df_x, skill_group == "Low"),
  aes(x = time, y = amplitude, group = interaction(Player_ID, curve), color = "blue")
) +
  geom_line(alpha = 0.15) +
  labs(
    title = "All X Curves (Low Skill Group)",
    x = "Normalized Time", y = "X Amplitude"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

