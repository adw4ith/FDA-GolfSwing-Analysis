
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(fda)

# Read the Excel file
full_data <- readxl::read_excel("Data/Player data.xlsx", sheet = "Study 2 data", col_names = FALSE)
group1 <- full_data[3:nrow(full_data), c(2,3,4,5,6)]
group2 <- full_data[3:nrow(full_data), c(8,9,10,11,12)]
colnames(group1) <- c("Player_ID", "Hcp", "Age", "Height", "Mass")
colnames(group2) <- c("Player_ID", "Hcp", "Age", "Height", "Mass")
group1$skill_group <- "High"
group2$skill_group <- "Low"
group1 <- group1[!is.na(group1$Player_ID), ]
group2 <- group2[!is.na(group2$Player_ID), ]
player_data <- bind_rows(group1, group2)
player_data <- player_data %>%
  mutate(
    Hcp = as.numeric(Hcp),
    Age = as.numeric(Age),
    Height = as.numeric(Height),
    Mass = as.numeric(Mass)
  )

# 2. Read Smoothed X Curves

all_smooth_x <- readRDS("smoothed_x_all_subjects.rds")
unitime <- seq(0, 1, length.out = 1498)


# 3. Match Subject IDs and Build Tidy Dataframe for All Curves
subjects <- c(4,5,7,9,10:27,29:49,52:57,59:63)
subjects <- paste0("00", subjects)
subjects[1:4] <- paste0("0", subjects[1:4])
subjects <- paste0("AMID", subjects)

all_x <- lapply(seq_along(all_smooth_x), function(i) {
  x_mat <- eval.fd(unitime, all_smooth_x[[i]])
  if (is.null(dim(x_mat))) x_mat <- matrix(x_mat, ncol = 1)
  data.frame(
    time = rep(unitime, ncol(x_mat)),
    amplitude = as.vector(x_mat),
    curve = rep(seq_len(ncol(x_mat)), each = length(unitime)),
    Player_ID = subjects[i]
  )
})
df_x <- bind_rows(all_x)
df_x <- left_join(df_x, player_data, by = "Player_ID")


# Mean & Quantile Curves by Group
df_group_summary <- df_x %>%
  group_by(time, skill_group) %>%
  summarise(
    mean_amp = mean(amplitude, na.rm = TRUE),
    q25 = quantile(amplitude, 0.25, na.rm = TRUE),
    q75 = quantile(amplitude, 0.75, na.rm = TRUE)
  )

ggplot(df_group_summary, aes(x = time, y = mean_amp, color = skill_group, fill = skill_group)) +
  geom_line(size = 1) +
  geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.2, color = NA) +
  labs(title = "Mean and Interquartile Range by Skill Group (X curves)") +
  theme_minimal()

# 5. Curve Features vs Player Attributes (e.g., max X vs Hcp)

curve_features <- df_x %>%
  group_by(Player_ID, curve, skill_group, Hcp, Age, Height, Mass) %>%
  summarise(
    max_x = max(amplitude, na.rm = TRUE),
    time_max_x = time[which.max(amplitude)],
    .groups = "drop"
  )

ggplot(curve_features, aes(x = Hcp, y = max_x, color = skill_group)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "Peak X Value vs Handicap", x = "Handicap", y = "Peak X Amplitude") +
  theme_minimal()

# 6. Subject Mean Curves by Group

subject_means_x <- df_x %>%
  group_by(Player_ID, skill_group, time) %>%
  summarise(mean_x = mean(amplitude, na.rm = TRUE), .groups = "drop")

ggplot(subject_means_x, aes(x = time, y = mean_x, group = Player_ID, color = skill_group)) +
  geom_line(alpha = 0.7) +
  labs(title = "Subject Mean X Curves by Skill Group") +
  theme_minimal()


# 7. Outlier Detection by Group (L2 distance)

group_means <- df_x %>%
  group_by(time, skill_group) %>%
  summarise(mean_amp = mean(amplitude, na.rm = TRUE), .groups = "drop")
df_x_with_mean <- left_join(df_x, group_means, by = c("time", "skill_group"))
curve_distance <- df_x_with_mean %>%
  group_by(Player_ID, curve, skill_group, Hcp) %>%
  summarise(
    l2_dist = sqrt(mean((amplitude - mean_amp)^2, na.rm = TRUE)),
    .groups = "drop"
  )

ggplot(curve_distance, aes(x = Hcp, y = l2_dist, color = skill_group)) +
  geom_point(alpha = 0.7, size = 2) +
  labs(title = "Curve Distance from Group Mean vs Handicap", x = "Handicap", y = "L2 Distance") +
  theme_minimal()


# 8. Density of Functional Features by Group

ggplot(curve_features, aes(x = max_x, fill = skill_group)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Peak X Amplitude by Skill Group", x = "Peak X Amplitude", y = "Density") +
  theme_minimal()



library(readxl)
library(dplyr)
library(ggplot2)

# Read impact sheets and map to Player_ID
impact_file <- "Data/Impact All.xlsx"
impact_sheets <- excel_sheets(impact_file)

# This function will convert '4i' -> 'AMID0004', etc.
sheet_to_player_id <- function(sheet) {
  num <- gsub("i", "", sheet)
  sprintf("AMID%04d", as.integer(num))
}
impact_player_ids <- sapply(impact_sheets, sheet_to_player_id)

# Read each impact sheet, compute mean metrics for that subject
impact_means <- lapply(seq_along(impact_sheets), function(i) {
  df <- read_excel(impact_file, sheet = impact_sheets[i])
  data.frame(
    Player_ID = impact_player_ids[i],
    mean_centre_speed = mean(df$`face centre speed (m/s)`, na.rm = TRUE),
    mean_dist_centre = mean(df$`Dist from Centre (mm)`, na.rm = TRUE)
  )
}) %>% bind_rows()

subject_peak_x <- df_x %>%
  group_by(Player_ID, skill_group) %>%
  summarise(mean_peak_x = mean(tapply(amplitude, curve, max)), .groups = "drop")

plot_df <- left_join(subject_peak_x, impact_means, by = "Player_ID")
plot_df_clean <- plot_df %>% filter(!is.na(mean_peak_x), !is.na(mean_dist_centre))

ggplot(plot_df_clean, aes(x = mean_peak_x, y = mean_dist_centre, color = skill_group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Swing Feature vs. Impact Metric",
    x = "Mean Peak X Amplitude (Swing Feature)",
    y = "Mean Distance from Centre (Impact Metric)",
    color = "Skill Group"
  ) +
  theme_minimal()


library(readxl)
library(dplyr)
library(ggplot2)

impact_file <- "Data/Impact All.xlsx"
impact_sheets <- excel_sheets(impact_file)

sheet_to_player_id <- function(sheet) {
  num <- gsub("i", "", sheet)
  sprintf("AMID%04d", as.integer(num))
}
impact_player_ids <- sapply(impact_sheets, sheet_to_player_id)

impact_means <- lapply(seq_along(impact_sheets), function(i) {
  df <- read_excel(impact_file, sheet = impact_sheets[i])
  data.frame(
    Player_ID = impact_player_ids[i],
    mean_centre_speed = mean(df$`face centre speed (m/s)`, na.rm = TRUE),
    mean_dist_centre = mean(df$`Dist from Centre (mm)`, na.rm = TRUE)
  )
}) %>% bind_rows()


plot_df <- left_join(subject_peak_x, impact_means, by = "Player_ID")
plot_df_clean <- plot_df %>% filter(!is.na(mean_peak_x), !is.na(mean_centre_speed))

ggplot(plot_df_clean, aes(x = mean_peak_x, y = mean_centre_speed, color = skill_group)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "Swing Feature vs. Face Centre Speed",
    x = "Mean Peak X Amplitude (Swing Feature)",
    y = "Mean Face Centre Speed (m/s)",
    color = "Skill Group"
  ) +
  theme_minimal()
