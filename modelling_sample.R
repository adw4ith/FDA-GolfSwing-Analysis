library(readxl)
library(dplyr)


# Group 1 data appears in columns 2 (Player ID) and 3 (Hcp)
group1 <- raw[, 2:3]
colnames(group1) <- c("Subject", "hcp")
# Remove rows that aren't valid player IDs
group1 <- group1 %>%
  filter(grepl("^AMID", Subject))


# Group 2 data appears in columns 8 (Player ID) and 9 (Hcp)
group2 <- raw[, 8:9]
colnames(group2) <- c("Subject", "hcp")
group2 <- group2 %>%
  filter(grepl("^AMID", Subject))

#  Combine both groups 
players <- bind_rows(group1, group2)
players$hcp <- as.numeric(players$hcp)

# Final player list
print(players)

# Extract FPC scores 
scores_array <- mfpca_res$scores  # [n_curves, n_harmonics, 3]
scores_df <- data.frame(
  FPC1_X = scores_array[, 1, 1],
  FPC1_Y = scores_array[, 1, 2],
  FPC1_Z = scores_array[, 1, 3],
  FPC2_X = scores_array[, 2, 1],
  FPC2_Y = scores_array[, 2, 2],
  FPC2_Z = scores_array[, 2, 3]
)

# Map each swing to a Subject ID
subjects <- c(4,5,7,9,10:27,29:49,52:57,59:63)
subjects <- paste0("00", subjects)
subjects[1:4] <- paste0("0", subjects[1:4])
subjects <- paste0("AMID", subjects)

subject_ids <- rep(subjects, each = 40)
subject_ids <- subject_ids[-91]  
scores_df$Subject <- subject_ids

final_df <- left_join(scores_df, players, by = "Subject")
# Check the merged data
str(final_df)
head(final_df)

#linear regression
model <- lm(hcp ~ FPC1_X + FPC1_Y + FPC1_Z + FPC2_X + FPC2_Y + FPC2_Z, data = final_df)
summary(model)

# Create skill group (median split)
median_hcp <- median(final_df$hcp, na.rm = TRUE)
final_df$skill_group <- ifelse(final_df$hcp <= median_hcp, "High", "Low")
final_df$skill_group <- factor(final_df$skill_group, levels = c("Low", "High"))

# Logistic regression
model_bin <- glm(skill_group ~ FPC1_X + FPC1_Y + FPC1_Z + FPC2_X + FPC2_Y + FPC2_Z,
                 data = final_df, family = "binomial")
summary(model_bin)

boxplot(FPC1_X ~ skill_group, data = final_df, main = "FPC1_X by Skill Group")
boxplot(FPC2_X ~ skill_group, data = final_df, main = "FPC2_X by Skill Group")

predicted <- ifelse(predict(model_bin, type = "response") > 0.5, "High", "Low")
table(Predicted = predicted, Actual = final_df$skill_group)

accuracy <- (1013 + 274) / sum(c(274, 1013, 646, 226))
print(accuracy)
cm <- table(Predicted = predicted, Actual = final_df$skill_group)
accuracy <- sum(diag(cm)) / sum(cm)
print(accuracy)

final_df$skill_group <- factor(final_df$skill_group, levels = c("Low", "High"))
predicted <- factor(predicted, levels = c("Low", "High"))

cm <- table(Predicted = predicted, Actual = final_df$skill_group)
accuracy <- sum(diag(cm)) / sum(cm)
print(accuracy)

