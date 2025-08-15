install.packages("future.apply") 
library(future.apply)
library(dplyr)
library(lme4)
library(future.apply)
library(readxl)

# 1. Prepare data
# ---- 1. Read player HCP data ----
raw <- read_excel("Data/Player data.xlsx") # <-- update with your actual file name

# Extract Group 1 and Group 2 as in your code
group1 <- raw[, 2:3]
colnames(group1) <- c("Subject_ID", "hcp")
group1 <- group1 %>% filter(grepl("^AMID", Subject_ID))

group2 <- raw[, 8:9]
colnames(group2) <- c("Subject_ID", "hcp")
group2 <- group2 %>% filter(grepl("^AMID", Subject_ID))

players <- bind_rows(group1, group2)
players$hcp <- as.numeric(players$hcp)

# ---- 2. Merge FPC means and HCP data ----
final_df <- left_join(subject_fpc_means, players, by = "Subject_ID")
colnames(players)
colnames(scores_df)
colnames(players)
# Check for missing values or mismatches
print(final_df)
scores_df <- data.frame(
  Subject_ID = subject_ids,
  mfpca_result$scores
)
colnames(scores_df)[-1] <- paste0("FPC", seq_len(ncol(mfpca_result$scores)))
scores_df$Curve_ID <- seq_len(nrow(scores_df))

df <- left_join(scores_df, players, by = "Subject_ID") %>%
  filter(!is.na(hcp), !if_any(starts_with("FPC"), is.na))

median_hcp <- median(df$hcp, na.rm = TRUE)
df$skill_group <- ifelse(df$hcp <= median_hcp, "High", "Low")
df$skill_group <- factor(df$skill_group, levels = c("Low", "High"))

# 2. Fit mixed-effects model (full data, for reference)
model <- glmer(skill_group ~ FPC1 + FPC2 + FPC3 + FPC4 + (1 | Subject_ID),
               data = df, family = binomial)
summary(model)

# Fit mixed-effects model using FPC1, FPC5, and FPC6
model_fpc1_5_6 <- glmer(skill_group ~ FPC1 + FPC5 + FPC6 + (1 | Subject_ID),
                        data = df, family = binomial)

# Display model summary
summary(model_fpc1_5_6)



# Odds ratios and CIs
or <- exp(fixef(model))
ci <- exp(confint(model, parm = "beta_", method = "Wald"))
odds_table <- cbind(Odds_Ratio = or, CI_lower = ci[,1], CI_upper = ci[,2])
print(round(odds_table, 3))

# In-sample predictions (not a valid measure of model performance, but shows model fit)
df$predicted_prob <- predict(model, type = "response")
df$predicted_class <- ifelse(df$predicted_prob > 0.5, "High", "Low")
df$predicted_class <- factor(df$predicted_class, levels = c("Low", "High"))
cm <- table(Predicted = df$predicted_class, Actual = df$skill_group)
print(cm)
accuracy <- sum(diag(cm)) / sum(cm)
cat("In-sample classification accuracy:", round(100 * accuracy, 2), "%\n")

# 3. Parallel LOSO Cross-Validation
plan(multisession, workers = parallel::detectCores() - 1) # Use all but 1 core

subjects <- unique(df$Subject_ID)
cat("Running LOSO CV on", length(subjects), "subjects in parallel...\n")

cv_list <- future_lapply(subjects, function(sid) {
  train_df <- df %>% filter(Subject_ID != sid)
  test_df  <- df %>% filter(Subject_ID == sid)
  if(n_distinct(train_df$skill_group) < 2) return(NULL)
  model_cv <- glmer(skill_group ~ FPC1 + FPC2 + FPC3 + FPC4 + (1 | Subject_ID),
                    data = train_df, family = binomial,
                    control = glmerControl(optimizer = "bobyqa")) # Slightly faster optimizer
  test_df$pred_prob <- predict(model_cv, newdata = test_df, type = "response", allow.new.levels = TRUE)
  test_df$pred_class <- ifelse(test_df$pred_prob > 0.5, "High", "Low")
  test_df$pred_class <- factor(test_df$pred_class, levels = c("Low", "High"))
  test_df
})

cv_results <- bind_rows(cv_list)

cm_cv <- table(Predicted = cv_results$pred_class, Actual = cv_results$skill_group)
print(cm_cv)
cv_accuracy <- sum(diag(cm_cv)) / sum(cm_cv)
cat("Parallel LOSO CV accuracy:", round(100 * cv_accuracy, 2), "%\n")


# Calculate sensitivity (recall), specificity, balanced accuracy
TP <- cm_cv["High", "High"]
TN <- cm_cv["Low", "Low"]
FP <- cm_cv["High", "Low"]
FN <- cm_cv["Low", "High"]

sensitivity <- TP / (TP + FN)  # True Positive Rate (High skill correctly classified)
specificity <- TN / (TN + FP)  # True Negative Rate (Low skill correctly classified)
balanced_acc <- (sensitivity + specificity) / 2

cat("Sensitivity:", round(100 * sensitivity, 2), "%\n")
cat("Specificity:", round(100 * specificity, 2), "%\n")
cat("Balanced Accuracy:", round(100 * balanced_acc, 2), "%\n")

install.packages("pROC")
library(pROC)
roc_obj <- roc(cv_results$skill_group, cv_results$pred_prob, levels = c("Low", "High"))
plot(roc_obj, main = "ROC Curve - LOSO CV")
auc(roc_obj)
 

library(ggplot2)
ggplot(cv_results, aes(x = predicted_prob, fill = skill_group)) +
  geom_histogram(position = "identity", alpha = 0.5, bins = 30) +
  labs(title = "Distribution of Predicted Probabilities by True Skill Group",
       x = "Predicted Probability of High Skill",
       y = "Count") +
  theme_minimal()

subject_acc <- cv_results %>%
  group_by(Subject_ID) %>%
  summarise(
    true_group = unique(skill_group),
    n_curves = n(),
    accuracy = mean(pred_class == skill_group)
  )
print(subject_acc)
# Plot accuracy per subject
ggplot(subject_acc, aes(x = reorder(Subject_ID, accuracy), y = accuracy, fill = true_group)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "Per-Subject Accuracy (LOSO CV)",
       x = "Subject",
       y = "Accuracy") +
  theme_minimal()



library(tidyr)
cv_results_long <- cv_results %>%
  pivot_longer(cols = starts_with("FPC"), names_to = "FPC", values_to = "Score")

ggplot(cv_results_long, aes(x = skill_group, y = Score, fill = skill_group)) +
  geom_boxplot() +
  facet_wrap(~ FPC, scales = "free_y") +
  labs(title = "Distribution of FPC Scores by Skill Group",
       x = "Skill Group", y = "FPC Score") +
  theme_minimal()

library(dplyr)
library(ggplot2)

# Bin the predicted probabilities
cv_results$prob_bin <- cut(cv_results$predicted_prob, breaks = seq(0, 1, 0.1), include.lowest = TRUE)

# For each bin, calculate the observed fraction of 'High'
calib_df <- cv_results %>%
  group_by(prob_bin) %>%
  summarise(
    mean_pred = mean(predicted_prob),
    frac_high = mean(skill_group == "High"),
    n = n()
  ) %>%
  filter(n > 0)  # Only plot bins with data

# Plot calibration
ggplot(calib_df, aes(x = mean_pred, y = frac_high)) +
  geom_line() +
  geom_point(size = 2) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
  labs(x = "Mean Predicted Probability (binned)", y = "Observed Fraction of High Skill",
       title = "Calibration Plot (LOSO CV)") +
  theme_minimal()


ggplot(cv_results, aes(x = skill_group, y = predicted_prob, fill = skill_group)) +
  geom_boxplot() +
  labs(title = "Model Predicted Probability by Skill Group",
       x = "True Skill Group", y = "Predicted Probability (High Skill)")






library(lme4)
library(caret)
library(pROC)

# 1. TRAIN-TEST SPLIT BY PLAYERS

# Get unique players and split them
unique_players <- unique(df$Subject_ID)
set.seed(123)  # For reproducibility

# 70% players for training, 30% for testing
n_train_players <- round(0.7 * length(unique_players))
train_players <- sample(unique_players, n_train_players)
test_players <- setdiff(unique_players, train_players)

# Split data
train_data <- df[df$Subject_ID %in% train_players, ]
test_data <- df[df$Subject_ID %in% test_players, ]

cat("Training players:", length(train_players), "\n")
cat("Testing players:", length(test_players), "\n")
cat("Training observations:", nrow(train_data), "\n")
cat("Testing observations:", nrow(test_data), "\n")


# 2. TRAIN MODEL

# Fit model on training data
model <- glmer(skill_group ~ FPC1 + FPC2 + FPC3 + FPC4 + (1 | Subject_ID),
               data = train_data, family = binomial)

# Fit model on training data using FPC1, FPC5, and FPC6
model_fpc1_5_6 <- glmer(skill_group ~ FPC1 + FPC5 + FPC6 + (1 | Subject_ID),
                        data = train_data, family = binomial)
print("Model Summary:")
summary(model)
scale(apply(train_data[, c("FPC1","FPC2","FPC3","FPC4")], 2, sd))


# 3. MAKE PREDICTIONS ON TEST DATA


# Predict probabilities (allow new levels for test players)
test_probs <- predict(model, newdata = test_data, type = "response", 
                      allow.new.levels = TRUE)

# Convert to classifications (threshold = 0.5)
test_preds <- ifelse(test_probs > 0.5, "High", "Low")
test_preds <- factor(test_preds, levels = c("Low", "High"))
test_actual <- factor(test_data$skill_group, levels = c("Low", "High"))


# 4. EVALUATE PERFORMANCE


# Confusion Matrix
conf_matrix <- confusionMatrix(test_preds, test_actual)
print(conf_matrix)

# Extract key metrics
accuracy <- conf_matrix$overall['Accuracy']
sensitivity <- conf_matrix$byClass['Sensitivity']
specificity <- conf_matrix$byClass['Specificity']
precision <- conf_matrix$byClass['Pos Pred Value']

# ROC and AUC
roc_obj <- roc(test_actual, test_probs)
auc_value <- auc(roc_obj)

# Print results
cat("\n=== PERFORMANCE RESULTS ===\n")
cat("Accuracy:", round(accuracy, 3), "\n")
cat("Sensitivity:", round(sensitivity, 3), "\n")
cat("Specificity:", round(specificity, 3), "\n")
cat("Precision:", round(precision, 3), "\n")
cat("AUC:", round(auc_value, 3), "\n")

# Plot ROC curve
plot(roc_obj, main = paste("ROC Curve (AUC =", round(auc_value, 3), ")"))


# 5. FEATURE IMPORTANCE


# Show which FPCs are significant
model_summary <- summary(model)
coefficients <- model_summary$coefficients

cat("\n=== SIGNIFICANT PREDICTORS ===\n")
for(i in 2:nrow(coefficients)) {  # Skip intercept
  predictor <- rownames(coefficients)[i]
  coef_val <- coefficients[i, "Estimate"]
  p_val <- coefficients[i, "Pr(>|z|)"]
  significance <- ifelse(p_val < 0.001, "***", 
                         ifelse(p_val < 0.01, "**", 
                                ifelse(p_val < 0.05, "*", "")))
  
  cat(predictor, ": coef =", round(coef_val, 3), 
      ", p =", round(p_val, 3), significance, "\n")
}

library(ggplot2)

# Convert confusion matrix table to data frame
cm_df <- as.data.frame(conf_matrix$table)
colnames(cm_df) <- c("Prediction", "Reference", "Freq")

# Create heatmap plot
ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", size = 6) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Confusion Matrix",
       x = "Actual",
       y = "Predicted") +
  theme_minimal() +
  theme(text = element_text(size = 14))


library(lme4)
library(caret)
library(pROC)
library(ggplot2)


# 1. TRAIN-TEST SPLIT BY PLAYERS

# Get unique players and split them
unique_players <- unique(df$Subject_ID)
set.seed(123)  # For reproducibility

# 70% players for training, 30% for testing
n_train_players <- round(0.7 * length(unique_players))
train_players <- sample(unique_players, n_train_players)
test_players <- setdiff(unique_players, train_players)

# Split data
train_data <- df[df$Subject_ID %in% train_players, ]
test_data <- df[df$Subject_ID %in% test_players, ]

cat("Training players:", length(train_players), "\n")
cat("Testing players:", length(test_players), "\n")
cat("Training observations:", nrow(train_data), "\n")
cat("Testing observations:", nrow(test_data), "\n")


# 2. TRAIN MODEL USING FPC1, FPC5, FPC6

# Fit model on training data using selected FPCs
model_fpc1_5_6 <- glmer(skill_group ~ FPC1 + FPC5 + FPC6 + (1 | Subject_ID),
                        data = train_data, family = binomial)

# Print model summary
cat("\n=== MODEL SUMMARY (FPC1, FPC5, FPC6) ===\n")
summary(model_fpc1_5_6)

# Standard deviation of predictors (optional scaling info)
scale(apply(train_data[, c("FPC1", "FPC5", "FPC6")], 2, sd))

# 3. MAKE PREDICTIONS ON TEST DATA


# Predict probabilities 
test_probs <- predict(model_fpc1_5_6, newdata = test_data, type = "response",
                      allow.new.levels = TRUE)

# Convert probabilities to predicted classes
test_preds <- ifelse(test_probs > 0.5, "High", "Low")
test_preds <- factor(test_preds, levels = c("Low", "High"))
test_actual <- factor(test_data$skill_group, levels = c("Low", "High"))

# 4. EVALUATE PERFORMANCE


# Confusion matrix
conf_matrix <- confusionMatrix(test_preds, test_actual)
print(conf_matrix)

# Extract performance metrics
accuracy <- conf_matrix$overall['Accuracy']
sensitivity <- conf_matrix$byClass['Sensitivity']
specificity <- conf_matrix$byClass['Specificity']
precision <- conf_matrix$byClass['Pos Pred Value']

# ROC and AUC
roc_obj <- roc(test_actual, test_probs)
auc_value <- auc(roc_obj)

# Print performance
cat("\n=== PERFORMANCE RESULTS ===\n")
cat("Accuracy:", round(accuracy, 3), "\n")
cat("Sensitivity:", round(sensitivity, 3), "\n")
cat("Specificity:", round(specificity, 3), "\n")
cat("Precision:", round(precision, 3), "\n")
cat("AUC:", round(auc_value, 3), "\n")

# Plot ROC curve
plot(roc_obj, main = paste("ROC Curve (AUC =", round(auc_value, 3), ")"))

# 5. FEATURE IMPORTANCE / SIGNIFICANCE

# Model summary coefficients
model_summary <- summary(model_fpc1_5_6)
coefficients <- model_summary$coefficients

cat("\n=== SIGNIFICANT PREDICTORS ===\n")
for(i in 2:nrow(coefficients)) {  # Skip intercept
  predictor <- rownames(coefficients)[i]
  coef_val <- coefficients[i, "Estimate"]
  p_val <- coefficients[i, "Pr(>|z|)"]
  significance <- ifelse(p_val < 0.001, "***", 
                         ifelse(p_val < 0.01, "**", 
                                ifelse(p_val < 0.05, "*", "")))
  
  cat(predictor, ": coef =", round(coef_val, 3), 
      ", p =", round(p_val, 3), significance, "\n")
}


# 6. VISUALIZE CONFUSION MATRIX

# Convert confusion matrix to data frame for plotting
cm_df <- as.data.frame(conf_matrix$table)
colnames(cm_df) <- c("Prediction", "Reference", "Freq")

# Plot heatmap
ggplot(cm_df, aes(x = Reference, y = Prediction, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), color = "black", size = 6) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "Confusion Matrix",
       x = "Actual",
       y = "Predicted") +
  theme_minimal() +
  theme(text = element_text(size = 14))


