# Load required libraries
library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
library(caret)
library(pROC)
library(gridExtra)

# Load dataset
colon <- read.csv("https://vincentarelbundock.github.io/Rdatasets/csv/survival/colon.csv")
colon <- na.omit(colon)  # Remove missing values

# Convert categorical variables
colon$rx <- factor(colon$rx, levels = c(0, 1, 2), labels = c("Obs", "Lev", "Lev+5FU"))
colon$status <- factor(colon$status, levels = c(0, 1), labels = c("Alive", "Dead"))

# Split the dataset into training and testing sets
set.seed(42)
train_index <- sample(1:nrow(colon), size = 0.7 * nrow(colon))
train <- colon[train_index, ]
test <- colon[-train_index, ]

### 1. Kaplan-Meier Analysis -----------------------------------------------------

# Survival curve by treatment (rx)
km_fit1 <- survfit(Surv(time, status == "Dead") ~ rx, data = colon)
ggsurvplot(km_fit1, data = colon, risk.table = TRUE, pval = TRUE,
           title = "Kaplan-Meier Curve for Treatment Groups", xlab = "Time (days)", ylab = "Survival Probability")

# Survival curve by node4
km_fit2 <- survfit(Surv(time, status == "Dead") ~ node4, data = colon)
ggsurvplot(km_fit2, data = colon, risk.table = TRUE, pval = TRUE,
           title = "Kaplan-Meier Curve for Node4 Groups", xlab = "Time (days)", ylab = "Survival Probability")

### 2. Cox Proportional Hazards Model --------------------------------------------

# Cox regression for treatment groups
cox_fit1 <- coxph(Surv(time, status == "Dead") ~ rx, data = colon)
summary(cox_fit1)

# Cox regression for node4
cox_fit2 <- coxph(Surv(time, status == "Dead") ~ node4, data = colon)
summary(cox_fit2)

# Cox regression with multiple predictors
cox_fit3 <- coxph(Surv(time, status == "Dead") ~ age + obstruct + perfor + nodes, data = colon)
summary(cox_fit3)

# Concordance Index for Cox models
c_index1 <- concordance.index(predict(cox_fit1), Surv(colon$time, colon$status == "Dead"))
print(paste("C-index for treatment model:", c_index1$c.index))

c_index2 <- concordance.index(predict(cox_fit2), Surv(colon$time, colon$status == "Dead"))
print(paste("C-index for node4 model:", c_index2$c.index))

# Schoenfeld test to check proportional hazard assumption
schoenfeld_test1 <- cox.zph(cox_fit1)
print(schoenfeld_test1)

schoenfeld_test2 <- cox.zph(cox_fit2)
print(schoenfeld_test2)

### 3. Logistic Regression -------------------------------------------------------

# Create binary outcome for logistic regression
colon$y <- ifelse(colon$time > median(colon$time), 1, 0)

# Fit logistic regression model
logit_model <- glm(y ~ rx + nodes, data = colon, family = binomial)
summary(logit_model)

# Predictions and ROC Curve
predictions <- predict(logit_model, type = "response")
roc_curve <- roc(colon$y, predictions)
plot(roc_curve, main = "ROC Curve for Logistic Regression Model")

# Model accuracy
accuracy <- sum((predictions > 0.5) == colon$y) / nrow(colon)
print(paste("Accuracy:", accuracy))

### 4. Model Diagnostics: QQ Plot, Residuals, and Scale Plot ----------------------

# QQ Plot
qqnorm(residuals(logit_model), main = "QQ Plot of Residuals")
qqline(residuals(logit_model))

# Residual Plot
plot(fitted(logit_model), residuals(logit_model), main = "Residual Plot",
     xlab = "Fitted Values", ylab = "Residuals", col = "blue", pch = 20)
abline(h = 0, col = "red")

# Scale-Location Plot (Variance Check)
plot(fitted(logit_model), sqrt(abs(residuals(logit_model))),
     main = "Scale-Location Plot", xlab = "Fitted Values", ylab = "sqrt(|Residuals|)", col = "blue", pch = 20)
abline(h = 0, col = "red")

### 5. Statistical Tests: T-test, F-test ------------------------------------------

# T-test to compare survival time between treatment groups
t_test_result <- t.test(time ~ rx, data = colon)
print(t_test_result)

# F-test to check variance homogeneity between treatment groups
f_test_result <- var.test(time ~ rx, data = colon)
print(f_test_result)

### 6. Visualizations: Boxplots and Bar Charts -----------------------------------

# Boxplot of survival time by treatment group
ggplot(colon, aes(x = rx, y = time, fill = rx)) +
  geom_boxplot() +
  labs(title = "Survival Time by Treatment Group", x = "Treatment", y = "Survival Time (days)") +
  theme_minimal()

# Bar chart of extent by treatment
ggplot(colon, aes(x = factor(extent), fill = rx)) +
  geom_bar(position = "dodge") +
  labs(title = "Extent of Disease by Treatment", x = "Extent", y = "Count") +
  theme_minimal()

# Scatter plot: Time vs Nodes by treatment
ggplot(colon, aes(x = nodes, y = time, color = rx)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Survival Time vs. Nodes by Treatment", x = "Number of Nodes", y = "Survival Time (days)") +
  theme_minimal()


