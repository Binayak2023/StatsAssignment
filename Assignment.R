# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   https://r-pkgs.org
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shif
#install.packages("tidyverse")
# Install necessary packages
#install.packages("ggplot2")
#install.packages("gridExtra")

# Load the necessary libraries
# Load the necessary libraries
library(ggplot2)
library(gridExtra)

# Load the dataset
data <- read.csv("C:\\Users\\user\\Downloads\\ds.csv")

# Inspect the column names
colnames(data)

# Rename columns for simplicity
colnames(data) <- c("Time", "x1", "x2", "x3", "x4", "x5")

# Check the first few rows of the dataset to ensure columns are renamed
head(data)

#Task 1
# Time series plots
p1 <- ggplot(data, aes(x = Time, y = x1)) + geom_line() + ggtitle("Gene x1 Time Series")
p2 <- ggplot(data, aes(x = Time, y = x2)) + geom_line() + ggtitle("Gene x2 Time Series")
p3 <- ggplot(data, aes(x = Time, y = x3)) + geom_line() + ggtitle("Gene x3 Time Series")
p4 <- ggplot(data, aes(x = Time, y = x4)) + geom_line() + ggtitle("Gene x4 Time Series")
p5 <- ggplot(data, aes(x = Time, y = x5)) + geom_line() + ggtitle("Gene x5 Time Series")

# Arrange plots in a grid
grid.arrange(p1, p2, p3, p4, p5, ncol = 2)

# Distribution plots
p1_dist <- ggplot(data, aes(x = x1)) + geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) + ggtitle("Distribution of Gene x1")
p2_dist <- ggplot(data, aes(x = x2)) + geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) + ggtitle("Distribution of Gene x2")
p3_dist <- ggplot(data, aes(x = x3)) + geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) + ggtitle("Distribution of Gene x3")
p4_dist <- ggplot(data, aes(x = x4)) + geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) + ggtitle("Distribution of Gene x4")
p5_dist <- ggplot(data, aes(x = x5)) + geom_histogram(binwidth = 0.1, fill = "blue", alpha = 0.7) + ggtitle("Distribution of Gene x5")

# Arrange distribution plots in a grid
grid.arrange(p1_dist, p2_dist, p3_dist, p4_dist, p5_dist, ncol = 2)

# Scatter plots and correlation coefficients
pairs(data[2:6], main = "Scatterplot Matrix")
correlation_matrix <- cor(data[2:6])
print(correlation_matrix)

# Split the data into training and testing sets (70% training, 30% testing)
set.seed(123)
train_indices <- sample(seq_len(nrow(data)), size = 0.7 * nrow(data))
train_data <- data[train_indices, ]
test_data <- data[-train_indices, ]

#  Task 2.1: Define the models
model1 <- lm(x2 ~ x4 + I(x3^2), data = train_data)
model2 <- lm(x2 ~ x4 + I(x3^2) + x5, data = train_data)
model3 <- lm(x2 ~ x3 + x4 + I(x5^3), data = train_data)
model4 <- lm(x2 ~ x4 + I(x3^2) + I(x5^3), data = train_data)
model5 <- lm(x2 ~ x4 + I(x1^2) + I(x3^2), data = train_data)

# Summarize the models to see the estimated parameters
summary(model1)
summary(model2)
summary(model3)
summary(model4)
summary(model5)

# Task 2.2: Compute RSS
rss1 <- sum(residuals(model1)^2)
rss2 <- sum(residuals(model2)^2)
rss3 <- sum(residuals(model3)^2)
rss4 <- sum(residuals(model4)^2)
rss5 <- sum(residuals(model5)^2)

# Task 2.3: Compute Log-Likelihood
n <- nrow(train_data)
sigma2_1 <- rss1 / (n - 1)
sigma2_2 <- rss2 / (n - 1)
sigma2_3 <- rss3 / (n - 1)
sigma2_4 <- rss4 / (n - 1)
sigma2_5 <- rss5 / (n - 1)

log_likelihood1 <- -n / 2 * log(2 * pi) - n / 2 * log(sigma2_1) - 1 / (2 * sigma2_1) * rss1
log_likelihood2 <- -n / 2 * log(2 * pi) - n / 2 * log(sigma2_2) - 1 / (2 * sigma2_2) * rss2
log_likelihood3 <- -n / 2 * log(2 * pi) - n / 2 * log(sigma2_3) - 1 / (2 * sigma2_3) * rss3
log_likelihood4 <- -n / 2 * log(2 * pi) - n / 2 * log(sigma2_4) - 1 / (2 * sigma2_4) * rss4
log_likelihood5 <- -n / 2 * log(2 * pi) - n / 2 * log(sigma2_5) - 1 / (2 * sigma2_5) * rss5

# Task 2.4: Compute AIC and BIC
k1 <- length(coef(model1))
k2 <- length(coef(model2))
k3 <- length(coef(model3))
k4 <- length(coef(model4))
k5 <- length(coef(model5))

aic1 <- 2 * k1 - 2 * log_likelihood1
aic2 <- 2 * k2 - 2 * log_likelihood2
aic3 <- 2 * k3 - 2 * log_likelihood3
aic4 <- 2 * k4 - 2 * log_likelihood4
aic5 <- 2 * k5 - 2 * log_likelihood5

bic1 <- k1 * log(n) - 2 * log_likelihood1
bic2 <- k2 * log(n) - 2 * log_likelihood2
bic3 <- k3 * log(n) - 2 * log_likelihood3
bic4 <- k4 * log(n) - 2 * log_likelihood4
bic5 <- k5 * log(n) - 2 * log_likelihood5

# Print RSS, log-likelihood, AIC, and BIC values
rss_values <- c(rss1, rss2, rss3, rss4, rss5)
log_likelihood_values <- c(log_likelihood1, log_likelihood2, log_likelihood3, log_likelihood4, log_likelihood5)
aic_values <- c(aic1, aic2, aic3, aic4, aic5)
bic_values <- c(bic1, bic2, bic3, bic4, bic5)

results <- data.frame(
  Model = paste("Model", 1:5),
  RSS = rss_values,
  Log_Likelihood = log_likelihood_values,
  AIC = aic_values,
  BIC = bic_values
)

print(results)

# Adjust the plotting parameters
par(mfrow = c(5, 2), mar = c(3, 3, 2, 1), cex.main = 1.5, cex.lab = 1.2, cex.axis = 1)  # Adjust the margins as necessary

# Plot residuals and Q-Q plots for each model (Task 2.5)

# Model 1
plot(model1, which = 1, main = "Model 1 Residuals")
qqnorm(residuals(model1), main = "Model 1 Q-Q Plot")
qqline(residuals(model1))

# Model 2
plot(model2, which = 1, main = "Model 2 Residuals")
qqnorm(residuals(model2), main = "Model 2 Q-Q Plot")
qqline(residuals(model2))

# Model 3
plot(model3, which = 1, main = "Model 3 Residuals")
qqnorm(residuals(model3), main = "Model 3 Q-Q Plot")
qqline(residuals(model3))

# Model 4
plot(model4, which = 1, main = "Model 4 Residuals")
qqnorm(residuals(model4), main = "Model 4 Q-Q Plot")
qqline(residuals(model4))

# Model 5
plot(model5, which = 1, main = "Model 5 Residuals")
qqnorm(residuals(model5), main = "Model 5 Q-Q Plot")
qqline(residuals(model5))

par(mfrow = c(1, 1))  # Reset the plotting parameters to default

#Task 2.7
# Assuming model5 is the best model based on previous tasks
best_model <- model5

# Predict on test data
predictions <- predict(best_model, newdata = test_data)

# Check if predictions object is created
head(predictions)

# Plot predictions vs actuals
ggplot() +
  geom_point(aes(x = test_data$x2, y = predictions)) +
  geom_abline(intercept = 0, slope = 1, col = "red") +
  ggtitle("Predictions vs Actuals") +
  xlab("Actual Values") +
  ylab("Predicted Values")

# Calculate 95% confidence intervals
conf_intervals <- predict(best_model, newdata = test_data, interval = "confidence", level = 0.95)

# Plot predictions with 95% confidence intervals
ggplot() +
  geom_point(aes(x = test_data$Time, y = test_data$x2), col = "blue") +
  geom_line(aes(x = test_data$Time, y = predictions), col = "red") +
  geom_ribbon(aes(x = test_data$Time, ymin = conf_intervals[, "lwr"], ymax = conf_intervals[, "upr"]), alpha = 0.2) +
  ggtitle("Predictions with 95% Confidence Intervals") +
  xlab("Time") +
  ylab("Gene x2")

# Extract and print the coefficients (theta) of Model 5
theta <- coef(model5)
print(theta)

#Task3
# Load necessary libraries
library(abc)
library(ggplot2)

# Define the observed data (response variable)
observed_data <- train_data$x2

# Define the number of simulations
num_simulations <- 10000

# Define prior distributions for the parameters using Uniform distributions
prior_intercept <- runif(num_simulations, min = 1.44619299 - 2, max = 1.44619299 + 2)
prior_x4 <- runif(num_simulations, min = -0.97499705 - 2, max = -0.97499705 + 2)

# Generate parameter samples
parameter_samples <- data.frame(intercept = prior_intercept, x4 = prior_x4)

# Simulate data using the generated parameter samples
simulated_data <- apply(parameter_samples, 1, function(params) {
  intercept <- params[1]
  x4 <- params[2]

  # Use fixed values for other parameters
  x1_squared <- 0.53906644 * train_data$x1^2
  x3_squared <- -0.06778866 * train_data$x3^2
  x4_values <- x4 * train_data$x4

  intercept + x4_values + x1_squared + x3_squared
})

# Calculate the distance metric (sum of squared errors)
distances <- apply(simulated_data, 1, function(sim) sum((sim - observed_data)^2))

# Set the tolerance level for acceptance
tolerance <- quantile(distances, 0.01)

# Accept samples with distances less than the tolerance
accepted_indices <- which(distances < tolerance)
accepted_samples <- parameter_samples[accepted_indices, ]

# Plot the posterior distributions
hist(accepted_samples$intercept, main = "Posterior Distribution of Intercept", xlab = "Intercept", breaks = 50, col = "blue")
hist(accepted_samples$x4, main = "Posterior Distribution of x4", xlab = "x4", breaks = 50, col = "blue")

# Summarize the accepted samples
summary(accepted_samples)

#Plot the joint and marginal posterior distributions

# Joint posterior distribution plot
joint_plot <- ggplot(accepted_samples, aes(x = intercept, y = x4)) +
  geom_point(alpha = 0.5) +
  labs(title = "Joint Posterior Distribution", x = "Intercept", y = "x4")

# Marginal posterior distributions
marginal_intercept <- ggplot(accepted_samples, aes(x = intercept)) +
  geom_histogram(fill = "blue", bins = 50, alpha = 0.7) +
  labs(title = "Marginal Posterior Distribution of Intercept", x = "Intercept", y = "Frequency")

marginal_x4 <- ggplot(accepted_samples, aes(x = x4)) +
  geom_histogram(fill = "blue", bins = 50, alpha = 0.7) +
  labs(title = "Marginal Posterior Distribution of x4", x = "x4", y = "Frequency")

# Arrange the plots
library(gridExtra)
grid.arrange(joint_plot, marginal_intercept, marginal_x4, ncol = 1)
