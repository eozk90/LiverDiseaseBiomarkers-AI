library(tidyr)
library(pROC)
library(formattable)
library(MASS)
library(bootStepAIC)
library(GGally)
library(verification)
library(ftExtra)
library(flextable)
library(openxlsx)
library(readxl)
library(knitr)
library(ciTools)
library(kableExtra)
library(survival)
library(survminer)
library(epitools)
library(dplyr)
library(broom)

excel_file <- "C:\\Users\\gulsa\\Downloads\\VCTE vs MRE FAST MRI April 9.xlsx"

# Read the Excel file
data <- read_excel(excel_file, sheet = "Sheet3")
head(data)
column_names <- colnames(data)
print(column_names)

perform_mann_whitney_roc <- function(data, group_column, parameter_indices, predictor_indices, response) {
  
  
  # Create an empty data frame to store results
  results_df <- data.frame(
    Parameter = character(),
    Test = character(),
    Mean_Group_0 = numeric(),
    Mean_Group_1 = numeric(),
    SD_Group_0 = numeric(),
    SD_Group_1 = numeric(),
    Count_Group_0 = numeric(),
    Count_Group_1 = numeric(),
    P_Value = numeric(),
    AUC = numeric(),
    CI_Lower = numeric(),
    CI_Upper = numeric(),
    Cutoff = numeric(),
    Sens = numeric(),
    Spec = numeric(),
    RuleIn = numeric(),
    RuleOut = numeric(),
    stringsAsFactors = FALSE
  )
  
  data$binary_response <- as.factor(ifelse(data[[response]] == 1, 1, 0))
  
  AUCs_for_DeLong_Test <- list()
  
  # Loop through each parameter index
  for (param_index in parameter_indices) {
    param_name <- colnames(data)[param_index]
    
    predictor <- colnames(data)[param_index]
    
    data_filtered <- na.omit(data[complete.cases(data[[predictor]], data[[response]]), 
                                  c(predictor, response, "binary_response")])
    
    # Perform ROC analysis
    ci_result <- ci.auc(data_filtered$binary_response,as.numeric(data_filtered[[predictor]]))
    roc_results <- roc(data_filtered$binary_response,as.numeric(data_filtered[[predictor]]))
    
    YoudenJindex = which.max(abs(as.array(unlist(roc_results$sensitivities))+as.array(unlist(roc_results$specificities))-1))
    Sensitivity = roc_results$sensitivities[YoudenJindex]
    Specificity = roc_results$specificities[YoudenJindex]
    Threshold = roc_results$thresholds[YoudenJindex]
    
    # Rule In cutoff (high specificity: few false positives)
    rule_in_index <- which.max(roc_results$specificities >= 0.95)
    rule_in_threshold <- roc_results$thresholds [rule_in_index]
    
    # Rule Out cutoff (high sensitivity: few false negatives)
    rule_out_index <- which.max(roc_results$sensitivities <= 0.95)
    rule_out_threshold <- roc_results$thresholds [rule_out_index]
    
    
    # Separate data into groups based on CPSH (0 and 1)
    group_0 <- data[data[[group_column]] == 0, ]
    group_1 <- data[data[[group_column]] == 1, ]
    
    # Calculate means and standard deviations for each group
    mean_group_0 <- mean(group_0[[param_name]], na.rm = TRUE)
    mean_group_1 <- mean(group_1[[param_name]], na.rm = TRUE)
    sd_group_0 <- sd(group_0[[param_name]], na.rm = TRUE)
    sd_group_1 <- sd(group_1[[param_name]], na.rm = TRUE)
    
    # Get data count for each group
    count_group_0 <- sum(!is.na(group_0[[param_name]]))
    count_group_1 <- sum(!is.na(group_1[[param_name]]))
    
    # Perform Mann-Whitney U test
    test_result <- wilcox.test(data[[param_name]] ~ data[[group_column]])
    
    # Store results in the data frame
    results_df <- rbind(results_df, data.frame(
      Parameter = param_name,
      Count_G_0 = count_group_0,
      Mean_G_0 = round(mean_group_0,3),
      SD_G_0 = round(sd_group_0,3),
      Count_G_1 = count_group_1,
      Mean_G_1 = round(mean_group_1,3),
      SD_G_1 = round(sd_group_1,3),
      P_Value = round(test_result$p.value,5),
      AUC = round(ci_result[2],digits=2),
      CI_L = round(ci_result[1],digits=2),
      CI_U = round(ci_result[3],digits=2),
      Cutoff = round(Threshold,3),
      Sens = round(Sensitivity,digits=3),
      Spec = round(Specificity,digits=3),
      RuleIn = round(rule_in_threshold,3),
      RuleOut = round(rule_out_threshold,3)
    ))
  }
  
  # Print the results table
  kable(results_df, format = "markdown")
}


### stats for prediction of liver decompensation
group_column <- "Decompensation" 
parameter_indices <- c(2,5:7,9,12,17,19,23:25,28:38,40,42,48,51,54:55)
response <- "Decompensation" 
predictor_indices <- c(2,5:7,9,12,17,19.23:25,28:38,40,42,48,51,54:55)
perform_mann_whitney_roc(data, group_column, parameter_indices, predictor_indices, response)

################################ combinations with logistic regression #####################################################
# Load data
excel_file <- "C:\\Users\\gulsa\\Downloads\\VCTE vs MRE FAST MRI April 9.xlsx"
data <- read_excel(excel_file, sheet = "Sheet3")

# Define combinations (each combination is a vector of column indices)
combinations_list <- list(
  c(60, 12, 9),   # VCTE + CAP
  c(60, 12, 48),  # VCTE + CP-Score
  c(60, 12, 51),  # VCTE + FIB-4
  c(60, 12, 54),  # VCTE + MELD 3.0
  c(60, 12, 55),  # VCTE + MELD 3.0-Na
  c(60, 19, 17),  # 2D MRE + MRI-PDFF
  c(60, 19, 48),  # 2D MRE + CP-Score
  c(60, 19, 51),  # 2D MRE + FIB-4
  c(60, 19, 54),  # 2D MRE + MELD 3.0
  c(60, 19, 55)   # 2D MRE + MELD 3.0-Na
)

# Iterate through each combination
for (columns_to_check in combinations_list) {
  data_filtered <- data[complete.cases(data[, columns_to_check]), ]
  
  # Construct formula as string
  formula_string <- paste0("data_filtered[[", columns_to_check[1], "]] ~ ",
                           paste0("data_filtered[[", columns_to_check[-1], "]]", collapse = " + "))
  
  # Fit logistic regression model
  model <- glm(formula_string, data = data_filtered, family = "binomial")
  
  # Stepwise selection
  mod_step <- stepAIC(model, direction = 'backward', trace = FALSE)
  
  # Predict probabilities
  predicted_probs <- predict(mod_step, newdata = data_filtered, type = "response")
  
  # ROC analysis
  roc_obj <- roc(data_filtered[[columns_to_check[1]]], predicted_probs)
  auc_value <- round(auc(roc_obj), 3)
  auc_ci <- ci(roc_obj)
  lower_ci <- round(auc_ci[1], 3)
  upper_ci <- round(auc_ci[3], 3)
  
  # Youden's J index
  YoudenJindex <- which.max(roc_obj$sensitivities + roc_obj$specificities - 1)
  Sensitivity <- round(roc_obj$sensitivities[YoudenJindex], 3)
  Specificity <- round(roc_obj$specificities[YoudenJindex], 3)
  
  # Print result
  cat("Predictors:", paste(colnames(data)[columns_to_check[-1]], collapse = " + "), "\n")
  cat("AUC:", auc_value, "CI:", lower_ci, "-", upper_ci,
      "Sensitivity:", Sensitivity, "Specificity:", Specificity, "\n\n")
}


#####################################################################  Cox proportional hazard analysis #########################################################################
excel_file <- "C:\\Users\\gulsa\\Downloads\\VCTE vs MRE FAST MRI April 9.xlsx"
data <- read_excel(excel_file, sheet = "Sheet3")
data$Gender <- ifelse(data$Gender == "Male", 1, 0)
column_names <- colnames(data)
# Define the parameter indices for univariate
#parameter_indices <- c(2,5:7,9,12,17,19,23:25,28:38,40,42,48,51,54:55)
# Define the parameter indices for multivariate
parameter_indices <- c(6:7,9,12,17,19,23:25,28:38,40,42,48,51,54:55)


# Create a list to store results
results_list <- list()

# Loop through each parameter
for (i in parameter_indices) {
  #For univariate analysis
  #subset_data <- data[, c(column_names[i], "Time To Follow Up or Decompensation", "Decompensation")]
  #survival_obj <- Surv(time = subset_data$"Time To Follow Up or Decompensation", event = subset_data$"Decompensation")
  #cox_model <- coxph(survival_obj ~ ., data=subset_data[, 1])
  #predicted_times <- predict(cox_model, newdata = subset_data, type = "expected")
  #results_list[[i]] <- summary(cox_model)
  #results_list[[i]]$conf.int
  
  #multivariate analysis adjusted for Age, Gender, BMI
  subset_data <- data[, c(column_names[i], "Age", "Gender", "BMI", "Time To Follow Up or Decompensation", "Decompensation")]
  survival_obj <- Surv(time = subset_data$"Time To Follow Up or Decompensation", event = subset_data$"Decompensation")
  ##Fit the Cox model using all other columns in data as predictors (excluding the survival object components)
  cox_model <- coxph(survival_obj ~ ., data=subset_data[, 1:4]) #for multivariate age, gender, bmi adjusted
  predicted_times <- predict(cox_model, newdata = subset_data, type = "expected")
  results_list[[i]] <- summary(cox_model)
  results_list[[i]]$conf.int
  
}


# Combine results into a table
results_table <- do.call(rbind, lapply(parameter_indices, function(i) {
  coef <- round(results_list[[i]]$coefficients[1,1], digits = 3) #coefficients[1] for univariate
  exp_coef <- round(results_list[[i]]$coefficients[1,2], digits = 3) #coefficients[2] for univariate
  se_coef <- round(results_list[[i]]$coefficients[1,3], digits = 3)  #coefficients[3] for univariate
  z_value <- round(results_list[[i]]$coefficients[1,4], digits = 3)  #coefficients[4] for univariate
  p_value <- round(results_list[[i]]$coefficients[1,5], digits = 3)  #coefficients[5] for univariate
  
  # Determine significance level
  if (p_value < 0.001) {
    significance <- "***"
  } else if (p_value < 0.01) {
    significance <- "**"
  } else if (p_value < 0.05) {
    significance <- "*"
  } else {
    significance <- ""
  }
  
  # Calculate confidence intervals
  ci_lower <- round(results_list[[i]]$conf.int[1,3], digits = 3)  #conf.int[3] for univariate
  ci_upper <- round(results_list[[i]]$conf.int[1,4], digits = 3) #conf.int[4] for univariate
  
  
  data.frame(
    "Parameter" = names(data)[i],
    "coef" = coef,
    "Hazard Ratio" = exp_coef,
    "CI Lower" = ci_lower,
    "CI Upper" = ci_upper,
    "SE Coef" = se_coef,
    "z-value" = z_value,
    "p-value" = p_value,
    "Significance" = significance
  )
}))

# Print the results table
print(results_table)