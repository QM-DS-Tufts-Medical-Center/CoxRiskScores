library(glmnet)
library(survival)
library(survminer)
library(ggplot2)
library(tidyr)
library(dplyr)


#' Fit the survival models for the phenotypes of interest
#'
#' This function constructs survival objects (`Surv`), where for each object 'time' is time spent with a condition/disease (referred as 'phenotype') and 'event' is the status of the phenotype.
#'
#' @param data A data frame containing survival time (follow-up time), and predictor variables.
#' @param pheno_cols A character vector specifying the names of phenotype columns with time (follow-up).
#'
#' @return A list containing:
#' - `surv_models`: A named list where each element is a `Surv` object for a corresponding phenotype.
#' @import survival
#' @export
create_surv_models <- function(data, pheno_names) {
  
  # Avoid modifying the original
  data_copy <- data  
  
  surv_models <- lapply(pheno_names, function(col_name) {
    event <- as.integer(data_copy[[paste0("Health_", col_name)]] == "Yes")
    
    time <- ifelse(
      data_copy[[paste0("Health_", col_name)]] == "Yes",
      data_copy[[paste0("Health_", col_name, "_Age")]],
      data_copy$Age_last
    )
    
    Surv(time = time, event = event)
  })
  
  names(surv_models) <- pheno_names
  return(surv_models = surv_models)
}


#' Optimize Cox Elastic Net Model
#'
#' This function performs cross-validated Cox proportional hazards modeling 
#' with elastic net regularization over a range of alpha values to find 
#' the best model for each survival response variable.
#'
#' @param survModels A named list of survival response variables, where each element 
#'        is a `Surv` object representing time-to-event data.
#' @param x A numeric matrix of predictor variables (omics features). Ensure that `x` is 
#'        preprocessed (e.g., scaled and log-transformed if necessary).
#' @param alpha_values A numeric vector of alpha values (default: `seq(0, 1, by = 0.1)`) 
#'        to search over. Alpha controls the balance between L1 (lasso) and L2 (ridge) penalties.
#'
#' @return A list of results for each survival model, including:
#' \describe{
#'   \item{`best_alpha`}{The optimal alpha value found via cross-validation.}
#'   \item{`best_lambda`}{The corresponding optimal lambda (regularization parameter).}
#'   \item{`min_cvm`}{The minimum cross-validated error achieved.}
#'   \item{`best_model`}{The fitted `glmnet` model with the best alpha and lambda.}
#'   \item{`pred`}{Predicted risk scores for the training data.}
#'   \item{`coefficients`}{A matrix of model coefficients at the best lambda.}
#' }
#'
#' @import glmnet survival
#' @export
#' @examples
#' # Example usage:
#' library(survival)
#' library(glmnet)
#'
#' # Simulated survival data
#' set.seed(123)
#' time <- rexp(100, rate = 0.1)
#' status <- sample(0:1, 100, replace = TRUE)
#' survModels <- list(outcome = Surv(time, status))
#' x <- matrix(rnorm(100 * 10), ncol = 10)
#'
#' # Run optimization
#' results <- optimize_cox_elastic_net(survModels, x)
#'
optimize_cox_elastic_net <- function(survModels, x, alpha_values = seq(0, 1, by = 0.1)) {
  library(glmnet)
  library(survival)
  
  results <- list()  # Initialize an empty list to store results
  
  for (col_name in names(survModels)) {
    cat("Processing", col_name, "\n")
    
    # Extract the current survival response variable (Surv object)
    y <- survModels[[col_name]]
    
    # Initialize variables to track the best model selection
    best_alpha <- NULL
    best_lambda <- NULL
    min_cvm <- Inf  # Set initial minimum cross-validation error to infinity
    best_model <- NULL
    
    # Loop through each alpha value to find the best one
    for (alpha in alpha_values) {
      cat("Testing alpha =", alpha, "\n")
      
      # Perform cross-validation with elastic net Cox regression
      cv_fit <- cv.glmnet(x, y, family = "cox", alpha = alpha)
      
      # Extract the best lambda and corresponding cross-validation error
      lambda_min <- cv_fit$lambda.min
      cvm_min <- min(cv_fit$cvm)
      
      # Update best model parameters if this is the lowest cvm found so far
      if (cvm_min < min_cvm) {
        min_cvm <- cvm_min
        best_alpha <- alpha
        best_lambda <- lambda_min
        best_model <- glmnet(x, y, family = "cox", alpha = alpha)
      }
    }
    
    # Compute prediction scores using the best model
    pred <- predict(best_model, newx = x, s = best_lambda, type = "link")
    
    # Extract coefficients as a named vector
    coefficients <- as.matrix(coef(best_model, s = best_lambda))
    
    # Store results for the current survival model
    results[[col_name]] <- list(
      best_alpha = best_alpha,
      best_lambda = best_lambda,
      min_cvm = min_cvm,
      best_model = best_model,
      pred = pred,
      coefficients = coefficients
    )
    
    cat("Best alpha for", col_name, ":", best_alpha, "\n")
    cat("Best lambda for", col_name, ":", best_lambda, "\n\n")
  }
  
  return(results) 
}



#' Plot Hazard Risk Scores
#'
#' This function plot computed risk scores
#'
#' @param RS_col A vector containing the Risk scores
#' @param col_name A name of the phenotype for which the Risk score is computed
#'
#'
#' @examples
#' \dontrun{
#'   visualize_HRS_df(HRS_df, meta_data = NULL, plot_type = "density")
#' }
#' 
#' @import Risk scores, 
#' @export 


visualize_HRS_df <- function(HRS_df, meta_data = NULL, plot_type = "density") {
  # Convert list to data frame
  df <- as.data.frame(HRS_df)
  df$CellID <- 1:nrow(df)  # optional index column
  
  # Reshape to long format
  df_long <- df %>%
    pivot_longer(cols = starts_with("HRS_"), names_to = "Variable", values_to = "HRS")
  
  if (!is.null(meta_data)) {
  if (is.vector(meta_data)) {
    meta_data <- data.frame(CellID = 1:length(meta_data), Meta = meta_data)
  } else if (is.data.frame(meta_data) && !"CellID" %in% colnames(meta_data)) {
    meta_data$CellID <- 1:nrow(meta_data)
  }
  
  # Join metadata
  df_long <- df_long %>%
    left_join(meta_data, by = "CellID")
}
  # Apply log to HRS
  df_long$log_HRS <- log(df_long$HRS)
  
  # Create plots
  if (plot_type == "histogram") {
    p <- ggplot(df_long, aes(x = log_HRS)) +
      geom_histogram(bins = 30, fill = "steelblue", color = "black") +
      facet_wrap(~Variable, scales = "free") +
      labs(title = "Histogram of log(HRS)", x = "log(HRS)", y = "Count") +
      theme_minimal()
    
  } else if (plot_type == "density") {
    p <- ggplot(df_long, aes(x = log_HRS, fill = Variable)) +
      geom_density(alpha = 0.5) +
      facet_wrap(~Variable, scales = "free") +
      labs(title = "Density of log(HRS)", x = "log(HRS)", y = "Density") +
      theme_minimal()
    
  } else if (plot_type == "scatter" && !is.null(meta_data) && "Age_at_Enrollment" %in% colnames(meta_data)) {
    p <- ggplot(df_long, aes(x = Meta, y = log_HRS)) +
      geom_point(alpha = 0.4) +
      geom_smooth(method = "loess", se = FALSE, color = "darkred") +
      facet_wrap(~Variable, scales = "free") +
      labs(title = "log(HRS) vs Age", x = "Age", y = "log(HRS)") +
      theme_minimal()
    
  } else {
    stop("Invalid plot_type or missing metadata for scatter plot")
  }
  
  return(p)
}

#' Compute and Plot Hazard Risk Scores (HRS)
#'
#' This function calculates hazard risk scores (HRS) for multiple models, stores them in a data frame, 
#' and generates histograms and scatter plots to visualize the scores in relation to age.
#'
#' @param x A matrix of omics variables.
#' @param results A named list containing model coefficients for each risk score calculation.
#'
#' @return A data frame containing the computed risk scores for each model along with the age variable.
#' @export
#'
#' @examples
#' # Example usage:
#' HRS_results <- compute_and_plot_HRS(x, cox_elastic_output)
compute_HRS <- function(x, cox_elastic_output) {
  # Initialize an empty data frame to store ProtRS scores
  HRS_df <- list()
  
  # Loop through each column name in results
  for (col_name in names(cox_elastic_output)) {
    cat("Processing HRS for", col_name, "\n")
    
    col_label <- as.character(col_name)
    
    # Compute HRS for the current column
    HRS_col_name <- as.vector(exp(as.matrix(x) %*% cox_elastic_output[[col_name]]$coefficients))
    
    # Add the HRS to the data frame
    HRS_df[[paste0("HRS_", col_label)]] <- HRS_col_name
    
  }
  
  return(HRS_df)
}



#' Prepare Data for Cox Regression
#'
#' This function preprocesses omics data by applying log transformation, and standardizing the values.
#'
#' @param data A data frame containing the omics data.
#'
#' @return A standardized matrix combining transformed omics data.
#' @importFrom stats scale
#' @export
prepare_data <- function(data) {
  x <- log(data)
  x <- scale(x, center = TRUE, scale = TRUE)
  
  return(x)
}



#' Perform Kaplan-Meier Survival Analysis for multiple models
#'
#' This function computes Kaplan-Meier survival curves based on hazard risk scores (HRS) and 
#' stratifies individuals into low, medium, and high-risk groups.
#'
#' @param results A named list containing survival models.
#' @param HRS_df A data frame containing computed Risk Scores.
#' @param data A data frame with follow up time and event status variables (data output of `create_surv_models()`)
#'
#' @return A list of Kaplan-Meier survival plots for each risk score model.
#' @export
#'
#' @examples
#' km_plots <- compute_km_analysis(results, HRS_df, data)
compute_km_analysis <- function(HRS_df, data) {
  km_plots <- list()  # Store KM plots for each model
  
  for (col_name in names(HRS_df)) {
    cat("Kaplan-Meier analysis for", col_name, "\n")
    
    # Extract the RS scores for the current column
    HRS_col_name <- HRS_df[[col_name]]
    
    # Ranking and Stratification
    ranked_indices <- order(HRS_col_name, decreasing = TRUE)
    ranked_scores <- HRS_col_name[ranked_indices]
    
    # Stratify into risk groups (low, medium, high risk)
    quantiles <- quantile(HRS_col_name, probs = c(0.33, 0.66), na.rm = TRUE)
    risk_groups <- cut(
      HRS_col_name,
      breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
      labels = c("Low Risk", "Medium Risk", "High Risk")
    )
    
    # Ensure time and event variables are numeric and non-missing
    event <- as.integer(data[[paste0("Health_", sub("HRS_", "", col_name))]] == "Yes")
    time <- ifelse(
      data[[paste0("Health_", sub("HRS_", "", col_name))]] == "Yes",
      data[[paste0("Health_", sub("HRS_", "", col_name), "_Age")]],
      data$Age_last
    )
    valid_indices <- !is.na(time) & !is.na(event)
    
    time <- time[valid_indices]
    event <- event[valid_indices]
    risk_groups <- risk_groups[valid_indices]
    
    # Kaplan-Meier Survival Analysis
    surv_obj <- Surv(time = time, event = event)
    fit <- survfit(surv_obj ~ risk_groups)
    
    # Plot Kaplan-Meier survival curves
    plot_obj <- ggsurvplot(
      fit,
      data = data.frame(
        HRS_col_name = HRS_col_name[valid_indices],
        time,
        event,
        risk_groups
      ),
      pval = TRUE,  # Log-rank test p-value
      conf.int = TRUE,  # Confidence intervals
      risk.table = TRUE,  # Risk table
      ggtheme = theme_minimal(),  # Theme
      risk.table.height = 0.3,
      title = paste("Kaplan-Meier Survival Analysis for", col_name),
      palette = c("#2E9FDF", "#E7B800", "#FC4E07")  # Colors for risk groups
    )
    
    # Store plot in the list
    km_plots[[col_name]] <- plot_obj
  }
  
  return(km_plots)
}



#' Plot All Kaplan-Meier Survival Results
#'
#' This function prints all Kaplan-Meier survival plots stored in a list.
#'
#' @param km_results A list of Kaplan-Meier survival plots returned from `compute_km_analysis()`.
#'
#' @export
#'
#' @examples
#' plot_all_km_results(km_results)
plot_all_km_results <- function(km_results) {
  for (col_name in names(km_results)) {
    cat("Displaying Kaplan-Meier plot for", col_name, "\n")
    print(km_results[[col_name]])
  }
}



#' Compute Hazard Risk Scores (HRS) and perform Kaplan-Meier analysis
#'
#' This function fits a Cox proportional hazards model using elastic net regularization, 
#' extracts model coefficients, computes hazard risk scores (HRS), and performs Kaplan-Meier analysis 
#' with visualization.
#'
#' @param x A matrix of predictor variables (omics features).
#' @param data A data frame containing phenotype columns with survival times.
#' @param pheno_cols A character vector of column names representing phenotype survival times in `data`.
#' @param Age A numeric vector representing the age of individuals in the dataset.
#' @param alpha_values A numeric vector of alpha values to search over during Cox elastic net regression. 
#'   Default is `seq(0, 1, by = 0.1)`.
#'
#' @return A data frame containing computed risk scores for each phenotype model along with the age variable.
#' @export
#'
#' @examples
#' \dontrun{
#'   HRS_results <- get_Cox_HRS(x, data, pheno_cols, Age)
#' }
#' 
#' @import glmnet survival survminer
get_Cox_HRS <- function(x, data, pheno_cols, Age, alpha_values = seq(0, 1, by = 0.1)) {
  
   # Scale, log-normalize omics data
  x<- prepare_data(x)
  
  # Create survival models
  surv_models <- create_surv_models(data, pheno_cols)
  
  # Perform Cox Elastic Net Regression
  results <- optimize_cox_elastic_net(surv_models, x, alpha_values)

  # Compute and visualize Risk Scores (RS)
  HRS_df <- compute_and_plot_RS(x, results, Age)
  
  # Perform Kaplan-Meier survival analysis
  km_plots <- compute_km_analysis(results, HRS_df, data)
  
  # Display Kaplan-Meier plots
  plot_all_km_results(km_plots)
  
  return(HRS_df)
}


