# NephroNet: Neural Network for Kidney Function Prediction

## Overview

NephroNet is a machine learning project that uses neural networks to predict kidney function decline in patients with Chronic Kidney Disease (CKD). The model predicts follow-up estimated Glomerular Filtration Rate (eGFR) values using clinical biomarkers and patient characteristics.

## Features

- **Neural Network Model**: Uses `caret` package with `nnet` method for regression
- **Feature Analysis**: Correlation analysis and feature importance using SHAP values
- **Model Validation**: Multiple validation approaches including external cohort validation
- **Visualization**: Comprehensive plotting for model performance, feature relationships, and predictions
- **ROC Analysis**: Binary classification for kidney function decline prediction
- **Synthetic Data Generation**: Gaussian synthetic data generation for model validation

## Dependencies

### Required R Packages

```r
# Core ML and Data Processing
library(caret)
library(neuralnet)
library(parallel)
library(doParallel)

# Data Manipulation
library(dplyr)
library(tidyverse)
library(tidyr)
library(readxl)
library(writexl)

# Visualization
library(ggplot2)
library(ggrepel)
library(ggforce)
library(ggstatsplot)
library(ggfortify)
library(pheatmap)
library(corrplot)

# Statistical Analysis
library(reshape2)
library(rstatix)
library(boot)
library(cluster)
library(fpc)
library(pROC)

# Feature Analysis
library(iml)          # For SHAP analysis
library(factoextra)
library(FactoMineR)

# Additional Utilities
library(Cairo)
library(magick)
library(caTools)
library(MASS)
```

## Data Requirements

### Input Data Structure

The main dataset should be an Excel file with the following columns:

- `Patient ID`: Unique identifier
- `Age (years)`: Patient age
- `Urea (mmol/l)`: Serum urea levels
- `Creatinine (umol/L)`: Serum creatinine levels
- `eGFR`: Baseline estimated Glomerular Filtration Rate
- `eGFR_Followup`: Follow-up eGFR (target variable)
- `Gender`: Patient sex (Male/Female)
- `Diabetes (E10-E14)`: Diabetes status
- `CVD`: Cardiovascular disease status
- `sTNFR1 (ng/ml)`: Soluble TNF Receptor 1 levels
- `sTNFR2 (ng/ml)`: Soluble TNF Receptor 2 levels

### File Paths

Update the following file paths in the script:
```r
# Main dataset
MainData <- read_excel("path/to/your/CKD_EPI.xlsx")

# Feature correlation data
Feats <- read_xlsx("path/to/your/Feature_Cross.xlsx")

# External validation cohort
CVDCohort <- read_xlsx("path/to/your/Anonymous_Final_Database.xlsx")
```

## Usage

### 1. Data Preparation

```r
# Set seed for reproducibility
set.seed(123)

# Load and clean data
MainData <- read_excel("your_data_file.xlsx")
MainData <- na.omit(MainData)

# Split data into train/test
splitR <- sample.split(DataFrame$eGFR_Followup, SplitRatio = 0.70)
trainData <- subset(DataFrame, splitR == TRUE)
testData <- subset(DataFrame, splitR == FALSE)
```

### 2. Model Training

```r
# Define training control
tr_control <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 3,
  allowParallel = TRUE,
  verboseIter = TRUE,
  search = "random"
)

# Train neural network
nn_model <- train(
  x = trainData[c(-1)],
  y = trainData$eGFR_Followup,
  method = "nnet",
  trControl = tr_control,
  linout = TRUE,
  trace = TRUE,
  maxit = 1000,
  tuneLength = 200,
  metric = "RMSE",
  preProcess = c("center", "scale")
)
```

### 3. Model Evaluation

```r
# Make predictions
predictions <- predict(nn_model, testData)

# Calculate performance metrics
CA <- lm(testData$eGFR_Followup ~ predictions)
summary(CA)

# Generate performance visualizations
ggplotRegressionCA(CA)
```

### 4. Feature Analysis

```r
# SHAP analysis for feature importance
predictor <- Predictor$new(nn_model, data = trainData[, -1])
shapley <- Shapley$new(predictor, x.interest = trainData[1, -1])

# Generate SHAP summary plot
SHAP_plot <- ggplot(SHAP_long, aes(x = `SHAP Value`, y = Feature, color = `Feature Value`)) +
  geom_jitter(width = 0.1, alpha = 1) +
  scale_color_gradient2(low = "dodgerblue3", mid = "gray90", high = "firebrick3")
```

## Key Functions

### Model Visualization
- `ggplotRegressionCA()`: Creates scatter plots with regression lines and statistics
- ROC curve generation for binary classification of kidney function decline
- Confusion matrix visualization for trend prediction

### Feature Analysis
- Correlation analysis between input features and target variable
- SHAP (SHapley Additive exPlanations) values for feature importance
- Individual feature relationship plots

### External Validation
- CVD cohort validation using separate dataset
- Synthetic data generation using multivariate Gaussian distribution
- Distribution comparison between original and synthetic data

## Output Files

The script generates several output files:

- **Model Performance**: Regression plots, ROC curves, confusion matrices
- **Feature Analysis**: Correlation plots, SHAP summary plots, individual feature relationships
- **Validation Results**: External cohort performance, synthetic data comparisons
- **Residual Analysis**: Residual plots for model diagnostics

All plots are saved as high-resolution TIFF files (350 DPI) suitable for publication.

## Model Architecture

- **Algorithm**: Neural Network (single hidden layer)
- **Activation**: Linear output layer for regression
- **Training**: 10-fold cross-validation with 3 repeats
- **Hyperparameter Optimization**: Random search with 200 iterations
- **Preprocessing**: Center and scale normalization
- **Early Stopping**: Absolute and relative tolerance thresholds

## Performance Metrics

- **Regression**: R², RMSE, slope, intercept
- **Classification**: AUC, sensitivity, specificity, accuracy
- **Feature Importance**: Mean absolute SHAP values
- **Model Validation**: External cohort performance, residual analysis

## Clinical Applications

This model can be used for:
- Predicting future kidney function in CKD patients
- Identifying patients at risk of rapid kidney function decline
- Supporting clinical decision-making in nephrology
- Research into biomarkers associated with kidney disease progression


## Contact

[mclarnon-t1@ulster.ac.uk]

## Version History

- v1.0: Initial release with basic neural network implementation
- v1.1: Added SHAP analysis and external validation
- v1.2: Implemented synthetic data generation and ROC analysis

## Troubleshooting

### Common Issues

1. **Memory Issues**: Reduce `tuneLength` parameter if running out of memory
2. **Convergence Problems**: Increase `maxit` or adjust tolerance parameters
3. **Missing Data**: Ensure all required columns are present in input data
4. **File Path Errors**: Update all file paths to match your directory structure

### Performance Optimization

- Use parallel processing for faster training
- Adjust cross-validation folds based on dataset size
- Consider feature selection to reduce model complexity

