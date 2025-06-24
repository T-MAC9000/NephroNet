# NephroNet: Neural Network for Chronic Kidney Disease Progression Prediction

## Overview

NephroNet is a machine learning model designed to predict chronic kidney disease (CKD) progression using neural networks. The model predicts follow-up estimated glomerular filtration rate (eGFR) values based on clinical biomarkers and patient characteristics.

## Features

- **Neural Network Prediction**: Uses neural networks to predict eGFR progression
- **Feature Analysis**: Correlation analysis and SHAP value interpretation
- **Cross-Validation**: 10-fold cross-validation for robust model evaluation
- **Synthetic Data Generation**: Gaussian synthetic data generation for model validation
- **Comprehensive Visualization**: Multiple plotting functions for model interpretation
- **ROC Analysis**: Receiver Operating Characteristic curves for classification performance

## Data Requirements

### Input Features
- **Urea** (mmol/L)
- **Creatinine** (µmol/L) 
- **sTNFR1** (ng/mL) - Soluble TNF Receptor 1
- **sTNFR2** (ng/mL) - Soluble TNF Receptor 2
- **Gender** (Male/Female)
- **Baseline eGFR** (mL/min per 1.73m²)

### Target Variable
- **Follow-up eGFR** (mL/min per 1.73m²)

## Installation

### Required R Packages

```r
# Core packages
install.packages(c(
  "writexl", "ggplot2", "parallel", "Cairo", "doParallel",
  "reshape2", "corrplot", "finalfit", "fmsb", "FactoMineR",
  "neuralnet", "matrixStats", "pathfindR", "caret", "readxl",
  "fpc", "ggforce", "clValid", "umap", "magick", "dplyr",
  "pheatmap", "ggrepel", "tidyverse", "rstatix", "caTools",
  "rlang", "tidyr", "gtools", "ggstatsplot", "ggfortify",
  "factoextra", "boot", "cluster", "MASS", "iml", "pROC"
))
```

## Usage

### 1. Data Preparation

```r
# Load your data
MainData <- read_excel("path/to/your/data.xlsx")

# The script expects columns:
# - Patient ID
# - Age (years) 
# - Urea (mmol/l)
# - Creatinine (umol/L)
# - eGFR (baseline)
# - eGFR_Followup
# - Gender
# - Diabetes (E10-E14)
# - CVD
# - Biomarker columns (sTNFR1, sTNFR2)
```

### 2. Model Training

```r
# Set seed for reproducibility
set.seed(123)

# Train the neural network model
nn_model <- train(
  x = trainData[c(-1)],           
  y = trainData$eGFR_Followup,                          
  method = "nnet",                  
  trControl = tr_control,           
  linout = TRUE,                    
  trace = TRUE,                    
  maxit = 1000,     
  preProcess = c("center", "scale"),
  abstol = 1e-6,                   
  reltol = 1e-6,
  tuneLength = 200
)
```

### 3. Model Evaluation

```r
# Make predictions
predictions <- predict(nn_model$finalModel, testData)

# Evaluate performance
CA <- lm(log2(testData$eGFR_Followup) ~ log2(predictions))
summary(CA)
```

### 4. Feature Importance Analysis

```r
# SHAP values for feature importance
predictor <- Predictor$new(nn_model, data = trainData[, -1])
shapley <- Shapley$new(predictor, x.interest = trainData[1, -1])
```

## Key Functions

### Data Processing
- **calc_eGFR()**: Calculates eGFR using CKD-EPI 2021 race-free formula
- **generate_gaussian_synthetic()**: Generates synthetic data for validation

### Visualization
- **ggplotRegressionCA()**: Creates regression plots with statistics
- **Correlation heatmaps**: Feature correlation analysis
- **SHAP summary plots**: Feature importance visualization
- **ROC curves**: Classification performance evaluation

### Model Validation
- **Confusion matrices**: Classification accuracy assessment
- **Bootstrap confidence intervals**: ROC curve confidence bounds
- **Kolmogorov-Smirnov tests**: Synthetic vs. original data comparison

## Output Files

The script generates multiple visualizations saved as TIFF files:
- Correlation plots
- Regression analyses
- Feature importance plots
- ROC curves
- Distribution comparisons
- Confusion matrices

## Model Performance Metrics

- **R² values**: Regression performance
- **AUC scores**: Classification performance
- **Sensitivity/Specificity**: Diagnostic accuracy
- **Bootstrap confidence intervals**: Statistical reliability

## File Structure

```
project/
├── data/
│   ├── CKD EPI.xlsx
│   ├── Feature Cross.xlsx
│   └── Anonymous_Final Database.xlsx
├── models/
│   └── saved_models.rds
├── figures/
│   ├── correlations/
│   ├── regressions/
│   └── validation/
└── README.md
```

## Notes

- **Parallel Processing**: The script uses parallel processing for faster model training
- **Cross-Validation**: 10-fold cross-validation ensures robust model evaluation
- **Data Transformation**: Log2 transformation applied to continuous variables
- **Synthetic Validation**: Gaussian synthetic data generation for external validation
- **Reproducibility**: Set seed values ensure reproducible results

## Clinical Relevance

NephroNet addresses the critical need for predicting CKD progression by:
- Integrating multiple biomarkers (sTNFR1, sTNFR2)
- Providing interpretable predictions through SHAP analysis
- Supporting clinical decision-making for patient monitoring
- Enabling early intervention strategies

## Limitations

- Model trained on specific patient populations
- Requires complete biomarker data
- Performance may vary across different clinical settings
- External validation recommended before clinical deployment
- Small sample sizes may capture untrue global patterns

## Citation

If you use this code or methodology, please cite the associated research paper and acknowledge the NephroNet framework; #UPDATE WHEN PUBLISHED#

## Support

For questions or issues:
1. Check the code comments for detailed explanations
2. Verify all required packages are installed
3. Ensure data format matches expected structure
4. Review file paths and adjust for your system

