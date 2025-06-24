##### Packages #####
library(writexl)
library(ggplot2)
library(parallel)
library(Cairo)
library(doParallel)
library(reshape2)
library(corrplot)
library(finalfit)
library(fmsb)
library(FactoMineR)
library(neuralnet)
library(matrixStats)
library(pathfindR)
library(caret)
library(readxl)
library(fpc)
library(ggforce)
library(clValid) # Retreieves CH argument for cluster.stats
library(umap)
library(magick)
#library(rgl)
library(dplyr)
library(pheatmap)
library(ggrepel)
library(tidyverse)
library(rstatix)
library(fpc)
library(caTools)
library(rlang)
library(tidyr)
library(gtools)
#library(hrbrthemes)
library(ggstatsplot)
library(ggfortify)
library(factoextra)
library(reshape2)
library(boot)
library(cluster)
set.seed(123)
##### ggsave #####
ggsave(
  'Influence.tiff',
  plot = ggplot2::last_plot(),
  device = NULL, ##### MAKE SURE TO FIX THIS CHANGE PLOT FROM RESIDUALS TO LASTPLOT
  path = "",
  scale = 1,
  width = 4000, #35
  height = 3000, #30
  units = c("px"),
  dpi = 350,
  limitsize = FALSE,
  bg = NULL,
)


ggplot2::last_plot()


residuals

dev.off()

last_plot()

##### Data entry #####
set.seed(123)

MainData<- as.data.frame(read_excel("C:\\Users\\Thomas\\Desktop\\CKD Sendotype Paper Figures High Quality\\NephroNet with Randox\\CKD EPI.xlsx"))
#Data Entry

MainData<-na.omit(MainData)
IDDD<-MainData$`Patient ID`
Age<-MainData$`Age (years)`
ID<-MainData[c(1)]
Urea<-MainData$`Urea (mmol/l)`
Creatinine1<-MainData$`Creatinine (umol/L)`
eGFR1<-MainData$eGFR
eGFR2<-MainData$eGFR_Followup
Gender<-MainData$Gender
Diabetes<-MainData$`Diabetes (E10-E14)`
CVD<-MainData$CVD
#Retrieval of Useful Information

Clinical<-data.frame(ID,Age,Urea,Creatinine1,eGFR1,eGFR2,Gender,MainData[c(4)],Diabetes,CVD)
#Merge of Useful Information

Data<-data.frame(Clinical,MainData[c(8,9)])
#Merge Randox Proteins

Data<-na.omit(Data)
ID<-Data[c(1)] # Keep ID
Data<-Data[c(-1)] # Remove ID


##### Wrangling #####

library(caret)
set.seed(123)


TestingCohort<-MainData

Feats<-read_xlsx("C:\\Users\\Thomas\\Desktop\\CKD Sendotype Paper Figures High Quality\\NephroNet with Randox\\Feature Corss.xlsx")
DataFrame<-MainData[,colnames(MainData) %in% Feats$Var1]

DataFrame<-data.frame(MainData[c(10)],MainData[c(11)],MainData[c(6)],DataFrame)

DataFrame$Creatinine..umol.L. <- DataFrame$Creatinine..umol.L. * 0.0113

splitR<-sample.split(DataFrame$eGFR_Followup, SplitRatio = 0.70)


trainData<-as.data.frame(subset(DataFrame, splitR == TRUE))
testData<-as.data.frame(subset(DataFrame, splitR == FALSE))

EGFR1Test<-testData$eGFR_Followup
TrainID<-as.data.frame(rep("Train",length(trainData$eGFR_Followup)))
TestID<-as.data.frame(rep("Test",length(testData$eGFR_Followup)))

colnames(TrainID)<-"ID"
colnames(TestID)<-"ID"

trainData<-cbind(TrainID,trainData)

testData<-cbind(TestID,testData)
WholeData<-rbind(trainData,testData)

WholeEGFR1<-WholeData[c(1,2,3)]

WholeData<-na.omit(WholeData)

NewWhole<-WholeData

IDs<-NewWhole[c(1)]

WholeD2<-NewWhole[c(-1)]

WholeD2Gender<-WholeD2[c(3)]

WholeD2<-WholeD2[c(-3)]

WholeD2<-log2(WholeD2)

WholeD2<-cbind(WholeD2,WholeD2Gender)

NewWholeD<-cbind(IDs,WholeD2)

BigBoySplit<-split(NewWholeD,NewWholeD$ID)

testData<-BigBoySplit$Test
trainData<-BigBoySplit$Train

testData<-testData[c(-1,-2)] # CHANGE BACK TO CREAT + PROTS
trainData<-trainData[c(-1,-2)]

formula <- eGFR2 ~.

##### Correlation Data #####

data <- MainData[c(-1)]  # Data minus ID

NewNames<-colnames(data)
NewNames[5]<-"Sex (Male Female)" # Change Gender to Sex per SAGER Guidelines

colnames(data)<-NewNames


correlation_matrix <- cor(data) #correlaton matrix

melted_corr <- melt(correlation_matrix) # format matrix
melted_corr<-as.data.frame(melted_corr)

Trial<-melted_corr[melted_corr[c(2)]=="eGFR_Followup",]
melted_Corr<-Trial

melted_Corr$value<-abs(melted_Corr$value)



melted_Corr<-arrange(melted_Corr,desc(melted_Corr$value))

Cols <- c("Significant Correlation"="#858585", "Insignificant Correlation" = "#90EE90")


melted_Corr<-melted_Corr[c(-1,-2,-4),]
melted_Corr<-melted_Corr %>%
  mutate(Colour=case_when(melted_Corr$value>0.5~"Significant Correlation",
                          TRUE~"Insignificant Correlation"))


melted_Corr %>%
  mutate(X1 = fct_reorder(Var1, value)) %>%
  ggplot(aes(x=value, y=X1, fill=Colour)) +
  geom_bar(stat="identity", alpha=.6, width=.8) +
  xlab("Absolute Pearson Correlation Coefficient")+
  ylab("Potential NephroNet Features")+
  ggtitle("Correlations between Input variables and follow-up eGFR")+
  theme(text = element_text(size=20),legend.position = "top", legend.justification.top = "left",
        axis.text.y = element_text( hjust=1,size=14))+
  scale_fill_manual(values = as.character(Cols))

corr_data<-split(melted_Corr,melted_Corr$Colour)
corr_data<-corr_data$`Significant Correlation`


#write_xlsx(corr_data,"")




##### NephroNet Model Building #####
# Define the training control with random search for hyperparameters
set.seed(87271)

tr_control <- trainControl(
  method = "cv",                  # Use cross-validation, can use repeatedcv too, add repeats=n,
  number = 10,                     # 10FCV
  allowParallel = TRUE,           # Enable parallel processing
  verboseIter = TRUE,         # Show progress of training
  search = "random",# Random search for hyperparameters
)


var<-trainData$eGFR_Followup

cores <- detectCores() - 1  # Detect number of available cores, using one less
cl <- makeCluster(cores)
registerDoParallel(cl)

nn_model <- train(
  x = trainData[c(-1)],           # Data - Var
  y = var,                          # Var
  method = "nnet",                  # caret nnet
  trControl = tr_control,           # tr_control
  linout = TRUE,                    # Linear activation func for back to continuous 
  trace = T,                    # Keep Trace toggled
  maxit = 1000,     
  preProcess = c("center", "scale"),
  abstol = 1e-6,                   # early stopping
  reltol = 1e-6, # Max iterations 
  tuneLength = 200    #tune length
)

stopCluster(cl) # stop Parallel
print(nn_model)

nnresults<-nn_model$finalModel
predictions<-predict(nnresults,testData)

testData$eGFR_Followup<-testData$eGFR_Followup
predictions<-predictions

CA <- lm(log2(testData$eGFR_Followup)~ log2(predictions))

ggplotRegressionCA <- function(CA) {
  
  ggplot(CA$model, aes(CA$model$`log2(testData$eGFR_Followup)`, CA$model$`log2(predictions)`)) + 
    geom_point(col="red",size=4) +
    geom_smooth(method = "lm", col = "black")+#,linetype = "dashed") +
    ylab("NephroNet Predicted log2(eGFR)") +
    xlab("Actual log2(eGFR)") +
    #scale_colour_manual(values = cols3)+
    ggtitle("Linear relationship of NephroNets Predictions of eGFR and Actual eGFR", 
            subtitle = paste("OLS R² =", signif(summary(CA)$adj.r.squared, 5),
                             "|  Intercept =", signif(CA$coef[[1]], 5),
                             "|  Slope =", signif(CA$coef[[2]], 5),
                             "|  p-value =", signif(summary(CA)$coef[2, 4], 5)))+
    guides(color = guide_legend(override.aes = list(size = 3)))+
    theme(
      plot.title = element_text(size = 20),  # Adjust title size
      plot.subtitle = element_text(size = 12),  # Adjust subtitle size
      axis.title.x = element_text(size = 18),  # Adjust x-axis label size
      axis.title.y = element_text(size = 18),  # Adjust y-axis label size
      axis.text.x = element_text(size = 10),  # Adjust x-axis tick label size
      axis.text.y = element_text(size = 10)  # Adjust y-axis tick label size
    )+
    theme_bw(base_size = 14)+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),  # Center title
          plot.subtitle = element_text(hjust = 0.5))  # Center subtitle
}

ggplotRegressionCA(CA)

#nn_model$trainingData

##### Baseline graph #####

predictions<-predict(nnresults,testData) # replace with testData, CVDTesting or synthetic_CVD


NewDf<-cbind(IDs,MainData)
NewDf<-split(NewDf,NewDf$ID)
NewDf<-NewDf$Test


Actual<-testData$eGFR_Followup
Pred<-predictions
Baseline<-NewDf$eGFR


Vis2<-as.data.frame(cbind(2^Actual,2^Pred,Baseline))
names(Vis2)<-c("Actual","Prediction","Baseline")

EGFRMetrics<-Vis2
names(EGFRMetrics)<-c("Actual","Prediction","Baseline")

EGFRMetrics <- EGFRMetrics %>%
  mutate(CASE = case_when(
    Prediction < Baseline & Actual < Baseline ~ "Correct Decline Prediction", # Both are declining
    Prediction >= Baseline & Actual >= Baseline ~ "Correct Stable Prediction", # Both are stable
    Prediction >= Baseline & Actual < Baseline ~ "Incorrect Stable Prediction", # Predicted stable but actually declining
    Prediction < Baseline & Actual >= Baseline ~ "Incorrect Decline Prediction",
    TRUE~"WHOOPSIES"# Predicted declining but actually stable
  ))

predictionss <- as.factor(rep(paste("Prediction", 1:length(Vis2$Actual)), length.out = n))


Vis2<-cbind(EGFRMetrics,predictionss)

Vis_long <- Vis2 %>%
  pivot_longer(cols = c("Actual", "Baseline", "Prediction"), 
               names_to = "Type", 
               values_to = "Value") %>%
  mutate(Predictions = predictionss)

data_long<-Vis_long


data_long<-arrange(data_long,data_long$Value)


data_long <- data_long %>%
  filter(Type == "Baseline") %>%
  arrange(Value) %>%
  mutate(predictionss = factor(predictionss, levels = unique(predictionss))) %>%
  dplyr::select(predictionss) %>%
  right_join(data_long, by = "predictionss")


ggplot(data_long, aes(y = predictionss)) +
  geom_segment(data = data_long %>% 
                 filter(Type %in% c("Baseline", "Prediction")) %>% 
                 pivot_wider(names_from = Type, values_from = Value),
               aes(x = Baseline, xend = Prediction, y = predictionss, yend = predictionss),
               color = "grey60", size = 0.5) +
  geom_segment(data = data_long %>% 
                 filter(Type %in% c("Prediction", "Actual")) %>% 
                 pivot_wider(names_from = Type, values_from = Value),
               aes(x = Prediction, xend = Actual, y = predictionss, yend = predictionss),
               color = "grey60", size = 0.5) +
  geom_point(aes(x = Value, color = Type, shape = Type), size = 3) +
  scale_color_manual(values = c("Baseline" = "#6a0dad", "Prediction" = "#ff7f0e", "Actual" = "#17becf")) +
  scale_shape_manual(values = c("Baseline" = 15, "Prediction" = 16, "Actual" = 17)) + # 15=square, 16=circle, 17=triangle
  theme_minimal() +
  labs(title = "Magnitudal Differences between eGFR",
       subtitle = "Baseline eGFR, 1Y Follow-Up, Predicted 1Y Follow-Up",
       x = "eGFR (mL/min per 1.73m²)",
       y = "Individual Patients") +
  theme_bw(base_size = 24) + 
  theme(legend.position = 'top', 
        legend.justification = 'left',
        legend.direction = 'horizontal') +
  theme(legend.title = element_text(size = 16)) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank()
  )




##### CVD Data Preprocessing #####

CVDCohort<-read_xlsx("C:\\Users\\Thomas\\Downloads\\Anonymous_Final Database.xlsx",sheet=1)
CVDID<-read_xlsx("C:\\Users\\Thomas\\Downloads\\CKD VHR(2).xlsx")
CVDCohort<-CVDCohort[CVDCohort$`Patient ID...2` %in% CVDID$`Patient ID`,]

sTNFR1<-CVDCohort$`sTNFR1 / TNFRSF1A (pg/ml) (R&D)`
sTNFR2<-CVDCohort$`sTNFR2 / TNFRSF1B (pg/ml) (R&D)`

Creat1CVD<-CVDID$`creatinine tp1`
Creat2CVD<-CVDID$`creatinine tp2`
AgeCVD<-CVDCohort$Age
GenderCVD<-CVDCohort$Sex

# Convert to mg/dL
Creat1CVD_mg_dl <- Creat1CVD / 88.4
Creat2CVD_mg_dl <- Creat2CVD / 88.4

# Normalize gender input


# CKD-EPI 2021 race-free function
calc_eGFR <- function(creatinine_mg_dl, age, sex){
  kappa <- ifelse(sex == "F", 0.7, 0.9)
  alpha <- ifelse(sex == "F", -0.241, -0.302)
  sex_factor <- ifelse(sex == "F", 1.012, 1)
  
  min_part <- pmin(creatinine_mg_dl / kappa, 1)
  max_part <- pmax(creatinine_mg_dl / kappa, 1)
  
  egfr <- 142 * (min_part ^ alpha) * (max_part ^ -1.2) * (0.9938 ^ age) * sex_factor
  return(egfr)
}

# Vectorized calculation
eGFR_CVD <- mapply(calc_eGFR, creatinine_mg_dl = Creat1CVD_mg_dl, 
                   age = AgeCVD, sex = GenderCVD)

eGFR_Followup <- mapply(calc_eGFR, creatinine_mg_dl = Creat2CVD_mg_dl, 
                        age = AgeCVD+1, sex = GenderCVD)

CVDTesting<-data.frame(eGFR_Followup,CVDID$`serum urea (mmol/L)`,Creat1CVD_mg_dl,sTNFR1,sTNFR2,eGFR_CVD)

CVDTesting<-log2(CVDTesting)

GenderCVD_num <- ifelse(GenderCVD == "F", 1, 0)
CVDTesting<-cbind(CVDTesting,GenderCVD_num)

g<-colnames(trainData)

g[6]<-"Baseline"
g[7]<-"Gender"
colnames(CVDTesting)<-g #replace testData for CVDTesting


##### CVD IMPUTATION #####

library(MASS)

generate_gaussian_synthetic <- function(data, n = 500) {
  numeric_data <- data[, sapply(data, is.numeric)]
  mu <- colMeans(numeric_data)
  sigma <- cov(numeric_data)
  
  synthetic_numeric <- mvrnorm(n = n, mu = mu, Sigma = sigma)
  synthetic_numeric <- as.data.frame(synthetic_numeric)
  
  # Add categorical (resampled) gender
  synthetic_numeric$Gender <- sample(data$Gender, n, replace = TRUE)
  
  return(synthetic_numeric)
}

# Usage
synthetic_CVD <- generate_gaussian_synthetic(CVDTesting, n = 500)

synthetic_CVD<-as.data.frame(synthetic_CVD) #replace testData for synthetic_CVD

Baseline<-synthetic_CVD$Baseline #replace CKD baseline in baseline graph for synth baseline graph


##### Linear Model #####

#rfe_results<-read_rds("C:\\Users\\Thomas\\Desktop\\CKD Sendotype Paper Figures High Quality\\Neural Network\\RFE Neural Network\\BEST MODEL\\RFE Model.rds")
#predictions<-predict(nn_model,testData)

#predictions<-predict(model$finalModel,testData)
predsCVD<-predict(nn_model$finalModel,synthetic_CVD) # replace with testData, CVDTesting or synthetic_CVD

CA <- lm(synthetic_CVD$eGFR_Followup~ predsCVD) # replace with testData, CVDTesting or synthetic_CVD


ggplotRegressionCA <- function(CA) {
  
  ggplot(CA$model, aes(CA$model$`synthetic_CVD$eGFR_Followup`, CA$model$predsCVD)) + 
    geom_point(col="red",size=4) +
    geom_smooth(method = "lm", col = "black")+#,linetype = "dashed") +
    ylab("NephroNet predicted log2(eGFRs)") +
    xlab("Actual log2(eGFRs)") +
    #scale_colour_manual(values = cols3)+
    ggtitle("Linear relationship of NephroNets Predictions of eGFR and Actual eGFR", 
            subtitle = paste("OLS R² =", signif(summary(CA)$adj.r.squared, 5),
                             "|  Intercept =", signif(CA$coef[[1]], 5),
                             "|  Slope =", signif(CA$coef[[2]], 5),
                             "|  p-value =", signif(summary(CA)$coef[2, 4], 5)))+
    guides(color = guide_legend(override.aes = list(size = 3)))+
    theme(
      plot.title = element_text(size = 20),  # Adjust title size
      plot.subtitle = element_text(size = 12),  # Adjust subtitle size
      axis.title.x = element_text(size = 14),  # Adjust x-axis label size
      axis.title.y = element_text(size = 14),  # Adjust y-axis label size
      axis.text.x = element_text(size = 10),  # Adjust x-axis tick label size
      axis.text.y = element_text(size = 10)  # Adjust y-axis tick label size
    )+
    theme_bw(base_size = 14)+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),  # Center title
          plot.subtitle = element_text(hjust = 0.5))  # Center subtitle
} # replace inside of function with testData, CVDTesting or synthetic_CVD

ggplotRegressionCA(CA)



##### Confusion Matrix for Classification #####
library(caret)
library(ggplot2)

# Copy dataset
tables <- Vis2

# Define directional labels
tables$ActualTrend <- ifelse(tables$Actual < tables$Baseline, "Declined", "Recovered")
tables$PredictedTrend <- ifelse(tables$Prediction < tables$Baseline, "Declined", "Recovered")

# Desired order: Declined = top-left
desired_order <- c("Declined", "Recovered")

# Set levels accordingly
tables$ActualTrend <- factor(tables$ActualTrend, levels = desired_order)
tables$PredictedTrend <- factor(tables$PredictedTrend, levels = desired_order)

# Confusion matrix (for stats)
cm <- confusionMatrix(tables$PredictedTrend, tables$ActualTrend)
print(cm)

# Create table for plot
cm_table <- table(Predicted = tables$PredictedTrend, Actual = tables$ActualTrend)
cm_df <- as.data.frame(cm_table)

# Reapply levels for ggplot control
cm_df$Predicted <- factor(cm_df$Predicted, levels = rev(desired_order))  # Flip y-axis to put Declined on top
cm_df$Actual <- factor(cm_df$Actual, levels = desired_order)

# Plot
ggplot(cm_df, aes(x = Actual, y = Predicted, fill = Freq)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Freq), size = 8) +
  scale_fill_gradientn(colors = rev(hcl.colors(25, "OrRd"))) +
  scale_y_discrete(limits = rev(desired_order)) + 
  theme_bw(base_size = 18) +
  labs(
    title = "Confusion Matrix: Predicted vs Actual Renal Trajectories",
    x = "Actual Trajectory (vs Baseline)",
    y = "Predicted Trajectory (vs Baseline)"
  )


##### Linear Models of Input Features (Continuous) #####

colnames(DataFrame)
Input<-DataFrame[c(-1)]


Names<-c("Follow-up eGFR   mL/min per 1.73m^2","Urea   mmol/L","Creatinine   µmol/L","Sex   Male/Female","sTNFR1   ng/mL","sTNFR2   ng/mL")


Names<-c("Follow up eGFR   mL/min per 1.73m^2","Sex   Male/Female","Urea    µmol/L","Creatinine   µmol/L","sTNFR1   ng/mL","sTNFR2   ng/mL")
DirNames<-c("Sex","Urea","Creatinine","sTNFR1","sTNFR2")

for (i in 1:length(DataFrame)){
  Var<-unlist(Input[c(1)])
  test<-unlist(Input[c(i+1)])
  CA <- lm(Var~test)
  
  
  
  ggplotRegressionCA <- function(CA) {
    
    ggplot(CA$model, aes(CA$model$Var, CA$model$test)) + 
      geom_point(col="red",size=4) +
      geom_smooth(method = "lm", col = "black")+#,linetype = "dashed") +
      ylab((paste0(Names[i+1]))) +
      xlab("Follow-up eGFR   mL/min per 1.73m^2") +
      #scale_colour_manual(values = cols3)+
      ggtitle("Linear relationship of Input Features and Follow-up eGFR", 
              subtitle = paste("OLS R² =", signif(summary(CA)$adj.r.squared, 5),
                               "|  Intercept =", signif(CA$coef[[1]], 5),
                               "|  Slope =", signif(CA$coef[[2]], 5),
                               "|  p-value =", signif(summary(CA)$coef[2, 4], 5)))+
      guides(color = guide_legend(override.aes = list(size = 3)))+
      theme(
        plot.title = element_text(size = 20),  # Adjust title size
        plot.subtitle = element_text(size = 12),  # Adjust subtitle size
        axis.title.x = element_text(size = 14),  # Adjust x-axis label size
        axis.title.y = element_text(size = 14),  # Adjust y-axis label size
        axis.text.x = element_text(size = 10),  # Adjust x-axis tick label size
        axis.text.y = element_text(size = 10)  # Adjust y-axis tick label size
      )+
      theme_bw(base_size = 14)+
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),  # Center title
            plot.subtitle = element_text(hjust = 0.5))  # Center subtitle
  }
  
  ggplotRegressionCA(CA)
  
  
  ggsave(
    paste0(DirNames[i],".tiff"),
    plot = ggplot2::last_plot(),
    device = NULL, ##### MAKE SURE TO FIX THIS CHANGE PLOT FROM RESIDUALS TO LASTPLOT
    path = "C:\\Users\\Thomas\\Desktop\\CKD Sendotype Paper Figures High Quality\\NephroNet with Randox\\Figures APRIL\\Input Feats",
    scale = 1,
    width = 3000, #35
    height = 3000, #30
    units = c("px"),
    dpi = 350,
    limitsize = FALSE,
    bg = NULL,
  )
  
  
}

##### Non-Linear Model of Input Features (Categorical) ##### 


# Ensure Gender is numeric (0 or 1)
Input$Gender <- as.numeric(Input$Gender)

# Fit logistic regression model
CA <- glm(Gender ~ Var, data = Input, family = binomial)

# Summary of the model
summary(CA)



ggplotRegressionCA <- function(CA) {
  
  # Compute McFadden's Pseudo R²
  pseudo_r2 <- 1 - (logLik(CA) / logLik(update(CA, . ~ 1)))
  
  ggplot(CA$model, aes(x = Var, y = Input$Gender)) + 
    geom_point(col = "red", size = 4) +
    geom_smooth(method = "glm", method.args = list(family = binomial), col = "black") +
    ylab("  Sex   0=Male 1=Female") +
    xlab("Follow-up eGFR   mL/min per 1.73m^2") +
    ggtitle("Non-Linear Relationship of Input Features and Follow-up eGFR", 
            subtitle = paste("Pseudo R² =", signif(pseudo_r2, 5),
                             "| Intercept =", signif(CA$coef[[1]], 5),
                             "| Slope =", signif(CA$coef[[2]], 5),
                             "| p-value =", signif(summary(CA)$coef[2, 4], 5))) +
    guides(color = guide_legend(override.aes = list(size = 3)))+
    theme(
      plot.title = element_text(size = 20),  # Adjust title size
      plot.subtitle = element_text(size = 12),  # Adjust subtitle size
      axis.title.x = element_text(size = 14),  # Adjust x-axis label size
      axis.title.y = element_text(size = 14),  # Adjust y-axis label size
      axis.text.x = element_text(size = 10),  # Adjust x-axis tick label size
      axis.text.y = element_text(size = 10)  # Adjust y-axis tick label size
    )+
    theme_bw(base_size = 14)+
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),  # Center title
          plot.subtitle = element_text(hjust = 0.5))
  
}

# Run the function
ggplotRegressionCA(CA)




write_rds(nn_model,"C:\\Users\\Thomas\\Desktop\\CKD Sendotype Paper Figures High Quality\\NephroNet with Randox\\Figures UPDATED\\CVD Validation\\CVD DIABETES\\CVD.rds")

##### ROC 2 #####
library(dplyr)
library(pROC)
library(ggplot2)

Df2 <- (Vis2)  # Assuming Vis2 contains eGFR Baseline, Predictions, and Actual values

# Step 1: Define Binary Outcome (Decline vs Stable/Increase)
Df2 <- Df2 %>%
  mutate(Decline = ifelse(Actual < Baseline, 1, 0))  # 1 = Decline, 0 = Stable/Increase

# Step 2: ROC Analysis using Predicted eGFR as a Threshold

#Df2$Actual <- Df2$Actual - Df2$Baseline
#Df2$Prediction <- Df2$Prediction - Df2$Baseline

# Normalize between 0 and 1
#normalize <- function(x) {
# return((x - min(x)) / (max(x) - min(x)))
#}

#Df2$Actual <- normalize(Df2$Actual)
#Df2$Prediction <- normalize(Df2$Prediction)

# Generate ROC Curve using normalized predicted values
roc_curve <- roc(Df2$Decline,Df2$Prediction)




# Step 3: Confidence Interval via Bootstrapping
x <- roc(Df2$Decline, Df2$Prediction, plot = TRUE,
         auc.polygon = TRUE, max.auc.polygon = TRUE, grid = TRUE,
         print.auc = TRUE, show.thres = TRUE, ci = TRUE,
         legacy.axes = TRUE)

dats <- ci.se(x, conf.level = 0.95, method = "bootstrap",
              boot.n = 5000, boot.stratified = TRUE, reuse.auc = TRUE,
              specificities = seq(0, 1, length.out = 25))

dat.ci <- data.frame(x = as.numeric(rownames(dats)),
                     lower = dats[, 1],
                     upper = dats[, 3])

# Step 4: Visualize ROC with Confidence Intervals
LRROC <- ggroc(x, legacy.axes = TRUE) +
  theme_minimal() +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", alpha = 0.7, color = "grey") +
  coord_equal() +
  geom_ribbon(data = dat.ci, aes(x = 1 - x, ymin = lower, ymax = upper), fill = "steelblue", alpha = 0.2) +
  xlab(expression("1-Specificity")) +
  ylab(expression("Sensitivity")) +
  ggtitle("NephroNet Renal Decline Prediction", 
          paste0("AUC: ", round(x$auc, 3), 
                 " (95% CI: ", round(x$ci[1], 3), " - ", round(x$ci[3], 3), ")"))+
  theme_bw(base_size = 16)

LRROC

##### SHAP 2 #####
library(iml)
library(ggplot2)
library(reshape2)
library(dplyr)

# Set up the predictor using the updated dataset
predictor <- Predictor$new(nn_model, data = trainData[, -1])
predictor$task <- "regression"

# Initialize empty list to collect SHAP values
shap_values_list <- list()

# Loop over observations to calculate SHAP values for each
for (i in 1:nrow(trainData)) {
  shapley <- Shapley$new(predictor, x.interest = trainData[i, -1])
  shap_row <- shapley$results[, c("feature", "phi")]
  shap_row_t <- setNames(as.data.frame(t(shap_row$phi)), shap_row$feature)
  shap_values_list[[i]] <- shap_row_t
  print(i)
}

# Combine SHAP values into a dataframe
SHAPDF <- do.call(rbind, shap_values_list)


# --- NEW: Extract the feature values (same shape) ---
FeatureVals <- trainData[, c("Urea..mmol.l.", "Creatinine..umol.L.", "sTNFR1..ng.ml.", "sTNFR2..ng.ml.", "Gender")]
FeatureVals <- FeatureVals[1:nrow(SHAPDF), ]  # Just in case

# Scale SHAP and Feature values individually by feature
SHAP_scaled <- SHAPDF
FeatureVals_scaled <- as.data.frame(lapply(FeatureVals, function(x) scale(x)))

# Melt both SHAP and feature values for ggplot
SHAP_long <- melt(SHAP_scaled, variable.name = "Feature", value.name = "SHAP Value")
FeatureVals_long <- melt(FeatureVals_scaled, variable.name = "Feature", value.name = "Feature Value")

# Combine the two data frames into one
SHAP_long$`Feature Value` <- FeatureVals_long$`Feature Value`

# Rename the features for clarity
SHAP_long$Feature <- recode(SHAP_long$Feature,
                            "Urea..mmol.l." = "Urea  µmol/L",
                            "Creatinine..umol.L." = "Creatinine  µmol/L",
                            "sTNFR1..ng.ml." = "sTNFR1  ng/mL",
                            "sTNFR2..ng.ml." = "sTNFR2  ng/mL",
                            "Gender" = "Sex  Male/Female")

# --- Final SHAP plot with renamed features and ordered by absolute SHAP importance ---
SHAP_plot <- ggplot(SHAP_long, aes(x = `SHAP Value`, y = Feature, color = `Feature Value`)) +
  geom_jitter(width = 0.1, alpha = 1) +
  scale_color_gradient2(low = "dodgerblue3", mid = "gray90", high = "firebrick3", midpoint = 0) +
  theme_bw(base_size = 16) +
  labs(
    title = "SHAP Summary Plot",
    subtitle = "(SHAP values calculated for all observations in train cohort)",
    x = "SHAP Value (mL/min per 1.73m^2)",
    y = ""
  ) +
  theme(axis.text.y = element_text(size = 16, family = "Arial")) +
  # Reorder features dynamically based on the mean absolute SHAP value (within ggplot call)
  scale_y_discrete(limits = names(sort(tapply(abs(SHAP_long$`SHAP Value`), SHAP_long$Feature, mean), decreasing = F)))

# Plot it
print(SHAP_plot)




varImp(nn_model$finalModel)
##### Importance Plot #####


Mean_Abs_SHAP <- tapply(abs(SHAP_long$`SHAP Value`), SHAP_long$Feature, mean)


feature_importance <- data.frame(
  Feature = names(Mean_Abs_SHAP),
  Mean_Abs_SHAP = as.numeric(Mean_Abs_SHAP)
)


feature_importance$Percent_Influence <- 100 * feature_importance$Mean_Abs_SHAP / sum(feature_importance$Mean_Abs_SHAP)


feature_importance$Feature <- factor(
  feature_importance$Feature,
  levels = feature_importance$Feature[order(feature_importance$Percent_Influence, decreasing = FALSE)]
)


ggplot(feature_importance, aes(x = Feature, y = Percent_Influence)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +
  coord_flip() +
  labs(
    title = "Feature Influence on Model Prediction",
    subtitle = "Percent Influence = Mean(|SHAP|) / Sum of Mean(|SHAP|) for all features",
    x = "",
    y = "Percentage Influence on Prediction (%)"
  ) +
  theme_bw(base_size = 18)



##### QC Histograms ####

Original<-CVDTesting
Synthetic<-synthetic_CVD



str(Original)
str(Synthetic)

library(ggplot2)
library(dplyr)
library(hrbrthemes)

ID1<-as.data.frame(rep("Original",length(Original$eGFR_Followup)))
colnames(ID1)<-"ID"
ID2<-as.data.frame(rep("Synthetic",length(Synthetic$eGFR_Followup)))
colnames(ID2)<-"ID"
ID<-rbind(ID1,ID2)
colnames(trainData)
Feats<-c("Follow up eGFR","Urea","Creatinine","sTNFR1","sTNFR2","Sex")
for(i in 1:length(Original)){
  
  FeatDf<-rbind(Original[c(i)],Synthetic[c(i)])
  IterDf<-data.frame(ID,FeatDf)
  str(IterDf)
  colnames(IterDf)<-c("Dataset","Feature")
  
  Org<-unlist(Original[c(i)])
  Syn<-unlist(Synthetic[c(i)])
  ksv<-ks.test(Org,Syn)
  ksv<-ksv$p.value
  IterDf$Dataset <- factor(IterDf$Dataset, levels = c("Synthetic", "Original"))  # Adjust the order here
  p <- IterDf %>%
    ggplot( aes(x=Feature, fill=Dataset)) +
    geom_histogram( color="#e9ecef", alpha=.5, position = 'stack',bins = 50) +
    scale_fill_manual(values=c("#404080", "#69b3a2")) +
    theme_ipsum() +
    labs(fill="")+
    ggtitle("Distribution of Variable differences between Original and Synthetic Cohorts",subtitle = paste0("Kolmogorov-Smirnov p-value: ",round(ksv,3)," Original (n=11), Synthetic (n=500)"))+
    theme_bw(base_size = 18)+
    xlab(paste0(Feats[i]))+
    ylab("Frequency")
  
  
  p
  
  
  ggsave(
    paste0(Feats[i],'.tiff'),
    plot = ggplot2::last_plot(),
    device = NULL, ##### MAKE SURE TO FIX THIS CHANGE PLOT FROM RESIDUALS TO LASTPLOT
    path = "C:\\Users\\Thomas\\Desktop\\CKD Sendotype Paper Figures High Quality\\NephroNet with Randox\\Figures APRIL\\CVD Validation\\Histograms",
    scale = 1,
    width = 4000, #35
    height = 3500, #30
    units = c("px"),
    dpi = 350,
    limitsize = FALSE,
    bg = NULL,
  )
  
}


# Example: comparing two continuous distributions
