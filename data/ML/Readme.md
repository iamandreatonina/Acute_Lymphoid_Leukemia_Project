# File di output from ML on HS

## Multiclass Classification Evaluation Metrics Explanation

## Overview

This document explains the evaluation metrics commonly used for multiclass classification problems. These metrics offer insights into different aspects of the model's performance, helping in understanding its strengths and weaknesses across various evaluation criteria.

## Metrics Explanation

### Accuracy
- **Description**: Provides an overall measure of how often the model's predictions are correct across all classes.
- **Focus**: Simple and intuitive but can be misleading in the presence of class imbalance.

### Kappa (KAP)
- **Description**: Measures the agreement between the model's predictions and the actual labels, considering the possibility of agreement occurring by chance.
- **Focus**: Provides a robust measure of classification accuracy, especially for imbalanced datasets.

### Sensitivity (SENS)
- **Description**: Quantifies the model's ability to correctly identify positive instances out of all actual positive instances.
- **Focus**: Minimizing false negatives.

### Specificity (SPEC)
- **Description**: Measures the model's ability to correctly identify negative instances out of all actual negative instances.
- **Focus**: Crucial in scenarios where false positives are costly.

### Positive Predictive Value (PPV)
- **Description**: Measures the accuracy of positive predictions made by the model.
- **Focus**: Minimizing false positives.

### Negative Predictive Value (NPV)
- **Description**: Measures the accuracy of negative predictions made by the model.
- **Focus**: Minimizing false negatives.

### Matthews Correlation Coefficient (MCC)
- **Description**: Measures the quality of binary classifications, considering both false positives and false negatives.
- **Range**: From -1 to 1 (1 indicates perfect prediction, 0 is no better than a random prediction, and -1 indicates total disagreement).

### Jaccard Index (J_INDEX)
- **Description**: Measures the similarity between two sets by calculating the ratio of the intersection to the union of the sets.
- **Focus**: Evaluate the overlap between the predicted and actual labels.

### Balanced Accuracy (BAL_ACCURACY)
- **Description**: Calculates the arithmetic mean of sensitivity and specificity.
- **Focus**: Provides a balanced measure that accounts for class imbalance.

### Detection Prevalence
- **Description**: Indicates the proportion of actual positive cases in the dataset.
- **Focus**: Helps in understanding the dataset's imbalance.

### Precision
- **Description**: Measures the accuracy of positive predictions focusing on the ratio of true positives to the sum of true and false positives.

### Recall
- **Description**: Measures the ratio of true positives to the sum of true positives and false negatives.
- **Focus**: Emphasizes the model's ability to identify all positive instances.

### F-measure (F_MEAS)
- **Description**: Harmonic mean of precision and recall.
- **Focus**: Provides a balanced evaluation metric considering both false positives and false negatives.

---
## Feature Importance in Random Forest (RF)

### Definition
Feature importance in RF quantifies each feature's contribution to the model's predictions based on the decrease in impurity when splitting the data.

### Calculation
- Permutation method, as implemented in the [Ranger library](https://github.com/imbs-hl/ranger) in R.
- Local permutation enables to recover the class-specific importance.
  
### Interpretation
- Higher importance indicates more influence on predictions.
- Lower importance suggests less relevance to the target variable.
- Negative importance suggests a negative impact on prediction.

### Visualization
- Importance scores can be visualized using bar plots in the related [image folder](https://github.com/iamandreatonina/Acute_Lymphoid_Leukemia_Project/tree/2f02a06eda1f30ab781c197093682e0671d7637c/Images/ML_importance).

### Importance in Practice
- **Insights**: Reveals relationships between features and the target.
- **Feature Selection**: Guides feature selection for improved model efficiency.
- **Interpretability**: Enhances model interpretability by highlighting influential features.



