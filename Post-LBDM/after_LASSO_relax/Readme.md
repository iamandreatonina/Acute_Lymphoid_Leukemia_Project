# Cross-Validated Lasso Model in R

This R code snippet demonstrates the construction of a cross-validated Relax Lasso model using the [`glmnet`](https://glmnet.stanford.edu/index.html) package. 
The Lasso model is a linear regression model with L1 regularization, commonly used for feature selection and regularization in predictive modeling.

## Code Explanation

The code is designed to fit a cross-validated Lasso model and then use the optimal lambda value for the final model fitting. Here's a breakdown of the key steps:

### 1. Cross-Validated Lasso Model

```R
# Fit cross-validated Lasso model
cv_lasso <- cv.glmnet(x = X,
                      y = Lasso_data$Cell_type, 
                      alpha = 0.5, 
                      grouped = FALSE,
                      parallel = TRUE,
                      relax = TRUE,
                      type.measure = "mse",
                      family = "gaussian",
                      nfolds = 15)
```
- The cv.glmnet function is used to perform cross-validated Lasso regression.
- X is the predictor matrix, and Lasso_data$Cell_type is the response variable.
- L1 regularization is controlled by the alpha parameter (set to 0.5 to do elastic-net).
- Other optional parameters include parallel processing, [`relaxation`](https://glmnet.stanford.edu/articles/relax.html), type of measurement (mean squared error), and the number of folds for cross-validation.

![mse](mse_cv_LASSO.png)

2. Fitting Lasso Model with Optimal Lambda

```R
optimal_lambda <- cv_lasso$lambda.min
```
The optimal regularization parameter (lambda) is selected based on the minimum mean squared error obtained during cross-validation.

```R
lasso_model <- glmnet(X, y, 
                      alpha = 0.5 
                      lambda = optimal_lambda,
                      family = "gaussian",
                      parallel = TRUE,
                      type.measure = "mse",
                      relax = TRUE)
```
- The glmnet function is used to fit the Lasso model with the optimal lambda value.
- X is the predictor matrix, and y is the response variable.
- Other parameters include the alpha value, family distribution, parallel processing, type of measurement (mean squared error), and relaxation.

[Coef](Post-LBDM/after_LASSO_relax/coef_lasso.png)
