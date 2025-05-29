# Changelog
All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/ ),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html ).

## [2.0.0] - 2025-04-05
### Changed
- Replaced ADASYN with **SMOGN** (a hybrid oversampling method combining SMOTE and Gaussian Noise based on distance between original samples) to improve compatibility and align with literature on handling highly imbalanced datasets.
- Reduced target class ratio from `1:1:1:1` to a more modest oversampling strategy that slightly increases minority class representation.
- Changed primary tuning metric from **MCC** to **KAPPA**, as KAPPA is more numerically stable and avoids issues with `NaN` values during optimization; final model selection still uses **MCC** for performance evaluation.

### Added
- Introduced **customized grid search** strategy aligned with relevant literature where applicable.
- Implemented **state-of-the-art probability calibration** — each model is calibrated using an isotonic regression to further reduce imbalance-related biases.

### Removed
- Removed the **bagged neural network model (`nnet` library)** due to outdated architecture, library incompatibility, and high computational cost. Deep learning remains a point of interest for future implementation.

### Experimental
- **MLP** and **TabPFN** models have been implemented and are functional, pending approval for integration into the main pipeline.
- **Improvements over matching flow algorithms**: A custom conditional matching flow implementation has been developed and is ready for validation and trial use in sample generation.
- Exploring a **new method for generating samples** via the **Forest Diffusion** library, which uses XGBoost internally to perform matching flow on tabular data — this approach shows promise and could be of interest for future oversampling strategies.