# 🏌️‍♂️ Functional Data Analysis of Golf Swing Motion

### Overview
This project applies **Functional Data Analysis (FDA)** techniques to **3D golf swing motion data** to understand and classify swing dynamics.  
It leverages smoothing, registration, and Functional Principal Component Analysis (FPCA) to capture continuous movement patterns and identify key sources of variation among golfers.

![Sample 3D Trajectory](images/YOUR_IMAGE_NAME_HERE.png)

---

## 📂 Project Structure

```
.
├── Data/                        # Raw and processed 3D swing data
├── Plots/                       # Output visualizations and diagnostics
├── 3D.R                         # Data import and initial processing
├── eda.R                         # Exploratory Data Analysis (EDA)
├── EDA2.R                        # Extended diagnostics and visualization
├── smoothing and registration.R  # Functional smoothing & curve alignment
├── Registration all.R            # Multi-axis registration process
├── fpca_xonly.R                  # FPCA on X-axis
├── mfpca.R                       # Multivariate FPCA on X, Y, Z axes
├── modelling.R                   # Final model training & classification
├── modelling_sample.R            # Example modeling pipeline
├── registered_xyz_fd.rds         # Cached registered functional data
├── smoothed_*_all_subjects.rds   # Smoothed curves for all subjects
└── FDA.Rproj                     # RStudio project file
```

---

## ⚙️ Analysis Workflow

1. **Data Loading & Preparation**  
   Load raw 3D trajectory data of golf swings, convert them into functional data objects using basis expansions.

2. **Exploratory Data Analysis (EDA)**  
   Visualize raw trajectories and compute summary statistics. Identify differences in swing shapes, timing, and amplitude.

   ![EDA Plot Example](images/YOUR_EDA_PLOT_IMAGE_HERE.png)

3. **Smoothing & Registration**  
   Apply smoothing to reduce noise and perform **time-warping alignment** to synchronize key swing events (e.g., impact point).

4. **Functional Principal Component Analysis (FPCA & mFPCA)**  
   - Perform FPCA for each axis (X, Y, Z).  
   - Conduct multivariate FPCA to capture correlated movements across axes.  
   - Extract principal components explaining swing variability.

   ![FPCA Plot Example](images/YOUR_FPCA_PLOT_IMAGE_HERE.png)

5. **Modeling & Classification**  
   - Use FPCA scores as input features for classification (e.g., expert vs novice).  
   - Evaluate models using accuracy, confusion matrices, and cross-validation.

6. **Visualization & Reporting**  
   - Visualize modes of variation, reconstructed curves, and performance plots.  
   - Output stored in the `Plots/` directory.

---

## 🧠 Key Insights
- **Functional representation** provides smoother, interpretable curves compared to discrete motion capture points.  
- **Registration** ensures swing phases align properly for comparison.  
- **FPCA** reveals dominant modes of variation — such as timing differences, backswing length, and rotation amplitude.  
- **Trip distance analogs** (movement amplitude) strongly predict skill or consistency.

---

## 🧰 Technologies Used
- **R**  
- **Packages:**  
  - `fda`, `fdapace`, `refund`, `tidyverse`, `mgcv`  
  - `ggplot2`, `gridExtra`, `reshape2`, `rstan` (optional for Bayesian modeling)  
- **RStudio** for workflow management and visualization.

---

## 🚀 How to Run

### Prerequisites
Ensure the following packages are installed:
```r
install.packages(c("fda", "fdapace", "tidyverse", "refund", "mgcv", "ggplot2", "reshape2", "gridExtra"))
```

### Steps
1. **Clone the repository**
   ```bash
   git clone https://github.com/adw4ith/FDA-GolfSwing-Analysis.git
   cd FDA-GolfSwing-Analysis
   ```

2. **Open the R project**
   Open `FDA.Rproj` in **RStudio**.

3. **Run the scripts sequentially:**
   ```r
   source("3D.R")
   source("eda.R")
   source("smoothing and registration.R")
   source("Registration all.R")
   source("fpca_xonly.R")
   source("mfpca.R")
   source("modelling.R")
   ```

4. **Check output**
   Processed `.rds` files (smoothed and registered data) and plots will be saved automatically in the `Data/` and `Plots/` folders.

---

## 📈 Future Work
- Integrate more complex classifiers (e.g., Random Forest, SVM).  
- Use dynamic time warping (DTW) for enhanced alignment.  
- Apply 3D visualization libraries to animate swing trajectories.  
- Evaluate swing consistency metrics across sessions.

---

## 🖼 Suggested Images to Add
You can include the following illustrative images in your `images/` folder:
1. `3d_swing_trajectory.png` — 3D trajectory of a sample swing  
2. `eda_summary.png` — Distribution or overlay of swings before/after registration  
3. `fpca_modes.png` — Visualization of first 2–3 FPCA modes  
4. `registration_alignment.png` — Before vs after alignment comparison  
5. `model_performance.png` — Model accuracy or confusion matrix  

Update their file names in the markdown under each placeholder.

---

## 👨‍💻 Author
**Adwaith T. S.**  
MSc Data Science Project — Functional Data Analysis of Golf Swing Motion  
