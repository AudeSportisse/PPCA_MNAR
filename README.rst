=======
Estimation and imputation in Probabilistic Principal Component Analysis with Missing Not At Random data
=======

This repository hosts the code to impute containing Missing Not At Random (MNAR) values by assuming a PPCA model and to estimate the PPCA parameters.

There are the following Notebooks: 

* **Notebook_synteticdata.Rmd** for comparing several methods on synthetic data in terms of estimation/imputation. 

* **Notebook_realdata.Rmd** for comparing several methods on real data in terms of imputation (by adding additional MNAR values and computing the prediction error for those values). 

The code to reproduce the simulations of the paper is available in: 

* **Simulations_Toy.R** for Figure 2 (simple setting n=1000, p=10, r=2, 7 self-masked MNAR missing variable). 

* **Simulations_noise.R** for Figures 5, 6 and 7 (testing the robustness to the noise). 

* **Simulations_percentageNA.R** for Figure 8 (varying the percentage of missing values). 

* **Simulations_fixedeffect.R** for Figure 9 (misspecification to the PPCA model, by considering a fixed effect model). 

* **Simulations_misspe_rank.R** for Figure 10 (misspecification to the rank).

* **Simulations_MNARgeneral.R** for Figure 11 (the missing variables are MNAR (and not self-masked MNAR)).

* **Simulations_HighDimension.R** for Figures 12, 13, 14, 15, 16 and 17 (higher dimension setting, n=1000, p=50).

* **Simulations_MNARparam_parallel.R** (if Method (iv) MNARparam is compared). 



