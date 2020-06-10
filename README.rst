=======
Estimation and imputation in Probabilistic Principal Component Analysis with Missing Not At Random data
=======

This repository hosts the code to impute containing Missing Not At Random (MNAR) values by assuming a PPCA model and to estimate the PPCA parameters.

There are the following Notebooks: 

* **Notebook_synteticdata.Rmd** for comparing several methods on synthetic data in terms of estimation/imputation. 

* **Notebook_realdata.Rmd** for comparing several methods on real data in terms of imputation (by adding additional MNAR values and computing the prediction error for those values). 

The code to reproduce the simulations of the paper is available in: 

* **Simulations_Toy.R** for Figure 2 (simple setting $$r=2$$, )

* **Simulations_noise.R** for Figures 5, 6 and 7 (testing the robustness to the noise)

* **Simulations_percentageNA.R**

* **Simulations_fixedeffect.R**

* **Simulations_misspe_rank.R**

* **Simulations_HighDimension.R**

* **Simulations_MNARparam_parallel.R**



