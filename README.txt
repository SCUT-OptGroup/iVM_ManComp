This code package contains implementations of several algorithms for related problems, along with demo and plotting scripts for reproducing the numerical experiments in the corresponding paper.

------------------------------------------------------------------------
Directory Structure
------------------------------------------------------------------------
The package includes five main folders:

1.  data/
    - Contains real-world datasets used for SSC experiments.

2.  functions/
    - Contains auxiliary/support functions used across all algorithms.

3.  SparSC/
    - Code for the SparSC algorithm (corresponds to Example 1.1 in the paper).

4.  r21PCA/
    - Code for the r21PCA algorithm (corresponds to Example 1.2 in the paper).

5.  Sdecom/
    - Code for the Sdecom algorithm (corresponds to Example 1.3 in the paper).

------------------------------------------------------------------------
Demo and Plotting Scripts
------------------------------------------------------------------------
The main directory contains several MATLAB scripts for running experiments and visualizing results:

- demo_realSSC.m
  Demo script for running SSC on real datasets.

- demo_r21PCA.m
  Demo script for running the r21PCA algorithm.

- demo_SDec_iter.m
  Demo script for for running SDec algorithm with $\lambda=0$.

- demo_SDecFL1_iter.m
  Demo script for for running SDec algorithm with $\lambda\neq 0$.

- plot_rsPCA_SNcg_APG.m
  Script for plotting results of the rsPCA with SNcg-APG method.

- plot_con_rate.m
  Script for plotting convergence rate curves.

- plot_synSSC.m
  Script for plotting SSC results on synthetic datasets.

- plot_gam_SDec.m
  Script for plotting the effect of the gamma parameter in SDec.

- plot_rsPCA_obj.m
  Script for plotting the objective function values and iterate times of rsPCA.

- plot_rsPCA_rho.m
  Script for plotting the effect of the $\rho$ parameter in rsPCA.

- plot_rsPCA_lambda.m
  Script for plotting the effect of the $\lambda$ parameter in rsPCA.

------------------------------------------------------------------------
Requirements
------------------------------------------------------------------------
- MATLAB (R2018b or later recommended)
- No additional toolboxes are required.

------------------------------------------------------------------------
Notes
------------------------------------------------------------------------
- Ensure the `data/` folder is in the correct path before running demos.
- Parameters such as regularization coefficients, maximum iterations, and convergence tolerances can be adjusted inside the demo scripts.
- Please send email to mahehao@mail.scut.edu.cn if you have any questions.

