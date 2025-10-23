Replication Code for Combining Matching and Synthetic Controls to Trade off Biases from Extrapolation and Interpolation
Kellogg, Mogstad, Pouliot, and Torgovitsky

This directory contains within it 6 main code files and a master file which runs them all.

The 6 files can be separated into two groups.

The first group contains code which runs the simulations and application exercises. The directory contains output from these codes already saved.
Running the first group requires gurobi.

The second group takes output from the first group, to make the plots for the paper. These make use of the 'textables' package (https://github.com/setzler/textables).

GROUP 1:
-SC_application.R: runs the comparative case study of Abadie and Gardeazabal (2003), along with placebo exercises based on that dataset.
-SC_EmpiricalARadjust_SNOW_byK.R: runs empirical monte carlo simulations, stores output

NOTE: The second program is set up to run on a linux cluster. Running it locally will require commenting out a few commands.

GROUP 2:
-plots_empiricalMC_descriptive.R: makes plots which describe the model used in the empirical monte carlo
-plots_empiricalMC_results.R: makes all plots based on the empirical monte carlo results
-plots_application.R: makes all plots for the application and associated empirical placebos
-plots_illustrative.R: makes plots of specific draws from the empirical monte carlo, for Section 2.

Additionally, support files and small datasets are found in the -auxfiles- folder:

-Estimator_Code.R: main code for solving all estimators given tuning parameters
-gurobiSC.R: code for a version of the SC estimator implemented using gurobi
-EmpiricalARoutcome_Code.R: code defining how we fit and draw from the DGP for the empirical monte carlo
-Cross-validation_Code_byK.R: code defining the cross-validation procedure by which we select tuning parameters for each estimator
-plot_rules.R: aesthetic rules applied to the plots.


