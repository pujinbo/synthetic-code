*******************************************************************************
* Article: Synthetic Control Estimator: A Generalzied Inference Procedure and
*          Confidence Sets
* Authors: Sergio Firpo and Vitor Possebom
* Code by: Vitor Possebom
* Goal: We use the California exercise to illustrate how to estimate confidence
*       sets for the synthetic control estimator using Stata.
*******************************************************************************
* Clean and organize the work environment
clear all
set more off
*******************************************************************************
* Load and organize the data
sysuse smoking
tsset state year
* Estimate the SC weights for each state
forvalues i = 1/39 {
	synth cigsale beer(1984(1)1988) lnincome retprice age15to24 cigsale(1988) cigsale(1980) cigsale(1975), trunit(`i') trperiod(1989)
	* Save the matrix of observed outcomes
	matrix temp1 = e(Y_treated)
	matrix Ymat = (nullmat(Ymat), temp1)
	* Save the matrix of weights
	matrix temp1 = e(W_weights)
	matrix temp2 = temp1[1...,2]
	matrix weightsmat = (nullmat(weightsmat), temp2)
}
********************************************************************************
* Define the vector v of best (worst) case scenario.
matrix v = J(1, colsof(Ymat), 0)
* Run the do-file with the function SCMCS. You have to change the next line to
* the path of the file function_SCM-CS_v05_stata.do.
do "C:\Users\VitorPossebom\Dropbox\artigos_proprios\published_articles\SCE_GIP_CS\function_confidence_sets\Stata_files\function_SCM-CS_v07_stata.do"
* Compute the confidence sets
SCMCS Ymat weightsmat 3 19 0 v 30 "constant" 4/39
