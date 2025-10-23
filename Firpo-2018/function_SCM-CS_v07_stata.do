********************************************************************************
* Article: Synthetic Control Method: Inference, Sensitivity Analysis and
*          Confidence Sets
* Authors: Sergio Firpo and Vitor Possebom
* Code by: Vitor Possebom
* Function: Implement Confidence Sets and the Sensitivity Analysis Mechanism
*           for the Synthetic Control Estimator.
* Version: 06 - This is a preliminary version of a function to implement the
*          Confidences Sets and the Sensitivity Analysis Mechanism for the
*          Synthetic Control Estimator proposed by Firpo and Possebom (2017).
*          We use the RMSPE test statistic to compute confidence sets in this
*          code. The user may want to manually change the test statistic. Any
*          questions, suggestions, comments, critics can be sent to Vitor
*          Possebom (vitoraugusto.possebom@yale.edu). In particular, the user
*          may want to manually change the graphical parameters at the end of
*          this code. In the code, there are comments stressing lines that
*		   the user may have to change depending on their dataset.
*******************************************************************************
capture program drop SCMCS
program define SCMCS, rclass
	args Ymat weightsmat treated T0 phi v precision etype significance
* 'Ymat' is a matrix that contains the realized outcome of interest. Each
* column represents a region and each row represents a time period. It was
* defined above
matrix Y = `Ymat'

* 'weightsmat' is a matrix that contains the estimated weights for the synthetic
* control unit of each region. Each columns represents a placebo region and
* each row, a comparison unit. Pay attention that the order of the regions in
* matrix 'weights' must be identical to the order of the regions in matrix
* 'Y'. 'weights' is a matrix with J-1 rows and J columns, where J is the
* number of observed regions. In order to construct this matrix, one can
* estimate a synthetic control unit using the command 'synth' and save the
* weights found in the vector W_weights for each region. Then, each vector
* W_weights will be a column in matrix 'weights'. It was defined above.
matrix weights = `weightsmat'

* 'treated' is the column in matrix 'Y' associated with the treated
* region. It must be a natural number less than or equal to the number of
* columns in matrix 'Y'.

* 'T0' is the row in matrix 'Y' associated with the last pre-intervention
* period. It must be a natural number less than or equal to the number of
* rows in matrix 'Y'.

* 'phi' is the sensitivity parameter defined either in step 6 or in step 7 of
* section 3 of Firpo and Possebom (2017). It has to be a positive real number.
* If you only want to construct the standard confidence interval under the
* assumption that each region is equally likely to receive treatment, set
* 'phi' to zero.

* 'v' is the worst (best) case scenario vector defined in step 6 (step 7) of
* section 3 of Firpo and Possebom (2017). It is row vector whose length is
* equal to the number of observed regions, i.e., J+1. The elements of this
* vector must be equal to 0 or 1. If you only want to construct the standard
* confidence interval under the assumption that each region is equally likely
* to receive treatment, let 'v' be a zero vector.

* 'precision' is a natural number. A larger value for 'precision' makes the
* estimation of the confidence sets more precise, requiring more computing
* time. A value between 20 and 30 is reasonable.

* 'etype' is a character vector equal to "linear" or "constant". 'etype' defines
* which confidence subset is implemented.

* 'significance' is the significance level in decimal form.

********************************************************************************
* Define the initial parameter for the intervention effect.
local erroru = "Upper bound was computed successfully"
local errorl = "Lower bound was computed successfully"
* Select the observed outcome for the treated unit.
matrix Y1 = Y[1..., `treated']
****************************************
* If the treated unit is the first or the last unit of your dataset, you
* must change the next two lines accordingly.
matrix Y0temp = Y[1..., 1..(`treated'-1)]
matrix Y0 = Y0temp,Y[1..., (`treated'+1)...]
****************************************
* Compute the estimated gaps for the treated unit.
matrix wtreated = weights[1...,`treated']
matrix gaps = Y1 - Y0 * wtreated
****************************************
* Compute the initial parameter
if "`etype'" == "constant" {
	local nrows = rowsof(Y1) - `T0'
	matrix U = J(1,`nrows',1)
	matrix gapsafter = gaps[(`T0' + 1)..., 1]
	matrix temp = U * gapsafter/`nrows'
	local param = temp[1,1]
}
else if "`etype'" == "linear" {
	local nrows = rowsof(Y1) - `T0'
	matrix temp = gaps[rowsof(Y1), 1]/`nrows'
	local param = temp[1,1]
}
****************************************
* Initialize the loop parameters
local s = abs(`param')/`param'
local ub = `param'
local lb = `param'
********************************************************************************
* We start our loop whose iterations increase the precision of our
* confidence sets.
forvalues p = 0(1)`precision' {
	local step = (1/2)^`p'
	local reject_l = 0
	local reject_u = 0
	* We start our loop whose interactions try to find the upper bound of our
	* confidence set.
	local attemptu1 = 1
		* The loop has to start by not rejecting the null hypothesis. If it 
		* starts by rejecting the null hypothesis, the confidence set for the
		* requested class of functions is empty. This parameter will help us to
		* guarantee that.
	local attemptl1 = 1
		* The loop has to start by not rejecting the null hypothesis. If it 
		* starts by rejecting the null hypothesis, the confidence set for the
		* requested class of functions is empty. This parameter will help us to
		* guarantee that.
	while `reject_u' == 0 {
		* Create the vector that represents the null hypothesis
		if "`etype'" == "constant" {
			matrix zero = J(`T0', 1, 0)
			matrix eff = J(`nrows', 1, `ub')
			matrix nh = zero\eff
		}
		else if "`etype'" == "linear" {
			matrix zero = J(`T0', 1, 0)
			matrix eff = J(`nrows', 1, .)
			forvalues i = 1(1)`nrows' {
				matrix eff[`i', 1] = `i'*`ub'
			}
			matrix nh = zero\eff
		}
		* Create a matrix to store the test statistics for each region.
		local nregion = colsof(Y)
		matrix ts = J(1,`nregion',.)
		* Estimate the gaps for each placebo unit.
		forvalues j = 1(1)`nregion' {
			* Create the matrices of observed outcomes under the null.
			if `j' == `treated' {
				matrix Y1 = Y[1...,`treated']
				if `treated' == 1 {
					matrix Y0 = Y[1..., (`treated'+1)...]
				}
				else if `treated' == `nregion' {
					matrix Y0 = Y[1..., 1..(`treated'-1)]
				}
				else {
					matrix Y0temp = Y[1..., 1..(`treated'-1)]
					matrix Y0 = Y0temp,Y[1..., (`treated'+1)...]
				}
			}
			else {
				matrix Y1 = Y[1..., `j'] + nh
				if `j' == 1 {
					matrix Y0 = Y[1...,2...]
				}
				else if `j' == `nregion' {
					matrix Y0 = Y[1...,1..`j'-1]
				}
				else {
					matrix Y0temp = Y[1..., 1..(`j'-1)]
					matrix Y0 = Y0temp,Y[1..., (`j'+1)...]
				}
				if `j' < `treated' {
					matrix temp = Y0[1..., `treated' - 1] - nh
					if `treated' == 2 {
						matrix Y0temp1 = Y0[1, 2...]
						matrix Y0 = temp,Y0temp1
					}
					else {
						matrix Y0temp1 = Y0[1...,1..`treated' - 2]
						matrix Y0temp2 = Y0[1..., `treated'...]
						matrix Y0 = Y0temp1,temp,Y0temp2
					}
				}
				else if `j' > `treated' {
					matrix temp = Y0[1..., `treated'] - nh
					if `treated' == `nregion' - 1 {
						matrix Y0temp1 = Y0[1..., 1..`treated' - 1]
						matrix Y0 = Y0temp1,temp
					}
					else {
						matrix Y0temp1 = Y0[1...,1..`treated' - 1]
						matrix Y0temp2 = Y0[1..., `treated' + 1...]
						matrix Y0 = Y0temp1,temp,Y0temp2
					}
				}
			}
			* Estimate the gaps
			matrix gaps = Y1 - Y0 * weights[1..., `j'] - nh
			* Estimate the test statistics (RMSPE)
			matrix temp = (gaps[`T0' + 1..., 1]' * gaps[`T0' + 1..., 1])/`nrows'
			local postp = temp[1,1]
			matrix temp = (gaps[1..`T0', 1]' * gaps[1..`T0', 1])/`T0'
			local prep = temp[1,1]
			matrix ts[1, `j'] = `postp'/`prep'
		}
		********************************
		* Test the null hypothesis
		
		* Create a vector that indicates which units have a RMSPE larger than
		* the RMSPE of the treated unit.
		matrix indicator = J(`nregion',1,0)
		forvalues j = 1(1)`nregion' {
			local tempc = ts[1, `j']
			local tempt = ts[1, `treated']
			if `tempc' >= `tempt' {
				matrix indicator[`j', 1] = 1
			}
		}
		
		* Create a vector with the assignment probabilities.
		matrix temp = J(1, `nregion', .)
		forvalues j = 1(1)`nregion' {
			local temp1 = v[1,`j']
			local temp2 = 2.7182818284590452^(`phi'*`temp1')
			matrix temp[1, `j'] = `temp2'
		}
		matrix U = J(`nregion', 1, 1)
		matrix temp1 = temp * U
		local temp2 = temp1[1,1]
		matrix probv = temp/`temp2'
		
		* Compute the p-value
		matrix temp = probv * indicator
		local pvalue = temp[1,1]
		if `pvalue' <= `significance' {
			local reject_u = 1
		}
		********************************
		* Update the upper bound parameter
		if `reject_u' == 0 {
			local temp = `param' * (`ub'/`param' + `s' * `step')
			local ub = `temp'
			local attemptu1 = 0
		}
		else {
			if `attemptu1' == 1 {
				local erroru = "The confidence set is empty for the requested class of functions."
				display "`erroru'"
				local reject_u = 1
			}
			local temp = `param' * (`ub'/`param' - `s' * `step')
			local ub = `temp'
		}
		if abs(`ub') > 100*abs(`param') {
			local erroru = "Upper bound was not found."
			local reject_u = 1
		}
	}
	****************************************************************************
	* We start our loop whose interactions try to find the lower bound of our
	* confidence set.
	while `reject_l' == 0 {
		* Create the vector that represents the null hypothesis
		if "`etype'" == "constant" {
			matrix zero = J(`T0', 1, 0)
			matrix eff = J(`nrows', 1, `lb')
			matrix nh = zero\eff
		}
		else if "`etype'" == "linear" {
			matrix zero = J(`T0', 1, 0)
			matrix eff = J(`nrows', 1, .)
			forvalues i = 1(1)`nrows' {
				matrix eff[`i', 1] = `i'*`lb'
			}
			matrix nh = zero\eff
		}
		* Create a matrix to store the test statistics for each region.
		local nregion = colsof(Y)
		matrix ts = J(1,`nregion',.)
		* Estimate the gaps for each placebo unit.
		forvalues j = 1(1)`nregion' {
			* Create the matrices of observed outcomes under the null.
			if `j' == `treated' {
				matrix Y1 = Y[1...,`treated']
				if `treated' == 1 {
					matrix Y0 = Y[1..., (`treated'+1)...]
				}
				else if `treated' == `nregion' {
					matrix Y0 = Y[1..., 1..(`treated'-1)]
				}
				else {
					matrix Y0temp = Y[1..., 1..(`treated'-1)]
					matrix Y0 = Y0temp,Y[1..., (`treated'+1)...]
				}
			}
			else {
				matrix Y1 = Y[1..., `j'] + nh
				if `j' == 1 {
					matrix Y0 = Y[1...,2...]
				}
				else if `j' == `nregion' {
					matrix Y0 = Y[1...,1..`j'-1]
				}
				else {
					matrix Y0temp = Y[1..., 1..(`j'-1)]
					matrix Y0 = Y0temp,Y[1..., (`j'+1)...]
				}
				if `j' < `treated' {
					matrix temp = Y0[1..., `treated' - 1] - nh
					if `treated' == 2 {
						matrix Y0temp1 = Y0[1, 2...]
						matrix Y0 = temp,Y0temp1
					}
					else {
						matrix Y0temp1 = Y0[1...,1..`treated' - 2]
						matrix Y0temp2 = Y0[1..., `treated'...]
						matrix Y0 = Y0temp1,temp,Y0temp2
					}
				}
				else if `j' > `treated' {
					matrix temp = Y0[1..., `treated'] - nh
					if `treated' == `nregion' - 1 {
						matrix Y0temp1 = Y0[1..., 1..`treated' - 1]
						matrix Y0 = Y0temp1,temp
					}
					else {
						matrix Y0temp1 = Y0[1...,1..`treated' - 1]
						matrix Y0temp2 = Y0[1..., `treated' + 1...]
						matrix Y0 = Y0temp1,temp,Y0temp2
					}
				}
			}
			* Estimate the gaps
			matrix gaps = Y1 - Y0 * weights[1..., `j'] - nh
			* Estimate the test statistics (RMSPE)
			matrix temp = (gaps[`T0' + 1..., 1]' * gaps[`T0' + 1..., 1])/`nrows'
			local postp = temp[1,1]
			matrix temp = (gaps[1..`T0', 1]' * gaps[1..`T0', 1])/`T0'
			local prep = temp[1,1]
			matrix ts[1, `j'] = `postp'/`prep'
		}
		********************************
		* Test the null hypothesis
		
		* Create a vector that indicates which units have a RMSPE larger than
		* the RMSPE of the treated unit.
		matrix indicator = J(`nregion',1,0)
		forvalues j = 1(1)`nregion' {
			local tempc = ts[1, `j']
			local tempt = ts[1, `treated']
			if `tempc' >= `tempt' {
				matrix indicator[`j', 1] = 1
			}
		}
		
		* Create a vector with the assignment probabilities.
		matrix temp = J(1, `nregion', .)
		forvalues j = 1(1)`nregion' {
			local temp1 = v[1,`j']
			local temp2 = 2.7182818284590452^(`phi'*`temp1')
			matrix temp[1, `j'] = `temp2'
		}
		matrix U = J(`nregion', 1, 1)
		matrix temp1 = temp * U
		local temp2 = temp1[1,1]
		matrix probv = temp/`temp2'
		
		* Compute the p-value
		matrix temp = probv * indicator
		local pvalue = temp[1,1]
		if `pvalue' <= `significance' {
			local reject_l = 1
		}
		********************************
		* Update the lower bound parameter
		if `reject_l' == 0 {
			local temp = `param' * (`lb'/`param' - `s' * `step')
			local lb = `temp'
			local attemptl1 = 0
		}
		else {
			local temp = `param' * (`lb'/`param' + `s' * `step')
			local lb = `temp'
			if `attemptl1' == 1 {
				local erroru = "The confidence set is empty for the requested class of functions."
				display "`errorl'"
				reject_l = 1
			}
		}
		if abs(`lb') > 100*abs(`param') {
			local errorl = "Lower bound was not found."
			local reject_l = 1
		}
	}
}
********************************************************************************
* Define the upper and lower bound of the confidence set
if "`etype'" == "constant" {
	* Upper bound
	matrix zero = J(`T0', 1, 0)
	matrix eff = J(`nrows', 1, `ub')
	matrix upper = zero\eff
	********************
	* Lower bound
	matrix zero = J(`T0', 1, 0)
	matrix eff = J(`nrows', 1, `lb')
	matrix lower = zero\eff
}
else if "`etype'" == "linear" {
* Upper bound
	matrix zero = J(`T0', 1, 0)
	matrix eff = J(`nrows', 1, .)
	forvalues i = 1(1)`nrows' {
		matrix eff[`i', 1] = `i'*`ub'
	}
	matrix upper = zero\eff
	********************
	* Upper bound
	matrix zero = J(`T0', 1, 0)
	matrix eff = J(`nrows', 1, .)
	forvalues i = 1(1)`nrows' {
		matrix eff[`i', 1] = `i'*`lb'
	}
	matrix lower = zero\eff
}
****************************************
* Report the final results
matrix results = upper,lower
matrix colnames results = "upper" "lower"
display "`erroru'"
display "`errorl'"
matrix list results

end
