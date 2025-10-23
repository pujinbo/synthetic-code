*****************************************************************************************
* Replication file for Cattaneo, Feng, Palomba, and Titiunik (2022)
*****************************************************************************************

*****************************************************************************************
** net install scpi, from(https://raw.githubusercontent.com/nppackages/scpi/master/stata) replace
*****************************************************************************************

cd "/Users/fpalomba/Dropbox/projects/scpi/CFPT_2022_application/code/stata/"

* Load dataset
use "final_data.dta", clear

keep if restricted2 == 1
gen lgdp = log(rgdppp)

global continents `" "Asia" "Africa" "Europe" "North America" "South America" "'
global sims = 200

foreach cont of global continents {

	qui levelsof countryname if continent == "`cont'" & treated == 1 & trDate <= 1992 & countryname != "Slovenia"
	global aux = r(levels)
	global units_cont
	local i = 1
	foreach unit of global aux {
		if `i' == 1 {
			global units_cont "`unit'"
		} 
		else {
			global units_cont "$units_cont, `unit'" 		
		}
		local ++i
	}


	*****************************************************
	* unit-time treatment effect
	*****************************************************
	scdatamulti lgdp lsc, dfname("python_scdatamulti") id(countryname) outcome(lgdp)       ///
				time(year) treatment(liberalization) cointegrated("True") constant("True") ///
				covadj("constant, trend") anticipation("1") post_est(10)                   ///
				units_est($units_cont) effect("unit-time")

	scpi, dfname("python_scdatamulti") name("L1-L2") sims($sims) 
	scplotmulti, uncertainty("gaussian") ptype("series")
	graph export "unit_time_`cont'.png", replace


	*****************************************************
	* average unit post-treatment effect
	*****************************************************
	scdatamulti lgdp lsc, dfname("python_scdatamulti") id(countryname) outcome(lgdp)       ///
				time(year) treatment(liberalization) cointegrated("True") constant("True") ///
				covadj("constant, trend") anticipation("1") post_est(10)                   ///
				units_est($units_cont) effect("unit")

	scpi, dfname("python_scdatamulti") name("L1-L2") sims($sims) 
	scplotmulti, uncertainty("gaussian") ptype("series")
	graph export "unit_`cont'.png", replace


	*****************************************************
	* average post-treatment effect on the treated
	*****************************************************
	scdatamulti lgdp lsc, dfname("python_scdatamulti") id(countryname) outcome(lgdp)       ///
				time(year) treatment(liberalization) cointegrated("True") constant("True") ///
				covadj("constant, trend") anticipation("1") post_est(10)                   ///
				units_est($units_cont) effect("time")

	scpi, dfname("python_scdatamulti") name("L1-L2") sims($sims) 
	scplotmulti, uncertainty("gaussian") ptype("series")
	graph export "time_`cont'.png", replace

}

