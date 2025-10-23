
clear

cd "...."

clear
import delimited "employment_peradded_nonstationary_FALSE_results.csv", clear 



egen pre = group(t0)
egen Rank_f1 = group(fac1)
egen Rank_fe = group(fe)




table fac3 t0, c(mean rec_sc_demean_fac3)



mat R = J(30,7,.)

****** FE


forvalues pre =1(1)4 {

local row=`pre'

summ sc_mean_eff_no_brk if pre==`pre' & Rank_fe==50 & t1==12
mat R[`row',1] = r(mean)

summ sc_demean_mean_eff_no_brk if pre==`pre' & Rank_fe==50 & t1==12
mat R[`row',2] = r(mean)

summ did_mean_eff_no_brk if pre==`pre' & Rank_fe==50 & t1==12
mat R[`row',3] = r(mean)

summ sc_sd_eff_no_brk if pre==`pre'  & Rank_fe==50 & t1==12
mat R[`row',5] = r(mean)

summ sc_demean_sd_eff_no_brk if pre==`pre' & Rank_fe==50 & t1==12
mat R[`row',6] = r(mean)

summ did_sd_eff_no_brk if pre==`pre' & Rank_fe==50 & t1==12
mat R[`row',7] = r(mean)


}



forvalues pre =1(1)4 {

local row=`pre'+5

summ sc_mean_eff_no_brk if pre==`pre' & Rank_fe==2 & t1==12
mat R[`row',1] = r(mean)

summ sc_demean_mean_eff_no_brk if pre==`pre' & Rank_fe==2 & t1==12
mat R[`row',2] = r(mean)

summ did_mean_eff_no_brk if pre==`pre' & Rank_fe==2 & t1==12
mat R[`row',3] = r(mean)

summ sc_sd_eff_no_brk if pre==`pre'  & Rank_fe==2 & t1==12
mat R[`row',5] = r(mean)

summ sc_demean_sd_eff_no_brk if pre==`pre' & Rank_fe==2 & t1==12
mat R[`row',6] = r(mean)

summ did_sd_eff_no_brk if pre==`pre' & Rank_fe==2 & t1==12
mat R[`row',7] = r(mean)


}




******* fac


forvalues pre =1(1)4 {

local row=`pre'+10

summ sc_mean_eff_brk_fac1 if pre==`pre' & Rank_f1==50 & t1==12
mat R[`row',1] = r(mean)

summ sc_demean_mean_eff_brk_fac1 if pre==`pre' & Rank_f1==50 & t1==12
mat R[`row',2] = r(mean)

summ did_mean_eff_brk_fac1 if pre==`pre' & Rank_f1==50 & t1==12
mat R[`row',3] = r(mean)

summ sc_sd_eff_brk_fac1 if pre==`pre'  & Rank_f1==50 & t1==12
mat R[`row',5] = r(mean)

summ sc_demean_sd_eff_brk_fac1 if pre==`pre' & Rank_f1==50 & t1==12
mat R[`row',6] = r(mean)

summ did_sd_eff_brk_fac1 if pre==`pre' & Rank_f1==50 & t1==12
mat R[`row',7] = r(mean)


}





forvalues pre =1(1)15 {

local row=`pre'+15

summ sc_mean_eff_brk_fac1 if pre==`pre' & Rank_f1==2 & t1==12
mat R[`row',1] = r(mean)

summ sc_demean_mean_eff_brk_fac1 if pre==`pre' & Rank_f1==2 & t1==12
mat R[`row',2] = r(mean)

summ did_mean_eff_brk_fac1 if pre==`pre' & Rank_f1==2 & t1==12
mat R[`row',3] = r(mean)

summ sc_sd_eff_brk_fac1 if pre==`pre'  & Rank_f1==2 & t1==12
mat R[`row',5] = r(mean)

summ sc_demean_sd_eff_brk_fac1 if pre==`pre' & Rank_f1==2 & t1==12
mat R[`row',6] = r(mean)

summ did_sd_eff_brk_fac1 if pre==`pre' & Rank_f1==2 & t1==12
mat R[`row',7] = r(mean)


}



*Table 1
mat li R


set scheme s2mono


la var sc_mean_eff_no_brk "Bias of original SC estimator" 
la var fe  "Fixed Effect" 


twoway (scatter  sc_mean_eff_no_brk fe if t1==12 & t0==120,  ms(d)) (scatter   sc_mean_eff_no_brk fe if t1==12 & t0==1200,  ms(0)), yscale(range(-1 1)) xlabel(-1.5(0.5)3.5)  ylabel(-1(0.2)1) yline(0) legend(label(1 "T{subscript:0} = 120") label(2 "T{subscript:0} = 1200") ) graphregion(color(white)) 

*Figure 1.A
graph export "Fig1A.pdf", replace


la var sc_demean_mean_eff_brk_fac1 "Bias of demeaned SC estimator" 
la var fac1  "Factor loading associated to 1st common factor" 


twoway (scatter  sc_demean_mean_eff_brk_fac1 fac1 if t1==12 & t0==120,  ms(d)) (scatter  sc_demean_mean_eff_brk_fac1 fac1 if t1==12 & t0==1200,  ms(0))   , yscale(range(-2.6 2)) xlabel(-3(0.5)2.2) ylabel(-2.5(0.5)2) yline(0) legend(label(1 "T{subscript:0} = 120") label(2 "T{subscript:0} = 1200") ) graphregion(color(white)) 

*Figure 1.B
graph export "Fig1B.pdf", replace


*Table A1
table fac1 t0, c(mean rej_05_q1_no_brk)
table fac1 t0, c(mean rej_05_q1_brk_fac1)







