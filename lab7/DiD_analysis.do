
use "Organized Crime and Political Quality.dta", clear

********************************************************************************
***** TWFE *****
********************************************************************************
** Baseline result in Table 1 column 1

* Table 1
* Results with pre-dissolution period always going back to 1985

* Original specification (xtreg)
//set matsize 10000
//xtset ID_municip year
//xi: xtreg MeanEduPol i.year i.ID_municip*trend mafiaben befcomgeneral, fe robust 

* Estimation with reghdfe (much faster)
reghdfe MeanEduPol befcomgeneral mafiaben, a(ID_municip year i.ID_municip#c.trend) cluster(ID_municip)

* Recoding of the treatment variable (T=1 after the dissolution)
gen befcomgeneral_recode = (1-befcomgeneral)*mafiaben

* Substantively equivalent results
reghdfe MeanEduPol befcomgeneral_recode mafiaben, a(ID_municip year i.ID_municip#c.trend) cluster(ID_municip)



********************************************************************************
***** De Chaisemartin & D'Haultfoeuille *****
********************************************************************************
// ssc install twowayfeweights, replace
// ssc install did_multiplegt, replace

** Create the group variable
bys desc_comune (year): gen change = (befcomgeneral_recode==1 & befcomgeneral_recode[_n-1]==0) 
by desc_comune: egen group = min(cond(change, year, .)) 
replace group = 0 if mafiaben==0



***** Negative weights diagnostics *****

** How many weights are negative?
** How much TE heterogeneity is needed to have beta_fe with the opposite sign of ATT?
twowayfeweights MeanEduPol ID_munic year befcomgeneral_recode, type(feTR) 


***** DC&DH estimator *****

** To speed up computation, we change the group from municipality to cohort

** Baseline call
did_multiplegt MeanEduPol group year befcomgeneral_recode, breps(100) cluster(ID_municip) count_switchers_tot

** Placebo test for pre-trends
did_multiplegt MeanEduPol group year befcomgeneral_recode, placebo(5) breps(100) cluster(ID_municip) count_switchers_tot

** Dynamic effect
did_multiplegt MeanEduPol group year befcomgeneral_recode, robust_dynamic dynamic(5) breps(100) cluster(ID_municip) count_switchers_tot


