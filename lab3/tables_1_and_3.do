
/* 2022-04-04

James G. MacKinnon, Morten O. Nielsen, and Matthew D. Webb,
"Cluster-robust inference: A guide to empirical practice", QED Working
Paper No. 1456, 2022.

This file generates all the CV3 results in Table 1 and the summary
statistics in Table 3 for the employment regression. */

*cd "C:\Users\mattw\Dropbox\research\JOE_survey\min_wage\data_csv"

*adopath + "C:\Users\mattw\Dropbox\research\influence\summcluster"

set more off
set seed 98789

clear all

cap log close
log using "tables_1_and_3.txt", text replace

timer clear 1
timer on 1 

global CVARS "styear state region year"   

clear 
insheet using "min_wage_teen_hours2.csv"


/*mata matrices for the tables*/

mata: tabone = J(3,24,.)
mata: tabthree = J(8,8,.)

/*counters*/

	*j for columns
	*i for rows
	
	local j = 1
	local i = 0

/*hours*/

    
/*one way CV1 and CV3 */
summclust mw, yvar(hours2) xvar(black female) fevar(educ age year state) ///
 cluster(styear) rho(0.5) jack
 
	mata:tabone[`j',4] =cvstuff[1,3]
	mata:tabone[`j',5] =cvstuff[1,4]
	
	mata:tabone[`j',7] =cvstuff[2,3]
	mata:tabone[`j',8] =cvstuff[2,4]
	
	mata:tabthree[1,1] = $G
	mata:tabthree[1,2] = gstarzero
	mata:tabthree[1,3] = clustsum[4,1]
	mata:tabthree[1,4] = clustsum[1,1]
	mata:tabthree[1,5] = clustsum[2,1]
	mata:tabthree[1,6] = clustsum[3,1]
	mata:tabthree[1,7] = clustsum[5,1]
	mata:tabthree[1,8] = clustsum[6,1]
	
 
 summclust mw, yvar(hours2) xvar(black female) fevar(educ age year state) ///
 cluster(state) rho(0.5) jack
 
	mata:tabone[`j',9] =cvstuff[1,3]
	mata:tabone[`j',10] =cvstuff[1,4]
	
	mata:tabone[`j',12] =cvstuff[2,3]
	mata:tabone[`j',13] =cvstuff[2,4]
	
	mata:tabthree[2,1] = $G
	mata:tabthree[2,2] = gstarzero
	mata:tabthree[2,3] = clustsum[4,1]
	mata:tabthree[2,4] = clustsum[1,1]
	mata:tabthree[2,5] = clustsum[2,1]
	mata:tabthree[2,6] = clustsum[3,1]
	mata:tabthree[2,7] = clustsum[5,1]
	mata:tabthree[2,8] = clustsum[6,1]
 
 summclust mw, yvar(hours2) xvar(black female) fevar(educ age year state) ///
 cluster(region) rho(0.5) jack
 
	mata:tabone[`j',14] =cvstuff[1,3]
	mata:tabone[`j',15] =cvstuff[1,4]
	
	mata:tabone[`j',17] =cvstuff[2,3]
	mata:tabone[`j',18] =cvstuff[2,4]
	
	mata:tabthree[3,1] = $G
	mata:tabthree[3,2] = gstarzero
	mata:tabthree[3,3] = clustsum[4,1]
	mata:tabthree[3,4] = clustsum[1,1]
	mata:tabthree[3,5] = clustsum[2,1]
	mata:tabthree[3,6] = clustsum[3,1]
	mata:tabthree[3,7] = clustsum[5,1]
	mata:tabthree[3,8] = clustsum[6,1]
 
 summclust mw, yvar(hours2) xvar(black female) fevar(educ age year state) ///
 cluster(year) rho(0.5) jack
 
	mata:tabthree[4,1] = $G
	mata:tabthree[4,2] = gstarzero
	mata:tabthree[4,3] = clustsum[4,1]
	mata:tabthree[4,4] = clustsum[1,1]
	mata:tabthree[4,5] = clustsum[2,1]
	mata:tabthree[4,6] = clustsum[3,1]
	mata:tabthree[4,7] = clustsum[5,1]
	mata:tabthree[4,8] = clustsum[6,1]
 


/*multiway*/
	qui reg hours2 mw black female i.educ i.age i.year i.state,  robust
		
		local beta = _b[mw]
		mata:tabone[`j',1] = `beta'
	
	test mw = 0
		local t = sqrt(r(F))
		local p = r(p)
		mata:tabone[`j',2] = `t'
		mata:tabone[`j',3] = `p'
	
	boottest mw, cluster(styear) noci reps(99999)
	
		local p = r(p)
		mata:tabone[`j',6] = `p'
	
	boottest mw, cluster(statefip) noci reps(99999)
	
		local p = r(p)
		mata:tabone[`j',11] = `p'
	
	boottest mw, cluster(region) noci weight(webb) reps(99999)
	
		local p = r(p)
		mata:tabone[`j',16] = `p'
	
	waldtest mw, cluster(state year) noci
	
		local t = r(t)
		local p = r(p)
		mata:tabone[`j',19] = `t'
		mata:tabone[`j',20] = `p'
		
	
	boottest mw, cluster(state year) bootcluster(year) noci  reps(99999)
	
		local p = r(p)
		mata:tabone[`j',21] = `p'

	waldtest mw, cluster(region year) noci
	
		local t = r(t)
		local p = r(p)
		mata:tabone[`j',22] = `t'
		mata:tabone[`j',23] = `p'
	
	
	boottest mw, cluster(region year) bootcluster(region) noci reps(99999)
	
		local p = r(p)
		mata:tabone[`j',24] = `p'

clear 
insheet using min_wage_teen_empop.csv

global YVARS "empop student "



/*empop and student*/
foreach yvar in $YVARS{

		local j = `j' + 1
		/*one way CV1 and CV3 */
	summclust mw, yvar(`yvar') xvar(black female) fevar(educ age year state) ///
	 cluster(styear) rho(0.5) jack
	 
		mata:tabone[`j',4] =cvstuff[1,3]
		mata:tabone[`j',5] =cvstuff[1,4]
		
		mata:tabone[`j',7] =cvstuff[2,3]
		mata:tabone[`j',8] =cvstuff[2,4]
		
		mata:tabthree[5,1] = $G
		mata:tabthree[5,2] = gstarzero
		mata:tabthree[5,3] = clustsum[4,1]
		mata:tabthree[5,4] = clustsum[1,1]
		mata:tabthree[5,5] = clustsum[2,1]
		mata:tabthree[5,6] = clustsum[3,1]
		mata:tabthree[5,7] = clustsum[5,1]
		mata:tabthree[5,8] = clustsum[6,1]
		
	 
	 summclust mw, yvar(`yvar') xvar(black female) fevar(educ age year state) ///
	 cluster(state) rho(0.5) jack
	 
		mata:tabone[`j',9] =cvstuff[1,3]
		mata:tabone[`j',10] =cvstuff[1,4]
		
		mata:tabone[`j',12] =cvstuff[2,3]
		mata:tabone[`j',13] =cvstuff[2,4]
		
		mata:tabthree[6,1] = $G
		mata:tabthree[6,2] = gstarzero
		mata:tabthree[6,3] = clustsum[4,1]
		mata:tabthree[6,4] = clustsum[1,1]
		mata:tabthree[6,5] = clustsum[2,1]
		mata:tabthree[6,6] = clustsum[3,1]
		mata:tabthree[6,7] = clustsum[5,1]
		mata:tabthree[6,8] = clustsum[6,1]
	 
	 summclust mw, yvar(`yvar') xvar(black female) fevar(educ age year state) ///
	 cluster(region) rho(0.5) jack
	 
		mata:tabone[`j',14] =cvstuff[1,3]
		mata:tabone[`j',15] =cvstuff[1,4]
		
		mata:tabone[`j',17] =cvstuff[2,3]
		mata:tabone[`j',18] =cvstuff[2,4]
		
		mata:tabthree[7,1] = $G
		mata:tabthree[7,2] = gstarzero
		mata:tabthree[7,3] = clustsum[4,1]
		mata:tabthree[7,4] = clustsum[1,1]
		mata:tabthree[7,5] = clustsum[2,1]
		mata:tabthree[7,6] = clustsum[3,1]
		mata:tabthree[7,7] = clustsum[5,1]
		mata:tabthree[7,8] = clustsum[6,1]
	 
	 summclust mw, yvar(`yvar') xvar(black female) fevar(educ age year state) ///
	 cluster(year) rho(0.5) jack
	 
		mata:tabthree[8,1] = $G
		mata:tabthree[8,2] = gstarzero
		mata:tabthree[8,3] = clustsum[4,1]
		mata:tabthree[8,4] = clustsum[1,1]
		mata:tabthree[8,5] = clustsum[2,1]
		mata:tabthree[8,6] = clustsum[3,1]
		mata:tabthree[8,7] = clustsum[5,1]
		mata:tabthree[8,8] = clustsum[6,1]
	 


	/*multiway*/
		qui reg `yvar' mw black female i.educ i.age i.year i.state,  robust
			
			local beta = _b[mw]
			mata:tabone[`j',1] = `beta'
		
		test mw = 0
			local t = sqrt(r(F))
			local p = r(p)
			mata:tabone[`j',2] = `t'
			mata:tabone[`j',3] = `p'
		
		boottest mw, cluster(styear) noci reps(99999)
		
			local p = r(p)
			mata:tabone[`j',6] = `p'
		
		boottest mw, cluster(statefip) noci reps(99999)
		
			local p = r(p)
			mata:tabone[`j',11] = `p'
		
		boottest mw, cluster(region) noci weight(webb) reps(99999)
		
			local p = r(p)
			mata:tabone[`j',16] = `p'
		
		waldtest mw, cluster(state year) noci
		
			local t = r(t)
			local p = r(p)
			mata:tabone[`j',19] = `t'
			mata:tabone[`j',20] = `p'
			
		
		boottest mw, cluster(state year) bootcluster(year) noci  reps(99999)
		
			local p = r(p)
			mata:tabone[`j',21] = `p'

		waldtest mw, cluster(region year) noci
		
			local t = r(t)
			local p = r(p)
			mata:tabone[`j',22] = `t'
			mata:tabone[`j',23] = `p'
		
		
		boottest mw, cluster(region year) bootcluster(region) noci reps(99999)
		
			local p = r(p)
			mata:tabone[`j',24] = `p'

}

mata: tabone'
mata: tabthree


timer off 1
timer list 1

cap log close
 
exit

