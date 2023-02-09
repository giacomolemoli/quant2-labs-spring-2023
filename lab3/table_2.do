/*number of bootstraps*/
	
	set seed 98789
	global B = 2

/*load the IM program and MNW programs into memory*/

	cap log close
	log using "table2.txt", text replace
	
/*set up matrix to create table*/
	mata: results = J(18,3,.)
	
/*load data*/
clear
insheet using "min_wage_teen_hours2.csv"

gen newid = _n

global XVARS black female i.educ i.age i.year i.state

/*----------------------------------------------*/
/*MNWTEST arguments*/

/* 1 - single outcome variable                  */
/* 2 - variable(s) of interest   (as global)    */
/* 3 - control variables (as global)            */
/* 4 - null clustering variable                 */
/* 5 - alternative clustering variable          */
/* 6 - number of bootstraps - enter 0 for none  */
/*----------------------------------------------*/

/*set up locals for the null and alt*/

		global null1 "newid"
		global null2 "newid"
		global null3 "newid"
		global null4 "styear"
		global null5 "styear"
		global null6 "statefip"
		
		global alt1 "styear"
		global alt2 "statefip"
		global alt3 "region"
		global alt4 "statefip"
		global alt5 "region"
		global alt6 "region"
	
/*loop over the six tests*/	
	
	forvalues i=1/6 {
		
		MNWTEST hours2 mw "$XVARS"  ${null`i'} ${alt`i'} $B
			local tau = r(tau)
			local pval = r(P)
			local bspval = r(BP)
			mata:results[3*(`i'-1) + 1,1] = `tau'
			mata:results[3*(`i'-1) + 2,1] = `pval'
			mata:results[3*(`i'-1) + 3,1] = `bspval'
	}
	
	
/*student and empop*/	

clear
insheet using "min_wage_teen_empop.csv"

gen newid = _n

global YVARS "empop student"

/*loop over the six tests*/	
	
	forvalues i=1/6 {
		
		/*loop over the two outcome variables*/
		local j = 1
		foreach yvar in $YVARS {
			
			local j = `j' + 1
			MNWTEST `yvar' mw "$XVARS"  ${null`i'} ${alt`i'} $B
			local tau = r(tau)
			local pval = r(P)
			local bspval = r(BP)
			
			mata:results[3*(`i'-1) + 1,`j'] = `tau'
			mata:results[3*(`i'-1) + 2,`j'] = `pval'
			mata:results[3*(`i'-1) + 3,`j'] = `bspval'
		}
		

	}
	
	mata: results2 = results[1::9,.],results[10::18,.]
	mata: results2


	cap log close