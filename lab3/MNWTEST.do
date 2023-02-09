/*-----------------------------*/
/* MNW cluster test */
/* implements procedures found in this paper */
/* https://ideas.repec.org/p/qed/wpaper/1428.html */ 
/* written by Matt Webb */
/* version 1.002 */
/* date 04/11/22 */
/*-----------------------------*/

/*program used by main program for MNWTEST*/
capture program drop VECHTEST
program define VECHTEST, rclass

	local y "`1'"
	 
	global XTEMP "`2'"
	 
	global CTEMP "`3'" 

	local null "`4'"
	 
	local alt1 "`5'"
		 
	local first `6'
		
	/*calculate residuals*/
	qui reg `y' $XTEMP $CTEMP	
	qui predict temp_uhat, res
	
	/*for the original data*/
	if (`first'==1) {
	
		/*fitted values and unrestricted residuals*/
		qui predict temp_xbr, xb
		qui gen temp_er = temp_uhat
		
		/*determine number of clusters under alternative*/
		qui egen temp_indexa = group(`alt1') if temp_uhat !=. 
		qui summ temp_indexa 
		global G = r(max)
		
		/*correct for possibility of non-unique naming of fine clusters*/
		/*determine number of clusters under null*/
		qui egen temp_indexn = group(`null') if temp_uhat !=. 
		qui egen temp_ints = group(temp_indexa temp_indexn) 
		qui summ temp_ints
		global H = r(max)
		
		/*put alternative cluster index into mata*/
		qui sort temp_indexa
		qui putmata cg = temp_indexa if temp_uhat!=., omitmissing replace
		
		/*put null cluster index into mata*/
		qui bysort temp_indexa temp_ints: gen temp_count = _n
		qui putmata mg = temp_indexa if temp_count==1 & temp_uhat!=., omitmissing replace
		
		/*values for m factors*/
			global NmK = e(df_r)
			global N = e(N)		
			global K = $N - $NmK
			
	/*creates table of subclusters index by G*/
		mata: infoM = panelsetup(mg,1)
			
	qui putmata ch = temp_ints if temp_uhat!=., omitmissing replace
	
	/*creates tables of obs by null cluster, or by alt cluster*/
	mata: infoH = panelsetup(ch,1)
	mata: infoG = panelsetup(cg, 1)
			
				
		/*calculate scores for each  XVAR*/
		local j = 0
		
		foreach xvar in $XTEMP {
			local j = `j' + 1
			qui reg `xvar' $CTEMP
			cap drop temp_z`j'
			qui predict temp_z`j', res
			global xnum = `j'

		} /*end of xvar*/
				
	/*m-factor under alternative*/
		global M_A = ($G/($G - 1)) * (($N-1)/${NmK})			

	/*m-factor under null*/
		global M_F = ($N-1)/(${NmK})*${H}/(${H}-1)	
		
		/*Gk and Hk matrik*/
		mata: Gk = J(${xnum}^2,${xnum}*(${xnum}+1)/2,0)
			
		forvalues i = 1/$xnum {	
			forvalues j = `i' /$xnum {
				local a = (`j'-1)*${xnum} + `i'			
				local b = (`i' - 1)*${xnum} + `j'			
				local c = (`j'-1)*(${xnum}-`j'/2)+`i'		
				mata: Gk[ `a' ,`c' ] =1
				mata: Gk[ `b',`c'] =1
			} /*end of j*/
		} /*end of i */
			
		mata: Hk = invsym(Gk'*Gk)*Gk'

		/*put G global into mata*/	
			mata: maxg = st_global("G")
			mata: maxg = strtoreal(maxg)
			
		/*put xnum variable into mata*/	
			mata: xnum = st_global("xnum")
			mata: xnum = strtoreal(xnum)
			
	} /*end of "first" if*/
	
	/*calculate scores for each  XVAR*/
	local j = 0
	
	foreach xvar in $XTEMP {
		local j = `j' + 1
		qui gen temp_s`j' = temp_uhat * temp_z`j'
		qui putmata s`j' = temp_s`j' if temp_uhat!=., omitmissing replace
	}
		
	local xnum = $xnum
	
	/*put all scores into one matrix*/
		mata: sh = J($N,$xnum,.)
		forvalues k = 1/$xnum{
				mata: sh[.,`k'] = s`k'
		}
			
		mata{
				
			/*initialize matrices*/
				 temp_sumg = J(xnum,xnum,0) 
				 *temp_sumh = J(xnum,xnum,0) 
				 temp_num_alt = J(xnum,xnum,0)
				 var_right = J(rows(Hk),rows(Hk),0)
				 var_left = J(rows(Hk),rows(Hk),0)
				 ALT = J(xnum ,xnum , 0)
				 NLL = J(xnum , xnum , 0)
				 theta = J(rows(Hk),1,0)
				 		
				/*sum scores by either alt or null clusters*/
				sh_h = panelsum(sh,infoH )
				sg = panelsum(sh,infoG)
				 					
			for (g=1; g<=maxg; g++){
										
				 temp_sumh = J(xnum,xnum,0) 
				 temp_var_left = J(xnum,xnum,0)
				 temp_var_right = J(xnum,xnum,0) 
				 
				 /*alt numerator*/
					temp_sg = sg[g,.]'
					temp_num_alt = temp_sg * temp_sg'
					ALT = ALT + temp_num_alt
					
				/*which obs are in cluster g*/
				strt = infoM[g,1]
				ed = infoM[g,2]

				for (i=strt; i<=ed; i++) {
				
					/*extract relevant row, i,  from score matrix*/
					sh1 = sh_h[i,.]'
					temp_cross = sh1 * sh1'
					
					/*var left*/
						temp_sumh = temp_sumh + temp_cross
					
					/*var right*/
						temp_var_right = Hk * (temp_cross # temp_cross) * Hk'
						var_right = var_right + temp_var_right

				} /*i loop */
				
				/*var left*/
					temp_var_left = Hk * (temp_sumh # temp_sumh) * Hk'
					var_left = var_left + temp_var_left
				
				/*null numerator*/
					temp_sumg = temp_sumg + temp_sumh
					
			} /*g loop*/
			
		/*theta*/
			NLL= temp_sumg
			NLL = ${M_F}:*NLL
			ALT = ${M_A}:*ALT
			theta = vech(ALT-NLL)
		
		/*variance*/
			var_left = 2*var_left
			var_right = 2*var_right
			var = var_left - var_right
		
		/*tau stat*/
		if (xnum == 1)  tau = theta/sqrt(var) 
			
		else if (xnum !=1) tau = theta'* invsym(var)*theta 
		
		} /*mata end*/
				
		return scalar xnum = $xnum
		
		cap drop temp_s*
		cap drop temp_uhat
		
end

/*--------------------------------------------------------*/


/*main program for cluster test*/
/*------------------------------*/
/* program stores the following in memory*/
/* r(P) - analytic P*/
/* r(BP) - bootstrap P*/
/* r(tau) - tau statistic*/
/*------------------------------*/
capture program drop MNWTEST
program define MNWTEST, rclass

	local y "`1'"
	 
	global x "`2'"
	 
	global CTRLVAR "`3'" 

	local null "`4'"
	 
	local alt1 "`5'"
		 	
	global BTEMP `6'
	
	/*estimate the test statistic for sample*/
		qui VECHTEST `y' "${x}" "${CTRLVAR}" `null' `alt1' 1
	
	local xnum = r(xnum)				

	/*store result as local and as tauhat for bootstrap P*/
		mata: tauhat = tau
		mata: st_numscalar("tau", tau)
		mata: st_local("tau", strofreal(tau))
		
	/*display summary statistics and results from analytic test*/
		disp "MNW test"
		disp "null is `null', "  "H is $H"
		disp "alt1 is `alt1', "  "G is $G"
		disp "theta is"
		mata: theta
		disp "tau is `tau' "
		
	/*calculate chi2 dof*/
		local dfchi = `xnum'*(`xnum'+1)/2
		disp "chi df is `dfchi' "
	
	/*calculate and display analytic P value*/
		if (`xnum' == 1) local MNW_P = min(normprob(`tau'),1-normprob(`tau'))
		else if (`xnum' >= 2) local MNW_P = chi2tail(`dfchi', `tau') 		

		disp "P value is `MNW_P' "
		
	if $BTEMP > 0 { 
		
		/*intitialize variables and bootstrap matrix*/
		qui gen temp_uni = .
		qui gen temp_pos = .
		qui gen temp_ernew = .
		qui gen temp_ywild = .
		mata: tau_star = J($BTEMP,1,.)
		
		sort `null'
		forvalues b = 1/$BTEMP {
				/*create bootstrap y - using wild cluster bootstrap*/
				qui by `null': replace temp_uni = runiform()	
				qui by `null': replace temp_pos = temp_uni[1]<.5   /*cluster level rademacher indicator*/
				qui replace temp_ernew = (2*temp_pos-1)*temp_er  /*transformed residuals */
				qui replace temp_ywild = temp_xbr + temp_ernew 
				
				/*calculate tau using bootstrap y*/
				qui VECHTEST temp_ywild "${x}" "${CTRLVAR}" `null' `alt1'  2
				mata: tau_star[`b',1]=tau
		
		}
			
		/*calculate and display bootstrap P value*/
		
			mata {
			
				if (`xnum'==1) temp_rej =  abs(tauhat[1,1]):<=abs(tau_star)
				else if (`xnum'>=2) temp_rej =  tauhat[1,1]:<=tau_star 
				
				temp_U = J(rows(temp_rej),1,1) 
				temp_sum = temp_U'*temp_rej
				boot_p = temp_sum / rows(temp_rej)
				
				st_numscalar("boot_p", boot_p)
				st_local("boot_p", strofreal(boot_p))
			}
			
			disp "Bootstrap P value `boot_p' "
			return scalar BP = `boot_p'
				
			cap drop temp_uni temp_pos temp_ernew temp_ywild

	} /*end of bootstrap "if" code*/
	
		return scalar tau = `tau'		
		return scalar xnum = $xnum
		return scalar P = `MNW_P'
		
		cap drop temp_ints temp_z* temp_er  temp_xbr temp_indexn temp_indexa temp_count
	
end  /*end of program*/	

/*--------------------------------------------------------*/

/*IM Test program, not used within MNWTEST*/
capture program drop IMTEST
program define IMTEST, rclass 

	local y "`1'"
	 
	local x "`2'"
	 
	global CTRLVAR "`3'" 

	local null "`4'"
	 
	local alt "`5'"
		 
/*determine number of clusters under alternative*/ 
	qui reg `y' `x' $CTRLVAR
	qui predict temp_uhat, resid
	
	qui egen temp_grp = group(`alt') if temp_uhat!=.
	qui summ temp_grp
	
/*store number of clusters as j in mata*/
	local j = r(max)
	mata: j = `j'
	
/*create empty matrices to store the j beta and s.e.*/
	mata: beta = J(j,1,.)
	mata: omega = J(j,1,.)
	
/*calculate the beta and s.e. per coarse cluster*/
	/*cluster s.e. under the null*/
	forvalues g=1/`j' {
		
		qui reg `y' `x' $CTRLVAR if temp_grp==`g', cluster(`null')
		
		/*store the beta and s.e. estimates in the resepective matrix*/
		local beta = _b[`x']
		mata: beta[`g',1]=`beta'
		local se = _se[`x']
		mata: omega[`g',1]=`se'
	}
	
/*Calculate the IM (2012) standard error*/	
	/*it is just the variance of the j beta */
	mata: S2 = variance(beta)
	/*display S2*/
	disp "IM TEST, variance \ P value"
	mata: S2
	
/*matrix for all the yj estimates*/	
mata: ybar = J(9999,1,.)

/*replication loop*/
forvalues k=1/9999{
	
	/*multiply the standard errors by standard normals*/
	mata: yj = omega:*invnormal(uniform(j,1))
	
	/*calculate the average of the yj*/
	mata avey = mean(yj[.,1])
	
	/*square the mean differences*/
	mata: sy2 = (yj:-avey):*(yj:-avey)
	
	/*sum the squares, divide by (j-1) */
	mata: ybk = (1/(j-1))* colsum(sy2[.,1])

	/*store in the matrix*/
	mata: ybar[`k',1]=ybk
}

/*calculate the p-value*/
	mata: temp_rej =  S2[1,1]:<ybar
	mata: temp_U = J(rows(temp_rej),1,1) 
	mata: temp_sum = temp_U'*temp_rej
	mata: IM_p = temp_sum / rows(temp_rej)

	/*display P value*/
	mata: IM_p
	
	cap drop temp_grp
	cap drop temp_uhat

end


/*--------------------------------------*/
/*change log*/
/*--------------------------------------*/

*1.002 change to one sided p-values 
