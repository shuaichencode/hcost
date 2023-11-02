*! version 1.4.2  14Nov2019
/*
HCOST calculates mean survival and costs, and the ICER with CIs, with censored data.

1. option for `method': 0 = BT estimator when cost history is not given, 1 = ZT estimator with cost history
2. option for time limit `l' (in days) is required 
3. For BT estimator without start and stop, assume the cost is accumulated before time limit L.
4. option for annual discounting rates of costs (with cost history available) and survival (`drcost' and `drsurv').
5. construct CI of ICER by Fieller method
6. Data must be stset before hcost is run
7. discount the cost according to its start date. Divide the interval into (user-specified) k-day subintervals,
 if the interval is larger than k days. 
*/

/*
Update log:
30Aug2018: add additional l[n+1]=orgs[n+1]=orgd[n+1]=orgr[n+1]=orgy[n+1]=0 to avoid error of "subscript invalid" 
16Jul2019: optimize program to increase computational speed: 
	 1. calculate KM estimator in advance and saving values for later use; 
	 2. do not run additional computation if BT method is used
17Jul2019: check potential data errors:
	1. check whether idvar is numeric; 
	2. check whether surv and delta are constant within id; 
	3. check whether same id appears in different groups
14Nov2019: change the way of truncating costs:
	1. Previous version only truncates cost if start and stop time are provided. If start and stop time are not 
		provided, we assume costs are evenly spread within 0 to follow-up time, and prorate costs if follow-up time 
		is truncated by l().
	2. Check data error: stop time is after the death time
*/


program hcost, eclass sortpreserve 

	version 11.2
	
	
	syntax varlist(numeric min=2 max=2) [if] [in] , l(real) [start(varname numeric) ///
		stop(varname numeric)  drsurv(real 0) drcost(real 0) icer(integer 0) ///
		k(integer 90) method(integer 1) group(varname numeric) level(real 95) *]
	
	st_is 2 analysis
	gettoken id cost : varlist
	marksample touse
	local surv  _t
	local delta  _d
	

	if (`l' <= 0) {
		di in red "L must be positive"
		exit 198
	}
	
	if (`k' <= 0) {
		di in red "k must be positive"
		exit 198
	}
	
	if (`drsurv' > 1 | `drsurv' < 0) {
		di in red "drsurv must be between 0 and 1"
		exit 198
	}
        
    if (`drcost' > 1 | `drcost' < 0) {
		di in red "drcost must be between 0 and 1"
		exit 198
	}
	
	if ((`method'!=0) & (`method'!=1)){
		di in red "method must be 0 or 1"
		exit 198
	}
	
	if ((`icer'!=0) & (`icer'!=1)){
		di in red "icer must be 0 or 1"
		exit 198
	}
	
	if ((`icer'==1) & ("`group'"=="")){
		di in red "group is required for calculating icer"
		exit 198
	}
	
	if ((`method'==1) & (("`start'"=="")|("`stop'"==""))){
		di in red "start and stop are required for method 1"
		exit 198
	}
	
	if ((`drcost'>0) & (("`start'"=="")|("`stop'"==""))){
		di in red "start and stop are required for discounting cost"
		exit 198
	}

	if(("`start'"!="")&("`stop'"!="")){
		if (`stop'>`surv' & `delta'==1){ 
			di in red "a stop time is larger than the death time"
			exit 459
		}
	}
		
	tempvar ix
	qui gen `ix' = .
	
	tempname b V b1 V1 b2 V2
	
	if "`group'" == "" {
	
		tempvar maxs 
		egen `maxs'=max(`surv')
		if(`l'>`maxs'){
			di in red "Warning: The time limit L is too large (greater than the last observation time)."
			exit 198
		}
	
		preserve
		qui keep if `touse'
		
		sort `id'
		qui by `id': replace `ix' = _n == 1 
		qui replace `ix' = `ix' & `touse'
		
		capture by `id':  /*
			*/ assert `surv'==`surv'[1] if `touse'
		if _rc { 
			di in red "Survival time not constant within `id'. Check survival time and precision of `id'."
			exit 459
		}
	
		capture by `id':  /*
			*/ assert `delta'==`delta'[1] if `touse'
		if _rc { 
			di in red "Death indicator not constant within `id'. Check indicator and precision of `id'."
			exit 459
		}
		

		
	
		if(("`start'"!="")&("`stop'"!="")){
			mata: CalCostSurv("`touse'","`ix'","`id'","`surv'","`delta'",	///
			"`start'","`stop'","`cost'", `drsurv', `drcost', `k', `l',`method')
		}
		else{
			mata: CalCostSurv0("`touse'","`ix'","`id'","`surv'","`delta'",	///
			"`cost'", `drsurv', `l')
		}
		
		mat `b' = J(1,2,0)
		mat `V' = J(2,2,0)

		mat `b'[1,1] = e(mean_cost)
		mat `b'[1,2] = e(mean_time)
		mat `V'[1,1] = e(v_cost)
		mat `V'[2,2] = e(v_time)
		mat `V'[1,2] = e(cov)
		mat `V'[2,1] = e(cov)
		
		local n=e(n_id)
		local nobs=e(n_obs)		
		local pcensor: di %8.2f e(perc)
		
		local names cost survival 

		mat colnames `b' = `names'
		mat rownames `V' = `names'
		mat colnames `V' = `names'
		ereturn post `b' `V', esample(`touse') 
		
		ereturn scalar n= `n'
		ereturn scalar nobs = `nobs'
		ereturn scalar pcensor = `pcensor'
		
		if (e(v_cost)<0){
			di in red "Warning: the estimated variance of the mean costs is negative. The time limit L may be too large."
		}
	
		if _caller() < 12 {
			ereturn local cmd "hcost"
			_table_ci, level(`level')			
		}			
		else {
			_coef_table, level(`level')
		}

		if (`method'==1){
			di "Method used: ZT (using cost history)"
		}
		else {
			di "Method used: BT (using total costs)"
		}
		display "Annual discounting rate for costs: `drcost'"
		display "Annual discounting rate for survival time: `drsurv'"	
		di "Time limit L (in days): `l'"

		
		restore
	}

	if "`group'" != "" {
		
		tempvar group2
		qui egen `group2' = group(`group')
		qui levelsof `group2', local(levels)
		local wdcnt : word count `levels'
		if (`wdcnt' > 2) {
			di in red "group can have at most two levels"
			exit 198
		} 
		
		tempfile one two
		
		save `one'
		qui keep if `touse'
		
		sort `id' 
		qui by `id': replace `ix' = _n == 1 
		qui replace `ix' = `ix' & `touse' 
			
		capture by `id':  /*
			*/ assert `group'==`group'[1] if `touse'
		if _rc { 
			di in red "same `id' appears in different groups"
			exit 459
		}
				
		capture by `id':  /*
			*/ assert `surv'==`surv'[1] if `touse'
		if _rc { 
			di in red "Survival time not constant within `id'. Check survival time and precision of `id'."
			exit 459
		}
	
		capture by `id':  /*
			*/ assert `delta'==`delta'[1] if `touse'
		if _rc { 
			di in red "Death indicator not constant within `id'. Check indicator and precision of `id'."
			exit 459
		}
		
		save `two'
		
		foreach lev of local levels {	
			
			qui keep if `touse' & `group2' == `lev'
			
			tempvar maxs 
			egen `maxs'=max(`surv')
			if(`l'>`maxs'){
				di in red "Warning: The time limit L is too large (greater than the last observation time)."
				exit 198
			}
	
			if(("`start'"!="")&("`stop'"!="")){
				mata: CalCostSurv("`touse'","`ix'","`id'","`surv'","`delta'",	///
				"`start'","`stop'","`cost'", `drsurv', `drcost', `k', `l',`method')
			}
			else{
				mata: CalCostSurv0("`touse'","`ix'","`id'","`surv'","`delta'",	///
				"`cost'", `drsurv', `l')
			}
			
			
			local mean_c`lev' = e(mean_cost)
			local mean_t`lev' = e(mean_time)
			local v_c`lev' = e(v_cost)
			local v_t`lev' = e(v_time)
			local cov`lev' = e(cov)
			
			mat `b' = J(1,2,0)
			mat `V' = J(2,2,0)
			
			mat `b'[1,1] = e(mean_cost)
			mat `b'[1,2] = e(mean_time)
			mat `V'[1,1] = e(v_cost)
			mat `V'[2,2] = e(v_time)
			mat `V'[2,1] = e(cov)
			mat `V'[1,2] = e(cov)

			local names cost survival 

			mat colnames `b' = `names'
			mat rownames `V' = `names'
			mat colnames `V' = `names'

			mat b`lev'=`b'
			mat V`lev'=`V'
			
			local n`lev'=e(n_id)
			local nobs`lev'=e(n_obs)
			local pcensor`lev': di %8.2f e(perc)

			ereturn post `b' `V'
			local gname=`group'[1]
					
			di
			display "Estimates for Group `lev' (`group'=`gname')"	
			
			if (e(v_cost)<0){
				di in red "Warning: the estimated variance of the mean costs is negative. The time limit L may be too large."
			}	
			
			if _caller() < 12 {
				ereturn local cmd "hcost"
				_table_ci, level(`level')			
		
			}		
			else {
				_coef_table, level(`level')
			}

		
			use `two', clear
			
		}
			mat `b' = J(1,2,0)
			mat `V' = J(2,2,0)
			
			mat `b'[1,1] = (`mean_c2' - `mean_c1')
			mat `b'[1,2] = (`mean_t2' - `mean_t1')
			mat `V'[1,1] = `v_c2' + `v_c1'
			mat `V'[2,2] = `v_t2' + `v_t1'
			mat `V'[2,1] = `cov2' + `cov1'
			mat `V'[1,2] = `cov2' + `cov1'

			local names cost survival

			mat colnames `b' = `names'
			mat rownames `V' = `names'
			mat colnames `V' = `names'

			ereturn post `b' `V'
	
	
			di
			display as result "Estimates for Difference Between Groups (Group 2 - Group 1)"
			_coef_table,  level(`level')

			if (`method'==1){
				di "Method used: ZT (using cost history)"
			}
			else {
				di "Method used: BT (using total costs)"
			}
					display "Annual discounting rate for costs: `drcost'"
					display "Annual discounting rate for survival time: `drsurv'"	
					di "Time limit L (in days): `l'"
						
			if (`icer'==1){
				mata: CalICER(`mean_c2',`mean_c1',`mean_t2',`mean_t1',	///
					`v_c2', `v_c1', `v_t2',`v_t1',`cov2',`cov1',`level')
					
				local icer: di %8.3f e(icer)
				local icer_low: di %8.3f e(icer_low) 
				local icer_up: di %8.3f e(icer_up)
				
				di
				if (`icer_low'<`icer_up'){
					display as result "Incremental cost-effectiveness ratio (`level'% CI): " `icer' " (" `icer_low' ", " `icer_up' ")" 
				}
				else{
					display as result "Incremental cost-effectiveness ratio (`level'% CI): " `icer' " (-Inf, " `icer_up' ")U(" `icer_low' ", Inf)" 
				}
			}

			ereturn clear
		
			ereturn mat V2 V2 
			ereturn mat b2 b2
			ereturn mat V1 V1
			ereturn mat b1 b1 
			
			ereturn scalar n1= `n1'
			ereturn scalar nobs1 = `nobs1'
			ereturn scalar pcensor1 = `pcensor1'
			ereturn scalar n2= `n2'
			ereturn scalar nobs2 = `nobs2'
			ereturn scalar pcensor2 = `pcensor2'
			
		use `one', clear
	}

end

program _table_ci
        if !c(noisily) {
                exit
        }
		if "`e(cmd)'" == "" {
			error 301
        }

        syntax [, Level(cilevel) COEFLegend SELEGEND cformat(passthru)]
        _get_diopts ignore, `cformat'
        local cfmt `"`s(cformat)'"'

        if "`coeflegend'`selegend'" != "" {
                if "`e(over)'" != "" {
                        local depname `"depname(Over)"'
                }
                else    local depname `"depname(" ")"'
                local coefttl `"coeftitle(`e(depvar)')"'
                _coef_table,    `depname'       ///
                                `coefttl'       ///
                                `coeflegend'    ///
                                `selegend'      ///
                                `cformat'
                exit
        }

        local ematrix "e(error)"
		
        mata: _st_sum_table(`level') 
		
end



 /****************************************************************/

version 11.2

mata:

mata set matastrict on

void CalCostSurv(
	string scalar touses,					///
	string scalar ixs,						///
	string scalar cids,						///
	string scalar survs,					///
	string scalar deltas,					///
	string scalar starts,					///
	string scalar stops,					///
	string scalar costs,					///
	real scalar drsurv,							///
    real scalar drcost,							///
	real scalar k,							///	
	real scalar nL,							///
	real scalar method							///
	)
{
	real scalar i, j, l, n, nobs, flag        
	real scalar CostMeanSW, CostVarSW, CostMeanIMP, CostVarIMP
	real scalar meanadd, varsub, covsub
	real scalar tmean, tvar, CovSW, CovIMP
	real scalar percens, constc
	
	real vector censor, tdelta, tsurv, tcost, s, kc
	real vector id, surv, delta, dsurv
	real vector logic, tmp, tmp1
	
	real matrix info, data, newdata
	
	id = st_data(.,cids, ixs)
	surv = st_data(.,survs, ixs)
	delta = st_data(.,deltas, ixs)

	pragma unset constc
	pragma unset meanadd
	pragma unset varsub
	pragma unset covsub


	data = st_data(., (cids, starts, stops, costs), touses)

	n 	 = rows(id)
	nobs = rows(data)
		
	st_numscalar("e(n_obs)", nobs)
		
	if(drcost>0)
	{
    //Divide the period into k-days 
	
		for(i=1;i<=nobs;i++){
			if(data[i,3]-data[i,2]>=k)
			{
				l=trunc((data[i,3]-data[i,2])/k)
				tmp1=data[i,4]
				tmp2=data[i,3]
				data[i,3]=data[i,2]+k-1
				data[i,4]=data[i,4]/(tmp2-data[i,2]+1)*k
			
				if(l>1)
				for(j=1;j<l;j++)
				{
					newdata=(data[i,1],data[i,2]+k*j,data[i,3]+k*j,data[i,4])
					data=(data \ newdata)
				}
				newdata=(data[i,1],data[i,2]+k*l,tmp2,tmp1/(tmp2-data[i,2]+1)*(tmp2-(data[i,2]+k*l)+1))
				data=(data \ newdata)
			}
		}
		nobs = rows(data)
		data=sort(data,1)
	}

	info = panelsetup(data,1)

	
	if(drcost>0)
	{	
	//Discount cost
		for(i=1;i<=n;i++){
			for (j=info[i,1]; j<=info[i,2]; j++) {
					if(data[j,2]<0)data[j,2]=0			//treat start<0 as 0
					if(drcost==1) data[j,4]=0
					else data[j,4] = data[j,4] * (1-drcost)^(data[j,2]/365.25)
			}
		}	
	}
	
	censor	= J(n,1,1)
	tcost	= J(n,1,0)
	s	= J(n,1,0)
	kc	= J(n,1,0)
	
	
	// TruncSurv()

	tdelta	= delta
	tsurv	= surv
	
	logic = (surv :>= nL)
	
	for (i=1; i<=n; i++) {
		if (logic[i]) {
			tsurv[i] = nL
			tdelta[i] = 1
		}
	}
	
	// CalCensor()
	// calculate percentage of observations censored

	tmp = 1 :- tdelta
	percens = mean(tmp)

	
	// Calculate total cost, and prorate cost if surv is truncated
	for (i=1; i<=n; i++) {
		for (j=info[i,1]; j<=info[i,2]; j++) {
			if (data[j,2]<=tsurv[i]) {
				if (data[j,3]>tsurv[i]) {
					tcost[i] = tcost[i] + data[j,4]*(tsurv[i]-data[j,2]+1) / (data[j,3]-data[j,2]+1)
				}
				else tcost[i] = tcost[i] + data[j,4]
			}
		}
	}
	censor = censor - tdelta
	KmCal(tsurv,censor,s)
	KmCal(tsurv,tdelta,kc)
	
	CostMeanSW = CalOurMean(tdelta, kc, tcost)
	CostVarSW = CalOurVar(tdelta, tsurv, s, kc, tcost, CostMeanSW)
	
	if (method==1)
	{
		CalMeanAdd(tsurv, tdelta, kc, s, id, data, info, tcost, meanadd, varsub, constc)
		CostMeanIMP = CostMeanSW + meanadd
		CostVarIMP = CostVarSW - varsub
	}
	
	// discount survival
	if(drsurv>0) dsurv = (365.25/drsurv) :* ( 1 :- exp(-drsurv :* (tsurv:/365.25)) )
	else dsurv=tsurv
	
	tmean = CalOurMean(tdelta, kc, dsurv)
	tvar = CalOurVar(tdelta, tsurv, s, kc, dsurv, tmean)
	CovSW = CalOurCov(tdelta, tsurv, s, kc, tcost, CostMeanSW, tmean, dsurv)
	
	if (method==1)
	{
		CalCovAdd(tsurv, tdelta, kc, s, id, data, info, covsub, dsurv)
		CovIMP = CovSW - covsub
	}
	
	if(method==1)st_numscalar("e(mean_cost)", CostMeanIMP)
	else st_numscalar("e(mean_cost)", CostMeanSW)
	
	if(method==1)st_numscalar("e(v_cost)", CostVarIMP)
	else st_numscalar("e(v_cost)", CostVarSW)
	
	if(method==1) st_numscalar("e(cov)", CovIMP)
	else st_numscalar("e(cov)", CovSW)
	
	st_numscalar("e(mean_time)", tmean)
	st_numscalar("e(v_time)", tvar)
	st_numscalar("e(n_id)", n)
	st_numscalar("e(perc)", percens)
}




void CalCostSurv0(
	string scalar touses,					///
	string scalar ixs,						///
	string scalar cids,						///
	string scalar survs,					///
	string scalar deltas,					///
	string scalar costs,					///
	real scalar drsurv,							///
	real scalar nL							///
	)
{
	real scalar i, j, l, n, nobs        
	real scalar CostMeanSW, CostVarSW
	real scalar meanadd, varsub, covsub
	real scalar tmean, tvar, CovSW
	real scalar percens, constc
	
	real vector censor, tdelta, tsurv, tcost, s, kc
	real vector id, surv, delta, dsurv
	real vector logic, tmp, tmp1
	
	real matrix info, data
	
	id = st_data(.,cids, ixs)
	surv = st_data(.,survs, ixs)
	delta = st_data(.,deltas, ixs)

	
	pragma unset constc
	pragma unset meanadd
	pragma unset varsub
	pragma unset covsub

	
	
	data = st_data(., (cids, costs), touses)
	info = panelsetup(data,1)

	n 	= rows(id)
	nobs	= rows(data)
	
	censor	= J(n,1,1)
	tcost	= J(n,1,0)
	s	= J(n,1,0)
	kc	= J(n,1,0)
	
	
	// TruncSurv()

	tdelta	= delta
	tsurv	= surv
	
	logic = (surv :>= nL)
	
	for (i=1; i<=n; i++) {
		if (logic[i]) {
			tsurv[i] = nL
			tdelta[i] = 1
		}
	}
	
	// CalCensor()
	// calculate percentage of observations censored

	tmp = 1 :- tdelta
	percens = mean(tmp)
	
	// CalTCost()
	for (i=1; i<=n; i++) {
		for (j=info[i,1]; j<=info[i,2]; j++) {
			if (surv[i]>tsurv[i]) {
				tcost[i] = tcost[i] + data[j,2]*tsurv[i] / surv[i]
			}
			else tcost[i] = tcost[i] + data[j,2]
		}
//		printf("%9.0g,%9.3f,%9.3f,%9.3f\n",i,tcost[i],tsurv[i],censor[i])
	}
	censor = censor - tdelta
	KmCal(tsurv,censor,s)
	KmCal(tsurv,tdelta,kc)
	
	CostMeanSW = CalOurMean(tdelta, kc, tcost)
	CostVarSW = CalOurVar(tdelta, tsurv, s, kc, tcost, CostMeanSW)

	if(drsurv>0) dsurv = (365.25/drsurv) :* ( 1 :- exp(-drsurv :* (tsurv:/365.25)) )
	else dsurv=tsurv
	tmean = CalOurMean(tdelta, kc, dsurv)
	tvar = CalOurVar(tdelta, tsurv, s, kc, dsurv, tmean)
	CovSW = CalOurCov(tdelta, tsurv, s, kc, tcost, CostMeanSW, tmean, dsurv)
	
	st_numscalar("e(mean_cost)", CostMeanSW)	
	st_numscalar("e(v_cost)", CostVarSW)
	st_numscalar("e(mean_time)", tmean)
	st_numscalar("e(v_time)", tvar)
	st_numscalar("e(cov)", CovSW)
	st_numscalar("e(n_id)", n)
	st_numscalar("e(n_obs)", nobs)
	st_numscalar("e(perc)", percens)
}

//K-M estimator 
void KmCal(		
	real vector tsurv,					///
	real vector tcensor,				///
	real vector kc						///
	)
{
	
	real scalar i, n, orgp, est, var 
	real vector orgy, orgs, orgr, orgd, kcorg, index
		
	n   = rows(tsurv)
	est = 1
	orgp = 0
	orgy = orgs = orgr = orgd = kcorg = index = J(n+1,1,0)
	
	datorg(tsurv,tcensor,orgy,orgs,orgr,orgd,orgp,index)
	
	for (i=1;i<=orgp;i++) {
		est = est * (orgr[i]-orgd[i])/orgr[i]
		kcorg[i]=est
		if (i>n) break
	}
	
	for (i=1; i<=n; i++) {
		kc[i] = kcorg[index[i]]	
//		printf("%9.0g, %9.0g, %9.3f, %9.3f\n", i,index[i],tsurv[i],kc[i])
	}

}

//sort data by follow-up time, and record numbers of events, risk set, etc 
void datorg(							///
	real vector time,					///
	real vector censor,					///
	real vector orgy,					/// 
	real vector orgs,					///
	real vector orgr,					/// 
	real vector orgd,					/// 
	real scalar orgp,					/// 
	real vector index
	)
{
	
	real scalar i, j, p, k, n, temp
	real vector l, ctime, ccensor
	
	n = rows(time)
	l = indexor= J(n+1,1,0)
	
	ctime = time
	ccensor = censor
	
	for (i=1; i<=n; i++) indexor[i]=i
	
	for (i=1; i<=n; i++) {

		p = i
		
		for (j=i+1; j<=n; j++) {
			if (ctime[p] > ctime[j]) {
				p = j		//record p for the smallest time since i
			}
		}
		if(p!=i){		//swap subject i & p if i is not smallest time
			temp = ctime[i]
			ctime[i] = ctime[p]
			ctime[p] = temp
			temp = ccensor[i]
			ccensor[i] = ccensor[p]
			ccensor[p] = temp
			temp=indexor[i]	//record index of pre-sorted time
			indexor[i]=indexor[p]
			indexor[p]=temp	
		}
	}
	
	p = 1
	for (i=1; i<n; i++) {
		if (ctime[i] != ctime[i+1]) {
			l[p] = i		//record only diff times
			index[indexor[i]]=p
			p++
		}
		else index[indexor[i]]=p		//record the index of pre-sorted times
	}
	
	l[p] = n	// !!
	index[indexor[n]]=p
	p++

	for (i=p; i<=n+1; i++) {
		l[i] = 0
		orgs[i] = 0
		orgd[i] = 0
		orgr[i] = 0
		orgy[i] = 0
	}
	
	for(i=1; i<p; i++) {
		orgy[i] = ctime[l[i]]
		orgs[i] = 0
		if (i==1) {
			k = 1
		}
		else k = l[i-1] + 1
			
		orgr[i] = n-k+1
		orgd[i] = 0

		for(j=k; j<=l[i]; j++) {
			if (ccensor[j]!=1) orgd[i] = orgd[i] + 1
			orgs[i] = orgs[i] + 1
		}
	}
	
	orgp = p-1
}


real scalar CalOurMean(					///
	real vector tdelta, 				///
	real vector kc,						///
	real vector tcost					///
	)
{
	
	real vector tmp
	
	tmp = (tdelta :== 1)
	tmp = (tcost :/ kc) :* tmp
	
	return(mean(tmp))
}


real scalar CalOurVar(					///
	real vector tdelta,					///
	real vector tsurv,					///
	real vector s,						///
	real vector kc,						///
	real vector tcost,					///
	real scalar mymean					///
	)
{

	real scalar i, j, n
	real scalar e, f, myvar, temp1, temp2
	
	n = rows(tdelta)
	temp1 = temp2 = 0
	
	for (i=1; i<=n; i++) {
		if (tdelta[i]==1) temp1 = temp1 + (tcost[i]-mymean)^2 / kc[i]
	}
  	temp1 = temp1 / n
  	
  	for (j=1; j<=n; j++) {
  		e = f = 0
  		if (tdelta[j]==0){
  			
			for (i=1; i<=n; i++){
				if(tdelta[i]==1 & tsurv[i]>=tsurv[j]) {
					e = e + tcost[i] / kc[i]
					f = f + tcost[i]^2 / kc[i]
				}
			}
			
			e = e / (s[j]*n)
			f = f / (s[j]*n)
			temp2 = temp2 + (f-e*e) / (kc[j]*kc[j])
		}

	}
	
	temp2 = temp2 / n
	
	myvar = temp1 + temp2
	myvar = myvar / n

	return(myvar)

}


void CalMeanAdd(						///
	real vector tsurv, 					///
	real vector tdelta,					///
	real vector kc,						///
	real vector s,						///
	real vector id,						///
	real matrix data,					///
	real matrix info,					///
	real vector tcost,					///
	real scalar meanadd,				///
	real scalar varsub,					///
	real scalar constc					///
	)
{
	
	real scalar n, nobs, i, j, k, y, beta
	real scalar gm, gmmu, gmu, gmustar, gmu2star
	real scalar part1, part2, part3, culcost
	
	
	n = rows(id)
	beta = 1
	part1 = part2 = part3 = 0


	// calculate ebar[j] and y[j] at censoring places
	
	for (j=1; j<=n; j++) {
		if (tdelta[j]==0) {
		gm=gmmu=gmu=gmustar=gmu2star=0
		y=0
			for (i=1; i<=n; i++) {
				if (tsurv[i]>=tsurv[j]) {
					culcost = 0
					for (k=info[i,1]; k<=info[i,2]; k++) {
						if (data[k,1]==id[i]&&data[k,2]<=tsurv[j]){
							if (data[k,3]>tsurv[j]) {
								culcost = culcost + data[k,4]*(tsurv[j]-data[k,2]+1) / (data[k,3]-data[k,2]+1)
							}
							else culcost = culcost + data[k,4]
						}
					}
					gmustar = gmustar + culcost
					gmu2star = gmu2star + culcost*culcost
					y++
					if(tdelta[i]==1) {
						gm = gm + tcost[i]/kc[i]
						gmmu = gmmu + tcost[i]*culcost/kc[i]
						gmu = gmu + culcost/kc[i]
					}
				}
			}
			gmustar = gmustar/y
			gmu2star = gmu2star/y
			gm = gm/(s[j]*n)
			gmmu= gmmu/(s[j]*n)
			gmu = gmu/(s[j]*n)
			part1 = part1 + (tcost[j]-gmustar)/kc[j]
			part2 = part2 + (gmmu-gm*gmu)/(kc[j]*kc[j])
			part3 = part3 + (gmu2star-gmustar*gmustar)/(kc[j]*kc[j])
		}
	}

	part1 = part1 / n

	meanadd = beta*part1
	varsub = (2*beta*part2-beta*beta*part3) / (n*n)
	constc = part2 / part3
	
}



real scalar CalOurCov(					///
	real vector tdelta,					///
	real vector tsurv,					///
	real vector s,						///
	real vector kc, 					///
	real vector tcost,					///
	real scalar mymean,					///
	real scalar tmean,					///
	real vector tstarsurv				///
	)
{
	
	real scalar i, j, n, gtc, gt, gc, mycov, temp1, temp2
	
	n = rows(tdelta)
	temp1 = temp2 = 0
	
	for (i=1; i<=n; i++) {
		if (tdelta[i]==1) temp1 = temp1 + tcost[i]*tstarsurv[i] / kc[i]
	}
	
	temp1 = temp1 / n
	temp1 = temp1 - mymean*tmean
	
	for (j=1; j<=n; j++) {
		gtc = gt = gc = 0
		if (tdelta[j]==0) {
			for (i=1; i<=n; i++) {
				if (tdelta[i]==1 & tsurv[i]>=tsurv[j]) {
					gtc = gtc + tcost[i]*tstarsurv[i] / kc[i]
					gt = gt + tstarsurv[i] / kc[i]
					gc = gc + tcost[i] / kc[i]
				}
			}
			gtc = gtc / (s[j]*n)
			gt = gt / (s[j]*n)
			gc = gc / (s[j]*n)
			temp2 = temp2 + (gtc-gt*gc) / (kc[j]*kc[j])
		}
	}
	temp2 = temp2 / n
	
	mycov = temp1 + temp2
	mycov = mycov / n
	
	return(mycov)
}


void CalCovAdd(							///
	real vector tsurv,					///
	real vector tdelta,					///
	real vector kc,						///
	real vector s,						///
	real vector id,						///
	real matrix data,					///
	real matrix info,					///
	real scalar covsub,					///
	real vector dsurv					///
	)
{

	real scalar i, j, k, n, nobs, beta
	real scalar gt, gmu, gtmu, part2, culcost
	real vector y
	
	n = rows(id)
	part2 = 0
	beta = 1
	y =  J(n,1,0)
	
	// First calculate ebar[j] and y[j] at censoring places
	
	for (j=1; j<=n; j++) {
if (tdelta[j]==0) {
	gt=gtmu=gmu=0
	for (i=1; i<=n; i++) {
		if (tdelta[i]==1&&tsurv[i]>=tsurv[j]) {
		culcost = 0
			for (k=info[i,1]; k<=info[i,2]; k++) {
				if (data[k,1]==id[i]&&data[k,2]<=tsurv[j]){
					if (data[k,3]>tsurv[j]) {
						culcost = culcost + data[k,4]*(tsurv[j]-data[k,2]+1) / (data[k,3]-data[k,2]+1)
					}
					else culcost = culcost + data[k,4]
					}
				}
				gt = gt + dsurv[i]/kc[i]
				gmu = gmu + culcost/kc[i]
				gtmu = gtmu + dsurv[i]*culcost/kc[i]			
			}
		}
		gt = gt/(s[j]*n)
		gmu = gmu/(s[j]*n)
		gtmu = gtmu/(s[j]*n)
		part2 = part2+ (gtmu-gmu*gt)/(kc[j]*kc[j])
		}
	}
	covsub = beta*part2/(n*n)
}

void CalICER(
	real scalar mean_c2,							///
	real scalar mean_c1,							///
	real scalar mean_t2,							///
	real scalar mean_t1,							///
	real scalar v_c2,							///
	real scalar v_c1,							///
	real scalar v_t2,							///
	real scalar v_t1,							///
	real scalar cov2,							///
	real scalar cov1,							///
	real scalar level
	)
{
	real scalar icer        
	real scalar icer_up, icer_low
	real scalar x, y, sxx, syy, sxy
	
	x = mean_c2 - mean_c1
	y = mean_t2 - mean_t1
	icer=x/y
	
	sxx=v_c2+v_c1
	syy=v_t2+v_t1
	sxy=cov2+cov1
	
	z=invnormal(0.5-level/200)^2

	f=(x*y-z*sxy)^2-(x*x-z*sxx)*(y*y-z*syy)
	icer_low=(x*y-z*sxy-sqrt(f))/(y*y-z*syy)
	icer_up=(x*y-z*sxy+sqrt(f))/(y*y-z*syy)
				
	st_numscalar("e(icer)", icer)	
	st_numscalar("e(icer_up)", icer_up)
	st_numscalar("e(icer_low)", icer_low)
}


void _st_sum_table(real scalar level)
{
	// info pulled from local macros
	string	scalar	bmat
	string	scalar	vmat
	string	scalar	emat
	// info pulled from e()
	real	scalar	k_eq
	string	matrix	stripe
	// work
	real	matrix	ms_info
	real	scalar	i
	// coefficient table class	
	class _sum_table scalar	T

	// get info from local macros
	bmat	= st_local("bmatrix")
	vmat	= st_local("vmatrix")
	emat	= st_local("ematrix")

	if (bmat == "") {
		bmat = "e(b)"
	}
	stripe	= st_matrixcolstripe(bmat)
	ms_info	= panelsetup(stripe, 1)
	k_eq	= rows(ms_info)
	if (bmat == "e(b)") {
		bmat = ""
	}
	T.set_b(bmat)
	T.set_V(vmat)
	T.set_errmat(emat)
	T.set_df(_ct_get_scalar("e(df_r)"))
	T.set_level(level)
	T.set_cformat(st_local("cfmt"))

	T.validate()
	T.di_titles()
	for (i=1; i<=k_eq; i++) {
		T.di_eq(i)
	}
	T.finish()
}

mata mlib create lhcost, replace
mata mlib add lhcost *()

end

exit
