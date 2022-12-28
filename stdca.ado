/*
PROGRAM: stdca.ado
PROGRAMMER: Daniel Sjoberg
DATE: 9/16/2014
DSECRIPTION:
stdca is an extension to dca for a survival-time endpoint. 
The program calculates the points on a decision curve and optionally
plots the decision curve, where the inputs are the probabilities of failure 
at time t associated with a single variable or multivariable
model. 
*/

capture program drop stdca
program stdca
	version 12.0

	syntax varlist(min=1 numeric) [if] [in], Timepoint(real) [xstart(numlist >0 <1 max=1 missingokay) xstop(numlist >0 <1 max=1 missingokay) ///
				xby(numlist >0 <1 max=1 missingokay) saving(string asis) smooth smoother(string) noGRAPH harm(string) INTERvention ///
				interventionper(real 1) interventionmin(real 0) ymin(real -0.05) ymax(real 1.0) PROBability(string) ///
				compet1(numlist missingokay) compet2(numlist missingokay) compet3(numlist missingokay) ///
				compet4(numlist missingokay) compet5(numlist missingokay) compet6(numlist missingokay) *]
	preserve
	
	*keeping selected observations
	capture keep `if' `in'
	
	*initializing default values for xstart xstop and xby and smoother
	if trim("`xstart'")=="" local xstart=0.01
	if trim("`xstop'")=="" local xstop=0.99
	if trim("`xby'")=="" local xby=0.01
	if trim("`smooth'")!="" & trim("`smoother'")=="" local smoother="3rssh"
	
	* confirming data is stset
	if "`_dta[_dta]'"!="st" {
		disp as error "data not st"
		exit 119
	}
	
	*checking that the same number of harms are specified as variables to check
	if wordcount("`varlist'")!=wordcount("`harm'") & trim("`harm'")!="" {
		disp as error "Number of harms specified must be the same as the number of predictors being checked."
		exit 198
	}

	*checking that the same number of probabilities are specified as variables to check
	if wordcount("`varlist'")!=wordcount("`probability'") & trim("`probability'")!="" {
		disp as error "Number of probabilities specified must be the same as the number of predictors being checked."
		exit 198
	}

	*asigning each predictor, harm, and probability an ID and default value if not specified.
	local varn=0
	foreach v of varlist `varlist' {
		local ++varn /*getting number of predictors*/
		local var`varn'="`v'"
		if trim("`harm'")!="" local harm`varn'="`=word("`harm'",`varn')'"
		else local harm`varn'=0
		if trim("`probability'")!="" local prob`varn'=upper("`=word("`probability'",`varn')'")
		else local prob`varn'="YES"
	}
	
	
	*model variable names being checked cannot be equal to "all" or "none"
	foreach v of varlist `varlist' {
		if inlist(trim("`v'"),"all","none") {
			disp as error `"Variable names cannot be equal to "all" or "none""'
			exit 198
		}
	}
	


	*only keeping observations that are stset and not missing model values
	qui keep if _st==1
	foreach v of varlist `varlist' {
		qui keep if !mi(`v')
	}
	
	*asserting all inputs are between 0 and 1 OR specified as non-probabilites.  If not a probability, then converting it to prob with logistic regression
	qui foreach i of numlist 1/`varn'{
		capture assert inrange(`var`i'',0,1)
		if "`prob`i''"=="YES" & _rc>0 {
			noi disp as error "`var`i'' must be between 0 and 1 OR sepcified as a non-probability in the probability option"
			exit 198
		}
		else if "`prob`i''"=="NO" {
			tempvar temppred`i' 
			tempvar s0`i'
			stcox `var`i'', basesurv(`s0`i'')
			predict `temppred`i'', xb
			qui sum `s0`i'' if _t<=`timepoint'
			replace `temppred`i''=1-`r(min)'^exp(`temppred`i'')

			replace `var`i''=`temppred`i''
			noi disp "`var`i'' converted to a probability with Cox regression. Due to linearity and proportional hazards assumption, miscalibration may occur."
		}
	}



	
	******************************************************************************************************************
	*** calculate net benefit for each threshold
	*** assume the data are already -stset-
	******************************************************************************************************************
	* getting the probability of the event for all subjects
	* this is used for the net benefit associated with treating all patients
		* if cmprsk option is not specified, then use Kaplan-Meier estimates
		if trim("`compet1'")=="" qui stdca_kmciest, timepoint(`timepoint') local(pd)
		*otherwise getting competing risk cum inc estimate
		else {
			*creating competing risk options for stcompet function if competing risks is being utilized.
			foreach i of numlist 1/6 {
				if trim("`compet`i''")!="" local cmprsk=trim("`cmprsk' compet`i'(`compet`i'')")
			}
			
			qui stdca_crciest, timepoint(`timepoint') local(pd) `cmprsk'
		}
	
	tempname stdcamemhold
    tempfile stdcaresults	
	postfile `stdcamemhold' threshold str100(model) nb using `stdcaresults'

	*looping over every threshold probability and calculating NB for all models
	local N=_N
	local tcount=0
	local threshold=`xstart'-`xby'
	while `threshold'+`xby'<=`xstop' {
		local threshold=`threshold'+`xby'
		local ++tcount
		
		* creating var to indicate if observation is at risk
		foreach model in all none `varlist' {
			
			*calculating NB for treat all and treat none
			if "`model'"=="all" post `stdcamemhold' (`threshold') ("`model'") (`pd' - (1 - `pd') * `threshold' / (1 - `threshold'))
			else if "`model'"=="none" post `stdcamemhold' (`threshold') ("`model'") (0) 
			*calculating NB for the rest of the xvars in varlist
			else  {
			
				*creating indicator if observation is above threshold risk
				capture drop `atrisk'
				tempvar atrisk
				g `atrisk'=(`model'>`threshold' & !mi(`model'))
				*determining avergage risk level among obs at risk
				qui summ `atrisk'
					local px=r(mean)
				
				*if no patients are at risk, NB cannot be calculated. 
				if `px'==0 {
					post `stdcamemhold' (`threshold') ("`model'") (.r) 
					continue  /*skipping the rest of the loop*/
				}
				
				*calculating risk among obs in the at risk set and saving out risk in local var pdgivenx.
				if trim("`cmprsk'")=="" qui stdca_kmciest if `atrisk'==1, timepoint(`timepoint') local(pdgivenx)
				else qui stdca_crciest if `atrisk'==1, timepoint(`timepoint') local(pdgivenx) `cmprsk'
								
				*if pdgivenx is missing it is because there is not enough followup in the atrisk group to get the risk estimate at that time.
				if `pdgivenx'==. post `stdcamemhold' (`threshold') ("`model'") (.f) 
				*posting result
				else post `stdcamemhold' (`threshold') ("`model'") (`pdgivenx'*`px' -  (1-`pdgivenx')*`px'*(`threshold'/(1-`threshold'))) 
				
				*grabbing variable label for otuput dataset
				local `model'label: variable label `model'
				if trim("``model'label'")=="" local `model'label `model'
			}
		}
	}
	postclose `stdcamemhold'
	*loading NB calculations
	use `stdcaresults', clear
	
	* applying harms if specified
	local vcount=0
	foreach v in `varlist' {
		local ++vcount
		*creating a local for the harms if specified
		if trim("`harm'")!="" {
			*creating a local for the harms if specified
			local `v'harm: word `vcount' of `harm'
			*modifying NB to account for harm
			qui replace nb=nb-``v'harm' if trim(model)==trim("`v'")
		}
	}
	
	****************************************************
	* creating warnings to print when NB not calculable.
		*no patients at risk
		sort model
		qui by model: egen atriskerrormin=min(threshold) if nb==.r
		sort model atriskerrormin
		qui by model: replace atriskerrormin=atriskerrormin[_n-1] if mi(atriskerrormin) & !mi(atriskerrormin[_n-1])
		*not enough followup among patients at risk
		sort model
		qui by model: egen followuperrormin=min(threshold) if nb==.f
		sort model followuperrormin
		qui by model: replace followuperrormin=followuperrormin[_n-1] if mi(followuperrormin) & !mi(followuperrormin[_n-1])
		
		qui g error=model+": No observations with risk greater than "+trim(string(atriskerrormin*100,"%9.1f"))+"%, and therefore net benefit not calculable in this range." if !mi(atriskerrormin) & atriskerrormin<followuperrormin
		qui replace error=model+": No observations with risk greater than "+trim(string(followuperrormin*100,"%9.1f"))+"% that have followup through the timepoint selected, and therefore net benefit not calculable in this range." if mi(error) & !mi(followuperrormin)
		
		local errcount=0
		qui levelsof error
		foreach err in `r(levels)' {
			local ++errcount
			noi disp "`err'"
		}
		drop followuperrormin atriskerrormin error
		
	* making dataset oneline per threshold value
	qui reshape wide nb, i(threshold) j(model) string
	
	
	* applying variable labels, and creating intervention vars if requested.
	sort threshold 
	foreach v in all none `varlist' {
		rename nb`v' `v'
		
		
		if "`v'"=="all" label variable `v' "Net Benefit: Treat All"
		else if "`v'"=="none" label variable `v' "Net Benefit: Treat None"
		else if trim("`harm'")=="" label variable `v' "Net Benefit: ``v'label'"
		else if trim("`harm'")!="" label variable `v' "Net Benefit: ``v'label' (``v'harm' harm applied)"
		
		*if intervention requested, then transforming variables to interventions avoided
		qui g `v'_i= (`v' - all)*`interventionper'/(threshold/(1-threshold))
			if "`v'"=="all" label variable `v'_i "Intervention: Treat All"
			else if "`v'"=="none" label variable `v'_i "Intervention: Treat None"
			else if trim("`harm'")=="" label variable `v'_i "Intervention: ``v'label'"
			else if trim("`harm'")!="" label variable `v'_i "Intervention: ``v'label' (``v'harm' harm applied)"
	}
	label variable threshold "Threshold Probability"
	
	*smoothing data if requested
	else if trim("`smooth'")=="smooth" {
		foreach v of varlist `varlist' {
			quietly smooth `smoother' `v', gen(sm_`v')
			if trim("`harm'")=="" label var sm_`v' "Smoothed Net Benefit: ``v'label'"
			else if trim("`harm'")!="" label var sm_`v' "Smoothed Net Benefit: ``v'label' (``v'harm' harm applied)"
			
			qui g sm_`v'_i= (sm_`v' - all)*`interventionper'/(threshold/(1-threshold))
			if trim("`harm'")=="" label var sm_`v'_i "Smoothed Intervention: ``v'label'"
			else if trim("`harm'")!="" label var sm_`v'_i "Smoothed Intervention: ``v'label' (``v'harm' harm applied)"		
		}
		
	}

	**  saving out permanent dataset if requested
	order threshold all none
	if trim(`"`saving'"')!="" save `saving'

	*if graph suppression option requested, then skipping the rest of the program
	if trim("`graph'")!="" exit
	
	**  creating variable list for plotting.  command will look like:   line `plotvarlist' threshold
		*if no smoothing plotting varlist is just the varlist
		if trim("`smooth'")=="" local plotvarlist="`varlist'"
		*otherwise, we need to add the sm_* prefix to the variable names
		else {
			foreach v of varlist `varlist' {
				local plotvarlist="`plotvarlist' sm_`v'"
			}
		}
		
		*if plotting NB rather than  interventions, add all none to plot
		if trim("`intervention'")=="" local plotvarlist="all none `plotvarlist'"
		*otherwise add *_i suffix to varlist
		else {
			local plotvarlist2="all_i none_i"
			foreach v of varlist `plotvarlist' {
				local plotvarlist2="`plotvarlist2' `v'_i"
			}
			local plotvarlist="`plotvarlist2'" 
		}
		
	***  deleting values outside ymin and ymax if plotting NB and lower than interventionmin for interventions 
	***  this is a cosmetic option for the graphing only
	foreach v of varlist `plotvarlist' {
		if trim("`intervention'")=="" qui replace `v'=. if !inrange(`v',`ymin',`ymax')
		else qui replace `v'=. if `v'<`interventionmin'
	}
	
	
	* creating title for NB or net beneift for figure
	if trim("`intervention'")=="" local plottitle=`""Net Benefit""'
	else if "`interventionper'"=="1" local plottitle=`""Net reduction in interventions""'
	else local plottitle=`""Net reduction in interventions" "per `interventionper' patients""'
	
	***********************************************
	***********  PLOTTING FIGURE  *****************
	***********************************************
	line `plotvarlist' threshold, cmissing(n) ytitle(`plottitle', margin(medium)) xtitle(, margin(medium)) `options'
	
end       
 

* we often need the Kaplan Meier estimate, this simple program returns it as a local variable
capture program drop stdca_kmciest
program stdca_kmciest, rclass
	syntax [if] , timepoint(real) [local(string)]
	
	preserve
	capture keep `if'
	tempfile kmall
	qui sts list, at(0 `timepoint') saving(`kmall') fail
	use `kmall', clear
	qui keep if time==`timepoint'
	assert _N==1
	if trim("`local'")!="" c_local `local'=failure
end

* we often need the cum inc with competing risk estimates, this simple program returns it as a local variable
* stcompet reference Stata Journal, 4(2):103-112
capture program drop stdca_crciest
program stdca_crciest, rclass
	syntax [if] , timepoint(real) [local(string) *] 
	
	preserve
	capture keep `if'	
	*checking that there are time points above _t
	qui count if _t>=`timepoint' & _d==1
	if `r(N)'==0 & trim("`local'")!="" c_local `local'=.
	
	else {
		tempvar ci
		/*computing cumulative incidence rates*/
		stcompet `ci'=ci, `options'
		qui sum `ci' if _t<=`timepoint' & _d==1
		/*if at least one observation has event before or on time t then calculate cum inc*/
		if `r(N)'>0 & trim("`local'")!="" c_local `local'=`r(max)' 
		/*if no events, then cum inc is 0 with missing confidence bounds */
		else if  trim("`local'")!="" c_local `local'=0
	}
end






