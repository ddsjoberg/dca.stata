/*
PROGRAM: dca.ado
PROGRAMMER: Daniel
DATE: 9/16/2014
DESCRIPTION:
dca calculates the points on a decision curve and optionally
plots the decision curve, where <yvar> is a binary response and
<xvars> are predictors of the binary response. 
The default is that all <xvars> specified are probabilities. If this
is not the case, then the user must specify the <prob> option. 

*/

capture program drop dca
program dca, rclass
	version 12.0

	syntax varlist(min=2 numeric) [if] [in], [xstart(numlist >0 <1 max=1 missingokay) xstop(numlist >0 <1 max=1 missingokay) ///
				xby(numlist >0 <1 max=1 missingokay) saving(string asis) smooth smoother(string) noGRAPH harm(string) ///
				PREValence(numlist max=1) INTERvention ///
				interventionper(real 100) interventionmin(real 0) ymin(real -0.05) ymax(real 1.0) PROBability(string) *]
	preserve
	
	*keeping observation in if/in
	capture keep `if' `in'
	
	*initializing default values for xstart xstop and xby and smoother
	if trim("`xstart'")=="" local xstart=0.01
	if trim("`xstop'")=="" local xstop=0.99
	if trim("`xby'")=="" local xby=0.01
	if trim("`smooth'")!="" & trim("`smoother'")=="" local smoother="3rssh"
	local probability=upper("`probability'")
	
	*extracting outcome varname
	local outcome=trim("`=word("`varlist'",1)'")
	*removing outcome variable from varlist
	local varlist=trim("`=subinstr("`varlist'","`outcome'","",1)'")
	

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

	*checking that the  probabilities are specified as YES or NO
	if trim("`probability'")!="" {
		foreach i of numlist 1/`=wordcount("`probability'")' {
			capture assert inlist("`=word("`probability'",`i')'","NO","YES")
			if _rc>0 {
				disp as error "Probabilities must be specified as Yes or No only."
				exit 198
			}
		}
	} 
	
	*asigning each predictor, harm, and probability an ID and default value if not specified.
	local varn=0
	foreach v of varlist `varlist' {
		local ++varn /*getting number of predictors*/
		local var`varn'="`v'"
		if trim("`harm'")!="" local harm`varn'="`=word("`harm'",`varn')'"
		else local harm`varn'=0
		if trim("`probability'")!="" local prob`varn'="`=word("`probability'",`varn')'"
		else local prob`varn'="YES"
	}

	
	*model variable names being checked cannot be equal to "all" or "none"
	foreach v of varlist `outcome' `varlist' {
		if inlist(trim("`v'"),"all","none") {
			disp as error `"Variable names cannot be equal to "all" or "none""'
			exit 198
		}
	}
	
	*only keeping observations that are not missing model values
	foreach v of varlist `outcome' `varlist' {
		qui keep if !mi(`v')
	}

	*asserting outcome coded as 0-1
	capture assert inrange(`outcome',0-0.000001,1+0.00001)
	if _rc>0 {
		disp as error "`outcome' cannot be less than 0 or greater than 1."
		exit 198
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
			logit `outcome' `var`i''
			predict `temppred`i''
			replace `var`i''=`temppred`i''
			noi disp "`var`i'' converted to a probability with logistic regression. Due to linearity assumption, miscalibration may occur."
		}
	}
	

	

	
	******************************************************************************************************************
	*** calculate net benefit for each threshold
	******************************************************************************************************************
	tempname dcamemhold
    tempfile dcaresults	
	postfile `dcamemhold' threshold str100(model) nb using `dcaresults'

	*looping over every threshold probability and calculating NB for all models
	if "`prevalence'" == "" {
		qui summ `outcome' 
			local prevalence=`r(mean)'
	}
	else  {
	    capture assert inrange(`prevalence',0,1)
		if "`prob`i''"=="YES" & _rc>0 {
			noi disp as error "prevalence argument must be between 0 and 1."
			exit 198
		}
	}
	 
	local N=_N
	return local N=_N
	local tcount=0
	local threshold=`xstart'-`xby'
	qui while `threshold'<`xstop' {
		local threshold=`threshold'+`xby'
		local ++tcount
		
			

		* creating var to indicate if observation is at risk
		qui foreach model in all none `varlist' {
		
			if "`model'"=="all" post `dcamemhold' (`threshold') ("`model'") (`prevalence'-(1-`prevalence')*`threshold'/(1-`threshold'))
			else if "`model'"=="none" post `dcamemhold' (`threshold') ("`model'") (0)
			else {
					
				*calculating TP and FP. 
				*TP=proportion of patients with disease given they are above the 
				*   threshold multiplied by number of men above threshold
				tempvar _testpos_
				gen `_testpos_' = `model'>=`threshold'
				sum `_testpos_' if `outcome'
				*if no patients above threshold then tp is 0
					if `r(N)'==0 {
						local tp_rate=0
					}
					else {
					*counting true positives
						local tp_rate=`r(mean)' * `prevalence'
					}	
				sum `_testpos_' if !`outcome'
					*if no patients above threshold then fp is 0
					if `r(N)'==0 {
						local fp_rate=0
					}
					else {
					*counting falseS positives
						local fp_rate=`r(mean)' * (1 - `prevalence')
					}		
			    drop `_testpos_'
				
				*posting the net benefit to results table
				post `dcamemhold' (`threshold') ("`model'") (`tp_rate' - (`fp_rate')*(`threshold'/(1-`threshold')))
				
				*grabbing variable label for otuput dataset
				local `model'label: variable label `model'
				if trim("``model'label'")=="" local `model'label `model'

			}
		} /*closing varlist loop*/
	} /*closing threshold loop*/
	postclose `dcamemhold'
		
		
	*loading NB calculations
	use `dcaresults', clear
	
	* applying harms if specified
	qui foreach i of numlist 1/`varn'{
		*modifying NB to account for harm
		qui replace nb=nb-`harm`i'' if trim(model)==trim("`var`i''") /*if harms not specified, then harm`i' was defaulted to zero earlier in the program*/
	}
	

	* making dataset oneline per threshold value
	qui reshape wide nb, i(threshold) j(model) string

	
	* applying variable lables, and creating intervention vars if requested.
	sort threshold 
	
	rename nball all
	rename nbnone none
	qui foreach i of numlist 1/`varn' {
		rename nb`var`i'' `var`i''
		
		if trim("`harm`i''")=="0" label variable `var`i'' "Net Benefit: ``var`i''label'"
		else label variable `var`i'' "Net Benefit: ``var`i''label' (`harm`i'' harm applied)"
		
		*transforming variables for interventions avoided
		qui g `var`i''_i= (`var`i'' - all)*`interventionper'/(threshold/(1-threshold))
		
		if trim("`harm`i''")=="0" label variable `var`i''_i "Intervention: ``var`i''label'"
		else label variable `var`i''_i "Intervention: ``var`i''label' (`harm`i'' harm applied)"
	}
	label variable threshold "Threshold Probability"
	label variable all "Net Benefit: Treat All"
	label variable none "Net Benefit: Treat None"
		
	*smoothing data if requested, and labelling new variables
	else if trim("`smooth'")=="smooth" {
		foreach i of numlist 1/`varn' {
			quietly smooth `smoother' `var`i'', gen(sm_`var`i'')
			
			if trim("`harm`i''")=="0" label var sm_`var`i'' "Smoothed Net Benefit: ``var`i''label'"
			else label var sm_`var`i'' "Smoothed Net Benefit: ``var`i''label' (`harm`i'' harm applied)"
			
			g sm_`var`i''_i= (sm_`var`i'' - all)*`interventionper'/(threshold/(1-threshold))
			if trim("`harm`i''")=="0" label var sm_`var`i''_i "Smoothed Intervention: ``var`i''label'"
			else label var sm_`var`i''_i "Smoothed Intervention: ``var`i''label' (`harm`i'' harm applied)"		
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
		
	* creating title for NB or net benefit for figure
	if trim("`intervention'")=="" local plottitle=`""Net Benefit""'
	else local plottitle=`""Net reduction in interventions" "per `interventionper' patients""'
		
		
	***********************************************
	***********  PLOTTING FIGURE  *****************
	***********************************************
	line `plotvarlist' threshold, cmissing(n) ytitle(`plottitle', margin(medium)) xtitle(, margin(medium)) `options'
	
end       
		
		
		
		
		
		
	
	


	
		
	

