*! version 0.2
capture program drop sfsd_p
program define sfsd_p
version 16

	syntax newvarname, [xb Residuals u su te ste mu]

	local spatialwvars `e(spatialwvars)'
	if `"`spatialwvars'"'=="" {
		if "`su'"!="" | "`ste'"!="" {
			di as error "su and ste cannot be specified in the non-spatial models"
			exit
		}
	}
	else{
		cap confirm var `spatialwvars'
		if _rc {
			di as error "For the postestimation, genwvars should be specified in the last estimate"
			exit 
		}
	}
	
	local cost "`e(function)'"
	local cost  =cond("`cost'" == "cost", -1,1)
	local nopts: word count `xb'  `residuals'  `u' `te' `mu' `su' `ste'
    if `nopts' >1 {
        display "{err}only one statistic may be specified"
        exit 498
    }
	if `nopts'==0 local xb xb


	if "`mu'"!=""{
        _predict `typlist' `varlist'  , xb eq(Mu)
        label var `varlist' "Prediction of Mu"
   }	   

	if "`su'"!=""{
		qui gen `typlist' `varlist' = __u_sp__ 
		label var `varlist'  "spatial-corrected inefficiency(u)"
		exit	
	}  

	if "`ste'"!=""{
		qui gen `typlist' `varlist' = exp(-__u_sp__)
		label var `varlist'  "spatial-corrected efficiency"
		exit	
	} 

	local yvar = e(depvar)

    if "`xb'" != "" {
        _predict `typlist' `varlist'  , xb
		local flagy = e(user)
		local flagy = (strpos("`flagy'","spsf"))
		if `flagy'>0{
			local rho = e(rho)
			qui replace `varlist' =`varlist' + `rho'*Wy_`yvar'
		}
    }
	else{
		tempvar yhat 
        qui _predict double `yhat'  , xb
		local flagy = e(user)
		local flagy = (strpos("`flagy'","spsf"))
		if `flagy'>0{
			local rho = e(rho)
			qui replace `yhat'=`yhat' + `rho'*Wy_`yvar'
		}			
	}
	
	if "`residuals'"!=""{
		qui gen `typlist' `varlist' = `yvar'-`yhat'
		label var `varlist' "prediction of residuals"
	}
	
	if "`omega'"!="" | "`u'"!="" | "`te'"!=""{

		local endovars  "`e(endvars)'"
		if "`endovars'"=="" & "`omega'"!=""{
			 display "{err}omega can only be specified in models with endogenous variables"
             exit 498
		}
		tempvar corterm
		qui gen double `corterm' = 0
		foreach v in `endovars'{
			tempvar r`v'
			qui _predict double `r`v''  , xb eq(`v')
			qui replace `r`v'' = (`v' - `r`v'' )
			qui replace `corterm' = `corterm' + `r`v''*_b["/eta_`v'"]
		}
		
	}
	
	if "`u'"!=""| "`te'"!="" | "`omega'"!=""{
		
		tempvar mu lnsigmav lnsigmau sigma2 mustar sigma2s
		qui _predict double `mu', xb eq(#2)
		qui _predict double `lnsigmav', xb eq(#3)
		qui _predict double `lnsigmau', xb eq(#4)
		if "`omega'"!=""{
			qui gen double `varlist' = `yvar'-`yhat' - `corterm'*exp(`lnsigmav')
			exit 
		}
		tempvar omega
		local endovars  "`e(endvars)'"
		if "`endovars'"==""{
			qui gen double `omega' = `yvar'-`yhat' 
		}
		else{
			qui gen double `omega' = `yvar'-`yhat' - `corterm'*exp(`lnsigmav')
		}

		qui replace `omega' = `cost'*`omega'	
		qui gen double `sigma2' = exp(2*`lnsigmav') + exp(2*`lnsigmau')
		qui gen double `mustar' = (exp(2*`lnsigmav') *`mu' - exp(2*`lnsigmau')*`omega')/`sigma2'
		qui gen double `sigma2s' = exp(2*`lnsigmav')* exp(2*`lnsigmau')/`sigma2'
		
		if "`u'"!=""{
			qui gen `typlist' `varlist' = `mustar' + sqrt(`sigma2s')*normalden(`mustar'/sqrt(`sigma2s'))/normal(`mustar'/sqrt(`sigma2s'))
			label var `varlist' "Prediction of inefficiency: u"
		}
		
		if "`te'"!=""{
			qui gen `typlist' `varlist' = exp(-`mustar' - sqrt(`sigma2s')*normalden(`mustar'/sqrt(`sigma2s'))/normal(`mustar'/sqrt(`sigma2s')))
			label var `varlist' "Prediction of efficiency: JLMS estimator"			
		}
	}
	
end


