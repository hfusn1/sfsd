* translog and xtsfsp packge should be installed in advance
capture log close
cap cd ./dta-files
adopath ++ ../ado_help_files

/////// Data generating process 1 //////
clear all
* The programme to generate the inefficiency term
capture mata: mata drop rtnorm3() 
mata
real function rtnorm3(
					 real vector mean, 
					 real vector sd, 
					 real vector lower, 
					 real vector upper) 
{   
	nrow=1
	ncol =1
	zlower = J(nrow,ncol,((lower :- mean) :/ sd))
	zupper = J(nrow,ncol,((upper :- mean) :/ sd))
	ind = select(1::length(zupper),zupper:==.)
	nzupper = 1
	u = normal(zlower) :+ (nzupper :- normal(zlower)) :* runiform(rows(mean), cols(mean))
	res = J(nrow,ncol,mean) :+ J(nrow,ncol,sd) :* invnormal(u)
	return(res)
}

end

* Setup
global N = 300
global T = 20
set seed 123456
mata mata matuse w_ex1,replace
mata: w1=wm
global beta 2
global rho  0.3
global theta 0.3
global alpha2 0.5
global gamma 0.2
global phi 0.3
global sigu2 0.4
global sigv2 0.2
global NT=$N * $T

clear 
set obs $NT 
qui{
gen e = rnormal()
gen v = rnormal()
gen x = rnormal()
gen g1 = rnormal()
gen g2 = rnormal()
gen z =1+0.2*g1+ 2*g2+ e
putmata x = x ,replace
putmata z = z ,replace
putmata v= v,replace	
}
 
* Generate the simulation data in Mata
mata{
id= J($T,1,(1::$N))
t = (1::$T)# J($N , 1,1)
us=J($NT,1,.)
y=us
wx=us
wz=wx
for(j=1; j<=$T ; j++){
	wm = w1[(1::$N),(1::$N)]
	ind = mm_which(t:==j)
	vs2 = sqrt( $sigv2 ) *v[ind]
	wx[ind] = wm*x[ind]
	wz[ind] = wm*z[ind]
	um = $alpha2 :+$gamma*z[ind]+$phi *wz[ind]
	us2 = rtnorm3(um,sqrt($sigu2 ),0,. )
	yyyy =  x[ind]*$beta +$theta *wx[ind]+ vs2 -us2
	y[ind] = matinv(I( $N )-$rho*wm)* ( yyyy )
	
}
}

* Get the data from Mata
clear
getmata id =id t=t y=y x=x z=z  
qui xtset id t
save sfsd_ex1,replace

/////// Data generating process 2 //////
clear all
* The programme to generate the inefficiency term
capture mata: mata drop rtnorm3() 
mata
real function rtnorm3(
					 real vector mean, 
					 real vector sd, 
					 real vector lower, 
					 real vector upper) 
{   
	nrow=1
	ncol =1
	zlower = J(nrow,ncol,((lower :- mean) :/ sd))
	zupper = J(nrow,ncol,((upper :- mean) :/ sd))
	ind = select(1::length(zupper),zupper:==.)
	nzupper = 1
	u = normal(zlower) :+ (nzupper :- normal(zlower)) :* runiform(rows(mean), cols(mean))
	res = J(nrow,ncol,mean) :+ J(nrow,ncol,sd) :* invnormal(u)
	return(res)
}

end

*Set up
global N = 300
global T = 10
set seed 123456
mata mata matuse w_ex2, replace
global beta 2
global rho  0.3
global theta 0.3
global alpha2 0.5
global gamma 0.2
global phi 0.3
global sigu2 0.4
global sigv2 0.2
global NT = $N * $T

clear 
set obs $NT 
qui {
        gen e = rnormal()
        gen v = rnormal()
        gen x = rnormal()
        gen g1 = rnormal()
        gen g2 = rnormal()
        gen z = 1 + 0.2 * g1 + 2 * g2 + e
        putmata x = x, replace
        putmata z = z, replace
        putmata v = v, replace	
    }

* Generate the simulation data in Mata
mata {
id = J($T, 1, (1 :: $N))
t = (1 :: $T) # J($N, 1, 1)
us = J($NT, 1, .)
y = us
wx = us
wz = wx

for (j = 1; j <= $T; j++) {
	if (j == 1) wm = w1
	else if (j == 2) wm = w2
	else if (j == 3) wm = w3
	else if (j == 4) wm = w4
	else if (j == 5) wm = w5
	else if (j == 6) wm = w6
	else if (j == 7) wm = w7
	else if (j == 8) wm = w8
	else if (j == 9) wm = w9
	else if (j == 10) wm = w10
	
	ind = mm_which(t :== j)
	vs2 = sqrt($sigv2) * v[ind]
	wx[ind] = wm * x[ind]
	wz[ind] = wm * z[ind]
	um = $alpha2 :+ $gamma * z[ind] + $phi * wz[ind]
	us2 = rtnorm3(um, sqrt($sigu2), 0, .)
	yyyy = x[ind] * $beta + $theta * wx[ind] + vs2 + us2
	y[ind] = matinv(I($N) - $rho * wm) * (yyyy)
}
}

* Get the data from Mata
clear
getmata id = id t = t y = y x = x z = z  
qui xtset id t
save sfsd_ex2,replace

///////4.1 SDF-STE model with time-constant spatial weight matrix  //////
clear all
sjlog using ../sfsd1, replace
* Setup
mata mata matuse w_ex1
use sfsd_ex1.dta

* Specify the constraint for multiple equations
constraint 1 [Frontier]x=10*[Mu]z

* Stochastic Durbin production model 
sfsd y x, id(id) time(t) noconstant wymat(wm, mata) wxmat(wm, mata) ///
wumat(wm, mata) mu(z) wxvars(x) wmuvars(z) constraints(1) genwvars nolog
sjlog close, replace

sjlog using ../sfsd2, replace
* Postestimation
predict uhat, u
predict te, te
predict ste, ste
* Summary statistics
summarize uhat te ste
sjlog close, replace

///////4.2 SDF-STE model with time-varying spatial weight matrixs //////
clear all
sjlog using ../sfsd3, replace
* Setup
mata mata matuse w_ex2
use sfsd_ex2.dta
local w w1 w2 w3 w4 w5 w6 w7 w8 w9 w10
* Stochastic Durbin cost model
sfsd y x, id(id) time(t) cost noconstant wmat(`w', mata) mu(z)     ///
wxvars(x) wmuvars(z) mex(x) meu(z) nolog
sjlog close, replace


///////4.3 Empirical example //////
clear all
sjlog using ../sfsd4
* Create a spatial-weighting matrix object (spmat object)
spmat use matrix using matrix.spmat,replace
* Copy matrix to mata
spmat getmatrix matrix wm

use example.dta
xtset city year
* Local the explained variable
local y lnGDP
* Local our key explanatory variable
local x lntour
* Local the input lnK_lnL lnK_2 lnL_2
local input lnK lnL lnK_lnL lnK_2 lnL_2
* Local other determinants
local control lnfin lnpop lnmkt

*Fit the SDF-STE model
sfsd `y' `x' `input', id(city) time(year) wmat(wm, mata) wxvars(`x') ///
mu(`x' `control') wmuvars(`x') mex(`x') meu(`x') genwvars delve mlplot

* Postestimation of economic efficiency
predict efficiency, ste
summarize efficiency
* Histograms and a kernel density estimate for economic efficiency 
hist efficiency, kdensity lcolor(black) fcolor(white) lwidth(thin)   ///
lpattern(solid) 
sjlog close, replace

graph export figure1.eps,replace
