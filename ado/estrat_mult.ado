*! estrat: version 1.0 Jeremy Ferwerda
*! estrat_mult: adaptation to fit multiple treatments, by Bernardo Ribeiro

cap program drop estrat_mult
program estrat_mult, eclass
	version 11
	syntax varlist(min=2 numeric ts) [if] [in] , potential(varlist min=1 max=1) PRed(varlist) CONtrol(varlist min=1 max=1) [cov(varlist) ///
			reps(numlist max=1 > 0) ///
			boot(numlist max=1 > 0) ///
			groups(numlist max=1 > 1 <= 10) ///
			loo_only ///
			rss_only ///
			savegroup ///
			]
			
	//local loo_only "loo_only"

	local depvar: word 1 of `varlist'
	local treatment: list varlist - depvar
	local narms: word count `treatment'
	local t=1
	tokenize `treatment'
	
	
	// Check that variables vary
	foreach i in `varlist' `potential' `pred' `control' {
		quietly summ `i'
		if `r(Var)' == 0 {
			display as error "All variables must vary"
			exit
		}	
	}

	// Check that control is binary
	qui: levelsof `control', local(tmp)
	foreach x of local tmp{
		if (`x' != 1 & `x'!=0){
			display as error "The control indicator must be a binary variable."
			exit
		}
	}
	
	// Check that treatments are binary
	foreach t in `treatment' {
		qui: levelsof `t', local(tmp)
		foreach x of local tmp {
			if (`x' != 1 & `x'!=0){
				display as error "The treatment indicators must be binary variables."
				exit
			}
		}
	}
	
	// Clean up if user previously terminated execution.
	capture drop looholder* 
	capture drop rssholder*
	capture drop betaa*
	capture drop estrat_loo_group
	
	// Defaults
	if ("`groups'" == "") local groups = 3
	if ("`reps'" == "") local reps = 100
	if ("`boot'" == "") local boot = 100
	
	// Counts
	qui: count if `control' == 0
	local ntreat = `r(N)'
	qui: count if `control' == 1
	local nuntreat = `r(N)'
		
	// Holders for reps
	forvalues p = 1/`groups' {
		forvalues t=1/`narms' {
			qui: gen looholder`p'`t'=.
			qui: gen rssholder`p'`t'=.
			qui: gen betaa`p'`t'=.
		}
		local p = `p' + 1
	}
	
	// Bootstrap Wrapper
	di as text "Boot Reps " _continue

	local maxrep = `boot' + 1 // Using the first rep for point estimates.
	if (`maxrep' > _N) set obs `maxrep'


	forvalues i = 1/`maxrep'{
		di as text "." _continue
		
		// Sample	
		preserve
		if (`i' != 1) bsample _N
		
		if ("`rss_only'" == ""){
		// LOO Estimation
			qui: reg `potential' `pred' if `control'==1, noheader notable 
			qui: predict fhat
			qui: predict lev, leverage
			qui: replace fhat = fhat - ( (lev/(1-lev))*(`potential'-fhat)  ) if `control'==1
			//qui: xtile fhat3 = fhat, nq(`groups')
		
		 	// Manual to match matlab
		 	if ("`savegroup'" != "") qui: gen estrat_oldorder = _n
			sort fhat
			qui: gen fhat_r = .
			qui: summ fhat
			local lastcap = `r(min)'
			local p = 1

			qui: count if fhat != .
			local newtotal = `r(N)'
		
			while `p' < `groups'{
				qui: summ fhat if _n == round(`newtotal'/(`groups'/`p')), meanonly
				local z = `p' + 1
				
				if (`p' == 1){
					qui: replace fhat_r = `p' if fhat <= `r(mean)'
				} 
				else {
					 qui: replace fhat_r = `p' if fhat <= `r(mean)' & fhat > `lastcap'
				}
				 qui: replace fhat_r = `z' if fhat > `r(mean)'
				
				local p = `p' + 1
				local lastcap = `r(mean)'
				
			}
			qui: replace fhat_r = . if fhat==.

			forvalues g = 1/`groups' {
				qui: capture reg `depvar' `treatment' `cov' if fhat_r==`g', noheader notable
				matrix b=e(b)
				if (`g' == 1){
					matrix looh = b[1,1..`narms']'
				}
				else {	
					if _rc!=0 {
					 matrix looh = looh, J(`narms',1,0)
					} 
					else { 
					matrix looh = looh, b[1,1..`narms']'
					}
				}
			}
		
			if (`i' != 1){
				matrix loo_se = looh
			} 
			else {
				matrix loo = looh
				if ("`savegroup'" != ""){
					qui: sort estrat_oldorder
					mata: m_preservegroup()
				}
			}
			
			}
		
		// RSS Estimation
			if ("`loo_only'" == ""){
			forvalues  m = 1/`reps' {
		
				qui: xtile vm = runiform() if `control'==1, nq(2) 
				qui: replace vm=2 if vm==.
				qui: reg `potential' `pred' if `control'==1 & vm==1, noheader notable
				qui: predict fhatrss
				//xtile fhat_r2 = fhatrss if vm!=1, nq(`groups')  
				
				qui: gen fhat_r2 = .
				qui: summ fhatrss if vm !=1
				local lastcap = `r(min)'
				local p = 1
							
				sort vm fhatrss
				qui: by vm: generate tempid = _n if vm != 1  & fhatrss != .
        
				qui: count if vm != 1 & fhatrss != .
				local newtotal = `r(N)'
        
				while `p' < `groups'{
					
					qui: summ fhatrss if tempid == round(`newtotal'/(`groups'/`p')), meanonly
					local z = `p' + 1
					
					if (`p' == 1){
						qui: replace fhat_r2 = `p' if fhatrss <= `r(mean)'
					} 
					else {
						qui: replace fhat_r2 = `p' if fhatrss <= `r(mean)' & fhatrss > `lastcap'
					}
					qui: replace fhat_r2 = `z' if fhatrss > `r(mean)'
        
					local p = `p' + 1
					local lastcap = `r(mean)'
					qui: replace fhat_r2 = . if vm == 1
				}
				
				qui: drop tempid			
				forvalues g = 1/`groups' {
					capture reg `depvar' `treatment' `cov' if fhat_r2==`g', noheader notable
					matrix b=e(b)
					forvalues t=1/`narms' {
						if _rc!=0 {
						qui: replace betaa`g'`t'= 0 in `m'
						} 
						else {
						qui: replace betaa`g'`t'= b[1,`t'] in `m'	
						}
					}
				}
				
				qui: drop vm fhatrss fhat_r2
			} 
			
			matrix rssh=J(`narms',`groups',0)
			forvalues g = 1/`groups' {
				forvalues t=1/`narms' {
					qui: capture summ betaa`g'`t', meanonly			
					if (`g' == 1){
						matrix rssh[`t',1] = `r(mean)'
					}
					else {	
						matrix rssh[`t',`g'] = `r(mean)'
					}
				}
			}
			
			
			// Matrix pass out of restore
			if (`i' != 1){
				matrix rss_se = rssh
			} 
			else {
				matrix rss = rssh
			}
			}
		restore
		
		// Save bootsrap iterations for later use
		if (`i' != 1){
			forvalues g = 1/`groups' {
				forvalues t=1/`narms' {
					if ("`rss_only'" == "") qui: replace looholder`g'`t'=loo_se[`t',`g'] in `i'
					if ("`loo_only'" == "") qui: replace rssholder`g'`t'=rss_se[`t',`g'] in `i'
				}
			}
		}
	} 
	qui: drop betaa*
	
	if ("`savegroup'" != "") mata: m_savegroup()
	
	preserve

	// Get Bootstrap SEs
	forvalues g = 1/`groups' {
		forvalues t=1/`narms'{
			if ("`rss_only'" == "") qui: egen sloo`g'`t' = sd(looholder`g'`t')
			if ("`loo_only'" == "") qui: egen srss`g'`t' = sd(rssholder`g'`t')
		}
	}


	// Clean up
	if ("`rss_only'" == "") qui: matrix sdloo=J(`narms',`groups', 0)
	if ("`loo_only'" == "") qui: matrix sdrss=J(`narms',`groups', 0)
	forvalues g=1/`groups' {
		forvalues t=1/`narms' {
			if ("`rss_only'" == "") {
				qui summ sloo`g'`t', meanonly
				matrix sdloo[`t',`g']=`r(mean)'
			}
			if ("`loo_only'" == "") {
				qui summ srss`g'`t', meanonly
				matrix sdrss[`t',`g']=`r(mean)'
			}
		}
	}
	
	restore 
	
	capture drop looholder*
	capture drop rssholder*  
	
	// Output table

	mata: m_estrat("`loo_only'","`rss_only'")
	
	// Ereturn
	ereturn post
	ereturn local depvar "`depvar'"
	ereturn local treatment "`treatment'"	
	ereturn local potential "`potential'"
	ereturn local predictors "`pred'"
	ereturn local covariates "`cov'"
	
	ereturn scalar groups=	`groups'
	ereturn scalar rssrep=		`reps'
	ereturn scalar bootrep=		`boot'
	ereturn scalar treated=		`ntreat'
	ereturn scalar untreated=	`nuntreat'

	if ("`rss_only'" == ""){	
	ereturn matrix LOO_C = loo
	ereturn matrix LOO_SE = sdloo
	}
	if ("`loo_only'" == ""){
	ereturn matrix RSS_C = rss
	ereturn matrix RSS_SE = sdrss	
	}

end

mata:
mata set matastrict on
mata set matafavor speed

void  m_preservegroup(){
	external real matrix groups
	real matrix X
	st_view(X, ., "fhat_r")
	groups = X
}

void m_savegroup(){
   real matrix T
   external real matrix groups
   (void) st_addvar("double","estrat_loo_group")

   st_view(T, ., "estrat_loo_group")
   T[.,.] = groups
}

void m_estrat(string scalar loo_only, string scalar rss_only){

real matrix loo, rss, sdloo, sdrss
real scalar i, quantile, nt, nu, reps, overflow, oq, boot, narms, t
string scalar depvar, covar, pred, potential

depvar = st_local("depvar")
covar = st_local("cov")
pred = st_local("pred")
potential = st_local("potential")

reps = strtoreal(st_local("reps"))
boot = strtoreal(st_local("boot"))
nt = strtoreal(st_local("ntreat"))
nu = strtoreal(st_local("nuntreat"))
quantile = strtoreal(st_local("groups"))
narms = strtoreal(st_local("narms"))

if (rss_only == ""){
loo = st_matrix("loo")
sdloo = st_matrix("sdloo")
}
if (loo_only == ""){
sdrss = st_matrix("sdrss")
rss = st_matrix("rss")
}

// Output results
printf("\n\n")
printf("{txt}Treated   = %11.0f\n", nt)
printf("{txt}Untreated = %11.0f\n", nu)
printf("{txt}Treatment arms = %6.0f\n", narms)
printf("{txt}RSS Reps  = %11.0f\n", reps)
printf("{txt}Boot Reps = %11.0f\n", boot)

for (t=1; t<=narms; t++) {

	printf("\n\n")
	printf("{txt}Treatment arm:   "+st_local(strofreal(t))+"\n")
	printf("{hline 9}{c TT}{hline %2.0f}\n",44)
	printf("         {c |}")

	if (rss_only == ""){
	printf("{space 7}LOO{space 9}SE{space 1}")
	}
	if (loo_only == ""){
	printf("{space 7}RSS{space 9}SE")
	}

	printf("\n{hline 9}{c +}{hline %2.0f}",44)

	for (i=1; i<= quantile; i++){
		printf("\n")
		printf("{txt}{space 1}Group %1.0f {c |}  ",i)
		
		if (rss_only == ""){
			printf("{res}%8.0g {space 2}", loo[t,i])
			printf("{res}%8.0g {space 2}", sdloo[t,i])	
		}
		
		if (loo_only == ""){
			printf("{res}%8.0g {space 2}", rss[t,i])
			printf("{res}%8.0g {space 2}", sdrss[t,i])
		}
	}		
	printf("\n{hline 9}{c BT}{hline %2.0f}\n",44)

}


printf("\n{txt}Predictors:          %s\n", pred)


 if (potential == depvar){
printf("{txt}Potential outcome:   %s\n", "Identical to dependent variable")
 }
 else {
printf("{txt}Potential outcome:   %s\n", potential)
 }
	
 if (covar == pred){
printf("{txt}Covariates: %s\n", "Identical to predictors")
 } 
 else  if (covar != ""){
printf("{txt}Covariates: %s\n", covar)
 } 
} 
end
