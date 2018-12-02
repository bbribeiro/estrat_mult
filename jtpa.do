clear
set more off

cd "C:\Users\Bernardo\Dropbox\WB\job\code\github\estrat_mult"

use jtpa, clear
set seed 1234
	
gen control=assignmt
	recode control 0=1 1=0
	
/* ------------------------ */
/* --- ORIGINAL PROGRAM --- */
/* ------------------------ */
mata: mata clear
cap program drop estrat
run "ado/estrat.ado"
set seed 1234

estrat earnings assignmt prevearn, reps(10) boot(10)

/* ------------------------ */
/* --- MODIFIED PROGRAM --- */
/* ------------------------ */
mata: mata clear
cap program drop estrat_mult
run "ado/estrat_mult.ado"
set seed 1234

/* Replicating original model */
estrat_mult earnings assignmt, potential(earnings) pred(prevearn) control(control) reps(10) boot(10)

/* Generate random multiple treatment arms */
set seed 1234
gen random = uniform()
gen ssize = _N
gen snumber = _n
gen treatment = 1
	replace treatment = 2 if snumber > (ssize/5) & snumber <= (ssize/5*2)
	replace treatment = 3 if snumber > (ssize/5*2) & snumber <= (ssize/5*3) 
	replace treatment = 4 if snumber > (ssize/5*3) & snumber <= (ssize/5*4) 
	replace treatment = 5 if snumber > (ssize/5*4)
	replace treatment = 0 if control == 1

tab treatment, gen(treatment)
	
/* Generate random potential outcome */
gen random_pot=runiform()


/* Run estrat_mul with multiple treatment */
estrat_mult earnings treatment2 treatment3 treatment4 treatment5 treatment6, potential(random_pot) pred(prevearn) control(control) reps(10) boot(10)
