clear
set more off

cd "C:\Users\Bernardo\Dropbox\WB\job\code\github\estrat"

use jtpa, clear
set seed 1234
	
gen control=assignmt
	recode control 0=1 1=0
	
gen treat1 = round(runiform())
gen treat2 = round(runiform())
gen treat3 = round(runiform())

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

estrat_mult earnings assignmt, potential(earnings) pred(prevearn) control(control) reps(10) boot(10)
