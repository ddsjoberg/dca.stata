* unit tests for the dca function

* load data
use "example-data\dca.dta", clear

tempfile dca_famhist
test_command dca cancer famhistory, nograph  xby(0.49) xstop(0.50) saving(`dca_famhist')
use `dca_famhist', clear

* check NB and interventions avoided are correct
test_assert "NB is correct: all[1]"  abs(all[1] -  0.13131313) < 0.0001
test_assert "NB is correct: all[2]"  abs(all[2] -  -.72000003) < 0.0001
test_assert "NB is correct: none[1]"  abs(none[1] -  0) < 0.0001
test_assert "NB is correct: none[2]"  abs(none[2] -  0) < 0.0001
test_assert "NB is correct: famhistory[1]"  abs(famhistory[1] -  .03077441) < 0.0001
test_assert "NB is correct: famhistory[2]"  abs(famhistory[2] -  -.08933333) < 0.0001
test_assert "NB is correct: famhistory_i[1]"  abs(famhistory_i[1] -  -995.33337) < 0.0001
test_assert "NB is correct: famhistory_i[2]"  abs(famhistory_i[2] -  63.066669) < 0.0001

* MA Checks
use "example-data\dca.dta", clear
logistic cancer marker
predict pred, pr

*net benefit = sensitivity × prevalence – (1 – specificity) × (1 – prevalence) × w 
*where w is the odds at the threshold probability
*thresholds: 0.1 and 0.2
g pred1=pred>0.1
g pred2=pred>0.2
roctab cancer pred1, detail
*nb = (.74285714*.14) - ((1-.56434109)*(.86)*(.1/.9)) = .06237037
roctab cancer pred2, detail
*nb = (.39047619*.14) - ((1-.88062016)*(.86)*(.2/.8)) = .029

tempfile dca_famhist2
test_command dca cancer pred1, prob( no) nograph prev(.14) xby(0.01) xstop(0.2) saving(`dca_famhist2')
use `dca_famhist2', clear

test_assert "NB is correct: marker[1]" abs(pred1[10]-((.74285714*.14) - ((1-.56434109)*(.86)*(.1/.9))) < 0.001
test_assert "NB is correct: marker[2]" abs(pred1[20]-((.39047619*.14) - ((1-.88062016)*(.86)*(.2/.8)))) < 0.001
