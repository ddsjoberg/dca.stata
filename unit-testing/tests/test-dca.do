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
