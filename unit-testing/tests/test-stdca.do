* unit tests for the dca function

* load data
use "example-data\dca.dta", clear
stset ttcancer, f(cancer)

tempfile dca_cancerpredmarker
test_command stdca cancerpredmarker, nograph timepoint(1.5) xby(0.49) xstop(0.50) saving(`dca_cancerpredmarker')
use `dca_cancerpredmarker', clear

* check NB and interventions avoided are correct
test_assert "NB is correct: all[1]"  abs(all[1] -  .21223511) < 0.0001
test_assert "NB is correct: all[2]"  abs(all[2] -  -.55977452) < 0.0001
test_assert "NB is correct: none[1]"  abs(none[1] -  0) < 0.0001
test_assert "NB is correct: none[2]"  abs(none[2] -  0) < 0.0001
test_assert "NB is correct: cancerpredmarker[1]"  abs(cancerpredmarker[1] -  .21282192) < 0.0001
test_assert "NB is correct: cancerpredmarker[2]"  abs(cancerpredmarker[2] -  .05783202) < 0.0001
test_assert "NB is correct: cancerpredmarker_i[1]"  abs(cancerpredmarker_i[1] -  0.058093967) < 0.0001
test_assert "NB is correct: cancerpredmarker_i[2]"  abs(cancerpredmarker_i[2] -  0.61760654) < 0.0001
