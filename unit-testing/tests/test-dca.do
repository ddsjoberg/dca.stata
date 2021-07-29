* unit tests for the dca function

* load data
use "example-data\dca.dta", clear

test_command dca cancer famhistory, nograph
