* set working directory
cd "C:\Users\sjobergd\GitHub\dca.stata"

* read in current dca and stdca functions
quietly include "dca.ado"
quietly include "stdca.ado"

* import unit testing functions test_command, test_assert, run_tests
quietly include "unit-testing\test_commands.do"

* run unit tests
run_tests, testfiledirectory("unit-testing\tests") stopiferror
