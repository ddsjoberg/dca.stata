* set working directory
cd "C:\Users\sjobergd\GitHub\dca.stata"

* import unit testing functions test_command, test_assert, run_tests
quietly include "https://raw.githubusercontent.com/btskinner/stata_unit_test/master/test_commands.do"

* run unit tests
run_tests, testfiledirectory("unit-testing\tests") stopiferror
