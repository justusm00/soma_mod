# CMake generated Testfile for 
# Source directory: /Users/justusmulthaup/soma_mod/testing
# Build directory: /Users/justusmulthaup/soma_mod/build/testing
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(quick "./quick.sh")
set_tests_properties(quick PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;73;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(walltime "./walltime.sh")
set_tests_properties(walltime PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;79;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(convert "./convert.sh")
set_tests_properties(convert PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;86;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(xml "./xml.sh")
set_tests_properties(xml PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;92;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(self-documentation "./self_documentation.sh")
set_tests_properties(self-documentation PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;98;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(async "/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10" "async.py" " -n 2")
set_tests_properties(async PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;105;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(AnaGen "./AnaGen.sh")
set_tests_properties(AnaGen PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;112;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(StructureFactor "/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10" "structure_factor.py" " -n 2")
set_tests_properties(StructureFactor PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;121;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(ExternalField "/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10" "external_field.py")
set_tests_properties(ExternalField PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;127;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(Restart "./RestartExact.sh")
set_tests_properties(Restart PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;140;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(Load-Balance "./LoadBalance.sh")
set_tests_properties(Load-Balance PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;153;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(tags-and-type-switching "./Tag.sh")
set_tests_properties(tags-and-type-switching PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;161;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(TestArea51 "./TestArea51.sh")
set_tests_properties(TestArea51 PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;168;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(PolytypeConversion "/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10" "conversion_test.py" "--additional-flags= -n 2")
set_tests_properties(PolytypeConversion PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;200;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(PolytypePartialConversion "/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10" "partial_conversion_test.py" "--additional-flags= -n 2")
set_tests_properties(PolytypePartialConversion PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;207;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(mobility-test "./mobility.sh")
set_tests_properties(mobility-test PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;216;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(weight-test "/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10" "weight.py" "--prefix=/opt/homebrew/bin/mpiexec -n 2" "--additional-flags= -n 2")
set_tests_properties(weight-test PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;218;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(default-stat "/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10" "statistics.py" "--prefix=/opt/homebrew/bin/mpiexec -n 2" "--additional-flags= -n 2")
set_tests_properties(default-stat PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;224;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(MPI-stat "/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10" "statistics.py" "--prefix=/opt/homebrew/bin/mpiexec -n 2" "--additional-flags= -n 2")
set_tests_properties(MPI-stat PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;226;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(MT-stat- "/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10" "statistics.py" "--prefix=/opt/homebrew/bin/mpiexec -n 2" "--additional-flags=-pMT -n 2")
set_tests_properties(MT-stat- PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;229;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(SET-stat- "/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10" "statistics.py" "--prefix=/opt/homebrew/bin/mpiexec -n 2" "--additional-flags=--iteration-alg=SET -n 2")
set_tests_properties(SET-stat- PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;237;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
add_test(SET-FIXED-stat- "/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10" "statistics.py" "--prefix=/opt/homebrew/bin/mpiexec -n 2" "--additional-flags=--iteration-alg=SET --set-generation-algorithm=FIXED -n 2")
set_tests_properties(SET-FIXED-stat- PROPERTIES  WORKING_DIRECTORY "/Users/justusmulthaup/soma_mod/build/testing" _BACKTRACE_TRIPLES "/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;239;add_test;/Users/justusmulthaup/soma_mod/testing/CMakeLists.txt;0;")
