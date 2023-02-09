#!/bin/bash

#   Copyright (C) 2020-2021 Ludwig Schneider
#
# This file is part of SOMA.
#
# SOMA is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# SOMA is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with SOMA.  If not, see <http://www.gnu.org/licenses/>.

/Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA --version


/Users/justusmulthaup/soma_mod/build/testing/../python_src/ConfGen.py -i coord.xml --no-ana-output

value=$?

if [ $value -eq 0 ]
then
    echo Test passed
else
    echo Test failed
    exit -1
fi


/opt/homebrew/bin/mpiexec  -n 1 /Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA -c coord.h5 -t 1 -f la.h5 --purpose="Test step 1"

value=$?

if [ $value -eq 0 ]
then
    echo Test passed
else
    echo Test failed
    exit -1
fi

if [ ! -f la.h5 ]; then
    echo "Test failed, output file 'la.h5' missing"
    exit -1
fi

/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10 -c "import h5py; f=h5py.File('la.h5'); exit(str(f['documentation'][0]).find('Test step 1')<0)"
value=$?
if [ $value -eq 0 ]
then
    echo Test passed
else
    echo Test failed
    exit -1
fi

/opt/homebrew/bin/mpiexec  -n 1 /Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA -c la.h5 -t 1 -f la.h5 --purpose="Test step 2"

value=$?

if [ $value -eq 0 ]
then
    echo Test passed
else
    echo Test failed
    exit -1
fi

if [ ! -f la.h5 ]; then
    echo "Test failed, output file 'la.h5' missing"
    exit -1
fi


/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10 -c "import h5py; f=h5py.File('la.h5'); exit(str(f['documentation'][0]).find('Test step 1')<0)"
value=$?
if [ $value -eq 0 ]
then
    echo Test passed
else
    echo Test failed
    exit -1
fi

/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10 -c "import h5py; f=h5py.File('la.h5'); exit(str(f['documentation'][1]).find('Test step 2')<0)"
value=$?
if [ $value -eq 0 ]
then
    echo Test passed
else
    echo Test failed
    exit -1
fi

rm la.h5
