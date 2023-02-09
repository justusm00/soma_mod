#!/bin/bash

#   Copyright (C) 2016-2021 Ludwig Schneider
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


rm -f mobility*h5

/Users/justusmulthaup/soma_mod/build/testing/../python_src/ConfGen.py -i mobility_low_concentration_eq.xml --no-ana-output

if [ ! $? -eq 0 ]; then
    exit -1
fi

/Users/justusmulthaup/soma_mod/build/testing/../python_src/ConfGen.py -i mobility_high_concentration_eq.xml --no-ana-output

if [ ! $? -eq 0 ]; then
    exit -2
fi

/opt/homebrew/bin/mpiexec -n 1 /Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA -c mobility_low_concentration_eq.h5 -f mobility_low_concentration.h5  -t 20000  -n 2

if [ ! $? -eq 0 ]; then
    exit -3
fi

/opt/homebrew/bin/mpiexec -n 1 /Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA -c mobility_high_concentration_eq.h5 -f mobility_high_concentration.h5 -t 20000  -n 2

if [ ! $? -eq 0 ]; then
    exit -3
fi

/Users/justusmulthaup/soma_mod/build/testing/../python_src/ConfGen.py -i mobility_high_concentration.xml --update

if [ ! $? -eq 0 ]; then
    exit -2
fi
/Users/justusmulthaup/soma_mod/build/testing/../python_src/ConfGen.py -i mobility_low_concentration.xml --update

if [ ! $? -eq 0 ]; then
    exit -2
fi


/opt/homebrew/bin/mpiexec -n 1 /Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA -c mobility_low_concentration.h5 -a mobility_low_concentration_ana.h5 -t 20000  -n 2

if [ ! $? -eq 0 ]; then
    exit -3
fi

/opt/homebrew/bin/mpiexec -n 1 /Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA -c mobility_high_concentration.h5 -a mobility_high_concentration_ana.h5 -t 20000  -n 2

if [ ! $? -eq 0 ]; then
    exit -4
fi

/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10 test_mobility.py

if [ ! $? -eq 0 ]; then
    exit -5
fi

rm -f mobility*h5
