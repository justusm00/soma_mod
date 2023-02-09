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

DELTA=0

/Users/justusmulthaup/soma_mod/build/testing/../python_src/ConfGen.py -i restartexact.xml --no-ana-output

if [ ! $? -eq 0 ]; then
    exit -1
fi

/opt/homebrew/bin/mpiexec  -n 1 /Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA -c restartexact.h5 -r 11 -t 1  -n 2
mv end.h5 step0.h5

if [ ! $? -eq 0 ]; then
    exit -5
fi

/opt/homebrew/bin/mpiexec  -n 1 /Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA -c restartexact.h5 -r 11 -t 1  -n 2
mv end.h5 norestart.h5

if [ ! $? -eq 0 ]; then
    exit -2
fi

/Users/justusmulthaup/soma_mod/build/testing/../python_src/ConfGen.py -i restartexact.xml -o step0.h5 --no-ana-output --update

if [ ! $? -eq 0 ]; then
    exit -6
fi
mv end.h5 restart.h5

/opt/homebrew/bin/mpiexec  -n 1 /Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA -c step0.h5 -r 11 -t 0  -n 2
if [ ! $? -eq 0 ]; then
    exit -3
fi
mv end.h5 restart.h5

/opt/homebrew/bin/h5diff -v -d $DELTA norestart.h5 restart.h5
if [ ! $? -eq 0 ]; then
    exit -4
fi
