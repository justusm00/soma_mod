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
SEED=51854
TIMESTEPS=136


/Users/justusmulthaup/soma_mod/build/testing/../python_src/ConfGen.py -i coord.xml --no-ana-output

if [ ! $? -eq 0 ]; then
    exit -1
fi

/opt/homebrew/bin/mpiexec -n 1 /Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA -c coord.h5 -t 0  -n 2
mv end.h5 start.h5


if [ ! $? -eq 0 ]; then
    exit -5
fi

echo A
/opt/homebrew/bin/mpiexec -n 2 /Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA -c start.h5 -r $SEED -t $TIMESTEPS --accepted-load-inbalance=0 -l 1  -n 2
if [ ! $? -eq 0 ]; then
    exit -2
fi
mv end.h5 load.h5

echo AA
/opt/homebrew/bin/mpiexec -n 2 /Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA -c start.h5 -r $SEED -t $TIMESTEPS  -n 2
if [ ! $? -eq 0 ]; then
    exit -3
fi
mv end.h5 noload.h5

echo B
/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10 compare_mixed_bead_data.py load.h5 noload.h5
echo C
if [ ! $? -eq 0 ]; then
    exit -4
fi
