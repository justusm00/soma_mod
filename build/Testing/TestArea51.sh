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
SEED=11
TIMESTEPS=42

/opt/homebrew/bin/mpiexec -n 1 /Users/justusmulthaup/soma_mod/build/testing/../python_src/ConfGen.py -i area51.xml --no-ana-output

if [ ! $? -eq 0 ]; then
    exit -1
fi

/opt/homebrew/bin/mpiexec -n 1 /Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA -c area51.h5 -r $SEED -t $TIMESTEPS  -n 2
if [ ! $? -eq 0 ]; then
    exit -2
fi
