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


/opt/homebrew/bin/mpiexec  -n 1 /Users/justusmulthaup/soma_mod/build/testing/../c_src/CONVERT coord.dat
value=$?

if [ $value -eq 0 ]
then
    echo -n ""
else
    echo Test failed on convert
    exit -1
fi

./compareDatH5.py coord.dat coord.h5

value=$?

if [ $value -eq 0 ]
then
    echo Test passed
else
    echo Test failed
    exit -1
fi