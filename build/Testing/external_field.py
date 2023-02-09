#!/opt/homebrew/Frameworks/Python.framework/Versions/3.10/bin/python3.10

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
import sys
sys.path.append("/Users/justusmulthaup/soma_mod/build/testing/../python_src")
import h5py
import argparse
import subprocess as sp
import numpy as np
SOMA = "/Users/justusmulthaup/soma_mod/build/testing/../c_src/SOMA"
CONVERT = "/Users/justusmulthaup/soma_mod/build/testing/../c_src/CONVERT"
CONFGEN = "/Users/justusmulthaup/soma_mod/build/testing/../python_src/ConfGen.py"
HANDLEANAH5 = "/Users/justusmulthaup/soma_mod/build/testing/../python_src/handleAnaH5.py"
epsilon = 2.

def convert_files():
    ret = sp.call([CONFGEN,"-i","external_field.xml"])
    if ret != 0: raise RuntimeError("Test failed")

def test():
    f = h5py.File("external_field.h5", 'r')

    external_field = f['/external_field']
    period = f['/external_field'].attrs['period']
    cos= np.array(f['/external_field'].attrs['cos'])
    value_field=external_field[0][0][0][0]
    if period==100 and cos[1]==1 and value_field==1:
       return 0

    else:
        print("The external field is not correctly defined")
        return -1


def main_test():
    convert_files()
    test()


if __name__ == "__main__":
    main_test()
