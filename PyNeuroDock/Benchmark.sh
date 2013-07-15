#!/bin/bash

# Copyright (C) 2013 by Eka A. Kurniawan
# eka.a.kurniawan(ta)gmail(tod)com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the
# Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

# To run: ./Benchmark.sh 2>&1 | tee Benchmark/Benchmark.txt

echo "========================================================================="
echo "Run Time Benchmark Python Sequentially"
echo "========================================================================="
echo ""

FILES=./Benchmark/Python/*
for f in $FILES
do
    echo "Running... $f"
    time /opt/local/bin/python2 NeuroDock.py -p $f
    echo "-------------------------------------------------------------------------"
done

echo "========================================================================="
echo "Run Time Benchmark Python-OpenCL GPU"
echo "========================================================================="
echo ""

FILES=./Benchmark/Python-OpenCL_GPU/*
for f in $FILES
do
    echo "Running... $f"
    time /opt/local/bin/python2 NeuroDock.py -p $f
    echo "-------------------------------------------------------------------------"
done

echo "========================================================================="
echo "Run Time Benchmark Python-OpenCL CPU"
echo "========================================================================="
echo ""

FILES=./Benchmark/Python-OpenCL_CPU/*
for f in $FILES
do
    echo "Running... $f"
    time /opt/local/bin/python2 NeuroDock.py -p $f
    echo "-------------------------------------------------------------------------"
done

