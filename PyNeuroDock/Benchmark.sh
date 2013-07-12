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

# To run: ./Benchmark.sh 2>&1 | tee ./Benchmark/run_time_benchmark.txt

echo "Run Time Benchmark"
echo ""

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_32_10.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_64_10.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_128_10.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_256_10.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_512_10.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_32_50.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_64_50.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_128_50.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_256_50.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_512_50.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_32_100.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_64_100.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_128_100.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_256_100.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_512_100.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf




echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_32_10.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_64_10.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_128_10.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_256_10.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_512_10.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_32_50.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_64_50.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_128_50.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_256_50.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_512_50.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_32_100.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_64_100.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_128_100.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_256_100.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_512_100.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf




echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_256_500.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_512_500.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_1024_500.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_256_1000.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_512_1000.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_1024_1000.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_2560_500.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_opencl_15360_500.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_256_500.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_512_500.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_1024_500.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_256_1000.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_512_1000.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

echo "-------------------------------------------------------------------------"
dpf="Benchmark/ind_runtime_seq_1024_1000.dpf"
echo "Running... $dpf"
time /opt/local/bin/python2 NeuroDock.py -p $dpf

