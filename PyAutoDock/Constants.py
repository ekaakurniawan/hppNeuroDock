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

# References:
#  - AutoDock 4.2.3 Source Code
#    http://autodock.scripps.edu
#  - IEEE-754 Analysis Calculator
#    http://babbage.cs.qc.edu/IEEE-754/
#  - numeric_limits::epsilon
#    http://msdn.microsoft.com/en-us/library/6x7575x3(v=vs.110).aspx


# Mathematical constant pi (22/7)
# The value of pi (based on IEEE-754 Analysis Calculator):
#  - Decimal                  : 3.14159265358979323846
#  - Single precision (32-bit): 3.1415925 or 0x40490FDA
#  - Double precision (64-bit): 3.1415926535897931 or 0x400921FB54442D18
PI = 3.1415926535897931

# Two times mathematical constant pi (2 * (22/7))
# The value of 2 x pi (based on IEEE-754 Analysis Calculator):
#  - Decimal                  : 6.28318530717958647692
#  - Single precision (32-bit): 6.2831850 or 0x40C90FDA
#  - Double precision (64-bit): 6.2831853071795862 or 0x401921FB54442D18
TWOPI = 6.2831853071795862

# To avoid divide-by-zero
# Based on C++ numeric_limits::epsilon for float (4bytes) is 1.19209E-7 and for
# double (8bytes) is 2.22045E-16. For FPGA implementation use following.
#  - Single precision (32-bit): 0.999999997E-6 or 0x358637BD
#  - Double precision (64-bit): 0.9999999999999998E-15 or 0x3CD203AF9EE75615
APPROX_ZERO = 0.9999999999999998E-15

DEG2RAD = PI / 180.0
