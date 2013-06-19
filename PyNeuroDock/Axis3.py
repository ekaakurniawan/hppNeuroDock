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


import math
import Constants as const

# Axis: x, y, z
class Axis3(object):
    @property
    def xyz(self):
        return [self.x, self.y, self.z]
    @xyz.setter
    def xyz(self, value):
        self.x = value[0]
        self.y = value[1]
        self.z = value[2]

    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return "Axis3: %6.2f, %6.2f, %6.2f" % (self.x, self.y, self.z)

    def __add__(self, other):
        return Axis3(self.x + other.x, self.y + other.y, self.z + other.z)
    
    def __sub__(self, other):
        return Axis3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __isub__(self, other):
        self.x -= other.x
        self.y -= other.y
        self.z -= other.z
        return self

    def sq_hypotenuse(self):
        return self.x * self.x + self.y * self.y + self.z * self.z

    def normalize(self):
        mag_xyz = math.sqrt(self.sq_hypotenuse())
        if mag_xyz > const.APPROX_ZERO:
            self.x /= mag_xyz
            self.y /= mag_xyz
            self.z /= mag_xyz
