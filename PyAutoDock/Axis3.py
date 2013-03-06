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
class Axis3:
    def __init__(self, x = 0.0, y = 0.0, z = 0.0):
        self.x = x
        self.y = y
        self.z = z

    def __repr__(self):
        return "Axis3: %.4f, %.4f, %.4f" % (self.x, self.y, self.z)

    def __add__(self, move):
        return Axis3(self.x+move.x, self.y+move.y, self.z+move.z)
    
    def __sub__(self, move):
        return Axis3(self.x-move.x, self.y-move.y, self.z-move.z)

    def transform(self, rot, move):
        tx=self.x*rot[0]+self.y*rot[3]+self.z*rot[6] + move.x
        ty=self.x*rot[1]+self.y*rot[4]+self.z*rot[7] + move.y
        tz=self.x*rot[2]+self.y*rot[5]+self.z*rot[8] + move.z
        self.x=tx
        self.y=ty
        self.z=tz
        return self

    def normalize(self):
        mag_xyz = math.sqrt(self.x * self.x + \
                            self.y * self.y + \
                            self.z * self.z)
        if mag_xyz > const.APPROX_ZERO:
            self.x /= mag_xyz
            self.y /= mag_xyz
            self.z /= mag_xyz

