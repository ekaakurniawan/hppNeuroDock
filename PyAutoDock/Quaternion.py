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
#  - http://en.wikipedia.org/wiki/Quaternion
#  - Lecture 9 on Quaternion from Computer Graphics, Fall 2009 by Kenneth Joy
#    https://itunes.apple.com/us/itunes-u/computer-graphics-fall-2009/id457893733
#  - http://www.euclideanspace.com/maths/algebra/realNormedAlgebra/quaternions
#  - AutoDock 4.2.3 Source Code
#    http://autodock.scripps.edu


import math
import Constants as const
from Axis3 import *

# Quaternion: a + bi + cj + dk
class Quaternion:
    def __init__(self, a = 1.0, b = 0.0, c = 0.0, d = 0.0):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    def __repr__(self):
        return "Quaternion: %.4f + %.4fi + %.4fj + %.4fk" % \
               (self.a, self.b, self.c, self.d)

    def __copy__(self):
        q = Quaternion()
        q.a = self.a
        q.b = self.b
        q.c = self.c
        q.d = self.d
        return q

    copy = __copy__

    # This should produce a uniformly distributed quaternion, according to
    # Shoemake, Graphics Gems III.6, pp.124-132, "Uniform Random Rotations",
    # published by Academic Press, Inc., (1992)
    # Generates a unformly-distributed random quaternion (UDQ)
    def uniform(self, rng):
        # strict Shoemake version
        x0 = rng.zero_to_one()
        r1 = math.sqrt(1.0 - x0)
        r2 = math.sqrt(x0)

        t1 = rng.zero_to_2pi()
        self.b = math.sin(t1) * r1
        self.c = math.cos(t1) * r1

        t2 = rng.zero_to_2pi()
        self.d = math.sin(t2) * r2
        self.a = math.cos(t2) * r2

    def __mul__(self, other):
        Aa = self.a
        Ab = self.b
        Ac = self.c
        Ad = self.d
        Ba = other.a
        Bb = other.b
        Bc = other.c
        Bd = other.d
        q = Quaternion()
        q.a = -Ab * Bb - Ac * Bc - Ad * Bd + Aa * Ba
        q.b =  Ab * Ba + Ac * Bd - Ad * Bc + Aa * Bb
        q.c = -Ab * Bd + Ac * Ba + Ad * Bb + Aa * Bc
        q.d =  Ab * Bc - Ac * Bb + Ad * Ba + Aa * Bd
        return q

    def __imul__(self, other):
        Aa = self.a
        Ab = self.b
        Ac = self.c
        Ad = self.d
        Ba = other.a
        Bb = other.b
        Bc = other.c
        Bd = other.d
        self.a = -Ab * Bb - Ac * Bc - Ad * Bd + Aa * Ba
        self.b =  Ab * Ba + Ac * Bd - Ad * Bc + Aa * Bb
        self.c = -Ab * Bd + Ac * Ba + Ad * Bb + Aa * Bc
        self.d =  Ab * Bc - Ac * Bb + Ad * Ba + Aa * Bd
        return self

    def __abs__(self):
        return math.sqrt(self.a * self.a + \
                         self.b * self.b + \
                         self.c * self.c + \
                         self.d * self.d)

    magnitude = __abs__

    def normalize(self):
        mag = self.magnitude()
        if mag > const.APPROX_ZERO:
            self.a /= mag
            self.b /= mag
            self.c /= mag
            self.d /= mag

    def conjugate(self):
        self.b = -self.b
        self.c = -self.c
        self.d = -self.d

    def identity(self):
        self.a = 1.0
        self.b = 0.0
        self.c = 0.0
        self.d = 0.0

    # Returns angle and axis
    # By convention, angle is in radians ranging from -pi to pi
    # Example:
    #  - Quaternion (a, b, c, d) = (0.707, -0.240, -0.665, 0.000)
    #    Angle-axis (ang, x, y, z) = (90.0deg (1.57rad), -0.340, -0.940, 0.000)
    def get_angle_axis(self):
        if (self.a < -1.0) & (self.a > 1.0):
            self.normalize()
        
        angle = 2.0 * math.acos(self.a)
        # Converting radians ranging from -pi to pi
        if angle > const.PI:
            angle -= const.TWOPI

        # Normalize axis
        axis = Axis3()
        s = math.sqrt(1.0 - self.a * self.a)
        if s < const.APPROX_ZERO:
            axis.x = 1.0
            axis.y = 0.0
            axis.z = 0.0
        else:
            axis.x = self.b / s
            axis.y = self.c / s
            axis.z = self.d / s
            axis.normalize()

        return angle, axis

    def set_angle_axis(self, angle, axis):
        half_angle = angle / 2
        self.a = math.cos(half_angle)

        axis.normalize()
        s = math.sin(half_angle)
        self.b = axis.x * s
        self.c = axis.y * s
        self.d = axis.z * s

