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
#  - AutoDock 4.2.3 Source Code (readfield.cc)
#    http://autodock.scripps.edu

DEBUG = False

from Axis3 import *

class Field:
    def __init__(self, filename = None):
        self.spacing = 0.0
        self.num_points = Axis3(0, 0, 0)
        self.num_points1 = Axis3(0, 0, 0)
        self.center = Axis3(0.0, 0.0, 0.0)
        self.lo = Axis3(0.0, 0.0, 0.0)
        self.hi = Axis3(0.0, 0.0, 0.0)
        if filename:
            self.read(filename)

    def read(self, filename):
        p_file = open(filename, 'r')

        # Spacing
        while True:
            line = p_file.readline()
            if ("#SPACING" in line):
                break
        self.spacing = float(line.split(" ")[1])

        # Number of Elements
        line = p_file.readline()
        self.num_points.axis = [int(element) for element in line.split(" ")[1:4]]
        # AutoGrid always adds number of elements with one for each dimension
        # as the center point
        self.num_points1.axis = [element + 1 for element in self.num_points.axis]
        
        # Center
        line = p_file.readline()
        self.center.axis = [float(center) for center in line.split(" ")[1:4]]

        # Macromolecule
        line = p_file.readline()

        # Grid Parameter File (gpf)
        line = p_file.readline()
        
        p_file.close()

        # Get minimum and maximum value of x, y, and z
        half_range = (self.num_points.x / 2) * self.spacing
        self.lo.x = self.center.x - half_range
        self.hi.x = self.center.x + half_range
        half_range = (self.num_points.y / 2) * self.spacing
        self.lo.y = self.center.y - half_range
        self.hi.y = self.center.y + half_range
        half_range = (self.num_points.z / 2) * self.spacing
        self.lo.z = self.center.z - half_range
        self.hi.z = self.center.z + half_range

    def test_print(self):
        print "Spacing     : %s" % self.spacing
        print "nPoints     : %s" % self.num_points.axis
        print "nPoints + 1 : %s" % self.num_points1.axis
        print "Center      : %s" % self.center.axis
        print "Lo          : %s" % self.lo.axis
        print "Hi          : %s" % self.hi.axis

if DEBUG:
    field = Field()
    field.read("./Parameters/hsg1_rigid.maps.fld")
    field.test_print()
