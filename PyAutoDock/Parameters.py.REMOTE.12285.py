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

class Field:
    def __init__(self):
        self.spacing = 0.0
        self.num_points = [0, 0, 0]

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
        self.num_points = [int(element) for element in line.split(" ")[1:4]]
        


        p_file.close()

    #bar - start
    def test_print(self):
        for i in xrange(-10, -0):
            print self.map[60][60][i]
    #bar - stop

map = ElectrostaticMap()
map.read("./Maps/hsg1_rigid.e.map")
map.test_print()

