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
#  - AutoDock 4.2.3 Source Code (readmap.cc)
#    http://autodock.scripps.edu

DEBUG = False

import numpy as np
import operator

class ElectrostaticMap:
    def __init__(self, filename = None, field = None):
        self.map = []
        if filename:
            self.read(filename, field)

    def read(self, filename, field):
        p_file = open(filename, 'r')

        # Parsing Header
        # --------------
        # 1. Grid Parameter File (GPF)
        #    Example: GRID_PARAMETER_FILE hsg1.gpf
        p_file.readline()
        # 2. Grid Data File
        #    Example: GRID_DATA_FILE hsg1_rigid.maps.fld
        p_file.readline()
        # 3. Macromolecule File
        #    Example: MACROMOLECULE hsg1_rigid.pdbqt
        p_file.readline()
        # 4. Spacing
        #    Example: SPACING 0.375
        p_file.readline()
        # 5. Map Size
        #    Example: NELEMENTS 60 60 60
        p_file.readline()
        # 6. Map Center
        #    Example: CENTER 2.500 6.500 -7.500
        p_file.readline()

        # Get Map Values
        # --------------
        volume = reduce(operator.mul, field.num_points1.xyz)
        # Axis allocation: map[z][y][x]
        self.map = (np.array(p_file.read().split('\n')[:volume], \
                             dtype = np.float)).reshape(field.num_points1.xyz)
        p_file.close()

    def test_print(self):
        for i in xrange(-10, -0):
            print self.map[60][60][i]

class DesolvationMap:
    def __init__(self, filename = None, field = None):
        self.map = []
        if filename:
            self.read(filename, field)

    def read(self, filename, field):
        p_file = open(filename, 'r')

        # Parsing Header
        # --------------
        # 1. Grid Parameter File (GPF)
        #    Example: GRID_PARAMETER_FILE hsg1.gpf
        p_file.readline()
        # 2. Grid Data File
        #    Example: GRID_DATA_FILE hsg1_rigid.maps.fld
        p_file.readline()
        # 3. Macromolecule File
        #    Example: MACROMOLECULE hsg1_rigid.pdbqt
        p_file.readline()
        # 4. Spacing
        #    Example: SPACING 0.375
        p_file.readline()
        # 5. Map Size
        #    Example: NELEMENTS 60 60 60
        p_file.readline()
        # 6. Map Center
        #    Example: CENTER 2.500 6.500 -7.500
        p_file.readline()

        # Get Map Values
        # --------------
        volume = reduce(operator.mul, field.num_points1.xyz)
        # Axis allocation: map[z][y][x]
        self.map = (np.array(p_file.read().split('\n')[:volume], \
                             dtype = np.float)).reshape(field.num_points1.xyz)
        p_file.close()

    def test_print(self):
        for i in xrange(-10, -0):
            print self.map[60][60][i]

class AtomTypeMap:
    def __init__(self, filename = None, field = None):
        self.map = []
        if filename:
            self.read(filename, field)

    def read(self, filename, field):
        p_file = open(filename, 'r')

        # Parsing Header
        # --------------
        # 1. Grid Parameter File (GPF)
        #    Example: GRID_PARAMETER_FILE hsg1.gpf
        p_file.readline()
        # 2. Grid Data File
        #    Example: GRID_DATA_FILE hsg1_rigid.maps.fld
        p_file.readline()
        # 3. Macromolecule File
        #    Example: MACROMOLECULE hsg1_rigid.pdbqt
        p_file.readline()
        # 4. Spacing
        #    Example: SPACING 0.375
        p_file.readline()
        # 5. Map Size
        #    Example: NELEMENTS 60 60 60
        p_file.readline()
        # 6. Map Center
        #    Example: CENTER 2.500 6.500 -7.500
        p_file.readline()

        # Get Map Values
        # --------------
        volume = reduce(operator.mul, field.num_points1.xyz)
        # Axis allocation: map[z][y][x]
        self.map = (np.array(p_file.read().split('\n')[:volume], \
                             dtype = np.float)).reshape(field.num_points1.xyz)
        p_file.close()

    def test_print(self):
        for i in xrange(-10, -0):
            print self.map[60][60][i]

if DEBUG:
    import Parameters
    field = Parameters.Field()
    field.read("./Parameters/hsg1_rigid.maps.fld")
    map = ElectrostaticMap()
    map.read("./Maps/hsg1_rigid.e.map", field)
    map.test_print()
