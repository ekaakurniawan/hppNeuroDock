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

import numpy as np
import operator

class AtomTypeMap:
    def __init__(self):
        self.map = None
    
    def read(self, filename):
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
        self.size = np.array(p_file.readline().split()[1:4], dtype = np.integer)
        # Map size being set is in nx, ny and nx but AutoGrid always adds one in
        # each dimension.
        self.size = self.size + [1, 1, 1]
        self.volume = reduce(operator.mul, self.size)
        # 6. Map Center
        #    Example: CENTER 2.500 6.500 -7.500
        p_file.readline()

        # Get Map Values
        # --------------
        self.map = (np.array(p_file.read().split('\n')[:self.volume], \
                             dtype = np.float)).reshape(self.size)
        p_file.close()

    #bar - start
    def test_print(self):
        for i in xrange(-10, -0):
            print self.map[60][60][i]
    #bar - stop


class DesolvationMap:
    def __init__(self):
        self.map = None

    def read(self, filename):
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
        self.size = np.array(p_file.readline().split()[1:4], dtype = np.integer)
        # Map size being set is in nx, ny and nx but AutoGrid always adds one in
        # each dimension.
        self.size = self.size + [1, 1, 1]
        self.volume = reduce(operator.mul, self.size)
        # 6. Map Center
        #    Example: CENTER 2.500 6.500 -7.500
        p_file.readline()

        # Get Map Values
        # --------------
        self.map = (np.array(p_file.read().split('\n')[:self.volume], \
                             dtype = np.float)).reshape(self.size)
        p_file.close()

    #bar - start
    def test_print(self):
        for i in xrange(-10, -0):
            print self.map[60][60][i]
    #bar - stop

class ElectrostaticMap:
    def __init__(self):
        self.map = None

    def read(self, filename):
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
        self.size = np.array(p_file.readline().split()[1:4], dtype = np.integer)
        # Map size being set is in nx, ny and nx but AutoGrid always adds one in
        # each dimension.
        self.size = self.size + [1, 1, 1]
        self.volume = reduce(operator.mul, self.size)
        # 6. Map Center
        #    Example: CENTER 2.500 6.500 -7.500
        p_file.readline()

        # Get Map Values
        # --------------
        self.map = (np.array(p_file.read().split('\n')[:self.volume], \
                             dtype = np.float)).reshape(self.size)
        p_file.close()

    #bar - start
    def test_print(self):
        for i in xrange(-10, -0):
            print self.map[60][60][i]
    #bar - stop

map = ElectrostaticMap()
map.read("./Maps/hsg1_rigid.e.map")
map.test_print()

