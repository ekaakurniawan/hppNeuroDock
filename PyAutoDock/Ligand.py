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
#  - AutoDock 4.2.3 Source Code (readPDBQT.cc)
#    http://autodock.scripps.edu

DEBUG = False

from Atom import *

class Ligand:
    def __init__(self):
        self.atoms = []

    def read_pdbqt(self, filename):
        with open(filename, 'r') as p_file:
            for line in p_file:
                # HETATM
                if "HETATM" in line:
                    data = line.split()
                    atom = Atom(data[12], \
                                [float(coord) for coord in data[6:9]], \
                                float(data[11]), abs(float(data[11])))
                    self.atoms.append(atom)

    def update_tcoord_model(self, filename):
        with open(filename, 'r') as p_file:
            for line in p_file:
                # Atom
                if "ATOM" in line:
                    data = line.split()
                    if data[3] != "IND": continue
                    self.atoms[int(data[1]) - 1].tcoord.axis = \
                        [float(coord) for coord in data[6:9]]

    def test_print(self):
        for i, atom in enumerate(self.atoms):
            print "%2s: %2s - %8.3f, %8.3f, %8.3f" % (i + 1, atom.type, \
                                                      atom.tcoord.x, \
                                                      atom.tcoord.y, \
                                                      atom.tcoord.z)

if DEBUG:
    ligand = Ligand()
    ligand.read_pdbqt("./Inputs/ind.pdbqt")
    ligand.test_print()
    ligand.update_tcoord_model("./Results/ind.model")
    ligand.test_print()
