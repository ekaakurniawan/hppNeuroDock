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

from Atom import Atom, Branch
from Axis3 import Axis3

class Protein:
    def __init__(self):
        self.flex_atoms = []        # flexible (rotatable) atoms
        self.flex_branches = []     # flexible (rotatable) branches

    def read_flex_pdbqt(self, filename):
        with open(filename, 'r') as p_file:
            branch_stack = []
            for line in p_file:
                # HETATM
                if line.startswith("ATOM"):
                    data = line.split()
                    atom_id = int(data[1])
                    tcoord = Axis3(float(data[6]), \
                                   float(data[7]), \
                                   float(data[8]))
                    atom = Atom(atom_id, data[12], tcoord, float(data[11]))
                    self.flex_atoms.append(atom)
                    # Update atom id into current active branches
                    for branch in branch_stack:
                        branch.atom_ids.append(atom_id)
                # BRANCH
                elif line.startswith("BRANCH"):
                    anchor_id, link_id = [int(x) for x in line.split()[1:]]
                    branch = Branch(anchor_id, link_id, [])
                    self.flex_branches.append(branch)
                    # Push active branch into branch_stack
                    branch_stack.append(branch)
                # ENDBRANCH
                elif line.startswith("ENDBRANCH"):
                    # Pop inactive branch from branch_stack
                    branch_stack.pop()

    def get_flex_atom_tcoords(self):
        tcoords = []
        for flex_atom in self.flex_atoms:
            tcoords.append(flex_atom.tcoord)
        return tcoords

    def set_flex_atom_tcoords(self, tcoords):
        for i, tcoord in enumerate(tcoords):
            self.flex_atoms[i].tcoord = tcoord

    def __repr__(self):
        ret = "Flexible (Rotatable) Atoms:\n"
        for flex_atom in self.flex_atoms:
            ret += "%2s: %2s - %8.3f, %8.3f, %8.3f\n" % (flex_atom.id, \
                                                         flex_atom.type, \
                                                         flex_atom.tcoord.x, \
                                                         flex_atom.tcoord.y, \
                                                         flex_atom.tcoord.z)
        ret += "\nBranch Information:\n"
        for flex_branch in self.flex_branches:
            ret += "%2s - %2s %s\n" % (flex_branch.anchor_id, \
                                       flex_branch.link_id, \
                                       flex_branch.atom_ids)
        return ret

#bar - start
#p = Protein()
#p.read_flex_pdbqt("./Inputs/hsg1_flex.pdbqt")
#print p
#bar - stop
