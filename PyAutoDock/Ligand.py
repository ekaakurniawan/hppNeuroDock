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
#  - AutoDock 4.2.3 Source Code (readPDBQT.cc, mkTorTree.cc)
#    http://autodock.scripps.edu

from Atom import Atom, Branch
from Axis3 import Axis3

class Ligand:
    def __init__(self):
        self.atoms = []
        self.branches = []
        self.root = None
        # Central rotation point of ligand
        self.about = Axis3(0.0, 0.0, 0.0)

    def read_pdbqt(self, filename):
        with open(filename, 'r') as p_file:
            branch_stack = []
            current_branch = None
            # Like atom ID, branch ID starts from 1.
            branch_id = 1
            for line in p_file:
                # HETATM
                if line.startswith("HETATM"):
                    data = line.split()
                    atom_id = int(data[1])
                    tcoord = Axis3(float(data[6]), \
                                   float(data[7]), \
                                   float(data[8]))
                    current_branch = branch_stack[-1]
                    atom = Atom(atom_id, data[12], tcoord, float(data[11]), \
                                current_branch)
                    self.atoms.append(atom)
                    # Update atom id into current active branches
                    # Note: First branch is root but when collecting atoms, it
                    # is treated as a branch)
                    for branch in branch_stack:
                        if atom_id != branch.link_id:
                            branch.all_atom_ids.append(atom_id)
                    if atom_id != branch_stack[-1].link_id:
                        branch_stack[-1].atom_ids.append(atom_id)
                # ROOT
                elif line.startswith("ROOT"):
                    # Ligand has only one root. Always set the ID to 1.
                    branch = Branch(branch_id, None, None, [], [], None, [])
                    self.root = branch
                    # Push root into branch_stack
                    branch_stack.append(branch)
                    branch_id += 1
                # ENDROOT
                elif line.startswith("ENDROOT"):
                    branch_stack.pop()
                # BRANCH
                elif line.startswith("BRANCH"):
                    anchor_id, link_id = [int(x) for x in line.split()[1:]]
                    # Get parent branch from the stack (which is the last one).
                    # If it empty, it means this branch directly linked to the
                    # root.
                    if not branch_stack:
                        parent_branch = self.root
                    else:
                        parent_branch = branch_stack[-1]
                    branch = Branch(branch_id, anchor_id, link_id, [], \
                                    [], parent_branch, [])
                    # Now the parent branch has this branch as the child
                    parent_branch.children.append(branch)
                    self.branches.append(branch)
                    # Push active branch into branch_stack
                    branch_stack.append(branch)
                    branch_id += 1
                # ENDBRANCH
                elif line.startswith("ENDBRANCH"):
                    # Pop inactive branch from branch_stack
                    branch_stack.pop()

    def update_tcoord_model(self, filename):
        with open(filename, 'r') as p_file:
            for line in p_file:
                # Atom
                if line.startswith("ATOM"):
                    data = line.split()
                    if data[3] != "IND": continue
                    self.atoms[int(data[1]) - 1].tcoord.xyz = \
                        [float(coord) for coord in data[6:9]]

    def get_atom_tcoords(self):
        tcoords = []
        for atom in self.atoms:
            tcoords.append(atom.tcoord)
        return tcoords

    def set_atom_tcoords(self, tcoords):
        for i, tcoord in enumerate(tcoords):
            self.atoms[i].tcoord = tcoord

    def __repr__(self):
        ret = "Atoms:\n"
        for atom in self.atoms:
            branch_parent_id = None
            try:
                branch_parent_id = atom.branch.parent.id
            except:
                pass
            ret += "%2s: %2s - %8.3f, %8.3f, %8.3f %4s %2s\n" % \
                   (atom.id, \
                    atom.type, \
                    atom.tcoord.x, \
                    atom.tcoord.y, \
                    atom.tcoord.z,
                    branch_parent_id, \
                    atom.branch.id)
        ret += "\nRoot Information:\n"
        ret += "%2s: %2s - %2s %s\n" % (self.root.id, \
                                        self.root.anchor_id, \
                                        self.root.link_id, \
                                        self.root.atom_ids)
        ret += "\nBranches Information:\n"
        for branch in self.branches:
            ret += "%2s: %2s - %2s %s\n" % (branch.id, \
                                            branch.anchor_id, \
                                            branch.link_id, \
                                            branch.atom_ids)
            ret += "    %s\n" % (branch.all_atom_ids)
        return ret

#bar - start
#l = Ligand()
#l.read_pdbqt("./Inputs/ind.pdbqt")
#print l
#bar - stop
