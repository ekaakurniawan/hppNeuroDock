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
#  - Root-mean-square deviation (bioinformatics)
#    http://en.wikipedia.org/wiki/Root-mean-square_deviation_(bioinformatics)

# Example:
# >>> from Analysis import Analysis
# >>> print Analysis.rmsd_pdbqts("./Results/ind_autodock_-15.66.pdbqt", "./Results/ind_pyneurodock_-14.27.pdbqt")

from Ligand import Ligand
from Protein import Protein
import numpy as np

class Analysis:
    # RMSD (Root-mean Square Deviation)
    @staticmethod
    def rmsd_pdbqts(pdbqt_file1, pdbqt_file2):
        ligand1 = Ligand()
        ligand1.read_pdbqt(pdbqt_file1)
        ligand2 = Ligand()
        ligand2.read_pdbqt(pdbqt_file2)
        return Analysis.rmsd_ligands(ligand1, ligand2)

    @staticmethod
    def rmsd_ligands(ligand1, ligand2):
        ttl_atoms = len(ligand1.atoms)
        tcoords1 = ligand1.get_atom_tcoords_in_numpy()
        tcoords2 = ligand2.get_atom_tcoords_in_numpy()
        return Analysis.rmsd_tcoords(tcoords1, tcoords2, ttl_atoms)

    @staticmethod
    def rmsd_tcoords(tcoords1, tcoords2, ttl_atoms):
        return np.sqrt(np.sum((tcoords1 - tcoords2) ** 2) / ttl_atoms)

    @staticmethod
    def rmsd_pdb(pdb_file1, pdb_file2):
        protein1 = Protein()
        protein1.read_pdb(pdb_file1)
        protein2 = Protein()
        protein2.read_pdb(pdb_file2)
        return Analysis.rmsd_proteins(protein1, protein2)

    @staticmethod
    def rmsd_proteins(protein1, protein2):
        ttl_atoms = len(protein1.rigid_atoms)
        tcoords1 = protein1.get_rigid_atom_tcoords_in_numpy()
        tcoords2 = protein2.get_rigid_atom_tcoords_in_numpy()
        return Analysis.rmsd_tcoords(tcoords1, tcoords2, ttl_atoms)

