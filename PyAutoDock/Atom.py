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
#  - AutoDock 4.2.3 Source Code (mdist.h, nonbonds.cc, weedbonds.cc,
#    parameters.h, read_parameter_library.cc)
#    http://autodock.scripps.edu
#  - Covalent Radius
#    http://en.wikipedia.org/wiki/Covalent_radius

import math
from Axis3 import Axis3

class Atom:
    # For atom type, please refer to bond_index column inside AD4.1_bound.dat
    # or AD4_parameters.dat file
    NUM_ATOM_TYPE = 7
    C, N, O, H, XX, P, S = range(NUM_ATOM_TYPE) # range starts from 0
    ATOM_TYPE = {'C' : C,
                 'A' : 0,
                 'N' : N,
                 'NA': 1,
                 'O' : O,
                 'OA': 2,
                 'H' : H,
                 'HD': 3,
                 'XX': XX,
                 'P' : P,
                 'S' : S,
                }

    def __init__(self, id = 0, type = '?', tcoord = Axis3(0.0, 0.0, 0.0), \
                 charge = 0.0, branch = None):
        self.id = id
        self.type = type
        self.tcoord = tcoord
        self.charge = charge
        self.branch = branch

class Branch:
    def __init__(self, id = 0, anchor_id = 0, link_id = 0, atom_ids = [], \
                 all_atom_ids = [], parent = None, children = []):
        self.id = id
        # Atom ID at parent branch
        self.anchor_id = anchor_id
        # Atom ID links to anchor atom
        self.link_id = link_id
        # Atom IDs in the branch except anchor and link atoms
        self.atom_ids = atom_ids
        # All atom IDs in the branch except anchor and link atoms; and atom IDs
        # of all sub-branches
        self.all_atom_ids = all_atom_ids
        # Parent and child branches of this branch
        self.parent = parent
        self.children = children

# Atomic Bond
class Bond:
    # Covalent bond distance tolerance
    DISTANCE_TOLERANCE = 0.1

    # Desolvation Parameters
    DESOLVATION_SIGMA = 3.6
    DESOLVATION = 0.01097

    # ELECSCALE converts between CGS units and SI units; see, e.g. p 254, 
    # "Molecular Modeling and Simulation", by Tamar Schlick, Springer.
    # Units of ELECSCALE are (Kcal/mol ) * (Angstrom / esu^2)
    # and this allows us to use distances in  Angstroms and charges in esu...
    ELECSCALE = 332.06363

    # Hydrogen Bonding Types
    NUM_H_BOND_TYPE = 6
    NON, DS, D1, AS, A1, A2 = range(NUM_H_BOND_TYPE) # range starts from 0
    H_BOND_TYPE = {'NON' : NON, # None
                   'DS'  : DS,  # Spherical donor
                   'D1'  : D1,  # Directional donor
                   'AS'  : AS,  # Spherical acceptor
                   'A1'  : A1,  # Acceptor of 1 directional hydrogen bond
                   'A2'  : A2,  # Acceptor of 2 directional hydrogen bond
                  }

    class EnergyParameter:
        def __init__(self, rij = 0.0, epsij = 0.0, vol = 0.0, solpar = 0.0, \
                     rij_hb = 0.0, epsij_hb = 0.0, hbond = '?', rec_index = 0, \
                     map_index = 0, bond_index = 0):
            # Lennard-Jones equilibrium separation
            self.rij = rij
            # Lennard-Jones energy well-depth
            self.epsij = epsij
            # Solvation volume
            self.vol = vol
            # Solvation parameter
            self.solpar = solpar
            # 12-10 Lennard-Jones equilibrium separation
            self.rij_hb = rij_hb
            # 12-10 Lennard-Jones energy well-depth
            self.epsij_hb = epsij_hb
            # Hydrogen bonding type
            self.hbond = hbond
            # Used to set up receptor atom_types
            self.rec_index = rec_index
            # Used to set up map atom_types
            self.map_index = map_index
            # Used to set up bonds; corresponds to the enum in mdist.h
            self.bond_index = bond_index

    class NonBond:
        def __init__(self, atom1 = 0, atom_type1 = '?', \
                     atom2 = 0, atom_type2 = '?', \
                     non_bond_type = 1, desolv = 0.0, q1q2 = 0.0):
            self.atom1 = atom1
            self.atom_type1 = atom_type1
            self.atom2 = atom2
            self.atom_type2 = atom_type2
            # Bond Types:
            #  1: Non-bond
            #  0: 1-1, 1-2, and 1-3 interactions
            #  4: 1-4 interaction
            self.non_bond_type = non_bond_type
            # Product of the partial atomic charges
            self.q1q2 = q1q2

    def __init__(self):
        # By default, 1-4 interactions is disabled
        self.include_1_4_interactions = False

        # Free energy coefficients for Van der Waals term
        self.fec_vdw = 0.0
        # Free energy coefficients for hydrogen bonding term
        self.fec_hbond = 0.0
        # Free energy coefficients for electrostatic term
        self.fec_estat = 0.0
        # Free energy coefficients for desolvation term
        self.fec_desolv = 0.0
        # Free energy coefficients for torsional term
        self.fec_tors = 0.0

        # Energy parameter of different atom types
        self.e_parms = {}

    # Read bonding parameter file
    def read(self, filename):
        with open("./Parameters/" + filename, 'r') as p_file:
            for line in p_file:
                if line.startswith("FE_coeff_vdW"):
                    self.fec_vdw = float(line.split()[1])

                if line.startswith("FE_coeff_hbond"):
                    self.fec_hbond = float(line.split()[1])

                if line.startswith("FE_coeff_estat"):
                    self.fec_estat = float(line.split()[1])

                if line.startswith("FE_coeff_desolv"):
                    self.fec_desolv = float(line.split()[1])

                if line.startswith("FE_coeff_tors"):
                    self.fec_tors = float(line.split()[1])

                if line.startswith("atom_par"):
                    data = line.split()
                    atom_type  = data[1]
                    rij        = float(data[2])
                    epsij      = float(data[3]) * self.fec_vdw
                    vol        = float(data[4])
                    solpar     = float(data[5])
                    rij_hb     = float(data[6])
                    epsij_hb   = float(data[7]) * self.fec_hbond
                    hbond_type = int(data[8])
                    rec_index  = int(data[9])
                    map_index  = int(data[10])
                    bond_index = int(data[11])

                    H_BOND_TYPE_AD4 = {0: 'NON',
                                       1: 'DS',
                                       2: 'D1',
                                       3: 'AS',
                                       4: 'A1',
                                       5: 'A2'}
                    hbond = H_BOND_TYPE_AD4.get(hbond_type, 'NON')
                    
                    e_parm = self.EnergyParameter(rij, epsij, vol, solpar, \
                                                  rij_hb, epsij_hb, hbond, \
                                                  rec_index, map_index, \
                                                  bond_index)
                    self.e_parms[atom_type] = e_parm

    # Get natural observation of the covalent bonding range of atom to atom
    # distance
    def calc_minmax_distance(self):
        minmax_distance = [[[0, 0] for x in xrange(Atom.NUM_ATOM_TYPE)] \
                           for x in xrange(Atom.NUM_ATOM_TYPE)]

        def set_minmax_distance(atom1, atom2, val):
            minmax_distance[atom1][atom2][0] = val[0] - self.DISTANCE_TOLERANCE
            minmax_distance[atom2][atom1][0] = minmax_distance[atom1][atom2][0]

            minmax_distance[atom1][atom2][1] = val[1] + self.DISTANCE_TOLERANCE
            minmax_distance[atom2][atom1][1] = minmax_distance[atom1][atom2][1]

        # Following values, unless otherwise stated, are taken from AutoDock,
        # taken from "Handbook of Chemistry and Physics" 44th edition
        set_minmax_distance(Atom.C,  Atom.C,  [1.20,    1.545])  # mindist[C][C] = 1.20, p. 3510 ; maxdist[C][C] = 1.545, p. 3511
        set_minmax_distance(Atom.C,  Atom.N,  [1.1,     1.479])  # mindist[C][N] = 1.1, p. 3510 ; maxdist[C][N] = 1.479, p. 3511
        set_minmax_distance(Atom.C,  Atom.O,  [1.15,    1.47])   # mindist[C][O] = 1.15, p. 3510 ; maxdist[C][O] = 1.47, p. 3512
        set_minmax_distance(Atom.C,  Atom.H,  [1.022,   1.12])   # p. 3518, p. 3517
        set_minmax_distance(Atom.C,  Atom.XX, [0.9,     1.545])  # mindist[C][XX] = 0.9, AutoDock 3 defaults ; maxdist[C][XX] = 1.545, p. 3511
        set_minmax_distance(Atom.C,  Atom.P,  [1.85,    1.89])   # mindist[C][P] = 1.85, p. 3510 ; maxdist[C][P] = 1.89, p. 3510
        set_minmax_distance(Atom.C,  Atom.S,  [1.55,    1.835])  # mindist[C][S] = 1.55, p. 3510 ; maxdist[C][S] = 1.835, p. 3512
        set_minmax_distance(Atom.N,  Atom.N,  [1.0974,  1.128])  # mindist[N][N] = 1.0974, p. 3513 ; maxdist[N][N] = 1.128, p. 3515
        set_minmax_distance(Atom.N,  Atom.O,  [1.0619,  1.25])   # mindist[N][O] = 1.0975, p. 3515 ; maxdist[N][O] = 1.128, p. 3515
        set_minmax_distance(Atom.N,  Atom.H,  [1.004,   1.041])  # mindist[N][H] = 1.004, p. 3516 ; maxdist[N][H] = 1.041, p. 3515
        set_minmax_distance(Atom.N,  Atom.XX, [0.9,     1.041])  # mindist[N][XX] = 0.9, AutoDock 3 defaults ; maxdist[N][XX] = 1.041, p. 3515
        set_minmax_distance(Atom.N,  Atom.P,  [1.4910,  1.4910]) # mindist[N][P] = 1.4910, p. 3515 ; maxdist[N][P] = 1.4910, p. 3515
        set_minmax_distance(Atom.N,  Atom.S,  [1.58,    1.672])  # mindist[N][S] = 1.58, 1czm.pdb sulfonamide ; maxdist[N][S] = 1.672, J. Chem. SOC., Dalton Trans., 1996, Pages 4063-4069
        set_minmax_distance(Atom.O,  Atom.O,  [1.208,   1.51])   # p.3513, p.3515
        set_minmax_distance(Atom.O,  Atom.H,  [0.955,   1.0289]) # mindist[O][H] = 0.955, p. 3515 ; maxdist[O][H] = 1.0289, p. 3515
        set_minmax_distance(Atom.O,  Atom.XX, [0.955,   2.1])    # AutoDock 3 defaults
        set_minmax_distance(Atom.O,  Atom.P,  [1.36,    1.67])   # mindist[O][P] = 1.36, p. 3516 ; maxdist[O][P] = 1.67, p. 3517
        set_minmax_distance(Atom.O,  Atom.S,  [1.41,    1.47])   # p. 3517, p. 3515
        set_minmax_distance(Atom.H,  Atom.H,  [100.0,  -100.0])  # impossible values to prevent such bonds from forming.
        set_minmax_distance(Atom.H,  Atom.XX, [0.9,     1.5])    # AutoDock 4 defaults
        set_minmax_distance(Atom.H,  Atom.P,  [1.40,    1.44])   # mindist[H][P] = 1.40, p. 3515 ; maxdist[H][P] = 1.44, p. 3515
        set_minmax_distance(Atom.H,  Atom.S,  [1.325,   1.3455]) # mindist[H][S] = 1.325, p. 3518 ; maxdist[H][S] = 1.3455, p. 3516
        set_minmax_distance(Atom.XX, Atom.XX, [0.9,     2.1])    # AutoDock 3 defaults
        set_minmax_distance(Atom.XX, Atom.P,  [0.9,     2.1])    # AutoDock 3 defaults
        set_minmax_distance(Atom.XX, Atom.S,  [1.325,   2.1])    # mindist[XX][S] = 1.325, p. 3518 ; maxdist[XX][S] = 2.1, AutoDock 3 defaults
        set_minmax_distance(Atom.P,  Atom.P,  [2.18,    2.23])   # mindist[P][P] = 2.18, p. 3513 ; maxdist[P][P] = 2.23, p. 3513
        set_minmax_distance(Atom.P,  Atom.S,  [1.83,    1.88])   # mindist[P][S] = 1.83, p. 3516 ; maxdist[P][S] = 1.88, p. 3515
        set_minmax_distance(Atom.S,  Atom.S,  [2.03,    2.05])   # mindist[S][S] = 2.03, p. 3515 ; maxdist[S][S] = 2.05, p. 3515

        return minmax_distance

    # Construct atom-to-atom covalent bonding matrix based on their distance for
    # ligand and protein individually
    def construct_bond_matrix(self, atoms, minmax_distance):
        total_atoms = len(atoms)
        # Construct an array of integer atom type from character atom type
        atom_types = []
        for i in xrange(total_atoms):
            atom_types.append(Atom.ATOM_TYPE[atoms[i].type])

        bond_matrix = [[] for i in xrange(total_atoms)]
        distances = [[] for i in xrange(total_atoms)]
        for i in xrange(total_atoms - 1):
            # Calculate atom-to-atom distance
            for j in xrange(i + 1, total_atoms):
                dx = atoms[i].tcoord.x - atoms[j].tcoord.x
                dy = atoms[i].tcoord.y - atoms[j].tcoord.y
                dz = atoms[i].tcoord.z - atoms[j].tcoord.z
                distances[j] = math.sqrt((dx * dx) + (dy * dy) + (dz * dz))
            sorted_idx = range(i + 1)
            sorted_idx += sorted(range(i + 1, total_atoms), \
                                 key=lambda k: distances[k])

            # Scan from the shortest atom-to-atom distance and check if it
            # within tolerable bonding range
            for j in xrange(i + 1, total_atoms):
                k = sorted_idx[j]
                distance = distances[k]
                min_distance = minmax_distance[atom_types[i]][atom_types[k]][0]
                max_distance = minmax_distance[atom_types[i]][atom_types[k]][1]
                if (distance >= min_distance) and (distance <= max_distance):
                    bond_matrix[i].append(k)
                    bond_matrix[k].append(i)

        return bond_matrix

    # Construct atom-to-atom non-bonding matrix
    def construct_non_bond_matrix(self, total_atoms):
        # Set initial non-bonding matrix values all to 1. As the detection
        # progressing, they will be set to:
        #  - 0 to be ignored (bonded), or to
        #  - 4 for another type of non-bonding atoms, 1-4 interactions
        non_bond_matrix = [[1 for i in xrange(total_atoms)] \
                           for j in xrange(total_atoms)]
        return non_bond_matrix

    # Weed out covalent bond by considering 1-1, 1-2, 1-3 and/or 1-4
    # interactions
    def weed_covalent_bond(self, bond_matrix, non_bond_matrix):
        total_atoms = len(bond_matrix)

        # If we include 1-4 interactions, mark them separately (4) otherwise
        # mark them as the rest, which is 0
        if self.include_1_4_interactions:
            mark_1_4 = 4
        else:
            mark_1_4 = 0
        # 1-1, 1-2, 1-3 and 1-4 interactions
        for i in xrange(total_atoms):
            non_bond_matrix[i][i] = 0                       # 1-1 interactions
            for j in bond_matrix[i]:
                non_bond_matrix[i][j] = 0                   # 1-2 interactions
                non_bond_matrix[j][i] = 0
                for k in bond_matrix[j]:
                    non_bond_matrix[i][k] = 0               # 1-3 interactions
                    non_bond_matrix[k][i] = 0
                    for l in bond_matrix[k]:
                        non_bond_matrix[i][l] = mark_1_4    # 1-4 interactions
                        non_bond_matrix[l][i] = mark_1_4

        return non_bond_matrix

    # For internal energy calculation, weed out:
    # - rigidly bonded root atoms
    # - anchor-link atoms
    # - link atoms in a same rigid body. Also, weed out link atom and the atoms
    #   in a same rigid body. These are considered 1-3 interactions.
    # Applicable for both ligand and protein.
    def weed_rigid_bond(self, non_bond_matrix, ligand, protein):
        # Starting index for protein atom IDs
        p_idx = len(ligand.atoms)
        
        # Weed out rigidly bonded root atoms
        for atom_id_i in ligand.root.atom_ids:
            for atom_id_j in ligand.root.atom_ids:
                non_bond_matrix[atom_id_i - 1][atom_id_j - 1] = 0
        for atom_id_i in protein.root.atom_ids:
            for atom_id_j in protein.root.atom_ids:
                non_bond_matrix[p_idx + atom_id_i - 1] \
                               [p_idx + atom_id_j - 1] = 0

        # Weed out rigidly bonded atoms in a same branch
        for branch in ligand.branches:
            for atom_id_i in branch.atom_ids:
                non_bond_matrix[atom_id_i - 1][branch.anchor_id - 1] = 0
                non_bond_matrix[branch.anchor_id - 1][atom_id_i - 1] = 0
                non_bond_matrix[atom_id_i - 1][branch.link_id - 1] = 0
                non_bond_matrix[branch.link_id - 1][atom_id_i - 1] = 0
                for atom_id_j in branch.atom_ids:
                    non_bond_matrix[atom_id_i - 1][atom_id_j - 1] = 0
        for branch in protein.flex_branches:
            for atom_id_i in branch.atom_ids:
                non_bond_matrix[p_idx + atom_id_i - 1] \
                               [p_idx + branch.anchor_id - 1] = 0
                non_bond_matrix[p_idx + branch.anchor_id - 1] \
                               [p_idx + atom_id_i - 1] = 0
                non_bond_matrix[p_idx + atom_id_i - 1] \
                               [p_idx + branch.link_id - 1] = 0
                non_bond_matrix[p_idx + branch.link_id - 1]\
                               [p_idx + atom_id_i - 1] = 0
                for atom_id_j in branch.atom_ids:
                    non_bond_matrix[p_idx + atom_id_i - 1] \
                                   [p_idx + atom_id_j - 1] = 0

        # Weed out anchor-link atoms
        for branch in ligand.branches:
            non_bond_matrix[branch.link_id - 1][branch.anchor_id - 1] = 0
        for branch in protein.flex_branches:
            non_bond_matrix[p_idx + branch.link_id - 1] \
                           [p_idx + branch.anchor_id - 1] = 0

        # Weed out link atoms in a same rigid body.
        # Also, weed out link atom and the atoms in a same rigid body.
        # These are considered 1-3 interactions.
        # Note: - Rigid body for a link atom is the parent branch as they are
        #         stick together
        #       - Rigid body of the rest of the atoms (except link atom) is the
        #         the branch they attach to
        for branch1 in ligand.branches:
            for branch2 in ligand.branches:
                if branch1.parent.id == branch2.parent.id:
                    non_bond_matrix[branch1.link_id - 1] \
                                   [branch2.link_id - 1] = 0
            for atom in ligand.atoms:
                if atom.branch.id == branch1.parent.id:
                    non_bond_matrix[atom.id - 1][branch1.link_id - 1] = 0
                    non_bond_matrix[branch1.link_id - 1][atom.id - 1] = 0
        for branch1 in protein.flex_branches:
            for branch2 in protein.flex_branches:
                if branch1.parent.id == branch2.parent.id:
                    non_bond_matrix[p_idx + branch1.link_id - 1] \
                                   [p_idx + branch2.link_id - 1] = 0
            for atom in protein.flex_atoms:
                if atom.branch.id == branch1.parent.id:
                    non_bond_matrix[p_idx + atom.id - 1] \
                                   [p_idx + branch1.link_id - 1] = 0
                    non_bond_matrix[p_idx + branch1.link_id - 1] \
                                   [p_idx + atom.id - 1] = 0

        return non_bond_matrix

    def convert_non_bond_matrix_to_list(self, non_bond_matrix, ligand, protein):
        ligand_len = len(ligand.atoms)
        protein_len = len(protein.flex_atoms)
        scale = self.ELECSCALE * self.fec_estat
        # Intramolecular non-bond for ligand
        non_bond_ligand = []
        for i in xrange(ligand_len):
            for j in xrange(i + 1, ligand_len):
                if (non_bond_matrix[i][j] != 1) and \
                   (non_bond_matrix[i][j] != 4):
                    continue
                atom_type1 = ligand.atoms[i].type
                atom_type2 = ligand.atoms[j].type
                desolv = (self.e_parms[atom_type2].vol * \
                             (self.e_parms[atom_type1].solpar + \
                                 (self.DESOLVATION * ligand.atoms[i].charge)) + \
                          self.e_parms[atom_type1].vol * \
                             (self.e_parms[atom_type2].solpar + \
                                 (self.DESOLVATION * ligand.atoms[j].charge)))
                q1q2 = scale * ligand.atoms[i].charge * ligand.atoms[j].charge

                nbi = self.NonBond(ligand.atoms[i].id, atom_type1, \
                                   ligand.atoms[j].id, atom_type2, \
                                   non_bond_matrix[i][j], desolv, q1q2)
                non_bond_ligand.append(nbi)
        
        # Intramolecular non-bond for ligand and receptor
        non_bond_ligand_receptor = []
        for i in xrange(ligand_len):
            for j in xrange(protein_len):
                if (non_bond_matrix[i][j + ligand_len] != 1) and \
                   (non_bond_matrix[i][j + ligand_len] != 4):
                        continue
                atom_type1 = ligand.atoms[i].type
                atom_type2 = protein.flex_atoms[j].type
                desolv = (self.e_parms[atom_type2].vol * \
                             (self.e_parms[atom_type1].solpar + \
                                 (self.DESOLVATION * ligand.atoms[i].charge)) + \
                          self.e_parms[atom_type1].vol * \
                             (self.e_parms[atom_type2].solpar + \
                                 (self.DESOLVATION * protein.flex_atoms[j].charge)))
                q1q2 = scale * ligand.atoms[i].charge * protein.flex_atoms[j].charge

                nbi = self.NonBond(ligand.atoms[i].id, atom_type1, \
                                   protein.flex_atoms[j].id, atom_type2, \
                                   non_bond_matrix[i][j + ligand_len], \
                                   desolv, q1q2)
                non_bond_ligand_receptor.append(nbi)

        # Intramolecular non-bond for receptor
        non_bond_receptor = []
        for i in xrange(protein_len):
            for j in xrange(i + 1, protein_len):
                if (non_bond_matrix[i + ligand_len][j + ligand_len] != 1) and \
                   (non_bond_matrix[i + ligand_len][j + ligand_len] != 4):
                        continue
                atom_type1 = protein.flex_atoms[i].type
                atom_type2 = protein.flex_atoms[j].type
                desolv = (self.e_parms[atom_type2].vol * \
                             (self.e_parms[atom_type1].solpar + \
                                 (self.DESOLVATION * protein.flex_atoms[i].charge)) + \
                          self.e_parms[atom_type1].vol * \
                             (self.e_parms[atom_type2].solpar + \
                                 (self.DESOLVATION * protein.flex_atoms[j].charge)))
                q1q2 = scale * protein.flex_atoms[i].charge * \
                       protein.flex_atoms[j].charge

                nbi = self.NonBond(protein.flex_atoms[i].id, atom_type1, \
                                   protein.flex_atoms[j].id, atom_type2, \
                                   non_bond_matrix[i + ligand_len] \
                                                  [j + ligand_len],
                                   desolv, q1q2)
                non_bond_receptor.append(nbi)

        return non_bond_ligand, non_bond_ligand_receptor, non_bond_receptor

    def __repr__(self):
        ret = "Bonding Parameters:\n"
        ret += "Free energy coefficient for Van der Waals term    = %6.3f\n" % \
                self.fec_vdw
        ret += "Free energy coefficient for hydrogen bonding term = %6.3f\n" % \
                self.fec_hbond
        ret += "Free energy coefficient for electrostatic term    = %6.3f\n" % \
                self.fec_estat
        ret += "Free energy coefficient for desolvation term      = %6.3f\n" % \
                self.fec_desolv
        ret += "Free energy coefficient for torsional term        = %6.3f\n" % \
                self.fec_tors
        ret += "\n"
        ret += "Energy Parameter of Different Atom Types:\n"
        for key in self.e_parms:
            ret += "  %-2s" % key
            ret += "  %4.2f" % self.e_parms[key].rij
            ret += "  %5.3f" % self.e_parms[key].epsij
            ret += "  %7.4f" % self.e_parms[key].vol
            ret += "  %8.5f" % self.e_parms[key].solpar
            ret += "  %3.1f" % self.e_parms[key].rij_hb
            ret += "  %3.1f" % self.e_parms[key].epsij_hb
            ret += "  %d"    % self.H_BOND_TYPE[self.e_parms[key].hbond]
            ret += "  %2d"   % self.e_parms[key].rec_index
            ret += "  %2d"   % self.e_parms[key].map_index
            ret += "  %d\n"  % self.e_parms[key].bond_index

        return ret

