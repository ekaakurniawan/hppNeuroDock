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
#  - AutoDock 4.2.3 Source Code (mkTorTree.cc, torsion.cc, trilinterp.cc,
#    torNorVec.cc, eintcal.cc)
#    http://autodock.scripps.edu

from Ligand import Ligand
from Protein import Protein
from Grid import Grid
from Quaternion import Quaternion
from Atom import Bond
import numpy as np
import pyopencl as cl

DEBUG = False

class DockingParameters:
    def __init__(self):
        # Calculate Internal Electrostatic Energies Flag
        self.calc_inter_elec_e = False

    def __repr__(self):
        ret = ""
        ret += "Calculate Internal Electrostatic Energies Flag         : %s" % \
            self.calc_inter_elec_e
        return ret

class Dock:
    SCALE_1_4_INTERACTIONS = 0.5

    def __init__(self):
        self.ligand = Ligand()
        self.protein = Protein()
        self.grid = Grid()
        self.bond = Bond()
        # Docking Parameters
        self.dps = DockingParameters()

        # Sorted ligand and protein branches ascendingly based on number of
        # atoms in the branch
        self.sorted_branches = []
        # Bonding lists
        self.non_bond_ligand = []
        self.non_bond_ligand_receptor = []
        self.non_bond_receptor = []

        # Binding torsional free energy
        self.torsional_energy = 0.0
        # Electrostatic
        self.elecs = []
        self.elec_total = 0.0
        # Van der Waals
        self.emaps = []
        self.emap_total = 0.0

    def get_total_torsions(self):
        ttl_torsions = 0
        ttl_torsions += len(self.ligand.branches)
        ttl_torsions += len(self.protein.flex_branches)
        return ttl_torsions

    def get_total_atoms(self):
        ttl_atoms = 0
        ttl_atoms += len(self.ligand.atoms)
        ttl_atoms += len(self.protein.flex_atoms)
        return ttl_atoms

    # Rotate rotatable branches/bonds for both ligand and protein.
    # Rotation is expected to be in radian.
    def rotate_branches(self, torsions):
        if not self.sorted_branches:
            for branch in self.ligand.branches:
                branch.molecule = 'l' # l for ligand
                self.sorted_branches.append(branch)
            for branch in self.protein.flex_branches:
                branch.molecule = 'p' # p for protein
                self.sorted_branches.append(branch)
            self.sorted_branches = sorted(self.sorted_branches, \
                                          key=lambda branch: len(branch.all_atom_ids))

        q_rotation = Quaternion()
        rot_i = 0
        for branch in self.sorted_branches:
            atoms = []
            atom_tcoords = []
            # Get the atoms from either ligand or protein based on branch
            # molecule information
            if branch.molecule == 'l':
                molecule_atoms = self.ligand.atoms
            else: # 'p'
                molecule_atoms = self.protein.flex_atoms
            # Get atom coordinates
            for atom in molecule_atoms:
                if atom.id in [branch.anchor_id] + [branch.link_id] + \
                              branch.all_atom_ids:
                    # Anchor and link atoms are not for rotation
                    if atom.id == branch.anchor_id:
                        anchor_tcoord = atom.tcoord
                        continue
                    if atom.id == branch.link_id:
                        link_tcoord = atom.tcoord
                        continue
                    atoms.append(atom)
                    atom_tcoords.append(atom.tcoord - link_tcoord)
            # Transform
            q_rotation.set_angle_axis(torsions[rot_i], \
                                      anchor_tcoord - link_tcoord)
            new_atom_tcoords = Quaternion.transform(link_tcoord, q_rotation, \
                                                    atom_tcoords)
            for i, atom in enumerate(atoms):
                atom.tcoord = new_atom_tcoords[i]

            rot_i += 1

    # Transform (translate and rotate) ligand root (whole body)
    def transform_ligand_root(self, translation, rotation):
        atom_tcoords = self.ligand.get_atom_tcoords()
        new_atom_tcoords = Quaternion.transform(translation, rotation, \
                                                atom_tcoords)
        self.ligand.set_atom_tcoords(new_atom_tcoords)

    # Set candidate binding mode
    # Return: - True: Success
    #         - False: New pose is out of grid
    def set_pose(self, translation, rotation, torsions):
        self.rotate_branches(torsions)
        self.transform_ligand_root(translation, rotation)
        if self.check_out_of_grid():
            return False
        else:
            return True

    # Set candidate binding mode from original pose
    # Return: - True: Success
    #         - False: New pose is out of grid
    def reset_pose(self, translation, rotation, torsions):
        self.ligand.reset_atoms()
        self.protein.reset_flex_atoms()
        return self.set_pose(translation, rotation, torsions)

    # Return true if either or both ligand or/and lexible parts of protein is
    # out of predefined grid space. Else, return false.
    def check_out_of_grid(self):
        lo_x, lo_y, lo_z = self.grid.field.lo.xyz
        hi_x, hi_y, hi_z = self.grid.field.hi.xyz
        for atom in self.ligand.atoms:
            if (atom.tcoord.x <= lo_x or \
                atom.tcoord.y <= lo_y or \
                atom.tcoord.z <= lo_z or \
                atom.tcoord.x >= hi_x or \
                atom.tcoord.y >= hi_y or \
                atom.tcoord.z >= hi_z): return True
        for atom in self.protein.flex_atoms:
            if (atom.tcoord.x <= lo_x or \
                atom.tcoord.y <= lo_y or \
                atom.tcoord.z <= lo_z or \
                atom.tcoord.x >= hi_x or \
                atom.tcoord.y >= hi_y or \
                atom.tcoord.z >= hi_z): return True
        return False

    def get_non_bond_list(self):
        minmax_distance = self.bond.calc_minmax_distance()
        ligand_bond_matrix = \
            self.bond.construct_bond_matrix(self.ligand.ori_atoms, minmax_distance)
        protein_bond_matrix = \
            self.bond.construct_bond_matrix(self.protein.ori_flex_atoms, \
                                            minmax_distance)

        # Before combining ligand and protein bond matrices, shift protein ids \
        # by protein start index (p_idx)
        p_idx = len(self.ligand.ori_atoms)
        for ids in protein_bond_matrix:
            for i, id in enumerate(ids):
                ids[i] = id + p_idx
        bond_matrix = ligand_bond_matrix
        bond_matrix += protein_bond_matrix

        non_bond_matrix = self.bond.construct_non_bond_matrix(len(bond_matrix))
        non_bond_matrix = self.bond.weed_covalent_bond(bond_matrix, \
                                                       non_bond_matrix)
        non_bond_matrix = self.bond.weed_rigid_bond(non_bond_matrix, \
                                                    self.ligand, self.protein)

        self.non_bond_ligand, self.non_bond_ligand_receptor, \
            self.non_bond_receptor = \
                self.bond.convert_non_bond_matrix_to_list(non_bond_matrix, \
                                                          self.ligand, \
                                                          self.protein)

    def print_non_bond_matrix(self, non_bond_matrix):
        print "non_bond_matrix:"
        for i in xrange(len(non_bond_matrix)):
            res = "%2s  " % (i + 1)
            for j in xrange(len(non_bond_matrix)):
                if non_bond_matrix[i][j]:
                    res += "|X"
                else:
                    res += "|_"
            print "%s" % res

    def print_non_bond_list(self, non_bond_list, title = ""):
        print title
        print " Atom1-Atom2    Scaled(q1xq2) "
        print "------------------------------"
        for nbi in non_bond_list:
            print " %5d-%-5d     %6.2f" % (nbi.atom1, nbi.atom2, nbi.q1q2)

    # 3D Linear Interpolation
    @staticmethod
    def calc_linInterp3(grid, ligand, protein):
        lo_x, lo_y, lo_z = grid.field.lo.xyz
        spacing = grid.field.spacing
        atom_len = len(ligand.atoms) + len(protein.flex_atoms)

        u = []
        v = []
        w = []
        for atom in ligand.atoms:
            u.append((atom.tcoord.x - lo_x) / spacing)
            v.append((atom.tcoord.y - lo_y) / spacing)
            w.append((atom.tcoord.z - lo_z) / spacing)
        for atom in protein.flex_atoms:
            u.append((atom.tcoord.x - lo_x) / spacing)
            v.append((atom.tcoord.y - lo_y) / spacing)
            w.append((atom.tcoord.z - lo_z) / spacing)

        u0 = [int(i) for i in u]
        v0 = [int(i) for i in v]
        w0 = [int(i) for i in w]
                
        u1 = [i + 1 for i in u0]
        v1 = [i + 1 for i in v0]
        w1 = [i + 1 for i in w0]
    
        p0u = [u[i] - float(u0[i]) for i in xrange(atom_len)]
        p0v = [v[i] - float(v0[i]) for i in xrange(atom_len)]
        p0w = [w[i] - float(w0[i]) for i in xrange(atom_len)]

        p1u = [float(u1[i]) - u[i] for i in xrange(atom_len)]
        p1v = [float(v1[i]) - v[i] for i in xrange(atom_len)]
        p1w = [float(w1[i]) - w[i] for i in xrange(atom_len)]

        p000 = [p0u[i] * p0v[i] * p0w[i] for i in xrange(atom_len)]
        p001 = [p0u[i] * p0v[i] * p1w[i] for i in xrange(atom_len)]
        p010 = [p0u[i] * p1v[i] * p0w[i] for i in xrange(atom_len)]
        p011 = [p0u[i] * p1v[i] * p1w[i] for i in xrange(atom_len)]
        p100 = [p1u[i] * p0v[i] * p0w[i] for i in xrange(atom_len)]
        p101 = [p1u[i] * p0v[i] * p1w[i] for i in xrange(atom_len)]
        p110 = [p1u[i] * p1v[i] * p0w[i] for i in xrange(atom_len)]
        p111 = [p1u[i] * p1v[i] * p1w[i] for i in xrange(atom_len)]

        return u0, v0, w0, u1, v1, w1, \
               p000, p001, p010, p011, p100, p101, p110, p111

    # Calculate free energy
    def calc_intermolecular_energy(self):
        u0, v0, w0, u1, v1, w1, \
            p000, p001, p010, p011, p100, p101, p110, p111 = \
            self.calc_linInterp3(self.grid, self.ligand, self.protein)

        atom_len = len(self.ligand.atoms) + len(self.protein.flex_atoms)
        protein_idx = len(self.ligand.atoms)

        es = [] # Electrostatic
        for i in xrange(atom_len):
            e = 0.0
            e += p000[i] * self.grid.maps['e'][w1[i]][v1[i]][u1[i]]
            e += p001[i] * self.grid.maps['e'][w1[i]][v1[i]][u0[i]]
            e += p010[i] * self.grid.maps['e'][w1[i]][v0[i]][u1[i]]
            e += p011[i] * self.grid.maps['e'][w1[i]][v0[i]][u0[i]]
            e += p100[i] * self.grid.maps['e'][w0[i]][v1[i]][u1[i]]
            e += p101[i] * self.grid.maps['e'][w0[i]][v1[i]][u0[i]]
            e += p110[i] * self.grid.maps['e'][w0[i]][v0[i]][u1[i]]
            e += p111[i] * self.grid.maps['e'][w0[i]][v0[i]][u0[i]]
            es.append(e)

        ds = [] # Desolvation
        for i in xrange(atom_len):
            d = 0.0
            d += p000[i] * self.grid.maps['d'][w1[i]][v1[i]][u1[i]]
            d += p001[i] * self.grid.maps['d'][w1[i]][v1[i]][u0[i]]
            d += p010[i] * self.grid.maps['d'][w1[i]][v0[i]][u1[i]]
            d += p011[i] * self.grid.maps['d'][w1[i]][v0[i]][u0[i]]
            d += p100[i] * self.grid.maps['d'][w0[i]][v1[i]][u1[i]]
            d += p101[i] * self.grid.maps['d'][w0[i]][v1[i]][u0[i]]
            d += p110[i] * self.grid.maps['d'][w0[i]][v0[i]][u1[i]]
            d += p111[i] * self.grid.maps['d'][w0[i]][v0[i]][u0[i]]
            ds.append(d)

        ms = [] # Atom Type
        for i, atom in enumerate(self.ligand.atoms):
            m = 0.0
            type = atom.type
            m += p000[i] * self.grid.maps[type][w1[i]][v1[i]][u1[i]]
            m += p001[i] * self.grid.maps[type][w1[i]][v1[i]][u0[i]]
            m += p010[i] * self.grid.maps[type][w1[i]][v0[i]][u1[i]]
            m += p011[i] * self.grid.maps[type][w1[i]][v0[i]][u0[i]]
            m += p100[i] * self.grid.maps[type][w0[i]][v1[i]][u1[i]]
            m += p101[i] * self.grid.maps[type][w0[i]][v1[i]][u0[i]]
            m += p110[i] * self.grid.maps[type][w0[i]][v0[i]][u1[i]]
            m += p111[i] * self.grid.maps[type][w0[i]][v0[i]][u0[i]]
            ms.append(m)
        for idx, atom in enumerate(self.protein.flex_atoms):
            m = 0.0
            type = atom.type
            i = protein_idx + idx
            m += p000[i] * self.grid.maps[type][w1[i]][v1[i]][u1[i]]
            m += p001[i] * self.grid.maps[type][w1[i]][v1[i]][u0[i]]
            m += p010[i] * self.grid.maps[type][w1[i]][v0[i]][u1[i]]
            m += p011[i] * self.grid.maps[type][w1[i]][v0[i]][u0[i]]
            m += p100[i] * self.grid.maps[type][w0[i]][v1[i]][u1[i]]
            m += p101[i] * self.grid.maps[type][w0[i]][v1[i]][u0[i]]
            m += p110[i] * self.grid.maps[type][w0[i]][v0[i]][u1[i]]
            m += p111[i] * self.grid.maps[type][w0[i]][v0[i]][u0[i]]
            ms.append(m)

        # Electrostatic
        self.elecs = []
        self.elec_total = 0.0
        for i, atom in enumerate(self.ligand.atoms):
            self.elec = es[i] * atom.charge
            self.elecs.append(self.elec)
            self.elec_total += self.elec
        for idx, atom in enumerate(self.protein.flex_atoms):
            if atom.id in self.protein.ignore_inter:
                self.elec = 0.0
            else:
                i = protein_idx + idx
                self.elec = es[i] * atom.charge
            self.elecs.append(self.elec)
            self.elec_total += self.elec

        # Van der Waals
        self.emaps = []
        self.emap_total = 0.0
        for i, atom in enumerate(self.ligand.atoms):
            self.emap = ms[i] + ds[i] * abs(atom.charge)
            self.emaps.append(self.emap)
            self.emap_total += self.emap
        for idx, atom in enumerate(self.protein.flex_atoms):
            if atom.id in self.protein.ignore_inter:
                self.emap = 0.0
            else:
                i = protein_idx + idx
                self.emap = ms[i] + ds[i] * abs(atom.charge)
            self.emaps.append(self.emap)
            self.emap_total += self.emap

        return self.elec_total + self.emap_total

    def calc_intramolecular_energy(self):
        ns_intl_1 = self.bond.EnergyTable.NS_INTL - 1
        ns_el_1 = self.bond.EnergyTable.NS_EL - 1

        total_e_internal = 0.0
        # Intramolecular in the ligand
        for nb in self.non_bond_ligand:
            atom1 = nb.atom1 - 1
            atom_type1 = nb.atom_type1
            atom2 = nb.atom2 - 1
            atom_type2 = nb.atom_type2
            non_bond_type = nb.non_bond_type
            desolv = nb.desolv
            q1q2 = nb.q1q2

            # Get atoms distance
            atom_tcoord1 = self.ligand.atoms[atom1].tcoord
            atom_tcoord2 = self.ligand.atoms[atom2].tcoord
            r_tcoord2 = atom_tcoord1 - atom_tcoord2
            r2 = r_tcoord2.sq_hypotenuse()
            r2 = max(self.bond.RMIN_ELEC2, r2)  # Clamp r2 at RMIN_ELEC2
            i = int(r2 * self.bond.EnergyTable.SQA_DIV)
            # Make sure the indexes are not greater than NS_INTL -1 and
            # NS_EL - 1 respectively
            i_ns_intl = min(i, ns_intl_1)
            i_ns_el = min(i, ns_el_1)

            e_internal = 0.0
            if self.dps.calc_inter_elec_e:
                # Calculate Electrostatic Energy
                e_elec = q1q2 * self.bond.bound_et.inv_r_epsilon[i_ns_el]
                e_internal += e_elec
            if r2 < self.bond.EnergyTable.NBC2:
                # Calculate Desolvation Energy
                e_desolv = desolv * self.bond.bound_et.solvation[i_ns_intl]
                # Calculate Van der Waals and Hydrogen Bond Energies
                if self.bond.include_1_4_interactions and non_bond_type == 4:
                    e_internal += self.SCALE_1_4_INTERACTIONS + \
                                  (self.bond.bound_et.vdw_hb[(atom_type1, atom_type2)][i_ns_intl] + e_desolv)
                else:
                    e_internal += self.bond.bound_et.vdw_hb[(atom_type1, atom_type2)][i_ns_intl] + e_desolv

            total_e_internal += e_internal

        # Intermolecular ligand-receptor
        for nb in self.non_bond_ligand_receptor:
            atom1 = nb.atom1 - 1
            atom_type1 = nb.atom_type1
            atom2 = nb.atom2 - 1
            atom_type2 = nb.atom_type2
            non_bond_type = nb.non_bond_type
            desolv = nb.desolv
            q1q2 = nb.q1q2

            # Get atoms distance
            atom_tcoord1 = self.ligand.atoms[atom1].tcoord
            atom_tcoord2 = self.protein.flex_atoms[atom2].tcoord
            r_tcoord2 = atom_tcoord1 - atom_tcoord2
            r2 = r_tcoord2.sq_hypotenuse()
            r2 = max(self.bond.RMIN_ELEC2, r2)  # Clamp r2 at RMIN_ELEC2
            i = int(r2 * self.bond.EnergyTable.SQA_DIV)
            # Make sure the indexes are not greater than NS_INTL -1 and
            # NS_EL - 1 respectively
            i_ns_intl = min(i, ns_intl_1)
            i_ns_el = min(i, ns_el_1)

            e_internal = 0.0
            if self.dps.calc_inter_elec_e:
                # Calculate Electrostatic Energy
                e_elec = q1q2 * self.bond.bound_et.inv_r_epsilon[i_ns_el]
                e_internal += e_elec
            if r2 < self.bond.EnergyTable.NBC2:
                # Calculate Desolvation Energy
                e_desolv = desolv * self.bond.bound_et.solvation[i_ns_intl]
                # Calculate Van der Waals and Hydrogen Bond Energies
                if self.bond.include_1_4_interactions and non_bond_type == 4:
                    e_internal += self.SCALE_1_4_INTERACTIONS + \
                        (self.bond.bound_et.vdw_hb[(atom_type1, atom_type2)][i_ns_intl] + e_desolv)
                else:
                    e_internal += self.bond.bound_et.vdw_hb[(atom_type1, atom_type2)][i_ns_intl] + e_desolv

            total_e_internal += e_internal

        # Intramolecular in the receptor
        for nb in self.non_bond_receptor:
            atom1 = nb.atom1 - 1
            atom2 = nb.atom2 - 1
            atom_type1 = nb.atom_type1
            atom_type2 = nb.atom_type2
            non_bond_type = nb.non_bond_type
            desolv = nb.desolv
            q1q2 = nb.q1q2

            # Get atoms distance
            atom_tcoord1 = self.protein.flex_atoms[atom1].tcoord
            atom_tcoord2 = self.protein.flex_atoms[atom2].tcoord
            r_tcoord2 = atom_tcoord1 - atom_tcoord2
            r2 = r_tcoord2.sq_hypotenuse()
            r2 = max(self.bond.RMIN_ELEC2, r2)  # Clamp r2 at RMIN_ELEC2
            i = int(r2 * self.bond.EnergyTable.SQA_DIV)
            # Make sure the indexes are not greater than NS_INTL -1 and
            # NS_EL - 1 respectively
            i_ns_intl = min(i, ns_intl_1)
            i_ns_el = min(i, ns_el_1)

            e_internal = 0.0
            if self.dps.calc_inter_elec_e:
                # Calculate Electrostatic Energy
                e_elec = q1q2 * self.bond.bound_et.inv_r_epsilon[i_ns_el]
                e_internal += e_elec
            if r2 < self.bond.EnergyTable.NBC2:
                # Calculate Desolvation Energy
                e_desolv = desolv * self.bond.bound_et.solvation[i_ns_intl]
                # Calculate Van der Waals and Hydrogen Bond Energies
                if self.bond.include_1_4_interactions and non_bond_type == 4:
                    e_internal += self.SCALE_1_4_INTERACTIONS + \
                        (self.bond.bound_et.vdw_hb[(atom_type1, atom_type2)][i_ns_intl] + e_desolv)
                else:
                    e_internal += self.bond.bound_et.vdw_hb[(atom_type1, atom_type2)][i_ns_intl] + e_desolv

            total_e_internal += e_internal

        return total_e_internal

    def calc_energy(self):
        intermolecular_energy = self.calc_intermolecular_energy()
        intramolecular_energy = self.calc_intramolecular_energy()
        return intermolecular_energy + intramolecular_energy

    def test_print(self):
        for i, atom in enumerate(self.ligand.atoms):
            print "%2s: %2s - %8.3f, %8.3f, %8.3f | %+9.2f | %+9.2f" % \
                (i + 1, atom.type, \
                 atom.tcoord.x, atom.tcoord.y, atom.tcoord.z, \
                 self.emaps[i], self.elecs[i])
        protein_idx = len(self.ligand.atoms)
        for idx, atom in enumerate(self.protein.flex_atoms):
            i = protein_idx + idx
            print "%2s: %2s - %8.3f, %8.3f, %8.3f | %+9.2f | %+9.2f" % \
                (idx + 1, atom.type, \
                 atom.tcoord.x, atom.tcoord.y, atom.tcoord.z, \
                 self.emaps[i], self.elecs[i])

class DockOpenCL(Dock):
    def __init__(self):
        Dock.__init__(self)
        # Keep information about branches rotation sequence that contains
        # atom IDs start from 1. Use 0 to denote a non-atom information.
        self.longest_branch = 0
        self.longest_branch_np = np.array([0], dtype = int)
        self.longest_branch_buf = None
        self.branches_rot_anchor_np = np.array([], dtype = int)
        self.branches_rot_anchor_buf = None
        self.branches_rot_link_np = np.array([], dtype = int)
        self.branches_rot_link_buf = None
        self.branches_rot_size_np = np.array([], dtype = int)
        self.branches_rot_size_buf = None
        self.branches_rot_seq_np = np.array([], dtype = int)
        self.branches_rot_seq_buf = None

        # OpenCL
        self.cl_ctx = None
        self.cl_queue = None
        self.cl_filename = "./OpenCL/Dock.cl"
        self.cl_prg = None

        # OpenCL device buffer
        self.num_points1_np = np.array([], dtype = int)
        self.num_points1_buf = None
        self.lo_grid_np = np.array([], dtype = float)
        self.lo_grid_buf = None
        self.hi_grid_np = np.array([], dtype = float)
        self.hi_grid_buf = None
        self.dist_grid_np = np.array([], dtype = float)
        self.dist_grid_buf = None
        self.field_spacing_np = np.array([], dtype = float)
        self.field_spacing_buf = None

        self.electrostatic_lut_np = np.array([], dtype = int)
        self.desolvation_lut_np = np.array([], dtype = int)
        self.atom_type_map_lut_np = np.array([], dtype = int)
        self.maps_np = np.array([], dtype = float)
        self.electrostatic_lut_buf = None
        self.desolvation_lut_buf = None
        self.atom_type_map_lut_buf = None
        self.maps_buf = None
        self.ttl_maps_np = np.array([], dtype = int)
        self.ttl_maps_buf = None

        self.ttl_atom_types_np = np.array([], dtype = int)
        self.ttl_atom_types_buf = None
        # Collection of atom properties for all atoms in ligand and protein.
        # It's 2D array i by j for atom ID and atom properties (atom type,
        # charge).
        # Indexed for atom ID starts from 1 (index 0 is not in use).
        self.ttl_atom_properties_np = np.array([], dtype = int)
        self.ttl_atom_properties_buf = None
        self.atoms_properties_np = np.array([], dtype = float)
        self.atoms_properties_buf = None
        self.ttl_torsions_np = np.array([], dtype = int)
        self.ttl_torsions_buf = None
        self.ttl_ligand_atoms_np = np.array([], dtype = int)
        self.ttl_ligand_atoms_buf = None
        self.ori_atom_tcoords_np = np.array([], dtype = float)
        self.ori_atom_tcoords_buf = None
        self.ttl_atoms_np = np.array([], dtype = int)
        self.ttl_atoms_buf = None

        # Non-bond properties
        self.ttl_non_bond_properties_np = np.array([], dtype = int)
        self.ttl_non_bond_properties_buf = None
        self.ttl_non_bond_list_np = np.array([], dtype = int)
        self.ttl_non_bond_list_buf = None
        # Matrix of non-bond list by non-bond properties
        self.non_bond_list_np = np.array([], dtype = float)
        self.non_bond_list_buf = None

        # Bond properties
        self.bond_properties_np = np.array([], dtype = float)
        self.bond_properties_buf = None
        self.calc_inter_elec_e_np = np.array([], dtype = int)
        self.calc_inter_elec_e_buf = None
        self.include_1_4_interactions_np = np.array([], dtype = int)
        self.include_1_4_interactions_buf = None
        # Energy tables
        self.et_inv_r_epsilon_np = np.array([], dtype = float)
        self.et_inv_r_epsilon_buf = None
        self.et_solvation_np = np.array([], dtype = float)
        self.et_solvation_buf = None
        self.et_vdw_hb_np = np.array([], dtype = float)
        self.et_vdw_hb_buf = None

        # Collection of atom coordinates for entire population. It's a 3D array
        # i by j by k for atom ID, individual and atom coordinate respectively.
        # Indexed for atom ID starts from 1 (index 0 is not in use).
        self.ttl_poses_np = np.array([], dtype = int)
        self.ttl_poses_buf = None
        self.ori_poses_np = None
        self.ori_poses_buf = None
        self.poses_np = None
        self.poses_buf = None
        self.protein_ignore_inter_np = np.array([], dtype = int)
        self.protein_ignore_inter_buf = None
        self.ttl_protein_ignore_inter_np = np.array([], dtype = int)
        self.ttl_protein_ignore_inter_buf = None

        # Intermolecular energy
        self.elecs_np = None
        self.elecs_buf = None
        self.emaps_np = None
        self.emaps_buf = None
        self.elec_totals_np = None
        self.elec_totals_buf = None
        self.emap_totals_np = None
        self.emap_totals_buf = None
        # Intramolecular energy
        self.e_internals_np = None
        self.e_internals_buf = None
        self.e_internal_totals_np = None
        self.e_internal_totals_buf = None
        # Total energy
        self.e_totals_np = None
        self.e_totals_buf = None

    def setup_opencl(self, cl_ctx = None, cl_queue = None):
        self.cl_ctx = cl_ctx
        self.cl_queue = cl_queue
        fh = open(self.cl_filename, 'r')
        cl_code = "".join(fh.readlines())
        self.cl_prg = cl.Program(cl_ctx, cl_code).build()

    #TODO: Use self class cl_ctx and cl_queue
    def setup_opencl_buffer(self, ttl_poses = 0, \
                            cl_ctx = None, cl_queue = None):
        mf = cl.mem_flags

        #TODO: Move following variables into respective class
        # Field information (OpenCL device buffer)
        self.num_points1_np = np.array(self.grid.field.num_points1.xyz, \
                                       dtype = int)
        self.num_points1_buf = cl.Buffer(cl_ctx, \
                                         mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                         hostbuf = self.num_points1_np)
        self.lo_grid_np = np.array(self.grid.field.lo.xyz, dtype = float)
        self.lo_grid_buf = cl.Buffer(cl_ctx, \
                                     mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                     hostbuf = self.lo_grid_np)
        self.hi_grid_np = np.array(self.grid.field.hi.xyz, dtype = float)
        self.hi_grid_buf = cl.Buffer(cl_ctx, \
                                     mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                     hostbuf = self.hi_grid_np)
        self.dist_grid_np = self.hi_grid_np - self.lo_grid_np
        self.dist_grid_buf = cl.Buffer(cl_ctx, \
                                       mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                       hostbuf = self.dist_grid_np)
        self.field_spacing_np = np.array(self.grid.field.spacing, dtype = float)
        self.field_spacing_buf = cl.Buffer(cl_ctx, \
                                           mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                           hostbuf = self.field_spacing_np)
        # Maps information (OpenCL device buffer)
        map_idx = 0
        self.electrostatic_lut_np = np.array([map_idx], dtype = int)
        self.maps_np = self.grid.maps['e'].ravel()
        map_idx += 1
        self.desolvation_lut_np = np.array([map_idx], dtype = int)
        self.maps_np = np.vstack([self.maps_np, self.grid.maps['d'].ravel()])
        self.atom_type_map_lut_np = np.array([-1 for i in xrange(len(self.ligand.atom_types))],
                                             dtype = int)
        for atom_type in self.ligand.atom_types:
            if atom_type in self.grid.maps:
                map_idx += 1
                self.atom_type_map_lut_np[self.ligand.atom_types.index(atom_type)] = map_idx
                self.maps_np = np.vstack([self.maps_np, \
                                                    self.grid.maps[atom_type].ravel()])
        self.maps_np = np.transpose(self.maps_np).ravel()
        self.ttl_maps_np = np.array([map_idx + 1], dtype = int)

        self.electrostatic_lut_buf = cl.Buffer(cl_ctx, \
                                               mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                               hostbuf = self.electrostatic_lut_np)
        self.desolvation_lut_buf = cl.Buffer(cl_ctx, \
                                             mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                             hostbuf = self.desolvation_lut_np)
        self.atom_type_map_lut_buf = cl.Buffer(cl_ctx, \
                                               mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                               hostbuf = self.atom_type_map_lut_np)
        self.maps_buf = cl.Buffer(cl_ctx, \
                                  mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                  hostbuf = self.maps_np)
        self.ttl_maps_buf = cl.Buffer(cl_ctx, \
                                      mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                      hostbuf = self.ttl_maps_np)

        # Atoms properties (OpenCL device buffer)
        self.ttl_atom_types_np = np.array([len(self.ligand.atom_types)], dtype = int)
        self.ttl_atom_types_buf = cl.Buffer(cl_ctx, \
                                            mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                            hostbuf = self.ttl_atom_types_np)
        ttl_atom_properties = 2 # Atom type, charge
        self.ttl_atom_properties_np = np.array(ttl_atom_properties, dtype = int)
        self.ttl_atom_properties_buf = cl.Buffer(cl_ctx, \
                                                 mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                                 hostbuf = self.ttl_atom_properties_np)
        atoms_properties = [0.0 for i in xrange(ttl_atom_properties)]
        for atom in self.ligand.atoms:
            atoms_properties.append(float(self.ligand.atom_types.index(atom.type)))
            atoms_properties.append(atom.charge)
        for atom in self.protein.flex_atoms:
            atoms_properties.append(float(self.ligand.atom_types.index(atom.type)))
            atoms_properties.append(atom.charge)
        self.atoms_properties_np = np.array(atoms_properties, dtype = float)
        self.atoms_properties_buf = cl.Buffer(cl_ctx, \
                                              mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                              hostbuf = self.atoms_properties_np)
        # Molecule information (OpenCL device buffer)
        protein_idx = len(self.ligand.atoms)
        self.ttl_torsions_np = np.array([self.get_total_torsions()], \
                                        dtype = int)
        self.ttl_torsions_buf = cl.Buffer(cl_ctx, \
                                          mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                          hostbuf = self.ttl_torsions_np)
        self.ttl_ligand_atoms_np = np.array([len(self.ligand.atoms)], \
                                            dtype = int)
        self.ttl_ligand_atoms_buf = cl.Buffer(cl_ctx, \
                                              mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                              hostbuf = self.ttl_ligand_atoms_np)
        self.ligand.reset_atoms()
        self.protein.reset_flex_atoms()
        self.ori_atom_tcoords_np = np.vstack([np.array([0., 0., 0.], dtype = float), \
                                              self.ligand.get_atom_tcoords_in_numpy(), \
                                              self.protein.get_flex_atom_tcoords_in_numpy()])
        self.ori_atom_tcoords_buf = cl.Buffer(cl_ctx, \
                                              mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                              hostbuf = self.ori_atom_tcoords_np)
        # Non-bond properties
        ttl_non_bond_properties = 7
        ttl_non_bond_list = 0
        non_bond_list = []
        for nb in self.non_bond_ligand:
            ttl_non_bond_list += 1
            non_bond_list.append(float(nb.atom1))
            non_bond_list.append(self.ligand.atom_types.index(nb.atom_type1))
            non_bond_list.append(float(nb.atom2))
            non_bond_list.append(self.ligand.atom_types.index(nb.atom_type2))
            non_bond_list.append(float(nb.non_bond_type))
            non_bond_list.append(nb.desolv)
            non_bond_list.append(nb.q1q2)
        for nb in self.non_bond_ligand_receptor:
            ttl_non_bond_list += 1
            non_bond_list.append(float(nb.atom1))
            non_bond_list.append(self.ligand.atom_types.index(nb.atom_type1))
            non_bond_list.append(float(protein_idx) + float(nb.atom2))
            non_bond_list.append(self.ligand.atom_types.index(nb.atom_type2))
            non_bond_list.append(float(nb.non_bond_type))
            non_bond_list.append(nb.desolv)
            non_bond_list.append(nb.q1q2)
        for nb in self.non_bond_receptor:
            ttl_non_bond_list += 1
            non_bond_list.append(float(protein_idx) + float(nb.atom1))
            non_bond_list.append(self.ligand.atom_types.index(nb.atom_type1))
            non_bond_list.append(float(protein_idx) + float(nb.atom2))
            non_bond_list.append(self.ligand.atom_types.index(nb.atom_type2))
            non_bond_list.append(float(nb.non_bond_type))
            non_bond_list.append(nb.desolv)
            non_bond_list.append(nb.q1q2)
        self.ttl_non_bond_properties_np = np.array([ttl_non_bond_properties], \
                                                   dtype = int)
        self.ttl_non_bond_properties_buf = cl.Buffer(cl_ctx, \
                                                     mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                                     hostbuf = self.ttl_non_bond_properties_np)
        self.ttl_non_bond_list_np = np.array([ttl_non_bond_list], dtype = int)
        self.ttl_non_bond_list_buf = cl.Buffer(cl_ctx, \
                                               mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                               hostbuf = self.ttl_non_bond_list_np)
        self.non_bond_list_np = np.array(non_bond_list, dtype = float)
        self.non_bond_list_buf = cl.Buffer(cl_ctx, \
                                           mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                           hostbuf = self.non_bond_list_np)
        # Bond properties
        bond_properties = []
        bond_properties.append(float(self.bond.EnergyTable.NS_INTL - 1))
        bond_properties.append(float(self.bond.EnergyTable.NS_EL - 1))
        bond_properties.append(self.bond.RMIN_ELEC2)
        bond_properties.append(self.bond.EnergyTable.SQA_DIV)
        bond_properties.append(self.bond.EnergyTable.NBC2)
        bond_properties.append(self.SCALE_1_4_INTERACTIONS)
        self.bond_properties_np = np.array(bond_properties, dtype = float)
        self.bond_properties_buf = cl.Buffer(cl_ctx, \
                                             mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                             hostbuf = self.bond_properties_np)

        if self.dps.calc_inter_elec_e:
            self.calc_inter_elec_e_np = np.array([1], dtype = int)
        else:
            self.calc_inter_elec_e_np = np.array([0], dtype = int)
        self.calc_inter_elec_e_buf = cl.Buffer(cl_ctx, \
                                               mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                               hostbuf = self.calc_inter_elec_e_np)
        if self.bond.include_1_4_interactions:
            self.include_1_4_interactions_np = np.array([1], dtype = int)
        else:
            self.include_1_4_interactions_np = np.array([0], dtype = int)
        self.include_1_4_interactions_buf = cl.Buffer(cl_ctx, \
                                                      mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                                      hostbuf = self.include_1_4_interactions_np)
        # Energy tables
        self.et_inv_r_epsilon_np = np.array(self.bond.bound_et.inv_r_epsilon, \
                                            dtype = float)
        self.et_solvation_np = np.array(self.bond.bound_et.solvation, \
                                        dtype = float)
        et_vdw_hb = []
        for i, at_i in enumerate(self.ligand.atom_types):
            for at_j in self.ligand.atom_types[i:]:
                et_vdw_hb.append(self.bond.bound_et.vdw_hb[(at_i, at_j)])
        self.et_vdw_hb_np = np.array(et_vdw_hb, dtype = float)

        self.et_inv_r_epsilon_buf = cl.Buffer(cl_ctx, \
                                              mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                              hostbuf = self.et_inv_r_epsilon_np)
        self.et_solvation_buf = cl.Buffer(cl_ctx, \
                                          mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                          hostbuf = self.et_solvation_np)
        self.et_vdw_hb_buf = cl.Buffer(cl_ctx, \
                                       mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                       hostbuf = self.et_vdw_hb_np)

        # Poses (OpenCL device buffer)
        ttl_atoms = self.get_total_atoms()
        self.ttl_atoms_np = np.array([ttl_atoms], dtype = int)
        self.ttl_atoms_buf = cl.Buffer(cl_ctx, \
                                       mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                       hostbuf = self.ttl_atoms_np)
        self.ttl_poses_np = np.array([ttl_poses], dtype = int)
        self.ttl_poses_buf = cl.Buffer(cl_ctx, \
                                       mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                       hostbuf = self.ttl_poses_np)
        self.ori_poses_np = np.hstack([self.ori_atom_tcoords_np] * ttl_poses).ravel()
        self.ori_poses_buf = cl.Buffer(cl_ctx, \
                                       mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                       hostbuf = self.ori_poses_np)
        # Poses holds total atoms + 1 due to starting index of 1
        self.poses_buf = cl.array.zeros(cl_queue, \
                                        ((ttl_atoms + 1) * ttl_poses * 3), \
                                        dtype = float)
        protein_ignore_inter = []
        for id in self.protein.ignore_inter:
            protein_ignore_inter.append(protein_idx + id)
        self.protein_ignore_inter_np = np.array([protein_ignore_inter], \
                                                dtype = int)
        self.protein_ignore_inter_buf = cl.Buffer(cl_ctx, \
                                                  mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                                  hostbuf = self.protein_ignore_inter_np)
        self.ttl_protein_ignore_inter_np = np.array([len(protein_ignore_inter)], \
                                                    dtype = int)
        self.ttl_protein_ignore_inter_buf = cl.Buffer(cl_ctx, \
                                                      mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                                      hostbuf = self.ttl_protein_ignore_inter_np)
        # Intermolecular energy
        self.elecs_buf = cl.array.zeros(cl_queue, \
                                        ((ttl_atoms + 1) * ttl_poses), \
                                        dtype = float)
        self.emaps_buf = cl.array.zeros(cl_queue, \
                                        ((ttl_atoms + 1) * ttl_poses), \
                                        dtype = float)
        self.elec_totals_buf = cl.array.zeros(cl_queue, (ttl_poses), \
                                              dtype = float)
        self.emap_totals_buf = cl.array.zeros(cl_queue, (ttl_poses), \
                                              dtype = float)
        # Intramolecular energy
        self.e_internals_buf = cl.array.zeros(cl_queue, (ttl_poses * ttl_non_bond_list), \
                                              dtype = float)
        self.e_internal_totals_buf = cl.array.zeros(cl_queue, (ttl_poses), \
                                                    dtype = float)
        # Total energy
        self.e_totals_buf = cl.array.zeros(cl_queue, (ttl_poses), dtype = float)
        # Protein and ligand orientations and comformations (OpenCL device
        # buffer)
        self.set_branches_rotation_sequence(cl_ctx)

    def get_longest_branch(self):
        longest_branch = 0
        for branch in self.ligand.branches:
            branch_len = len(branch.all_atom_ids)
            if branch_len > longest_branch:
                longest_branch = branch_len
        for branch in self.protein.flex_branches:
            branch_len = len(branch.all_atom_ids)
            if branch_len > longest_branch:
                longest_branch = branch_len
        return longest_branch
        
    # Branches rotation sequence starts from the closest to leaf
    def set_branches_rotation_sequence(self, cl_ctx = None):
        for branch in self.ligand.branches:
            branch.molecule = 'l' # l for ligand
            self.sorted_branches.append(branch)
        for branch in self.protein.flex_branches:
            branch.molecule = 'p' # p for protein
            self.sorted_branches.append(branch)
        self.sorted_branches = sorted(self.sorted_branches, \
                                      key=lambda branch: len(branch.all_atom_ids))

        self.longest_branch = self.get_longest_branch()
        protein_idx = len(self.ligand.atoms)
        branches_rot_anchor = []
        branches_rot_link = []
        branches_rot_size = []
        branches_rot_seq = []
        for branch in self.sorted_branches:
            branch_rot_seq = [0 for i in xrange(self.longest_branch)]
            # Get the atoms from either ligand or protein based on branch
            # molecule information
            if branch.molecule == 'l':
                molecule_atoms = self.ligand.atoms
                shift_atom_id = 0
            else: # 'p'
                molecule_atoms = self.protein.flex_atoms
                shift_atom_id = protein_idx
            # Get atom coordinates
            seq_idx = 0
            for atom in molecule_atoms:
                if atom.id in [branch.anchor_id] + [branch.link_id] + \
                    branch.all_atom_ids:
                        # Anchor and link atoms are not for rotation
                        if atom.id == branch.anchor_id:
                            branches_rot_anchor.append(atom.id + shift_atom_id)
                            continue
                        if atom.id == branch.link_id:
                            branches_rot_link.append(atom.id + shift_atom_id)
                            continue
                        branch_rot_seq[seq_idx] = atom.id + shift_atom_id
                        seq_idx += 1
            branches_rot_size.append(seq_idx)
            branches_rot_seq.append(branch_rot_seq)

        mf = cl.mem_flags
        self.longest_branch_np = np.array([self.longest_branch], dtype = int)
        self.longest_branch_buf = cl.Buffer(cl_ctx, \
                                            mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                            hostbuf = self.longest_branch_np)
        self.branches_rot_anchor_np = np.array(branches_rot_anchor, dtype = int)
        self.branches_rot_anchor_buf = cl.Buffer(cl_ctx, \
                                                 mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                                 hostbuf = self.branches_rot_anchor_np)
        self.branches_rot_link_np = np.array(branches_rot_link, dtype = int)
        self.branches_rot_link_buf = cl.Buffer(cl_ctx, \
                                               mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                               hostbuf = self.branches_rot_link_np)
        self.branches_rot_size_np = np.array(branches_rot_size, dtype = int)
        self.branches_rot_size_buf = cl.Buffer(cl_ctx, \
                                               mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                               hostbuf = self.branches_rot_size_np)
        self.branches_rot_seq_np = np.array(branches_rot_seq, dtype = int).ravel()
        self.branches_rot_seq_buf = cl.Buffer(cl_ctx, \
                                              mf.READ_ONLY | mf.COPY_HOST_PTR, \
                                              hostbuf = self.branches_rot_seq_np)

        if DEBUG:
            print "Branches Rotation Sequence Table:"
            for i in xrange(self.get_total_torsions()):
                print "%2d %2d (%2d) %s" % (self.branches_rot_anchor_np[i], \
                                            self.branches_rot_link_np[i], \
                                            self.branches_rot_size_np[i], \
                                            self.branches_rot_seq_np[(i * self.longest_branch): \
                                                                     ((i + 1) * self.longest_branch)])

    #TODO: Use self class cl_queue
    def set_poses(self, ttl_poses = 0, individuals_buf = None, \
                  cl_queue = None):
        print self.ttl_torsions_np
        # Rotate rotatable branches/bonds for both ligand and protein.
        # Rotation is expected to be in radian.
        self.cl_prg.rotate_branches(cl_queue, \
                                    (ttl_poses * self.longest_branch,), None, \
                                    self.ttl_torsions_buf, \
                                    individuals_buf.data, \

                                    self.longest_branch_buf, \
                                    self.branches_rot_anchor_buf, \
                                    self.branches_rot_link_buf, \
                                    self.branches_rot_size_buf, \
                                    self.branches_rot_seq_buf, \

                                    self.ttl_poses_buf, \
                                    self.poses_buf.data)
        ttl_ligand_atoms = int(self.ttl_ligand_atoms_np[0])
        # Transform (translate and rotate) ligand root (whole body)
        self.cl_prg.transform_ligand_root(cl_queue, \
                                          (ttl_ligand_atoms * ttl_poses,), \
                                          None, \
                                          self.ttl_ligand_atoms_buf, \
                                          individuals_buf.data, \
                                          self.ttl_torsions_buf, \
                                          self.ttl_poses_buf, \
                                          self.poses_buf.data)

    #TODO: Use self class cl_ctx and cl_queue
    def reset_poses(self, ttl_poses = 0, individuals_buf = None, \
                    cl_ctx = None, cl_queue = None):
        cl.enqueue_copy(cl_queue, self.poses_buf.data, self.ori_poses_buf)

        if DEBUG:
            self.poses_np = self.poses_buf.get()
            print "Ori:"
            for i in xrange(self.ttl_atoms_np[0] + 1):
                print self.poses_np[(i * ttl_poses * 3):((i * ttl_poses * 3) + 3)]
        
        self.set_poses(ttl_poses, individuals_buf, cl_queue)

        if DEBUG:
            self.poses_np = self.poses_buf.get()
            print "New:"
            for i in xrange(self.ttl_atoms_np[0] + 1):
                print self.poses_np[((i * ttl_poses * 3) + 447):((i * ttl_poses * 3) + 450)]

    # Calculate free energy
    def calc_intermolecular_energy(self):
        ttl_atoms = int(self.ttl_atoms_np[0])
        ttl_poses = int(self.ttl_poses_np[0])
        self.cl_prg.calc_inter_energy(self.cl_queue, \
                                      (ttl_atoms * ttl_poses,), None, \
                                      self.ttl_poses_buf, \
                                      self.lo_grid_buf, \
                                      self.hi_grid_buf, \
                                      self.field_spacing_buf, \
                                      self.poses_buf.data, \

                                      self.num_points1_buf, \
                                      self.ttl_maps_buf, \
                                      self.electrostatic_lut_buf, \
                                      self.desolvation_lut_buf, \
                                      self.atom_type_map_lut_buf, \
                                      self.maps_buf, \

                                      self.ttl_atom_properties_buf, \
                                      self.atoms_properties_buf, \
                                      self.protein_ignore_inter_buf, \
                                      self.ttl_protein_ignore_inter_buf, \

                                      self.elecs_buf.data, \
                                      self.emaps_buf.data)
        self.cl_prg.calc_total_inter_energy(self.cl_queue, \
                                            (2 * ttl_poses,), None, \
                                            self.ttl_atoms_buf, \
                                            self.ttl_poses_buf, \
                                            self.elecs_buf.data, \
                                            self.emaps_buf.data, \

                                            self.elec_totals_buf.data, \
                                            self.emap_totals_buf.data)

        if DEBUG:
            self.elec_totals_np = self.elec_totals_buf.get()
            print self.elec_totals_np
            self.emap_totals_np = self.emap_totals_buf.get()
            print self.emap_totals_np
    
    def calc_intramolecular_energy(self):
        ttl_poses = int(self.ttl_poses_np[0])
        ttl_non_bond_list = int(self.ttl_non_bond_list_np[0])
        self.cl_prg.calc_intra_energy(self.cl_queue, \
                                      (ttl_poses * ttl_non_bond_list,), None, \
                                      self.ttl_poses_buf, \
                                      self.poses_buf.data, \
                                      self.lo_grid_buf, \
                                      self.hi_grid_buf, \

                                      self.ttl_non_bond_list_buf, \
                                      self.ttl_non_bond_properties_buf, \
                                      self.non_bond_list_buf, \

                                      self.ttl_atom_types_buf, \
                                      self.bond_properties_buf, \

                                      self.calc_inter_elec_e_buf, \
                                      self.include_1_4_interactions_buf, \

                                      self.et_inv_r_epsilon_buf, \
                                      self.et_solvation_buf, \
                                      self.et_vdw_hb_buf, \

                                      self.e_internals_buf.data)

        if DEBUG:
            self.e_internals_np = self.e_internals_buf.get()
            print self.e_internals_np

        self.cl_prg.calc_total_intra_energy(self.cl_queue, \
                                            (ttl_poses,), None, \
                                            self.ttl_poses_buf, \
                                            self.ttl_non_bond_list_buf, \
                                            self.e_internals_buf.data, \

                                            self.e_internal_totals_buf.data)

        if DEBUG:
            self.e_internal_totals_np = self.e_internal_totals_buf.get()
            print self.e_internal_totals_np

    def calc_total_energy(self):
        ttl_poses = int(self.ttl_poses_np[0])
        self.cl_prg.calc_total_energy(self.cl_queue, (ttl_poses,), None, \
                                      self.elec_totals_buf.data, \
                                      self.emap_totals_buf.data, \
                                      self.e_internal_totals_buf.data, \
                                      self.e_totals_buf.data)

    def calc_energy(self):
        self.calc_intermolecular_energy()
        self.calc_intramolecular_energy()
        self.calc_total_energy()

