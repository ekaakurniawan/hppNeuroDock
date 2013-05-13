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
#  - AutoDock 4.2.3 Source Code (main.cc, parse_dpf_line.cc, dpftoken.h)
#    http://autodock.scripps.edu

from Grid import Field
from Dock import Dock, DockingParameters
from Map import ElectrostaticMap, DesolvationMap, AtomTypeMap
from Axis3 import Axis3
from Quaternion import Quaternion
from Constants import DEG2RAD
from Atom import Bond

class NeuroDock:
    def __init__(self, docking_parameter_file = None):
        self.docking_parameter_file = docking_parameter_file
        # Instantiate Dock
        self.dock = Dock()
        
        self.grid_field_file = ""
        self.atom_type_map_files = {}
        self.electrostatic_map_file = ""
        self.desolvation_map_file = ""
        self.ligand_file = ""
        self.protein_file = ""

    def run(self):
        with open(self.docking_parameter_file, 'r') as p_file:
            for line in p_file:
                if line.startswith("intelec"):
                    self.dock.dps.calc_inter_elec_e = True

                # Atomic bonding parameter file
                if line.startswith("b_prm"):
                    self.dock.bond.read(line.split()[1])
                    print self.dock.bond #bar

                if line.startswith("fld"):
                    self.grid_field_file = "./Parameters/" + line.split()[1]
                    self.dock.grid.field = Field(self.grid_field_file)

                if line.startswith("ligand_types"):
                    for type in line.split('#')[0].split()[1:]:
                        self.atom_type_map_files[type] = ""

                if line.startswith("map"):
                    filename = line.split()[1]
                    type = filename.split('.')[1]
                    self.atom_type_map_files[type] = "./Maps/" + filename
                    self.dock.grid.maps[type] = AtomTypeMap(self.atom_type_map_files[type], self.dock.grid.field).map

                if line.startswith("elecmap"):
                    self.electrostatic_map_file = "./Maps/" + line.split()[1]
                    self.dock.grid.maps['e'] = ElectrostaticMap(self.electrostatic_map_file, self.dock.grid.field).map

                if line.startswith("desolvmap"):
                    self.desolvation_map_file = "./Maps/" + line.split()[1]
                    self.dock.grid.maps['d'] = DesolvationMap(self.desolvation_map_file, self.dock.grid.field).map

                # Set movable molecule (ligand)
                if line.startswith("move"):
                    self.ligand_file = "./Inputs/" + line.split()[1]
                    self.dock.ligand.read_pdbqt(self.ligand_file)

                # Set flexible portion of molecule (normally protein receptor)
                if line.startswith("flexres"):
                    self.protein_file = "./Inputs/" + line.split()[1]
                    self.dock.protein.read_flex_pdbqt(self.protein_file)

                if line.startswith("about"):
                    about = Axis3(0.0, 0.0, 0.0)
                    about.xyz = [float(axis) for axis in line.split()[1:4]]
                    self.dock.ligand.about = about
                    # Zero-out central point. Applicable only for all ligand
                    # atoms (not protein atoms)
                    for i in xrange(len(self.dock.ligand.atoms)):
                        self.dock.ligand.atoms[i].tcoord -= about

                if line.startswith("include_1_4_interactions"):
                    # By default, 1-4 interactions is disabled. Only 1-1, 1-2,
                    # and 1-3 interactions are considered.
                    self.dock.bond.include_1_4_interactions = True

                if line.startswith("ga_run"):
                    print self.dock.ligand #bar
                    print self.dock.protein #bar

                    non_bond_ligand, non_bond_ligand_receptor, \
                        non_bond_receptor = self.dock.get_non_bond_list()

                    self.dock.print_non_bond_list(non_bond_ligand, \
                                                  "Non-bonded Pair Ligand-Ligand:") #bar
                    self.dock.print_non_bond_list(non_bond_ligand_receptor, \
                                                  "Non-bonded Pair Ligand-Receptor:") #bar
                    self.dock.print_non_bond_list(non_bond_receptor, \
                                                  "Non-bonded Pair Receptor-Receptor:") #bar

                    translation = Axis3(2.056477, 5.846611, -7.245407)
                    rotation = Quaternion(0.532211, 0.379383, 0.612442, 0.444674)
                    torsion = [DEG2RAD * x for x in [-122.13, -179.41, \
                                                     -141.59,  177.29, \
                                                     -179.46,   -9.31, \
                                                      132.37,  -89.19, \
                                                       78.43,   22.22, \
                                                       71.37,   59.52]]
                    if self.dock.set_pose(translation, rotation, torsion):
                        self.dock.calc_energy()
                        self.dock.test_print() #bar


#bar - start
docking_parameter_file = "./Parameters/ind.dpf"
neuroDock = NeuroDock(docking_parameter_file)
neuroDock.run()
#bar - stop
