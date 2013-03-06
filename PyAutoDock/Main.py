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

from Ligand import Ligand
from Grid import Grid
from Parameters import Field
from Map import ElectrostaticMap, DesolvationMap, AtomTypeMap
from Dock import Dock

ligand = Ligand()
ligand.read_pdbqt("./Inputs/ind.pdbqt")
ligand.update_tcoord_model("./Results/ind.model")

grid = Grid()
grid.field = Field("./Parameters/hsg1_rigid.maps.fld")
grid.maps['e'] = ElectrostaticMap("./Maps/hsg1_rigid.e.map", grid.field).map
grid.maps['d'] = DesolvationMap("./Maps/hsg1_rigid.d.map", grid.field).map
grid.maps['A'] = AtomTypeMap("./Maps/hsg1_rigid.A.map", grid.field).map
grid.maps['C'] = AtomTypeMap("./Maps/hsg1_rigid.C.map", grid.field).map
grid.maps['HD'] = AtomTypeMap("./Maps/hsg1_rigid.HD.map", grid.field).map
grid.maps['N'] = AtomTypeMap("./Maps/hsg1_rigid.N.map", grid.field).map
grid.maps['NA'] = AtomTypeMap("./Maps/hsg1_rigid.NA.map", grid.field).map
grid.maps['OA'] = AtomTypeMap("./Maps/hsg1_rigid.OA.map", grid.field).map

dock = Dock(grid, ligand)
dock.calc_energy()
dock.test_print() #bar
