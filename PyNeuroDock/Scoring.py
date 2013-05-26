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

from NeuroDock import NeuroDock
from Axis3 import Axis3
from Quaternion import Quaternion
from Constants import DEG2RAD

class Scoring:
    def __init__(self):
        docking_parameter_file = "./Parameters/ind_scoring.dpf"
        self.neuroDock = NeuroDock(docking_parameter_file)
        self.neuroDock.run()
        self.neuroDock.dock.get_non_bond_list()

    def __call__(self, pose):
        translation = Axis3(pose[0], pose[1], pose[2])
        rotation = Quaternion(pose[3], pose[4], pose[5], pose[6])
        torsion = [DEG2RAD * x for x in pose[7:]]
        return self.neuroDock.dock.energy(translation, rotation, torsion)


#bar - start
sc = Scoring()
pose = [2.056477, 5.846611, -7.245407, \
        0.532211, 0.379383, 0.612442, 0.444674, \
        -122.13, -179.41, \
        -141.59,  177.29, \
        -179.46,   -9.31, \
         132.37,  -89.19, \
          78.43,   22.22, \
          71.37,   59.52]

print "Energy:\n%s" % str(sc(pose))
print "Energy:\n%s" % str(sc(pose))
print "Energy:\n%s" % str(sc(pose))
#bar - stop
