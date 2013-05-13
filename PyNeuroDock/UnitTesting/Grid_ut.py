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

import unittest
from Grid import Grid, Field
from Axis3 import Axis3

class ParametersField(unittest.TestCase):
    def testRead(self):
        field = Field()
        field.read("./Parameters/hsg1_rigid.maps.fld")
        
        self.assertEquals(field.spacing, 0.375)
        self.assertEquals(str(field.num_points), str(Axis3(60, 60, 60)))
        self.assertEquals(str(field.num_points1), str(Axis3(61, 61, 61)))
        self.assertEquals(str(field.center), str(Axis3(2.500, 6.500, -7.500)))
        self.assertEquals(str(field.lo), str(Axis3(-8.750, -4.750, -18.750)))
        self.assertEquals(str(field.hi), str(Axis3(13.750, 17.750, 3.750)))

def suite():
    suite1 = unittest.makeSuite(ParametersField)
    return unittest.TestSuite((suite1, ))

if __name__ == '__main__':
    unittest.main()
