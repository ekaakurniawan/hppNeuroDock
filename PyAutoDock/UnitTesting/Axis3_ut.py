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
from Axis3 import *

class ListTestInit(unittest.TestCase):
    def testInit(self):
        axis = Axis3()
        self.assertEquals(axis.xyz, [0.0, 0.0, 0.0])

        axis = Axis3(1.2, 2.3, 3.4)
        self.assertEquals(axis.xyz, [1.2, 2.3, 3.4])

class ListTestNormalize(unittest.TestCase):
    def testNormalize(self):
        axis = Axis3(1.2, 2.3, 3.4)
        axis.normalize()
        self.assertLessEqual(axis.x - 0.2806, 0.0001)
        self.assertGreaterEqual(axis.x - 0.2806, -0.0001)
        self.assertLessEqual(axis.y - 0.5378, 0.0001)
        self.assertGreaterEqual(axis.y - 0.5378, -0.0001)
        self.assertLessEqual(axis.z - 0.7950, 0.0001)
        self.assertGreaterEqual(axis.z - 0.7950, -0.0001)

def suite():
    suite1 = unittest.makeSuite(ListTestInit)
    suite2 = unittest.makeSuite(ListTestNormalize)
    return unittest.TestSuite((suite1, suite2))

if __name__ == '__main__':
    unittest.main()
