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
from Map import ElectrostaticMap, DesolvationMap, AtomTypeMap
from Grid import Grid, Field

class MapElectrostatic(unittest.TestCase):
    def testRead(self):
        field = Field()
        field.read("./Parameters/hsg1_rigid.maps.fld")
        em = ElectrostaticMap()
        em.read("./Maps/hsg1_rigid.e.map", field)

        # 1
        exp = [-0.486,
               -0.436,
               -0.410,
               -0.408,
               -0.426,
               -0.459,
               -0.499,
               -0.555,
               -0.620,
               -0.680]
        for i in xrange(10):
            self.assertEquals(em.map[0][0][i], exp[i])

        # 316
        exp = [-0.368,
               -0.175,
                0.001,
                0.053,
               -0.019,
               -0.117,
               -0.187,
               -0.249,
               -0.312,
               -0.395]
        for i in xrange(10):
            self.assertEquals(em.map[0][5][10 + i], exp[i])

        # 183110
        exp = [-0.622,
               -0.730,
               -0.411,
                0.367,
                1.223,
                1.625,
                1.909,
                2.677,
                2.379,
               -0.580]
        for i in xrange(10):
            self.assertEquals(em.map[49][12][48 + i], exp[i])

        # 226972
        exp = [-0.005,
               -0.006,
               -0.007,
               -0.007,
               -0.007,
               -0.008,
               -0.008,
               -0.007,
               -0.007,
               -0.007]
        for i in xrange(10):
            self.assertEquals(em.map[60][60][51 + i], exp[i])

class MapDesolvation(unittest.TestCase):
    def testRead(self):
        field = Field()
        field.read("./Parameters/hsg1_rigid.maps.fld")
        em = ElectrostaticMap()
        em.read("./Maps/hsg1_rigid.d.map", field)

        # 1
        exp = [0.704,
               0.701,
               0.693,
               0.691,
               0.685,
               0.689,
               0.689,
               0.689,
               0.694,
               0.690]
        for i in xrange(10):
            self.assertEquals(em.map[0][0][i], exp[i])

        # 316
        exp = [0.766,
               0.764,
               0.751,
               0.753,
               0.741,
               0.731,
               0.730,
               0.728,
               0.722,
               0.721]
        for i in xrange(10):
            self.assertEquals(em.map[0][5][10 + i], exp[i])

        # 183110
        exp = [1.218,
               1.232,
               1.238,
               1.259,
               1.259,
               1.255,
               1.242,
               1.240,
               1.244,
               1.243]
        for i in xrange(10):
            self.assertEquals(em.map[49][12][48 + i], exp[i])

        # 226972
        exp = [0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0]
        for i in xrange(10):
            self.assertEquals(em.map[60][60][51 + i], exp[i])

class MapAtomType(unittest.TestCase):
    def testRead(self):
        field = Field()
        field.read("./Parameters/hsg1_rigid.maps.fld")
        em = ElectrostaticMap()
        em.read("./Maps/hsg1_rigid.C.map", field)

        # 1
        exp = [1.092,
               0.806,
               0.553,
               0.263,
               0.190,
               0.598,
               2.408,
               8.181,
               23.834,
               54.212]
        for i in xrange(10):
            self.assertEquals(em.map[0][0][i], exp[i])

        # 316
        exp = [2165.672,
               16349.299,
               70409.602,
               78148.484,
               21661.604,
               2812.866,
               387.952,
               75.797,
               21.270,
               10.186]
        for i in xrange(10):
            self.assertEquals(em.map[0][5][10 + i], exp[i])

        # 183110
        exp = [425.622,
               1197.815,
               4611.195,
               22100.568,
               57932.355,
               44180.094,
               108190.430,
               102640.172,
               116771.336,
               199047.891]
        for i in xrange(10):
            self.assertEquals(em.map[49][12][48 + i], exp[i])

        # 226972
        exp = [0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0,
               0.0]
        for i in xrange(10):
            self.assertEquals(em.map[60][60][51 + i], exp[i])

def suite():
    suite1 = unittest.makeSuite(MapElectrostatic)
    suite2 = unittest.makeSuite(MapDesolvation)
    suite3 = unittest.makeSuite(MapAtomType)
    return unittest.TestSuite((suite1, suite2, suite3))

if __name__ == '__main__':
    unittest.main()
