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
from Quaternion import *
from LFSR import *

class QuaternionTestInit(unittest.TestCase):
    def testInit(self):
        q = Quaternion()
        self.assertEquals(str(q), "Quaternion: 1.0000 + 0.0000i + 0.0000j + 0.0000k")
        q = Quaternion(1.2, 2.3, 3.4, 4.5)
        self.assertEquals(str(q), "Quaternion: 1.2000 + 2.3000i + 3.4000j + 4.5000k")

class QuaternionTestCopy(unittest.TestCase):
    def testCopy(self):
        q = Quaternion(1.2, 2.3, 3.4, 4.5)
        q1 = q.copy()
        q1.a = 100.002
        q1.b = 200.003
        q1.c = 300.004
        q1.d = 400.005
        self.assertEquals(str(q), "Quaternion: 1.2000 + 2.3000i + 3.4000j + 4.5000k")
        self.assertEquals(str(q1), "Quaternion: 100.0020 + 200.0030i + 300.0040j + 400.0050k")

class QuaternionTestUniform(unittest.TestCase):
    def testUniform(self):
        exp_q = ["Quaternion: -0.5105 + -0.7011i + 0.0180j + 0.4976k",
                 "Quaternion: 0.1267 + -0.5028i + -0.7474j + -0.4153k",
                 "Quaternion: 0.1504 + 0.0998i + -0.3025j + 0.9359k",
                 "Quaternion: 0.6418 + -0.5843i + 0.2151j + -0.4476k",
                 "Quaternion: 0.9726 + -0.0334i + 0.2176j + -0.0742k",
                 "Quaternion: 0.0096 + 0.0015i + -0.0778j + 0.9969k",
                 "Quaternion: -0.4398 + 0.5669i + -0.2332j + -0.6564k",
                 "Quaternion: 0.2429 + 0.2133i + -0.3556j + 0.8769k",
                 "Quaternion: -0.0521 + -0.3024i + -0.8972j + 0.3175k"]
                 
        exp_angle = [-2.0702,
                     2.8876,
                     2.8396,
                     1.7478,
                     0.4689,
                     3.1225,
                     -2.2309,
                     2.6510,
                     -3.0374]

        exp_axis = ["Axis3: -0.8153, 0.0209, 0.5787",
                    "Axis3: -0.5069, -0.7535, -0.4187",
                    "Axis3: 0.1009, -0.3059, 0.9467",
                    "Axis3: -0.7619, 0.2806, -0.5837",
                    "Axis3: -0.1437, 0.9367, -0.3193",
                    "Axis3: 0.0015, -0.0778, 0.9970",
                    "Axis3: 0.6312, -0.2597, -0.7309",
                    "Axis3: 0.2199, -0.3666, 0.9040",
                    "Axis3: -0.3028, -0.8985, 0.3179"]

        rng = LFSR(lfsr = 1070, bit_len = 16)
        q = Quaternion()
        for i in xrange(9):
            q.uniform(rng)
            angle, axis = q.get_angle_axis()
            self.assertEquals(str(q), exp_q[i])
            self.assertLessEqual(angle - exp_angle[i], 0.0001)
            self.assertGreaterEqual(angle - exp_angle[i], -0.0001)
            self.assertEquals(str(axis), exp_axis[i])

class QuaternionTestMultiplication(unittest.TestCase):
    def testMultiplication(self):
        q = Quaternion(0.3627, 0.3898, 0.8427, 0.0798)
        q1 = Quaternion(0.5221, -0.6938, 0.0376, -0.4946)
        q *= q1
        self.assertEquals(str(q), "Quaternion: 0.4676 + -0.4679i + 0.5910j + 0.4616k")

        # Magnitude checking after multiplication
        self.assertLessEqual(1 - q.magnitude(), 0.00001)
        self.assertGreaterEqual(1 - q.magnitude(), -0.00001)
        self.assertLessEqual(1 - q1.magnitude(), 0.00001)
        self.assertGreaterEqual(1 - q1.magnitude(), -0.00001)

class QuaternionTestNormalize(unittest.TestCase):
    def testNormalize(self):
        q = Quaternion(0.4676, -0.4679, 0.5910, 0.4616)
        q.normalize()
        self.assertLessEqual(1 - q.magnitude(), 0.00001)
        self.assertGreaterEqual(1 - q.magnitude(), -0.00001)

        q = Quaternion(-0.0918, 0.0463, -0.2413, -0.9650)
        q.normalize()
        self.assertLessEqual(1 - q.magnitude(), 0.00001)
        self.assertGreaterEqual(1 - q.magnitude(), -0.00001)

        q = Quaternion(a = 1000.0002, b = 2.03, c = 0.04, d = 40004.5)
        q.normalize()
        self.assertLessEqual(1 - q.magnitude(), 0.00001)
        self.assertGreaterEqual(1 - q.magnitude(), -0.00001)

class QuaternionTestConjugate(unittest.TestCase):
    def testConjugate(self):
        q = Quaternion(-0.3308, -0.7273, -0.3217, 0.5080)
        q.conjugate()
        self.assertEquals(str(q), "Quaternion: -0.3308 + 0.7273i + 0.3217j + -0.5080k")

class QuaternionTestIdentity(unittest.TestCase):
    def testInit(self):
        q = Quaternion(-0.3308, 0.7273, 0.3217, -0.5080)
        q.identity()
        self.assertEquals(str(q), "Quaternion: 1.0000 + 0.0000i + 0.0000j + 0.0000k")
        self.assertLessEqual(1 - q.magnitude(), 0.00001)
        self.assertGreaterEqual(1 - q.magnitude(), -0.00001)

class QuaternionTestConversion(unittest.TestCase):
    def testQuaternionToAngleAxis(self):
        q = Quaternion(a = 0.707, b = -0.240, c = -0.665, d = 0.000)
        angle, axis = q.get_angle_axis()
        self.assertLessEqual(angle - 1.5710983268, 0.0000000001)
        self.assertGreaterEqual(angle - 1.5710983268, -0.0000000001)
        self.assertEquals(str(axis), "Axis3: -0.3395, -0.9406, 0.0000")

    def testAngleAxisToQuaternion(self):
        q = Quaternion()
        q.identity()
        angle = 1.57
        axis = Axis3(x = -0.340, y = -0.940, z = 0.000)
        q.set_angle_axis(angle, axis)
        self.assertEquals(str(q), "Quaternion: 0.7074 + -0.2404i + -0.6647j + 0.0000k")

def suite():
    suite1 = unittest.makeSuite(QuaternionTestInit)
    suite2 = unittest.makeSuite(QuaternionTestCopy)
    suite3 = unittest.makeSuite(QuaternionTestUniform)
    suite4 = unittest.makeSuite(QuaternionTestMultiplication)
    suite5 = unittest.makeSuite(QuaternionTestNormalize)
    suite6 = unittest.makeSuite(QuaternionTestConjugate)
    suite7 = unittest.makeSuite(QuaternionTestIdentity)
    suite8 = unittest.makeSuite(QuaternionTestConversion)
    return unittest.TestSuite((suite1, suite2, suite3, suite4, suite5, suite6, \
                               suite7, suite8, ))

if __name__ == '__main__':
    unittest.main()
