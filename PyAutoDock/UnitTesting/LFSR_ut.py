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
from LFSR import *

class ListTestGenerate(unittest.TestCase):
    def testGenerate(self):
        exp = [33303,
               49419,
               24709,
               12354,
               38945,
               52240,
               58888,
               29444,
               14722]
        lfsr = LFSR(lfsr = 1070, bit_len = 16)
        for i in xrange(9):
            self.assertEquals(lfsr.generate(), exp[i])

        exp = [0.612319946289,
               0.80615234375,
               0.903076171875,
               0.951538085938,
               0.975769042969,
               0.987884521484,
               0.993942260742,
               0.496963500977,
               0.248474121094]
        for i in xrange(9):
            res = lfsr.zero_to_one()
            self.assertLessEqual(res - exp[i], 0.00000000001)
            self.assertGreaterEqual(res - exp[i], -0.00000000001)

        exp = [3.92219712703,
               1.96109856351,
               4.12209399845,
               5.20259171591,
               2.60124792106,
               1.30062396053,
               0.650311980264,
               3.46670070682,
               1.73330241651]
        for i in xrange(9):
            res = lfsr.zero_to_2pi()
            self.assertLessEqual(res - exp[i], 0.00000000001)
            self.assertGreaterEqual(res - exp[i], -0.00000000001)

        exp = [-1.0,
               -1.0,
               -1.0,
               1.0,
               1.0,
               1.0,
               1.0,
               -1.0,
               1.0]
        for i in xrange(9):
            self.assertEquals(lfsr.sign(), exp[i])

def suite():
    suite1 = unittest.makeSuite(ListTestGenerate)
    return unittest.TestSuite((suite1, ))

if __name__ == '__main__':
    unittest.main()
