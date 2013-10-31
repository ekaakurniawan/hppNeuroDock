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
#  - http://en.wikipedia.org/wiki/Linear_feedback_shift_register
#  - Efficient Shift Registers, LFSR Counters, and Long Pseudo-Random Sequence
#    Generators
#    http://www.xilinx.com/support/documentation/application_notes/xapp052.pdf
#  - IEEE-754 Analysis Calculator
#    http://babbage.cs.qc.edu/IEEE-754/
#  - Xilinx LogiCORE IP Floating-Point Operator
#    http://www.xilinx.com/support/documentation/ip_documentation/floating_point_ds335.pdf

# Linear Feedback Shift Register (LFSR)
# -------------------------------------
# Generates number in between 0 and 2^n where n is the bit length. Note that 0
# and 2^n are exclusive.


import Constants as const

class LFSR:
    xnor_shifts_list = {
        4:  [0, 1],             # 4, 3
        8:  [0, 2, 3, 4],       # 8, 6, 5, 4
        16: [0, 1, 3, 12],      # 16, 15, 13, 4
        32: [0, 10, 30, 31],    # 32, 22, 2, 1
        48: [0, 1, 27, 28],     # 48, 47, 21, 20
        64: [0, 1, 3, 4],       # 64, 63, 61, 60
    }

    def __init__(self, lfsr = 1070, bit_len = 16):
        self.lfsr = lfsr
        self.bit_len = bit_len
        self.bit_len_dec = self.bit_len - 1
        self.xnor_shifts = self.xnor_shifts_list[bit_len]
        self.denominator = 2.0 ** self.bit_len

    def generate(self):
        bit = self.lfsr >> self.xnor_shifts[0]
        for xnor_shift in self.xnor_shifts[1:]:
            bit = ~(bit ^ (self.lfsr >> xnor_shift))
        bit &= 1
        self.lfsr = (self.lfsr >> 1) | (bit << self.bit_len_dec)
        return self.lfsr

    # This implementation follows genunf function from AutoDock that takes two
    # parameters (low and high) and returns uniform distribution in between
    # them. Note that the low and high boundaries are exclusive. In this case,
    # the low value is 0.0 and the high value is 1.0.
    def zero_to_one(self):
        # For fixed bit length, the denominator is implemented as constant
        return (self.generate() / self.denominator)

    # This implementation follows genunf function from AutoDock that takes two
    # parameters (low and high) and returns uniform distribution in between
    # them. Note that the low and high boundaries are exclusive. In this case,
    # the low value is 0.0 and the high value is 2 x pi.
    def zero_to_2pi(self):
        # For fixed bit length, the denominator is implemented as constant
        return (self.zero_to_one() * const.TWOPI)

    # Randomly returns floating point number in between negative PI and
    # positive PI
    def neg_pi_to_pi(self):
        # For fixed bit length, the denominator is implemented as constant
        return ((self.zero_to_one() - 0.5) * const.TWOPI)

    # Randomly returns -1.0 or 1.0
    # Based on IEEE-754 Analysis Calculator,
    #  - Single precision (32-bit) for -1.0 is 0xBF800000
    #                              for  1.0 is 0x3F800000
    #  - Double precision (64-bit) for -1.0 is 0xBFF0000000000000
    #                              for  1.0 is 0x3FF0000000000000
    # Based on Xilinx LogiCORE IP Floating-Point Operator, sign bit (s) is
    # placed at the MSB of each floating point number. s = 0 determines positive
    # number whereas s = 1 determines the negative.
    def sign(self):
        # Can be implemented by checking the MSB of LFSR
        if self.generate() < ((2 ** self.bit_len) / 2):
            return -1.0
        else:
            return 1.0

