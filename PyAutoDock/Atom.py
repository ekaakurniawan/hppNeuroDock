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
#  - AutoDock 4.2.3 Source Code (mdist.h)
#    http://autodock.scripps.edu

from Axis3 import Axis3

class Atom:
    NUM_SUPPORTED_ATOM = 7
    C, N, O, H, XX, P, S = range(NUM_SUPPORTED_ATOM)

    def __init__(self, id = 0, type = "?", tcoord = Axis3(0.0, 0.0, 0.0), \
                 charge = 0.0):
        self.id = id
        self.type = type
        self.tcoord = tcoord
        self.charge = charge

class Branch:
    def __init__(self, anchor_id = 0, link_id = 0, atom_ids = []):
        self.anchor_id = anchor_id
        self.link_id = link_id
        self.atom_ids = atom_ids

class Bond:
    DISTANCE_TOLERANCE = 0.1

    def __init__(self):
        pass

    # Get minimum to maximum range of atom to atom distance
    @staticmethod
    def get_minmax_distance():
        minmax_distance = [[[0, 0] for x in xrange(Atom.NUM_SUPPORTED_ATOM)] \
                           for x in xrange(Atom.NUM_SUPPORTED_ATOM)]

        def set_minmax_distance(atom1, atom2, val):
            minmax_distance[atom1][atom2][0] = val[0] - Bond.DISTANCE_TOLERANCE
            minmax_distance[atom2][atom1][0] = minmax_distance[atom1][atom2][0]

            minmax_distance[atom1][atom2][1] = val[1] + Bond.DISTANCE_TOLERANCE
            minmax_distance[atom2][atom1][1] = minmax_distance[atom1][atom2][1]

        # Following values, unless otherwise stated, are taken from AutoDock,
        # taken from "Handbook of Chemistry and Physics" 44th edition
        set_minmax_distance(Atom.C,  Atom.C,  [1.20,    1.545])  # mindist[C][C] = 1.20, p. 3510 ; maxdist[C][C] = 1.545, p. 3511
        set_minmax_distance(Atom.C,  Atom.N,  [1.1,     1.479])  # mindist[C][N] = 1.1, p. 3510 ; maxdist[C][N] = 1.479, p. 3511
        set_minmax_distance(Atom.C,  Atom.O,  [1.15,    1.47])   # mindist[C][O] = 1.15, p. 3510 ; maxdist[C][O] = 1.47, p. 3512
        set_minmax_distance(Atom.C,  Atom.H,  [1.022,   1.12])   # p. 3518, p. 3517
        set_minmax_distance(Atom.C,  Atom.XX, [0.9,     1.545])  # mindist[C][XX] = 0.9, AutoDock 3 defaults ; maxdist[C][XX] = 1.545, p. 3511
        set_minmax_distance(Atom.C,  Atom.P,  [1.85,    1.89])   # mindist[C][P] = 1.85, p. 3510 ; maxdist[C][P] = 1.89, p. 3510
        set_minmax_distance(Atom.C,  Atom.S,  [1.55,    1.835])  # mindist[C][S] = 1.55, p. 3510 ; maxdist[C][S] = 1.835, p. 3512
        set_minmax_distance(Atom.N,  Atom.N,  [1.0974,  1.128])  # mindist[N][N] = 1.0974, p. 3513 ; maxdist[N][N] = 1.128, p. 3515
        set_minmax_distance(Atom.N,  Atom.O,  [1.0619,  1.25])   # mindist[N][O] = 1.0975, p. 3515 ; maxdist[N][O] = 1.128, p. 3515
        set_minmax_distance(Atom.N,  Atom.H,  [1.004,   1.041])  # mindist[N][H] = 1.004, p. 3516 ; maxdist[N][H] = 1.041, p. 3515
        set_minmax_distance(Atom.N,  Atom.XX, [0.9,     1.041])  # mindist[N][XX] = 0.9, AutoDock 3 defaults ; maxdist[N][XX] = 1.041, p. 3515
        set_minmax_distance(Atom.N,  Atom.P,  [1.4910,  1.4910]) # mindist[N][P] = 1.4910, p. 3515 ; maxdist[N][P] = 1.4910, p. 3515
        set_minmax_distance(Atom.N,  Atom.S,  [1.58,    1.672])  # mindist[N][S] = 1.58, 1czm.pdb sulfonamide ; maxdist[N][S] = 1.672, J. Chem. SOC., Dalton Trans., 1996, Pages 4063-4069
        set_minmax_distance(Atom.O,  Atom.O,  [1.208,   1.51])   # p.3513, p.3515
        set_minmax_distance(Atom.O,  Atom.H,  [0.955,   1.0289]) # mindist[O][H] = 0.955, p. 3515 ; maxdist[O][H] = 1.0289, p. 3515
        set_minmax_distance(Atom.O,  Atom.XX, [0.955,   2.1])    # AutoDock 3 defaults
        set_minmax_distance(Atom.O,  Atom.P,  [1.36,    1.67])   # mindist[O][P] = 1.36, p. 3516 ; maxdist[O][P] = 1.67, p. 3517
        set_minmax_distance(Atom.O,  Atom.S,  [1.41,    1.47])   # p. 3517, p. 3515
        set_minmax_distance(Atom.H,  Atom.H,  [100.0,  -100.0])  # impossible values to prevent such bonds from forming.
        set_minmax_distance(Atom.H,  Atom.XX, [0.9,     1.5])    # AutoDock 4 defaults
        set_minmax_distance(Atom.H,  Atom.P,  [1.40,    1.44])   # mindist[H][P] = 1.40, p. 3515 ; maxdist[H][P] = 1.44, p. 3515
        set_minmax_distance(Atom.H,  Atom.S,  [1.325,   1.3455]) # mindist[H][S] = 1.325, p. 3518 ; maxdist[H][S] = 1.3455, p. 3516
        set_minmax_distance(Atom.XX, Atom.XX, [0.9,     2.1])    # AutoDock 3 defaults
        set_minmax_distance(Atom.XX, Atom.P,  [0.9,     2.1])    # AutoDock 3 defaults
        set_minmax_distance(Atom.XX, Atom.S,  [1.325,   2.1])    # mindist[XX][S] = 1.325, p. 3518 ; maxdist[XX][S] = 2.1, AutoDock 3 defaults
        set_minmax_distance(Atom.P,  Atom.P,  [2.18,    2.23])   # mindist[P][P] = 2.18, p. 3513 ; maxdist[P][P] = 2.23, p. 3513
        set_minmax_distance(Atom.P,  Atom.S,  [1.83,    1.88])   # mindist[P][S] = 1.83, p. 3516 ; maxdist[P][S] = 1.88, p. 3515
        set_minmax_distance(Atom.S,  Atom.S,  [2.03,    2.05])   # mindist[S][S] = 2.03, p. 3515 ; maxdist[S][S] = 2.05, p. 3515

        return minmax_distance

    # Construct atom-to-atom covalent bond matrix
    @staticmethod
    def construct_bond_matrix():
        pass

    # Construct atom-to-atom non-covalent bond matrix by considering 1-1, 1-2,
    # 1-3 and 1-4 interactions
    @staticmethod
    def construct_non_bond_matrix():
        pass

#bar - start
#print Bond.get_minmax_distance()
#bar - stop
