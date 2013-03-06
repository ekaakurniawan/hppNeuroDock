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
#  - AutoDock 4.2.3 Source Code (trilinterp.cc)
#    http://autodock.scripps.edu

DEBUG = False

class Dock:
    def __init__(self, grid, ligand):
        self.grid = grid
        self.ligand = ligand
        self.elecs = []
        self.elec_total = 0.0
        self.emaps = []
        self.emap_total = 0.0

    def interp3(self):
        lo_x, lo_y, lo_z = self.grid.field.lo.axis
        spacing = self.grid.field.spacing
        
        u = []
        v = []
        w = []
        for atom in self.ligand.atoms:
            u.append((atom.tcoord.x - lo_x) / spacing)
            v.append((atom.tcoord.y - lo_y) / spacing)
            w.append((atom.tcoord.z - lo_z) / spacing)

        u0 = [int(i) for i in u]
        v0 = [int(i) for i in v]
        w0 = [int(i) for i in w]
                
        u1 = [i + 1 for i in u0]
        v1 = [i + 1 for i in v0]
        w1 = [i + 1 for i in w0]

        atom_len = len(self.ligand.atoms)
    
        p0u = [u[i] - float(u0[i]) for i in xrange(atom_len)]
        p0v = [v[i] - float(v0[i]) for i in xrange(atom_len)]
        p0w = [w[i] - float(w0[i]) for i in xrange(atom_len)]

        p1u = [float(u1[i]) - u[i] for i in xrange(atom_len)]
        p1v = [float(v1[i]) - v[i] for i in xrange(atom_len)]
        p1w = [float(w1[i]) - w[i] for i in xrange(atom_len)]

        p000 = [p0u[i] * p0v[i] * p0w[i] for i in xrange(atom_len)]
        p001 = [p0u[i] * p0v[i] * p1w[i] for i in xrange(atom_len)]
        p010 = [p0u[i] * p1v[i] * p0w[i] for i in xrange(atom_len)]
        p011 = [p0u[i] * p1v[i] * p1w[i] for i in xrange(atom_len)]
        p100 = [p1u[i] * p0v[i] * p0w[i] for i in xrange(atom_len)]
        p101 = [p1u[i] * p0v[i] * p1w[i] for i in xrange(atom_len)]
        p110 = [p1u[i] * p1v[i] * p0w[i] for i in xrange(atom_len)]
        p111 = [p1u[i] * p1v[i] * p1w[i] for i in xrange(atom_len)]

        return u0, v0, w0, u1, v1, w1, \
               p000, p001, p010, p011, p100, p101, p110, p111

    def calc_energy(self):
        u0, v0, w0, u1, v1, w1, \
        p000, p001, p010, p011, p100, p101, p110, p111 = self.interp3()
                
        atom_len = len(self.ligand.atoms)
                
        es = []
        for i in xrange(atom_len):
            e = 0.0
            e += p000[i] * self.grid.maps['e'][w1[i]][v1[i]][u1[i]]
            e += p001[i] * self.grid.maps['e'][w1[i]][v1[i]][u0[i]]
            e += p010[i] * self.grid.maps['e'][w1[i]][v0[i]][u1[i]]
            e += p011[i] * self.grid.maps['e'][w1[i]][v0[i]][u0[i]]
            e += p100[i] * self.grid.maps['e'][w0[i]][v1[i]][u1[i]]
            e += p101[i] * self.grid.maps['e'][w0[i]][v1[i]][u0[i]]
            e += p110[i] * self.grid.maps['e'][w0[i]][v0[i]][u1[i]]
            e += p111[i] * self.grid.maps['e'][w0[i]][v0[i]][u0[i]]
            es.append(e)

        ds = []
        for i in xrange(atom_len):
            d = 0.0
            d += p000[i] * self.grid.maps['d'][w1[i]][v1[i]][u1[i]]
            d += p001[i] * self.grid.maps['d'][w1[i]][v1[i]][u0[i]]
            d += p010[i] * self.grid.maps['d'][w1[i]][v0[i]][u1[i]]
            d += p011[i] * self.grid.maps['d'][w1[i]][v0[i]][u0[i]]
            d += p100[i] * self.grid.maps['d'][w0[i]][v1[i]][u1[i]]
            d += p101[i] * self.grid.maps['d'][w0[i]][v1[i]][u0[i]]
            d += p110[i] * self.grid.maps['d'][w0[i]][v0[i]][u1[i]]
            d += p111[i] * self.grid.maps['d'][w0[i]][v0[i]][u0[i]]
            ds.append(d)

        ms = []
        for i, atom in enumerate(self.ligand.atoms):
            m = 0.0
            type = atom.type
            m += p000[i] * self.grid.maps[type][w1[i]][v1[i]][u1[i]]
            m += p001[i] * self.grid.maps[type][w1[i]][v1[i]][u0[i]]
            m += p010[i] * self.grid.maps[type][w1[i]][v0[i]][u1[i]]
            m += p011[i] * self.grid.maps[type][w1[i]][v0[i]][u0[i]]
            m += p100[i] * self.grid.maps[type][w0[i]][v1[i]][u1[i]]
            m += p101[i] * self.grid.maps[type][w0[i]][v1[i]][u0[i]]
            m += p110[i] * self.grid.maps[type][w0[i]][v0[i]][u1[i]]
            m += p111[i] * self.grid.maps[type][w0[i]][v0[i]][u0[i]]
            ms.append(m)

        self.elecs = []
        self.elec_total = 0.0
        for i, atom in enumerate(self.ligand.atoms):
            self.elec = es[i] * atom.charge
            self.elecs.append(self.elec)
            self.elec_total += self.elec

        self.emaps = []
        self.emap_total = 0.0
        for i, atom in enumerate(self.ligand.atoms):
            self.emap = ms[i] + ds[i] * atom.abs_charge
            self.emaps.append(self.emap)
            self.emap_total += self.emap

    def test_print(self):
        for i, atom in enumerate(self.ligand.atoms):
            print "%2s: %2s - %8.3f, %8.3f, %8.3f | %+7.2f | %+7.2f" % \
                (i + 1, atom.type, \
                 atom.tcoord.x, atom.tcoord.y, atom.tcoord.z, \
                 self.emaps[i], self.elecs[i])
