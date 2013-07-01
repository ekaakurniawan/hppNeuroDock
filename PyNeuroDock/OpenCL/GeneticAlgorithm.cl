// Copyright (C) 2013 by Eka A. Kurniawan
// eka.a.kurniawan(ta)gmail(tod)com
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; version 2 of the License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the
// Free Software Foundation, Inc.,
// 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

// Construct individuals from previously generated random numbers
__kernel void construct_individuals(__global double *individuals,
                                    __global const double *lo_grid,
                                    __global const double *dist_grid,
                                    __global const long *ttl_torsions)
{
    long i = get_global_id(0);          // individuals
    long j = 3 + 4 + ttl_torsions[0];   // DNA (genes)
    long st_idx = i * j;
    double two_pi = 2 * M_PI;

    // Translation genes. Starts from index 0 to 2
    individuals[st_idx + 0] = lo_grid[0] + 
                              (individuals[st_idx + 0] * dist_grid[0]);
    individuals[st_idx + 1] = lo_grid[1] +
                              (individuals[st_idx + 1] * dist_grid[1]);
    individuals[st_idx + 2] = lo_grid[2] +
                              (individuals[st_idx + 2] * dist_grid[2]);
    // Rotation genes. Starts from index 3 to 6
    //x0 = rng.zero_to_one()
    double x0 = individuals[st_idx + 3];
    // t1 = rng.zero_to_2pi()
    double t1 = individuals[st_idx + 4] * two_pi;
    // t2 = rng.zero_to_2pi()
    double t2 = individuals[st_idx + 5] * two_pi;
    double r1 = sqrt(1.0 - x0);
    double r2 = sqrt(x0);
    individuals[st_idx + 4] = sin(t1) * r1;
    individuals[st_idx + 5] = cos(t1) * r1;
    individuals[st_idx + 6] = sin(t2) * r2;
    individuals[st_idx + 3] = cos(t2) * r2;

    // Torsion genes. Starts from index 7 to the last index
    for (long tor_idx = 7; tor_idx < j; tor_idx++) {
        // rng.neg_pi_to_pi()
        individuals[st_idx + tor_idx] = (individuals[st_idx + tor_idx] - 0.5) *
                                        two_pi;
    }

    //bar - start
    if (i == 0 || i == 20 || i == 149) {
        individuals[i * 19 + 0]  =  2.056477;
        individuals[i * 19 + 1]  =  5.846611;
        individuals[i * 19 + 2]  =  -7.245407;
        individuals[i * 19 + 3]  =  0.532211;
        individuals[i * 19 + 4]  =  0.379383;
        individuals[i * 19 + 5]  =  0.612442;
        individuals[i * 19 + 6]  =  0.444674;
        individuals[i * 19 + 7]  =  radians(-122.13);
        individuals[i * 19 + 8]  =  radians(-179.41);
        individuals[i * 19 + 9]  =  radians(-141.59);
        individuals[i * 19 + 10] =  radians(177.29);
        individuals[i * 19 + 11] =  radians(-179.46);
        individuals[i * 19 + 12] =  radians(-9.31);
        individuals[i * 19 + 13] =  radians(132.37);
        individuals[i * 19 + 14] =  radians(-89.19);
        individuals[i * 19 + 15] =  radians(78.43);
        individuals[i * 19 + 16] =  radians(22.22);
        individuals[i * 19 + 17] =  radians(71.37);
        individuals[i * 19 + 18] =  radians(59.52);
    }
    //bar - stop
}

