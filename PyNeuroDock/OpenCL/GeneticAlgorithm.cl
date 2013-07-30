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

#define DEBUG 0

// Individual
#define I_TRANS_X_IDX   0
#define I_TRANS_Y_IDX   1
#define I_TRANS_Z_IDX   2
#define I_ROT_A_IDX     3
#define I_ROT_B_IDX     4
#define I_ROT_C_IDX     5
#define I_ROT_D_IDX     6
#define I_TOR_START_IDX 7

// Population random number
#define PRN_P1_IDX                  0   // Parent1
#define PRN_P2_IDX                  1   // Parent2
#define PRN_MUTATION_CHANCE_IDX     2   // Overall mutation probability
#define PRN_TRANS_IDX               3   // Translation mutation probability
#define PRN_ROT_IDX                 4   // Rotation mutation probability
#define PRN_TOR_START_IDX           5   // Torsion mutation probabilities
#define PRN_POSE_TRANS_X_IDX        (PRN_TOR_START_IDX + ttl_torsions[0] + 0)
#define PRN_POSE_TRANS_Y_IDX        (PRN_TOR_START_IDX + ttl_torsions[0] + 1)
#define PRN_POSE_TRANS_Z_IDX        (PRN_TOR_START_IDX + ttl_torsions[0] + 2)
#define PRN_POSE_ROT_A_IDX          (PRN_TOR_START_IDX + ttl_torsions[0] + 3)
#define PRN_POSE_ROT_B_IDX          (PRN_TOR_START_IDX + ttl_torsions[0] + 4)
#define PRN_POSE_ROT_C_IDX          (PRN_TOR_START_IDX + ttl_torsions[0] + 5)
#define PRN_POSE_ROT_D_IDX          (PRN_TOR_START_IDX + ttl_torsions[0] + 6)
#define PRN_POSE_TOR_START_IDX      (PRN_TOR_START_IDX + ttl_torsions[0] + 7)

// Crossover modes
// Separate probabilities for translation/rotation genes
#define CM_SEPARATE (long)0
// Combined probability for translation/rotation genes
#define CM_COMBINE  (long)1

// Construct individuals from previously generated random numbers
__kernel void construct_individuals(__global const double *lo_grid,
                                    __global const double *dist_grid,
                                    __global const long *dna_size,

                                    __global double *individuals)
{
    // Individual ID
    long i_id = get_global_id(0);
    // DNA (gene) ID
    long st_idx = i_id * dna_size[0];
    double two_pi = 2 * M_PI;

    // Translation genes
    individuals[st_idx + I_TRANS_X_IDX] = lo_grid[0] +
                                          (individuals[st_idx + I_TRANS_X_IDX] * dist_grid[0]);
    individuals[st_idx + I_TRANS_Y_IDX] = lo_grid[1] +
                                          (individuals[st_idx + I_TRANS_Y_IDX] * dist_grid[1]);
    individuals[st_idx + I_TRANS_Z_IDX] = lo_grid[2] +
                                          (individuals[st_idx + I_TRANS_Z_IDX] * dist_grid[2]);
    // Rotation genes
    // x0 = rng.zero_to_one()
    double x0 = individuals[st_idx + I_ROT_A_IDX];
    // t1 = rng.zero_to_2pi()
    double t1 = individuals[st_idx + I_ROT_B_IDX] * two_pi;
    // t2 = rng.zero_to_2pi()
    double t2 = individuals[st_idx + I_ROT_C_IDX] * two_pi;
    double r1 = sqrt(1.0 - x0);
    double r2 = sqrt(x0);
    individuals[st_idx + I_ROT_B_IDX] = sin(t1) * r1;
    individuals[st_idx + I_ROT_C_IDX] = cos(t1) * r1;
    individuals[st_idx + I_ROT_D_IDX] = sin(t2) * r2;
    individuals[st_idx + I_ROT_A_IDX] = cos(t2) * r2;

    // Torsion genes
    for (long tor_idx = I_TOR_START_IDX; tor_idx < dna_size[0]; tor_idx++) {
        // rng.neg_pi_to_pi()
        individuals[st_idx + tor_idx] = (individuals[st_idx + tor_idx] - 0.5) *
                                        two_pi;
    }

    if (DEBUG == 1) {
        if (i_id == 0 || i_id == 20 || i_id == 149) {
            individuals[i_id * 19 + 0]  =  2.056477;
            individuals[i_id * 19 + 1]  =  5.846611;
            individuals[i_id * 19 + 2]  =  -7.245407;
            individuals[i_id * 19 + 3]  =  0.532211;
            individuals[i_id * 19 + 4]  =  0.379383;
            individuals[i_id * 19 + 5]  =  0.612442;
            individuals[i_id * 19 + 6]  =  0.444674;
            individuals[i_id * 19 + 7]  =  radians(-122.13);
            individuals[i_id * 19 + 8]  =  radians(-179.41);
            individuals[i_id * 19 + 9]  =  radians(-141.59);
            individuals[i_id * 19 + 10] =  radians(177.29);
            individuals[i_id * 19 + 11] =  radians(-179.46);
            individuals[i_id * 19 + 12] =  radians(-9.31);
            individuals[i_id * 19 + 13] =  radians(132.37);
            individuals[i_id * 19 + 14] =  radians(-89.19);
            individuals[i_id * 19 + 15] =  radians(78.43);
            individuals[i_id * 19 + 16] =  radians(22.22);
            individuals[i_id * 19 + 17] =  radians(71.37);
            individuals[i_id * 19 + 18] =  radians(59.52);
        }
    }
}

__kernel void calc_chances(__global const double *e_totals,
                           __global const long *normalizer,
                           __global const long *max_inherited_prob,

                           __global long *chances)
{
    // Individual ID
    long i_id = get_global_id(0);
    double score = e_totals[i_id] / (double)normalizer[0];

    // To handle different interpretation of Python infinity by different
    // processors.
    // - NVIDIA GeForce GT 650M catches it as 9223372036854775807
    // - Intel Core i7 catches it as         -9223372036854775808
    if (e_totals[i_id] == INFINITY) {
        chances[i_id] = 1;
        return;
    }
    if (score < 0.0) {
        chances[i_id] = max_inherited_prob[0];
        return;
    }
    double power = log(score);
    if ((long)power < max_inherited_prob[0]) {
        chances[i_id] = (long)((double)max_inherited_prob[0] - power);
    } else {
        chances[i_id] = 1;
    }
}

//TODO: Update the new individuals into individials array instead of
//      new_individuals array
__kernel void reproduce(__global const long *population_size,
                        __global const long *chances,
                        __global const long *ttl_reproduction_rns,
                        __global const double *reproduction_rns,

                        __global const long *dna_size,
                        __global const double *individuals,

                        __global const long *crossover_translation_mode,
                        __global const long *crossover_rotation_mode,
                        __global const double *crossover_probability,

                        __global const double *mutation_chance,
                        __global const double *mutation_probability,
                        __global const long *ttl_torsions,
                        __global const double *lo_grid,
                        __global const double *dist_grid,

                        __global long *chances_sum,
                        __global double *dna1,
                        __global double *dna2,

                        __global double *new_individuals)
{
    // Individual ID
    long i_id = get_global_id(0);

    // Pick parents
    long ttl_chances = 0;
    for (long i = 0; i < population_size[0]; i++) {
        ttl_chances += chances[i];
        // Prefix-sum
        chances_sum[i] = ttl_chances;
    }
    ttl_chances += 1;
    // Start index of reproduction random number
    long start_r_rns_idx = i_id * ttl_reproduction_rns[0];
    long p1_loc = (long)(reproduction_rns[start_r_rns_idx + PRN_P1_IDX] *
                         ttl_chances);
    long p2_loc = (long)(reproduction_rns[start_r_rns_idx + PRN_P2_IDX] *
                         ttl_chances);
    // p1 location has to be smaller than or equal to p2
    if (p2_loc < p1_loc) {
        long ptmp_loc = p1_loc;
        p1_loc = p2_loc;
        p2_loc = ptmp_loc;
    }
    // Get parent 1 and 2 IDs
    long p1_id, p2_id;
    for (long i = 0; i < population_size[0]; i++) {
        if (p1_loc <= chances_sum[i]) {
            p1_id = i;
            break;
        }
    }
    for (long i = p1_id; i < population_size[0]; i++) {
        if (p2_loc <= chances_sum[i]) {
            p2_id = i;
            break;
        }
    }
    // Get parent 1 and 2 DNAs
    long start_src_idx;
    start_src_idx = p1_id * dna_size[0];
    long start_dna_idx = i_id * dna_size[0];
    for (long i = 0; i < dna_size[0]; i++) {
        dna1[start_dna_idx + i] = individuals[start_src_idx + i];
    }
    start_src_idx = p2_id * dna_size[0];
    for (long i = 0; i < dna_size[0]; i++) {
        dna2[start_dna_idx + i] = individuals[start_src_idx + i];
    }
    long start_dst_idx = i_id * dna_size[0];
    // Crossover
    // If parent 1 and 2 IDs point to a same DNA, skip crossover
    if (p1_id != p2_id) {
        // Translation genes
        if (crossover_translation_mode[0] == CM_SEPARATE) {
            if (new_individuals[start_dst_idx + I_TRANS_X_IDX] > crossover_probability[0]) {
                dna1[start_dna_idx + I_TRANS_X_IDX] = dna2[start_dna_idx + I_TRANS_X_IDX];
            }
            if (new_individuals[start_dst_idx + I_TRANS_Y_IDX] > crossover_probability[0]) {
                dna1[start_dna_idx + I_TRANS_Y_IDX] = dna2[start_dna_idx + I_TRANS_Y_IDX];
            }
            if (new_individuals[start_dst_idx + I_TRANS_Z_IDX] > crossover_probability[0]) {
                dna1[start_dna_idx + I_TRANS_Z_IDX] = dna2[start_dna_idx + I_TRANS_Z_IDX];
            }
        } else { // CM_COMBINE
            if (new_individuals[start_dst_idx + I_TRANS_X_IDX] > crossover_probability[0]) {
                for (long i = 0; i < 3; i++) {
                    dna1[start_dna_idx + I_TRANS_X_IDX + i] = dna2[start_dna_idx + I_TRANS_X_IDX + i];
                }
            }
        }
        // Rotation genes
        if (crossover_rotation_mode[0] == CM_SEPARATE) {
            if (new_individuals[start_dst_idx + I_ROT_A_IDX] > crossover_probability[0]) {
                dna1[start_dna_idx + I_ROT_A_IDX] = dna2[start_dna_idx + I_ROT_A_IDX];
            }
            if (new_individuals[start_dst_idx + I_ROT_B_IDX] > crossover_probability[0]) {
                dna1[start_dna_idx + I_ROT_B_IDX] = dna2[start_dna_idx + I_ROT_B_IDX];
            }
            if (new_individuals[start_dst_idx + I_ROT_C_IDX] > crossover_probability[0]) {
                dna1[start_dna_idx + I_ROT_C_IDX] = dna2[start_dna_idx + I_ROT_C_IDX];
            }
            if (new_individuals[start_dst_idx + I_ROT_D_IDX] > crossover_probability[0]) {
                dna1[start_dna_idx + I_ROT_D_IDX] = dna2[start_dna_idx + I_ROT_D_IDX];
            }
        } else { // CM_COMBINE
            if (new_individuals[start_dst_idx + I_ROT_A_IDX] > crossover_probability[0]) {
                for (long i = 0; i < 4; i++) {
                    dna1[start_dna_idx + I_ROT_A_IDX + i] = dna2[start_dna_idx + I_ROT_A_IDX + i];
                }
            }
        }
        // Torsion genes
        for (long i = 0; i < ttl_torsions[0]; i++) {
            if (new_individuals[start_dst_idx + I_TOR_START_IDX + i] > crossover_probability[0]) {
                dna1[start_dna_idx + I_TOR_START_IDX + i] = dna2[start_dna_idx + I_TOR_START_IDX + i];
            }
        }
    }
    // Mutation
    if (reproduction_rns[start_r_rns_idx + PRN_MUTATION_CHANCE_IDX] < mutation_chance[0]) {
        double two_pi = 2 * M_PI;
        // Translation genes
        if (reproduction_rns[start_r_rns_idx + PRN_TRANS_IDX] < mutation_probability[0]) {
            // Translation genes. Starts from index 0 to 2
            dna1[start_dna_idx + I_TRANS_X_IDX] = lo_grid[0] +
                                                  (reproduction_rns[start_r_rns_idx + PRN_POSE_TRANS_X_IDX] * dist_grid[0]);
            dna1[start_dna_idx + I_TRANS_Y_IDX] = lo_grid[1] +
                                                  (reproduction_rns[start_r_rns_idx + PRN_POSE_TRANS_Y_IDX] * dist_grid[1]);
            dna1[start_dna_idx + I_TRANS_Z_IDX] = lo_grid[2] +
                                                  (reproduction_rns[start_r_rns_idx + PRN_POSE_TRANS_Z_IDX] * dist_grid[2]);
        }
        // Rotation genes
        if (reproduction_rns[start_r_rns_idx + PRN_ROT_IDX] < mutation_probability[0]) {
            // x0 = rng.zero_to_one()
            double x0 = reproduction_rns[start_r_rns_idx + PRN_POSE_ROT_A_IDX];
            // t1 = rng.zero_to_2pi()
            double t1 = reproduction_rns[start_r_rns_idx + PRN_POSE_ROT_B_IDX] * two_pi;
            // t2 = rng.zero_to_2pi()
            double t2 = reproduction_rns[start_r_rns_idx + PRN_POSE_ROT_C_IDX] * two_pi;
            double r1 = sqrt(1.0 - x0);
            double r2 = sqrt(x0);
            dna1[start_dna_idx + I_ROT_B_IDX] = sin(t1) * r1;
            dna1[start_dna_idx + I_ROT_C_IDX] = cos(t1) * r1;
            dna1[start_dna_idx + I_ROT_D_IDX] = sin(t2) * r2;
            dna1[start_dna_idx + I_ROT_A_IDX] = cos(t2) * r2;
        }
        // Tortion genes
        for (long i = 0; i < ttl_torsions[0]; i++) {
            if (reproduction_rns[start_r_rns_idx + PRN_TOR_START_IDX + i] < mutation_probability[0]) {
                // rng.neg_pi_to_pi()
                dna1[start_dna_idx + I_TOR_START_IDX + i] = (reproduction_rns[start_r_rns_idx + PRN_POSE_TOR_START_IDX + i] - 0.5) *
                                                            two_pi;
            }
        }
    }
    // Commit the new individual
    for (long i = 0; i < dna_size[0]; i++) {
        new_individuals[start_dst_idx + i] = dna1[start_dna_idx + i];
    }
}

