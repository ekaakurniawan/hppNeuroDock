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

// References:
// - Converting a Triangular Matrix to an Array
//   http://jamesmccaffrey.wordpress.com/2010/05/14/converting-a-triangular-matrix-to-an-array/

#pragma OPENCL EXTENSION cl_khr_fp64 : enable

//TODO: Add few characters in front of definition names to act like namespace
// Geometry index
#define TRANSLATION_START_IDX   0
#define ROTATION_START_IDX      3
#define TORSION_START_IDX       7

// Atom properties
#define ATOM_TYPE_IDX       0
#define ATOM_CHARGE_IDX     1

// Non-bond properties
#define ATOM_ID1_IDX        0
#define ATOM_TYPE1_IDX      1
#define ATOM_ID2_IDX        2
#define ATOM_TYPE2_IDX      3
#define NON_BOND_TYPE_IDX   4
#define DESOLV_IDX          5
#define Q1Q2_IDX            6

// Bond properties
#define NS_INTL_1_IDX               0
#define NS_EL_1_IDX                 1
#define RMIN_ELEC2_IDX              2
#define SQA_DIV_IDX                 3
#define NBC2_IDX                    4
#define SCALE_1_4_INTERACTIONS_IDX  5


__kernel void rotate_branches(__global const long *ttl_torsions,
                              __global const double *individuals,

                              __global const long *longest_branch,
                              __global const long *branches_rot_anchor_buf,
                              __global const long *branches_rot_link_buf,
                              __global const long *branches_rot_size_buf,
                              __global const long *branches_rot_seq_buf,

                              __global const long *ttl_poses,
                              __global double *poses)
{
    // Thread ID
    long thread_id = get_global_id(0);
    // Pose ID or individual ID
    long pose_id = thread_id / longest_branch[0];
    // Branch rotation, j axis
    long br_j = thread_id % longest_branch[0];
    // Branch rotation, i axis or torsion index
    for (long br_i = 0; br_i < ttl_torsions[0]; br_i++) {
        // If atom ID equals 0, skip the conversion
        long atom_tcoord_id = branches_rot_seq_buf[(br_i * longest_branch[0]) +
                                                           br_j];
        if (atom_tcoord_id == 0) {
            barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
            continue;
        }
        // Torsion angle
        double tor_angle = individuals[(pose_id * (3 + 4 + ttl_torsions[0])) +
                                       TORSION_START_IDX + br_i];
        // Torsion axis
        double anchor_tcoord[3];
        double link_tcoord[3];
        double3 tor_axis;
        long anchor_tcoord_id = branches_rot_anchor_buf[br_i];
        long link_tcoord_id =  branches_rot_link_buf[br_i];
        for (long i = 0; i < 3; i++) {
            anchor_tcoord[i] = poses[(anchor_tcoord_id * ttl_poses[0] * 3) +
                                     (pose_id * 3) + i];
            link_tcoord[i] = poses[(link_tcoord_id * ttl_poses[0] * 3) +
                                   (pose_id * 3) + i];
            tor_axis[i] = anchor_tcoord[i] - link_tcoord[i];
        }
        // Rotation in quaternion
        double rotation[4];
        double half_tor_angle = tor_angle / 2;
        rotation[0] = cos(half_tor_angle);
        tor_axis = normalize(tor_axis);
        double s = sin(half_tor_angle);
        rotation[1] = tor_axis[0] * s;
        rotation[2] = tor_axis[1] * s;
        rotation[3] = tor_axis[2] * s;
        // Transform
        double atom_tcoord[3];
        for (long i = 0; i < 3; i++) {
            atom_tcoord[i] = poses[(atom_tcoord_id * ttl_poses[0] * 3) +
                                   (pose_id * 3) + i] -
                             link_tcoord[i];
        }
        double a = rotation[0];
        double b = rotation[1];
        double c = rotation[2];
        double d = rotation[3];

        double db = b + b;
        double dc = c + c;
        double dd = d + d;

        double a_db = a * db;
        double a_dc = a * dc;
        double a_dd = a * dd;

        double b_db_i = 1.0 - (b * db);

        double c_db = c * db;
        double c_dc = c * dc;
        double c_dc_i = 1.0 - c_dc;

        double d_db = d * db;
        double d_dc = d * dc;
        double d_dd = d * dd;

        double r_xx = c_dc_i - d_dd;
        double r_xy = c_db   + a_dd;
        double r_xz = d_db   - a_dc;
        double r_yx = c_db   - a_dd;
        double r_yy = b_db_i - d_dd;
        double r_yz = d_dc   + a_db;
        double r_zx = d_db   + a_dc;
        double r_zy = d_dc   - a_db;
        double r_zz = b_db_i - c_dc;

        double new_atom_tcoord[3];
        new_atom_tcoord[0] =  atom_tcoord[0] * r_xx;
        new_atom_tcoord[0] += atom_tcoord[1] * r_xy;
        new_atom_tcoord[0] += atom_tcoord[2] * r_xz;

        new_atom_tcoord[1] =  atom_tcoord[0] * r_yx;
        new_atom_tcoord[1] += atom_tcoord[1] * r_yy;
        new_atom_tcoord[1] += atom_tcoord[2] * r_yz;

        new_atom_tcoord[2] =  atom_tcoord[0] * r_zx;
        new_atom_tcoord[2] += atom_tcoord[1] * r_zy;
        new_atom_tcoord[2] += atom_tcoord[2] * r_zz;

        for (long i = 0; i < 3; i++) {
            poses[(atom_tcoord_id * ttl_poses[0] * 3) + (pose_id * 3) + i] =
                new_atom_tcoord[i] + link_tcoord[i];
        }

        barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
    }
}

__kernel void transform_ligand_root(__global const long *ttl_ligand_atoms,
                                    __global const double *individuals,
                                    __global const long *ttl_torsions,
                                    __global const long *ttl_poses,
                                    __global double *poses)
{
    long thread_id = get_global_id(0);
    // Atom ID starts from index 1
    long atom_id = (thread_id / ttl_poses[0]) + 1;
    // Pose ID or individual ID
    long pose_id = thread_id % ttl_poses[0];
    // Translation
    double translation[3];
    for (long i = 0; i < 3; i++) {
        translation[i] = individuals[(pose_id * (3 + 4 + ttl_torsions[0])) +
                                     TRANSLATION_START_IDX + i];
    }
    // Rotation
    double rotation[4];
    for (long i = 0; i < 4; i++) {
        rotation[i] = individuals[(pose_id * (3 + 4 + ttl_torsions[0])) +
                                  ROTATION_START_IDX + i];
    }
    // Atom coordinate
    double atom_tcoord[3];
    for (long i = 0; i < 3; i++) {
        atom_tcoord[i] = poses[(atom_id * ttl_poses[0] * 3) +
                               (pose_id * 3) + i];
    }
    // Transform
    double a = rotation[0];
    double b = rotation[1];
    double c = rotation[2];
    double d = rotation[3];

    double db = b + b;
    double dc = c + c;
    double dd = d + d;

    double a_db = a * db;
    double a_dc = a * dc;
    double a_dd = a * dd;

    double b_db_i = 1.0 - (b * db);

    double c_db = c * db;
    double c_dc = c * dc;
    double c_dc_i = 1.0 - c_dc;

    double d_db = d * db;
    double d_dc = d * dc;
    double d_dd = d * dd;

    double r_xx = c_dc_i - d_dd;
    double r_xy = c_db   + a_dd;
    double r_xz = d_db   - a_dc;
    double r_yx = c_db   - a_dd;
    double r_yy = b_db_i - d_dd;
    double r_yz = d_dc   + a_db;
    double r_zx = d_db   + a_dc;
    double r_zy = d_dc   - a_db;
    double r_zz = b_db_i - c_dc;

    double new_atom_tcoord[3];
    new_atom_tcoord[0] =  atom_tcoord[0] * r_xx;
    new_atom_tcoord[0] += atom_tcoord[1] * r_xy;
    new_atom_tcoord[0] += atom_tcoord[2] * r_xz;

    new_atom_tcoord[1] =  atom_tcoord[0] * r_yx;
    new_atom_tcoord[1] += atom_tcoord[1] * r_yy;
    new_atom_tcoord[1] += atom_tcoord[2] * r_yz;

    new_atom_tcoord[2] =  atom_tcoord[0] * r_zx;
    new_atom_tcoord[2] += atom_tcoord[1] * r_zy;
    new_atom_tcoord[2] += atom_tcoord[2] * r_zz;

    for (long i = 0; i < 3; i++) {
        poses[(atom_id * ttl_poses[0] * 3) + (pose_id * 3) + i] =
            new_atom_tcoord[i] + translation[i];
    }
}

__kernel void calc_inter_energy(__global const long *ttl_poses,
                                __global const double *lo_grid,
                                __global const double *hi_grid,
                                __global const double *field_spacing,
                                __global const double *poses,

                                __global const long *num_points1,
                                __global const long *ttl_maps,
                                __global const long *electrostatic_lut,
                                __global const long *desolvation_lut,
                                __global const long *atom_type_map_lut,
                                __global const double *maps,

                                __global const long *ttl_atom_properties,
                                __global const double *atoms_properties,
                                __global const long *protein_ignore_inter,
                                __global const long *ttl_protein_ignore_inter,

                                __global double *elecs,
                                __global double *emaps)
{
    long thread_id = get_global_id(0);
    // Atom ID starts from index 1
    long atom_id = (thread_id / ttl_poses[0]) + 1;
    // Pose ID or individual ID
    long pose_id = thread_id % ttl_poses[0];
    // Exclude the atom and the first atom branching out of root from
    // intermolecular energy calculation
    for (long i = 0; i < ttl_protein_ignore_inter[0]; i++) {
        if (atom_id == protein_ignore_inter[i]) {
            elecs[(atom_id * ttl_poses[0]) + pose_id] = 0.0;
            emaps[(atom_id * ttl_poses[0]) + pose_id] = 0.0;
            return;
        }
    }
    // Atom coordinate
    double atom_tcoord[3];
    for (long i = 0; i < 3; i++) {
        atom_tcoord[i] = poses[(atom_id * ttl_poses[0] * 3) +
                               (pose_id * 3) + i];
    }
    // Check out of grid
    if (atom_tcoord[0] <= lo_grid[0] ||
        atom_tcoord[1] <= lo_grid[1] ||
        atom_tcoord[2] <= lo_grid[2] ||
        atom_tcoord[0] >= hi_grid[0] ||
        atom_tcoord[1] >= hi_grid[1] ||
        atom_tcoord[2] >= hi_grid[2]) {

        elecs[(atom_id * ttl_poses[0]) + pose_id] = INFINITY;
        emaps[(atom_id * ttl_poses[0]) + pose_id] = INFINITY;
        return;
    }
    // 3D Linear Interpolation
    double u = (atom_tcoord[0] - lo_grid[0]) / field_spacing[0];
    double v = (atom_tcoord[1] - lo_grid[1]) / field_spacing[0];
    double w = (atom_tcoord[2] - lo_grid[2]) / field_spacing[0];

    long u0 = (long)u;
    long v0 = (long)v;
    long w0 = (long)w;

    long u1 = u0 + 1;
    long v1 = v0 + 1;
    long w1 = w0 + 1;

    double p0u = u - (double)u0;
    double p0v = v - (double)v0;
    double p0w = w - (double)w0;

    double p1u = (double)u1 - u;
    double p1v = (double)v1 - v;
    double p1w = (double)w1 - w;

    double p000 = p0u * p0v * p0w;
    double p001 = p0u * p0v * p1w;
    double p010 = p0u * p1v * p0w;
    double p011 = p0u * p1v * p1w;
    double p100 = p1u * p0v * p0w;
    double p101 = p1u * p0v * p1w;
    double p110 = p1u * p1v * p0w;
    double p111 = p1u * p1v * p1w;
    
    // AutoDock uses z, y, x axis order for maps
    long num_points1_2 = num_points1[2] * num_points1[1];
    long num_points1_1 = num_points1[2];
    // Atom type
    long atom_type_id = (long)atoms_properties[(atom_id * ttl_atom_properties[0]) +
                                               ATOM_TYPE_IDX];

    // Energy calculation
    double e = 0.0; // Electrostatic
    double d = 0.0; // Desolvation
    double m = 0.0; // Aton type

    long w1v1u1_idx = ttl_maps[0] *
                      ((w1 * num_points1_2) + (v1 * num_points1_1) + u1);
    e += p000 * maps[w1v1u1_idx + electrostatic_lut[0]];
    d += p000 * maps[w1v1u1_idx + desolvation_lut[0]];
    m += p000 * maps[w1v1u1_idx + atom_type_map_lut[atom_type_id]];
    long w1v1u0_idx = ttl_maps[0] *
                      ((w1 * num_points1_2) + (v1 * num_points1_1) + u0);
    e += p001 * maps[w1v1u0_idx + electrostatic_lut[0]];
    d += p001 * maps[w1v1u0_idx + desolvation_lut[0]];
    m += p001 * maps[w1v1u0_idx + atom_type_map_lut[atom_type_id]];
    long w1v0u1_idx = ttl_maps[0] *
                      ((w1 * num_points1_2) + (v0 * num_points1_1) + u1);
    e += p010 * maps[w1v0u1_idx + electrostatic_lut[0]];
    d += p010 * maps[w1v0u1_idx + desolvation_lut[0]];
    m += p010 * maps[w1v0u1_idx + atom_type_map_lut[atom_type_id]];
    long w1v0u0_idx = ttl_maps[0] *
                      ((w1 * num_points1_2) + (v0 * num_points1_1) + u0);
    e += p011 * maps[w1v0u0_idx + electrostatic_lut[0]];
    d += p011 * maps[w1v0u0_idx + desolvation_lut[0]];
    m += p011 * maps[w1v0u0_idx + atom_type_map_lut[atom_type_id]];
    long w0v1u1_idx = ttl_maps[0] *
                      ((w0 * num_points1_2) + (v1 * num_points1_1) + u1);
    e += p100 * maps[w0v1u1_idx + electrostatic_lut[0]];
    d += p100 * maps[w0v1u1_idx + desolvation_lut[0]];
    m += p100 * maps[w0v1u1_idx + atom_type_map_lut[atom_type_id]];
    long w0v1u0_idx = ttl_maps[0] *
                      ((w0 * num_points1_2) + (v1 * num_points1_1) + u0);
    e += p101 * maps[w0v1u0_idx + electrostatic_lut[0]];
    d += p101 * maps[w0v1u0_idx + desolvation_lut[0]];
    m += p101 * maps[w0v1u0_idx + atom_type_map_lut[atom_type_id]];
    long w0v0u1_idx = ttl_maps[0] *
                      ((w0 * num_points1_2) + (v0 * num_points1_1) + u1);
    e += p110 * maps[w0v0u1_idx + electrostatic_lut[0]];
    d += p110 * maps[w0v0u1_idx + desolvation_lut[0]];
    m += p110 * maps[w0v0u1_idx + atom_type_map_lut[atom_type_id]];
    long w0v0u0_idx = ttl_maps[0] *
                      ((w0 * num_points1_2) + (v0 * num_points1_1) + u0);
    e += p111 * maps[w0v0u0_idx + electrostatic_lut[0]];
    d += p111 * maps[w0v0u0_idx + desolvation_lut[0]];
    m += p111 * maps[w0v0u0_idx + atom_type_map_lut[atom_type_id]];

    double atom_charge = atoms_properties[(atom_id * ttl_atom_properties[0]) +
                                          ATOM_CHARGE_IDX];
    double abs_atom_charge = fabs(atom_charge);

    elecs[(atom_id * ttl_poses[0]) + pose_id] = e * atom_charge;
    emaps[(atom_id * ttl_poses[0]) + pose_id] = m + (d * abs_atom_charge);
}

__kernel void calc_total_inter_energy(__global const long *ttl_atoms,
                                      __global const long *ttl_poses,
                                      __global const double *elecs,
                                      __global const double *emaps,

                                      __global double *elec_totals,
                                      __global double *emap_totals)
{
    long thread_id = get_global_id(0);
    // Atom ID starts from index 1
    long e_id = (thread_id / ttl_poses[0]);
    // Pose ID or individual ID
    long pose_id = thread_id % ttl_poses[0];
    // elec totals
    if (e_id == 0) {
        double elec_total = 0.0;
        for (long i = 1; i < ttl_atoms[0] + 1; i++) {
            elec_total += elecs[(i * ttl_poses[0]) + pose_id];
        }
        elec_totals[pose_id] = elec_total;
    }
    // emap totals
    if (e_id == 1) {
        double emap_total = 0.0;
        for (long i = 1; i < ttl_atoms[0] + 1; i++) {
            emap_total += emaps[(i * ttl_poses[0]) + pose_id];
        }
        emap_totals[pose_id] = emap_total;
    }
}

__kernel void calc_intra_energy(__global const long *ttl_poses,
                                __global const double *poses,
                                __global const double *lo_grid,
                                __global const double *hi_grid,

                                __global const long *ttl_non_bond_list,
                                __global const long *ttl_non_bond_properties,
                                __global const double *non_bond_list,

                                __global const long *ttl_atom_types,
                                __global const double *bond_properties,

                                __global const long *calc_inter_elec_e,
                                __global const long *include_1_4_interactions,

                                __global const double *et_inv_r_epsilon,
                                __global const double *et_solvation,
                                __global const double *et_vdw_hb,
                                
                                __global double *e_internals)
{
    long thread_id = get_global_id(0);
    // Pose ID or individual ID
    long pose_id = (thread_id / ttl_non_bond_list[0]);
    // Non-bond list ID
    long nb_id = thread_id % ttl_non_bond_list[0];
    // Get non-bond properties
    long non_bond_start_idx = nb_id * ttl_non_bond_properties[0];
    long atom_id1 = non_bond_list[non_bond_start_idx + ATOM_ID1_IDX];
    long atom_type1 = non_bond_list[non_bond_start_idx + ATOM_TYPE1_IDX];
    long atom_id2 = non_bond_list[non_bond_start_idx + ATOM_ID2_IDX];
    long atom_type2 = non_bond_list[non_bond_start_idx + ATOM_TYPE2_IDX];
    long non_bond_type = non_bond_list[non_bond_start_idx + NON_BOND_TYPE_IDX];
    double desolv = non_bond_list[non_bond_start_idx + DESOLV_IDX];
    double q1q2 = non_bond_list[non_bond_start_idx + Q1Q2_IDX];

    // Atom coordinate
    double3 atom_tcoord1;
    for (long i = 0; i < 3; i++) {
        atom_tcoord1[i] = poses[(atom_id1 * ttl_poses[0] * 3) +
                                (pose_id * 3) + i];
    }
    double3 atom_tcoord2;
    for (long i = 0; i < 3; i++) {
        atom_tcoord2[i] = poses[(atom_id2 * ttl_poses[0] * 3) +
                                (pose_id * 3) + i];
    }
    // Check out of grid
    if (atom_tcoord1[0] <= lo_grid[0] ||
        atom_tcoord1[1] <= lo_grid[1] ||
        atom_tcoord1[2] <= lo_grid[2] ||
        atom_tcoord1[0] >= hi_grid[0] ||
        atom_tcoord1[1] >= hi_grid[1] ||
        atom_tcoord1[2] >= hi_grid[2] ||
        atom_tcoord2[0] <= lo_grid[0] ||
        atom_tcoord2[1] <= lo_grid[1] ||
        atom_tcoord2[2] <= lo_grid[2] ||
        atom_tcoord2[0] >= hi_grid[0] ||
        atom_tcoord2[1] >= hi_grid[1] ||
        atom_tcoord2[2] >= hi_grid[2]) {

        e_internals[(pose_id * ttl_non_bond_list[0]) + nb_id] = INFINITY;
        return;
    }
    // Calculate distance square
    double3 r_tcoord2 = atom_tcoord1 - atom_tcoord2;
    double r2 = (r_tcoord2[0] * r_tcoord2[0]) +
                (r_tcoord2[1] * r_tcoord2[1]) +
                (r_tcoord2[2] * r_tcoord2[2]);
    r2 = max(bond_properties[RMIN_ELEC2_IDX], r2);  // Clamp r2 at RMIN_ELEC2
    long index = (long)(r2 * bond_properties[SQA_DIV_IDX]);
    // Make sure the indexes are not greater than NS_INTL -1 and NS_EL - 1
    // respectively
    long i_ns_intl = min(index, (long)bond_properties[NS_INTL_1_IDX]);
    long i_ns_el = min(index, (long)bond_properties[NS_EL_1_IDX]);

    double e_internal = 0.0;
    if (calc_inter_elec_e[0] == 1) {
        // Calculate Electrostatic Energy
        e_internal += q1q2 * et_inv_r_epsilon[i_ns_el];
    }
    if (r2 < bond_properties[NBC2_IDX]) {
        // Calculate Desolvation Energy
        double e_desolv = desolv * et_solvation[i_ns_intl];
        // Calculate Van der Waals and Hydrogen Bond Energies
        long a1 = atom_type1;
        long a2 = atom_type2;
        if (atom_type1 > atom_type2) {
            a1 = atom_type2;
            a2 = atom_type1;
        }
        long n = ttl_atom_types[0];
        long col = (n * a1) + a2 - ((a1 * (a1 + 1)) / 2);
        long ns_intl = (long)bond_properties[NS_INTL_1_IDX] + 1;
        if (include_1_4_interactions[0] == 1 && non_bond_type == 4) {
            e_internal += bond_properties[SCALE_1_4_INTERACTIONS_IDX] +
                          (et_vdw_hb[(col * ns_intl) + i_ns_intl] + e_desolv);
        } else {
            e_internal += et_vdw_hb[(col * ns_intl) + i_ns_intl] + e_desolv;
        }
    }
    e_internals[(pose_id * ttl_non_bond_list[0]) + nb_id] = e_internal;
}

__kernel void calc_total_intra_energy(__global const long *ttl_poses,
                                      __global const long *ttl_non_bond_list,
                                      __global const double *e_internals,

                                      __global double *e_internal_totals)
{
    // Pose ID or individual ID
    long pose_id = get_global_id(0);

    // Internal energy totals
    double total_e_internal = 0.0;
    for (long i = 0; i < ttl_non_bond_list[0]; i++) {
        total_e_internal += e_internals[(pose_id * ttl_non_bond_list[0]) + i];
    }
    e_internal_totals[pose_id] = total_e_internal;
}

__kernel void calc_total_energy(__global const double *elec_totals,
                                __global const double *emap_totals,
                                __global const double *e_internal_totals,

                                __global double *e_totals)
{
    // Pose ID or individual ID
    long pose_id = get_global_id(0);
    e_totals[pose_id] = elec_totals[pose_id] +
                        emap_totals[pose_id] +
                        e_internal_totals[pose_id];
}

