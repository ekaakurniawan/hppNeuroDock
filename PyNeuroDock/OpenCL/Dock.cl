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

#define TRANSLATION_START_IDX 0
#define ROTATION_START_IDX 3
#define TORSION_START_IDX 7

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
        long atom_tcoord_id = branches_rot_seq_buf[(br_i * longest_branch[0]) + \
                                                           br_j];
        if (atom_tcoord_id == 0) {
            barrier(CLK_LOCAL_MEM_FENCE | CLK_GLOBAL_MEM_FENCE);
            continue;
        }
        // Torsion angle
        double tor_angle = individuals[(pose_id * (3 + 4 + ttl_torsions[0])) + \
                                       TORSION_START_IDX + br_i];
        // Torsion axis
        double anchor_tcoord[3];
        double link_tcoord[3];
        double tor_axis[3];
        long anchor_tcoord_id = branches_rot_anchor_buf[br_i];
        long link_tcoord_id =  branches_rot_link_buf[br_i];
        for (long i = 0; i < 3; i++) {
            anchor_tcoord[i] = poses[(anchor_tcoord_id * ttl_poses[0] * 3) + \
                                     (pose_id * 3) + i];
            link_tcoord[i] = poses[(link_tcoord_id * ttl_poses[0] * 3) + \
                                   (pose_id * 3) + i];
            tor_axis[i] = anchor_tcoord[i] - link_tcoord[i];
        }
        // Rotation in quaternion
        double rotation[4];
        double half_tor_angle = tor_angle / 2;
        rotation[0] = cos(half_tor_angle);
        double tor_axis_mag = sqrt(pow(tor_axis[0], 2.0) +
                                   pow(tor_axis[1], 2.0) +
                                   pow(tor_axis[2], 2.0));
        if (tor_axis_mag > 0.0) {
            tor_axis[0] /= tor_axis_mag;
            tor_axis[1] /= tor_axis_mag;
            tor_axis[2] /= tor_axis_mag;
        }
        double s = sin(half_tor_angle);
        rotation[1] = tor_axis[0] * s;
        rotation[2] = tor_axis[1] * s;
        rotation[3] = tor_axis[2] * s;
        // Transform
        double atom_tcoord[3];
        for (long i = 0; i < 3; i++) {
            atom_tcoord[i] = poses[(atom_tcoord_id * ttl_poses[0] * 3) + \
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
    // Pose ID or individual ID
    long pose_id = thread_id / ttl_ligand_atoms[0];
    // Atom ID starts from index 1
    long atom_id = (thread_id % ttl_ligand_atoms[0]) + 1;
    if (atom_id == 0) return;
    // Translation
    double translation[3];
    for (long i = 0; i < 3; i++) {
        translation[i] = individuals[(pose_id * (3 + 4 + ttl_torsions[0])) + \
                                  TRANSLATION_START_IDX + i];
    }
    // Rotation
    double rotation[4];
    for (long i = 0; i < 4; i++) {
        rotation[i] = individuals[(pose_id * (3 + 4 + ttl_torsions[0])) + \
                               ROTATION_START_IDX + i];
    }
    // Atom coordinate
    double atom_tcoord[3];
    for (long i = 0; i < 3; i++) {
        atom_tcoord[i] = poses[(atom_id * ttl_poses[0] * 3) + \
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

