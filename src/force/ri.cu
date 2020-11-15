/*
    Copyright 2017 Zheyong Fan, Ville Vierimaa, Mikko Ervasti, and Ari Harju
    This file is part of GPUMD.
    GPUMD is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    GPUMD is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with GPUMD.  If not, see <http://www.gnu.org/licenses/>.
*/

/*----------------------------------------------------------------------------80
The class dealing with the rigid-ion potential.

Reference for the method of evaluating the Coulomb force in the rigid-ion
potential:

[1] C. J. Fennell and J. D. Gezelter. Is the Ewald summation still necessary?
Pairwise alternatives to the accepted standard for long-range electrostatics,
J. Chem. Phys. 124, 234104 (2006).
------------------------------------------------------------------------------*/

#include "ri.cuh"
#include "utilities/common.cuh"
#include "utilities/error.cuh"

#define BLOCK_SIZE_FORCE 128
#define RI_ALPHA 0.2
#define RI_ALPHA_SQ 0.04
#define RI_PI_FACTOR 0.225675833419103 // ALPHA * 2 / SQRT(PI)

RI::RI(FILE* fid, int num_types)
{
  printf("Use the rigid-ion potential.\n");

  int N = num_types*(num_types+1)/2;

  fscanf(fid, "%lf", &ri_para.cutoff);

  double q[num_types];
  for (int i = 0; i < num_types; i++) {
    fscanf(fid, "%lf", &q[i]);
  }

  double x[N][3];
  for (int i = 0; i < N; i++) {
    int count = fscanf(fid, "%lf%lf%lf", &x[i][0], &x[i][1], &x[i][2]);
    PRINT_SCANF_ERROR(count, 3, "Reading error for rigid-ion potential.");
  }

  if (num_types == 2){
    ri_para.a11 = x[0][0];
    ri_para.b11 = 1.0 / x[0][1];
    ri_para.c11 = x[0][2];
    ri_para.a12 = x[1][0];
    ri_para.b12 = 1.0 / x[1][1];
    ri_para.c12 = x[1][2];
    ri_para.a22 = x[2][0];
    ri_para.b22 = 1.0 / x[2][1];
    ri_para.c22 = x[2][2];
    ri_para.qq11 = q[0] * q[0] * K_C;
    ri_para.qq12 = q[0] * q[1] * K_C;
    ri_para.qq22 = q[1] * q[1] * K_C;
  }
  else if (num_types == 3){
    ri_para.a11 = x[0][0];
    ri_para.b11 = 1.0 / x[0][1];
    ri_para.c11 = x[0][2];
    ri_para.a12 = x[1][0];
    ri_para.b12 = 1.0 / x[1][1];
    ri_para.c12 = x[1][2];
    ri_para.a13 = x[2][0];
    ri_para.b13 = 1.0 / x[2][1];
    ri_para.c13 = x[2][2];
    ri_para.a22 = x[3][0];
    ri_para.b22 = 1.0 / x[3][1];
    ri_para.c22 = x[3][2];
    ri_para.a23 = x[4][0];
    ri_para.b23 = 1.0 / x[4][1];
    ri_para.c23 = x[4][2];
    ri_para.a33 = x[5][0];
    ri_para.b33 = 1.0 / x[5][1];
    ri_para.c33 = x[5][2];
    ri_para.qq11 = q[0] * q[0] * K_C;
    ri_para.qq12 = q[0] * q[1] * K_C;
    ri_para.qq13 = q[0] * q[2] * K_C;
    ri_para.qq22 = q[1] * q[1] * K_C;
    ri_para.qq23 = q[1] * q[2] * K_C;
    ri_para.qq33 = q[2] * q[2] * K_C;
  }
  else if (num_types == 4){
    ri_para.a11 = x[0][0];
    ri_para.b11 = 1.0 / x[0][1];
    ri_para.c11 = x[0][2];
    ri_para.a12 = x[1][0];
    ri_para.b12 = 1.0 / x[1][1];
    ri_para.c12 = x[1][2];
    ri_para.a13 = x[2][0];
    ri_para.b13 = 1.0 / x[2][1];
    ri_para.c13 = x[2][2];
    ri_para.a14 = x[3][0];
    ri_para.b14 = 1.0 / x[3][1];
    ri_para.c14 = x[3][2];
    ri_para.a22 = x[4][0];
    ri_para.b22 = 1.0 / x[4][1];
    ri_para.c22 = x[4][2];
    ri_para.a23 = x[5][0];
    ri_para.b23 = 1.0 / x[5][1];
    ri_para.c23 = x[5][2];
    ri_para.a24 = x[6][0];
    ri_para.b24 = 1.0 / x[6][1];
    ri_para.c24 = x[6][2];
    ri_para.a33 = x[7][0];
    ri_para.b33 = 1.0 / x[7][1];
    ri_para.c33 = x[7][2];
    ri_para.a34 = x[8][0];
    ri_para.b34 = 1.0 / x[8][1];
    ri_para.c34 = x[8][2];
    ri_para.a44 = x[9][0];
    ri_para.b44 = 1.0 / x[9][1];
    ri_para.c44 = x[9][2];
    ri_para.qq11 = q[0] * q[0] * K_C;
    ri_para.qq12 = q[0] * q[1] * K_C;
    ri_para.qq13 = q[0] * q[2] * K_C;
    ri_para.qq14 = q[0] * q[3] * K_C;
    ri_para.qq22 = q[1] * q[1] * K_C;
    ri_para.qq23 = q[1] * q[2] * K_C;
    ri_para.qq24 = q[1] * q[3] * K_C;
    ri_para.qq33 = q[2] * q[2] * K_C;
    ri_para.qq34 = q[2] * q[3] * K_C;
    ri_para.qq44 = q[3] * q[3] * K_C;
  }
  else if (num_types == 5){
    ri_para.a11 = x[0][0];
    ri_para.b11 = 1.0 / x[0][1];
    ri_para.c11 = x[0][2];
    ri_para.a12 = x[1][0];
    ri_para.b12 = 1.0 / x[1][1];
    ri_para.c12 = x[1][2];
    ri_para.a13 = x[2][0];
    ri_para.b13 = 1.0 / x[2][1];
    ri_para.c13 = x[2][2];
    ri_para.a14 = x[3][0];
    ri_para.b14 = 1.0 / x[3][1];
    ri_para.c14 = x[3][2];
    ri_para.a15 = x[4][0];
    ri_para.b15 = 1.0 / x[4][1];
    ri_para.c15 = x[4][2];
    ri_para.a22 = x[5][0];
    ri_para.b22 = 1.0 / x[5][1];
    ri_para.c22 = x[5][2];
    ri_para.a23 = x[6][0];
    ri_para.b23 = 1.0 / x[6][1];
    ri_para.c23 = x[6][2];
    ri_para.a24 = x[7][0];
    ri_para.b24 = 1.0 / x[7][1];
    ri_para.c24 = x[7][2];
    ri_para.a25 = x[8][0];
    ri_para.b25 = 1.0 / x[8][1];
    ri_para.c25 = x[8][2];
    ri_para.a33 = x[9][0];
    ri_para.b33 = 1.0 / x[9][1];
    ri_para.c33 = x[9][2];
    ri_para.a34 = x[10][0];
    ri_para.b34 = 1.0 / x[10][1];
    ri_para.c34 = x[10][2];
    ri_para.a35 = x[11][0];
    ri_para.b35 = 1.0 / x[11][1];
    ri_para.c35 = x[11][2];
    ri_para.a44 = x[12][0];
    ri_para.b44 = 1.0 / x[12][1];
    ri_para.c44 = x[12][2];
    ri_para.a45 = x[13][0];
    ri_para.b45 = 1.0 / x[13][1];
    ri_para.c45 = x[13][2];
    ri_para.a55 = x[14][0];
    ri_para.b55 = 1.0 / x[14][1];
    ri_para.c55 = x[14][2];
    ri_para.qq11 = q[0] * q[0] * K_C;
    ri_para.qq12 = q[0] * q[1] * K_C;
    ri_para.qq13 = q[0] * q[2] * K_C;
    ri_para.qq14 = q[0] * q[3] * K_C;
    ri_para.qq15 = q[0] * q[4] * K_C;
    ri_para.qq22 = q[1] * q[1] * K_C;
    ri_para.qq23 = q[1] * q[2] * K_C;
    ri_para.qq24 = q[1] * q[3] * K_C;
    ri_para.qq25 = q[1] * q[4] * K_C;
    ri_para.qq33 = q[2] * q[2] * K_C;
    ri_para.qq34 = q[2] * q[3] * K_C;
    ri_para.qq35 = q[2] * q[4] * K_C;
    ri_para.qq44 = q[3] * q[3] * K_C;
    ri_para.qq45 = q[3] * q[4] * K_C;
    ri_para.qq55 = q[4] * q[4] * K_C;
  }
  else if (num_types == 6){
    ri_para.a11 = x[0][0];
    ri_para.b11 = 1.0 / x[0][1];
    ri_para.c11 = x[0][2];
    ri_para.a12 = x[1][0];
    ri_para.b12 = 1.0 / x[1][1];
    ri_para.c12 = x[1][2];
    ri_para.a13 = x[2][0];
    ri_para.b13 = 1.0 / x[2][1];
    ri_para.c13 = x[2][2];
    ri_para.a14 = x[3][0];
    ri_para.b14 = 1.0 / x[3][1];
    ri_para.c14 = x[3][2];
    ri_para.a15 = x[4][0];
    ri_para.b15 = 1.0 / x[4][1];
    ri_para.c15 = x[4][2];
    ri_para.a16 = x[5][0];
    ri_para.b16 = 1.0 / x[5][1];
    ri_para.c16 = x[5][2];
    ri_para.a22 = x[6][0];
    ri_para.b22 = 1.0 / x[6][1];
    ri_para.c22 = x[6][2];
    ri_para.a23 = x[7][0];
    ri_para.b23 = 1.0 / x[7][1];
    ri_para.c23 = x[7][2];
    ri_para.a24 = x[8][0];
    ri_para.b24 = 1.0 / x[8][1];
    ri_para.c24 = x[8][2];
    ri_para.a25 = x[9][0];
    ri_para.b25 = 1.0 / x[9][1];
    ri_para.c25 = x[9][2];
    ri_para.a26 = x[10][0];
    ri_para.b26 = 1.0 / x[10][1];
    ri_para.c26 = x[10][2];
    ri_para.a33 = x[11][0];
    ri_para.b33 = 1.0 / x[11][1];
    ri_para.c33 = x[11][2];
    ri_para.a34 = x[12][0];
    ri_para.b34 = 1.0 / x[12][1];
    ri_para.c34 = x[12][2];
    ri_para.a35 = x[13][0];
    ri_para.b35 = 1.0 / x[13][1];
    ri_para.c35 = x[13][2];
    ri_para.a36 = x[14][0];
    ri_para.b36 = 1.0 / x[14][1];
    ri_para.c36 = x[14][2];
    ri_para.a44 = x[15][0];
    ri_para.b44 = 1.0 / x[15][1];
    ri_para.c44 = x[15][2];
    ri_para.a45 = x[16][0];
    ri_para.b45 = 1.0 / x[16][1];
    ri_para.c45 = x[16][2];
    ri_para.a46 = x[17][0];
    ri_para.b46 = 1.0 / x[17][1];
    ri_para.c46 = x[17][2];
    ri_para.a55 = x[18][0];
    ri_para.b55 = 1.0 / x[18][1];
    ri_para.c55 = x[18][2];
    ri_para.a56 = x[19][0];
    ri_para.b56 = 1.0 / x[19][1];
    ri_para.c56 = x[19][2];
    ri_para.a66 = x[20][0];
    ri_para.b66 = 1.0 / x[20][1];
    ri_para.c66 = x[20][2];
    ri_para.qq11 = q[0] * q[0] * K_C;
    ri_para.qq12 = q[0] * q[1] * K_C;
    ri_para.qq13 = q[0] * q[2] * K_C;
    ri_para.qq14 = q[0] * q[3] * K_C;
    ri_para.qq15 = q[0] * q[4] * K_C;
    ri_para.qq16 = q[0] * q[5] * K_C;
    ri_para.qq22 = q[1] * q[1] * K_C;
    ri_para.qq23 = q[1] * q[2] * K_C;
    ri_para.qq24 = q[1] * q[3] * K_C;
    ri_para.qq25 = q[1] * q[4] * K_C;
    ri_para.qq26 = q[1] * q[5] * K_C;
    ri_para.qq33 = q[2] * q[2] * K_C;
    ri_para.qq34 = q[2] * q[3] * K_C;
    ri_para.qq35 = q[2] * q[4] * K_C;
    ri_para.qq36 = q[2] * q[5] * K_C;
    ri_para.qq44 = q[3] * q[3] * K_C;
    ri_para.qq45 = q[3] * q[4] * K_C;
    ri_para.qq46 = q[3] * q[5] * K_C;
    ri_para.qq55 = q[4] * q[4] * K_C;
    ri_para.qq56 = q[4] * q[5] * K_C;
    ri_para.qq55 = q[5] * q[5] * K_C;
  }

  rc = ri_para.cutoff; // force cutoff

  ri_para.v_rc = erfc(RI_ALPHA * rc) / rc;
  ri_para.dv_rc = -erfc(RI_ALPHA * rc) / (rc * rc);
  ri_para.dv_rc -= RI_PI_FACTOR * exp(-RI_ALPHA_SQ * rc * rc) / rc;
}

RI::~RI(void)
{
  // nothing
}

// get U_ij and (d U_ij / d r_ij) / r_ij (the RI potential)
static __device__ void
find_p2_and_f2(int type1, int type2, RI_Para ri, double d12sq, double& p2, double& f2)
{
  double a, b, c, qq;

  if (type1 == 0 && type2 == 0) {
    a = ri.a11;
    b = ri.b11;
    c = ri.c11;
    qq = ri.qq11;
  } else if  ((type1 == 0 && type2 == 1) || (type1 == 1 && type2 == 0))  {
    a = ri.a12;
    b = ri.b12;
    c = ri.c12;
    qq = ri.qq12;
  } else if (type1 == 1 && type2 == 1) {
    a = ri.a22;
    b = ri.b22;
    c = ri.c22;
    qq = ri.qq22;
  } else if  ((type1 == 0 && type2 == 2) || (type1 == 2 && type2 == 0))  {
    a = ri.a13;
    b = ri.b13;
    c = ri.c13;
    qq = ri.qq13;
  } else if ((type1 == 1 && type2 == 2) || (type1 == 2 && type2 == 1))  {
    a = ri.a23;
    b = ri.b23;
    c = ri.c23;
    qq = ri.qq23;
  } else if (type1 == 2 && type2 == 2) {
    a = ri.a33;
    b = ri.b33;
    c = ri.c33;
    qq = ri.qq33;
  } else if  ((type1 == 0 && type2 == 3) || (type1 == 3 && type2 == 0))  {
    a = ri.a14;
    b = ri.b14;
    c = ri.c14;
    qq = ri.qq14;
  } else if ((type1 == 1 && type2 == 3) || (type1 == 3 && type2 == 1))  {
    a = ri.a24;
    b = ri.b24;
    c = ri.c24;
    qq = ri.qq24;
  } else if ((type1 == 2 && type2 == 3) || (type1 == 3 && type2 == 2))  {
    a = ri.a34;
    b = ri.b34;
    c = ri.c34;
    qq = ri.qq34;
  } else if (type1 == 3 && type2 == 3) {
    a = ri.a44;
    b = ri.b44;
    c = ri.c44;
    qq = ri.qq44;
  } else if  ((type1 == 0 && type2 == 4) || (type1 == 4 && type2 == 0))  {
    a = ri.a15;
    b = ri.b15;
    c = ri.c15;
    qq = ri.qq15;
  } else if ((type1 == 1 && type2 == 4) || (type1 == 4 && type2 == 1))  {
    a = ri.a25;
    b = ri.b25;
    c = ri.c25;
    qq = ri.qq25;
  } else if ((type1 == 2 && type2 == 4) || (type1 == 4 && type2 == 2))  {
    a = ri.a35;
    b = ri.b35;
    c = ri.c35;
    qq = ri.qq35;
  } else if ((type1 == 3 && type2 == 4) || (type1 == 4 && type2 == 3))  {
    a = ri.a45;
    b = ri.b45;
    c = ri.c45;
    qq = ri.qq45;
  } else if (type1 == 4 && type2 == 4) {
    a = ri.a55;
    b = ri.b55;
    c = ri.c55;
    qq = ri.qq55;
  } else if  ((type1 == 0 && type2 == 5) || (type1 == 5 && type2 == 0))  {
    a = ri.a16;
    b = ri.b16;
    c = ri.c16;
    qq = ri.qq16;
  } else if ((type1 == 1 && type2 == 5) || (type1 == 5 && type2 == 1))  {
    a = ri.a26;
    b = ri.b26;
    c = ri.c26;
    qq = ri.qq26;
  } else if ((type1 == 2 && type2 == 5) || (type1 == 5 && type2 == 2))  {
    a = ri.a36;
    b = ri.b36;
    c = ri.c36;
    qq = ri.qq36;
  } else if ((type1 == 3 && type2 == 5) || (type1 == 5 && type2 == 3))  {
    a = ri.a46;
    b = ri.b46;
    c = ri.c46;
    qq = ri.qq46;
  } else if ((type1 == 4 && type2 == 5) || (type1 == 5 && type2 == 4))  {
    a = ri.a56;
    b = ri.b56;
    c = ri.c56;
    qq = ri.qq56;
  } else if (type1 == 5 && type2 == 5) {
    a = ri.a66;
    b = ri.b66;
    c = ri.c66;
    qq = ri.qq66;
  }

  double d12 = sqrt(d12sq);
  double d12inv = 1.0 / d12;
  double d12inv3 = d12inv * d12inv * d12inv;
  double exponential = exp(-d12 * b); // b = 1/rho
  double erfc_r = erfc(RI_ALPHA * d12) * d12inv;
  p2 = a * exponential - c * d12inv3 * d12inv3;
  p2 += qq * (erfc_r - ri.v_rc - ri.dv_rc * (d12 - ri.cutoff));
  f2 = 6.0 * c * (d12inv3 * d12inv3 * d12inv) - a * exponential * b;
  f2 -= qq * (erfc_r * d12inv + RI_PI_FACTOR * d12inv * exp(-RI_ALPHA_SQ * d12sq) + ri.dv_rc);
  f2 *= d12inv;
}

// force evaluation kernel
static __global__ void gpu_find_force(
  RI_Para ri,
  const int number_of_particles,
  const int N1,
  const int N2,
  const Box box,
  const int* g_neighbor_number,
  const int* g_neighbor_list,
  const int* g_type,
  const int shift,
  const double* __restrict__ g_x,
  const double* __restrict__ g_y,
  const double* __restrict__ g_z,
  double* g_fx,
  double* g_fy,
  double* g_fz,
  double* g_virial,
  double* g_potential)
{
  int n1 = blockIdx.x * blockDim.x + threadIdx.x + N1; // particle index
  double s_fx = 0.0;                                   // force_x
  double s_fy = 0.0;                                   // force_y
  double s_fz = 0.0;                                   // force_z
  double s_pe = 0.0;                                   // potential energy
  double s_sxx = 0.0;                                  // virial_stress_xx
  double s_sxy = 0.0;                                  // virial_stress_xy
  double s_sxz = 0.0;                                  // virial_stress_xz
  double s_syx = 0.0;                                  // virial_stress_yx
  double s_syy = 0.0;                                  // virial_stress_yy
  double s_syz = 0.0;                                  // virial_stress_yz
  double s_szx = 0.0;                                  // virial_stress_zx
  double s_szy = 0.0;                                  // virial_stress_zy
  double s_szz = 0.0;                                  // virial_stress_zz

  if (n1 >= N1 && n1 < N2) {
    int neighbor_number = g_neighbor_number[n1];
    int type1 = g_type[n1] - shift;
    double x1 = g_x[n1];
    double y1 = g_y[n1];
    double z1 = g_z[n1];

    for (int i1 = 0; i1 < neighbor_number; ++i1) {
      int n2 = g_neighbor_list[n1 + number_of_particles * i1];
      int type2 = g_type[n2] - shift;

      double x12 = g_x[n2] - x1;
      double y12 = g_y[n2] - y1;
      double z12 = g_z[n2] - z1;
      apply_mic(box, x12, y12, z12);
      double d12sq = x12 * x12 + y12 * y12 + z12 * z12;

      double p2, f2;

      // RI
      if (d12sq >= ri.cutoff * ri.cutoff) {
        continue;
      }
      find_p2_and_f2(type1, type2, ri, d12sq, p2, f2);

      // treat two-body potential in the same way as many-body potential
      double f12x = f2 * x12 * 0.5;
      double f12y = f2 * y12 * 0.5;
      double f12z = f2 * z12 * 0.5;
      double f21x = -f12x;
      double f21y = -f12y;
      double f21z = -f12z;

      // accumulate force
      s_fx += f12x - f21x;
      s_fy += f12y - f21y;
      s_fz += f12z - f21z;

      // accumulate potential energy and virial
      s_pe += p2 * 0.5; // two-body potential
      s_sxx += x12 * f21x;
      s_sxy += x12 * f21y;
      s_sxz += x12 * f21z;
      s_syx += y12 * f21x;
      s_syy += y12 * f21y;
      s_syz += y12 * f21z;
      s_szx += z12 * f21x;
      s_szy += z12 * f21y;
      s_szz += z12 * f21z;
    }

    // save force
    g_fx[n1] += s_fx;
    g_fy[n1] += s_fy;
    g_fz[n1] += s_fz;

    // save virial
    // xx xy xz    0 3 4
    // yx yy yz    6 1 5
    // zx zy zz    7 8 2
    g_virial[n1 + 0 * number_of_particles] += s_sxx;
    g_virial[n1 + 1 * number_of_particles] += s_syy;
    g_virial[n1 + 2 * number_of_particles] += s_szz;
    g_virial[n1 + 3 * number_of_particles] += s_sxy;
    g_virial[n1 + 4 * number_of_particles] += s_sxz;
    g_virial[n1 + 5 * number_of_particles] += s_syz;
    g_virial[n1 + 6 * number_of_particles] += s_syx;
    g_virial[n1 + 7 * number_of_particles] += s_szx;
    g_virial[n1 + 8 * number_of_particles] += s_szy;

    // save potential
    g_potential[n1] += s_pe;
  }
}

// Find force and related quantities for pair potentials (A wrapper)
void RI::compute(
  const int type_shift,
  const Box& box,
  const Neighbor& neighbor,
  const GPU_Vector<int>& type,
  const GPU_Vector<double>& position_per_atom,
  GPU_Vector<double>& potential_per_atom,
  GPU_Vector<double>& force_per_atom,
  GPU_Vector<double>& virial_per_atom)
{
  const int number_of_atoms = type.size();
  int grid_size = (N2 - N1 - 1) / BLOCK_SIZE_FORCE + 1;

  gpu_find_force<<<grid_size, BLOCK_SIZE_FORCE>>>(
    ri_para, number_of_atoms, N1, N2, box, neighbor.NN_local.data(), neighbor.NL_local.data(),
    type.data(), type_shift, position_per_atom.data(), position_per_atom.data() + number_of_atoms,
    position_per_atom.data() + number_of_atoms * 2, force_per_atom.data(),
    force_per_atom.data() + number_of_atoms, force_per_atom.data() + 2 * number_of_atoms,
    virial_per_atom.data(), potential_per_atom.data());
  CUDA_CHECK_KERNEL
}
