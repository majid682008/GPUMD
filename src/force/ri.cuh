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

#pragma once
#include "potential.cuh"
#include <stdio.h>
#include <vector>

struct RI_Para {
  double a11, b11, c11, qq11;
  double a12, b12, c12, qq12;
  double a13, b13, c13, qq13;
  double a14, b14, c14, qq14;
  double a15, b15, c15, qq15;
  double a16, b16, c16, qq16;
  double a22, b22, c22, qq22;
  double a23, b23, c23, qq23;
  double a24, b24, c24, qq24;
  double a25, b25, c25, qq25;
  double a26, b26, c26, qq26;
  double a33, b33, c33, qq33;
  double a34, b34, c34, qq34;
  double a35, b35, c35, qq35;
  double a36, b36, c36, qq36;
  double a44, b44, c44, qq44;
  double a45, b45, c45, qq45;
  double a46, b46, c46, qq46;
  double a55, b55, c55, qq55;
  double a56, b56, c56, qq56;
  double a66, b66, c66, qq66;
  double v_rc, dv_rc; // potential and its derivative at the cutoff distance
  double cutoff;
};

class RI : public Potential
{
public:
  RI(FILE*, int);
  virtual ~RI(void);

  virtual void compute(
    const int type_shift,
    const Box& box,
    const Neighbor& neighbor,
    const GPU_Vector<int>& type,
    const GPU_Vector<double>& position,
    GPU_Vector<double>& potential,
    GPU_Vector<double>& force,
    GPU_Vector<double>& virial);

  void initialize_ri(FILE* fid);

protected:
  RI_Para ri_para;
};
