/*
 * Human_in_Loop.c
 *
 * Code generation for model "Human_in_Loop".
 *
 * Model version              : 1.987
 * Simulink Coder version : 8.13 (R2017b) 24-Jul-2017
 * C source code generated on : Wed May  8 14:27:10 2019
 *
 * Target selection: rti1202.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Custom Processor->Custom
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#include "Human_in_Loop_trc_ptr.h"
#include "Human_in_Loop.h"
#include "Human_in_Loop_private.h"

/* Block signals (auto storage) */
B_Human_in_Loop_T Human_in_Loop_B;

/* Block states (auto storage) */
DW_Human_in_Loop_T Human_in_Loop_DW;

/* Real-time model */
RT_MODEL_Human_in_Loop_T Human_in_Loop_M_;
RT_MODEL_Human_in_Loop_T *const Human_in_Loop_M = &Human_in_Loop_M_;

/* Forward declaration for local functions */
static void Human_in_Loop_mldivide(const real_T A[16], real_T B[4]);
static void Human_in_Loop_power(const real_T a_data[], const int32_T a_size[2],
  real_T y_data[], int32_T y_size[2]);
static void Human_in_Loop_power_a(const real_T a_data[], const int32_T a_size[2],
  real_T y_data[], int32_T y_size[2]);

/* Forward declaration for local functions */
static real_T Human_in_Loop_xnrm2(const real_T x[30], int32_T ix0);
static real_T Human_in_Loop_xnrm2_o(int32_T n, const real_T x[30], int32_T ix0);
static void Human_in_Loop_xgeqp3(real_T A[30], real_T tau[2], int32_T jpvt[2]);

/*
 * Output and update for atomic system:
 *    '<S19>/Mux'
 *    '<S21>/Mux'
 */
void Human_in_Loop_Mux(real_T rtu_x1, real_T rtu_x2, B_Mux_Human_in_Loop_T
  *localB)
{
  /* MATLAB Function 'Sensor Data/Encoder module/Mux': '<S31>:1' */
  /* '<S31>:1:3' */
  localB->x[0] = rtu_x1;
  localB->x[1] = rtu_x2;
}

real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T tmp;
  real_T tmp_0;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = (rtNaN);
  } else {
    tmp = fabs(u0);
    tmp_0 = fabs(u1);
    if (rtIsInf(u1)) {
      if (tmp == 1.0) {
        y = 1.0;
      } else if (tmp > 1.0) {
        if (u1 > 0.0) {
          y = (rtInf);
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = (rtInf);
      }
    } else if (tmp_0 == 0.0) {
      y = 1.0;
    } else if (tmp_0 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = (rtNaN);
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S1>/Torque track' */
static void Human_in_Loop_mldivide(const real_T A[16], real_T B[4])
{
  real_T b_A[16];
  int8_T ipiv[4];
  int32_T j;
  int32_T ix;
  real_T smax;
  real_T s;
  int32_T k;
  int32_T iy;
  int32_T c_ix;
  int32_T d;
  int32_T ijA;
  int32_T kAcol;
  memcpy(&b_A[0], &A[0], sizeof(real_T) << 4U);
  ipiv[0] = 1;
  ipiv[1] = 2;
  ipiv[2] = 3;
  for (j = 0; j < 3; j++) {
    kAcol = j * 5;
    iy = 0;
    ix = kAcol;
    smax = fabs(b_A[kAcol]);
    for (k = 2; k <= 4 - j; k++) {
      ix++;
      s = fabs(b_A[ix]);
      if (s > smax) {
        iy = k - 1;
        smax = s;
      }
    }

    if (b_A[kAcol + iy] != 0.0) {
      if (iy != 0) {
        ipiv[j] = (int8_T)((j + iy) + 1);
        ix = j;
        iy += j;
        smax = b_A[ix];
        b_A[ix] = b_A[iy];
        b_A[iy] = smax;
        ix += 4;
        iy += 4;
        smax = b_A[ix];
        b_A[ix] = b_A[iy];
        b_A[iy] = smax;
        ix += 4;
        iy += 4;
        smax = b_A[ix];
        b_A[ix] = b_A[iy];
        b_A[iy] = smax;
        ix += 4;
        iy += 4;
        smax = b_A[ix];
        b_A[ix] = b_A[iy];
        b_A[iy] = smax;
      }

      iy = (kAcol - j) + 4;
      for (ix = kAcol + 1; ix + 1 <= iy; ix++) {
        b_A[ix] /= b_A[kAcol];
      }
    }

    iy = kAcol;
    ix = kAcol + 4;
    for (k = 1; k <= 3 - j; k++) {
      smax = b_A[ix];
      if (b_A[ix] != 0.0) {
        c_ix = kAcol + 1;
        d = (iy - j) + 8;
        for (ijA = 5 + iy; ijA + 1 <= d; ijA++) {
          b_A[ijA] += b_A[c_ix] * -smax;
          c_ix++;
        }
      }

      ix += 4;
      iy += 4;
    }
  }

  if (ipiv[0] != 1) {
    smax = B[0];
    B[0] = B[ipiv[0] - 1];
    B[ipiv[0] - 1] = smax;
  }

  if (ipiv[1] != 2) {
    smax = B[1];
    B[1] = B[ipiv[1] - 1];
    B[ipiv[1] - 1] = smax;
  }

  if (ipiv[2] != 3) {
    smax = B[2];
    B[2] = B[ipiv[2] - 1];
    B[ipiv[2] - 1] = smax;
  }

  if (B[0] != 0.0) {
    for (iy = 1; iy + 1 < 5; iy++) {
      B[iy] -= B[0] * b_A[iy];
    }
  }

  if (B[1] != 0.0) {
    for (iy = 2; iy + 1 < 5; iy++) {
      B[iy] -= b_A[iy + 4] * B[1];
    }
  }

  if (B[2] != 0.0) {
    for (iy = 3; iy + 1 < 5; iy++) {
      B[iy] -= b_A[iy + 8] * B[2];
    }
  }

  if (B[3] != 0.0) {
    B[3] /= b_A[15];
    for (iy = 0; iy + 1 < 4; iy++) {
      B[iy] -= b_A[iy + 12] * B[3];
    }
  }

  if (B[2] != 0.0) {
    B[2] /= b_A[10];
    for (iy = 0; iy + 1 < 3; iy++) {
      B[iy] -= b_A[iy + 8] * B[2];
    }
  }

  if (B[1] != 0.0) {
    B[1] /= b_A[5];
    for (iy = 0; iy + 1 < 2; iy++) {
      B[iy] -= b_A[iy + 4] * B[1];
    }
  }

  if (B[0] != 0.0) {
    B[0] /= b_A[0];
  }
}

/* Function for MATLAB Function: '<S1>/Torque track' */
static void Human_in_Loop_power(const real_T a_data[], const int32_T a_size[2],
  real_T y_data[], int32_T y_size[2])
{
  int32_T loop_ub;
  int32_T z1_size_idx_1;
  y_size[1] = (int16_T)a_size[1];
  z1_size_idx_1 = y_size[1];
  loop_ub = y_size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&Human_in_Loop_B.z1_data[0], &y_data[0], loop_ub * sizeof(real_T));
  }

  for (loop_ub = 0; loop_ub + 1 <= y_size[1]; loop_ub++) {
    Human_in_Loop_B.z1_data[loop_ub] = a_data[loop_ub] * a_data[loop_ub];
  }

  y_size[0] = 1;
  y_size[1] = z1_size_idx_1;
  if (0 <= z1_size_idx_1 - 1) {
    memcpy(&y_data[0], &Human_in_Loop_B.z1_data[0], z1_size_idx_1 * sizeof
           (real_T));
  }
}

/* Function for MATLAB Function: '<S1>/Torque track' */
static void Human_in_Loop_power_a(const real_T a_data[], const int32_T a_size[2],
  real_T y_data[], int32_T y_size[2])
{
  int32_T loop_ub;
  int32_T z1_size_idx_1;
  y_size[1] = (int16_T)a_size[1];
  z1_size_idx_1 = y_size[1];
  loop_ub = y_size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&Human_in_Loop_B.z1_data_b[0], &y_data[0], loop_ub * sizeof(real_T));
  }

  for (loop_ub = 0; loop_ub + 1 <= y_size[1]; loop_ub++) {
    Human_in_Loop_B.z1_data_b[loop_ub] = rt_powd_snf(a_data[loop_ub], 3.0);
  }

  y_size[0] = 1;
  y_size[1] = z1_size_idx_1;
  if (0 <= z1_size_idx_1 - 1) {
    memcpy(&y_data[0], &Human_in_Loop_B.z1_data_b[0], z1_size_idx_1 * sizeof
           (real_T));
  }
}

/* System initialize for function-call system: '<Root>/Control Module' */
void Human_in_Loo_ControlModule_Init(void)
{
  /* SystemInitialize for MATLAB Function: '<S1>/Torque track' */
  Human_in_Loop_DW.last_footstate = 0.0;

  /* SystemInitialize for MATLAB Function: '<S1>/Controller' */
  Human_in_Loop_DW.calib_state = 0.0;
}

/* System reset for function-call system: '<Root>/Control Module' */
void Human_in_Lo_ControlModule_Reset(void)
{
  /* SystemReset for MATLAB Function: '<S1>/Torque track' */
  Human_in_Loop_DW.last_footstate = 0.0;

  /* SystemReset for MATLAB Function: '<S1>/Controller' */
  Human_in_Loop_DW.calib_state = 0.0;
}

/* Output and update for function-call system: '<Root>/Control Module' */
void Human_in_Loop_ControlModule(void)
{
  real_T mode;
  real_T footstate;
  real_T peak_torque;
  real_T torque_measure;
  real_T troque_delta;
  real_T stride_index;
  real_T index_peak;
  real_T index_rise;
  real_T index_fall;
  real_T parm1[4];
  real_T b;
  int32_T d;
  int32_T f;
  int32_T h;
  int32_T y;
  int32_T b_n;
  int32_T br;
  int32_T ar;
  int32_T ia;
  real_T tmp[16];
  int16_T tmp_data[750];
  int32_T loop_ub;
  int32_T tmp_size[2];
  int32_T tmp_size_0[2];
  int32_T tmp_size_1[2];
  int32_T tmp_size_2[2];
  int32_T tmp_size_3[2];
  int32_T tmp_size_4[2];
  int32_T tmp_size_idx_1;
  int16_T tmp_0;

  /* MATLAB Function: '<S1>/Torque track' */
  /* MATLAB Function 'Control Module/Torque track': '<S9>:1' */
  /* '<S9>:1:42' */
  /* '<S9>:1:19' */
  mode = Human_in_Loop_B.RT1_n[0];

  /* '<S9>:1:20' */
  footstate = Human_in_Loop_B.RT1_n[1];

  /* '<S9>:1:21' */
  /* '<S9>:1:22' */
  /* '<S9>:1:24' */
  peak_torque = Human_in_Loop_B.RT2_h[0];

  /* '<S9>:1:25' */
  /* '<S9>:1:26' */
  /* '<S9>:1:27' */
  /* '<S9>:1:29' */
  torque_measure = Human_in_Loop_B.RT1[0];

  /* '<S9>:1:30' */
  troque_delta = Human_in_Loop_B.RT1[1];

  /* '<S9>:1:31' */
  stride_index = Human_in_Loop_B.RT1_n[3] * 500.0 + 1.0;
  if (stride_index > 750.0) {
    /* '<S9>:1:32' */
    /* '<S9>:1:33' */
    stride_index = 750.0;
  }

  if ((Human_in_Loop_DW.last_footstate == 0.0) && (Human_in_Loop_B.RT1_n[1] ==
       1.0) && ((Human_in_Loop_B.RT1_n[0] == 2.0) || (Human_in_Loop_B.RT1_n[0] ==
        1.0))) {
    /* '<S9>:1:37' */
    /* '<S9>:1:39' */
    /* '<S9>:1:40' */
    memset(&Human_in_Loop_B.torque_track[0], 0, 750U * sizeof(real_T));
    memset(&Human_in_Loop_B.torque_delta_track[0], 0, 750U * sizeof(real_T));

    /* '<S9>:1:41' */
    /* '<S9>:1:42' */
    /* '<S9>:1:43' */
    index_peak = floor(Human_in_Loop_B.RT2_h[2] / 100.0 * Human_in_Loop_B.RT1_n
                       [2] * 500.0);

    /* '<S9>:1:44' */
    index_rise = index_peak - floor(Human_in_Loop_B.RT2_h[1] / 100.0 *
      Human_in_Loop_B.RT1_n[2] * 500.0);

    /* '<S9>:1:45' */
    index_fall = floor(Human_in_Loop_B.RT2_h[3] / 100.0 * Human_in_Loop_B.RT1_n
                       [2] * 500.0) + index_peak;

    /* '<S9>:1:48' */
    /* '<S9>:1:52' */
    /* '<S9>:1:53' */
    parm1[0] = 0.0;
    parm1[1] = Human_in_Loop_B.RT2_h[0];
    parm1[2] = 0.0;
    parm1[3] = 0.0;
    tmp[0] = 1.0;
    tmp[4] = index_rise;
    tmp[8] = index_rise * index_rise;
    tmp[12] = rt_powd_snf(index_rise, 3.0);
    tmp[1] = 1.0;
    tmp[5] = index_peak;
    tmp[9] = index_peak * index_peak;
    tmp[13] = rt_powd_snf(index_peak, 3.0);
    tmp[2] = 0.0;
    tmp[6] = 1.0;
    tmp[10] = 2.0 * index_rise;
    tmp[14] = index_rise * index_rise * 3.0;
    tmp[3] = 0.0;
    tmp[7] = 1.0;
    tmp[11] = 2.0 * index_peak;
    tmp[15] = index_peak * index_peak * 3.0;
    Human_in_Loop_mldivide(tmp, parm1);
    b = index_peak - index_rise;
    if (1.0 > b) {
      b_n = 0;
    } else {
      b_n = (int32_T)b;
    }

    if (index_rise > index_peak - 1.0) {
      br = 1;
      d = 0;
      ar = 0;
      f = 0;
      ia = 0;
      h = 0;
      y = 0;
    } else {
      br = (int32_T)index_rise;
      d = (int32_T)(index_peak - 1.0);
      ar = (int32_T)index_rise - 1;
      f = (int32_T)(index_peak - 1.0);
      ia = (int32_T)index_rise - 1;
      h = (int32_T)(index_peak - 1.0);
      y = (int32_T)index_rise - 1;
    }

    tmp_size_idx_1 = f - ar;
    loop_ub = f - ar;
    for (f = 0; f < loop_ub; f++) {
      tmp_data[f] = (int16_T)((int16_T)(ar + f) + 1);
    }

    tmp_size_2[0] = 1;
    tmp_size_2[1] = tmp_size_idx_1;
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.tmp_data_c[f] = tmp_data[f];
    }

    Human_in_Loop_power(Human_in_Loop_B.tmp_data_c, tmp_size_2,
                        Human_in_Loop_B.tmp_data_k, tmp_size_3);
    tmp_size_idx_1 = h - ia;
    loop_ub = h - ia;
    for (f = 0; f < loop_ub; f++) {
      tmp_data[f] = (int16_T)((int16_T)(ia + f) + 1);
    }

    tmp_size_4[0] = 1;
    tmp_size_4[1] = tmp_size_idx_1;
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.tmp_data_c[f] = tmp_data[f];
    }

    Human_in_Loop_power_a(Human_in_Loop_B.tmp_data_c, tmp_size_4,
                          Human_in_Loop_B.tmp_data_cx, tmp_size_2);
    for (f = 0; f < b_n; f++) {
      Human_in_Loop_B.tmp_data[f << 2] = 1.0;
    }

    loop_ub = d - br;
    for (f = 0; f <= loop_ub; f++) {
      Human_in_Loop_B.tmp_data[1 + (f << 2)] = (int16_T)((int16_T)((br + f) - 1)
        + 1);
    }

    loop_ub = tmp_size_3[1];
    for (f = 0; f < loop_ub; f++) {
      Human_in_Loop_B.tmp_data[2 + (f << 2)] =
        Human_in_Loop_B.tmp_data_k[tmp_size_3[0] * f];
    }

    loop_ub = tmp_size_2[1];
    for (f = 0; f < loop_ub; f++) {
      Human_in_Loop_B.tmp_data[3 + (f << 2)] =
        Human_in_Loop_B.tmp_data_cx[tmp_size_2[0] * f];
    }

    br = b_n;
    for (f = 0; f < b_n; f++) {
      Human_in_Loop_B.result_data[f << 2] = Human_in_Loop_B.tmp_data[f << 2];
      Human_in_Loop_B.result_data[1 + (f << 2)] = Human_in_Loop_B.tmp_data[(f <<
        2) + 1];
      Human_in_Loop_B.result_data[2 + (f << 2)] = Human_in_Loop_B.tmp_data[(f <<
        2) + 2];
      Human_in_Loop_B.result_data[3 + (f << 2)] = Human_in_Loop_B.tmp_data[(f <<
        2) + 3];
    }

    tmp_0 = (int16_T)br;
    tmp_size_idx_1 = tmp_0;
    b_n = br - 1;
    if (0 <= tmp_size_idx_1 - 1) {
      memset(&Human_in_Loop_B.tmp_data_c[0], 0, tmp_size_idx_1 * sizeof(real_T));
    }

    if (br != 0) {
      for (br = 1; br - 1 <= b_n; br++) {
        for (d = br; d <= br; d++) {
          Human_in_Loop_B.tmp_data_c[d - 1] = 0.0;
        }
      }

      br = 0;
      for (d = 0; d <= b_n; d++) {
        ar = -1;
        for (f = br; f + 1 <= br + 4; f++) {
          if (Human_in_Loop_B.result_data[f] != 0.0) {
            ia = ar;
            for (h = d; h + 1 <= d + 1; h++) {
              ia++;
              Human_in_Loop_B.tmp_data_c[h] += Human_in_Loop_B.result_data[f] *
                parm1[ia];
            }
          }

          ar++;
        }

        br += 4;
      }
    }

    /* '<S9>:1:54' */
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.torque_track[y + f] = Human_in_Loop_B.tmp_data_c[f];
    }

    b = index_peak - index_rise;
    if (1.0 > b) {
      b_n = 0;
    } else {
      b_n = (int32_T)b;
    }

    if (index_rise > index_peak - 1.0) {
      br = 1;
      d = 0;
      ar = 0;
      f = 0;
      y = 0;
    } else {
      br = (int32_T)index_rise;
      d = (int32_T)(index_peak - 1.0);
      ar = (int32_T)index_rise - 1;
      f = (int32_T)(index_peak - 1.0);
      y = (int32_T)index_rise - 1;
    }

    tmp_size_idx_1 = f - ar;
    loop_ub = f - ar;
    for (f = 0; f < loop_ub; f++) {
      tmp_data[f] = (int16_T)((int16_T)(ar + f) + 1);
    }

    tmp_size_1[0] = 1;
    tmp_size_1[1] = tmp_size_idx_1;
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.tmp_data_c[f] = tmp_data[f];
    }

    Human_in_Loop_power(Human_in_Loop_B.tmp_data_c, tmp_size_1,
                        Human_in_Loop_B.tmp_data_k, tmp_size_2);
    for (f = 0; f < b_n; f++) {
      Human_in_Loop_B.tmp_data_m[3 * f] = 1.0;
    }

    loop_ub = d - br;
    for (f = 0; f <= loop_ub; f++) {
      Human_in_Loop_B.tmp_data_m[1 + 3 * f] = (real_T)(int16_T)((int16_T)((br +
        f) - 1) + 1) * 2.0;
    }

    loop_ub = tmp_size_2[1];
    for (f = 0; f < loop_ub; f++) {
      Human_in_Loop_B.tmp_data_m[2 + 3 * f] =
        Human_in_Loop_B.tmp_data_k[tmp_size_2[0] * f] * 3.0;
    }

    br = b_n;
    for (f = 0; f < b_n; f++) {
      Human_in_Loop_B.b_result_data[3 * f] = Human_in_Loop_B.tmp_data_m[3 * f];
      Human_in_Loop_B.b_result_data[1 + 3 * f] = Human_in_Loop_B.tmp_data_m[3 *
        f + 1];
      Human_in_Loop_B.b_result_data[2 + 3 * f] = Human_in_Loop_B.tmp_data_m[3 *
        f + 2];
    }

    tmp_0 = (int16_T)br;
    tmp_size_idx_1 = tmp_0;
    b_n = br - 1;
    if (0 <= tmp_size_idx_1 - 1) {
      memset(&Human_in_Loop_B.tmp_data_c[0], 0, tmp_size_idx_1 * sizeof(real_T));
    }

    if (br != 0) {
      for (br = 1; br - 1 <= b_n; br++) {
        for (d = br; d <= br; d++) {
          Human_in_Loop_B.tmp_data_c[d - 1] = 0.0;
        }
      }

      br = 0;
      for (d = 0; d <= b_n; d++) {
        ar = -1;
        for (f = br; f + 1 <= br + 3; f++) {
          if (Human_in_Loop_B.b_result_data[f] != 0.0) {
            ia = ar;
            for (h = d; h + 1 <= d + 1; h++) {
              ia++;
              Human_in_Loop_B.tmp_data_c[h] += parm1[1 + ia] *
                Human_in_Loop_B.b_result_data[f];
            }
          }

          ar++;
        }

        br += 3;
      }
    }

    /* '<S9>:1:55' */
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.torque_delta_track[y + f] = Human_in_Loop_B.tmp_data_c[f] *
        500.0;
    }

    /* '<S9>:1:57' */
    /* '<S9>:1:61' */
    /* '<S9>:1:62' */
    parm1[0] = peak_torque;
    parm1[1] = 0.0;
    parm1[2] = 0.0;
    parm1[3] = 0.0;
    tmp[0] = 1.0;
    tmp[4] = index_peak;
    tmp[8] = index_peak * index_peak;
    tmp[12] = rt_powd_snf(index_peak, 3.0);
    tmp[1] = 1.0;
    tmp[5] = index_fall;
    tmp[9] = index_fall * index_fall;
    tmp[13] = rt_powd_snf(index_fall, 3.0);
    tmp[2] = 0.0;
    tmp[6] = 1.0;
    tmp[10] = 2.0 * index_peak;
    tmp[14] = index_peak * index_peak * 3.0;
    tmp[3] = 0.0;
    tmp[7] = 1.0;
    tmp[11] = 2.0 * index_fall;
    tmp[15] = index_fall * index_fall * 3.0;
    Human_in_Loop_mldivide(tmp, parm1);
    peak_torque = (index_fall + 1.0) - index_peak;
    if (1.0 > peak_torque) {
      b_n = 0;
    } else {
      b_n = (int32_T)peak_torque;
    }

    if (index_peak > index_fall) {
      br = 1;
      d = 0;
      ar = 0;
      f = 0;
      ia = 0;
      h = 0;
      y = 0;
    } else {
      br = (int32_T)index_peak;
      d = (int32_T)index_fall;
      ar = (int32_T)index_peak - 1;
      f = (int32_T)index_fall;
      ia = (int32_T)index_peak - 1;
      h = (int32_T)index_fall;
      y = (int32_T)index_peak - 1;
    }

    tmp_size_idx_1 = f - ar;
    loop_ub = f - ar;
    for (f = 0; f < loop_ub; f++) {
      tmp_data[f] = (int16_T)((int16_T)(ar + f) + 1);
    }

    tmp_size[0] = 1;
    tmp_size[1] = tmp_size_idx_1;
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.tmp_data_c[f] = tmp_data[f];
    }

    Human_in_Loop_power(Human_in_Loop_B.tmp_data_c, tmp_size,
                        Human_in_Loop_B.tmp_data_k, tmp_size_2);
    tmp_size_idx_1 = h - ia;
    loop_ub = h - ia;
    for (f = 0; f < loop_ub; f++) {
      tmp_data[f] = (int16_T)((int16_T)(ia + f) + 1);
    }

    tmp_size_0[0] = 1;
    tmp_size_0[1] = tmp_size_idx_1;
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.tmp_data_c[f] = tmp_data[f];
    }

    Human_in_Loop_power_a(Human_in_Loop_B.tmp_data_c, tmp_size_0,
                          Human_in_Loop_B.tmp_data_cx, tmp_size_3);
    for (f = 0; f < b_n; f++) {
      Human_in_Loop_B.tmp_data[f << 2] = 1.0;
    }

    loop_ub = d - br;
    for (f = 0; f <= loop_ub; f++) {
      Human_in_Loop_B.tmp_data[1 + (f << 2)] = (int16_T)((int16_T)((br + f) - 1)
        + 1);
    }

    loop_ub = tmp_size_2[1];
    for (f = 0; f < loop_ub; f++) {
      Human_in_Loop_B.tmp_data[2 + (f << 2)] =
        Human_in_Loop_B.tmp_data_k[tmp_size_2[0] * f];
    }

    loop_ub = tmp_size_3[1];
    for (f = 0; f < loop_ub; f++) {
      Human_in_Loop_B.tmp_data[3 + (f << 2)] =
        Human_in_Loop_B.tmp_data_cx[tmp_size_3[0] * f];
    }

    br = b_n;
    for (f = 0; f < b_n; f++) {
      Human_in_Loop_B.result_data[f << 2] = Human_in_Loop_B.tmp_data[f << 2];
      Human_in_Loop_B.result_data[1 + (f << 2)] = Human_in_Loop_B.tmp_data[(f <<
        2) + 1];
      Human_in_Loop_B.result_data[2 + (f << 2)] = Human_in_Loop_B.tmp_data[(f <<
        2) + 2];
      Human_in_Loop_B.result_data[3 + (f << 2)] = Human_in_Loop_B.tmp_data[(f <<
        2) + 3];
    }

    tmp_0 = (int16_T)br;
    tmp_size_idx_1 = tmp_0;
    b_n = br - 1;
    if (0 <= tmp_size_idx_1 - 1) {
      memset(&Human_in_Loop_B.tmp_data_c[0], 0, tmp_size_idx_1 * sizeof(real_T));
    }

    if (br != 0) {
      for (br = 1; br - 1 <= b_n; br++) {
        for (d = br; d <= br; d++) {
          Human_in_Loop_B.tmp_data_c[d - 1] = 0.0;
        }
      }

      br = 0;
      for (d = 0; d <= b_n; d++) {
        ar = -1;
        for (f = br; f + 1 <= br + 4; f++) {
          if (Human_in_Loop_B.result_data[f] != 0.0) {
            ia = ar;
            for (h = d; h + 1 <= d + 1; h++) {
              ia++;
              Human_in_Loop_B.tmp_data_c[h] += Human_in_Loop_B.result_data[f] *
                parm1[ia];
            }
          }

          ar++;
        }

        br += 4;
      }
    }

    /* '<S9>:1:63' */
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.torque_track[y + f] = Human_in_Loop_B.tmp_data_c[f];
    }

    /* '<S9>:1:67' */
    /* '<S9>:1:68' */
    /* '<S9>:1:69' */
    /* '<S9>:1:70' */
    for (f = 0; f < 750; f++) {
      Human_in_Loop_DW.TorqueMem[f << 2] = Human_in_Loop_B.torque_track[f];
      Human_in_Loop_DW.TorqueMem[1 + (f << 2)] = 0.0;
      Human_in_Loop_DW.TorqueMem[2 + (f << 2)] =
        Human_in_Loop_B.torque_delta_track[f];
      Human_in_Loop_DW.TorqueMem[3 + (f << 2)] = 0.0;
    }
  }

  /* '<S9>:1:73' */
  Human_in_Loop_DW.last_footstate = footstate;

  /* '<S9>:1:74' */
  Human_in_Loop_DW.TorqueMem[1 + (((int32_T)stride_index - 1) << 2)] =
    torque_measure;

  /* '<S9>:1:75' */
  Human_in_Loop_DW.TorqueMem[3 + (((int32_T)stride_index - 1) << 2)] =
    troque_delta;
  if (mode == 2.0) {
    /* '<S9>:1:77' */
    /* '<S9>:1:78' */
    mode = Human_in_Loop_DW.TorqueMem[((int32_T)stride_index - 1) << 2];

    /* '<S9>:1:79' */
    stride_index = Human_in_Loop_DW.TorqueMem[(((int32_T)stride_index - 1) << 2)
      + 2];
  } else {
    /* '<S9>:1:81' */
    mode = 0.0;

    /* '<S9>:1:82' */
    stride_index = 0.0;
  }

  /* '<S9>:1:85' */
  /* '<S9>:1:86' */
  Human_in_Loop_B.torque_des = mode;
  Human_in_Loop_B.torque_delta_des = stride_index;
  for (f = 0; f < 750; f++) {
    Human_in_Loop_B.torque_trace[f << 1] = Human_in_Loop_DW.TorqueMem[f << 2];
    Human_in_Loop_B.torque_trace[1 + (f << 1)] = Human_in_Loop_DW.TorqueMem[(f <<
      2) + 1];
  }

  for (f = 0; f < 750; f++) {
    Human_in_Loop_B.torque_delta_trace[f << 1] = Human_in_Loop_DW.TorqueMem[(f <<
      2) + 2];
    Human_in_Loop_B.torque_delta_trace[1 + (f << 1)] =
      Human_in_Loop_DW.TorqueMem[(f << 2) + 3];
  }

  /* End of MATLAB Function: '<S1>/Torque track' */

  /* MATLAB Function: '<S1>/Controller' */
  /* MATLAB Function 'Control Module/Controller': '<S7>:1' */
  /* '<S7>:1:21' */
  /* '<S7>:1:22' */
  /* '<S7>:1:26' */
  /* '<S7>:1:27' */
  /* '<S7>:1:31' */
  /* '<S7>:1:33' */
  /* '<S7>:1:36' */
  /* '<S7>:1:37' */
  /* '<S7>:1:41' */
  /* '<S7>:1:42' */
  /* '<S7>:1:44' */
  /* '<S7>:1:45' */
  switch ((int32_T)Human_in_Loop_B.RT1_n[0]) {
   case 1:
    /* '<S7>:1:49' */
    Human_in_Loop_B.motor_vel_cmd = 0.0;

    /* '<S7>:1:50' */
    Human_in_Loop_DW.calib_state = 0.0;
    break;

   case 3:
    /* '<S7>:1:53' */
    Human_in_Loop_B.motor_vel_cmd = -Human_in_Loop_P.Controller_SLACK_SPEED *
      5.0;
    break;

   case 4:
    if ((Human_in_Loop_DW.calib_state == 0.0) && (Human_in_Loop_B.RT1[0] <
         Human_in_Loop_P.Controller_CALIB_TORQUE)) {
      /* '<S7>:1:56' */
      /* '<S7>:1:57' */
      Human_in_Loop_B.motor_vel_cmd = Human_in_Loop_P.Controller_CALIB_SPEED *
        5.0;
    } else if ((Human_in_Loop_DW.calib_state == 0.0) && (Human_in_Loop_B.RT1[0] >
                Human_in_Loop_P.Controller_CALIB_TORQUE)) {
      /* '<S7>:1:58' */
      /* '<S7>:1:59' */
      Human_in_Loop_DW.calib_state = 1.0;

      /* '<S7>:1:60' */
      Human_in_Loop_B.motor_vel_cmd = 0.0;
    } else {
      /* '<S7>:1:62' */
      Human_in_Loop_B.motor_vel_cmd = 0.0;
    }
    break;

   case 2:
    /* '<S7>:1:66' */
    switch (Human_in_Loop_P.Controller_run_mode) {
     case 1:
      if (Human_in_Loop_B.RT2[0] > 0.0) {
        /* '<S7>:1:69' */
        /* '<S7>:1:70' */
        stride_index = (Human_in_Loop_B.RT2[0] - Human_in_Loop_B.RT2[0] / 50.0 *
                        Human_in_Loop_P.Controller_FOLLOW_SLACK_ANGLE) -
          Human_in_Loop_B.RT3[0] * 0.33333333333333331;
      } else {
        /* '<S7>:1:72' */
        stride_index = Human_in_Loop_B.RT2[0] - Human_in_Loop_B.RT3[0] *
          0.33333333333333331;
      }

      /* '<S7>:1:74' */
      Human_in_Loop_B.motor_vel_cmd = Human_in_Loop_B.RT1_a[3] * stride_index *
        5.0 / 0.05;
      break;

     case 2:
      if (Human_in_Loop_B.RT1_n[1] == 1.0) {
        /* '<S7>:1:77' */
        /* '<S7>:1:78' */
        /* '<S7>:1:79' */
        /* '<S7>:1:80' */
        /* '<S7>:1:81' */
        Human_in_Loop_B.motor_vel_cmd = (((Human_in_Loop_B.torque_des -
          Human_in_Loop_B.RT1[0]) * Human_in_Loop_B.RT1_a[0] +
          (Human_in_Loop_B.torque_delta_des - Human_in_Loop_B.RT1[1]) *
          Human_in_Loop_B.RT1_a[1]) + Human_in_Loop_B.RT1_a[4] *
          Human_in_Loop_B.torque_delta_des) * 5.0 / 0.05;
      } else {
        /* '<S7>:1:83' */
        /* '<S7>:1:84' */
        Human_in_Loop_B.motor_vel_cmd = (0.0 - Human_in_Loop_B.RT3[0] *
          0.33333333333333331) * Human_in_Loop_B.RT1_a[3] * 5.0 / 0.05;
      }
      break;

     default:
      /* '<S7>:1:87' */
      Human_in_Loop_B.motor_vel_cmd = 0.0;
      break;
    }
    break;

   case 0:
    /* '<S7>:1:91' */
    Human_in_Loop_B.motor_vel_cmd = 0.0;
    break;

   default:
    /* '<S7>:1:94' */
    Human_in_Loop_B.motor_vel_cmd = 0.0;
    break;
  }

  /* End of MATLAB Function: '<S1>/Controller' */

  /* MATLAB Function: '<S8>/MATLAB Function' */
  /* MATLAB Function 'Control Module/Motor/MATLAB Function': '<S11>:1' */
  /* '<S11>:1:3' */
  /* '<S11>:1:4' */
  if (Human_in_Loop_B.RT1[0] > Human_in_Loop_P.MATLABFunction_MAX_TORQUE) {
    /* '<S11>:1:6' */
    /* '<S11>:1:7' */
    Human_in_Loop_B.vel = 0.0;
  } else if (Human_in_Loop_B.RT3[0] >
             Human_in_Loop_P.MATLABFunction_MAX_MOTOR_ANGLE) {
    /* '<S11>:1:8' */
    /* '<S11>:1:9' */
    Human_in_Loop_B.vel = 0.0;
  } else if (Human_in_Loop_B.RT3[0] <
             Human_in_Loop_P.MATLABFunction_MIN_MOTOR_ANGLE) {
    /* '<S11>:1:10' */
    /* '<S11>:1:11' */
    Human_in_Loop_B.vel = 0.0;
  } else if (Human_in_Loop_B.motor_vel_cmd >
             Human_in_Loop_P.MATLABFunction_MAX_SPEED) {
    /* '<S11>:1:12' */
    /* '<S11>:1:13' */
    Human_in_Loop_B.vel = Human_in_Loop_P.MATLABFunction_MAX_SPEED;
  } else if (Human_in_Loop_B.motor_vel_cmd <
             -Human_in_Loop_P.MATLABFunction_MAX_SPEED) {
    /* '<S11>:1:14' */
    /* '<S11>:1:15' */
    Human_in_Loop_B.vel = -Human_in_Loop_P.MATLABFunction_MAX_SPEED;
  } else {
    /* '<S11>:1:17' */
    Human_in_Loop_B.vel = Human_in_Loop_B.motor_vel_cmd;
  }

  /* End of MATLAB Function: '<S8>/MATLAB Function' */

  /* Gain: '<S8>/Gain2' */
  Human_in_Loop_B.Gain2_j = Human_in_Loop_P.Gain2_Gain * Human_in_Loop_B.vel;

  /* Gain: '<S8>/Gain1' */
  Human_in_Loop_B.Gain1_p = Human_in_Loop_P.Gain1_Gain * Human_in_Loop_B.Gain2_j;

  /* S-Function (rti_commonblock): '<S10>/S-Function1' */
  /* This comment workarounds a code generation problem */

  /* --- Human_in_Loop/Control Module/Motor/DAC_CLASS1_BL1 --- */
  /* --- [RTI120X, DAC C1] - Channel: 16 --- */
  {
    /* define variables required for DAC realtime functions */
    Float64 inportDacData= 0.0;
    inportDacData = (real_T) Human_in_Loop_B.Gain1_p;

    /* write value of CL1 DAC for output channel 16 */
    DacCl1AnalogOut_setOutputValue(pRTIDacC1AnalogOut_Ch_16,
      DAC_CLASS1_CHANNEL_16, inportDacData);
    DacCl1AnalogOut_write(pRTIDacC1AnalogOut_Ch_16);
  }
}

/* Termination for function-call system: '<Root>/Control Module' */
void Human_in_Loo_ControlModule_Term(void)
{
  /* Terminate for S-Function (rti_commonblock): '<S10>/S-Function1' */

  /* --- Human_in_Loop/Control Module/Motor/DAC_CLASS1_BL1 --- */
  /* --- [RTI120X, DAC C1] - Channel: 16 --- */

  /* All channel outputs are set to high impedance state */
  DacCl1AnalogOut_setOutputHighZ(pRTIDacC1AnalogOut_Ch_16, DAC_CLASS1_HIGH_Z_ON);

  /* Deactivates AnalogOut functionality */
  DacCl1AnalogOut_stop(pRTIDacC1AnalogOut_Ch_16);
}

/* Function for MATLAB Function: '<S21>/MATLAB Function' */
static real_T Human_in_Loop_xnrm2(const real_T x[30], int32_T ix0)
{
  real_T y;
  real_T scale;
  real_T absxk;
  real_T t;
  int32_T k;
  y = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = ix0; k <= ix0 + 14; k++) {
    absxk = fabs(x[k - 1]);
    if (absxk > scale) {
      t = scale / absxk;
      y = y * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      y += t * t;
    }
  }

  return scale * sqrt(y);
}

/* Function for MATLAB Function: '<S21>/MATLAB Function' */
static real_T Human_in_Loop_xnrm2_o(int32_T n, const real_T x[30], int32_T ix0)
{
  real_T y;
  real_T scale;
  int32_T kend;
  real_T absxk;
  real_T t;
  int32_T k;
  y = 0.0;
  if (!(n < 1)) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

real_T rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T a;
  a = fabs(u0);
  y = fabs(u1);
  if (a < y) {
    a /= y;
    y *= sqrt(a * a + 1.0);
  } else if (a > y) {
    y /= a;
    y = sqrt(y * y + 1.0) * a;
  } else {
    if (!rtIsNaN(y)) {
      y = a * 1.4142135623730951;
    }
  }

  return y;
}

/* Function for MATLAB Function: '<S21>/MATLAB Function' */
static void Human_in_Loop_xgeqp3(real_T A[30], real_T tau[2], int32_T jpvt[2])
{
  real_T work[2];
  real_T vn1[2];
  real_T vn2[2];
  int32_T k;
  int32_T i_i;
  int32_T pvt;
  int32_T ix;
  real_T smax;
  int32_T iy;
  real_T xnorm;
  int32_T c_ix;
  real_T b_c;
  int32_T d_ix;
  int32_T f;
  int32_T ijA;
  int32_T exitg1;
  boolean_T exitg2;
  jpvt[0] = 1;
  work[0] = 0.0;
  smax = Human_in_Loop_xnrm2(A, 1);
  vn2[0] = smax;
  vn1[0] = smax;
  jpvt[1] = 2;
  work[1] = 0.0;
  smax = Human_in_Loop_xnrm2(A, 16);
  vn2[1] = smax;
  vn1[1] = smax;
  for (k = 0; k < 2; k++) {
    i_i = k * 15 + k;
    pvt = 0;
    if (2 - k > 1) {
      ix = k;
      smax = fabs(vn1[k]);
      iy = 2;
      while (iy <= 2 - k) {
        ix++;
        b_c = fabs(vn1[ix]);
        if (b_c > smax) {
          pvt = 1;
          smax = b_c;
        }

        iy = 3;
      }
    }

    pvt += k;
    if (pvt + 1 != k + 1) {
      ix = 15 * pvt;
      iy = 15 * k;
      for (c_ix = 0; c_ix < 15; c_ix++) {
        smax = A[ix];
        A[ix] = A[iy];
        A[iy] = smax;
        ix++;
        iy++;
      }

      ix = jpvt[pvt];
      jpvt[pvt] = jpvt[k];
      jpvt[k] = ix;
      vn1[pvt] = vn1[k];
      vn2[pvt] = vn2[k];
    }

    smax = A[i_i];
    b_c = 0.0;
    xnorm = Human_in_Loop_xnrm2_o(14 - k, A, i_i + 2);
    if (xnorm != 0.0) {
      xnorm = rt_hypotd_snf(A[i_i], xnorm);
      if (A[i_i] >= 0.0) {
        xnorm = -xnorm;
      }

      if (fabs(xnorm) < 1.0020841800044864E-292) {
        pvt = 0;
        ix = (i_i - k) + 15;
        do {
          pvt++;
          for (iy = i_i + 1; iy + 1 <= ix; iy++) {
            A[iy] *= 9.9792015476736E+291;
          }

          xnorm *= 9.9792015476736E+291;
          smax *= 9.9792015476736E+291;
        } while (!(fabs(xnorm) >= 1.0020841800044864E-292));

        xnorm = rt_hypotd_snf(smax, Human_in_Loop_xnrm2_o(14 - k, A, i_i + 2));
        if (smax >= 0.0) {
          xnorm = -xnorm;
        }

        b_c = (xnorm - smax) / xnorm;
        smax = 1.0 / (smax - xnorm);
        ix = (i_i - k) + 15;
        for (iy = i_i + 1; iy + 1 <= ix; iy++) {
          A[iy] *= smax;
        }

        for (ix = 1; ix <= pvt; ix++) {
          xnorm *= 1.0020841800044864E-292;
        }

        smax = xnorm;
      } else {
        b_c = (xnorm - A[i_i]) / xnorm;
        smax = 1.0 / (A[i_i] - xnorm);
        pvt = (i_i - k) + 15;
        for (ix = i_i + 1; ix + 1 <= pvt; ix++) {
          A[ix] *= smax;
        }

        smax = xnorm;
      }
    }

    tau[k] = b_c;
    A[i_i] = smax;
    if (k + 1 < 2) {
      smax = A[i_i];
      A[i_i] = 1.0;
      if (tau[0] != 0.0) {
        pvt = 30;
        ix = i_i + 14;
        while ((pvt - 15 > 0) && (A[ix] == 0.0)) {
          pvt--;
          ix--;
        }

        ix = 1;
        exitg2 = false;
        while ((!exitg2) && (ix > 0)) {
          iy = 16;
          do {
            exitg1 = 0;
            if (iy <= pvt) {
              if (A[iy - 1] != 0.0) {
                exitg1 = 1;
              } else {
                iy++;
              }
            } else {
              ix = 0;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        pvt = 15;
        ix = 0;
      }

      if (pvt - 15 > 0) {
        if (ix != 0) {
          work[0] = 0.0;
          iy = 0;
          c_ix = 16;
          while (c_ix <= 16) {
            c_ix = i_i;
            b_c = 0.0;
            for (d_ix = 16; d_ix <= pvt; d_ix++) {
              b_c += A[d_ix - 1] * A[c_ix];
              c_ix++;
            }

            work[iy] += b_c;
            iy++;
            c_ix = 31;
          }
        }

        if (!(-tau[0] == 0.0)) {
          iy = 0;
          c_ix = 0;
          d_ix = 1;
          while (d_ix <= ix) {
            if (work[c_ix] != 0.0) {
              b_c = work[c_ix] * -tau[0];
              d_ix = i_i;
              f = pvt + iy;
              for (ijA = iy + 15; ijA + 1 <= f; ijA++) {
                A[ijA] += A[d_ix] * b_c;
                d_ix++;
              }
            }

            c_ix++;
            iy += 15;
            d_ix = 2;
          }
        }
      }

      A[i_i] = smax;
    }

    i_i = k + 2;
    while (i_i < 3) {
      if (vn1[1] != 0.0) {
        smax = fabs(A[15 + k]) / vn1[1];
        smax = 1.0 - smax * smax;
        if (smax < 0.0) {
          smax = 0.0;
        }

        b_c = vn1[1] / vn2[1];
        b_c = b_c * b_c * smax;
        if (b_c <= 1.4901161193847656E-8) {
          vn1[1] = Human_in_Loop_xnrm2_o(14 - k, A, k + 17);
          vn2[1] = vn1[1];
        } else {
          vn1[1] *= sqrt(smax);
        }
      }

      i_i = 3;
    }
  }
}

/* Model output function */
void Human_in_Loop_output(void)
{
  real_T x[2];
  real_T A[30];
  real_T tau[2];
  int32_T jpvt[2];
  real_T tol;
  real_T B[15];
  int32_T c_i;
  static const int8_T b_A[30] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

  real_T B_0;
  int32_T i;

  /* S-Function (rti_commonblock): '<S38>/S-Function1' */
  /* This comment workarounds a code generation problem */

  /* --- Human_in_Loop/Sensor Data/Torque module/ADC_CLASS1_BL6 --- */
  /* --- [RTI120X, ADC C1] - Channel: 6 --- */
  {
    UInt32 readStateFlag[] = { 1 };

    /* define variable required for adc cl1 realtime functions */
    UInt32 IsNew = 0;

    /* wait until conversion result is available */
    while (IsNew == 0) {
      AdcCl1AnalogIn_isDataReady(pRTIAdcC1AnalogIn_Ch_6, &IsNew);
    }

    /* read conversion result from hardware */
    AdcCl1AnalogIn_read(pRTIAdcC1AnalogIn_Ch_6);

    /* retrieve conversion result */
    AdcCl1AnalogIn_getSingleValue(pRTIAdcC1AnalogIn_Ch_6, readStateFlag, (real_T*)
      &Human_in_Loop_B.SFunction1);
  }

  /* Gain: '<S21>/Gain' */
  Human_in_Loop_B.Gain = Human_in_Loop_P.Gain_Gain * Human_in_Loop_B.SFunction1;

  /* DiscreteFilter: '<S21>/0.4low2' */
  tol = Human_in_Loop_B.Gain;
  tol -= Human_in_Loop_P.u4low2_DenCoef[1] * Human_in_Loop_DW.u4low2_states[0];
  tol -= Human_in_Loop_P.u4low2_DenCoef[2] * Human_in_Loop_DW.u4low2_states[1];
  tol -= Human_in_Loop_P.u4low2_DenCoef[3] * Human_in_Loop_DW.u4low2_states[2];
  tol /= Human_in_Loop_P.u4low2_DenCoef[0];
  Human_in_Loop_DW.u4low2_tmp = tol;
  tol = Human_in_Loop_P.u4low2_NumCoef[0] * Human_in_Loop_DW.u4low2_tmp;
  tol += Human_in_Loop_P.u4low2_NumCoef[1] * Human_in_Loop_DW.u4low2_states[0];
  tol += Human_in_Loop_P.u4low2_NumCoef[2] * Human_in_Loop_DW.u4low2_states[1];
  tol += Human_in_Loop_P.u4low2_NumCoef[3] * Human_in_Loop_DW.u4low2_states[2];
  Human_in_Loop_B.u4low2 = tol;

  /* MATLAB Function: '<S21>/Data process' */
  /* MATLAB Function 'Sensor Data/Torque module/Data process': '<S39>:1' */
  if (Human_in_Loop_P.Dataprocess_BT_RESET_TORQUE) {
    /* '<S39>:1:10' */
    /* '<S39>:1:11' */
    Human_in_Loop_DW.torque_zero = Human_in_Loop_B.u4low2 *
      Human_in_Loop_P.Dataprocess_load_vol_gain +
      Human_in_Loop_P.Dataprocess_load_vol_offset;
  }

  /* '<S39>:1:14' */
  Human_in_Loop_B.torque = (Human_in_Loop_B.u4low2 *
    Human_in_Loop_P.Dataprocess_load_vol_gain +
    Human_in_Loop_P.Dataprocess_load_vol_offset) - Human_in_Loop_DW.torque_zero;

  /* End of MATLAB Function: '<S21>/Data process' */

  /* UnitDelay: '<S37>/Unit Delay1' */
  Human_in_Loop_B.x2k1 = Human_in_Loop_DW.UnitDelay1_DSTATE;

  /* UnitDelay: '<S37>/Unit Delay' */
  Human_in_Loop_B.x1k1 = Human_in_Loop_DW.UnitDelay_DSTATE;

  /* Gain: '<S37>/Gain1' */
  B_0 = Human_in_Loop_P.uOrderTD_T1 * Human_in_Loop_P.uOrderTD_T2;
  tol = 1.0 / B_0;
  Human_in_Loop_B.Gain1 = tol * Human_in_Loop_B.x1k1;

  /* Gain: '<S37>/Gain2' */
  tol = Human_in_Loop_P.uOrderTD_T1 + Human_in_Loop_P.uOrderTD_T2;
  B_0 = Human_in_Loop_P.uOrderTD_T1 * Human_in_Loop_P.uOrderTD_T2;
  tol /= B_0;
  Human_in_Loop_B.Gain2 = tol * Human_in_Loop_B.x2k1;

  /* UnitDelay: '<S37>/Unit Delay2' */
  Human_in_Loop_B.UnitDelay2 = Human_in_Loop_DW.UnitDelay2_DSTATE;

  /* Gain: '<S37>/Gain4' */
  B_0 = Human_in_Loop_P.uOrderTD_T1 * Human_in_Loop_P.uOrderTD_T2;
  tol = 1.0 / B_0;
  Human_in_Loop_B.Gain4 = tol * Human_in_Loop_B.UnitDelay2;

  /* Sum: '<S37>/Add2' */
  Human_in_Loop_B.Add2 = (Human_in_Loop_B.Gain1 + Human_in_Loop_B.Gain2) -
    Human_in_Loop_B.Gain4;

  /* Gain: '<S37>/Gain3' */
  Human_in_Loop_B.Gain3 = Human_in_Loop_P.uOrderTD_Ts * Human_in_Loop_B.Add2;

  /* Sum: '<S37>/Add1' */
  Human_in_Loop_B.x2k = Human_in_Loop_B.x2k1 - Human_in_Loop_B.Gain3;

  /* MATLAB Function: '<S21>/Mux' */
  Human_in_Loop_Mux(Human_in_Loop_B.torque, Human_in_Loop_B.x2k,
                    &Human_in_Loop_B.sf_Mux);

  /* RateTransition: '<S4>/RT1' */
  switch (Human_in_Loop_DW.RT1_read_buf) {
   case 0:
    Human_in_Loop_DW.RT1_write_buf = 1;
    break;

   case 1:
    Human_in_Loop_DW.RT1_write_buf = 0;
    break;

   default:
    Human_in_Loop_DW.RT1_write_buf = (int8_T)(Human_in_Loop_DW.RT1_last_buf_wr ==
      0);
    break;
  }

  if (Human_in_Loop_DW.RT1_write_buf != 0) {
    Human_in_Loop_DW.RT1_Buffer1[0] = Human_in_Loop_B.sf_Mux.x[0];
    Human_in_Loop_DW.RT1_Buffer1[1] = Human_in_Loop_B.sf_Mux.x[1];
  } else {
    Human_in_Loop_DW.RT1_Buffer0[0] = Human_in_Loop_B.sf_Mux.x[0];
    Human_in_Loop_DW.RT1_Buffer0[1] = Human_in_Loop_B.sf_Mux.x[1];
  }

  Human_in_Loop_DW.RT1_last_buf_wr = Human_in_Loop_DW.RT1_write_buf;
  Human_in_Loop_DW.RT1_write_buf = -1;

  /* End of RateTransition: '<S4>/RT1' */

  /* S-Function (rti_commonblock): '<S27>/S-Function1' */
  /* This comment workarounds a code generation problem */

  /* Gain: '<S19>/Gain' */
  Human_in_Loop_B.Gain_l = Human_in_Loop_P.Gain_Gain_c *
    Human_in_Loop_B.SFunction1_o1;

  /* MATLAB Function: '<S19>/Data process' */
  /* MATLAB Function 'Sensor Data/Encoder module/Data process': '<S24>:1' */
  if (Human_in_Loop_P.Dataprocess_BT_RESET_ANKLE) {
    /* '<S24>:1:8' */
    /* '<S24>:1:9' */
    Human_in_Loop_DW.angle_zero_f = Human_in_Loop_B.Gain_l;
  }

  /* '<S24>:1:12' */
  Human_in_Loop_B.angle_m = Human_in_Loop_B.Gain_l -
    Human_in_Loop_DW.angle_zero_f;

  /* End of MATLAB Function: '<S19>/Data process' */

  /* Gain: '<S19>/Gain1' */
  Human_in_Loop_B.Gain1_b = Human_in_Loop_P.Gain1_Gain_e *
    Human_in_Loop_B.SFunction1_o2;

  /* MATLAB Function: '<S19>/Mux' */
  Human_in_Loop_Mux(Human_in_Loop_B.angle_m, Human_in_Loop_B.Gain1_b,
                    &Human_in_Loop_B.sf_Mux_p);

  /* RateTransition: '<S4>/RT2' */
  switch (Human_in_Loop_DW.RT2_read_buf) {
   case 0:
    Human_in_Loop_DW.RT2_write_buf = 1;
    break;

   case 1:
    Human_in_Loop_DW.RT2_write_buf = 0;
    break;

   default:
    Human_in_Loop_DW.RT2_write_buf = (int8_T)(Human_in_Loop_DW.RT2_last_buf_wr ==
      0);
    break;
  }

  if (Human_in_Loop_DW.RT2_write_buf != 0) {
    Human_in_Loop_DW.RT2_Buffer1[0] = Human_in_Loop_B.sf_Mux_p.x[0];
    Human_in_Loop_DW.RT2_Buffer1[1] = Human_in_Loop_B.sf_Mux_p.x[1];
  } else {
    Human_in_Loop_DW.RT2_Buffer0[0] = Human_in_Loop_B.sf_Mux_p.x[0];
    Human_in_Loop_DW.RT2_Buffer0[1] = Human_in_Loop_B.sf_Mux_p.x[1];
  }

  Human_in_Loop_DW.RT2_last_buf_wr = Human_in_Loop_DW.RT2_write_buf;
  Human_in_Loop_DW.RT2_write_buf = -1;

  /* End of RateTransition: '<S4>/RT2' */

  /* S-Function (rti_commonblock): '<S28>/S-Function1' */
  /* This comment workarounds a code generation problem */

  /* Gain: '<S19>/Gain2' */
  Human_in_Loop_B.Gain2_h = Human_in_Loop_P.Gain2_Gain_e *
    Human_in_Loop_B.SFunction1_o1_k;

  /* MATLAB Function: '<S19>/Data process1' */
  /* MATLAB Function 'Sensor Data/Encoder module/Data process1': '<S25>:1' */
  if (Human_in_Loop_P.Dataprocess1_BT_RESET_MOTOR) {
    /* '<S25>:1:8' */
    /* '<S25>:1:9' */
    Human_in_Loop_DW.angle_zero = Human_in_Loop_B.Gain2_h;
  }

  /* '<S25>:1:12' */
  Human_in_Loop_B.angle = Human_in_Loop_B.Gain2_h - Human_in_Loop_DW.angle_zero;

  /* End of MATLAB Function: '<S19>/Data process1' */

  /* Gain: '<S19>/Gain3' */
  Human_in_Loop_B.Gain3_m = Human_in_Loop_P.Gain3_Gain *
    Human_in_Loop_B.SFunction1_o2_k;

  /* UnitDelay: '<S23>/Unit Delay1' */
  Human_in_Loop_B.x2k1_k = Human_in_Loop_DW.UnitDelay1_DSTATE_h;

  /* UnitDelay: '<S23>/Unit Delay' */
  Human_in_Loop_B.x1k1_m = Human_in_Loop_DW.UnitDelay_DSTATE_e;

  /* Gain: '<S23>/Gain1' */
  B_0 = Human_in_Loop_P.uOrderTD_T1_l * Human_in_Loop_P.uOrderTD_T2_h;
  tol = 1.0 / B_0;
  Human_in_Loop_B.Gain1_m = tol * Human_in_Loop_B.x1k1_m;

  /* Gain: '<S23>/Gain2' */
  tol = Human_in_Loop_P.uOrderTD_T1_l + Human_in_Loop_P.uOrderTD_T2_h;
  B_0 = Human_in_Loop_P.uOrderTD_T1_l * Human_in_Loop_P.uOrderTD_T2_h;
  tol /= B_0;
  Human_in_Loop_B.Gain2_c = tol * Human_in_Loop_B.x2k1_k;

  /* UnitDelay: '<S23>/Unit Delay2' */
  Human_in_Loop_B.UnitDelay2_b = Human_in_Loop_DW.UnitDelay2_DSTATE_a;

  /* Gain: '<S23>/Gain4' */
  B_0 = Human_in_Loop_P.uOrderTD_T1_l * Human_in_Loop_P.uOrderTD_T2_h;
  tol = 1.0 / B_0;
  Human_in_Loop_B.Gain4_g = tol * Human_in_Loop_B.UnitDelay2_b;

  /* Sum: '<S23>/Add2' */
  Human_in_Loop_B.Add2_m = (Human_in_Loop_B.Gain1_m + Human_in_Loop_B.Gain2_c) -
    Human_in_Loop_B.Gain4_g;

  /* Gain: '<S23>/Gain3' */
  Human_in_Loop_B.Gain3_k = Human_in_Loop_P.uOrderTD_Ts_n *
    Human_in_Loop_B.Add2_m;

  /* Sum: '<S23>/Add1' */
  Human_in_Loop_B.x2k_i = Human_in_Loop_B.x2k1_k - Human_in_Loop_B.Gain3_k;

  /* MATLAB Function: '<S19>/Mux1' */
  /* MATLAB Function 'Sensor Data/Encoder module/Mux1': '<S32>:1' */
  /* '<S32>:1:3' */
  Human_in_Loop_B.x_i[0] = Human_in_Loop_B.angle;
  Human_in_Loop_B.x_i[1] = Human_in_Loop_B.Gain3_m;
  Human_in_Loop_B.x_i[2] = Human_in_Loop_B.x2k_i;

  /* RateTransition: '<S4>/RT3' */
  switch (Human_in_Loop_DW.RT3_read_buf) {
   case 0:
    Human_in_Loop_DW.RT3_write_buf = 1;
    break;

   case 1:
    Human_in_Loop_DW.RT3_write_buf = 0;
    break;

   default:
    Human_in_Loop_DW.RT3_write_buf = (int8_T)(Human_in_Loop_DW.RT3_last_buf_wr ==
      0);
    break;
  }

  if (Human_in_Loop_DW.RT3_write_buf != 0) {
    Human_in_Loop_DW.RT3_Buffer1[0] = Human_in_Loop_B.x_i[0];
    Human_in_Loop_DW.RT3_Buffer1[1] = Human_in_Loop_B.x_i[1];
    Human_in_Loop_DW.RT3_Buffer1[2] = Human_in_Loop_B.x_i[2];
  } else {
    Human_in_Loop_DW.RT3_Buffer0[0] = Human_in_Loop_B.x_i[0];
    Human_in_Loop_DW.RT3_Buffer0[1] = Human_in_Loop_B.x_i[1];
    Human_in_Loop_DW.RT3_Buffer0[2] = Human_in_Loop_B.x_i[2];
  }

  Human_in_Loop_DW.RT3_last_buf_wr = Human_in_Loop_DW.RT3_write_buf;
  Human_in_Loop_DW.RT3_write_buf = -1;

  /* End of RateTransition: '<S4>/RT3' */

  /* S-Function (rti_commonblock): '<S34>/S-Function1' */
  /* This comment workarounds a code generation problem */

  /* MATLAB Function: '<S20>/FootSwitch Filter' */
  /* MATLAB Function 'Sensor Data/FootSwitch module/FootSwitch Filter': '<S36>:1' */
  /* '<S36>:1:6' */
  if (Human_in_Loop_DW.foot_state == 0.0) {
    /* '<S36>:1:14' */
    if (Human_in_Loop_B.SFunction1_a) {
      /* '<S36>:1:15' */
      /* '<S36>:1:16' */
      Human_in_Loop_DW.foot_state = 1.0;
    } else {
      /* '<S36>:1:18' */
      Human_in_Loop_DW.foot_state = 0.0;
    }
  } else if (Human_in_Loop_B.SFunction1_a) {
    /* '<S36>:1:21' */
    /* '<S36>:1:22' */
    Human_in_Loop_DW.foot_state = 1.0;
  } else {
    /* '<S36>:1:24' */
    Human_in_Loop_DW.filter_time += 0.0002;
    if (Human_in_Loop_DW.filter_time > 0.1) {
      /* '<S36>:1:25' */
      /* '<S36>:1:26' */
      Human_in_Loop_DW.filter_time = 0.0;

      /* '<S36>:1:27' */
      Human_in_Loop_DW.foot_state = 0.0;
    }
  }

  /* '<S36>:1:32' */
  Human_in_Loop_B.state_c = Human_in_Loop_DW.foot_state;

  /* End of MATLAB Function: '<S20>/FootSwitch Filter' */

  /* MATLAB Function: '<S5>/State Machine' */
  /* MATLAB Function 'State Module/State Machine': '<S43>:1' */
  /* '<S43>:1:20' */
  /* '<S43>:1:21' */
  /* '<S43>:1:22' */
  /* '<S43>:1:23' */
  /* '<S43>:1:26' */
  /* '<S43>:1:27' */
  /* '<S43>:1:28' */
  /* '<S43>:1:29' */
  /* '<S43>:1:30' */
  if (Human_in_Loop_P.StateMachine_BT_RUN) {
    /* '<S43>:1:33' */
    /* '<S43>:1:34' */
    Human_in_Loop_DW.bt_run = 1.0;
  }

  if (Human_in_Loop_P.StateMachine_BT_CALIB) {
    /* '<S43>:1:36' */
    /* '<S43>:1:37' */
    Human_in_Loop_DW.reg_mode = 4.0;
  }

  if (Human_in_Loop_P.StateMachine_BT_SLACK) {
    /* '<S43>:1:39' */
    /* '<S43>:1:40' */
    Human_in_Loop_DW.reg_mode = 3.0;
  }

  if (Human_in_Loop_P.StateMachine_BT_IDLE) {
    /* '<S43>:1:42' */
    /* '<S43>:1:43' */
    Human_in_Loop_DW.reg_mode = 1.0;
  }

  if (Human_in_Loop_P.StateMachine_BT_ERROR) {
    /* '<S43>:1:45' */
    /* '<S43>:1:46' */
    Human_in_Loop_DW.reg_mode = 0.0;
  }

  if ((Human_in_Loop_DW.bt_run == 1.0) && (Human_in_Loop_DW.reg_last_switch ==
       0.0) && (Human_in_Loop_B.state_c == 1.0)) {
    /* '<S43>:1:49' */
    /* '<S43>:1:50' */
    Human_in_Loop_DW.reg_mode = 2.0;

    /* '<S43>:1:51' */
    Human_in_Loop_DW.bt_run = 0.0;
  }

  if ((Human_in_Loop_DW.reg_mode == 2.0) || (Human_in_Loop_DW.reg_mode == 1.0))
  {
    /* '<S43>:1:54' */
    if ((Human_in_Loop_DW.reg_last_switch == 0.0) && (Human_in_Loop_B.state_c ==
         1.0)) {
      /* '<S43>:1:55' */
      /* '<S43>:1:56' */
      Human_in_Loop_DW.reg_state = 1.0;

      /* '<S43>:1:57' */
      Human_in_Loop_DW.reg_stride_time = 0.618 *
        Human_in_Loop_DW.reg_stride_time + 0.382 *
        Human_in_Loop_DW.reg_stride_time_count;

      /* '<S43>:1:58' */
      Human_in_Loop_DW.reg_stride_time_count = 0.0;
    } else if ((Human_in_Loop_DW.reg_state == 1.0) &&
               (Human_in_Loop_DW.reg_stride_time_count > 0.6 *
                Human_in_Loop_DW.reg_stride_time)) {
      /* '<S43>:1:59' */
      /* '<S43>:1:60' */
      Human_in_Loop_DW.reg_state = 0.0;

      /* '<S43>:1:61' */
      Human_in_Loop_DW.reg_stride_time_count += 0.0002;
    } else {
      /* '<S43>:1:63' */
      Human_in_Loop_DW.reg_stride_time_count += 0.0002;
    }
  }

  /* '<S43>:1:67' */
  Human_in_Loop_DW.reg_last_switch = Human_in_Loop_B.state_c;
  if (Human_in_Loop_DW.reg_stride_time > 1.5) {
    /* '<S43>:1:68' */
    /* '<S43>:1:69' */
    Human_in_Loop_DW.reg_stride_time = 1.5;
  } else {
    if (Human_in_Loop_DW.reg_stride_time < 0.5) {
      /* '<S43>:1:70' */
      /* '<S43>:1:71' */
      Human_in_Loop_DW.reg_stride_time = 0.5;
    }
  }

  /* '<S43>:1:74' */
  Human_in_Loop_B.mode = Human_in_Loop_DW.reg_mode;

  /* '<S43>:1:75' */
  Human_in_Loop_B.state = Human_in_Loop_DW.reg_state;

  /* '<S43>:1:76' */
  Human_in_Loop_B.stride_time = Human_in_Loop_DW.reg_stride_time;

  /* '<S43>:1:77' */
  Human_in_Loop_B.stride_timer = Human_in_Loop_DW.reg_stride_time_count;

  /* End of MATLAB Function: '<S5>/State Machine' */

  /* MATLAB Function: '<S5>/Mux1' */
  /* MATLAB Function 'State Module/Mux1': '<S42>:1' */
  /* '<S42>:1:3' */
  Human_in_Loop_B.x[0] = Human_in_Loop_B.mode;
  Human_in_Loop_B.x[1] = Human_in_Loop_B.state;
  Human_in_Loop_B.x[2] = Human_in_Loop_B.stride_time;
  Human_in_Loop_B.x[3] = Human_in_Loop_B.stride_timer;

  /* RateTransition: '<S5>/RT1' */
  switch (Human_in_Loop_DW.RT1_read_buf_l) {
   case 0:
    Human_in_Loop_DW.RT1_write_buf_n = 1;
    break;

   case 1:
    Human_in_Loop_DW.RT1_write_buf_n = 0;
    break;

   default:
    Human_in_Loop_DW.RT1_write_buf_n = (int8_T)
      (Human_in_Loop_DW.RT1_last_buf_wr_m == 0);
    break;
  }

  if (Human_in_Loop_DW.RT1_write_buf_n != 0) {
    Human_in_Loop_DW.RT1_Buffer1_f[0] = Human_in_Loop_B.x[0];
    Human_in_Loop_DW.RT1_Buffer1_f[1] = Human_in_Loop_B.x[1];
    Human_in_Loop_DW.RT1_Buffer1_f[2] = Human_in_Loop_B.x[2];
    Human_in_Loop_DW.RT1_Buffer1_f[3] = Human_in_Loop_B.x[3];
  } else {
    Human_in_Loop_DW.RT1_Buffer0_a[0] = Human_in_Loop_B.x[0];
    Human_in_Loop_DW.RT1_Buffer0_a[1] = Human_in_Loop_B.x[1];
    Human_in_Loop_DW.RT1_Buffer0_a[2] = Human_in_Loop_B.x[2];
    Human_in_Loop_DW.RT1_Buffer0_a[3] = Human_in_Loop_B.x[3];
  }

  Human_in_Loop_DW.RT1_last_buf_wr_m = Human_in_Loop_DW.RT1_write_buf_n;
  Human_in_Loop_DW.RT1_write_buf_n = -1;

  /* End of RateTransition: '<S5>/RT1' */

  /* MATLAB Function: '<S12>/Mux1' incorporates:
   *  Constant: '<S12>/Kd'
   *  Constant: '<S12>/Kl'
   *  Constant: '<S12>/Ko'
   *  Constant: '<S12>/Kp'
   *  Constant: '<S12>/Ks'
   */
  /* MATLAB Function 'Parameter Module/Control Parameter/Mux1': '<S14>:1' */
  /* '<S14>:1:3' */
  Human_in_Loop_B.x_o[0] = Human_in_Loop_P.Kp_Value;
  Human_in_Loop_B.x_o[1] = Human_in_Loop_P.Kd_Value;
  Human_in_Loop_B.x_o[2] = Human_in_Loop_P.Kl_Value;
  Human_in_Loop_B.x_o[3] = Human_in_Loop_P.Ks_Value;
  Human_in_Loop_B.x_o[4] = Human_in_Loop_P.Ko_Value;

  /* RateTransition: '<S2>/RT1' */
  switch (Human_in_Loop_DW.RT1_read_buf_p) {
   case 0:
    Human_in_Loop_DW.RT1_write_buf_e = 1;
    break;

   case 1:
    Human_in_Loop_DW.RT1_write_buf_e = 0;
    break;

   default:
    Human_in_Loop_DW.RT1_write_buf_e = (int8_T)
      (Human_in_Loop_DW.RT1_last_buf_wr_a == 0);
    break;
  }

  if (Human_in_Loop_DW.RT1_write_buf_e != 0) {
    for (i = 0; i < 5; i++) {
      Human_in_Loop_DW.RT1_Buffer1_e[i] = Human_in_Loop_B.x_o[i];
    }
  } else {
    for (i = 0; i < 5; i++) {
      Human_in_Loop_DW.RT1_Buffer0_o[i] = Human_in_Loop_B.x_o[i];
    }
  }

  Human_in_Loop_DW.RT1_last_buf_wr_a = Human_in_Loop_DW.RT1_write_buf_e;
  Human_in_Loop_DW.RT1_write_buf_e = -1;

  /* End of RateTransition: '<S2>/RT1' */

  /* MATLAB Function: '<S13>/Mux1' incorporates:
   *  Constant: '<S13>/fall_time'
   *  Constant: '<S13>/peak_time'
   *  Constant: '<S13>/peak_torque'
   *  Constant: '<S13>/rise_time'
   */
  /* MATLAB Function 'Parameter Module/Torque Parameter/Mux1': '<S15>:1' */
  /* '<S15>:1:3' */
  Human_in_Loop_B.x_k[0] = Human_in_Loop_P.peak_torque_Value;
  Human_in_Loop_B.x_k[1] = Human_in_Loop_P.rise_time_Value;
  Human_in_Loop_B.x_k[2] = Human_in_Loop_P.peak_time_Value;
  Human_in_Loop_B.x_k[3] = Human_in_Loop_P.fall_time_Value;

  /* RateTransition: '<S2>/RT2' */
  switch (Human_in_Loop_DW.RT2_read_buf_k) {
   case 0:
    Human_in_Loop_DW.RT2_write_buf_i = 1;
    break;

   case 1:
    Human_in_Loop_DW.RT2_write_buf_i = 0;
    break;

   default:
    Human_in_Loop_DW.RT2_write_buf_i = (int8_T)
      (Human_in_Loop_DW.RT2_last_buf_wr_d == 0);
    break;
  }

  if (Human_in_Loop_DW.RT2_write_buf_i != 0) {
    Human_in_Loop_DW.RT2_Buffer1_i[0] = Human_in_Loop_B.x_k[0];
    Human_in_Loop_DW.RT2_Buffer1_i[1] = Human_in_Loop_B.x_k[1];
    Human_in_Loop_DW.RT2_Buffer1_i[2] = Human_in_Loop_B.x_k[2];
    Human_in_Loop_DW.RT2_Buffer1_i[3] = Human_in_Loop_B.x_k[3];
  } else {
    Human_in_Loop_DW.RT2_Buffer0_k[0] = Human_in_Loop_B.x_k[0];
    Human_in_Loop_DW.RT2_Buffer0_k[1] = Human_in_Loop_B.x_k[1];
    Human_in_Loop_DW.RT2_Buffer0_k[2] = Human_in_Loop_B.x_k[2];
    Human_in_Loop_DW.RT2_Buffer0_k[3] = Human_in_Loop_B.x_k[3];
  }

  Human_in_Loop_DW.RT2_last_buf_wr_d = Human_in_Loop_DW.RT2_write_buf_i;
  Human_in_Loop_DW.RT2_write_buf_i = -1;

  /* End of RateTransition: '<S2>/RT2' */

  /* S-Function (rti_commonblock): '<S6>/S-Function1' */

  /* This comment workarounds a code generation problem */

  /* End of Outputs for S-Function (rti_commonblock): '<S6>/S-Function1' */

  /* Gain: '<S37>/Gain' */
  Human_in_Loop_B.Gain_g = Human_in_Loop_P.uOrderTD_Ts * Human_in_Loop_B.x2k1;

  /* Sum: '<S37>/Add' */
  Human_in_Loop_B.x1k = Human_in_Loop_B.Gain_g + Human_in_Loop_B.x1k1;

  /* MATLAB Function: '<S21>/MATLAB Function' */
  /* MATLAB Function 'Sensor Data/Torque module/MATLAB Function': '<S40>:1' */
  /* '<S40>:1:13' */
  /* '<S40>:1:8' */
  memcpy(&B[0], &Human_in_Loop_DW.data[1], 14U * sizeof(real_T));
  B[14] = Human_in_Loop_B.torque;
  memcpy(&Human_in_Loop_DW.data[0], &B[0], 15U * sizeof(real_T));

  /* '<S40>:1:13' */
  for (i = 0; i < 30; i++) {
    A[i] = b_A[i];
  }

  Human_in_Loop_xgeqp3(A, tau, jpvt);
  i = 0;
  tol = 15.0 * fabs(A[0]) * 2.2204460492503131E-16;
  while ((i < 2) && (!(fabs(A[15 * i + i]) <= tol))) {
    i++;
  }

  x[0] = 0.0;
  x[1] = 0.0;
  memcpy(&B[0], &Human_in_Loop_DW.data[0], 15U * sizeof(real_T));
  if (tau[0] != 0.0) {
    tol = B[0];
    for (c_i = 1; c_i + 1 < 16; c_i++) {
      tol += A[c_i] * B[c_i];
    }

    tol *= tau[0];
    if (tol != 0.0) {
      B[0] -= tol;
      for (c_i = 1; c_i + 1 < 16; c_i++) {
        B[c_i] -= A[c_i] * tol;
      }
    }
  }

  if (tau[1] != 0.0) {
    tol = B[1];
    for (c_i = 2; c_i + 1 < 16; c_i++) {
      tol += A[c_i + 15] * B[c_i];
    }

    tol *= tau[1];
    if (tol != 0.0) {
      B[1] -= tol;
      for (c_i = 2; c_i + 1 < 16; c_i++) {
        B[c_i] -= A[c_i + 15] * tol;
      }
    }
  }

  for (c_i = 0; c_i + 1 <= i; c_i++) {
    x[jpvt[c_i] - 1] = B[c_i];
  }

  for (i--; i + 1 > 0; i--) {
    x[jpvt[i] - 1] /= A[15 * i + i];
    c_i = 1;
    while (c_i <= i) {
      x[jpvt[0] - 1] -= x[jpvt[i] - 1] * A[15 * i];
      c_i = 2;
    }
  }

  /* '<S40>:1:15' */
  Human_in_Loop_B.torque_dot = x[1] * 5000.0;

  /* End of MATLAB Function: '<S21>/MATLAB Function' */

  /* SampleTimeMath: '<S26>/TSamp'
   *
   * About '<S26>/TSamp':
   *  y = u * K where K = 1 / ( w * Ts )
   */
  Human_in_Loop_B.TSamp = Human_in_Loop_B.Gain3_m * Human_in_Loop_P.TSamp_WtEt;

  /* UnitDelay: '<S26>/UD' */
  Human_in_Loop_B.Uk1 = Human_in_Loop_DW.UD_DSTATE;

  /* Sum: '<S26>/Diff' */
  Human_in_Loop_B.Diff = Human_in_Loop_B.TSamp - Human_in_Loop_B.Uk1;

  /* DiscreteFilter: '<S19>/0.4low1' */
  tol = Human_in_Loop_B.Diff;
  tol -= Human_in_Loop_P.u4low1_DenCoef[1] * Human_in_Loop_DW.u4low1_states[0];
  tol -= Human_in_Loop_P.u4low1_DenCoef[2] * Human_in_Loop_DW.u4low1_states[1];
  tol -= Human_in_Loop_P.u4low1_DenCoef[3] * Human_in_Loop_DW.u4low1_states[2];
  tol /= Human_in_Loop_P.u4low1_DenCoef[0];
  Human_in_Loop_DW.u4low1_tmp = tol;
  tol = Human_in_Loop_P.u4low1_NumCoef[0] * Human_in_Loop_DW.u4low1_tmp;
  tol += Human_in_Loop_P.u4low1_NumCoef[1] * Human_in_Loop_DW.u4low1_states[0];
  tol += Human_in_Loop_P.u4low1_NumCoef[2] * Human_in_Loop_DW.u4low1_states[1];
  tol += Human_in_Loop_P.u4low1_NumCoef[3] * Human_in_Loop_DW.u4low1_states[2];
  Human_in_Loop_B.u4low1 = tol;

  /* UnitDelay: '<S22>/Unit Delay' */
  Human_in_Loop_B.UnitDelay = Human_in_Loop_DW.UnitDelay_DSTATE_j;

  /* UnitDelay: '<S22>/Unit Delay1' */
  Human_in_Loop_B.UnitDelay1 = Human_in_Loop_DW.UnitDelay1_DSTATE_a;

  /* Sum: '<S22>/Add1' */
  Human_in_Loop_B.Add1 = Human_in_Loop_B.UnitDelay - Human_in_Loop_B.UnitDelay1;

  /* Gain: '<S22>/Gain' */
  Human_in_Loop_B.Gain_p = Human_in_Loop_P.uOrderTD_Ts_g * Human_in_Loop_B.Add1;

  /* Gain: '<S22>/Gain1' */
  B_0 = Human_in_Loop_P.uOrderTD_T;
  tol = 1.0 / B_0;
  Human_in_Loop_B.Gain1_d = tol * Human_in_Loop_B.Gain_p;

  /* Sum: '<S22>/Add2' */
  Human_in_Loop_B.Add2_i = Human_in_Loop_B.UnitDelay1 + Human_in_Loop_B.Gain1_d;

  /* Sum: '<S22>/Add3' */
  Human_in_Loop_B.Add3 = Human_in_Loop_B.Gain3_m - Human_in_Loop_B.Add2_i;

  /* Gain: '<S22>/Gain2' */
  B_0 = Human_in_Loop_P.uOrderTD_T;
  tol = 1.0 / B_0;
  Human_in_Loop_B.Gain2_i = tol * Human_in_Loop_B.Add3;

  /* Gain: '<S23>/Gain' */
  Human_in_Loop_B.Gain_i = Human_in_Loop_P.uOrderTD_Ts_n *
    Human_in_Loop_B.x2k1_k;

  /* Sum: '<S23>/Add' */
  Human_in_Loop_B.x1k_i = Human_in_Loop_B.Gain_i + Human_in_Loop_B.x1k1_m;

  /* MATLAB Function: '<S19>/Mux2' */
  /* MATLAB Function 'Sensor Data/Encoder module/Mux2': '<S33>:1' */
  /* '<S33>:1:3' */
  Human_in_Loop_B.x_d[0] = Human_in_Loop_B.u4low1;
  Human_in_Loop_B.x_d[1] = Human_in_Loop_B.Gain2_i;
  Human_in_Loop_B.x_d[2] = Human_in_Loop_B.x2k_i;

  /* S-Function (rti_commonblock): '<S29>/S-Function1' incorporates:
   *  Constant: '<S19>/VCC1'
   */
  /* This comment workarounds a code generation problem */

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup1 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 11 --- */
  {
    /* define variables required for BitOut realtime functions */
    UInt32 outputData = 0;

    /* write output state value to digital output channel 11-11 on port 3 */
    outputData = ( ( ( (UInt32)Human_in_Loop_P.VCC1_Value) << (11 - 1)) |
                  outputData);
    DioCl1DigOut_setChMaskOutData(pRTIDioC1DigOut_Port_3_Ch_11, outputData);
    DioCl1DigOut_write(pRTIDioC1DigOut_Port_3_Ch_11);
  }

  /* S-Function (rti_commonblock): '<S30>/S-Function1' incorporates:
   *  Constant: '<S19>/VCC3'
   */
  /* This comment workarounds a code generation problem */

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup3 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 13 --- */
  {
    /* define variables required for BitOut realtime functions */
    UInt32 outputData = 0;

    /* write output state value to digital output channel 13-13 on port 3 */
    outputData = ( ( ( (UInt32)Human_in_Loop_P.VCC3_Value) << (13 - 1)) |
                  outputData);
    DioCl1DigOut_setChMaskOutData(pRTIDioC1DigOut_Port_3_Ch_13, outputData);
    DioCl1DigOut_write(pRTIDioC1DigOut_Port_3_Ch_13);
  }

  /* S-Function (rti_commonblock): '<S35>/S-Function1' incorporates:
   *  Constant: '<S20>/Constant'
   */
  /* This comment workarounds a code generation problem */

  /* --- Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_OUT_BL1 --- */
  /* --- [RTI120X, BITOUT] - Port: 1 - Channel: 1 --- */
  {
    /* define variables required for BitOut realtime functions */
    UInt32 outputData = 0;

    /* write output state value to digital output channel 1-1 on port 1 */
    outputData = ( ( ( (UInt32)Human_in_Loop_P.Constant_Value) << (1 - 1)) |
                  outputData);
    DioCl1DigOut_setChMaskOutData(pRTIDioC1DigOut_Port_1_Ch_1, outputData);
    DioCl1DigOut_write(pRTIDioC1DigOut_Port_1_Ch_1);
  }
}

/* Model update function */
void Human_in_Loop_update(void)
{
  /* Update for DiscreteFilter: '<S21>/0.4low2' */
  Human_in_Loop_DW.u4low2_states[2] = Human_in_Loop_DW.u4low2_states[1];
  Human_in_Loop_DW.u4low2_states[1] = Human_in_Loop_DW.u4low2_states[0];
  Human_in_Loop_DW.u4low2_states[0] = Human_in_Loop_DW.u4low2_tmp;

  /* Update for UnitDelay: '<S37>/Unit Delay1' */
  Human_in_Loop_DW.UnitDelay1_DSTATE = Human_in_Loop_B.x2k;

  /* Update for UnitDelay: '<S37>/Unit Delay' */
  Human_in_Loop_DW.UnitDelay_DSTATE = Human_in_Loop_B.x1k;

  /* Update for UnitDelay: '<S37>/Unit Delay2' */
  Human_in_Loop_DW.UnitDelay2_DSTATE = Human_in_Loop_B.torque;

  /* Update for UnitDelay: '<S23>/Unit Delay1' */
  Human_in_Loop_DW.UnitDelay1_DSTATE_h = Human_in_Loop_B.x2k_i;

  /* Update for UnitDelay: '<S23>/Unit Delay' */
  Human_in_Loop_DW.UnitDelay_DSTATE_e = Human_in_Loop_B.x1k_i;

  /* Update for UnitDelay: '<S23>/Unit Delay2' */
  Human_in_Loop_DW.UnitDelay2_DSTATE_a = Human_in_Loop_B.Gain3_m;

  /* Update for UnitDelay: '<S26>/UD' */
  Human_in_Loop_DW.UD_DSTATE = Human_in_Loop_B.TSamp;

  /* Update for DiscreteFilter: '<S19>/0.4low1' */
  Human_in_Loop_DW.u4low1_states[2] = Human_in_Loop_DW.u4low1_states[1];
  Human_in_Loop_DW.u4low1_states[1] = Human_in_Loop_DW.u4low1_states[0];
  Human_in_Loop_DW.u4low1_states[0] = Human_in_Loop_DW.u4low1_tmp;

  /* Update for UnitDelay: '<S22>/Unit Delay' */
  Human_in_Loop_DW.UnitDelay_DSTATE_j = Human_in_Loop_B.Gain3_m;

  /* Update for UnitDelay: '<S22>/Unit Delay1' */
  Human_in_Loop_DW.UnitDelay1_DSTATE_a = Human_in_Loop_B.Add2_i;

  /* Update absolute time for base rate */
  /* The "clockTick0" counts the number of times the code of this task has
   * been executed. The absolute time is the multiplication of "clockTick0"
   * and "Timing.stepSize0". Size of "clockTick0" ensures timer will not
   * overflow during the application lifespan selected.
   * Timer of this task consists of two 32 bit unsigned integers.
   * The two integers represent the low bits Timing.clockTick0 and the high bits
   * Timing.clockTickH0. When the low bit overflows to 0, the high bits increment.
   */
  if (!(++Human_in_Loop_M->Timing.clockTick0)) {
    ++Human_in_Loop_M->Timing.clockTickH0;
  }

  Human_in_Loop_M->Timing.taskTime0 = Human_in_Loop_M->Timing.clockTick0 *
    Human_in_Loop_M->Timing.stepSize0 + Human_in_Loop_M->Timing.clockTickH0 *
    Human_in_Loop_M->Timing.stepSize0 * 4294967296.0;
}

/* Model initialize function */
void Human_in_Loop_initialize(void)
{
  /* Registration code */

  /* initialize non-finites */
  rt_InitInfAndNaN(sizeof(real_T));

  /* initialize real-time model */
  (void) memset((void *)Human_in_Loop_M, 0,
                sizeof(RT_MODEL_Human_in_Loop_T));
  Human_in_Loop_M->Timing.stepSize0 = 0.0002;

  /* block I/O */
  (void) memset(((void *) &Human_in_Loop_B), 0,
                sizeof(B_Human_in_Loop_T));

  /* states (dwork) */
  (void) memset((void *)&Human_in_Loop_DW, 0,
                sizeof(DW_Human_in_Loop_T));

  {
    /* user code (registration function declaration) */
    /*Initialize global TRC pointers. */
    Human_in_Loop_rti_init_trc_pointers();
  }

  {
    int32_T i;

    /* Start for RateTransition: '<S4>/RT1' */
    Human_in_Loop_B.RT1[0] = Human_in_Loop_P.RT1_InitialCondition;
    Human_in_Loop_B.RT1[1] = Human_in_Loop_P.RT1_InitialCondition;

    /* Start for RateTransition: '<S4>/RT2' */
    Human_in_Loop_B.RT2[0] = Human_in_Loop_P.RT2_InitialCondition;
    Human_in_Loop_B.RT2[1] = Human_in_Loop_P.RT2_InitialCondition;

    /* Start for RateTransition: '<S4>/RT3' */
    Human_in_Loop_B.RT3[0] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_B.RT3[1] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_B.RT3[2] = Human_in_Loop_P.RT3_InitialCondition;

    /* Start for RateTransition: '<S5>/RT1' */
    Human_in_Loop_B.RT1_n[0] = Human_in_Loop_P.RT1_InitialCondition_i;
    Human_in_Loop_B.RT1_n[1] = Human_in_Loop_P.RT1_InitialCondition_i;
    Human_in_Loop_B.RT1_n[2] = Human_in_Loop_P.RT1_InitialCondition_i;
    Human_in_Loop_B.RT1_n[3] = Human_in_Loop_P.RT1_InitialCondition_i;

    /* Start for RateTransition: '<S2>/RT1' */
    for (i = 0; i < 5; i++) {
      Human_in_Loop_B.RT1_a[i] = Human_in_Loop_P.RT1_InitialCondition_e;
    }

    /* End of Start for RateTransition: '<S2>/RT1' */

    /* Start for RateTransition: '<S2>/RT2' */
    Human_in_Loop_B.RT2_h[0] = Human_in_Loop_P.RT2_InitialCondition_e;
    Human_in_Loop_B.RT2_h[1] = Human_in_Loop_P.RT2_InitialCondition_e;
    Human_in_Loop_B.RT2_h[2] = Human_in_Loop_P.RT2_InitialCondition_e;
    Human_in_Loop_B.RT2_h[3] = Human_in_Loop_P.RT2_InitialCondition_e;

    /* Start for DataStoreMemory: '<Root>/Data Store Memory' */
    memcpy(&Human_in_Loop_DW.TorqueMem[0],
           &Human_in_Loop_P.DataStoreMemory_InitialValue[0], 4400U * sizeof
           (real_T));
  }

  {
    int32_T i;

    /* InitializeConditions for DiscreteFilter: '<S21>/0.4low2' */
    Human_in_Loop_DW.u4low2_states[0] = Human_in_Loop_P.u4low2_InitialStates;
    Human_in_Loop_DW.u4low2_states[1] = Human_in_Loop_P.u4low2_InitialStates;
    Human_in_Loop_DW.u4low2_states[2] = Human_in_Loop_P.u4low2_InitialStates;

    /* InitializeConditions for UnitDelay: '<S37>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE =
      Human_in_Loop_P.UnitDelay1_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S37>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE =
      Human_in_Loop_P.UnitDelay_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S37>/Unit Delay2' */
    Human_in_Loop_DW.UnitDelay2_DSTATE =
      Human_in_Loop_P.UnitDelay2_InitialCondition;

    /* InitializeConditions for RateTransition: '<S4>/RT1' */
    Human_in_Loop_DW.RT1_Buffer0[0] = Human_in_Loop_P.RT1_InitialCondition;
    Human_in_Loop_DW.RT1_Buffer0[1] = Human_in_Loop_P.RT1_InitialCondition;
    Human_in_Loop_DW.RT1_write_buf = -1;
    Human_in_Loop_DW.RT1_read_buf = -1;

    /* InitializeConditions for RateTransition: '<S4>/RT2' */
    Human_in_Loop_DW.RT2_Buffer0[0] = Human_in_Loop_P.RT2_InitialCondition;
    Human_in_Loop_DW.RT2_Buffer0[1] = Human_in_Loop_P.RT2_InitialCondition;
    Human_in_Loop_DW.RT2_write_buf = -1;
    Human_in_Loop_DW.RT2_read_buf = -1;

    /* InitializeConditions for UnitDelay: '<S23>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE_h =
      Human_in_Loop_P.UnitDelay1_InitialCondition_j;

    /* InitializeConditions for UnitDelay: '<S23>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_e =
      Human_in_Loop_P.UnitDelay_InitialCondition_j;

    /* InitializeConditions for UnitDelay: '<S23>/Unit Delay2' */
    Human_in_Loop_DW.UnitDelay2_DSTATE_a =
      Human_in_Loop_P.UnitDelay2_InitialCondition_h;

    /* InitializeConditions for RateTransition: '<S4>/RT3' */
    Human_in_Loop_DW.RT3_Buffer0[0] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_DW.RT3_Buffer0[1] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_DW.RT3_Buffer0[2] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_DW.RT3_write_buf = -1;
    Human_in_Loop_DW.RT3_read_buf = -1;

    /* InitializeConditions for RateTransition: '<S5>/RT1' */
    Human_in_Loop_DW.RT1_Buffer0_a[0] = Human_in_Loop_P.RT1_InitialCondition_i;
    Human_in_Loop_DW.RT1_Buffer0_a[1] = Human_in_Loop_P.RT1_InitialCondition_i;
    Human_in_Loop_DW.RT1_Buffer0_a[2] = Human_in_Loop_P.RT1_InitialCondition_i;
    Human_in_Loop_DW.RT1_Buffer0_a[3] = Human_in_Loop_P.RT1_InitialCondition_i;
    Human_in_Loop_DW.RT1_write_buf_n = -1;
    Human_in_Loop_DW.RT1_read_buf_l = -1;

    /* InitializeConditions for RateTransition: '<S2>/RT1' */
    for (i = 0; i < 5; i++) {
      Human_in_Loop_DW.RT1_Buffer0_o[i] = Human_in_Loop_P.RT1_InitialCondition_e;
    }

    Human_in_Loop_DW.RT1_write_buf_e = -1;
    Human_in_Loop_DW.RT1_read_buf_p = -1;

    /* End of InitializeConditions for RateTransition: '<S2>/RT1' */

    /* InitializeConditions for RateTransition: '<S2>/RT2' */
    Human_in_Loop_DW.RT2_Buffer0_k[0] = Human_in_Loop_P.RT2_InitialCondition_e;
    Human_in_Loop_DW.RT2_Buffer0_k[1] = Human_in_Loop_P.RT2_InitialCondition_e;
    Human_in_Loop_DW.RT2_Buffer0_k[2] = Human_in_Loop_P.RT2_InitialCondition_e;
    Human_in_Loop_DW.RT2_Buffer0_k[3] = Human_in_Loop_P.RT2_InitialCondition_e;
    Human_in_Loop_DW.RT2_write_buf_i = -1;
    Human_in_Loop_DW.RT2_read_buf_k = -1;

    /* InitializeConditions for UnitDelay: '<S26>/UD' */
    Human_in_Loop_DW.UD_DSTATE = Human_in_Loop_P.DiscreteDerivative_ICPrevScaled;

    /* InitializeConditions for DiscreteFilter: '<S19>/0.4low1' */
    Human_in_Loop_DW.u4low1_states[0] = Human_in_Loop_P.u4low1_InitialStates;
    Human_in_Loop_DW.u4low1_states[1] = Human_in_Loop_P.u4low1_InitialStates;
    Human_in_Loop_DW.u4low1_states[2] = Human_in_Loop_P.u4low1_InitialStates;

    /* InitializeConditions for UnitDelay: '<S22>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_j =
      Human_in_Loop_P.UnitDelay_InitialCondition_g;

    /* InitializeConditions for UnitDelay: '<S22>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE_a =
      Human_in_Loop_P.UnitDelay1_InitialCondition_i;

    /* SystemInitialize for MATLAB Function: '<S21>/Data process' */
    Human_in_Loop_DW.torque_zero = 0.0;

    /* SystemInitialize for MATLAB Function: '<S19>/Data process' */
    Human_in_Loop_DW.angle_zero_f = 0.0;

    /* SystemInitialize for MATLAB Function: '<S19>/Data process1' */
    Human_in_Loop_DW.angle_zero = 0.0;

    /* SystemInitialize for MATLAB Function: '<S20>/FootSwitch Filter' */
    Human_in_Loop_DW.foot_state = 0.0;
    Human_in_Loop_DW.filter_time = 0.0;

    /* SystemInitialize for MATLAB Function: '<S5>/State Machine' */
    Human_in_Loop_DW.reg_stride_time = 1.0;
    Human_in_Loop_DW.reg_stride_time_count = 0.0;
    Human_in_Loop_DW.reg_mode = 1.0;
    Human_in_Loop_DW.reg_state = 1.0;
    Human_in_Loop_DW.bt_run = 0.0;
    Human_in_Loop_DW.reg_last_switch = 1.0;

    /* SystemInitialize for S-Function (rti_commonblock): '<S6>/S-Function1' incorporates:
     *  SubSystem: '<Root>/Control Module'
     */
    Human_in_Loo_ControlModule_Init();

    /* End of SystemInitialize for S-Function (rti_commonblock): '<S6>/S-Function1' */

    /* SystemInitialize for MATLAB Function: '<S21>/MATLAB Function' */
    for (i = 0; i < 15; i++) {
      Human_in_Loop_DW.data[i] = 1.0;
    }

    /* End of SystemInitialize for MATLAB Function: '<S21>/MATLAB Function' */
  }
}

/* Model terminate function */
void Human_in_Loop_terminate(void)
{
  /* Terminate for S-Function (rti_commonblock): '<S27>/S-Function1' */

  /* --- Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL1 --- */
  /* --- [RTIEMC, Encoder] - DIO class: 2 - Unit: 1 - Port: 1 - Channel: 1 --- */
  {
    /* Deactivates encoder interface functionality */
    DioCl2EncoderIn_stop(pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1);
  }

  /* Terminate for S-Function (rti_commonblock): '<S28>/S-Function1' */

  /* --- Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL3 --- */
  /* --- [RTIEMC, Encoder] - DIO class: 2 - Unit: 3 - Port: 1 - Channel: 5 --- */
  {
    /* Deactivates encoder interface functionality */
    DioCl2EncoderIn_stop(pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5);
  }

  /* Terminate for S-Function (rti_commonblock): '<S6>/S-Function1' incorporates:
   *  SubSystem: '<Root>/Control Module'
   */
  Human_in_Loo_ControlModule_Term();

  /* End of Terminate for S-Function (rti_commonblock): '<S6>/S-Function1' */

  /* Terminate for S-Function (rti_commonblock): '<S29>/S-Function1' incorporates:
   *  Constant: '<S19>/VCC1'
   */

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup1 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 11 --- */

  /* disable digital output channel 11-11 on port 3 *
   * (set to high-impedance), when the simulation terminates       */
  DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_3_Ch_11,
    DIO_CLASS1_HIGH_Z_ON);
  DioCl1DigOut_write(pRTIDioC1DigOut_Port_3_Ch_11);

  /* Terminate for S-Function (rti_commonblock): '<S30>/S-Function1' incorporates:
   *  Constant: '<S19>/VCC3'
   */

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup3 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 13 --- */

  /* disable digital output channel 13-13 on port 3 *
   * (set to high-impedance), when the simulation terminates       */
  DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_3_Ch_13,
    DIO_CLASS1_HIGH_Z_ON);
  DioCl1DigOut_write(pRTIDioC1DigOut_Port_3_Ch_13);

  /* Terminate for S-Function (rti_commonblock): '<S35>/S-Function1' incorporates:
   *  Constant: '<S20>/Constant'
   */

  /* --- Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_OUT_BL1 --- */
  /* --- [RTI120X, BITOUT] - Port: 1 - Channel: 1 --- */

  /* disable digital output channel 1-1 on port 1 *
   * (set to high-impedance), when the simulation terminates       */
  DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_1_Ch_1,
    DIO_CLASS1_HIGH_Z_ON);
  DioCl1DigOut_write(pRTIDioC1DigOut_Port_1_Ch_1);
}
