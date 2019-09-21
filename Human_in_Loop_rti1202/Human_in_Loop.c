/*
 * Human_in_Loop.c
 *
 * Code generation for model "Human_in_Loop".
 *
 * Model version              : 1.1241
 * Simulink Coder version : 8.13 (R2017b) 24-Jul-2017
 * C source code generated on : Sat Sep 21 19:43:10 2019
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

/* Continuous states */
X_Human_in_Loop_T Human_in_Loop_X;

/* Block states (auto storage) */
DW_Human_in_Loop_T Human_in_Loop_DW;

/* Previous zero-crossings (trigger) states */
PrevZCX_Human_in_Loop_T Human_in_Loop_PrevZCX;

/* Real-time model */
RT_MODEL_Human_in_Loop_T Human_in_Loop_M_;
RT_MODEL_Human_in_Loop_T *const Human_in_Loop_M = &Human_in_Loop_M_;

/* Forward declaration for local functions */
static void Human_in_Loop_mldivide(const real_T A[16], real_T B[4]);
static void Human_in_Loop_power(const real_T a_data[], const int32_T a_size[2],
  real_T y_data[], int32_T y_size[2]);
static void Human_in_Loop_power_j(const real_T a_data[], const int32_T a_size[2],
  real_T y_data[], int32_T y_size[2]);

/* Forward declaration for local functions */
static real_T Human_in_Loop_norm(const real_T x_data[], const int32_T x_size[2]);
static real_T Human_in_Loop_mean(const real_T x[10]);
static real_T Human_in_Loop_xnrm2(const real_T x[30], int32_T ix0);
static real_T Human_in_Loop_xnrm2_o(int32_T n, const real_T x[30], int32_T ix0);
static void Human_in_Loop_xgeqp3(real_T A[30], real_T tau[2], int32_T jpvt[2]);
static void rate_scheduler(void);

/*
 *   This function updates active task flag for each subrate.
 * The function is called at model base rate, hence the
 * generated code self-manages all its subrates.
 */
static void rate_scheduler(void)
{
  /* Compute which subrates run during the next base time step.  Subrates
   * are an integer multiple of the base rate counter.  Therefore, the subtask
   * counter is reset when it reaches its limit (zero means run).
   */
  (Human_in_Loop_M->Timing.TaskCounters.TID[2])++;
  if ((Human_in_Loop_M->Timing.TaskCounters.TID[2]) > 9) {/* Sample time: [0.002s, 0.0s] */
    Human_in_Loop_M->Timing.TaskCounters.TID[2] = 0;
  }
}

/*
 * This function updates continuous states using the ODE1 fixed-step
 * solver algorithm
 */
static void rt_ertODEUpdateContinuousStates(RTWSolverInfo *si )
{
  time_T tnew = rtsiGetSolverStopTime(si);
  time_T h = rtsiGetStepSize(si);
  real_T *x = rtsiGetContStates(si);
  ODE1_IntgData *id = (ODE1_IntgData *)rtsiGetSolverData(si);
  real_T *f0 = id->f[0];
  int_T i;
  int_T nXc = 24;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  Human_in_Loop_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/*
 * Output and update for atomic system:
 *    '<S29>/Mux'
 *    '<S32>/Mux'
 */
void Human_in_Loop_Mux(real_T rtu_x1, real_T rtu_x2, B_Mux_Human_in_Loop_T
  *localB)
{
  /* MATLAB Function 'Sensor Data/Encoder module/Mux': '<S60>:1' */
  /* '<S60>:1:3' */
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
static void Human_in_Loop_power_j(const real_T a_data[], const int32_T a_size[2],
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

real_T rt_roundd_snf(real_T u)
{
  real_T y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

/* System initialize for function-call system: '<Root>/Control Module' */
void Human_in_Loo_ControlModule_Init(void)
{
  /* SystemInitialize for MATLAB Function: '<S1>/Torque track' */
  Human_in_Loop_DW.last_footstate = 0.0;

  /* SystemInitialize for MATLAB Function: '<S1>/LRN' */
  Human_in_Loop_DW.last_footstate_not_empty = false;
  Human_in_Loop_DW.last_footstate_a = 0.0;
  memset(&Human_in_Loop_DW.torque_error_memory[0], 0, 1000U * sizeof(real_T));
  memset(&Human_in_Loop_DW.lrn_cmd_memory[0], 0, 1000U * sizeof(real_T));

  /* SystemInitialize for MATLAB Function: '<S1>/Controller' */
  Human_in_Loop_DW.calib_state = 0.0;
}

/* System reset for function-call system: '<Root>/Control Module' */
void Human_in_Lo_ControlModule_Reset(void)
{
  /* SystemReset for MATLAB Function: '<S1>/Torque track' */
  Human_in_Loop_DW.last_footstate = 0.0;

  /* SystemReset for MATLAB Function: '<S1>/LRN' */
  Human_in_Loop_DW.last_footstate_not_empty = false;
  Human_in_Loop_DW.last_footstate_a = 0.0;
  memset(&Human_in_Loop_DW.torque_error_memory[0], 0, 1000U * sizeof(real_T));
  memset(&Human_in_Loop_DW.lrn_cmd_memory[0], 0, 1000U * sizeof(real_T));

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
  int32_T b;
  real_T c;
  int32_T e;
  int32_T g;
  int32_T i;
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

  /* MATLAB Function: '<S1>/Torque track' incorporates:
   *  Constant: '<S1>/torque_offset'
   */
  /* MATLAB Function 'Control Module/Torque track': '<S12>:1' */
  /* '<S12>:1:43' */
  /* '<S12>:1:19' */
  mode = Human_in_Loop_B.RT1[0];

  /* '<S12>:1:20' */
  footstate = Human_in_Loop_B.RT1[1];

  /* '<S12>:1:21' */
  /* '<S12>:1:22' */
  /* '<S12>:1:24' */
  peak_torque = Human_in_Loop_B.RT3[0];

  /* '<S12>:1:25' */
  /* '<S12>:1:26' */
  /* '<S12>:1:27' */
  /* '<S12>:1:29' */
  torque_measure = Human_in_Loop_B.RT4[0];

  /* '<S12>:1:30' */
  troque_delta = Human_in_Loop_B.RT4[1];

  /* '<S12>:1:31' */
  stride_index = Human_in_Loop_B.RT1[3] * 500.0 + 1.0;
  if (stride_index > 750.0) {
    /* '<S12>:1:32' */
    /* '<S12>:1:33' */
    stride_index = 750.0;
  }

  if ((Human_in_Loop_DW.last_footstate == 0.0) && (Human_in_Loop_B.RT1[1] == 1.0)
      && ((Human_in_Loop_B.RT1[0] == 2.0) || (Human_in_Loop_B.RT1[0] == 1.0))) {
    /* '<S12>:1:37' */
    /* '<S12>:1:39' */
    /* '<S12>:1:41' */
    memset(&Human_in_Loop_B.torque_track[0], 0, 750U * sizeof(real_T));
    memset(&Human_in_Loop_B.torque_delta_track[0], 0, 750U * sizeof(real_T));

    /* '<S12>:1:42' */
    /* '<S12>:1:43' */
    /* '<S12>:1:44' */
    index_peak = floor(Human_in_Loop_B.RT3[2] / 100.0 * Human_in_Loop_B.RT1[2] *
                       500.0);

    /* '<S12>:1:45' */
    index_rise = index_peak - floor(Human_in_Loop_B.RT3[1] / 100.0 *
      Human_in_Loop_B.RT1[2] * 500.0);

    /* '<S12>:1:46' */
    index_fall = floor(Human_in_Loop_B.RT3[3] / 100.0 * Human_in_Loop_B.RT1[2] *
                       500.0) + index_peak;
    if (1.0 > index_rise - 1.0) {
      b = 0;
    } else {
      b = (int32_T)(index_rise - 1.0);
    }

    /* '<S12>:1:49' */
    for (g = 0; g < b; g++) {
      Human_in_Loop_B.tmp_data_cx[g] = (real_T)(int16_T)(1 + (int16_T)g) /
        index_rise * Human_in_Loop_P.torque_offset_Value;
    }

    if (0 <= b - 1) {
      memcpy(&Human_in_Loop_B.torque_track[0], &Human_in_Loop_B.tmp_data_cx[0],
             b * sizeof(real_T));
    }

    /* '<S12>:1:53' */
    /* '<S12>:1:57' */
    /* '<S12>:1:58' */
    parm1[0] = Human_in_Loop_P.torque_offset_Value;
    parm1[1] = Human_in_Loop_B.RT3[0];
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
    c = index_peak - index_rise;
    if (1.0 > c) {
      b_n = 0;
    } else {
      b_n = (int32_T)c;
    }

    if (index_rise > index_peak - 1.0) {
      br = 1;
      e = 0;
      ar = 0;
      g = 0;
      ia = 0;
      i = 0;
      b = 0;
    } else {
      br = (int32_T)index_rise;
      e = (int32_T)(index_peak - 1.0);
      ar = (int32_T)index_rise - 1;
      g = (int32_T)(index_peak - 1.0);
      ia = (int32_T)index_rise - 1;
      i = (int32_T)(index_peak - 1.0);
      b = (int32_T)index_rise - 1;
    }

    tmp_size_idx_1 = g - ar;
    loop_ub = g - ar;
    for (g = 0; g < loop_ub; g++) {
      tmp_data[g] = (int16_T)((int16_T)(ar + g) + 1);
    }

    tmp_size_2[0] = 1;
    tmp_size_2[1] = tmp_size_idx_1;
    for (g = 0; g < tmp_size_idx_1; g++) {
      Human_in_Loop_B.tmp_data_cx[g] = tmp_data[g];
    }

    Human_in_Loop_power(Human_in_Loop_B.tmp_data_cx, tmp_size_2,
                        Human_in_Loop_B.tmp_data_c, tmp_size_3);
    tmp_size_idx_1 = i - ia;
    loop_ub = i - ia;
    for (g = 0; g < loop_ub; g++) {
      tmp_data[g] = (int16_T)((int16_T)(ia + g) + 1);
    }

    tmp_size_4[0] = 1;
    tmp_size_4[1] = tmp_size_idx_1;
    for (g = 0; g < tmp_size_idx_1; g++) {
      Human_in_Loop_B.tmp_data_cx[g] = tmp_data[g];
    }

    Human_in_Loop_power_j(Human_in_Loop_B.tmp_data_cx, tmp_size_4,
                          Human_in_Loop_B.tmp_data_k, tmp_size_2);
    for (g = 0; g < b_n; g++) {
      Human_in_Loop_B.tmp_data[g << 2] = 1.0;
    }

    loop_ub = e - br;
    for (g = 0; g <= loop_ub; g++) {
      Human_in_Loop_B.tmp_data[1 + (g << 2)] = (int16_T)((int16_T)((br + g) - 1)
        + 1);
    }

    loop_ub = tmp_size_3[1];
    for (g = 0; g < loop_ub; g++) {
      Human_in_Loop_B.tmp_data[2 + (g << 2)] =
        Human_in_Loop_B.tmp_data_c[tmp_size_3[0] * g];
    }

    loop_ub = tmp_size_2[1];
    for (g = 0; g < loop_ub; g++) {
      Human_in_Loop_B.tmp_data[3 + (g << 2)] =
        Human_in_Loop_B.tmp_data_k[tmp_size_2[0] * g];
    }

    br = b_n;
    for (g = 0; g < b_n; g++) {
      Human_in_Loop_B.result_data[g << 2] = Human_in_Loop_B.tmp_data[g << 2];
      Human_in_Loop_B.result_data[1 + (g << 2)] = Human_in_Loop_B.tmp_data[(g <<
        2) + 1];
      Human_in_Loop_B.result_data[2 + (g << 2)] = Human_in_Loop_B.tmp_data[(g <<
        2) + 2];
      Human_in_Loop_B.result_data[3 + (g << 2)] = Human_in_Loop_B.tmp_data[(g <<
        2) + 3];
    }

    tmp_0 = (int16_T)br;
    tmp_size_idx_1 = tmp_0;
    b_n = br - 1;
    if (0 <= tmp_size_idx_1 - 1) {
      memset(&Human_in_Loop_B.tmp_data_cx[0], 0, tmp_size_idx_1 * sizeof(real_T));
    }

    if (br != 0) {
      for (br = 1; br - 1 <= b_n; br++) {
        for (e = br; e <= br; e++) {
          Human_in_Loop_B.tmp_data_cx[e - 1] = 0.0;
        }
      }

      br = 0;
      for (e = 0; e <= b_n; e++) {
        ar = -1;
        for (g = br; g + 1 <= br + 4; g++) {
          if (Human_in_Loop_B.result_data[g] != 0.0) {
            ia = ar;
            for (i = e; i + 1 <= e + 1; i++) {
              ia++;
              Human_in_Loop_B.tmp_data_cx[i] += Human_in_Loop_B.result_data[g] *
                parm1[ia];
            }
          }

          ar++;
        }

        br += 4;
      }
    }

    /* '<S12>:1:59' */
    for (g = 0; g < tmp_size_idx_1; g++) {
      Human_in_Loop_B.torque_track[b + g] = Human_in_Loop_B.tmp_data_cx[g];
    }

    c = index_peak - index_rise;
    if (1.0 > c) {
      b_n = 0;
    } else {
      b_n = (int32_T)c;
    }

    if (index_rise > index_peak - 1.0) {
      br = 1;
      e = 0;
      ar = 0;
      g = 0;
      b = 0;
    } else {
      br = (int32_T)index_rise;
      e = (int32_T)(index_peak - 1.0);
      ar = (int32_T)index_rise - 1;
      g = (int32_T)(index_peak - 1.0);
      b = (int32_T)index_rise - 1;
    }

    tmp_size_idx_1 = g - ar;
    loop_ub = g - ar;
    for (g = 0; g < loop_ub; g++) {
      tmp_data[g] = (int16_T)((int16_T)(ar + g) + 1);
    }

    tmp_size_1[0] = 1;
    tmp_size_1[1] = tmp_size_idx_1;
    for (g = 0; g < tmp_size_idx_1; g++) {
      Human_in_Loop_B.tmp_data_cx[g] = tmp_data[g];
    }

    Human_in_Loop_power(Human_in_Loop_B.tmp_data_cx, tmp_size_1,
                        Human_in_Loop_B.tmp_data_c, tmp_size_2);
    for (g = 0; g < b_n; g++) {
      Human_in_Loop_B.tmp_data_m[3 * g] = 1.0;
    }

    loop_ub = e - br;
    for (g = 0; g <= loop_ub; g++) {
      Human_in_Loop_B.tmp_data_m[1 + 3 * g] = (real_T)(int16_T)((int16_T)((br +
        g) - 1) + 1) * 2.0;
    }

    loop_ub = tmp_size_2[1];
    for (g = 0; g < loop_ub; g++) {
      Human_in_Loop_B.tmp_data_m[2 + 3 * g] =
        Human_in_Loop_B.tmp_data_c[tmp_size_2[0] * g] * 3.0;
    }

    br = b_n;
    for (g = 0; g < b_n; g++) {
      Human_in_Loop_B.b_result_data[3 * g] = Human_in_Loop_B.tmp_data_m[3 * g];
      Human_in_Loop_B.b_result_data[1 + 3 * g] = Human_in_Loop_B.tmp_data_m[3 *
        g + 1];
      Human_in_Loop_B.b_result_data[2 + 3 * g] = Human_in_Loop_B.tmp_data_m[3 *
        g + 2];
    }

    tmp_0 = (int16_T)br;
    tmp_size_idx_1 = tmp_0;
    b_n = br - 1;
    if (0 <= tmp_size_idx_1 - 1) {
      memset(&Human_in_Loop_B.tmp_data_cx[0], 0, tmp_size_idx_1 * sizeof(real_T));
    }

    if (br != 0) {
      for (br = 1; br - 1 <= b_n; br++) {
        for (e = br; e <= br; e++) {
          Human_in_Loop_B.tmp_data_cx[e - 1] = 0.0;
        }
      }

      br = 0;
      for (e = 0; e <= b_n; e++) {
        ar = -1;
        for (g = br; g + 1 <= br + 3; g++) {
          if (Human_in_Loop_B.b_result_data[g] != 0.0) {
            ia = ar;
            for (i = e; i + 1 <= e + 1; i++) {
              ia++;
              Human_in_Loop_B.tmp_data_cx[i] += parm1[1 + ia] *
                Human_in_Loop_B.b_result_data[g];
            }
          }

          ar++;
        }

        br += 3;
      }
    }

    /* '<S12>:1:60' */
    for (g = 0; g < tmp_size_idx_1; g++) {
      Human_in_Loop_B.torque_delta_track[b + g] = Human_in_Loop_B.tmp_data_cx[g]
        * 500.0;
    }

    /* '<S12>:1:63' */
    /* '<S12>:1:67' */
    /* '<S12>:1:68' */
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
      e = 0;
      ar = 0;
      g = 0;
      ia = 0;
      i = 0;
      b = 0;
    } else {
      br = (int32_T)index_peak;
      e = (int32_T)index_fall;
      ar = (int32_T)index_peak - 1;
      g = (int32_T)index_fall;
      ia = (int32_T)index_peak - 1;
      i = (int32_T)index_fall;
      b = (int32_T)index_peak - 1;
    }

    tmp_size_idx_1 = g - ar;
    loop_ub = g - ar;
    for (g = 0; g < loop_ub; g++) {
      tmp_data[g] = (int16_T)((int16_T)(ar + g) + 1);
    }

    tmp_size[0] = 1;
    tmp_size[1] = tmp_size_idx_1;
    for (g = 0; g < tmp_size_idx_1; g++) {
      Human_in_Loop_B.tmp_data_cx[g] = tmp_data[g];
    }

    Human_in_Loop_power(Human_in_Loop_B.tmp_data_cx, tmp_size,
                        Human_in_Loop_B.tmp_data_c, tmp_size_2);
    tmp_size_idx_1 = i - ia;
    loop_ub = i - ia;
    for (g = 0; g < loop_ub; g++) {
      tmp_data[g] = (int16_T)((int16_T)(ia + g) + 1);
    }

    tmp_size_0[0] = 1;
    tmp_size_0[1] = tmp_size_idx_1;
    for (g = 0; g < tmp_size_idx_1; g++) {
      Human_in_Loop_B.tmp_data_cx[g] = tmp_data[g];
    }

    Human_in_Loop_power_j(Human_in_Loop_B.tmp_data_cx, tmp_size_0,
                          Human_in_Loop_B.tmp_data_k, tmp_size_3);
    for (g = 0; g < b_n; g++) {
      Human_in_Loop_B.tmp_data[g << 2] = 1.0;
    }

    loop_ub = e - br;
    for (g = 0; g <= loop_ub; g++) {
      Human_in_Loop_B.tmp_data[1 + (g << 2)] = (int16_T)((int16_T)((br + g) - 1)
        + 1);
    }

    loop_ub = tmp_size_2[1];
    for (g = 0; g < loop_ub; g++) {
      Human_in_Loop_B.tmp_data[2 + (g << 2)] =
        Human_in_Loop_B.tmp_data_c[tmp_size_2[0] * g];
    }

    loop_ub = tmp_size_3[1];
    for (g = 0; g < loop_ub; g++) {
      Human_in_Loop_B.tmp_data[3 + (g << 2)] =
        Human_in_Loop_B.tmp_data_k[tmp_size_3[0] * g];
    }

    br = b_n;
    for (g = 0; g < b_n; g++) {
      Human_in_Loop_B.result_data[g << 2] = Human_in_Loop_B.tmp_data[g << 2];
      Human_in_Loop_B.result_data[1 + (g << 2)] = Human_in_Loop_B.tmp_data[(g <<
        2) + 1];
      Human_in_Loop_B.result_data[2 + (g << 2)] = Human_in_Loop_B.tmp_data[(g <<
        2) + 2];
      Human_in_Loop_B.result_data[3 + (g << 2)] = Human_in_Loop_B.tmp_data[(g <<
        2) + 3];
    }

    tmp_0 = (int16_T)br;
    tmp_size_idx_1 = tmp_0;
    b_n = br - 1;
    if (0 <= tmp_size_idx_1 - 1) {
      memset(&Human_in_Loop_B.tmp_data_cx[0], 0, tmp_size_idx_1 * sizeof(real_T));
    }

    if (br != 0) {
      for (br = 1; br - 1 <= b_n; br++) {
        for (e = br; e <= br; e++) {
          Human_in_Loop_B.tmp_data_cx[e - 1] = 0.0;
        }
      }

      br = 0;
      for (e = 0; e <= b_n; e++) {
        ar = -1;
        for (g = br; g + 1 <= br + 4; g++) {
          if (Human_in_Loop_B.result_data[g] != 0.0) {
            ia = ar;
            for (i = e; i + 1 <= e + 1; i++) {
              ia++;
              Human_in_Loop_B.tmp_data_cx[i] += Human_in_Loop_B.result_data[g] *
                parm1[ia];
            }
          }

          ar++;
        }

        br += 4;
      }
    }

    /* '<S12>:1:69' */
    for (g = 0; g < tmp_size_idx_1; g++) {
      Human_in_Loop_B.torque_track[b + g] = Human_in_Loop_B.tmp_data_cx[g];
    }

    /* '<S12>:1:72' */
    /* '<S12>:1:73' */
    /* '<S12>:1:74' */
    /* '<S12>:1:75' */
    for (g = 0; g < 750; g++) {
      Human_in_Loop_DW.TorqueMem[g << 2] = Human_in_Loop_B.torque_track[g];
      Human_in_Loop_DW.TorqueMem[1 + (g << 2)] = 0.0;
      Human_in_Loop_DW.TorqueMem[2 + (g << 2)] =
        Human_in_Loop_B.torque_delta_track[g];
      Human_in_Loop_DW.TorqueMem[3 + (g << 2)] = 0.0;
    }
  }

  /* '<S12>:1:78' */
  Human_in_Loop_DW.last_footstate = footstate;

  /* '<S12>:1:79' */
  Human_in_Loop_DW.TorqueMem[1 + (((int32_T)stride_index - 1) << 2)] =
    torque_measure;

  /* '<S12>:1:80' */
  Human_in_Loop_DW.TorqueMem[3 + (((int32_T)stride_index - 1) << 2)] =
    troque_delta;
  if (mode == 2.0) {
    /* '<S12>:1:82' */
    /* '<S12>:1:83' */
    mode = Human_in_Loop_DW.TorqueMem[((int32_T)stride_index - 1) << 2];

    /* '<S12>:1:84' */
    stride_index = Human_in_Loop_DW.TorqueMem[(((int32_T)stride_index - 1) << 2)
      + 2];
  } else {
    /* '<S12>:1:86' */
    mode = 0.0;

    /* '<S12>:1:87' */
    stride_index = 0.0;
  }

  /* '<S12>:1:90' */
  /* '<S12>:1:91' */
  Human_in_Loop_B.torque_des = mode;
  Human_in_Loop_B.torque_delta_des = stride_index;
  for (g = 0; g < 750; g++) {
    Human_in_Loop_B.torque_trace[g << 1] = Human_in_Loop_DW.TorqueMem[g << 2];
    Human_in_Loop_B.torque_trace[1 + (g << 1)] = Human_in_Loop_DW.TorqueMem[(g <<
      2) + 1];
  }

  for (g = 0; g < 750; g++) {
    Human_in_Loop_B.torque_delta_trace[g << 1] = Human_in_Loop_DW.TorqueMem[(g <<
      2) + 2];
    Human_in_Loop_B.torque_delta_trace[1 + (g << 1)] =
      Human_in_Loop_DW.TorqueMem[(g << 2) + 3];
  }

  /* End of MATLAB Function: '<S1>/Torque track' */

  /* MATLAB Function: '<S1>/LRN' */
  /* MATLAB Function 'Control Module/LRN': '<S10>:1' */
  if (!Human_in_Loop_DW.last_footstate_not_empty) {
    /* '<S10>:1:10' */
    Human_in_Loop_DW.last_footstate_not_empty = true;

    /* '<S10>:1:14' */
    Human_in_Loop_DW.last_torque_parm[0] = Human_in_Loop_B.RT3[0];
    Human_in_Loop_DW.last_torque_parm[1] = Human_in_Loop_B.RT3[2];
  }

  /* '<S10>:1:18' */
  /* '<S10>:1:22' */
  /* '<S10>:1:25' */
  /* '<S10>:1:26' */
  /* '<S10>:1:28' */
  /* '<S10>:1:30' */
  /* '<S10>:1:32' */
  stride_index = Human_in_Loop_B.RT1[3] * 500.0 + 1.0;
  if (stride_index > 750.0) {
    /* '<S10>:1:33' */
    /* '<S10>:1:34' */
    stride_index = 750.0;
  }

  if ((Human_in_Loop_DW.last_footstate_a == 0.0) && (Human_in_Loop_B.RT1[1] ==
       1.0) && (Human_in_Loop_B.RT1[0] == 2.0) && Human_in_Loop_P.LRN_BT_LRN_ON)
  {
    /* '<S10>:1:37' */
    /* '<S10>:1:38' */
    /* '<S10>:1:39' */
    mode = 1.0 - Human_in_Loop_P.LRN_error_filter_k;

    /* '<S10>:1:40' */
    footstate = Human_in_Loop_B.RT2[2];
    for (g = 0; g < 750; g++) {
      peak_torque = (Human_in_Loop_DW.TorqueMem[g << 2] -
                     Human_in_Loop_DW.TorqueMem[(g << 2) + 1]) *
        Human_in_Loop_P.LRN_error_filter_k + mode *
        Human_in_Loop_DW.torque_error_memory[g];
      Human_in_Loop_DW.torque_error_memory[g] = peak_torque;
      Human_in_Loop_DW.lrn_cmd_memory[g] = Human_in_Loop_P.LRN_lrn_shrink *
        Human_in_Loop_DW.lrn_cmd_memory[g] + footstate *
        Human_in_Loop_DW.torque_error_memory[g];
    }
  }

  if (Human_in_Loop_P.LRN_BT_LRN_CLEAR || (Human_in_Loop_DW.last_torque_parm[0]
       != Human_in_Loop_B.RT3[0]) || (Human_in_Loop_DW.last_torque_parm[1] !=
       Human_in_Loop_B.RT3[2])) {
    /* '<S10>:1:43' */
    /* '<S10>:1:44' */
    /* '<S10>:1:45' */
    memset(&Human_in_Loop_DW.torque_error_memory[0], 0, 1000U * sizeof(real_T));
    memset(&Human_in_Loop_DW.lrn_cmd_memory[0], 0, 1000U * sizeof(real_T));
  }

  if (Human_in_Loop_B.RT1[0] == 2.0) {
    /* '<S10>:1:48' */
    b = Human_in_Loop_P.LRN_time_delay;
    if (b < -2147482897) {
      b = MAX_int32_T;
    } else {
      b = 750 - b;
    }

    if (stride_index >= b) {
      /* '<S10>:1:49' */
      /* '<S10>:1:50' */
      Human_in_Loop_B.lrn_cmd = 0.0;
    } else {
      /* '<S10>:1:52' */
      mode = rt_roundd_snf(stride_index + (real_T)Human_in_Loop_P.LRN_time_delay);
      if (mode < 2.147483648E+9) {
        if (mode >= -2.147483648E+9) {
          g = (int32_T)mode;
        } else {
          g = MIN_int32_T;
        }
      } else {
        g = MAX_int32_T;
      }

      Human_in_Loop_B.lrn_cmd = Human_in_Loop_DW.lrn_cmd_memory[g - 1];
    }
  } else {
    /* '<S10>:1:55' */
    Human_in_Loop_B.lrn_cmd = 0.0;
  }

  /* '<S10>:1:58' */
  Human_in_Loop_DW.last_footstate_a = Human_in_Loop_B.RT1[1];

  /* '<S10>:1:59' */
  Human_in_Loop_DW.last_torque_parm[0] = Human_in_Loop_B.RT3[0];
  Human_in_Loop_DW.last_torque_parm[1] = Human_in_Loop_B.RT3[2];

  /* '<S10>:1:60' */
  memcpy(&Human_in_Loop_B.lrn_mem[0], &Human_in_Loop_DW.lrn_cmd_memory[0], 750U *
         sizeof(real_T));

  /* End of MATLAB Function: '<S1>/LRN' */

  /* MATLAB Function: '<S1>/Controller' */
  /* MATLAB Function 'Control Module/Controller': '<S9>:1' */
  /* '<S9>:1:23' */
  /* '<S9>:1:24' */
  /* '<S9>:1:26' */
  /* '<S9>:1:27' */
  /* '<S9>:1:28' */
  /* '<S9>:1:29' */
  /* '<S9>:1:32' */
  /* '<S9>:1:33' */
  /* '<S9>:1:37' */
  /* '<S9>:1:38' */
  /* '<S9>:1:39' */
  /* '<S9>:1:42' */
  /* '<S9>:1:43' */
  /* '<S9>:1:47' */
  /* '<S9>:1:48' */
  /* '<S9>:1:49' */
  /* '<S9>:1:50' */
  /* '<S9>:1:51' */
  /* '<S9>:1:52' */
  switch ((int32_T)Human_in_Loop_B.RT1[0]) {
   case 1:
    /* '<S9>:1:56' */
    Human_in_Loop_B.motor_vel_cmd = 0.0;

    /* '<S9>:1:57' */
    Human_in_Loop_DW.calib_state = 0.0;
    break;

   case 3:
    /* '<S9>:1:60' */
    Human_in_Loop_B.motor_vel_cmd = -150.0;
    break;

   case 4:
    if ((Human_in_Loop_DW.calib_state == 0.0) && (Human_in_Loop_B.RT4[0] < 5.0))
    {
      /* '<S9>:1:63' */
      /* '<S9>:1:64' */
      Human_in_Loop_B.motor_vel_cmd = 150.0;
    } else if ((Human_in_Loop_DW.calib_state == 0.0) && (Human_in_Loop_B.RT4[0] >
                5.0)) {
      /* '<S9>:1:65' */
      /* '<S9>:1:66' */
      Human_in_Loop_DW.calib_state = 1.0;

      /* '<S9>:1:67' */
      Human_in_Loop_B.motor_vel_cmd = 0.0;

      /* '<S9>:1:68' */
      /* '<S9>:1:69' */
      Human_in_Loop_DW.ParmReg[1] = (Human_in_Loop_DW.ParmReg[1] +
        Human_in_Loop_B.RT6[0]) - 20.0;
    } else {
      /* '<S9>:1:71' */
      Human_in_Loop_B.motor_vel_cmd = 0.0;
    }
    break;

   case 2:
    /* '<S9>:1:75' */
    switch (Human_in_Loop_P.Controller_MODE) {
     case 1:
      /* '<S9>:1:78' */
      /* '<S9>:1:79' */
      /* '<S9>:1:80' */
      Human_in_Loop_B.motor_vel_cmd = (((Human_in_Loop_B.RT5[0] +
        Human_in_Loop_P.Controller_FOLLOW_SLACK_ANGLE) - Human_in_Loop_B.RT6[0] *
        0.37037037037037035) * Human_in_Loop_B.RT2[3] + Human_in_Loop_B.RT5[1] *
        0.37037037037037035 * Human_in_Loop_B.RT2[4]) * 5.0 / 0.05;
      break;

     case 2:
      if (Human_in_Loop_B.RT1[1] == 1.0) {
        /* '<S9>:1:83' */
        /* '<S9>:1:84' */
        /* '<S9>:1:85' */
        /* '<S9>:1:86' */
        /* '<S9>:1:87' */
        Human_in_Loop_B.motor_vel_cmd = ((((Human_in_Loop_B.torque_des -
          Human_in_Loop_B.RT4[0]) * Human_in_Loop_B.RT2[0] +
          (Human_in_Loop_B.torque_delta_des - Human_in_Loop_B.RT4[1]) *
          Human_in_Loop_B.RT2[1]) + Human_in_Loop_B.lrn_cmd) +
          Human_in_Loop_B.RT2[5] * Human_in_Loop_B.torque_delta_des) * 5.0 /
          0.05;
      } else {
        /* '<S9>:1:89' */
        /* '<S9>:1:90' */
        Human_in_Loop_B.motor_vel_cmd = (0.0 - Human_in_Loop_B.RT6[0] *
          0.37037037037037035) * 2.0 * 5.0 / 0.05;
      }
      break;

     default:
      /* '<S9>:1:93' */
      Human_in_Loop_B.motor_vel_cmd = 0.0;
      break;
    }
    break;

   case 0:
    /* '<S9>:1:97' */
    Human_in_Loop_B.motor_vel_cmd = 0.0;
    break;

   default:
    /* '<S9>:1:100' */
    Human_in_Loop_B.motor_vel_cmd = 0.0;
    break;
  }

  /* End of MATLAB Function: '<S1>/Controller' */

  /* MATLAB Function: '<S11>/MATLAB Function' */
  /* MATLAB Function 'Control Module/Motor/MATLAB Function': '<S14>:1' */
  /* '<S14>:1:3' */
  /* '<S14>:1:4' */
  if (Human_in_Loop_B.RT4[0] > Human_in_Loop_P.MATLABFunction_MAX_TORQUE) {
    /* '<S14>:1:6' */
    /* '<S14>:1:7' */
    Human_in_Loop_B.vel = 0.0;
  } else if (Human_in_Loop_B.RT6[0] >
             Human_in_Loop_P.MATLABFunction_MAX_MOTOR_ANGLE) {
    /* '<S14>:1:8' */
    /* '<S14>:1:9' */
    Human_in_Loop_B.vel = 0.0;
  } else if (Human_in_Loop_B.RT6[0] <
             Human_in_Loop_P.MATLABFunction_MIN_MOTOR_ANGLE) {
    /* '<S14>:1:10' */
    /* '<S14>:1:11' */
    Human_in_Loop_B.vel = 0.0;
  } else if (Human_in_Loop_B.motor_vel_cmd >
             Human_in_Loop_P.MATLABFunction_MAX_SPEED) {
    /* '<S14>:1:12' */
    /* '<S14>:1:13' */
    Human_in_Loop_B.vel = Human_in_Loop_P.MATLABFunction_MAX_SPEED;
  } else if (Human_in_Loop_B.motor_vel_cmd <
             -Human_in_Loop_P.MATLABFunction_MAX_SPEED) {
    /* '<S14>:1:14' */
    /* '<S14>:1:15' */
    Human_in_Loop_B.vel = -Human_in_Loop_P.MATLABFunction_MAX_SPEED;
  } else {
    /* '<S14>:1:17' */
    Human_in_Loop_B.vel = Human_in_Loop_B.motor_vel_cmd;
  }

  /* End of MATLAB Function: '<S11>/MATLAB Function' */

  /* Gain: '<S11>/Gain2' */
  Human_in_Loop_B.Gain2_o = Human_in_Loop_P.Gain2_Gain * Human_in_Loop_B.vel;

  /* Gain: '<S11>/Gain1' */
  Human_in_Loop_B.Gain1_h = Human_in_Loop_P.Gain1_Gain * Human_in_Loop_B.Gain2_o;

  /* S-Function (rti_commonblock): '<S13>/S-Function1' */
  /* This comment workarounds a code generation problem */

  /* --- Human_in_Loop/Control Module/Motor/DAC_CLASS1_BL1 --- */
  /* --- [RTI120X, DAC C1] - Channel: 16 --- */
  {
    /* define variables required for DAC realtime functions */
    Float64 inportDacData= 0.0;
    inportDacData = (real_T) Human_in_Loop_B.Gain1_h;

    /* write value of CL1 DAC for output channel 16 */
    DacCl1AnalogOut_setOutputValue(pRTIDacC1AnalogOut_Ch_16,
      DAC_CLASS1_CHANNEL_16, inportDacData);
    DacCl1AnalogOut_write(pRTIDacC1AnalogOut_Ch_16);
  }
}

/* Termination for function-call system: '<Root>/Control Module' */
void Human_in_Loo_ControlModule_Term(void)
{
  /* Terminate for S-Function (rti_commonblock): '<S13>/S-Function1' */

  /* --- Human_in_Loop/Control Module/Motor/DAC_CLASS1_BL1 --- */
  /* --- [RTI120X, DAC C1] - Channel: 16 --- */

  /* All channel outputs are set to high impedance state */
  DacCl1AnalogOut_setOutputHighZ(pRTIDacC1AnalogOut_Ch_16, DAC_CLASS1_HIGH_Z_ON);

  /* Deactivates AnalogOut functionality */
  DacCl1AnalogOut_stop(pRTIDacC1AnalogOut_Ch_16);
}

/* Function for MATLAB Function: '<S3>/torque_track_loss' */
static real_T Human_in_Loop_norm(const real_T x_data[], const int32_T x_size[2])
{
  real_T y;
  real_T scale;
  real_T absxk;
  real_T t;
  int32_T k;
  if (x_size[1] == 0) {
    y = 0.0;
  } else {
    y = 0.0;
    if (x_size[1] == 1) {
      y = fabs(x_data[0]);
    } else {
      scale = 3.3121686421112381E-170;
      for (k = 1; k <= x_size[1]; k++) {
        absxk = fabs(x_data[k - 1]);
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

/* Function for MATLAB Function: '<S3>/torque_track_loss' */
static real_T Human_in_Loop_mean(const real_T x[10])
{
  real_T b_y;
  int32_T k;
  b_y = x[0];
  for (k = 0; k < 9; k++) {
    b_y += x[k + 1];
  }

  return b_y / 10.0;
}

/* Function for MATLAB Function: '<S32>/MATLAB Function' */
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

/* Function for MATLAB Function: '<S32>/MATLAB Function' */
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

/* Function for MATLAB Function: '<S32>/MATLAB Function' */
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
  int_T tid = 0;
  real_T torque_zero;
  real_T x[6];
  real_T varargin_1[2];
  int32_T ixstart;
  real_T x_0[2];
  real_T A[30];
  int32_T rankR;
  real_T B[15];
  static const int8_T b_A[30] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

  int32_T xoffset;
  int32_T c;
  real_T B_0;
  ZCEventType zcEvent;
  int32_T i;
  real_T tmp[10];
  int32_T tmp_size[2];
  boolean_T exitg1;
  if (rtmIsMajorTimeStep(Human_in_Loop_M)) {
    /* set solver stop time */
    if (!(Human_in_Loop_M->Timing.clockTick0+1)) {
      rtsiSetSolverStopTime(&Human_in_Loop_M->solverInfo,
                            ((Human_in_Loop_M->Timing.clockTickH0 + 1) *
        Human_in_Loop_M->Timing.stepSize0 * 4294967296.0));
    } else {
      rtsiSetSolverStopTime(&Human_in_Loop_M->solverInfo,
                            ((Human_in_Loop_M->Timing.clockTick0 + 1) *
        Human_in_Loop_M->Timing.stepSize0 + Human_in_Loop_M->Timing.clockTickH0 *
        Human_in_Loop_M->Timing.stepSize0 * 4294967296.0));
    }
  }                                    /* end MajorTimeStep */

  /* Update absolute time of base rate at minor time step */
  if (rtmIsMinorTimeStep(Human_in_Loop_M)) {
    Human_in_Loop_M->Timing.t[0] = rtsiGetT(&Human_in_Loop_M->solverInfo);
  }

  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S77>/S-Function1' */
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
      AdcCl1AnalogIn_getSingleValue(pRTIAdcC1AnalogIn_Ch_6, readStateFlag,
        (real_T*) &Human_in_Loop_B.SFunction1);
    }

    /* Gain: '<S32>/Gain' */
    Human_in_Loop_B.Gain = Human_in_Loop_P.Gain_Gain *
      Human_in_Loop_B.SFunction1;

    /* DiscreteFilter: '<S32>/0.4low2' */
    torque_zero = Human_in_Loop_B.Gain;
    torque_zero -= Human_in_Loop_P.u4low2_DenCoef[1] *
      Human_in_Loop_DW.u4low2_states[0];
    torque_zero -= Human_in_Loop_P.u4low2_DenCoef[2] *
      Human_in_Loop_DW.u4low2_states[1];
    torque_zero -= Human_in_Loop_P.u4low2_DenCoef[3] *
      Human_in_Loop_DW.u4low2_states[2];
    torque_zero /= Human_in_Loop_P.u4low2_DenCoef[0];
    Human_in_Loop_DW.u4low2_tmp = torque_zero;
    torque_zero = Human_in_Loop_P.u4low2_NumCoef[0] *
      Human_in_Loop_DW.u4low2_tmp;
    torque_zero += Human_in_Loop_P.u4low2_NumCoef[1] *
      Human_in_Loop_DW.u4low2_states[0];
    torque_zero += Human_in_Loop_P.u4low2_NumCoef[2] *
      Human_in_Loop_DW.u4low2_states[1];
    torque_zero += Human_in_Loop_P.u4low2_NumCoef[3] *
      Human_in_Loop_DW.u4low2_states[2];
    Human_in_Loop_B.u4low2 = torque_zero;

    /* MATLAB Function: '<S32>/Data process' */
    /* MATLAB Function 'Sensor Data/Torque module/Data process': '<S79>:1' */
    /* '<S79>:1:4' */
    torque_zero = Human_in_Loop_DW.ParmReg[0];
    if (Human_in_Loop_P.Dataprocess_BT_RESET_TORQUE) {
      /* '<S79>:1:7' */
      /* '<S79>:1:8' */
      torque_zero = Human_in_Loop_B.u4low2 *
        Human_in_Loop_P.Dataprocess_load_vol_gain +
        Human_in_Loop_P.Dataprocess_load_vol_offset;

      /* '<S79>:1:9' */
      Human_in_Loop_DW.ParmReg[0] = torque_zero;
    }

    /* '<S79>:1:12' */
    Human_in_Loop_B.torque = (Human_in_Loop_B.u4low2 *
      Human_in_Loop_P.Dataprocess_load_vol_gain +
      Human_in_Loop_P.Dataprocess_load_vol_offset) - torque_zero;

    /* End of MATLAB Function: '<S32>/Data process' */

    /* UnitDelay: '<S76>/Unit Delay1' */
    Human_in_Loop_B.x2k1 = Human_in_Loop_DW.UnitDelay1_DSTATE;

    /* UnitDelay: '<S76>/Unit Delay' */
    Human_in_Loop_B.x1k1 = Human_in_Loop_DW.UnitDelay_DSTATE;

    /* Gain: '<S76>/Gain1' */
    B_0 = Human_in_Loop_P.uOrderTD_T1 * Human_in_Loop_P.uOrderTD_T2;
    torque_zero = 1.0 / B_0;
    Human_in_Loop_B.Gain1 = torque_zero * Human_in_Loop_B.x1k1;

    /* Gain: '<S76>/Gain2' */
    torque_zero = Human_in_Loop_P.uOrderTD_T1 + Human_in_Loop_P.uOrderTD_T2;
    B_0 = Human_in_Loop_P.uOrderTD_T1 * Human_in_Loop_P.uOrderTD_T2;
    torque_zero /= B_0;
    Human_in_Loop_B.Gain2 = torque_zero * Human_in_Loop_B.x2k1;

    /* UnitDelay: '<S76>/Unit Delay2' */
    Human_in_Loop_B.UnitDelay2 = Human_in_Loop_DW.UnitDelay2_DSTATE;

    /* Gain: '<S76>/Gain4' */
    B_0 = Human_in_Loop_P.uOrderTD_T1 * Human_in_Loop_P.uOrderTD_T2;
    torque_zero = 1.0 / B_0;
    Human_in_Loop_B.Gain4 = torque_zero * Human_in_Loop_B.UnitDelay2;

    /* Sum: '<S76>/Add2' */
    Human_in_Loop_B.Add2 = (Human_in_Loop_B.Gain1 + Human_in_Loop_B.Gain2) -
      Human_in_Loop_B.Gain4;

    /* Gain: '<S76>/Gain3' */
    Human_in_Loop_B.Gain3 = Human_in_Loop_P.uOrderTD_Ts * Human_in_Loop_B.Add2;

    /* Sum: '<S76>/Add1' */
    Human_in_Loop_B.x2k = Human_in_Loop_B.x2k1 - Human_in_Loop_B.Gain3;

    /* MATLAB Function: '<S32>/Mux' */
    Human_in_Loop_Mux(Human_in_Loop_B.torque, Human_in_Loop_B.x2k,
                      &Human_in_Loop_B.sf_Mux);

    /* RateTransition: '<Root>/RT4' */
    switch (Human_in_Loop_DW.RT4_read_buf) {
     case 0:
      Human_in_Loop_DW.RT4_write_buf = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT4_write_buf = 0;
      break;

     default:
      Human_in_Loop_DW.RT4_write_buf = (int8_T)(Human_in_Loop_DW.RT4_last_buf_wr
        == 0);
      break;
    }

    if (Human_in_Loop_DW.RT4_write_buf != 0) {
      Human_in_Loop_DW.RT4_Buffer1[0] = Human_in_Loop_B.sf_Mux.x[0];
      Human_in_Loop_DW.RT4_Buffer1[1] = Human_in_Loop_B.sf_Mux.x[1];
    } else {
      Human_in_Loop_DW.RT4_Buffer0[0] = Human_in_Loop_B.sf_Mux.x[0];
      Human_in_Loop_DW.RT4_Buffer0[1] = Human_in_Loop_B.sf_Mux.x[1];
    }

    Human_in_Loop_DW.RT4_last_buf_wr = Human_in_Loop_DW.RT4_write_buf;
    Human_in_Loop_DW.RT4_write_buf = -1;

    /* End of RateTransition: '<Root>/RT4' */

    /* S-Function (rti_commonblock): '<S56>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* Gain: '<S29>/Gain' */
    Human_in_Loop_B.Gain_g = Human_in_Loop_P.Gain_Gain_o *
      Human_in_Loop_B.SFunction1_o1;

    /* MATLAB Function: '<S29>/Data process' */
    /* MATLAB Function 'Sensor Data/Encoder module/Data process': '<S53>:1' */
    /* '<S53>:1:3' */
    torque_zero = Human_in_Loop_DW.ParmReg[2];
    if (Human_in_Loop_P.Dataprocess_BT_RESET_ANKLE) {
      /* '<S53>:1:6' */
      /* '<S53>:1:7' */
      torque_zero = Human_in_Loop_B.Gain_g;

      /* '<S53>:1:8' */
      Human_in_Loop_DW.ParmReg[2] = Human_in_Loop_B.Gain_g;
    }

    /* '<S53>:1:11' */
    Human_in_Loop_B.angle_k = Human_in_Loop_B.Gain_g - torque_zero;

    /* End of MATLAB Function: '<S29>/Data process' */

    /* Gain: '<S29>/Gain1' */
    Human_in_Loop_B.Gain1_o = Human_in_Loop_P.Gain1_Gain_m *
      Human_in_Loop_B.SFunction1_o2;

    /* MATLAB Function: '<S29>/Mux' */
    Human_in_Loop_Mux(Human_in_Loop_B.angle_k, Human_in_Loop_B.Gain1_o,
                      &Human_in_Loop_B.sf_Mux_p);

    /* RateTransition: '<Root>/RT5' */
    switch (Human_in_Loop_DW.RT5_read_buf) {
     case 0:
      Human_in_Loop_DW.RT5_write_buf = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT5_write_buf = 0;
      break;

     default:
      Human_in_Loop_DW.RT5_write_buf = (int8_T)(Human_in_Loop_DW.RT5_last_buf_wr
        == 0);
      break;
    }

    if (Human_in_Loop_DW.RT5_write_buf != 0) {
      Human_in_Loop_DW.RT5_Buffer1[0] = Human_in_Loop_B.sf_Mux_p.x[0];
      Human_in_Loop_DW.RT5_Buffer1[1] = Human_in_Loop_B.sf_Mux_p.x[1];
    } else {
      Human_in_Loop_DW.RT5_Buffer0[0] = Human_in_Loop_B.sf_Mux_p.x[0];
      Human_in_Loop_DW.RT5_Buffer0[1] = Human_in_Loop_B.sf_Mux_p.x[1];
    }

    Human_in_Loop_DW.RT5_last_buf_wr = Human_in_Loop_DW.RT5_write_buf;
    Human_in_Loop_DW.RT5_write_buf = -1;

    /* End of RateTransition: '<Root>/RT5' */

    /* S-Function (rti_commonblock): '<S57>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* Gain: '<S29>/Gain2' */
    Human_in_Loop_B.Gain2_c = Human_in_Loop_P.Gain2_Gain_h *
      Human_in_Loop_B.SFunction1_o1_j;

    /* MATLAB Function: '<S29>/Data process1' */
    /* MATLAB Function 'Sensor Data/Encoder module/Data process1': '<S54>:1' */
    /* '<S54>:1:3' */
    torque_zero = Human_in_Loop_DW.ParmReg[1];
    if (Human_in_Loop_P.Dataprocess1_BT_RESET_MOTOR) {
      /* '<S54>:1:5' */
      /* '<S54>:1:6' */
      torque_zero = Human_in_Loop_B.Gain2_c;

      /* '<S54>:1:7' */
      Human_in_Loop_DW.ParmReg[1] = Human_in_Loop_B.Gain2_c;
    }

    /* '<S54>:1:10' */
    Human_in_Loop_B.angle = Human_in_Loop_B.Gain2_c - torque_zero;

    /* End of MATLAB Function: '<S29>/Data process1' */

    /* Gain: '<S29>/Gain3' */
    Human_in_Loop_B.Gain3_h = Human_in_Loop_P.Gain3_Gain *
      Human_in_Loop_B.SFunction1_o2_j;

    /* UnitDelay: '<S52>/Unit Delay1' */
    Human_in_Loop_B.x2k1_d = Human_in_Loop_DW.UnitDelay1_DSTATE_f;

    /* UnitDelay: '<S52>/Unit Delay' */
    Human_in_Loop_B.x1k1_e = Human_in_Loop_DW.UnitDelay_DSTATE_a;

    /* Gain: '<S52>/Gain1' */
    B_0 = Human_in_Loop_P.uOrderTD_T1_p * Human_in_Loop_P.uOrderTD_T2_k;
    torque_zero = 1.0 / B_0;
    Human_in_Loop_B.Gain1_p = torque_zero * Human_in_Loop_B.x1k1_e;

    /* Gain: '<S52>/Gain2' */
    torque_zero = Human_in_Loop_P.uOrderTD_T1_p + Human_in_Loop_P.uOrderTD_T2_k;
    B_0 = Human_in_Loop_P.uOrderTD_T1_p * Human_in_Loop_P.uOrderTD_T2_k;
    torque_zero /= B_0;
    Human_in_Loop_B.Gain2_n = torque_zero * Human_in_Loop_B.x2k1_d;

    /* UnitDelay: '<S52>/Unit Delay2' */
    Human_in_Loop_B.UnitDelay2_i = Human_in_Loop_DW.UnitDelay2_DSTATE_e;

    /* Gain: '<S52>/Gain4' */
    B_0 = Human_in_Loop_P.uOrderTD_T1_p * Human_in_Loop_P.uOrderTD_T2_k;
    torque_zero = 1.0 / B_0;
    Human_in_Loop_B.Gain4_o = torque_zero * Human_in_Loop_B.UnitDelay2_i;

    /* Sum: '<S52>/Add2' */
    Human_in_Loop_B.Add2_f = (Human_in_Loop_B.Gain1_p + Human_in_Loop_B.Gain2_n)
      - Human_in_Loop_B.Gain4_o;

    /* Gain: '<S52>/Gain3' */
    Human_in_Loop_B.Gain3_hf = Human_in_Loop_P.uOrderTD_Ts_n *
      Human_in_Loop_B.Add2_f;

    /* Sum: '<S52>/Add1' */
    Human_in_Loop_B.x2k_a = Human_in_Loop_B.x2k1_d - Human_in_Loop_B.Gain3_hf;

    /* MATLAB Function: '<S29>/Mux1' */
    /* MATLAB Function 'Sensor Data/Encoder module/Mux1': '<S61>:1' */
    /* '<S61>:1:3' */
    Human_in_Loop_B.x_f[0] = Human_in_Loop_B.angle;
    Human_in_Loop_B.x_f[1] = Human_in_Loop_B.Gain3_h;
    Human_in_Loop_B.x_f[2] = Human_in_Loop_B.x2k_a;

    /* RateTransition: '<Root>/RT6' */
    switch (Human_in_Loop_DW.RT6_read_buf) {
     case 0:
      Human_in_Loop_DW.RT6_write_buf = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT6_write_buf = 0;
      break;

     default:
      Human_in_Loop_DW.RT6_write_buf = (int8_T)(Human_in_Loop_DW.RT6_last_buf_wr
        == 0);
      break;
    }

    if (Human_in_Loop_DW.RT6_write_buf != 0) {
      Human_in_Loop_DW.RT6_Buffer1[0] = Human_in_Loop_B.x_f[0];
      Human_in_Loop_DW.RT6_Buffer1[1] = Human_in_Loop_B.x_f[1];
      Human_in_Loop_DW.RT6_Buffer1[2] = Human_in_Loop_B.x_f[2];
    } else {
      Human_in_Loop_DW.RT6_Buffer0[0] = Human_in_Loop_B.x_f[0];
      Human_in_Loop_DW.RT6_Buffer0[1] = Human_in_Loop_B.x_f[1];
      Human_in_Loop_DW.RT6_Buffer0[2] = Human_in_Loop_B.x_f[2];
    }

    Human_in_Loop_DW.RT6_last_buf_wr = Human_in_Loop_DW.RT6_write_buf;
    Human_in_Loop_DW.RT6_write_buf = -1;

    /* End of RateTransition: '<Root>/RT6' */

    /* S-Function (rti_commonblock): '<S63>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* MATLAB Function: '<S30>/FootSwitch Filter' */
    /* MATLAB Function 'Sensor Data/FootSwitch module/FootSwitch Filter': '<S65>:1' */
    /* '<S65>:1:6' */
    if (Human_in_Loop_DW.foot_state == 0.0) {
      /* '<S65>:1:15' */
      if (Human_in_Loop_B.SFunction1_it) {
        /* '<S65>:1:16' */
        /* '<S65>:1:17' */
        Human_in_Loop_DW.foot_state = 1.0;
      } else {
        /* '<S65>:1:19' */
        Human_in_Loop_DW.foot_state = 0.0;
      }
    } else if (Human_in_Loop_B.SFunction1_it) {
      /* '<S65>:1:22' */
      /* '<S65>:1:23' */
      Human_in_Loop_DW.foot_state = 1.0;
    } else {
      /* '<S65>:1:25' */
      Human_in_Loop_DW.filter_time += 0.0002;
      if (Human_in_Loop_DW.filter_time > 0.1) {
        /* '<S65>:1:26' */
        /* '<S65>:1:27' */
        Human_in_Loop_DW.filter_time = 0.0;

        /* '<S65>:1:28' */
        Human_in_Loop_DW.foot_state = 0.0;
      }
    }

    /* '<S65>:1:33' */
    Human_in_Loop_B.state_m = Human_in_Loop_DW.foot_state;

    /* End of MATLAB Function: '<S30>/FootSwitch Filter' */

    /* MATLAB Function: '<S7>/State Machine' */
    /* MATLAB Function 'State Module/State Machine': '<S83>:1' */
    /* '<S83>:1:21' */
    /* '<S83>:1:22' */
    /* '<S83>:1:23' */
    /* '<S83>:1:24' */
    /* '<S83>:1:27' */
    /* '<S83>:1:28' */
    /* '<S83>:1:29' */
    /* '<S83>:1:30' */
    /* '<S83>:1:31' */
    if (Human_in_Loop_P.StateMachine_BT_RUN) {
      /* '<S83>:1:34' */
      /* '<S83>:1:35' */
      Human_in_Loop_DW.bt_run = 1.0;
    }

    if (Human_in_Loop_P.StateMachine_BT_CALIB) {
      /* '<S83>:1:37' */
      /* '<S83>:1:38' */
      Human_in_Loop_DW.reg_mode = 4.0;
    }

    if (Human_in_Loop_P.StateMachine_BT_SLACK) {
      /* '<S83>:1:40' */
      /* '<S83>:1:41' */
      Human_in_Loop_DW.reg_mode = 3.0;
    }

    if (Human_in_Loop_P.StateMachine_BT_IDLE) {
      /* '<S83>:1:43' */
      /* '<S83>:1:44' */
      Human_in_Loop_DW.reg_mode = 1.0;
    }

    if (Human_in_Loop_P.StateMachine_BT_ERROR) {
      /* '<S83>:1:46' */
      /* '<S83>:1:47' */
      Human_in_Loop_DW.reg_mode = 0.0;
    }

    if ((Human_in_Loop_DW.bt_run == 1.0) && (Human_in_Loop_DW.reg_last_switch ==
         0.0) && (Human_in_Loop_B.state_m == 1.0)) {
      /* '<S83>:1:50' */
      /* '<S83>:1:51' */
      Human_in_Loop_DW.reg_mode = 2.0;

      /* '<S83>:1:52' */
      Human_in_Loop_DW.bt_run = 0.0;
    }

    if ((Human_in_Loop_DW.reg_mode == 2.0) || (Human_in_Loop_DW.reg_mode == 1.0))
    {
      /* '<S83>:1:55' */
      if ((Human_in_Loop_DW.reg_last_switch == 0.0) && (Human_in_Loop_B.state_m ==
           1.0)) {
        /* '<S83>:1:56' */
        /* '<S83>:1:57' */
        Human_in_Loop_DW.reg_state = 1.0;
        if ((Human_in_Loop_DW.reg_stride_time_count < 1.5) &&
            (Human_in_Loop_DW.reg_stride_time_count > 0.5)) {
          /* '<S83>:1:59' */
          /* '<S83>:1:60' */
          Human_in_Loop_DW.reg_stride_time = 0.618 *
            Human_in_Loop_DW.reg_stride_time + 0.382 *
            Human_in_Loop_DW.reg_stride_time_count;
        }

        /* '<S83>:1:63' */
        Human_in_Loop_DW.reg_stride_time_count = 0.0;
      } else if ((Human_in_Loop_DW.reg_state == 1.0) &&
                 (Human_in_Loop_DW.reg_stride_time_count > 0.65 *
                  Human_in_Loop_DW.reg_stride_time)) {
        /* '<S83>:1:65' */
        /* '<S83>:1:66' */
        Human_in_Loop_DW.reg_state = 0.0;

        /* '<S83>:1:67' */
        Human_in_Loop_DW.reg_stride_time_count += 0.0002;
      } else {
        /* '<S83>:1:69' */
        Human_in_Loop_DW.reg_stride_time_count += 0.0002;
      }
    }

    /* '<S83>:1:73' */
    Human_in_Loop_DW.reg_last_switch = Human_in_Loop_B.state_m;
    if (Human_in_Loop_DW.reg_stride_time > 1.5) {
      /* '<S83>:1:74' */
      /* '<S83>:1:75' */
      Human_in_Loop_DW.reg_stride_time = 1.5;
    } else {
      if (Human_in_Loop_DW.reg_stride_time < 0.5) {
        /* '<S83>:1:76' */
        /* '<S83>:1:77' */
        Human_in_Loop_DW.reg_stride_time = 0.5;
      }
    }

    /* '<S83>:1:80' */
    Human_in_Loop_B.mode = Human_in_Loop_DW.reg_mode;

    /* '<S83>:1:81' */
    Human_in_Loop_B.state = Human_in_Loop_DW.reg_state;

    /* '<S83>:1:82' */
    Human_in_Loop_B.stride_time = Human_in_Loop_DW.reg_stride_time;

    /* '<S83>:1:83' */
    Human_in_Loop_B.stride_timer = Human_in_Loop_DW.reg_stride_time_count;

    /* End of MATLAB Function: '<S7>/State Machine' */

    /* MATLAB Function: '<S7>/Mux1' */
    /* MATLAB Function 'State Module/Mux1': '<S82>:1' */
    /* '<S82>:1:3' */
    Human_in_Loop_B.x[0] = Human_in_Loop_B.mode;
    Human_in_Loop_B.x[1] = Human_in_Loop_B.state;
    Human_in_Loop_B.x[2] = Human_in_Loop_B.stride_time;
    Human_in_Loop_B.x[3] = Human_in_Loop_B.stride_timer;

    /* RateTransition: '<Root>/RT1' */
    switch (Human_in_Loop_DW.RT1_read_buf) {
     case 0:
      Human_in_Loop_DW.RT1_write_buf = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT1_write_buf = 0;
      break;

     default:
      Human_in_Loop_DW.RT1_write_buf = (int8_T)(Human_in_Loop_DW.RT1_last_buf_wr
        == 0);
      break;
    }

    if (Human_in_Loop_DW.RT1_write_buf != 0) {
      Human_in_Loop_DW.RT1_Buffer1[0] = Human_in_Loop_B.x[0];
      Human_in_Loop_DW.RT1_Buffer1[1] = Human_in_Loop_B.x[1];
      Human_in_Loop_DW.RT1_Buffer1[2] = Human_in_Loop_B.x[2];
      Human_in_Loop_DW.RT1_Buffer1[3] = Human_in_Loop_B.x[3];
    } else {
      Human_in_Loop_DW.RT1_Buffer0[0] = Human_in_Loop_B.x[0];
      Human_in_Loop_DW.RT1_Buffer0[1] = Human_in_Loop_B.x[1];
      Human_in_Loop_DW.RT1_Buffer0[2] = Human_in_Loop_B.x[2];
      Human_in_Loop_DW.RT1_Buffer0[3] = Human_in_Loop_B.x[3];
    }

    Human_in_Loop_DW.RT1_last_buf_wr = Human_in_Loop_DW.RT1_write_buf;
    Human_in_Loop_DW.RT1_write_buf = -1;

    /* End of RateTransition: '<Root>/RT1' */

    /* MATLAB Function: '<S21>/Mux1' incorporates:
     *  Constant: '<S21>/Kd'
     *  Constant: '<S21>/Kl'
     *  Constant: '<S21>/Ko'
     *  Constant: '<S21>/Kp'
     *  Constant: '<S21>/Ksd'
     *  Constant: '<S21>/Ksp'
     *  Constant: '<S21>/control reset'
     */
    /* MATLAB Function 'Parameter Module/Control Parameter/Mux1': '<S23>:1' */
    if (Human_in_Loop_P.controlreset_Value == 1.0) {
      /* '<S23>:1:3' */
      /* '<S23>:1:4' */
      for (i = 0; i < 6; i++) {
        Human_in_Loop_B.x_h[i] = 0.0;
      }
    } else {
      /* '<S23>:1:6' */
      Human_in_Loop_B.x_h[0] = Human_in_Loop_P.Kp_Value;
      Human_in_Loop_B.x_h[1] = Human_in_Loop_P.Kd_Value;
      Human_in_Loop_B.x_h[2] = Human_in_Loop_P.Kl_Value;
      Human_in_Loop_B.x_h[3] = Human_in_Loop_P.Ksp_Value;
      Human_in_Loop_B.x_h[4] = Human_in_Loop_P.Ksd_Value;
      Human_in_Loop_B.x_h[5] = Human_in_Loop_P.Ko_Value;
    }

    /* End of MATLAB Function: '<S21>/Mux1' */

    /* RateTransition: '<Root>/RT2' */
    switch (Human_in_Loop_DW.RT2_read_buf) {
     case 0:
      Human_in_Loop_DW.RT2_write_buf = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT2_write_buf = 0;
      break;

     default:
      Human_in_Loop_DW.RT2_write_buf = (int8_T)(Human_in_Loop_DW.RT2_last_buf_wr
        == 0);
      break;
    }

    if (Human_in_Loop_DW.RT2_write_buf != 0) {
      for (i = 0; i < 6; i++) {
        Human_in_Loop_DW.RT2_Buffer1[i] = Human_in_Loop_B.x_h[i];
      }
    } else {
      for (i = 0; i < 6; i++) {
        Human_in_Loop_DW.RT2_Buffer0[i] = Human_in_Loop_B.x_h[i];
      }
    }

    Human_in_Loop_DW.RT2_last_buf_wr = Human_in_Loop_DW.RT2_write_buf;
    Human_in_Loop_DW.RT2_write_buf = -1;

    /* End of RateTransition: '<Root>/RT2' */

    /* MATLAB Function: '<S22>/Mux1' incorporates:
     *  Constant: '<S22>/fall_time'
     *  Constant: '<S22>/peak_time'
     *  Constant: '<S22>/peak_torque'
     *  Constant: '<S22>/rise_time'
     */
    /* MATLAB Function 'Parameter Module/Torque Parameter/Mux1': '<S24>:1' */
    /* '<S24>:1:3' */
    Human_in_Loop_B.x_k[0] = Human_in_Loop_P.peak_torque_Value;
    Human_in_Loop_B.x_k[1] = Human_in_Loop_P.rise_time_Value;
    Human_in_Loop_B.x_k[2] = Human_in_Loop_P.peak_time_Value;
    Human_in_Loop_B.x_k[3] = Human_in_Loop_P.fall_time_Value;

    /* RateTransition: '<Root>/RT3' */
    switch (Human_in_Loop_DW.RT3_read_buf) {
     case 0:
      Human_in_Loop_DW.RT3_write_buf = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT3_write_buf = 0;
      break;

     default:
      Human_in_Loop_DW.RT3_write_buf = (int8_T)(Human_in_Loop_DW.RT3_last_buf_wr
        == 0);
      break;
    }

    if (Human_in_Loop_DW.RT3_write_buf != 0) {
      Human_in_Loop_DW.RT3_Buffer1[0] = Human_in_Loop_B.x_k[0];
      Human_in_Loop_DW.RT3_Buffer1[1] = Human_in_Loop_B.x_k[1];
      Human_in_Loop_DW.RT3_Buffer1[2] = Human_in_Loop_B.x_k[2];
      Human_in_Loop_DW.RT3_Buffer1[3] = Human_in_Loop_B.x_k[3];
    } else {
      Human_in_Loop_DW.RT3_Buffer0[0] = Human_in_Loop_B.x_k[0];
      Human_in_Loop_DW.RT3_Buffer0[1] = Human_in_Loop_B.x_k[1];
      Human_in_Loop_DW.RT3_Buffer0[2] = Human_in_Loop_B.x_k[2];
      Human_in_Loop_DW.RT3_Buffer0[3] = Human_in_Loop_B.x_k[3];
    }

    Human_in_Loop_DW.RT3_last_buf_wr = Human_in_Loop_DW.RT3_write_buf;
    Human_in_Loop_DW.RT3_write_buf = -1;

    /* End of RateTransition: '<Root>/RT3' */

    /* S-Function (rti_commonblock): '<S8>/S-Function1' */

    /* This comment workarounds a code generation problem */

    /* End of Outputs for S-Function (rti_commonblock): '<S8>/S-Function1' */

    /* Gain: '<S3>/Gain' */
    Human_in_Loop_B.Gain_i = Human_in_Loop_P.Gain_Gain_j * Human_in_Loop_B.x[3];

    /* Delay: '<S16>/Delay3' */
    Human_in_Loop_B.Delay3 = Human_in_Loop_DW.Delay3_DSTATE;

    /* Delay: '<S16>/Delay2' */
    Human_in_Loop_B.Delay2 = Human_in_Loop_DW.Delay2_DSTATE;

    /* Delay: '<S16>/Delay7' */
    Human_in_Loop_B.Delay7 = Human_in_Loop_DW.Delay7_DSTATE;

    /* Delay: '<S16>/Delay6' */
    Human_in_Loop_B.Delay6 = Human_in_Loop_DW.Delay6_DSTATE;
  }

  /* StateSpace: '<S40>/low_pass' */
  Human_in_Loop_B.low_pass = 0.0;
  Human_in_Loop_B.low_pass += Human_in_Loop_P.low_pass_C *
    Human_in_Loop_X.low_pass_CSTATE[1];

  /* MATLAB Function: '<S40>/MVC' */
  /* MATLAB Function 'Sensor Data/EMG module/Preprocessing/MVC': '<S46>:1' */
  /* '<S46>:1:3' */
  Human_in_Loop_B.y_b = (Human_in_Loop_B.low_pass -
    Human_in_Loop_P.MVC_r_SOL_MIN) / (Human_in_Loop_P.MVC_r_SOL_MAX -
    Human_in_Loop_P.MVC_r_SOL_MIN);

  /* StateSpace: '<S41>/low_pass' */
  Human_in_Loop_B.low_pass_e = 0.0;
  Human_in_Loop_B.low_pass_e += Human_in_Loop_P.low_pass_C_a *
    Human_in_Loop_X.low_pass_CSTATE_b[1];

  /* MATLAB Function: '<S41>/MVC' */
  /* MATLAB Function 'Sensor Data/EMG module/Preprocessing1/MVC': '<S47>:1' */
  /* '<S47>:1:3' */
  Human_in_Loop_B.y_g = (Human_in_Loop_B.low_pass_e -
    Human_in_Loop_P.MVC_r_M_GAS_MIN) / (Human_in_Loop_P.MVC_r_M_GAS_MAX -
    Human_in_Loop_P.MVC_r_M_GAS_MIN);

  /* StateSpace: '<S42>/low_pass' */
  Human_in_Loop_B.low_pass_m = 0.0;
  Human_in_Loop_B.low_pass_m += Human_in_Loop_P.low_pass_C_g *
    Human_in_Loop_X.low_pass_CSTATE_m[1];

  /* MATLAB Function: '<S42>/MVC' */
  /* MATLAB Function 'Sensor Data/EMG module/Preprocessing2/MVC': '<S48>:1' */
  /* '<S48>:1:3' */
  Human_in_Loop_B.y_j = (Human_in_Loop_B.low_pass_m -
    Human_in_Loop_P.MVC_r_L_GAS_MIN) / (Human_in_Loop_P.MVC_r_L_GAS_MAX -
    Human_in_Loop_P.MVC_r_L_GAS_MIN);

  /* StateSpace: '<S43>/low_pass' */
  Human_in_Loop_B.low_pass_d = 0.0;
  Human_in_Loop_B.low_pass_d += Human_in_Loop_P.low_pass_C_e *
    Human_in_Loop_X.low_pass_CSTATE_d[1];

  /* MATLAB Function: '<S43>/MVC' */
  /* MATLAB Function 'Sensor Data/EMG module/Preprocessing3/MVC': '<S49>:1' */
  /* '<S49>:1:3' */
  Human_in_Loop_B.y_e = (Human_in_Loop_B.low_pass_d -
    Human_in_Loop_P.MVC_l_M_GAS_MIN) / (Human_in_Loop_P.MVC_l_M_GAS_MAX -
    Human_in_Loop_P.MVC_l_M_GAS_MIN);

  /* StateSpace: '<S45>/low_pass' */
  Human_in_Loop_B.low_pass_p = 0.0;
  Human_in_Loop_B.low_pass_p += Human_in_Loop_P.low_pass_C_o *
    Human_in_Loop_X.low_pass_CSTATE_j[1];

  /* MATLAB Function: '<S45>/MVC' */
  /* MATLAB Function 'Sensor Data/EMG module/Preprocessing9/MVC': '<S50>:1' */
  /* '<S50>:1:3' */
  Human_in_Loop_B.y = (Human_in_Loop_B.low_pass_p -
                       Human_in_Loop_P.MVC_l_L_GAS_MIN) /
    (Human_in_Loop_P.MVC_l_L_GAS_MAX - Human_in_Loop_P.MVC_l_L_GAS_MIN);

  /* StateSpace: '<S44>/low_pass' */
  Human_in_Loop_B.low_pass_i = 0.0;
  Human_in_Loop_B.low_pass_i += Human_in_Loop_P.low_pass_C_p *
    Human_in_Loop_X.low_pass_CSTATE_a[1];

  /* MATLAB Function: '<S28>/Mux1' */
  /* MATLAB Function 'Sensor Data/EMG module/Mux1': '<S39>:1' */
  /* '<S39>:1:3' */
  Human_in_Loop_B.x_c[0] = Human_in_Loop_B.y_b;
  Human_in_Loop_B.x_c[1] = Human_in_Loop_B.y_g;
  Human_in_Loop_B.x_c[2] = Human_in_Loop_B.y_j;
  Human_in_Loop_B.x_c[3] = Human_in_Loop_B.y_e;
  Human_in_Loop_B.x_c[4] = Human_in_Loop_B.y;
  Human_in_Loop_B.x_c[5] = Human_in_Loop_B.low_pass_i;
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    for (i = 0; i < 6; i++) {
      /* Delay: '<S16>/Delay1' */
      Human_in_Loop_B.Delay1[i] = Human_in_Loop_DW.Delay1_DSTATE[i];

      /* Delay: '<S16>/Delay5' */
      Human_in_Loop_B.Delay5[i] = Human_in_Loop_DW.Delay5_DSTATE[i];

      /* Delay: '<S16>/Delay4' */
      Human_in_Loop_B.Delay4[i] = Human_in_Loop_DW.Delay4_DSTATE[i];

      /* Delay: '<S16>/Delay8' */
      Human_in_Loop_B.Delay8[i] = Human_in_Loop_DW.Delay8_DSTATE[i];
    }

    /* MATLAB Function: '<S16>/MATLAB Function' incorporates:
     *  Constant: '<S16>/Sample Time'
     */
    /* MATLAB Function 'ObjectiveCala/Single Cycle Analysis1/MATLAB Function': '<S20>:1' */
    /* '<S20>:1:28' */
    /* '<S20>:1:11' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_B.Max[i] = Human_in_Loop_B.Delay1[i];
    }

    /* '<S20>:1:12' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_B.Mean[i] = Human_in_Loop_B.Delay5[i];
    }

    /* '<S20>:1:13' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_B.RMS[i] = Human_in_Loop_B.Delay4[i];
    }

    /* '<S20>:1:14' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_B.iEMG[i] = Human_in_Loop_B.Delay8[i];
    }

    /* '<S20>:1:16' */
    torque_zero = Human_in_Loop_B.Delay7;

    /* '<S20>:1:17' */
    Human_in_Loop_B.Cycle_Frequency = Human_in_Loop_B.Delay6;

    /* '<S20>:1:19' */
    for (i = 0; i < 6; i++) {
      /* '<S20>:1:19' */
      /* '<S20>:1:20' */
      varargin_1[0] = Human_in_Loop_DW.EMG_Memory[i];
      varargin_1[1] = Human_in_Loop_B.x_c[i];
      ixstart = 1;
      B_0 = varargin_1[0];
      if (rtIsNaN(varargin_1[0])) {
        rankR = 2;
        exitg1 = false;
        while ((!exitg1) && (rankR < 3)) {
          ixstart = 2;
          if (!rtIsNaN(varargin_1[1])) {
            B_0 = varargin_1[1];
            exitg1 = true;
          } else {
            rankR = 3;
          }
        }
      }

      if ((ixstart < 2) && (varargin_1[1] > B_0)) {
        B_0 = varargin_1[1];
      }

      Human_in_Loop_DW.EMG_Memory[i] = B_0;

      /* '<S20>:1:21' */
      Human_in_Loop_DW.EMG_Memory[6 + i] += Human_in_Loop_B.x_c[i];

      /* '<S20>:1:22' */
      Human_in_Loop_DW.EMG_Memory[12 + i] += Human_in_Loop_B.x_c[i] *
        Human_in_Loop_B.x_c[i];

      /* '<S20>:1:23' */
      Human_in_Loop_DW.EMG_Memory[18 + i] += fabs(Human_in_Loop_B.x_c[i]);
    }

    if ((Human_in_Loop_B.x[1] == 1.0) && (Human_in_Loop_B.Delay2 == 0.0)) {
      /* '<S20>:1:26' */
      /* '<S20>:1:27' */
      for (ixstart = 0; ixstart < 6; ixstart++) {
        Human_in_Loop_B.Max[ixstart] = Human_in_Loop_DW.EMG_Memory[ixstart];
      }

      /* '<S20>:1:28' */
      for (ixstart = 0; ixstart < 6; ixstart++) {
        Human_in_Loop_B.Mean[ixstart] = Human_in_Loop_DW.EMG_Memory[6 + ixstart]
          / Human_in_Loop_B.Delay3;
      }

      /* '<S20>:1:29' */
      for (i = 0; i < 6; i++) {
        torque_zero = Human_in_Loop_DW.EMG_Memory[12 + i] /
          Human_in_Loop_B.Delay3;
        torque_zero = sqrt(torque_zero);
        x[i] = torque_zero;
      }

      for (ixstart = 0; ixstart < 6; ixstart++) {
        Human_in_Loop_B.RMS[ixstart] = x[ixstart];
      }

      /* '<S20>:1:30' */
      for (ixstart = 0; ixstart < 6; ixstart++) {
        Human_in_Loop_B.iEMG[ixstart] = Human_in_Loop_DW.EMG_Memory[18 + ixstart]
          / Human_in_Loop_B.Delay3;
      }

      /* '<S20>:1:32' */
      /* '<S20>:1:33' */
      /* '<S20>:1:34' */
      /* '<S20>:1:35' */
      for (ixstart = 0; ixstart < 6; ixstart++) {
        Human_in_Loop_DW.EMG_Memory[ixstart] = 0.0;
        Human_in_Loop_DW.EMG_Memory[6 + ixstart] = 0.0;
        Human_in_Loop_DW.EMG_Memory[12 + ixstart] = 0.0;
        Human_in_Loop_DW.EMG_Memory[18 + ixstart] = 0.0;
      }

      /* '<S20>:1:37' */
      torque_zero = Human_in_Loop_B.Delay3 *
        Human_in_Loop_P.SingleCycleAnalysis1_Ts;

      /* '<S20>:1:38' */
      Human_in_Loop_B.Cycle_Frequency = 1.0 / torque_zero;
    }

    Human_in_Loop_B.Cycle_Time = torque_zero;

    /* End of MATLAB Function: '<S16>/MATLAB Function' */

    /* Outputs for Triggered SubSystem: '<S15>/Mean Calculate' incorporates:
     *  TriggerPort: '<S18>/Trigger'
     */
    if (rtmIsMajorTimeStep(Human_in_Loop_M)) {
      zcEvent = rt_ZCFcn(RISING_ZERO_CROSSING,
                         &Human_in_Loop_PrevZCX.MeanCalculate_Trig_ZCE,
                         (Human_in_Loop_B.x[1]));
      if (zcEvent != NO_ZCEVENT) {
        /* UnitDelay: '<S18>/Unit Delay' */
        Human_in_Loop_B.UnitDelay_a = Human_in_Loop_DW.UnitDelay_DSTATE_g;

        /* UnitDelay: '<S18>/Unit Delay1' */
        memcpy(&Human_in_Loop_B.UnitDelay1_k[0],
               &Human_in_Loop_DW.UnitDelay1_DSTATE_e[0], 26U * sizeof(real_T));

        /* SignalConversion: '<S19>/TmpSignal ConversionAt SFunction Inport1' incorporates:
         *  MATLAB Function: '<S18>/MATLAB Function1'
         */
        Human_in_Loop_B.TmpSignalConversionAtSFunctionI[0] =
          Human_in_Loop_B.Cycle_Time;
        Human_in_Loop_B.TmpSignalConversionAtSFunctionI[1] =
          Human_in_Loop_B.Cycle_Frequency;
        for (i = 0; i < 6; i++) {
          Human_in_Loop_B.TmpSignalConversionAtSFunctionI[i + 2] =
            Human_in_Loop_B.Max[i];
        }

        for (i = 0; i < 6; i++) {
          Human_in_Loop_B.TmpSignalConversionAtSFunctionI[i + 8] =
            Human_in_Loop_B.Mean[i];
        }

        for (i = 0; i < 6; i++) {
          Human_in_Loop_B.TmpSignalConversionAtSFunctionI[i + 14] =
            Human_in_Loop_B.RMS[i];
        }

        for (i = 0; i < 6; i++) {
          Human_in_Loop_B.TmpSignalConversionAtSFunctionI[i + 20] =
            Human_in_Loop_B.iEMG[i];
        }

        /* End of SignalConversion: '<S19>/TmpSignal ConversionAt SFunction Inport1' */

        /* MATLAB Function: '<S18>/MATLAB Function1' incorporates:
         *  Constant: '<S18>/Constant'
         */
        /* MATLAB Function 'ObjectiveCala/Multi Cycle Analysis1/Mean Calculate/MATLAB Function1': '<S19>:1' */
        /* '<S19>:1:12' */
        /* '<S19>:1:6' */
        /* '<S19>:1:8' */
        torque_zero = Human_in_Loop_B.UnitDelay_a + 1.0;

        /* '<S19>:1:9' */
        ixstart = (int32_T)(Human_in_Loop_B.UnitDelay_a + 1.0);
        memcpy(&Human_in_Loop_B.Mean_c[0], &Human_in_Loop_B.UnitDelay1_k[0], 26U
               * sizeof(real_T));
        memcpy(&Human_in_Loop_DW.SingleCycleData[ixstart * 26 + -26],
               &Human_in_Loop_B.TmpSignalConversionAtSFunctionI[0], 26U * sizeof
               (real_T));
        if (Human_in_Loop_B.UnitDelay_a + 1.0 >=
            Human_in_Loop_P.MultiCycleAnalysis1_N) {
          /* '<S19>:1:11' */
          if (1.0 > Human_in_Loop_P.MultiCycleAnalysis1_N) {
            rankR = 0;
          } else {
            rankR = (int32_T)Human_in_Loop_P.MultiCycleAnalysis1_N;
          }

          /* '<S19>:1:12' */
          if (rankR == 0) {
            memset(&Human_in_Loop_B.Mean_c[0], 0, 26U * sizeof(real_T));
          } else {
            for (ixstart = 0; ixstart < 26; ixstart++) {
              Human_in_Loop_B.Mean_c[ixstart] =
                Human_in_Loop_DW.SingleCycleData[ixstart / 26 * 26 + ixstart %
                26];
            }

            for (i = 2; i <= rankR; i++) {
              xoffset = (i - 1) * 26;
              for (ixstart = 0; ixstart < 26; ixstart++) {
                c = xoffset + ixstart;
                Human_in_Loop_B.Mean_c[ixstart] +=
                  Human_in_Loop_DW.SingleCycleData[c / 26 * 26 + c % 26];
              }
            }
          }

          i = rankR;
          for (ixstart = 0; ixstart < 26; ixstart++) {
            Human_in_Loop_B.Mean_c[ixstart] /= (real_T)i;
          }

          /* '<S19>:1:13' */
          torque_zero = 0.0;
        }

        Human_in_Loop_B.count = torque_zero;

        /* Update for UnitDelay: '<S18>/Unit Delay' */
        Human_in_Loop_DW.UnitDelay_DSTATE_g = Human_in_Loop_B.count;

        /* Update for UnitDelay: '<S18>/Unit Delay1' */
        memcpy(&Human_in_Loop_DW.UnitDelay1_DSTATE_e[0],
               &Human_in_Loop_B.Mean_c[0], 26U * sizeof(real_T));
      }
    }

    /* End of Outputs for SubSystem: '<S15>/Mean Calculate' */

    /* MATLAB Function: '<S3>/torque_track_loss' */
    /* MATLAB Function 'ObjectiveCala/torque_track_loss': '<S17>:1' */
    if ((Human_in_Loop_DW.last_footstate_p == 0.0) && (Human_in_Loop_B.x[1] ==
         1.0)) {
      /* '<S17>:1:18' */
      torque_zero = floor(Human_in_Loop_B.x[2] * 500.0);
      if (1.0 > torque_zero) {
        c = 0;
      } else {
        c = (int32_T)torque_zero;
      }

      /* '<S17>:1:19' */
      tmp_size[0] = 1;
      tmp_size[1] = c;
      for (ixstart = 0; ixstart < c; ixstart++) {
        Human_in_Loop_B.tmp_data_mb[ixstart] =
          Human_in_Loop_DW.TorqueMem[ixstart << 2] - Human_in_Loop_DW.TorqueMem
          [(ixstart << 2) + 1];
      }

      Human_in_Loop_DW.loss_reg = Human_in_Loop_norm(Human_in_Loop_B.tmp_data_mb,
        tmp_size) / sqrt(Human_in_Loop_B.x[2] * 500.0);

      /* '<S17>:1:20' */
      memcpy(&tmp[0], &Human_in_Loop_DW.loss_mem[1], 9U * sizeof(real_T));
      tmp[9] = Human_in_Loop_DW.loss_reg;
      memcpy(&Human_in_Loop_DW.loss_mem[0], &tmp[0], 10U * sizeof(real_T));
    }

    /* '<S17>:1:23' */
    Human_in_Loop_DW.last_footstate_p = Human_in_Loop_B.x[1];

    /* '<S17>:1:24' */
    /* '<S17>:1:25' */
    Human_in_Loop_B.torque_track_rmse = Human_in_Loop_DW.loss_reg;
    Human_in_Loop_B.torque_track_rmse_mean = Human_in_Loop_mean
      (Human_in_Loop_DW.loss_mem);

    /* End of MATLAB Function: '<S3>/torque_track_loss' */

    /* Gain: '<S76>/Gain' */
    Human_in_Loop_B.Gain_h = Human_in_Loop_P.uOrderTD_Ts * Human_in_Loop_B.x2k1;

    /* Sum: '<S76>/Add' */
    Human_in_Loop_B.x1k = Human_in_Loop_B.Gain_h + Human_in_Loop_B.x1k1;

    /* MATLAB Function: '<S32>/MATLAB Function' */
    /* MATLAB Function 'Sensor Data/Torque module/MATLAB Function': '<S80>:1' */
    /* '<S80>:1:13' */
    /* '<S80>:1:8' */
    memcpy(&B[0], &Human_in_Loop_DW.data[1], 14U * sizeof(real_T));
    B[14] = Human_in_Loop_B.torque;
    memcpy(&Human_in_Loop_DW.data[0], &B[0], 15U * sizeof(real_T));

    /* '<S80>:1:13' */
    for (ixstart = 0; ixstart < 30; ixstart++) {
      A[ixstart] = b_A[ixstart];
    }

    Human_in_Loop_xgeqp3(A, varargin_1, tmp_size);
    rankR = 0;
    torque_zero = 15.0 * fabs(A[0]) * 2.2204460492503131E-16;
    while ((rankR < 2) && (!(fabs(A[15 * rankR + rankR]) <= torque_zero))) {
      rankR++;
    }

    x_0[0] = 0.0;
    x_0[1] = 0.0;
    memcpy(&B[0], &Human_in_Loop_DW.data[0], 15U * sizeof(real_T));
    if (varargin_1[0] != 0.0) {
      torque_zero = B[0];
      for (i = 1; i + 1 < 16; i++) {
        torque_zero += A[i] * B[i];
      }

      torque_zero *= varargin_1[0];
      if (torque_zero != 0.0) {
        B[0] -= torque_zero;
        for (i = 1; i + 1 < 16; i++) {
          B[i] -= A[i] * torque_zero;
        }
      }
    }

    if (varargin_1[1] != 0.0) {
      torque_zero = B[1];
      for (i = 2; i + 1 < 16; i++) {
        torque_zero += A[i + 15] * B[i];
      }

      torque_zero *= varargin_1[1];
      if (torque_zero != 0.0) {
        B[1] -= torque_zero;
        for (i = 2; i + 1 < 16; i++) {
          B[i] -= A[i + 15] * torque_zero;
        }
      }
    }

    for (i = 0; i + 1 <= rankR; i++) {
      x_0[tmp_size[i] - 1] = B[i];
    }

    for (ixstart = rankR - 1; ixstart + 1 > 0; ixstart--) {
      x_0[tmp_size[ixstart] - 1] /= A[15 * ixstart + ixstart];
      i = 1;
      while (i <= ixstart) {
        x_0[tmp_size[0] - 1] -= x_0[tmp_size[ixstart] - 1] * A[15 * ixstart];
        i = 2;
      }
    }

    /* '<S80>:1:15' */
    Human_in_Loop_B.torque_dot = x_0[1] * 5000.0;

    /* End of MATLAB Function: '<S32>/MATLAB Function' */

    /* S-Function (rti_commonblock): '<S78>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* SampleTimeMath: '<S55>/TSamp'
     *
     * About '<S55>/TSamp':
     *  y = u * K where K = 1 / ( w * Ts )
     */
    Human_in_Loop_B.TSamp = Human_in_Loop_B.Gain3_h * Human_in_Loop_P.TSamp_WtEt;

    /* UnitDelay: '<S55>/UD' */
    Human_in_Loop_B.Uk1 = Human_in_Loop_DW.UD_DSTATE;

    /* Sum: '<S55>/Diff' */
    Human_in_Loop_B.Diff = Human_in_Loop_B.TSamp - Human_in_Loop_B.Uk1;

    /* DiscreteFilter: '<S29>/0.4low1' */
    torque_zero = Human_in_Loop_B.Diff;
    torque_zero -= Human_in_Loop_P.u4low1_DenCoef[1] *
      Human_in_Loop_DW.u4low1_states[0];
    torque_zero -= Human_in_Loop_P.u4low1_DenCoef[2] *
      Human_in_Loop_DW.u4low1_states[1];
    torque_zero -= Human_in_Loop_P.u4low1_DenCoef[3] *
      Human_in_Loop_DW.u4low1_states[2];
    torque_zero /= Human_in_Loop_P.u4low1_DenCoef[0];
    Human_in_Loop_DW.u4low1_tmp = torque_zero;
    torque_zero = Human_in_Loop_P.u4low1_NumCoef[0] *
      Human_in_Loop_DW.u4low1_tmp;
    torque_zero += Human_in_Loop_P.u4low1_NumCoef[1] *
      Human_in_Loop_DW.u4low1_states[0];
    torque_zero += Human_in_Loop_P.u4low1_NumCoef[2] *
      Human_in_Loop_DW.u4low1_states[1];
    torque_zero += Human_in_Loop_P.u4low1_NumCoef[3] *
      Human_in_Loop_DW.u4low1_states[2];
    Human_in_Loop_B.u4low1 = torque_zero;

    /* UnitDelay: '<S51>/Unit Delay' */
    Human_in_Loop_B.UnitDelay = Human_in_Loop_DW.UnitDelay_DSTATE_d;

    /* UnitDelay: '<S51>/Unit Delay1' */
    Human_in_Loop_B.UnitDelay1 = Human_in_Loop_DW.UnitDelay1_DSTATE_f2;

    /* Sum: '<S51>/Add1' */
    Human_in_Loop_B.Add1 = Human_in_Loop_B.UnitDelay -
      Human_in_Loop_B.UnitDelay1;

    /* Gain: '<S51>/Gain' */
    Human_in_Loop_B.Gain_m = Human_in_Loop_P.uOrderTD_Ts_d *
      Human_in_Loop_B.Add1;

    /* Gain: '<S51>/Gain1' */
    B_0 = Human_in_Loop_P.uOrderTD_T;
    torque_zero = 1.0 / B_0;
    Human_in_Loop_B.Gain1_b = torque_zero * Human_in_Loop_B.Gain_m;

    /* Sum: '<S51>/Add2' */
    Human_in_Loop_B.Add2_h = Human_in_Loop_B.UnitDelay1 +
      Human_in_Loop_B.Gain1_b;

    /* Sum: '<S51>/Add3' */
    Human_in_Loop_B.Add3 = Human_in_Loop_B.Gain3_h - Human_in_Loop_B.Add2_h;

    /* Gain: '<S51>/Gain2' */
    B_0 = Human_in_Loop_P.uOrderTD_T;
    torque_zero = 1.0 / B_0;
    Human_in_Loop_B.Gain2_i = torque_zero * Human_in_Loop_B.Add3;

    /* Gain: '<S52>/Gain' */
    Human_in_Loop_B.Gain_o = Human_in_Loop_P.uOrderTD_Ts_n *
      Human_in_Loop_B.x2k1_d;

    /* Sum: '<S52>/Add' */
    Human_in_Loop_B.x1k_k = Human_in_Loop_B.Gain_o + Human_in_Loop_B.x1k1_e;

    /* MATLAB Function: '<S29>/Mux2' */
    /* MATLAB Function 'Sensor Data/Encoder module/Mux2': '<S62>:1' */
    /* '<S62>:1:3' */
    Human_in_Loop_B.x_i[0] = Human_in_Loop_B.u4low1;
    Human_in_Loop_B.x_i[1] = Human_in_Loop_B.Gain2_i;
    Human_in_Loop_B.x_i[2] = Human_in_Loop_B.x2k_a;

    /* S-Function (rti_commonblock): '<S58>/S-Function1' incorporates:
     *  Constant: '<S29>/VCC1'
     */
    /* This comment workarounds a code generation problem */

    /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup1 --- */
    /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 15 --- */
    {
      /* define variables required for BitOut realtime functions */
      UInt32 outputData = 0;

      /* write output state value to digital output channel 15-15 on port 3 */
      outputData = ( ( ( (UInt32)Human_in_Loop_P.VCC1_Value) << (15 - 1)) |
                    outputData);
      DioCl1DigOut_setChMaskOutData(pRTIDioC1DigOut_Port_3_Ch_15, outputData);
      DioCl1DigOut_write(pRTIDioC1DigOut_Port_3_Ch_15);
    }

    /* S-Function (rti_commonblock): '<S59>/S-Function1' incorporates:
     *  Constant: '<S29>/VCC3'
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

    /* S-Function (rti_commonblock): '<S64>/S-Function1' incorporates:
     *  Constant: '<S30>/Constant'
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

    /* S-Function (rti_commonblock): '<S33>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* --- Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL1 --- */
    /* --- [RTI120X, ADC C1] - Channel: 1 --- */
    {
      UInt32 readStateFlag[] = { 1 };

      /* define variable required for adc cl1 realtime functions */
      UInt32 IsNew = 0;

      /* wait until conversion result is available */
      while (IsNew == 0) {
        AdcCl1AnalogIn_isDataReady(pRTIAdcC1AnalogIn_Ch_1, &IsNew);
      }

      /* read conversion result from hardware */
      AdcCl1AnalogIn_read(pRTIAdcC1AnalogIn_Ch_1);

      /* retrieve conversion result */
      AdcCl1AnalogIn_getSingleValue(pRTIAdcC1AnalogIn_Ch_1, readStateFlag,
        (real_T*) &Human_in_Loop_B.SFunction1_g);
    }

    /* S-Function (rti_commonblock): '<S34>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* --- Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL2 --- */
    /* --- [RTI120X, ADC C1] - Channel: 2 --- */
    {
      UInt32 readStateFlag[] = { 1 };

      /* define variable required for adc cl1 realtime functions */
      UInt32 IsNew = 0;

      /* wait until conversion result is available */
      while (IsNew == 0) {
        AdcCl1AnalogIn_isDataReady(pRTIAdcC1AnalogIn_Ch_2, &IsNew);
      }

      /* read conversion result from hardware */
      AdcCl1AnalogIn_read(pRTIAdcC1AnalogIn_Ch_2);

      /* retrieve conversion result */
      AdcCl1AnalogIn_getSingleValue(pRTIAdcC1AnalogIn_Ch_2, readStateFlag,
        (real_T*) &Human_in_Loop_B.SFunction1_i);
    }

    /* S-Function (rti_commonblock): '<S35>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* --- Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL3 --- */
    /* --- [RTI120X, ADC C1] - Channel: 3 --- */
    {
      UInt32 readStateFlag[] = { 1 };

      /* define variable required for adc cl1 realtime functions */
      UInt32 IsNew = 0;

      /* wait until conversion result is available */
      while (IsNew == 0) {
        AdcCl1AnalogIn_isDataReady(pRTIAdcC1AnalogIn_Ch_3, &IsNew);
      }

      /* read conversion result from hardware */
      AdcCl1AnalogIn_read(pRTIAdcC1AnalogIn_Ch_3);

      /* retrieve conversion result */
      AdcCl1AnalogIn_getSingleValue(pRTIAdcC1AnalogIn_Ch_3, readStateFlag,
        (real_T*) &Human_in_Loop_B.SFunction1_m);
    }

    /* S-Function (rti_commonblock): '<S36>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* --- Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL4 --- */
    /* --- [RTI120X, ADC C1] - Channel: 4 --- */
    {
      UInt32 readStateFlag[] = { 1 };

      /* define variable required for adc cl1 realtime functions */
      UInt32 IsNew = 0;

      /* wait until conversion result is available */
      while (IsNew == 0) {
        AdcCl1AnalogIn_isDataReady(pRTIAdcC1AnalogIn_Ch_4, &IsNew);
      }

      /* read conversion result from hardware */
      AdcCl1AnalogIn_read(pRTIAdcC1AnalogIn_Ch_4);

      /* retrieve conversion result */
      AdcCl1AnalogIn_getSingleValue(pRTIAdcC1AnalogIn_Ch_4, readStateFlag,
        (real_T*) &Human_in_Loop_B.SFunction1_j);
    }

    /* S-Function (rti_commonblock): '<S37>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* --- Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL5 --- */
    /* --- [RTI120X, ADC C1] - Channel: 9 --- */
    {
      UInt32 readStateFlag[] = { 1 };

      /* define variable required for adc cl1 realtime functions */
      UInt32 IsNew = 0;

      /* wait until conversion result is available */
      while (IsNew == 0) {
        AdcCl1AnalogIn_isDataReady(pRTIAdcC1AnalogIn_Ch_9, &IsNew);
      }

      /* read conversion result from hardware */
      AdcCl1AnalogIn_read(pRTIAdcC1AnalogIn_Ch_9);

      /* retrieve conversion result */
      AdcCl1AnalogIn_getSingleValue(pRTIAdcC1AnalogIn_Ch_9, readStateFlag,
        (real_T*) &Human_in_Loop_B.SFunction1_l);
    }

    /* S-Function (rti_commonblock): '<S38>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* --- Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL6 --- */
    /* --- [RTI120X, ADC C1] - Channel: 10 --- */
    {
      UInt32 readStateFlag[] = { 1 };

      /* define variable required for adc cl1 realtime functions */
      UInt32 IsNew = 0;

      /* wait until conversion result is available */
      while (IsNew == 0) {
        AdcCl1AnalogIn_isDataReady(pRTIAdcC1AnalogIn_Ch_10, &IsNew);
      }

      /* read conversion result from hardware */
      AdcCl1AnalogIn_read(pRTIAdcC1AnalogIn_Ch_10);

      /* retrieve conversion result */
      AdcCl1AnalogIn_getSingleValue(pRTIAdcC1AnalogIn_Ch_10, readStateFlag,
        (real_T*) &Human_in_Loop_B.SFunction1_k);
    }

    /* Gain: '<S28>/Gain' */
    Human_in_Loop_B.Gain_mi = Human_in_Loop_P.EMGmodule_Kemg *
      Human_in_Loop_B.SFunction1_g;

    /* Gain: '<S28>/Gain1' */
    Human_in_Loop_B.Gain1_a = Human_in_Loop_P.EMGmodule_Kemg *
      Human_in_Loop_B.SFunction1_i;

    /* Gain: '<S28>/Gain2' */
    Human_in_Loop_B.Gain2_p = Human_in_Loop_P.EMGmodule_Kemg *
      Human_in_Loop_B.SFunction1_m;

    /* Gain: '<S28>/Gain3' */
    Human_in_Loop_B.Gain3_n = Human_in_Loop_P.EMGmodule_Kemg *
      Human_in_Loop_B.SFunction1_j;

    /* Gain: '<S28>/Gain4' */
    Human_in_Loop_B.Gain4_m = Human_in_Loop_P.EMGmodule_Kemg *
      Human_in_Loop_B.SFunction1_l;

    /* Gain: '<S28>/Gain5' */
    Human_in_Loop_B.Gain5 = Human_in_Loop_P.EMGmodule_Kemg *
      Human_in_Loop_B.SFunction1_k;
  }

  /* StateSpace: '<S40>/high_pass' */
  Human_in_Loop_B.high_pass = 0.0;
  Human_in_Loop_B.high_pass += Human_in_Loop_P.high_pass_C[0] *
    Human_in_Loop_X.high_pass_CSTATE[0];
  Human_in_Loop_B.high_pass += Human_in_Loop_P.high_pass_C[1] *
    Human_in_Loop_X.high_pass_CSTATE[1];
  Human_in_Loop_B.high_pass += Human_in_Loop_P.high_pass_D *
    Human_in_Loop_B.Gain_mi;

  /* Abs: '<S40>/Abs' */
  Human_in_Loop_B.Abs = fabs(Human_in_Loop_B.high_pass);

  /* StateSpace: '<S41>/high_pass' */
  Human_in_Loop_B.high_pass_h = 0.0;
  Human_in_Loop_B.high_pass_h += Human_in_Loop_P.high_pass_C_b[0] *
    Human_in_Loop_X.high_pass_CSTATE_p[0];
  Human_in_Loop_B.high_pass_h += Human_in_Loop_P.high_pass_C_b[1] *
    Human_in_Loop_X.high_pass_CSTATE_p[1];
  Human_in_Loop_B.high_pass_h += Human_in_Loop_P.high_pass_D_d *
    Human_in_Loop_B.Gain1_a;

  /* Abs: '<S41>/Abs' */
  Human_in_Loop_B.Abs_i = fabs(Human_in_Loop_B.high_pass_h);

  /* StateSpace: '<S42>/high_pass' */
  Human_in_Loop_B.high_pass_hy = 0.0;
  Human_in_Loop_B.high_pass_hy += Human_in_Loop_P.high_pass_C_by[0] *
    Human_in_Loop_X.high_pass_CSTATE_f[0];
  Human_in_Loop_B.high_pass_hy += Human_in_Loop_P.high_pass_C_by[1] *
    Human_in_Loop_X.high_pass_CSTATE_f[1];
  Human_in_Loop_B.high_pass_hy += Human_in_Loop_P.high_pass_D_o *
    Human_in_Loop_B.Gain2_p;

  /* Abs: '<S42>/Abs' */
  Human_in_Loop_B.Abs_l = fabs(Human_in_Loop_B.high_pass_hy);

  /* StateSpace: '<S43>/high_pass' */
  Human_in_Loop_B.high_pass_j = 0.0;
  Human_in_Loop_B.high_pass_j += Human_in_Loop_P.high_pass_C_p[0] *
    Human_in_Loop_X.high_pass_CSTATE_a[0];
  Human_in_Loop_B.high_pass_j += Human_in_Loop_P.high_pass_C_p[1] *
    Human_in_Loop_X.high_pass_CSTATE_a[1];
  Human_in_Loop_B.high_pass_j += Human_in_Loop_P.high_pass_D_n *
    Human_in_Loop_B.Gain3_n;

  /* Abs: '<S43>/Abs' */
  Human_in_Loop_B.Abs_b = fabs(Human_in_Loop_B.high_pass_j);

  /* StateSpace: '<S44>/high_pass' */
  Human_in_Loop_B.high_pass_e = 0.0;
  Human_in_Loop_B.high_pass_e += Human_in_Loop_P.high_pass_C_e[0] *
    Human_in_Loop_X.high_pass_CSTATE_o[0];
  Human_in_Loop_B.high_pass_e += Human_in_Loop_P.high_pass_C_e[1] *
    Human_in_Loop_X.high_pass_CSTATE_o[1];
  Human_in_Loop_B.high_pass_e += Human_in_Loop_P.high_pass_D_l *
    Human_in_Loop_B.Gain5;

  /* Abs: '<S44>/Abs' */
  Human_in_Loop_B.Abs_lh = fabs(Human_in_Loop_B.high_pass_e);

  /* StateSpace: '<S45>/high_pass' */
  Human_in_Loop_B.high_pass_b = 0.0;
  Human_in_Loop_B.high_pass_b += Human_in_Loop_P.high_pass_C_h[0] *
    Human_in_Loop_X.high_pass_CSTATE_ff[0];
  Human_in_Loop_B.high_pass_b += Human_in_Loop_P.high_pass_C_h[1] *
    Human_in_Loop_X.high_pass_CSTATE_ff[1];
  Human_in_Loop_B.high_pass_b += Human_in_Loop_P.high_pass_D_e *
    Human_in_Loop_B.Gain4_m;

  /* Abs: '<S45>/Abs' */
  Human_in_Loop_B.Abs_p = fabs(Human_in_Loop_B.high_pass_b);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S72>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:3 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X3]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X3]->processed) {
        can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X3]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X3]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X3]->timestamp > 0.0) {
        if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X3]->processed) {
          Human_in_Loop_B.SFunction1_o4 = (real_T)
            can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X3]->processed;
          Human_in_Loop_B.SFunction1_o5 = (real_T)
            can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X3]->timestamp;
          CAN_Msg = can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X3]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "Accel_x" (0|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            Human_in_Loop_B.SFunction1_o1_c = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "Accel_y" (16|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            Human_in_Loop_B.SFunction1_o2_l = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "Accel_z" (32|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            Human_in_Loop_B.SFunction1_o3 = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }

      if (!can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X3]->processed) {
        /* ... set RX status to 0 because no new message has arrived */
        Human_in_Loop_B.SFunction1_o4 = 0.0;
      }
    }

    /* DataTypeConversion: '<S73>/Data Type Conversion' */
    Human_in_Loop_B.DataTypeConversion = Human_in_Loop_B.SFunction1_o1_c;

    /* Gain: '<S73>/Gain' */
    Human_in_Loop_B.Gain_gd = Human_in_Loop_P.Gain_Gain_a *
      Human_in_Loop_B.DataTypeConversion;

    /* DataTypeConversion: '<S74>/Data Type Conversion' */
    Human_in_Loop_B.DataTypeConversion_j = Human_in_Loop_B.SFunction1_o2_l;

    /* Gain: '<S74>/Gain' */
    Human_in_Loop_B.Gain_b = Human_in_Loop_P.Gain_Gain_h *
      Human_in_Loop_B.DataTypeConversion_j;

    /* DataTypeConversion: '<S75>/Data Type Conversion' */
    Human_in_Loop_B.DataTypeConversion_m = Human_in_Loop_B.SFunction1_o3;

    /* Gain: '<S75>/Gain' */
    Human_in_Loop_B.Gain_n = Human_in_Loop_P.Gain_Gain_jd *
      Human_in_Loop_B.DataTypeConversion_m;

    /* S-Function (rti_commonblock): '<S66>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* S-Function (rti_commonblock): '<S67>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* S-Function (rti_commonblock): '<S70>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:100 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64]->processed) {
        can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64]->timestamp > 0.0) {
        if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64]->processed) {
          Human_in_Loop_B.SFunction1_o2_f = (real_T)
            can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64]->processed;
          Human_in_Loop_B.SFunction1_o3_g = (real_T)
            can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64]->timestamp;
          CAN_Msg = can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "RX signal" (0|32, standard signal, IEEE 32, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.SgnBytes.Byte2 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte3 = CAN_Msg[3];
            Human_in_Loop_B.SFunction1_o1_h = ((real_T) CAN_Sgn.IeeeVal32);
          }
        }
      }

      if (!can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64]->processed) {
        /* ... set RX status to 0 because no new message has arrived */
        Human_in_Loop_B.SFunction1_o2_f = 0.0;
      }
    }
  }

  /* Sin: '<S68>/Sine Wave' */
  Human_in_Loop_B.SineWave = sin(Human_in_Loop_P.SineWave_Freq *
    Human_in_Loop_M->Timing.t[0] + Human_in_Loop_P.SineWave_Phase) *
    Human_in_Loop_P.SineWave_Amp + Human_in_Loop_P.SineWave_Bias;
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S71>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN TX Message Block: "TX Message" Id:100 */
    {
      UInt32 CAN_Msg[8] = { 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };

      Float32 delayTime = 0.0;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->processed) {
        can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->timestamp > 0.0) {
        Human_in_Loop_B.SFunction1_o1_l = (real_T)
          can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->processed;
        Human_in_Loop_B.SFunction1_o2_jb = (real_T)
          can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->timestamp;
        Human_in_Loop_B.SFunction1_o3_p = (real_T)
          can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->deltatime;
        Human_in_Loop_B.SFunction1_o4_k = (real_T)
          can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->delaytime;
      }

      /* ... Encode Simulink signals of TX and RM blocks*/
      {
        rtican_Signal_t CAN_Sgn;

        /* ...... "TX signal" (0|32, standard signal, IEEE 32, little endian) */
        /* Add or substract 0.5 in order to round to nearest integer */
        CAN_Sgn.IeeeVal32 = (Float32) (( Human_in_Loop_B.SineWave ));
        CAN_Msg[0] |= CAN_Sgn.SgnBytes.Byte0;
        CAN_Msg[1] |= CAN_Sgn.SgnBytes.Byte1;
        CAN_Msg[2] |= CAN_Sgn.SgnBytes.Byte2;
        CAN_Msg[3] |= CAN_Sgn.SgnBytes.Byte3;
      }

      /* ... Write the data to the CAN microcontroller and trigger the sending of the message */
      can_tp1_msg_send(can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64], 8, &(CAN_Msg
        [0]), delayTime);
    }
  }

  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[2] == 0) {
    /* S-Function (rti_commonblock): '<S2>/S-Function1' */
    /* This comment workarounds a code generation problem */
  }

  /* user code (Output function Trailer for TID0) */

  /* dSPACE Data Capture block: Human_in_Loop/Data Capture1 */
  /* ... Service number: 1 */
  /* ... Service name:   500Hz */
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[2] == 0) {
    DsDaq_Service(0, 0, 1, (DsDaqSTimestampStruct *)
                  rtk_current_task_absolute_time_ptr_get());
  }
}

/* Model update function */
void Human_in_Loop_update(void)
{
  int32_T i;
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* Update for DiscreteFilter: '<S32>/0.4low2' */
    Human_in_Loop_DW.u4low2_states[2] = Human_in_Loop_DW.u4low2_states[1];
    Human_in_Loop_DW.u4low2_states[1] = Human_in_Loop_DW.u4low2_states[0];
    Human_in_Loop_DW.u4low2_states[0] = Human_in_Loop_DW.u4low2_tmp;

    /* Update for UnitDelay: '<S76>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE = Human_in_Loop_B.x2k;

    /* Update for UnitDelay: '<S76>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE = Human_in_Loop_B.x1k;

    /* Update for UnitDelay: '<S76>/Unit Delay2' */
    Human_in_Loop_DW.UnitDelay2_DSTATE = Human_in_Loop_B.torque;

    /* Update for UnitDelay: '<S52>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE_f = Human_in_Loop_B.x2k_a;

    /* Update for UnitDelay: '<S52>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_a = Human_in_Loop_B.x1k_k;

    /* Update for UnitDelay: '<S52>/Unit Delay2' */
    Human_in_Loop_DW.UnitDelay2_DSTATE_e = Human_in_Loop_B.Gain3_h;

    /* Update for Delay: '<S16>/Delay3' */
    Human_in_Loop_DW.Delay3_DSTATE = Human_in_Loop_B.Gain_i;

    /* Update for Delay: '<S16>/Delay2' */
    Human_in_Loop_DW.Delay2_DSTATE = Human_in_Loop_B.x[1];

    /* Update for Delay: '<S16>/Delay7' */
    Human_in_Loop_DW.Delay7_DSTATE = Human_in_Loop_B.Cycle_Time;

    /* Update for Delay: '<S16>/Delay6' */
    Human_in_Loop_DW.Delay6_DSTATE = Human_in_Loop_B.Cycle_Frequency;

    /* Update for Delay: '<S16>/Delay1' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_DW.Delay1_DSTATE[i] = Human_in_Loop_B.Max[i];
    }

    /* End of Update for Delay: '<S16>/Delay1' */

    /* Update for Delay: '<S16>/Delay5' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_DW.Delay5_DSTATE[i] = Human_in_Loop_B.Mean[i];
    }

    /* End of Update for Delay: '<S16>/Delay5' */

    /* Update for Delay: '<S16>/Delay4' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_DW.Delay4_DSTATE[i] = Human_in_Loop_B.RMS[i];
    }

    /* End of Update for Delay: '<S16>/Delay4' */

    /* Update for Delay: '<S16>/Delay8' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_DW.Delay8_DSTATE[i] = Human_in_Loop_B.iEMG[i];
    }

    /* End of Update for Delay: '<S16>/Delay8' */

    /* Update for UnitDelay: '<S55>/UD' */
    Human_in_Loop_DW.UD_DSTATE = Human_in_Loop_B.TSamp;

    /* Update for DiscreteFilter: '<S29>/0.4low1' */
    Human_in_Loop_DW.u4low1_states[2] = Human_in_Loop_DW.u4low1_states[1];
    Human_in_Loop_DW.u4low1_states[1] = Human_in_Loop_DW.u4low1_states[0];
    Human_in_Loop_DW.u4low1_states[0] = Human_in_Loop_DW.u4low1_tmp;

    /* Update for UnitDelay: '<S51>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_d = Human_in_Loop_B.Gain3_h;

    /* Update for UnitDelay: '<S51>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE_f2 = Human_in_Loop_B.Add2_h;
  }

  if (rtmIsMajorTimeStep(Human_in_Loop_M)) {
    rt_ertODEUpdateContinuousStates(&Human_in_Loop_M->solverInfo);
  }

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

  Human_in_Loop_M->Timing.t[0] = rtsiGetSolverStopTime
    (&Human_in_Loop_M->solverInfo);

  {
    /* Update absolute timer for sample time: [0.0002s, 0.0s] */
    /* The "clockTick1" counts the number of times the code of this task has
     * been executed. The resolution of this integer timer is 0.0002, which is the step size
     * of the task. Size of "clockTick1" ensures timer will not overflow during the
     * application lifespan selected.
     * Timer of this task consists of two 32 bit unsigned integers.
     * The two integers represent the low bits Timing.clockTick1 and the high bits
     * Timing.clockTickH1. When the low bit overflows to 0, the high bits increment.
     */
    Human_in_Loop_M->Timing.clockTick1++;
    if (!Human_in_Loop_M->Timing.clockTick1) {
      Human_in_Loop_M->Timing.clockTickH1++;
    }
  }

  rate_scheduler();
}

/* Derivatives for root system: '<Root>' */
void Human_in_Loop_derivatives(void)
{
  XDot_Human_in_Loop_T *_rtXdot;
  _rtXdot = ((XDot_Human_in_Loop_T *) Human_in_Loop_M->derivs);

  /* Derivatives for StateSpace: '<S40>/low_pass' */
  _rtXdot->low_pass_CSTATE[0] = 0.0;
  _rtXdot->low_pass_CSTATE[1] = 0.0;
  _rtXdot->low_pass_CSTATE[0] += Human_in_Loop_P.low_pass_A[0] *
    Human_in_Loop_X.low_pass_CSTATE[0];
  _rtXdot->low_pass_CSTATE[0] += Human_in_Loop_P.low_pass_A[1] *
    Human_in_Loop_X.low_pass_CSTATE[1];
  _rtXdot->low_pass_CSTATE[1] += Human_in_Loop_P.low_pass_A[2] *
    Human_in_Loop_X.low_pass_CSTATE[0];
  _rtXdot->low_pass_CSTATE[0] += Human_in_Loop_P.low_pass_B *
    Human_in_Loop_B.Abs;

  /* Derivatives for StateSpace: '<S41>/low_pass' */
  _rtXdot->low_pass_CSTATE_b[0] = 0.0;
  _rtXdot->low_pass_CSTATE_b[1] = 0.0;
  _rtXdot->low_pass_CSTATE_b[0] += Human_in_Loop_P.low_pass_A_k[0] *
    Human_in_Loop_X.low_pass_CSTATE_b[0];
  _rtXdot->low_pass_CSTATE_b[0] += Human_in_Loop_P.low_pass_A_k[1] *
    Human_in_Loop_X.low_pass_CSTATE_b[1];
  _rtXdot->low_pass_CSTATE_b[1] += Human_in_Loop_P.low_pass_A_k[2] *
    Human_in_Loop_X.low_pass_CSTATE_b[0];
  _rtXdot->low_pass_CSTATE_b[0] += Human_in_Loop_P.low_pass_B_j *
    Human_in_Loop_B.Abs_i;

  /* Derivatives for StateSpace: '<S42>/low_pass' */
  _rtXdot->low_pass_CSTATE_m[0] = 0.0;
  _rtXdot->low_pass_CSTATE_m[1] = 0.0;
  _rtXdot->low_pass_CSTATE_m[0] += Human_in_Loop_P.low_pass_A_e[0] *
    Human_in_Loop_X.low_pass_CSTATE_m[0];
  _rtXdot->low_pass_CSTATE_m[0] += Human_in_Loop_P.low_pass_A_e[1] *
    Human_in_Loop_X.low_pass_CSTATE_m[1];
  _rtXdot->low_pass_CSTATE_m[1] += Human_in_Loop_P.low_pass_A_e[2] *
    Human_in_Loop_X.low_pass_CSTATE_m[0];
  _rtXdot->low_pass_CSTATE_m[0] += Human_in_Loop_P.low_pass_B_f *
    Human_in_Loop_B.Abs_l;

  /* Derivatives for StateSpace: '<S43>/low_pass' */
  _rtXdot->low_pass_CSTATE_d[0] = 0.0;
  _rtXdot->low_pass_CSTATE_d[1] = 0.0;
  _rtXdot->low_pass_CSTATE_d[0] += Human_in_Loop_P.low_pass_A_l[0] *
    Human_in_Loop_X.low_pass_CSTATE_d[0];
  _rtXdot->low_pass_CSTATE_d[0] += Human_in_Loop_P.low_pass_A_l[1] *
    Human_in_Loop_X.low_pass_CSTATE_d[1];
  _rtXdot->low_pass_CSTATE_d[1] += Human_in_Loop_P.low_pass_A_l[2] *
    Human_in_Loop_X.low_pass_CSTATE_d[0];
  _rtXdot->low_pass_CSTATE_d[0] += Human_in_Loop_P.low_pass_B_n *
    Human_in_Loop_B.Abs_b;

  /* Derivatives for StateSpace: '<S45>/low_pass' */
  _rtXdot->low_pass_CSTATE_j[0] = 0.0;
  _rtXdot->low_pass_CSTATE_j[1] = 0.0;
  _rtXdot->low_pass_CSTATE_j[0] += Human_in_Loop_P.low_pass_A_c[0] *
    Human_in_Loop_X.low_pass_CSTATE_j[0];
  _rtXdot->low_pass_CSTATE_j[0] += Human_in_Loop_P.low_pass_A_c[1] *
    Human_in_Loop_X.low_pass_CSTATE_j[1];
  _rtXdot->low_pass_CSTATE_j[1] += Human_in_Loop_P.low_pass_A_c[2] *
    Human_in_Loop_X.low_pass_CSTATE_j[0];
  _rtXdot->low_pass_CSTATE_j[0] += Human_in_Loop_P.low_pass_B_m *
    Human_in_Loop_B.Abs_p;

  /* Derivatives for StateSpace: '<S44>/low_pass' */
  _rtXdot->low_pass_CSTATE_a[0] = 0.0;
  _rtXdot->low_pass_CSTATE_a[1] = 0.0;
  _rtXdot->low_pass_CSTATE_a[0] += Human_in_Loop_P.low_pass_A_ec[0] *
    Human_in_Loop_X.low_pass_CSTATE_a[0];
  _rtXdot->low_pass_CSTATE_a[0] += Human_in_Loop_P.low_pass_A_ec[1] *
    Human_in_Loop_X.low_pass_CSTATE_a[1];
  _rtXdot->low_pass_CSTATE_a[1] += Human_in_Loop_P.low_pass_A_ec[2] *
    Human_in_Loop_X.low_pass_CSTATE_a[0];
  _rtXdot->low_pass_CSTATE_a[0] += Human_in_Loop_P.low_pass_B_h *
    Human_in_Loop_B.Abs_lh;

  /* Derivatives for StateSpace: '<S40>/high_pass' */
  _rtXdot->high_pass_CSTATE[0] = 0.0;

  /* Derivatives for StateSpace: '<S41>/high_pass' */
  _rtXdot->high_pass_CSTATE_p[0] = 0.0;

  /* Derivatives for StateSpace: '<S42>/high_pass' */
  _rtXdot->high_pass_CSTATE_f[0] = 0.0;

  /* Derivatives for StateSpace: '<S43>/high_pass' */
  _rtXdot->high_pass_CSTATE_a[0] = 0.0;

  /* Derivatives for StateSpace: '<S44>/high_pass' */
  _rtXdot->high_pass_CSTATE_o[0] = 0.0;

  /* Derivatives for StateSpace: '<S45>/high_pass' */
  _rtXdot->high_pass_CSTATE_ff[0] = 0.0;

  /* Derivatives for StateSpace: '<S40>/high_pass' */
  _rtXdot->high_pass_CSTATE[1] = 0.0;

  /* Derivatives for StateSpace: '<S41>/high_pass' */
  _rtXdot->high_pass_CSTATE_p[1] = 0.0;

  /* Derivatives for StateSpace: '<S42>/high_pass' */
  _rtXdot->high_pass_CSTATE_f[1] = 0.0;

  /* Derivatives for StateSpace: '<S43>/high_pass' */
  _rtXdot->high_pass_CSTATE_a[1] = 0.0;

  /* Derivatives for StateSpace: '<S44>/high_pass' */
  _rtXdot->high_pass_CSTATE_o[1] = 0.0;

  /* Derivatives for StateSpace: '<S45>/high_pass' */
  _rtXdot->high_pass_CSTATE_ff[1] = 0.0;

  /* Derivatives for StateSpace: '<S40>/high_pass' */
  _rtXdot->high_pass_CSTATE[0] += Human_in_Loop_P.high_pass_A[0] *
    Human_in_Loop_X.high_pass_CSTATE[1];
  _rtXdot->high_pass_CSTATE[1] += Human_in_Loop_P.high_pass_A[1] *
    Human_in_Loop_X.high_pass_CSTATE[0];
  _rtXdot->high_pass_CSTATE[1] += Human_in_Loop_P.high_pass_A[2] *
    Human_in_Loop_X.high_pass_CSTATE[1];
  _rtXdot->high_pass_CSTATE[1] += Human_in_Loop_P.high_pass_B *
    Human_in_Loop_B.Gain_mi;

  /* Derivatives for StateSpace: '<S41>/high_pass' */
  _rtXdot->high_pass_CSTATE_p[0] += Human_in_Loop_P.high_pass_A_g[0] *
    Human_in_Loop_X.high_pass_CSTATE_p[1];
  _rtXdot->high_pass_CSTATE_p[1] += Human_in_Loop_P.high_pass_A_g[1] *
    Human_in_Loop_X.high_pass_CSTATE_p[0];
  _rtXdot->high_pass_CSTATE_p[1] += Human_in_Loop_P.high_pass_A_g[2] *
    Human_in_Loop_X.high_pass_CSTATE_p[1];
  _rtXdot->high_pass_CSTATE_p[1] += Human_in_Loop_P.high_pass_B_i *
    Human_in_Loop_B.Gain1_a;

  /* Derivatives for StateSpace: '<S42>/high_pass' */
  _rtXdot->high_pass_CSTATE_f[0] += Human_in_Loop_P.high_pass_A_p[0] *
    Human_in_Loop_X.high_pass_CSTATE_f[1];
  _rtXdot->high_pass_CSTATE_f[1] += Human_in_Loop_P.high_pass_A_p[1] *
    Human_in_Loop_X.high_pass_CSTATE_f[0];
  _rtXdot->high_pass_CSTATE_f[1] += Human_in_Loop_P.high_pass_A_p[2] *
    Human_in_Loop_X.high_pass_CSTATE_f[1];
  _rtXdot->high_pass_CSTATE_f[1] += Human_in_Loop_P.high_pass_B_h *
    Human_in_Loop_B.Gain2_p;

  /* Derivatives for StateSpace: '<S43>/high_pass' */
  _rtXdot->high_pass_CSTATE_a[0] += Human_in_Loop_P.high_pass_A_a[0] *
    Human_in_Loop_X.high_pass_CSTATE_a[1];
  _rtXdot->high_pass_CSTATE_a[1] += Human_in_Loop_P.high_pass_A_a[1] *
    Human_in_Loop_X.high_pass_CSTATE_a[0];
  _rtXdot->high_pass_CSTATE_a[1] += Human_in_Loop_P.high_pass_A_a[2] *
    Human_in_Loop_X.high_pass_CSTATE_a[1];
  _rtXdot->high_pass_CSTATE_a[1] += Human_in_Loop_P.high_pass_B_n *
    Human_in_Loop_B.Gain3_n;

  /* Derivatives for StateSpace: '<S44>/high_pass' */
  _rtXdot->high_pass_CSTATE_o[0] += Human_in_Loop_P.high_pass_A_a0[0] *
    Human_in_Loop_X.high_pass_CSTATE_o[1];
  _rtXdot->high_pass_CSTATE_o[1] += Human_in_Loop_P.high_pass_A_a0[1] *
    Human_in_Loop_X.high_pass_CSTATE_o[0];
  _rtXdot->high_pass_CSTATE_o[1] += Human_in_Loop_P.high_pass_A_a0[2] *
    Human_in_Loop_X.high_pass_CSTATE_o[1];
  _rtXdot->high_pass_CSTATE_o[1] += Human_in_Loop_P.high_pass_B_n3 *
    Human_in_Loop_B.Gain5;

  /* Derivatives for StateSpace: '<S45>/high_pass' */
  _rtXdot->high_pass_CSTATE_ff[0] += Human_in_Loop_P.high_pass_A_f[0] *
    Human_in_Loop_X.high_pass_CSTATE_ff[1];
  _rtXdot->high_pass_CSTATE_ff[1] += Human_in_Loop_P.high_pass_A_f[1] *
    Human_in_Loop_X.high_pass_CSTATE_ff[0];
  _rtXdot->high_pass_CSTATE_ff[1] += Human_in_Loop_P.high_pass_A_f[2] *
    Human_in_Loop_X.high_pass_CSTATE_ff[1];
  _rtXdot->high_pass_CSTATE_ff[1] += Human_in_Loop_P.high_pass_B_g *
    Human_in_Loop_B.Gain4_m;
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

  {
    /* Setup solver object */
    rtsiSetSimTimeStepPtr(&Human_in_Loop_M->solverInfo,
                          &Human_in_Loop_M->Timing.simTimeStep);
    rtsiSetTPtr(&Human_in_Loop_M->solverInfo, &rtmGetTPtr(Human_in_Loop_M));
    rtsiSetStepSizePtr(&Human_in_Loop_M->solverInfo,
                       &Human_in_Loop_M->Timing.stepSize0);
    rtsiSetdXPtr(&Human_in_Loop_M->solverInfo, &Human_in_Loop_M->derivs);
    rtsiSetContStatesPtr(&Human_in_Loop_M->solverInfo, (real_T **)
                         &Human_in_Loop_M->contStates);
    rtsiSetNumContStatesPtr(&Human_in_Loop_M->solverInfo,
      &Human_in_Loop_M->Sizes.numContStates);
    rtsiSetNumPeriodicContStatesPtr(&Human_in_Loop_M->solverInfo,
      &Human_in_Loop_M->Sizes.numPeriodicContStates);
    rtsiSetPeriodicContStateIndicesPtr(&Human_in_Loop_M->solverInfo,
      &Human_in_Loop_M->periodicContStateIndices);
    rtsiSetPeriodicContStateRangesPtr(&Human_in_Loop_M->solverInfo,
      &Human_in_Loop_M->periodicContStateRanges);
    rtsiSetErrorStatusPtr(&Human_in_Loop_M->solverInfo, (&rtmGetErrorStatus
      (Human_in_Loop_M)));
    rtsiSetRTModelPtr(&Human_in_Loop_M->solverInfo, Human_in_Loop_M);
  }

  rtsiSetSimTimeStep(&Human_in_Loop_M->solverInfo, MAJOR_TIME_STEP);
  Human_in_Loop_M->intgData.f[0] = Human_in_Loop_M->odeF[0];
  Human_in_Loop_M->contStates = ((X_Human_in_Loop_T *) &Human_in_Loop_X);
  rtsiSetSolverData(&Human_in_Loop_M->solverInfo, (void *)
                    &Human_in_Loop_M->intgData);
  rtsiSetSolverName(&Human_in_Loop_M->solverInfo,"ode1");
  rtmSetTPtr(Human_in_Loop_M, &Human_in_Loop_M->Timing.tArray[0]);
  Human_in_Loop_M->Timing.stepSize0 = 0.0002;

  /* block I/O */
  (void) memset(((void *) &Human_in_Loop_B), 0,
                sizeof(B_Human_in_Loop_T));

  /* states (continuous) */
  {
    (void) memset((void *)&Human_in_Loop_X, 0,
                  sizeof(X_Human_in_Loop_T));
  }

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

    /* Start for RateTransition: '<Root>/RT4' */
    Human_in_Loop_B.RT4[0] = Human_in_Loop_P.RT4_InitialCondition;
    Human_in_Loop_B.RT4[1] = Human_in_Loop_P.RT4_InitialCondition;

    /* Start for RateTransition: '<Root>/RT5' */
    Human_in_Loop_B.RT5[0] = Human_in_Loop_P.RT5_InitialCondition;
    Human_in_Loop_B.RT5[1] = Human_in_Loop_P.RT5_InitialCondition;

    /* Start for RateTransition: '<Root>/RT6' */
    Human_in_Loop_B.RT6[0] = Human_in_Loop_P.RT6_InitialCondition;
    Human_in_Loop_B.RT6[1] = Human_in_Loop_P.RT6_InitialCondition;
    Human_in_Loop_B.RT6[2] = Human_in_Loop_P.RT6_InitialCondition;

    /* Start for RateTransition: '<Root>/RT1' */
    Human_in_Loop_B.RT1[0] = Human_in_Loop_P.RT1_InitialCondition;
    Human_in_Loop_B.RT1[1] = Human_in_Loop_P.RT1_InitialCondition;
    Human_in_Loop_B.RT1[2] = Human_in_Loop_P.RT1_InitialCondition;
    Human_in_Loop_B.RT1[3] = Human_in_Loop_P.RT1_InitialCondition;

    /* Start for RateTransition: '<Root>/RT2' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_B.RT2[i] = Human_in_Loop_P.RT2_InitialCondition;
    }

    /* End of Start for RateTransition: '<Root>/RT2' */

    /* Start for RateTransition: '<Root>/RT3' */
    Human_in_Loop_B.RT3[0] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_B.RT3[1] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_B.RT3[2] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_B.RT3[3] = Human_in_Loop_P.RT3_InitialCondition;

    /* Start for Triggered SubSystem: '<S15>/Mean Calculate' */
    /* Start for DataStoreMemory: '<S18>/Data Store Memory' */
    memcpy(&Human_in_Loop_DW.SingleCycleData[0],
           &Human_in_Loop_P.DataStoreMemory_InitialValue[0], 2600U * sizeof
           (real_T));

    /* End of Start for SubSystem: '<S15>/Mean Calculate' */

    /* Start for DataStoreMemory: '<S16>/Data Store Memory' */
    memcpy(&Human_in_Loop_DW.EMG_Memory[0],
           &Human_in_Loop_P.DataStoreMemory_InitialValue_o[0], 24U * sizeof
           (real_T));

    /* Start for DataStoreMemory: '<Root>/Data Store Memory' */
    memcpy(&Human_in_Loop_DW.TorqueMem[0],
           &Human_in_Loop_P.DataStoreMemory_InitialValue_l[0], 4400U * sizeof
           (real_T));

    /* Start for DataStoreMemory: '<Root>/Data Store Memory1' */
    memcpy(&Human_in_Loop_DW.ParmReg[0],
           &Human_in_Loop_P.DataStoreMemory1_InitialValue[0], 10U * sizeof
           (real_T));
  }

  Human_in_Loop_PrevZCX.MeanCalculate_Trig_ZCE = UNINITIALIZED_ZCSIG;

  {
    int32_T i;

    /* InitializeConditions for DiscreteFilter: '<S32>/0.4low2' */
    Human_in_Loop_DW.u4low2_states[0] = Human_in_Loop_P.u4low2_InitialStates;
    Human_in_Loop_DW.u4low2_states[1] = Human_in_Loop_P.u4low2_InitialStates;
    Human_in_Loop_DW.u4low2_states[2] = Human_in_Loop_P.u4low2_InitialStates;

    /* InitializeConditions for UnitDelay: '<S76>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE =
      Human_in_Loop_P.UnitDelay1_InitialCondition_b;

    /* InitializeConditions for UnitDelay: '<S76>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE =
      Human_in_Loop_P.UnitDelay_InitialCondition_b;

    /* InitializeConditions for UnitDelay: '<S76>/Unit Delay2' */
    Human_in_Loop_DW.UnitDelay2_DSTATE =
      Human_in_Loop_P.UnitDelay2_InitialCondition;

    /* InitializeConditions for RateTransition: '<Root>/RT4' */
    Human_in_Loop_DW.RT4_Buffer0[0] = Human_in_Loop_P.RT4_InitialCondition;
    Human_in_Loop_DW.RT4_Buffer0[1] = Human_in_Loop_P.RT4_InitialCondition;
    Human_in_Loop_DW.RT4_write_buf = -1;
    Human_in_Loop_DW.RT4_read_buf = -1;

    /* InitializeConditions for RateTransition: '<Root>/RT5' */
    Human_in_Loop_DW.RT5_Buffer0[0] = Human_in_Loop_P.RT5_InitialCondition;
    Human_in_Loop_DW.RT5_Buffer0[1] = Human_in_Loop_P.RT5_InitialCondition;
    Human_in_Loop_DW.RT5_write_buf = -1;
    Human_in_Loop_DW.RT5_read_buf = -1;

    /* InitializeConditions for UnitDelay: '<S52>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE_f =
      Human_in_Loop_P.UnitDelay1_InitialCondition_bo;

    /* InitializeConditions for UnitDelay: '<S52>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_a =
      Human_in_Loop_P.UnitDelay_InitialCondition_c;

    /* InitializeConditions for UnitDelay: '<S52>/Unit Delay2' */
    Human_in_Loop_DW.UnitDelay2_DSTATE_e =
      Human_in_Loop_P.UnitDelay2_InitialCondition_c;

    /* InitializeConditions for RateTransition: '<Root>/RT6' */
    Human_in_Loop_DW.RT6_Buffer0[0] = Human_in_Loop_P.RT6_InitialCondition;
    Human_in_Loop_DW.RT6_Buffer0[1] = Human_in_Loop_P.RT6_InitialCondition;
    Human_in_Loop_DW.RT6_Buffer0[2] = Human_in_Loop_P.RT6_InitialCondition;
    Human_in_Loop_DW.RT6_write_buf = -1;
    Human_in_Loop_DW.RT6_read_buf = -1;

    /* InitializeConditions for RateTransition: '<Root>/RT1' */
    Human_in_Loop_DW.RT1_Buffer0[0] = Human_in_Loop_P.RT1_InitialCondition;
    Human_in_Loop_DW.RT1_Buffer0[1] = Human_in_Loop_P.RT1_InitialCondition;
    Human_in_Loop_DW.RT1_Buffer0[2] = Human_in_Loop_P.RT1_InitialCondition;
    Human_in_Loop_DW.RT1_Buffer0[3] = Human_in_Loop_P.RT1_InitialCondition;
    Human_in_Loop_DW.RT1_write_buf = -1;
    Human_in_Loop_DW.RT1_read_buf = -1;

    /* InitializeConditions for RateTransition: '<Root>/RT2' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_DW.RT2_Buffer0[i] = Human_in_Loop_P.RT2_InitialCondition;
    }

    Human_in_Loop_DW.RT2_write_buf = -1;
    Human_in_Loop_DW.RT2_read_buf = -1;

    /* End of InitializeConditions for RateTransition: '<Root>/RT2' */

    /* InitializeConditions for RateTransition: '<Root>/RT3' */
    Human_in_Loop_DW.RT3_Buffer0[0] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_DW.RT3_Buffer0[1] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_DW.RT3_Buffer0[2] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_DW.RT3_Buffer0[3] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_DW.RT3_write_buf = -1;
    Human_in_Loop_DW.RT3_read_buf = -1;

    /* InitializeConditions for Delay: '<S16>/Delay3' */
    Human_in_Loop_DW.Delay3_DSTATE = Human_in_Loop_P.Delay3_InitialCondition;

    /* InitializeConditions for Delay: '<S16>/Delay2' */
    Human_in_Loop_DW.Delay2_DSTATE = Human_in_Loop_P.Delay2_InitialCondition;

    /* InitializeConditions for Delay: '<S16>/Delay7' */
    Human_in_Loop_DW.Delay7_DSTATE = Human_in_Loop_P.Delay7_InitialCondition;

    /* InitializeConditions for Delay: '<S16>/Delay6' */
    Human_in_Loop_DW.Delay6_DSTATE = Human_in_Loop_P.Delay6_InitialCondition;

    /* InitializeConditions for StateSpace: '<S40>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE[0] =
      Human_in_Loop_P.low_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S41>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_b[0] =
      Human_in_Loop_P.low_pass_InitialCondition_i;

    /* InitializeConditions for StateSpace: '<S42>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_m[0] =
      Human_in_Loop_P.low_pass_InitialCondition_n;

    /* InitializeConditions for StateSpace: '<S43>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_d[0] =
      Human_in_Loop_P.low_pass_InitialCondition_o;

    /* InitializeConditions for StateSpace: '<S45>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_j[0] =
      Human_in_Loop_P.low_pass_InitialCondition_ir;

    /* InitializeConditions for StateSpace: '<S44>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_a[0] =
      Human_in_Loop_P.low_pass_InitialCondition_h;

    /* InitializeConditions for StateSpace: '<S40>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE[1] =
      Human_in_Loop_P.low_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S41>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_b[1] =
      Human_in_Loop_P.low_pass_InitialCondition_i;

    /* InitializeConditions for StateSpace: '<S42>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_m[1] =
      Human_in_Loop_P.low_pass_InitialCondition_n;

    /* InitializeConditions for StateSpace: '<S43>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_d[1] =
      Human_in_Loop_P.low_pass_InitialCondition_o;

    /* InitializeConditions for StateSpace: '<S45>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_j[1] =
      Human_in_Loop_P.low_pass_InitialCondition_ir;

    /* InitializeConditions for StateSpace: '<S44>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_a[1] =
      Human_in_Loop_P.low_pass_InitialCondition_h;
    for (i = 0; i < 6; i++) {
      /* InitializeConditions for Delay: '<S16>/Delay1' */
      Human_in_Loop_DW.Delay1_DSTATE[i] =
        Human_in_Loop_P.Delay1_InitialCondition;

      /* InitializeConditions for Delay: '<S16>/Delay5' */
      Human_in_Loop_DW.Delay5_DSTATE[i] =
        Human_in_Loop_P.Delay5_InitialCondition;

      /* InitializeConditions for Delay: '<S16>/Delay4' */
      Human_in_Loop_DW.Delay4_DSTATE[i] =
        Human_in_Loop_P.Delay4_InitialCondition;

      /* InitializeConditions for Delay: '<S16>/Delay8' */
      Human_in_Loop_DW.Delay8_DSTATE[i] =
        Human_in_Loop_P.Delay8_InitialCondition;
    }

    /* InitializeConditions for UnitDelay: '<S55>/UD' */
    Human_in_Loop_DW.UD_DSTATE = Human_in_Loop_P.DiscreteDerivative_ICPrevScaled;

    /* InitializeConditions for DiscreteFilter: '<S29>/0.4low1' */
    Human_in_Loop_DW.u4low1_states[0] = Human_in_Loop_P.u4low1_InitialStates;
    Human_in_Loop_DW.u4low1_states[1] = Human_in_Loop_P.u4low1_InitialStates;
    Human_in_Loop_DW.u4low1_states[2] = Human_in_Loop_P.u4low1_InitialStates;

    /* InitializeConditions for UnitDelay: '<S51>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_d =
      Human_in_Loop_P.UnitDelay_InitialCondition_i;

    /* InitializeConditions for UnitDelay: '<S51>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE_f2 =
      Human_in_Loop_P.UnitDelay1_InitialCondition_e;

    /* InitializeConditions for StateSpace: '<S40>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE[0] =
      Human_in_Loop_P.high_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S41>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_p[0] =
      Human_in_Loop_P.high_pass_InitialCondition_a;

    /* InitializeConditions for StateSpace: '<S42>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_f[0] =
      Human_in_Loop_P.high_pass_InitialCondition_d;

    /* InitializeConditions for StateSpace: '<S43>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_a[0] =
      Human_in_Loop_P.high_pass_InitialCondition_m;

    /* InitializeConditions for StateSpace: '<S44>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_o[0] =
      Human_in_Loop_P.high_pass_InitialCondition_h;

    /* InitializeConditions for StateSpace: '<S45>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_ff[0] =
      Human_in_Loop_P.high_pass_InitialCondition_e;

    /* InitializeConditions for StateSpace: '<S40>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE[1] =
      Human_in_Loop_P.high_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S41>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_p[1] =
      Human_in_Loop_P.high_pass_InitialCondition_a;

    /* InitializeConditions for StateSpace: '<S42>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_f[1] =
      Human_in_Loop_P.high_pass_InitialCondition_d;

    /* InitializeConditions for StateSpace: '<S43>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_a[1] =
      Human_in_Loop_P.high_pass_InitialCondition_m;

    /* InitializeConditions for StateSpace: '<S44>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_o[1] =
      Human_in_Loop_P.high_pass_InitialCondition_h;

    /* InitializeConditions for StateSpace: '<S45>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_ff[1] =
      Human_in_Loop_P.high_pass_InitialCondition_e;

    /* SystemInitialize for MATLAB Function: '<S30>/FootSwitch Filter' */
    Human_in_Loop_DW.foot_state = 0.0;
    Human_in_Loop_DW.filter_time = 0.0;

    /* SystemInitialize for MATLAB Function: '<S7>/State Machine' */
    Human_in_Loop_DW.reg_stride_time = 1.0;
    Human_in_Loop_DW.reg_stride_time_count = 0.0;
    Human_in_Loop_DW.reg_mode = 1.0;
    Human_in_Loop_DW.reg_state = 1.0;
    Human_in_Loop_DW.bt_run = 0.0;
    Human_in_Loop_DW.reg_last_switch = 1.0;

    /* SystemInitialize for S-Function (rti_commonblock): '<S8>/S-Function1' incorporates:
     *  SubSystem: '<Root>/Control Module'
     */
    Human_in_Loo_ControlModule_Init();

    /* End of SystemInitialize for S-Function (rti_commonblock): '<S8>/S-Function1' */

    /* SystemInitialize for Triggered SubSystem: '<S15>/Mean Calculate' */
    /* InitializeConditions for UnitDelay: '<S18>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_g =
      Human_in_Loop_P.UnitDelay_InitialCondition;
    for (i = 0; i < 26; i++) {
      /* InitializeConditions for UnitDelay: '<S18>/Unit Delay1' */
      Human_in_Loop_DW.UnitDelay1_DSTATE_e[i] =
        Human_in_Loop_P.UnitDelay1_InitialCondition;

      /* SystemInitialize for Outport: '<S18>/Mean' */
      Human_in_Loop_B.Mean_c[i] = Human_in_Loop_P.Mean_Y0;
    }

    /* SystemInitialize for Outport: '<S18>/Count' */
    Human_in_Loop_B.count = Human_in_Loop_P.Count_Y0;

    /* End of SystemInitialize for SubSystem: '<S15>/Mean Calculate' */

    /* SystemInitialize for MATLAB Function: '<S3>/torque_track_loss' */
    Human_in_Loop_DW.last_footstate_p = 0.0;
    Human_in_Loop_DW.loss_reg = 0.0;
    memset(&Human_in_Loop_DW.loss_mem[0], 0, 10U * sizeof(real_T));

    /* SystemInitialize for MATLAB Function: '<S32>/MATLAB Function' */
    for (i = 0; i < 15; i++) {
      Human_in_Loop_DW.data[i] = 1.0;
    }

    /* End of SystemInitialize for MATLAB Function: '<S32>/MATLAB Function' */
  }
}

/* Model terminate function */
void Human_in_Loop_terminate(void)
{
  /* Terminate for S-Function (rti_commonblock): '<S56>/S-Function1' */

  /* --- Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL1 --- */
  /* --- [RTIEMC, Encoder] - DIO class: 2 - Unit: 5 - Port: 1 - Channel: 10 --- */
  {
    /* Deactivates encoder interface functionality */
    DioCl2EncoderIn_stop(pRTIEmcEncoder_Unit_5_DioCl_2_Port_1_Ch10);
  }

  /* Terminate for S-Function (rti_commonblock): '<S57>/S-Function1' */

  /* --- Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL3 --- */
  /* --- [RTIEMC, Encoder] - DIO class: 2 - Unit: 3 - Port: 1 - Channel: 5 --- */
  {
    /* Deactivates encoder interface functionality */
    DioCl2EncoderIn_stop(pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5);
  }

  /* Terminate for S-Function (rti_commonblock): '<S8>/S-Function1' incorporates:
   *  SubSystem: '<Root>/Control Module'
   */
  Human_in_Loo_ControlModule_Term();

  /* End of Terminate for S-Function (rti_commonblock): '<S8>/S-Function1' */

  /* Terminate for S-Function (rti_commonblock): '<S58>/S-Function1' incorporates:
   *  Constant: '<S29>/VCC1'
   */

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup1 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 15 --- */

  /* disable digital output channel 15-15 on port 3 *
   * (set to high-impedance), when the simulation terminates       */
  DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_3_Ch_15,
    DIO_CLASS1_HIGH_Z_ON);
  DioCl1DigOut_write(pRTIDioC1DigOut_Port_3_Ch_15);

  /* Terminate for S-Function (rti_commonblock): '<S59>/S-Function1' incorporates:
   *  Constant: '<S29>/VCC3'
   */

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup3 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 13 --- */

  /* disable digital output channel 13-13 on port 3 *
   * (set to high-impedance), when the simulation terminates       */
  DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_3_Ch_13,
    DIO_CLASS1_HIGH_Z_ON);
  DioCl1DigOut_write(pRTIDioC1DigOut_Port_3_Ch_13);

  /* Terminate for S-Function (rti_commonblock): '<S64>/S-Function1' incorporates:
   *  Constant: '<S30>/Constant'
   */

  /* --- Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_OUT_BL1 --- */
  /* --- [RTI120X, BITOUT] - Port: 1 - Channel: 1 --- */

  /* disable digital output channel 1-1 on port 1 *
   * (set to high-impedance), when the simulation terminates       */
  DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_1_Ch_1,
    DIO_CLASS1_HIGH_Z_ON);
  DioCl1DigOut_write(pRTIDioC1DigOut_Port_1_Ch_1);

  /* Terminate for S-Function (rti_commonblock): '<S72>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:3 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X3])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S66>/S-Function1' */

  /* dSPACE RTICAN STD Srvc-Message Block */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (CANTP1_RX_SPMSG_M1_C1_STD)) == DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S67>/S-Function1' */

  /* dSPACE RTICAN STD Srvc-Message Block */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (CANTP1_RX_SPMSG_M1_C2_STD)) == DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S70>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:100 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][1] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S71>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:100 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][1] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }
}
