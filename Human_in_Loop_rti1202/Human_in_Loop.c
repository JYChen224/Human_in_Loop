/*
 * Human_in_Loop.c
 *
 * Code generation for model "Human_in_Loop".
 *
 * Model version              : 1.1163
 * Simulink Coder version : 8.13 (R2017b) 24-Jul-2017
 * C source code generated on : Wed May 22 18:33:58 2019
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
static real_T Human_in_Loop_xnrm2(int32_T n, const real_T x_data[], int32_T ix0);
static void Human_in_Loop_xgeqp3(real_T A_data[], int32_T A_size[2], real_T
  tau_data[], int32_T *tau_size, int32_T jpvt[2]);
static void Human_in_Loop_mldivide(const real_T A_data[], const int32_T A_size[2],
  const real_T B_data[], const int32_T *B_size, real_T Y[2]);

/* Forward declaration for local functions */
static void Human_in_Loop_mldivide_g(const real_T A[16], real_T B[4]);
static void Human_in_Loop_power(const real_T a_data[], const int32_T a_size[2],
  real_T y_data[], int32_T y_size[2]);
static void Human_in_Loop_power_a(const real_T a_data[], const int32_T a_size[2],
  real_T y_data[], int32_T y_size[2]);

/* Forward declaration for local functions */
static real_T Human_in_Loop_norm(const real_T x_data[], const int32_T x_size[2]);
static real_T Human_in_Loop_mean(const real_T x[10]);
static real_T Human_in_Loop_xnrm2_o(const real_T x[30], int32_T ix0);
static real_T Human_in_Loop_xnrm2_ov(int32_T n, const real_T x[30], int32_T ix0);
static void Human_in_Loop_xgeqp3_n(real_T A[30], real_T tau[2], int32_T jpvt[2]);
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
  int_T nXc = 48;
  rtsiSetSimTimeStep(si,MINOR_TIME_STEP);
  rtsiSetdX(si, f0);
  Human_in_Loop_derivatives();
  rtsiSetT(si, tnew);
  for (i = 0; i < nXc; ++i) {
    x[i] += h * f0[i];
  }

  rtsiSetSimTimeStep(si,MAJOR_TIME_STEP);
}

/* Function for MATLAB Function: '<S38>/Estimation' */
static real_T Human_in_Loop_xnrm2(int32_T n, const real_T x_data[], int32_T ix0)
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
      y = fabs(x_data[ix0 - 1]);
    } else {
      scale = 3.3121686421112381E-170;
      kend = (ix0 + n) - 1;
      for (k = ix0; k <= kend; k++) {
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

/* Function for MATLAB Function: '<S38>/Estimation' */
static void Human_in_Loop_xgeqp3(real_T A_data[], int32_T A_size[2], real_T
  tau_data[], int32_T *tau_size, int32_T jpvt[2])
{
  int32_T m;
  real_T work[2];
  real_T vn1[2];
  int32_T k;
  int32_T ix;
  real_T smax;
  int32_T iy;
  real_T xnorm;
  int32_T c_ix;
  real_T b_c;
  int32_T f;
  int32_T b_ia;
  int32_T jy;
  int32_T exitg1;
  boolean_T exitg2;
  m = A_size[0];
  *tau_size = 2;
  jpvt[0] = 1;
  work[0] = 0.0;
  smax = Human_in_Loop_xnrm2(m, A_data, 1);
  k = 1 + m;
  vn1[0] = smax;
  jpvt[1] = 2;
  work[1] = 0.0;
  smax = Human_in_Loop_xnrm2(m, A_data, k);
  vn1[1] = smax;
  k = 0;
  ix = 0;
  smax = fabs(vn1[0]);
  iy = 2;
  while (iy <= 2) {
    ix++;
    b_c = fabs(vn1[ix]);
    if (b_c > smax) {
      k = 1;
      smax = b_c;
    }

    iy = 3;
  }

  if (k + 1 != 1) {
    ix = m * k;
    iy = 0;
    for (jy = 1; jy <= m; jy++) {
      smax = A_data[ix];
      A_data[ix] = A_data[iy];
      A_data[iy] = smax;
      ix++;
      iy++;
    }

    ix = jpvt[k];
    jpvt[k] = 1;
    jpvt[0] = ix;
  }

  smax = A_data[0];
  b_c = 0.0;
  if (!(m <= 0)) {
    xnorm = Human_in_Loop_xnrm2(m - 1, A_data, 2);
    if (xnorm != 0.0) {
      xnorm = rt_hypotd_snf(A_data[0], xnorm);
      if (A_data[0] >= 0.0) {
        xnorm = -xnorm;
      }

      if (fabs(xnorm) < 1.0020841800044864E-292) {
        k = 0;
        ix = m;
        do {
          k++;
          for (iy = 1; iy + 1 <= ix; iy++) {
            A_data[iy] *= 9.9792015476736E+291;
          }

          xnorm *= 9.9792015476736E+291;
          smax *= 9.9792015476736E+291;
        } while (!(fabs(xnorm) >= 1.0020841800044864E-292));

        xnorm = rt_hypotd_snf(smax, Human_in_Loop_xnrm2(m - 1, A_data, 2));
        if (smax >= 0.0) {
          xnorm = -xnorm;
        }

        b_c = (xnorm - smax) / xnorm;
        smax = 1.0 / (smax - xnorm);
        for (iy = 1; iy + 1 <= ix; iy++) {
          A_data[iy] *= smax;
        }

        for (ix = 1; ix <= k; ix++) {
          xnorm *= 1.0020841800044864E-292;
        }

        smax = xnorm;
      } else {
        b_c = (xnorm - A_data[0]) / xnorm;
        smax = 1.0 / (A_data[0] - xnorm);
        k = m;
        for (ix = 1; ix + 1 <= k; ix++) {
          A_data[ix] *= smax;
        }

        smax = xnorm;
      }
    }
  }

  tau_data[0] = b_c;
  A_data[0] = smax;
  smax = A_data[0];
  A_data[0] = 1.0;
  if (tau_data[0] != 0.0) {
    k = m;
    ix = m - 1;
    while ((k > 0) && (A_data[ix] == 0.0)) {
      k--;
      ix--;
    }

    ix = 1;
    exitg2 = false;
    while ((!exitg2) && (ix > 0)) {
      iy = m;
      do {
        exitg1 = 0;
        if (iy + 1 <= m + k) {
          if (A_data[iy] != 0.0) {
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
    k = 0;
    ix = 0;
  }

  if (k > 0) {
    if (ix != 0) {
      work[0] = 0.0;
      iy = 0;
      for (jy = 1 + m; jy <= m + 1; jy += m) {
        c_ix = 0;
        b_c = 0.0;
        f = (jy + k) - 1;
        for (b_ia = jy; b_ia <= f; b_ia++) {
          b_c += A_data[b_ia - 1] * A_data[c_ix];
          c_ix++;
        }

        work[iy] += b_c;
        iy++;
      }
    }

    if (!(-tau_data[0] == 0.0)) {
      iy = m;
      jy = 0;
      c_ix = 1;
      while (c_ix <= ix) {
        if (work[jy] != 0.0) {
          b_c = work[jy] * -tau_data[0];
          c_ix = 0;
          f = k + iy;
          for (b_ia = iy; b_ia + 1 <= f; b_ia++) {
            A_data[b_ia] += A_data[c_ix] * b_c;
            c_ix++;
          }
        }

        jy++;
        iy += m;
        c_ix = 2;
      }
    }
  }

  A_data[0] = smax;
  jy = 1 + m;
  m--;
  smax = A_data[jy];
  b_c = 0.0;
  if (!(m <= 0)) {
    xnorm = Human_in_Loop_xnrm2(m - 1, A_data, jy + 2);
    if (xnorm != 0.0) {
      xnorm = rt_hypotd_snf(A_data[jy], xnorm);
      if (A_data[jy] >= 0.0) {
        xnorm = -xnorm;
      }

      if (fabs(xnorm) < 1.0020841800044864E-292) {
        k = 0;
        ix = jy + m;
        do {
          k++;
          for (iy = jy + 1; iy + 1 <= ix; iy++) {
            A_data[iy] *= 9.9792015476736E+291;
          }

          xnorm *= 9.9792015476736E+291;
          smax *= 9.9792015476736E+291;
        } while (!(fabs(xnorm) >= 1.0020841800044864E-292));

        xnorm = rt_hypotd_snf(smax, Human_in_Loop_xnrm2(m - 1, A_data, jy + 2));
        if (smax >= 0.0) {
          xnorm = -xnorm;
        }

        b_c = (xnorm - smax) / xnorm;
        smax = 1.0 / (smax - xnorm);
        ix = jy + m;
        for (iy = jy + 1; iy + 1 <= ix; iy++) {
          A_data[iy] *= smax;
        }

        for (ix = 1; ix <= k; ix++) {
          xnorm *= 1.0020841800044864E-292;
        }

        smax = xnorm;
      } else {
        b_c = (xnorm - A_data[jy]) / xnorm;
        smax = 1.0 / (A_data[jy] - xnorm);
        k = jy + m;
        for (ix = jy + 1; ix + 1 <= k; ix++) {
          A_data[ix] *= smax;
        }

        smax = xnorm;
      }
    }
  }

  tau_data[1] = b_c;
  A_data[jy] = smax;
}

/* Function for MATLAB Function: '<S38>/Estimation' */
static void Human_in_Loop_mldivide(const real_T A_data[], const int32_T A_size[2],
  const real_T B_data[], const int32_T *B_size, real_T Y[2])
{
  real_T b_A_data[40];
  real_T tau_data[2];
  int32_T jpvt[2];
  int32_T rankR;
  real_T tol;
  real_T b_B_data[20];
  int32_T m;
  int32_T c_i;
  int32_T b_A_size[2];
  real_T x_data_idx_1;
  real_T x_data_idx_0;
  real_T x_data_idx_2;
  real_T x_data_idx_3;
  if (A_size[0] == 2) {
    x_data_idx_0 = A_data[0];
    x_data_idx_1 = A_data[1];
    x_data_idx_2 = A_data[2];
    x_data_idx_3 = A_data[3];
    jpvt[0] = 1;
    rankR = 0;
    if (fabs(A_data[1]) > fabs(A_data[0])) {
      rankR = 1;
    }

    if (A_data[rankR] != 0.0) {
      if (rankR != 0) {
        jpvt[0] = 2;
        x_data_idx_0 = A_data[0];
        x_data_idx_1 = A_data[1];
        x_data_idx_2 = A_data[2];
        x_data_idx_3 = A_data[3];
        tol = x_data_idx_0;
        x_data_idx_0 = x_data_idx_1;
        x_data_idx_1 = tol;
        tol = x_data_idx_2;
        x_data_idx_2 = x_data_idx_3;
        x_data_idx_3 = tol;
      }

      x_data_idx_1 /= x_data_idx_0;
    }

    if (x_data_idx_2 != 0.0) {
      x_data_idx_3 += x_data_idx_1 * -x_data_idx_2;
    }

    Y[0] = B_data[0];
    Y[1] = B_data[1];
    if (jpvt[0] != 1) {
      Y[0] = B_data[1];
      Y[1] = B_data[0];
    }

    if (Y[0] != 0.0) {
      rankR = 2;
      while (rankR < 3) {
        Y[1] -= Y[0] * x_data_idx_1;
        rankR = 3;
      }
    }

    if (Y[1] != 0.0) {
      Y[1] /= x_data_idx_3;
      rankR = 1;
      while (rankR <= 1) {
        Y[0] -= Y[1] * x_data_idx_2;
        rankR = 2;
      }
    }

    if (Y[0] != 0.0) {
      Y[0] /= x_data_idx_0;
    }
  } else {
    b_A_size[0] = A_size[0];
    b_A_size[1] = 2;
    m = A_size[0] * A_size[1];
    if (0 <= m - 1) {
      memcpy(&b_A_data[0], &A_data[0], m * sizeof(real_T));
    }

    Human_in_Loop_xgeqp3(b_A_data, b_A_size, tau_data, &m, jpvt);
    rankR = 0;
    tol = (real_T)b_A_size[0] * fabs(b_A_data[0]) * 2.2204460492503131E-16;
    while ((rankR < 2) && (!(fabs(b_A_data[b_A_size[0] * rankR + rankR]) <= tol)))
    {
      rankR++;
    }

    Y[0] = 0.0;
    Y[1] = 0.0;
    m = *B_size;
    if (0 <= m - 1) {
      memcpy(&b_B_data[0], &B_data[0], m * sizeof(real_T));
    }

    m = b_A_size[0];
    if (tau_data[0] != 0.0) {
      tol = b_B_data[0];
      for (c_i = 1; c_i + 1 <= m; c_i++) {
        tol += b_A_data[c_i] * b_B_data[c_i];
      }

      tol *= tau_data[0];
      if (tol != 0.0) {
        b_B_data[0] -= tol;
        for (c_i = 1; c_i + 1 <= m; c_i++) {
          b_B_data[c_i] -= b_A_data[c_i] * tol;
        }
      }
    }

    if (tau_data[1] != 0.0) {
      tol = b_B_data[1];
      for (c_i = 2; c_i + 1 <= m; c_i++) {
        tol += b_A_data[c_i + b_A_size[0]] * b_B_data[c_i];
      }

      tol *= tau_data[1];
      if (tol != 0.0) {
        b_B_data[1] -= tol;
        for (c_i = 2; c_i + 1 <= m; c_i++) {
          b_B_data[c_i] -= b_A_data[c_i + b_A_size[0]] * tol;
        }
      }
    }

    for (m = 0; m + 1 <= rankR; m++) {
      Y[jpvt[m] - 1] = b_B_data[m];
    }

    for (rankR--; rankR + 1 > 0; rankR--) {
      Y[jpvt[rankR] - 1] /= b_A_data[b_A_size[0] * rankR + rankR];
      m = 1;
      while (m <= rankR) {
        Y[jpvt[0] - 1] -= Y[jpvt[rankR] - 1] * b_A_data[b_A_size[0] * rankR];
        m = 2;
      }
    }
  }
}

/* System initialize for function-call system: '<S34>/Serial Decoding System' */
void Human_SerialDecodingSystem_Init(void)
{
  int32_T i;

  /* SystemInitialize for MATLAB Function: '<S38>/Estimation' */
  Human_in_Loop_DW.last_count = -1.0;
  for (i = 0; i < 20; i++) {
    Human_in_Loop_DW.x[i] = 10.0 * (real_T)i + 10.0;
    Human_in_Loop_DW.y[i] = 0.0;
  }

  Human_in_Loop_DW.theta[0] = 0.0;
  Human_in_Loop_DW.theta[1] = 0.0;
  Human_in_Loop_DW.num = 1.0;

  /* End of SystemInitialize for MATLAB Function: '<S38>/Estimation' */

  /* SystemInitialize for Outport: '<S38>/E' */
  Human_in_Loop_B.E = Human_in_Loop_P.E_Y0;
}

/* System reset for function-call system: '<S34>/Serial Decoding System' */
void Huma_SerialDecodingSystem_Reset(void)
{
  int32_T i;

  /* SystemReset for MATLAB Function: '<S38>/Estimation' */
  Human_in_Loop_DW.last_count = -1.0;
  for (i = 0; i < 20; i++) {
    Human_in_Loop_DW.x[i] = 10.0 * (real_T)i + 10.0;
    Human_in_Loop_DW.y[i] = 0.0;
  }

  Human_in_Loop_DW.theta[0] = 0.0;
  Human_in_Loop_DW.theta[1] = 0.0;
  Human_in_Loop_DW.num = 1.0;

  /* End of SystemReset for MATLAB Function: '<S38>/Estimation' */
}

/* Output and update for function-call system: '<S34>/Serial Decoding System' */
void Human_in_L_SerialDecodingSystem(void)
{
  real_T Data1;
  real_T Data2;
  real_T i;
  real_T y_data[20];
  real_T x_data[20];
  real_T b_x[12];
  real_T a[24];
  real_T tmp_data[40];
  int32_T i_0;
  int32_T i_1;
  int32_T loop_ub;
  int32_T tmp_size[2];

  /* S-Function (rti_commonblock): '<S39>/S-Function1' */
  /* This comment workarounds a code generation problem */

  /* dSPACE I/O Board DS1202SER #1 Unit:GENSER Group:RECEIVE */
  {
    Human_in_Loop_B.SFunction1_o3_l = dsser_receive(rtiDS1202SER_B1_Ser[0], 16U,
      (UInt8 *) &Human_in_Loop_B.SFunction1_o1_l[0], (UInt32 *)
      &Human_in_Loop_B.SFunction1_o2_b);
  }

  /* MATLAB Function: '<S38>/MATLAB Function1' */
  /* MATLAB Function 'Sensor Data/Cosmed/COSMED/Serial Decoding System/MATLAB Function1': '<S41>:1' */
  /* '<S41>:1:3' */
  Data1 = 0.0;

  /* '<S41>:1:4' */
  Data2 = 0.0;

  /* '<S41>:1:6' */
  for (i = 1.0; Human_in_Loop_B.SFunction1_o1_l[(int32_T)i - 1] != 44; i++) {
    /* '<S41>:1:8' */
    /* '<S41>:1:9' */
    Data1 = ((real_T)Human_in_Loop_B.SFunction1_o1_l[(int32_T)i - 1] - 48.0) +
      Data1 * 10.0;

    /* '<S41>:1:10' */
  }

  /* '<S41>:1:13' */
  for (i++; i < Human_in_Loop_B.SFunction1_o2_b; i++) {
    /* '<S41>:1:15' */
    /* '<S41>:1:16' */
    Data2 = ((real_T)Human_in_Loop_B.SFunction1_o1_l[(int32_T)i - 1] - 48.0) +
      Data2 * 10.0;

    /* '<S41>:1:17' */
  }

  Human_in_Loop_B.Data1 = Data1;
  Human_in_Loop_B.Data2 = Data2;

  /* End of MATLAB Function: '<S38>/MATLAB Function1' */

  /* MATLAB Function: '<S38>/Estimation' */
  /* MATLAB Function 'Sensor Data/Cosmed/COSMED/Serial Decoding System/Estimation': '<S40>:1' */
  /* '<S40>:1:49' */
  if ((Human_in_Loop_DW.last_count != Human_in_Loop_B.RT2_g) ||
      (Human_in_Loop_B.RT2_g == 0.0)) {
    /* '<S40>:1:21' */
    /* '<S40>:1:22' */
    Human_in_Loop_DW.last_count = Human_in_Loop_B.RT2_g;

    /* '<S40>:1:23' */
    /* '<S40>:1:24' */
    for (i_0 = 0; i_0 < 20; i_0++) {
      Human_in_Loop_DW.x[i_0] = 10.0 * (real_T)i_0 + 10.0;
      Human_in_Loop_DW.y[i_0] = 0.0;
    }

    /* '<S40>:1:25' */
    Human_in_Loop_DW.theta[0] = 0.0;
    Human_in_Loop_DW.theta[1] = 0.0;

    /* '<S40>:1:26' */
    Human_in_Loop_DW.num = 1.0;
  } else {
    /* '<S40>:1:28' */
    Human_in_Loop_DW.num++;
  }

  if (Human_in_Loop_B.RT2_g != 0.0) {
    /* '<S40>:1:32' */
    /* '<S40>:1:33' */
    Human_in_Loop_DW.x[(int32_T)Human_in_Loop_DW.num - 1] =
      Human_in_Loop_B.RT1_c;

    /* '<S40>:1:34' */
    Human_in_Loop_DW.y[(int32_T)Human_in_Loop_DW.num - 1] = 0.278 *
      Human_in_Loop_B.Data1 + 0.075 * Human_in_Loop_B.Data2;
    if (Human_in_Loop_DW.num >= 2.0) {
      /* '<S40>:1:37' */
      /* '<S40>:1:38' */
      i_0 = (int32_T)Human_in_Loop_DW.num;
      loop_ub = (int32_T)Human_in_Loop_DW.num;
      for (i_1 = 0; i_1 < loop_ub; i_1++) {
        y_data[i_1] = -Human_in_Loop_DW.x[i_1] / 42.0;
      }

      if (0 <= i_0 - 1) {
        memcpy(&x_data[0], &y_data[0], i_0 * sizeof(real_T));
      }

      for (loop_ub = 0; loop_ub + 1 <= i_0; loop_ub++) {
        x_data[loop_ub] = exp(x_data[loop_ub]);
      }

      /* '<S40>:1:39' */
      /* '<S40>:1:40' */
      tmp_size[0] = i_0;
      tmp_size[1] = 2;
      for (i_1 = 0; i_1 < i_0; i_1++) {
        tmp_data[i_1] = 1.0 - x_data[i_1];
      }

      loop_ub = (int32_T)Human_in_Loop_DW.num;
      for (i_1 = 0; i_1 < loop_ub; i_1++) {
        tmp_data[i_1 + i_0] = 1.0;
      }

      i_0 = (int32_T)Human_in_Loop_DW.num;
      loop_ub = (int32_T)Human_in_Loop_DW.num;
      if (0 <= loop_ub - 1) {
        memcpy(&y_data[0], &Human_in_Loop_DW.y[0], loop_ub * sizeof(real_T));
      }

      Human_in_Loop_mldivide(tmp_data, tmp_size, y_data, &i_0,
        Human_in_Loop_DW.theta);
    }
  }

  /* '<S40>:1:45' */
  /* '<S40>:1:47' */
  memcpy(&Human_in_Loop_B.Time_p[0], &Human_in_Loop_DW.x[0], 12U * sizeof(real_T));

  /* '<S40>:1:48' */
  memcpy(&Human_in_Loop_B.Y[0], &Human_in_Loop_DW.y[0], 12U * sizeof(real_T));

  /* '<S40>:1:49' */
  for (i_0 = 0; i_0 < 12; i_0++) {
    b_x[i_0] = -Human_in_Loop_B.Time_p[i_0] / 42.0;
  }

  for (loop_ub = 0; loop_ub < 12; loop_ub++) {
    Data1 = b_x[loop_ub];
    Data1 = exp(Data1);
    a[loop_ub] = 1.0 - Data1;
    a[loop_ub + 12] = 1.0;
  }

  /* '<S40>:1:50' */
  Human_in_Loop_B.E = Human_in_Loop_DW.theta[0] + Human_in_Loop_DW.theta[1];
  for (i_1 = 0; i_1 < 12; i_1++) {
    Human_in_Loop_B.Y_E[i_1] = 0.0;
    Human_in_Loop_B.Y_E[i_1] += a[i_1] * Human_in_Loop_DW.theta[0];
    Human_in_Loop_B.Y_E[i_1] += a[i_1 + 12] * Human_in_Loop_DW.theta[1];
  }

  for (i_1 = 0; i_1 < 12; i_1++) {
    Data1 = a[i_1] * Human_in_Loop_DW.theta[0];
    Data1 += a[i_1 + 12] * Human_in_Loop_DW.theta[1];
    b_x[i_1] = Data1;
  }

  for (i_1 = 0; i_1 < 12; i_1++) {
    Human_in_Loop_B.est_curve[i_1 << 1] = Human_in_Loop_B.Y[i_1];
  }

  for (i_1 = 0; i_1 < 12; i_1++) {
    Human_in_Loop_B.est_curve[1 + (i_1 << 1)] = b_x[i_1];
  }

  /* End of MATLAB Function: '<S38>/Estimation' */
}

/*
 * Output and update for atomic system:
 *    '<S29>/Mux1'
 *    '<S29>/Mux3'
 */
void Human_in_Loop_Mux1(real_T rtu_x1, real_T rtu_x2, real_T rtu_x3, real_T
  rtu_x4, real_T rtu_x5, real_T rtu_x6, B_Mux1_Human_in_Loop_T *localB)
{
  /* MATLAB Function 'Sensor Data/EMG module/Mux1': '<S48>:1' */
  /* '<S48>:1:3' */
  localB->x[0] = rtu_x1;
  localB->x[1] = rtu_x2;
  localB->x[2] = rtu_x3;
  localB->x[3] = rtu_x4;
  localB->x[4] = rtu_x5;
  localB->x[5] = rtu_x6;
}

/*
 * System initialize for atomic system:
 *    '<S52>/MATLAB Function'
 *    '<S53>/MATLAB Function'
 *    '<S56>/MATLAB Function'
 *    '<S57>/MATLAB Function'
 *    '<S58>/MATLAB Function'
 *    '<S59>/MATLAB Function'
 */
void Human_in_Lo_MATLABFunction_Init(DW_MATLABFunction_Human_in_Lo_T *localDW)
{
  memset(&localDW->data_mem[0], 0, 500U * sizeof(real_T));
}

/*
 * Output and update for atomic system:
 *    '<S52>/MATLAB Function'
 *    '<S53>/MATLAB Function'
 *    '<S56>/MATLAB Function'
 *    '<S57>/MATLAB Function'
 *    '<S58>/MATLAB Function'
 *    '<S59>/MATLAB Function'
 */
void Human_in_Loop_MATLABFunction(real_T rtu_data,
  B_MATLABFunction_Human_in_Loo_T *localB, DW_MATLABFunction_Human_in_Lo_T
  *localDW)
{
  real_T y;
  real_T scale;
  real_T absxk;
  real_T t;
  real_T tmp[500];
  int32_T k;

  /* MATLAB Function 'Sensor Data/EMG module/Preprocessing10/MATLAB Function': '<S62>:1' */
  /* '<S62>:1:8' */
  memcpy(&tmp[0], &localDW->data_mem[1], 499U * sizeof(real_T));
  tmp[499] = rtu_data;

  /* '<S62>:1:10' */
  y = 0.0;
  scale = 3.3121686421112381E-170;
  for (k = 0; k < 500; k++) {
    localDW->data_mem[k] = tmp[k];
    absxk = localDW->data_mem[k] * localDW->data_mem[k];
    absxk = fabs(absxk);
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
  localB->y = y / 22.360679774997898;
}

/*
 * Output and update for atomic system:
 *    '<S30>/Mux'
 *    '<S33>/Mux'
 */
void Human_in_Loop_Mux(real_T rtu_x1, real_T rtu_x2, B_Mux_Human_in_Loop_T
  *localB)
{
  /* MATLAB Function 'Sensor Data/Encoder module/Mux': '<S77>:1' */
  /* '<S77>:1:3' */
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
static void Human_in_Loop_mldivide_g(const real_T A[16], real_T B[4])
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
  /* MATLAB Function 'Control Module/Torque track': '<S12>:1' */
  /* '<S12>:1:42' */
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
    /* '<S12>:1:40' */
    memset(&Human_in_Loop_B.torque_track[0], 0, 750U * sizeof(real_T));
    memset(&Human_in_Loop_B.torque_delta_track[0], 0, 750U * sizeof(real_T));

    /* '<S12>:1:41' */
    /* '<S12>:1:42' */
    /* '<S12>:1:43' */
    index_peak = floor(Human_in_Loop_B.RT3[2] / 100.0 * Human_in_Loop_B.RT1[2] *
                       500.0);

    /* '<S12>:1:44' */
    index_rise = index_peak - floor(Human_in_Loop_B.RT3[1] / 100.0 *
      Human_in_Loop_B.RT1[2] * 500.0);

    /* '<S12>:1:45' */
    index_fall = floor(Human_in_Loop_B.RT3[3] / 100.0 * Human_in_Loop_B.RT1[2] *
                       500.0) + index_peak;

    /* '<S12>:1:48' */
    /* '<S12>:1:52' */
    /* '<S12>:1:53' */
    parm1[0] = 0.0;
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
    Human_in_Loop_mldivide_g(tmp, parm1);
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

    /* '<S12>:1:54' */
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

    /* '<S12>:1:55' */
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.torque_delta_track[y + f] = Human_in_Loop_B.tmp_data_c[f] *
        500.0;
    }

    /* '<S12>:1:57' */
    /* '<S12>:1:61' */
    /* '<S12>:1:62' */
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
    Human_in_Loop_mldivide_g(tmp, parm1);
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

    /* '<S12>:1:63' */
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.torque_track[y + f] = Human_in_Loop_B.tmp_data_c[f];
    }

    /* '<S12>:1:67' */
    /* '<S12>:1:69' */
    for (f = 0; f < 750; f++) {
      Human_in_Loop_DW.TorqueMem[f << 2] = Human_in_Loop_B.torque_track[f];
      Human_in_Loop_DW.TorqueMem[2 + (f << 2)] =
        Human_in_Loop_B.torque_delta_track[f];
    }
  }

  /* '<S12>:1:73' */
  Human_in_Loop_DW.last_footstate = footstate;

  /* '<S12>:1:74' */
  Human_in_Loop_DW.TorqueMem[1 + (((int32_T)stride_index - 1) << 2)] =
    torque_measure;

  /* '<S12>:1:75' */
  Human_in_Loop_DW.TorqueMem[3 + (((int32_T)stride_index - 1) << 2)] =
    troque_delta;
  if (mode == 2.0) {
    /* '<S12>:1:77' */
    /* '<S12>:1:78' */
    mode = Human_in_Loop_DW.TorqueMem[((int32_T)stride_index - 1) << 2];

    /* '<S12>:1:79' */
    stride_index = Human_in_Loop_DW.TorqueMem[(((int32_T)stride_index - 1) << 2)
      + 2];
  } else {
    /* '<S12>:1:81' */
    mode = 0.0;

    /* '<S12>:1:82' */
    stride_index = 0.0;
  }

  /* '<S12>:1:85' */
  /* '<S12>:1:86' */
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

  /* MATLAB Function: '<S1>/LRN' */
  /* MATLAB Function 'Control Module/LRN': '<S10>:1' */
  /* '<S10>:1:15' */
  /* '<S10>:1:19' */
  /* '<S10>:1:22' */
  /* '<S10>:1:23' */
  /* '<S10>:1:25' */
  /* '<S10>:1:27' */
  /* '<S10>:1:29' */
  stride_index = Human_in_Loop_B.RT1[3] * 500.0 + 1.0;
  if (stride_index > 750.0) {
    /* '<S10>:1:30' */
    /* '<S10>:1:31' */
    stride_index = 750.0;
  }

  if ((Human_in_Loop_DW.last_footstate_a == 0.0) && (Human_in_Loop_B.RT1[1] ==
       1.0) && (Human_in_Loop_B.RT1[0] == 2.0) && Human_in_Loop_P.LRN_BT_LRN_ON)
  {
    /* '<S10>:1:34' */
    /* '<S10>:1:35' */
    /* '<S10>:1:36' */
    mode = 1.0 - Human_in_Loop_P.LRN_error_filter_k;

    /* '<S10>:1:37' */
    footstate = Human_in_Loop_B.RT2[2];
    for (f = 0; f < 750; f++) {
      peak_torque = (Human_in_Loop_DW.TorqueMem[f << 2] -
                     Human_in_Loop_DW.TorqueMem[(f << 2) + 1]) *
        Human_in_Loop_P.LRN_error_filter_k + mode *
        Human_in_Loop_DW.torque_error_memory[f];
      Human_in_Loop_DW.torque_error_memory[f] = peak_torque;
      Human_in_Loop_DW.lrn_cmd_memory[f] = Human_in_Loop_P.LRN_lrn_shrink *
        Human_in_Loop_DW.lrn_cmd_memory[f] + footstate *
        Human_in_Loop_DW.torque_error_memory[f];
    }
  }

  if (Human_in_Loop_P.LRN_BT_LRN_CLEAR) {
    /* '<S10>:1:40' */
    /* '<S10>:1:41' */
    /* '<S10>:1:42' */
    memset(&Human_in_Loop_DW.torque_error_memory[0], 0, 1000U * sizeof(real_T));
    memset(&Human_in_Loop_DW.lrn_cmd_memory[0], 0, 1000U * sizeof(real_T));
  }

  if (Human_in_Loop_B.RT1[0] == 2.0) {
    /* '<S10>:1:45' */
    b_n = Human_in_Loop_P.LRN_time_delay;
    if (b_n < -2147482897) {
      b_n = MAX_int32_T;
    } else {
      b_n = 750 - b_n;
    }

    if (stride_index >= b_n) {
      /* '<S10>:1:46' */
      /* '<S10>:1:47' */
      Human_in_Loop_B.lrn_cmd = 0.0;
    } else {
      /* '<S10>:1:49' */
      mode = rt_roundd_snf(stride_index + (real_T)Human_in_Loop_P.LRN_time_delay);
      if (mode < 2.147483648E+9) {
        if (mode >= -2.147483648E+9) {
          f = (int32_T)mode;
        } else {
          f = MIN_int32_T;
        }
      } else {
        f = MAX_int32_T;
      }

      Human_in_Loop_B.lrn_cmd = Human_in_Loop_DW.lrn_cmd_memory[f - 1];
    }
  } else {
    /* '<S10>:1:52' */
    Human_in_Loop_B.lrn_cmd = 0.0;
  }

  /* '<S10>:1:55' */
  Human_in_Loop_DW.last_footstate_a = Human_in_Loop_B.RT1[1];

  /* '<S10>:1:56' */
  memcpy(&Human_in_Loop_B.lrn_mem[0], &Human_in_Loop_DW.lrn_cmd_memory[0], 750U *
         sizeof(real_T));

  /* End of MATLAB Function: '<S1>/LRN' */

  /* MATLAB Function: '<S1>/Controller' */
  /* MATLAB Function 'Control Module/Controller': '<S9>:1' */
  /* '<S9>:1:21' */
  /* '<S9>:1:22' */
  /* '<S9>:1:26' */
  /* '<S9>:1:27' */
  /* '<S9>:1:31' */
  /* '<S9>:1:33' */
  /* '<S9>:1:36' */
  /* '<S9>:1:37' */
  /* '<S9>:1:41' */
  /* '<S9>:1:42' */
  /* '<S9>:1:43' */
  /* '<S9>:1:44' */
  switch ((int32_T)Human_in_Loop_B.RT1[0]) {
   case 1:
    /* '<S9>:1:48' */
    Human_in_Loop_B.motor_vel_cmd = 0.0;

    /* '<S9>:1:49' */
    Human_in_Loop_DW.calib_state = 0.0;
    break;

   case 3:
    /* '<S9>:1:52' */
    Human_in_Loop_B.motor_vel_cmd = -Human_in_Loop_P.Controller_SLACK_SPEED *
      5.0;
    break;

   case 4:
    if ((Human_in_Loop_DW.calib_state == 0.0) && (Human_in_Loop_B.RT4[0] <
         Human_in_Loop_P.Controller_CALIB_TORQUE)) {
      /* '<S9>:1:55' */
      /* '<S9>:1:56' */
      Human_in_Loop_B.motor_vel_cmd = Human_in_Loop_P.Controller_CALIB_SPEED *
        5.0;
    } else if ((Human_in_Loop_DW.calib_state == 0.0) && (Human_in_Loop_B.RT4[0] >
                Human_in_Loop_P.Controller_CALIB_TORQUE)) {
      /* '<S9>:1:57' */
      /* '<S9>:1:58' */
      Human_in_Loop_DW.calib_state = 1.0;

      /* '<S9>:1:59' */
      Human_in_Loop_B.motor_vel_cmd = 0.0;
    } else {
      /* '<S9>:1:61' */
      Human_in_Loop_B.motor_vel_cmd = 0.0;
    }
    break;

   case 2:
    /* '<S9>:1:65' */
    switch (Human_in_Loop_P.Controller_run_mode) {
     case 1:
      if (Human_in_Loop_B.RT5[0] > 0.0) {
        /* '<S9>:1:68' */
        /* '<S9>:1:69' */
        stride_index = (Human_in_Loop_B.RT5[0] - Human_in_Loop_B.RT5[0] / 50.0 *
                        Human_in_Loop_P.Controller_FOLLOW_SLACK_ANGLE) -
          Human_in_Loop_B.RT6[0] * 0.33333333333333331;
      } else {
        /* '<S9>:1:71' */
        stride_index = Human_in_Loop_B.RT5[0] - Human_in_Loop_B.RT6[0] *
          0.33333333333333331;
      }

      /* '<S9>:1:73' */
      Human_in_Loop_B.motor_vel_cmd = Human_in_Loop_B.RT2[3] * stride_index *
        5.0 / 0.05;
      break;

     case 2:
      if (Human_in_Loop_B.RT1[1] == 1.0) {
        /* '<S9>:1:76' */
        /* '<S9>:1:77' */
        /* '<S9>:1:78' */
        /* '<S9>:1:79' */
        /* '<S9>:1:80' */
        Human_in_Loop_B.motor_vel_cmd = ((((Human_in_Loop_B.torque_des -
          Human_in_Loop_B.RT4[0]) * Human_in_Loop_B.RT2[0] +
          (Human_in_Loop_B.torque_delta_des - Human_in_Loop_B.RT4[1]) *
          Human_in_Loop_B.RT2[1]) + Human_in_Loop_B.lrn_cmd) +
          Human_in_Loop_B.RT2[4] * Human_in_Loop_B.torque_delta_des) * 5.0 /
          0.05;
      } else {
        /* '<S9>:1:82' */
        /* '<S9>:1:83' */
        Human_in_Loop_B.motor_vel_cmd = (0.0 - Human_in_Loop_B.RT6[0] *
          0.33333333333333331) * Human_in_Loop_B.RT2[3] * 5.0 / 0.05;
      }
      break;

     default:
      /* '<S9>:1:86' */
      Human_in_Loop_B.motor_vel_cmd = 0.0;
      break;
    }
    break;

   case 0:
    /* '<S9>:1:90' */
    Human_in_Loop_B.motor_vel_cmd = 0.0;
    break;

   default:
    /* '<S9>:1:93' */
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
  Human_in_Loop_B.Gain2_j = Human_in_Loop_P.Gain2_Gain * Human_in_Loop_B.vel;

  /* Gain: '<S11>/Gain1' */
  Human_in_Loop_B.Gain1_pi = Human_in_Loop_P.Gain1_Gain *
    Human_in_Loop_B.Gain2_j;

  /* S-Function (rti_commonblock): '<S13>/S-Function1' */
  /* This comment workarounds a code generation problem */

  /* --- Human_in_Loop/Control Module/Motor/DAC_CLASS1_BL1 --- */
  /* --- [RTI120X, DAC C1] - Channel: 16 --- */
  {
    /* define variables required for DAC realtime functions */
    Float64 inportDacData= 0.0;
    inportDacData = (real_T) Human_in_Loop_B.Gain1_pi;

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

/* Function for MATLAB Function: '<S33>/MATLAB Function' */
static real_T Human_in_Loop_xnrm2_o(const real_T x[30], int32_T ix0)
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

/* Function for MATLAB Function: '<S33>/MATLAB Function' */
static real_T Human_in_Loop_xnrm2_ov(int32_T n, const real_T x[30], int32_T ix0)
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

/* Function for MATLAB Function: '<S33>/MATLAB Function' */
static void Human_in_Loop_xgeqp3_n(real_T A[30], real_T tau[2], int32_T jpvt[2])
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
  smax = Human_in_Loop_xnrm2_o(A, 1);
  vn2[0] = smax;
  vn1[0] = smax;
  jpvt[1] = 2;
  work[1] = 0.0;
  smax = Human_in_Loop_xnrm2_o(A, 16);
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
    xnorm = Human_in_Loop_xnrm2_ov(14 - k, A, i_i + 2);
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

        xnorm = rt_hypotd_snf(smax, Human_in_Loop_xnrm2_ov(14 - k, A, i_i + 2));
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
          vn1[1] = Human_in_Loop_xnrm2_ov(14 - k, A, k + 17);
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
  real_T x[6];
  real_T varargin_1[2];
  int32_T ixstart;
  real_T x_0[2];
  real_T A[30];
  int32_T rankR;
  real_T tol;
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
    /* S-Function (rti_commonblock): '<S96>/S-Function1' */
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

    /* Gain: '<S33>/Gain' */
    Human_in_Loop_B.Gain = Human_in_Loop_P.Gain_Gain *
      Human_in_Loop_B.SFunction1;

    /* DiscreteFilter: '<S33>/0.4low2' */
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

    /* MATLAB Function: '<S33>/Data process' */
    /* MATLAB Function 'Sensor Data/Torque module/Data process': '<S97>:1' */
    if (Human_in_Loop_P.Dataprocess_BT_RESET_TORQUE) {
      /* '<S97>:1:10' */
      /* '<S97>:1:11' */
      Human_in_Loop_DW.torque_zero = Human_in_Loop_B.u4low2 *
        Human_in_Loop_P.Dataprocess_load_vol_gain +
        Human_in_Loop_P.Dataprocess_load_vol_offset;
    }

    /* '<S97>:1:14' */
    Human_in_Loop_B.torque = (Human_in_Loop_B.u4low2 *
      Human_in_Loop_P.Dataprocess_load_vol_gain +
      Human_in_Loop_P.Dataprocess_load_vol_offset) -
      Human_in_Loop_DW.torque_zero;

    /* End of MATLAB Function: '<S33>/Data process' */

    /* UnitDelay: '<S95>/Unit Delay1' */
    Human_in_Loop_B.x2k1 = Human_in_Loop_DW.UnitDelay1_DSTATE;

    /* UnitDelay: '<S95>/Unit Delay' */
    Human_in_Loop_B.x1k1 = Human_in_Loop_DW.UnitDelay_DSTATE;

    /* Gain: '<S95>/Gain1' */
    B_0 = Human_in_Loop_P.uOrderTD_T1 * Human_in_Loop_P.uOrderTD_T2;
    tol = 1.0 / B_0;
    Human_in_Loop_B.Gain1 = tol * Human_in_Loop_B.x1k1;

    /* Gain: '<S95>/Gain2' */
    tol = Human_in_Loop_P.uOrderTD_T1 + Human_in_Loop_P.uOrderTD_T2;
    B_0 = Human_in_Loop_P.uOrderTD_T1 * Human_in_Loop_P.uOrderTD_T2;
    tol /= B_0;
    Human_in_Loop_B.Gain2 = tol * Human_in_Loop_B.x2k1;

    /* UnitDelay: '<S95>/Unit Delay2' */
    Human_in_Loop_B.UnitDelay2 = Human_in_Loop_DW.UnitDelay2_DSTATE;

    /* Gain: '<S95>/Gain4' */
    B_0 = Human_in_Loop_P.uOrderTD_T1 * Human_in_Loop_P.uOrderTD_T2;
    tol = 1.0 / B_0;
    Human_in_Loop_B.Gain4 = tol * Human_in_Loop_B.UnitDelay2;

    /* Sum: '<S95>/Add2' */
    Human_in_Loop_B.Add2 = (Human_in_Loop_B.Gain1 + Human_in_Loop_B.Gain2) -
      Human_in_Loop_B.Gain4;

    /* Gain: '<S95>/Gain3' */
    Human_in_Loop_B.Gain3 = Human_in_Loop_P.uOrderTD_Ts * Human_in_Loop_B.Add2;

    /* Sum: '<S95>/Add1' */
    Human_in_Loop_B.x2k = Human_in_Loop_B.x2k1 - Human_in_Loop_B.Gain3;

    /* MATLAB Function: '<S33>/Mux' */
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

    /* S-Function (rti_commonblock): '<S73>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* Gain: '<S30>/Gain' */
    Human_in_Loop_B.Gain_l = Human_in_Loop_P.Gain_Gain_c *
      Human_in_Loop_B.SFunction1_o1;

    /* MATLAB Function: '<S30>/Data process' */
    /* MATLAB Function 'Sensor Data/Encoder module/Data process': '<S70>:1' */
    if (Human_in_Loop_P.Dataprocess_BT_RESET_ANKLE) {
      /* '<S70>:1:8' */
      /* '<S70>:1:9' */
      Human_in_Loop_DW.angle_zero_f = Human_in_Loop_B.Gain_l;
    }

    /* '<S70>:1:12' */
    Human_in_Loop_B.angle_m = Human_in_Loop_B.Gain_l -
      Human_in_Loop_DW.angle_zero_f;

    /* End of MATLAB Function: '<S30>/Data process' */

    /* Gain: '<S30>/Gain1' */
    Human_in_Loop_B.Gain1_b = Human_in_Loop_P.Gain1_Gain_e *
      Human_in_Loop_B.SFunction1_o2;

    /* MATLAB Function: '<S30>/Mux' */
    Human_in_Loop_Mux(Human_in_Loop_B.angle_m, Human_in_Loop_B.Gain1_b,
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

    /* S-Function (rti_commonblock): '<S74>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* Gain: '<S30>/Gain2' */
    Human_in_Loop_B.Gain2_h = Human_in_Loop_P.Gain2_Gain_e *
      Human_in_Loop_B.SFunction1_o1_k;

    /* MATLAB Function: '<S30>/Data process1' */
    /* MATLAB Function 'Sensor Data/Encoder module/Data process1': '<S71>:1' */
    if (Human_in_Loop_P.Dataprocess1_BT_RESET_MOTOR) {
      /* '<S71>:1:8' */
      /* '<S71>:1:9' */
      Human_in_Loop_DW.angle_zero = Human_in_Loop_B.Gain2_h;
    }

    /* '<S71>:1:12' */
    Human_in_Loop_B.angle = Human_in_Loop_B.Gain2_h -
      Human_in_Loop_DW.angle_zero;

    /* End of MATLAB Function: '<S30>/Data process1' */

    /* Gain: '<S30>/Gain3' */
    Human_in_Loop_B.Gain3_m = Human_in_Loop_P.Gain3_Gain *
      Human_in_Loop_B.SFunction1_o2_k;

    /* UnitDelay: '<S69>/Unit Delay1' */
    Human_in_Loop_B.x2k1_k = Human_in_Loop_DW.UnitDelay1_DSTATE_h;

    /* UnitDelay: '<S69>/Unit Delay' */
    Human_in_Loop_B.x1k1_m = Human_in_Loop_DW.UnitDelay_DSTATE_e;

    /* Gain: '<S69>/Gain1' */
    B_0 = Human_in_Loop_P.uOrderTD_T1_l * Human_in_Loop_P.uOrderTD_T2_h;
    tol = 1.0 / B_0;
    Human_in_Loop_B.Gain1_m = tol * Human_in_Loop_B.x1k1_m;

    /* Gain: '<S69>/Gain2' */
    tol = Human_in_Loop_P.uOrderTD_T1_l + Human_in_Loop_P.uOrderTD_T2_h;
    B_0 = Human_in_Loop_P.uOrderTD_T1_l * Human_in_Loop_P.uOrderTD_T2_h;
    tol /= B_0;
    Human_in_Loop_B.Gain2_c = tol * Human_in_Loop_B.x2k1_k;

    /* UnitDelay: '<S69>/Unit Delay2' */
    Human_in_Loop_B.UnitDelay2_b = Human_in_Loop_DW.UnitDelay2_DSTATE_a;

    /* Gain: '<S69>/Gain4' */
    B_0 = Human_in_Loop_P.uOrderTD_T1_l * Human_in_Loop_P.uOrderTD_T2_h;
    tol = 1.0 / B_0;
    Human_in_Loop_B.Gain4_g = tol * Human_in_Loop_B.UnitDelay2_b;

    /* Sum: '<S69>/Add2' */
    Human_in_Loop_B.Add2_m = (Human_in_Loop_B.Gain1_m + Human_in_Loop_B.Gain2_c)
      - Human_in_Loop_B.Gain4_g;

    /* Gain: '<S69>/Gain3' */
    Human_in_Loop_B.Gain3_k = Human_in_Loop_P.uOrderTD_Ts_n *
      Human_in_Loop_B.Add2_m;

    /* Sum: '<S69>/Add1' */
    Human_in_Loop_B.x2k_i = Human_in_Loop_B.x2k1_k - Human_in_Loop_B.Gain3_k;

    /* MATLAB Function: '<S30>/Mux1' */
    /* MATLAB Function 'Sensor Data/Encoder module/Mux1': '<S78>:1' */
    /* '<S78>:1:3' */
    Human_in_Loop_B.x_i[0] = Human_in_Loop_B.angle;
    Human_in_Loop_B.x_i[1] = Human_in_Loop_B.Gain3_m;
    Human_in_Loop_B.x_i[2] = Human_in_Loop_B.x2k_i;

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
      Human_in_Loop_DW.RT6_Buffer1[0] = Human_in_Loop_B.x_i[0];
      Human_in_Loop_DW.RT6_Buffer1[1] = Human_in_Loop_B.x_i[1];
      Human_in_Loop_DW.RT6_Buffer1[2] = Human_in_Loop_B.x_i[2];
    } else {
      Human_in_Loop_DW.RT6_Buffer0[0] = Human_in_Loop_B.x_i[0];
      Human_in_Loop_DW.RT6_Buffer0[1] = Human_in_Loop_B.x_i[1];
      Human_in_Loop_DW.RT6_Buffer0[2] = Human_in_Loop_B.x_i[2];
    }

    Human_in_Loop_DW.RT6_last_buf_wr = Human_in_Loop_DW.RT6_write_buf;
    Human_in_Loop_DW.RT6_write_buf = -1;

    /* End of RateTransition: '<Root>/RT6' */

    /* S-Function (rti_commonblock): '<S80>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* MATLAB Function: '<S31>/FootSwitch Filter' */
    /* MATLAB Function 'Sensor Data/FootSwitch module/FootSwitch Filter': '<S82>:1' */
    /* '<S82>:1:6' */
    if (Human_in_Loop_DW.foot_state == 0.0) {
      /* '<S82>:1:14' */
      if (Human_in_Loop_B.SFunction1_a) {
        /* '<S82>:1:15' */
        /* '<S82>:1:16' */
        Human_in_Loop_DW.foot_state = 1.0;
      } else {
        /* '<S82>:1:18' */
        Human_in_Loop_DW.foot_state = 0.0;
      }
    } else if (Human_in_Loop_B.SFunction1_a) {
      /* '<S82>:1:21' */
      /* '<S82>:1:22' */
      Human_in_Loop_DW.foot_state = 1.0;
    } else {
      /* '<S82>:1:24' */
      Human_in_Loop_DW.filter_time += 0.0002;
      if (Human_in_Loop_DW.filter_time > 0.4) {
        /* '<S82>:1:25' */
        /* '<S82>:1:26' */
        Human_in_Loop_DW.filter_time = 0.0;

        /* '<S82>:1:27' */
        Human_in_Loop_DW.foot_state = 0.0;
      }
    }

    /* '<S82>:1:32' */
    Human_in_Loop_B.state_c = Human_in_Loop_DW.foot_state;

    /* End of MATLAB Function: '<S31>/FootSwitch Filter' */

    /* MATLAB Function: '<S7>/State Machine' */
    /* MATLAB Function 'State Module/State Machine': '<S101>:1' */
    /* '<S101>:1:20' */
    /* '<S101>:1:21' */
    /* '<S101>:1:22' */
    /* '<S101>:1:23' */
    /* '<S101>:1:26' */
    /* '<S101>:1:27' */
    /* '<S101>:1:28' */
    /* '<S101>:1:29' */
    /* '<S101>:1:30' */
    if (Human_in_Loop_P.StateMachine_BT_RUN) {
      /* '<S101>:1:33' */
      /* '<S101>:1:34' */
      Human_in_Loop_DW.bt_run = 1.0;
    }

    if (Human_in_Loop_P.StateMachine_BT_CALIB) {
      /* '<S101>:1:36' */
      /* '<S101>:1:37' */
      Human_in_Loop_DW.reg_mode = 4.0;
    }

    if (Human_in_Loop_P.StateMachine_BT_SLACK) {
      /* '<S101>:1:39' */
      /* '<S101>:1:40' */
      Human_in_Loop_DW.reg_mode = 3.0;
    }

    if (Human_in_Loop_P.StateMachine_BT_IDLE) {
      /* '<S101>:1:42' */
      /* '<S101>:1:43' */
      Human_in_Loop_DW.reg_mode = 1.0;
    }

    if (Human_in_Loop_P.StateMachine_BT_ERROR) {
      /* '<S101>:1:45' */
      /* '<S101>:1:46' */
      Human_in_Loop_DW.reg_mode = 0.0;
    }

    if ((Human_in_Loop_DW.bt_run == 1.0) && (Human_in_Loop_DW.reg_last_switch ==
         0.0) && (Human_in_Loop_B.state_c == 1.0)) {
      /* '<S101>:1:49' */
      /* '<S101>:1:50' */
      Human_in_Loop_DW.reg_mode = 2.0;

      /* '<S101>:1:51' */
      Human_in_Loop_DW.bt_run = 0.0;
    }

    if ((Human_in_Loop_DW.reg_mode == 2.0) || (Human_in_Loop_DW.reg_mode == 1.0))
    {
      /* '<S101>:1:54' */
      if ((Human_in_Loop_DW.reg_last_switch == 0.0) && (Human_in_Loop_B.state_c ==
           1.0)) {
        /* '<S101>:1:55' */
        /* '<S101>:1:56' */
        Human_in_Loop_DW.reg_state = 1.0;

        /* '<S101>:1:57' */
        Human_in_Loop_DW.reg_stride_time = 0.618 *
          Human_in_Loop_DW.reg_stride_time + 0.382 *
          Human_in_Loop_DW.reg_stride_time_count;

        /* '<S101>:1:58' */
        Human_in_Loop_DW.reg_stride_time_count = 0.0;
      } else if ((Human_in_Loop_DW.reg_state == 1.0) &&
                 (Human_in_Loop_DW.reg_stride_time_count > 0.65 *
                  Human_in_Loop_DW.reg_stride_time)) {
        /* '<S101>:1:59' */
        /* '<S101>:1:60' */
        Human_in_Loop_DW.reg_state = 0.0;

        /* '<S101>:1:61' */
        Human_in_Loop_DW.reg_stride_time_count += 0.0002;
      } else {
        /* '<S101>:1:63' */
        Human_in_Loop_DW.reg_stride_time_count += 0.0002;
      }
    }

    /* '<S101>:1:67' */
    Human_in_Loop_DW.reg_last_switch = Human_in_Loop_B.state_c;
    if (Human_in_Loop_DW.reg_stride_time > 1.5) {
      /* '<S101>:1:68' */
      /* '<S101>:1:69' */
      Human_in_Loop_DW.reg_stride_time = 1.5;
    } else {
      if (Human_in_Loop_DW.reg_stride_time < 0.5) {
        /* '<S101>:1:70' */
        /* '<S101>:1:71' */
        Human_in_Loop_DW.reg_stride_time = 0.5;
      }
    }

    /* '<S101>:1:74' */
    Human_in_Loop_B.mode = Human_in_Loop_DW.reg_mode;

    /* '<S101>:1:75' */
    Human_in_Loop_B.state = Human_in_Loop_DW.reg_state;

    /* '<S101>:1:76' */
    Human_in_Loop_B.stride_time = Human_in_Loop_DW.reg_stride_time;

    /* '<S101>:1:77' */
    Human_in_Loop_B.stride_timer = Human_in_Loop_DW.reg_stride_time_count;

    /* End of MATLAB Function: '<S7>/State Machine' */

    /* MATLAB Function: '<S7>/Mux1' */
    /* MATLAB Function 'State Module/Mux1': '<S100>:1' */
    /* '<S100>:1:3' */
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
     *  Constant: '<S21>/Ks'
     */
    /* MATLAB Function 'Parameter Module/Control Parameter/Mux1': '<S23>:1' */
    /* '<S23>:1:3' */
    Human_in_Loop_B.x_o[0] = Human_in_Loop_P.Kp_Value;
    Human_in_Loop_B.x_o[1] = Human_in_Loop_P.Kd_Value;
    Human_in_Loop_B.x_o[2] = Human_in_Loop_P.Kl_Value;
    Human_in_Loop_B.x_o[3] = Human_in_Loop_P.Ks_Value;
    Human_in_Loop_B.x_o[4] = Human_in_Loop_P.Ko_Value;

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
      for (i = 0; i < 5; i++) {
        Human_in_Loop_DW.RT2_Buffer1[i] = Human_in_Loop_B.x_o[i];
      }
    } else {
      for (i = 0; i < 5; i++) {
        Human_in_Loop_DW.RT2_Buffer0[i] = Human_in_Loop_B.x_o[i];
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

  /* StateSpace: '<S50>/low_pass' */
  Human_in_Loop_B.low_pass = 0.0;
  Human_in_Loop_B.low_pass += Human_in_Loop_P.low_pass_C *
    Human_in_Loop_X.low_pass_CSTATE[1];

  /* StateSpace: '<S51>/low_pass' */
  Human_in_Loop_B.low_pass_f = 0.0;
  Human_in_Loop_B.low_pass_f += Human_in_Loop_P.low_pass_C_g *
    Human_in_Loop_X.low_pass_CSTATE_o[1];

  /* StateSpace: '<S54>/low_pass' */
  Human_in_Loop_B.low_pass_l = 0.0;
  Human_in_Loop_B.low_pass_l += Human_in_Loop_P.low_pass_C_d *
    Human_in_Loop_X.low_pass_CSTATE_n[1];

  /* StateSpace: '<S55>/low_pass' */
  Human_in_Loop_B.low_pass_d = 0.0;
  Human_in_Loop_B.low_pass_d += Human_in_Loop_P.low_pass_C_i *
    Human_in_Loop_X.low_pass_CSTATE_p[1];

  /* StateSpace: '<S61>/low_pass' */
  Human_in_Loop_B.low_pass_k = 0.0;
  Human_in_Loop_B.low_pass_k += Human_in_Loop_P.low_pass_C_a *
    Human_in_Loop_X.low_pass_CSTATE_on[1];

  /* StateSpace: '<S60>/low_pass' */
  Human_in_Loop_B.low_pass_fc = 0.0;
  Human_in_Loop_B.low_pass_fc += Human_in_Loop_P.low_pass_C_c *
    Human_in_Loop_X.low_pass_CSTATE_f[1];

  /* MATLAB Function: '<S29>/Mux1' */
  Human_in_Loop_Mux1(Human_in_Loop_B.low_pass, Human_in_Loop_B.low_pass_f,
                     Human_in_Loop_B.low_pass_l, Human_in_Loop_B.low_pass_d,
                     Human_in_Loop_B.low_pass_k, Human_in_Loop_B.low_pass_fc,
                     &Human_in_Loop_B.sf_Mux1_c);
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
    tol = Human_in_Loop_B.Delay7;

    /* '<S20>:1:17' */
    Human_in_Loop_B.Cycle_Frequency = Human_in_Loop_B.Delay6;

    /* '<S20>:1:19' */
    for (i = 0; i < 6; i++) {
      /* '<S20>:1:19' */
      /* '<S20>:1:20' */
      varargin_1[0] = Human_in_Loop_DW.EMG_Memory[i];
      varargin_1[1] = Human_in_Loop_B.sf_Mux1_c.x[i];
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
      Human_in_Loop_DW.EMG_Memory[6 + i] += Human_in_Loop_B.sf_Mux1_c.x[i];

      /* '<S20>:1:22' */
      Human_in_Loop_DW.EMG_Memory[12 + i] += Human_in_Loop_B.sf_Mux1_c.x[i] *
        Human_in_Loop_B.sf_Mux1_c.x[i];

      /* '<S20>:1:23' */
      Human_in_Loop_DW.EMG_Memory[18 + i] += fabs(Human_in_Loop_B.sf_Mux1_c.x[i]);
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
        tol = Human_in_Loop_DW.EMG_Memory[12 + i] / Human_in_Loop_B.Delay3;
        tol = sqrt(tol);
        x[i] = tol;
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
      tol = Human_in_Loop_B.Delay3 * Human_in_Loop_P.SingleCycleAnalysis1_Ts;

      /* '<S20>:1:38' */
      Human_in_Loop_B.Cycle_Frequency = 1.0 / tol;
    }

    Human_in_Loop_B.Cycle_Time = tol;

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
        tol = Human_in_Loop_B.UnitDelay_a + 1.0;

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
          tol = 0.0;
        }

        Human_in_Loop_B.count = tol;

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
      tol = floor(Human_in_Loop_B.x[2] * 500.0);
      if (1.0 > tol) {
        c = 0;
      } else {
        c = (int32_T)tol;
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

    /* Gain: '<S95>/Gain' */
    Human_in_Loop_B.Gain_g = Human_in_Loop_P.uOrderTD_Ts * Human_in_Loop_B.x2k1;

    /* Sum: '<S95>/Add' */
    Human_in_Loop_B.x1k = Human_in_Loop_B.Gain_g + Human_in_Loop_B.x1k1;

    /* MATLAB Function: '<S33>/MATLAB Function' */
    /* MATLAB Function 'Sensor Data/Torque module/MATLAB Function': '<S98>:1' */
    /* '<S98>:1:13' */
    /* '<S98>:1:8' */
    memcpy(&B[0], &Human_in_Loop_DW.data[1], 14U * sizeof(real_T));
    B[14] = Human_in_Loop_B.torque;
    memcpy(&Human_in_Loop_DW.data[0], &B[0], 15U * sizeof(real_T));

    /* '<S98>:1:13' */
    for (ixstart = 0; ixstart < 30; ixstart++) {
      A[ixstart] = b_A[ixstart];
    }

    Human_in_Loop_xgeqp3_n(A, varargin_1, tmp_size);
    rankR = 0;
    tol = 15.0 * fabs(A[0]) * 2.2204460492503131E-16;
    while ((rankR < 2) && (!(fabs(A[15 * rankR + rankR]) <= tol))) {
      rankR++;
    }

    x_0[0] = 0.0;
    x_0[1] = 0.0;
    memcpy(&B[0], &Human_in_Loop_DW.data[0], 15U * sizeof(real_T));
    if (varargin_1[0] != 0.0) {
      tol = B[0];
      for (i = 1; i + 1 < 16; i++) {
        tol += A[i] * B[i];
      }

      tol *= varargin_1[0];
      if (tol != 0.0) {
        B[0] -= tol;
        for (i = 1; i + 1 < 16; i++) {
          B[i] -= A[i] * tol;
        }
      }
    }

    if (varargin_1[1] != 0.0) {
      tol = B[1];
      for (i = 2; i + 1 < 16; i++) {
        tol += A[i + 15] * B[i];
      }

      tol *= varargin_1[1];
      if (tol != 0.0) {
        B[1] -= tol;
        for (i = 2; i + 1 < 16; i++) {
          B[i] -= A[i + 15] * tol;
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

    /* '<S98>:1:15' */
    Human_in_Loop_B.torque_dot = x_0[1] * 5000.0;

    /* End of MATLAB Function: '<S33>/MATLAB Function' */

    /* SampleTimeMath: '<S72>/TSamp'
     *
     * About '<S72>/TSamp':
     *  y = u * K where K = 1 / ( w * Ts )
     */
    Human_in_Loop_B.TSamp = Human_in_Loop_B.Gain3_m * Human_in_Loop_P.TSamp_WtEt;

    /* UnitDelay: '<S72>/UD' */
    Human_in_Loop_B.Uk1 = Human_in_Loop_DW.UD_DSTATE;

    /* Sum: '<S72>/Diff' */
    Human_in_Loop_B.Diff = Human_in_Loop_B.TSamp - Human_in_Loop_B.Uk1;

    /* DiscreteFilter: '<S30>/0.4low1' */
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

    /* UnitDelay: '<S68>/Unit Delay' */
    Human_in_Loop_B.UnitDelay = Human_in_Loop_DW.UnitDelay_DSTATE_j;

    /* UnitDelay: '<S68>/Unit Delay1' */
    Human_in_Loop_B.UnitDelay1 = Human_in_Loop_DW.UnitDelay1_DSTATE_a;

    /* Sum: '<S68>/Add1' */
    Human_in_Loop_B.Add1 = Human_in_Loop_B.UnitDelay -
      Human_in_Loop_B.UnitDelay1;

    /* Gain: '<S68>/Gain' */
    Human_in_Loop_B.Gain_p = Human_in_Loop_P.uOrderTD_Ts_g *
      Human_in_Loop_B.Add1;

    /* Gain: '<S68>/Gain1' */
    B_0 = Human_in_Loop_P.uOrderTD_T;
    tol = 1.0 / B_0;
    Human_in_Loop_B.Gain1_d = tol * Human_in_Loop_B.Gain_p;

    /* Sum: '<S68>/Add2' */
    Human_in_Loop_B.Add2_i = Human_in_Loop_B.UnitDelay1 +
      Human_in_Loop_B.Gain1_d;

    /* Sum: '<S68>/Add3' */
    Human_in_Loop_B.Add3 = Human_in_Loop_B.Gain3_m - Human_in_Loop_B.Add2_i;

    /* Gain: '<S68>/Gain2' */
    B_0 = Human_in_Loop_P.uOrderTD_T;
    tol = 1.0 / B_0;
    Human_in_Loop_B.Gain2_i = tol * Human_in_Loop_B.Add3;

    /* Gain: '<S69>/Gain' */
    Human_in_Loop_B.Gain_io = Human_in_Loop_P.uOrderTD_Ts_n *
      Human_in_Loop_B.x2k1_k;

    /* Sum: '<S69>/Add' */
    Human_in_Loop_B.x1k_i = Human_in_Loop_B.Gain_io + Human_in_Loop_B.x1k1_m;

    /* MATLAB Function: '<S30>/Mux2' */
    /* MATLAB Function 'Sensor Data/Encoder module/Mux2': '<S79>:1' */
    /* '<S79>:1:3' */
    Human_in_Loop_B.x_d[0] = Human_in_Loop_B.u4low1;
    Human_in_Loop_B.x_d[1] = Human_in_Loop_B.Gain2_i;
    Human_in_Loop_B.x_d[2] = Human_in_Loop_B.x2k_i;

    /* S-Function (rti_commonblock): '<S75>/S-Function1' incorporates:
     *  Constant: '<S30>/VCC1'
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

    /* S-Function (rti_commonblock): '<S76>/S-Function1' incorporates:
     *  Constant: '<S30>/VCC3'
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

    /* S-Function (rti_commonblock): '<S81>/S-Function1' incorporates:
     *  Constant: '<S31>/Constant'
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

    /* S-Function (rti_commonblock): '<S42>/S-Function1' */
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
        (real_T*) &Human_in_Loop_B.SFunction1_o);
    }

    /* S-Function (rti_commonblock): '<S43>/S-Function1' */
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
        (real_T*) &Human_in_Loop_B.SFunction1_k);
    }

    /* S-Function (rti_commonblock): '<S44>/S-Function1' */
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
        (real_T*) &Human_in_Loop_B.SFunction1_j);
    }

    /* S-Function (rti_commonblock): '<S45>/S-Function1' */
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
        (real_T*) &Human_in_Loop_B.SFunction1_d);
    }

    /* S-Function (rti_commonblock): '<S46>/S-Function1' */
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

    /* S-Function (rti_commonblock): '<S47>/S-Function1' */
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
        (real_T*) &Human_in_Loop_B.SFunction1_i);
    }

    /* Gain: '<S29>/Gain' */
    Human_in_Loop_B.Gain_k = Human_in_Loop_P.Gain_Gain_h *
      Human_in_Loop_B.SFunction1_o;

    /* Gain: '<S29>/Gain1' */
    Human_in_Loop_B.Gain1_p = Human_in_Loop_P.Gain1_Gain_b *
      Human_in_Loop_B.SFunction1_k;

    /* Gain: '<S29>/Gain2' */
    Human_in_Loop_B.Gain2_cz = Human_in_Loop_P.Gain2_Gain_k *
      Human_in_Loop_B.SFunction1_j;

    /* Gain: '<S29>/Gain3' */
    Human_in_Loop_B.Gain3_f = Human_in_Loop_P.Gain3_Gain_b *
      Human_in_Loop_B.SFunction1_d;

    /* Gain: '<S29>/Gain4' */
    Human_in_Loop_B.Gain4_h = Human_in_Loop_P.Gain4_Gain *
      Human_in_Loop_B.SFunction1_l;

    /* Gain: '<S29>/Gain5' */
    Human_in_Loop_B.Gain5 = Human_in_Loop_P.Gain5_Gain *
      Human_in_Loop_B.SFunction1_i;
  }

  /* StateSpace: '<S56>/high_pass' */
  Human_in_Loop_B.high_pass = 0.0;
  Human_in_Loop_B.high_pass += Human_in_Loop_P.high_pass_C *
    Human_in_Loop_X.high_pass_CSTATE[1];

  /* Abs: '<S56>/Abs' */
  Human_in_Loop_B.Abs = fabs(Human_in_Loop_B.high_pass);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* MATLAB Function: '<S56>/MATLAB Function' */
    Human_in_Loop_MATLABFunction(Human_in_Loop_B.Abs,
      &Human_in_Loop_B.sf_MATLABFunction_g,
      &Human_in_Loop_DW.sf_MATLABFunction_g);
  }

  /* StateSpace: '<S57>/high_pass' */
  Human_in_Loop_B.high_pass_b = 0.0;
  Human_in_Loop_B.high_pass_b += Human_in_Loop_P.high_pass_C_i *
    Human_in_Loop_X.high_pass_CSTATE_b[1];

  /* Abs: '<S57>/Abs' */
  Human_in_Loop_B.Abs_k = fabs(Human_in_Loop_B.high_pass_b);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* MATLAB Function: '<S57>/MATLAB Function' */
    Human_in_Loop_MATLABFunction(Human_in_Loop_B.Abs_k,
      &Human_in_Loop_B.sf_MATLABFunction_m,
      &Human_in_Loop_DW.sf_MATLABFunction_m);
  }

  /* StateSpace: '<S58>/high_pass' */
  Human_in_Loop_B.high_pass_d = 0.0;
  Human_in_Loop_B.high_pass_d += Human_in_Loop_P.high_pass_C_m *
    Human_in_Loop_X.high_pass_CSTATE_c[1];

  /* Abs: '<S58>/Abs' */
  Human_in_Loop_B.Abs_b = fabs(Human_in_Loop_B.high_pass_d);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* MATLAB Function: '<S58>/MATLAB Function' */
    Human_in_Loop_MATLABFunction(Human_in_Loop_B.Abs_b,
      &Human_in_Loop_B.sf_MATLABFunction_n,
      &Human_in_Loop_DW.sf_MATLABFunction_n);
  }

  /* StateSpace: '<S59>/high_pass' */
  Human_in_Loop_B.high_pass_o = 0.0;
  Human_in_Loop_B.high_pass_o += Human_in_Loop_P.high_pass_C_j *
    Human_in_Loop_X.high_pass_CSTATE_m[1];

  /* Abs: '<S59>/Abs' */
  Human_in_Loop_B.Abs_n = fabs(Human_in_Loop_B.high_pass_o);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* MATLAB Function: '<S59>/MATLAB Function' */
    Human_in_Loop_MATLABFunction(Human_in_Loop_B.Abs_n,
      &Human_in_Loop_B.sf_MATLABFunction_l,
      &Human_in_Loop_DW.sf_MATLABFunction_l);
  }

  /* StateSpace: '<S52>/high_pass' */
  Human_in_Loop_B.high_pass_l = 0.0;
  Human_in_Loop_B.high_pass_l += Human_in_Loop_P.high_pass_C_mg *
    Human_in_Loop_X.high_pass_CSTATE_a[1];

  /* Abs: '<S52>/Abs' */
  Human_in_Loop_B.Abs_m = fabs(Human_in_Loop_B.high_pass_l);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* MATLAB Function: '<S52>/MATLAB Function' */
    Human_in_Loop_MATLABFunction(Human_in_Loop_B.Abs_m,
      &Human_in_Loop_B.sf_MATLABFunction_g5,
      &Human_in_Loop_DW.sf_MATLABFunction_g5);
  }

  /* StateSpace: '<S53>/high_pass' */
  Human_in_Loop_B.high_pass_g = 0.0;
  Human_in_Loop_B.high_pass_g += Human_in_Loop_P.high_pass_C_a *
    Human_in_Loop_X.high_pass_CSTATE_my[1];

  /* Abs: '<S53>/Abs' */
  Human_in_Loop_B.Abs_f = fabs(Human_in_Loop_B.high_pass_g);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* MATLAB Function: '<S53>/MATLAB Function' */
    Human_in_Loop_MATLABFunction(Human_in_Loop_B.Abs_f,
      &Human_in_Loop_B.sf_MATLABFunction_p,
      &Human_in_Loop_DW.sf_MATLABFunction_p);

    /* MATLAB Function: '<S29>/Mux3' */
    Human_in_Loop_Mux1(Human_in_Loop_B.sf_MATLABFunction_g.y,
                       Human_in_Loop_B.sf_MATLABFunction_m.y,
                       Human_in_Loop_B.sf_MATLABFunction_n.y,
                       Human_in_Loop_B.sf_MATLABFunction_l.y,
                       Human_in_Loop_B.sf_MATLABFunction_g5.y,
                       Human_in_Loop_B.sf_MATLABFunction_p.y,
                       &Human_in_Loop_B.sf_Mux3);
  }

  /* StateSpace: '<S50>/high_pass' */
  Human_in_Loop_B.high_pass_h = 0.0;
  Human_in_Loop_B.high_pass_h += Human_in_Loop_P.high_pass_C_ja[0] *
    Human_in_Loop_X.high_pass_CSTATE_p[0];
  Human_in_Loop_B.high_pass_h += Human_in_Loop_P.high_pass_C_ja[1] *
    Human_in_Loop_X.high_pass_CSTATE_p[1];
  Human_in_Loop_B.high_pass_h += Human_in_Loop_P.high_pass_D *
    Human_in_Loop_B.Gain_k;

  /* Abs: '<S50>/Abs' */
  Human_in_Loop_B.Abs_h = fabs(Human_in_Loop_B.high_pass_h);

  /* StateSpace: '<S51>/high_pass' */
  Human_in_Loop_B.high_pass_oe = 0.0;
  Human_in_Loop_B.high_pass_oe += Human_in_Loop_P.high_pass_C_d[0] *
    Human_in_Loop_X.high_pass_CSTATE_n[0];
  Human_in_Loop_B.high_pass_oe += Human_in_Loop_P.high_pass_C_d[1] *
    Human_in_Loop_X.high_pass_CSTATE_n[1];
  Human_in_Loop_B.high_pass_oe += Human_in_Loop_P.high_pass_D_m *
    Human_in_Loop_B.Gain1_p;

  /* Abs: '<S51>/Abs' */
  Human_in_Loop_B.Abs_d = fabs(Human_in_Loop_B.high_pass_oe);

  /* StateSpace: '<S54>/high_pass' */
  Human_in_Loop_B.high_pass_ba = 0.0;
  Human_in_Loop_B.high_pass_ba += Human_in_Loop_P.high_pass_C_g[0] *
    Human_in_Loop_X.high_pass_CSTATE_pw[0];
  Human_in_Loop_B.high_pass_ba += Human_in_Loop_P.high_pass_C_g[1] *
    Human_in_Loop_X.high_pass_CSTATE_pw[1];
  Human_in_Loop_B.high_pass_ba += Human_in_Loop_P.high_pass_D_j *
    Human_in_Loop_B.Gain2_cz;

  /* Abs: '<S54>/Abs' */
  Human_in_Loop_B.Abs_j = fabs(Human_in_Loop_B.high_pass_ba);

  /* StateSpace: '<S55>/high_pass' */
  Human_in_Loop_B.high_pass_hl = 0.0;
  Human_in_Loop_B.high_pass_hl += Human_in_Loop_P.high_pass_C_h[0] *
    Human_in_Loop_X.high_pass_CSTATE_k[0];
  Human_in_Loop_B.high_pass_hl += Human_in_Loop_P.high_pass_C_h[1] *
    Human_in_Loop_X.high_pass_CSTATE_k[1];
  Human_in_Loop_B.high_pass_hl += Human_in_Loop_P.high_pass_D_a *
    Human_in_Loop_B.Gain3_f;

  /* Abs: '<S55>/Abs' */
  Human_in_Loop_B.Abs_a = fabs(Human_in_Loop_B.high_pass_hl);

  /* StateSpace: '<S60>/high_pass' */
  Human_in_Loop_B.high_pass_a = 0.0;
  Human_in_Loop_B.high_pass_a += Human_in_Loop_P.high_pass_C_dd[0] *
    Human_in_Loop_X.high_pass_CSTATE_l[0];
  Human_in_Loop_B.high_pass_a += Human_in_Loop_P.high_pass_C_dd[1] *
    Human_in_Loop_X.high_pass_CSTATE_l[1];
  Human_in_Loop_B.high_pass_a += Human_in_Loop_P.high_pass_D_n *
    Human_in_Loop_B.Gain5;

  /* Abs: '<S60>/Abs' */
  Human_in_Loop_B.Abs_dj = fabs(Human_in_Loop_B.high_pass_a);

  /* StateSpace: '<S61>/high_pass' */
  Human_in_Loop_B.high_pass_e = 0.0;
  Human_in_Loop_B.high_pass_e += Human_in_Loop_P.high_pass_C_dc[0] *
    Human_in_Loop_X.high_pass_CSTATE_h[0];
  Human_in_Loop_B.high_pass_e += Human_in_Loop_P.high_pass_C_dc[1] *
    Human_in_Loop_X.high_pass_CSTATE_h[1];
  Human_in_Loop_B.high_pass_e += Human_in_Loop_P.high_pass_D_mh *
    Human_in_Loop_B.Gain4_h;

  /* Abs: '<S61>/Abs' */
  Human_in_Loop_B.Abs_fy = fabs(Human_in_Loop_B.high_pass_e);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* MATLAB Function: '<S28>/Timer' */
    /* MATLAB Function 'Sensor Data/Cosmed/Timer': '<S35>:1' */
    /* '<S35>:1:7' */
    /* '<S35>:1:8' */
    if (Human_in_Loop_P.Timer_BT_METAB_EST) {
      /* '<S35>:1:16' */
      if (Human_in_Loop_DW.time > 120.0) {
        /* '<S35>:1:17' */
        /* '<S35>:1:18' */
        Human_in_Loop_DW.time = 0.0;

        /* '<S35>:1:19' */
        Human_in_Loop_DW.count++;

        /* '<S35>:1:20' */
        i = 1;
      } else {
        /* '<S35>:1:22' */
        Human_in_Loop_DW.time += 0.0002;

        /* '<S35>:1:23' */
        i = 0;
      }
    } else {
      /* '<S35>:1:26' */
      i = 0;

      /* '<S35>:1:27' */
      Human_in_Loop_DW.count = 0.0;

      /* '<S35>:1:28' */
      Human_in_Loop_DW.time += 0.0002;
    }

    /* '<S35>:1:31' */
    Human_in_Loop_B.Trigger = i;

    /* '<S35>:1:32' */
    Human_in_Loop_B.Time = Human_in_Loop_DW.time;

    /* '<S35>:1:33' */
    Human_in_Loop_B.Count = Human_in_Loop_DW.count;

    /* End of MATLAB Function: '<S28>/Timer' */

    /* RateTransition: '<S28>/RT1' */
    switch (Human_in_Loop_DW.RT1_read_buf_p) {
     case 0:
      Human_in_Loop_DW.RT1_write_buf_f = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT1_write_buf_f = 0;
      break;

     default:
      Human_in_Loop_DW.RT1_write_buf_f = (int8_T)
        (Human_in_Loop_DW.RT1_last_buf_wr_j == 0);
      break;
    }

    if (Human_in_Loop_DW.RT1_write_buf_f != 0) {
      Human_in_Loop_DW.RT1_Buffer1_f = Human_in_Loop_B.Time;
    } else {
      Human_in_Loop_DW.RT1_Buffer0_e = Human_in_Loop_B.Time;
    }

    Human_in_Loop_DW.RT1_last_buf_wr_j = Human_in_Loop_DW.RT1_write_buf_f;
    Human_in_Loop_DW.RT1_write_buf_f = -1;

    /* End of RateTransition: '<S28>/RT1' */

    /* RateTransition: '<S28>/RT2' */
    switch (Human_in_Loop_DW.RT2_read_buf_a) {
     case 0:
      Human_in_Loop_DW.RT2_write_buf_g = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT2_write_buf_g = 0;
      break;

     default:
      Human_in_Loop_DW.RT2_write_buf_g = (int8_T)
        (Human_in_Loop_DW.RT2_last_buf_wr_a == 0);
      break;
    }

    if (Human_in_Loop_DW.RT2_write_buf_g != 0) {
      Human_in_Loop_DW.RT2_Buffer1_o = Human_in_Loop_B.Count;
    } else {
      Human_in_Loop_DW.RT2_Buffer0_i = Human_in_Loop_B.Count;
    }

    Human_in_Loop_DW.RT2_last_buf_wr_a = Human_in_Loop_DW.RT2_write_buf_g;
    Human_in_Loop_DW.RT2_write_buf_g = -1;

    /* End of RateTransition: '<S28>/RT2' */

    /* S-Function (rti_commonblock): '<S36>/S-Function1' */

    /* This comment workarounds a code generation problem */

    /* End of Outputs for S-Function (rti_commonblock): '<S36>/S-Function1' */

    /* S-Function (rti_commonblock): '<S37>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* S-Function (rti_commonblock): '<S90>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* dSPACE RTICAN RX Message Block: "RX Message" Id:1 */
    {
      UInt32 *CAN_Msg;

      /* ... Read status and timestamp info (previous message) */
      can_tp1_msg_read_from_mem(can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X1]);

      /* Convert timestamp */
      if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X1]->processed) {
        can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X1]->timestamp =
          rtk_dsts_time_to_simtime_convert
          (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X1]->timestamp);
      }

      /* Messages with timestamp zero have been received in pause/stop state
         and must not be handled.
       */
      if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X1]->timestamp > 0.0) {
        if (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X1]->processed) {
          Human_in_Loop_B.SFunction1_o5 = (real_T)
            can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X1]->processed;
          Human_in_Loop_B.SFunction1_o6 = (real_T)
            can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X1]->timestamp;
          CAN_Msg = can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X1]->data;

          /* ... Decode CAN Message */
          {
            rtican_Signal_t CAN_Sgn;

            /* ...... "Angle1" (0|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[0];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[1];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            Human_in_Loop_B.SFunction1_o1_d = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "Angle2" (16|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            Human_in_Loop_B.SFunction1_o2_d = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "Gyro1" (32|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[4];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[5];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            Human_in_Loop_B.SFunction1_o3 = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "Gyro2" (48|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[6];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[7];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            Human_in_Loop_B.SFunction1_o4 = ((real_T) CAN_Sgn.SignedSgn);
          }
        }
      }

      if (!can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X1]->processed) {
        /* ... set RX status to 0 because no new message has arrived */
        Human_in_Loop_B.SFunction1_o5 = 0.0;
      }
    }

    /* DataTypeConversion: '<S93>/Data Type Conversion' */
    Human_in_Loop_B.DataTypeConversion = Human_in_Loop_B.SFunction1_o3;

    /* Gain: '<S93>/Gain' */
    Human_in_Loop_B.Gain_h = Human_in_Loop_P.Gain_Gain_cd *
      Human_in_Loop_B.DataTypeConversion;

    /* DataTypeConversion: '<S94>/Data Type Conversion' */
    Human_in_Loop_B.DataTypeConversion_f = Human_in_Loop_B.SFunction1_o4;

    /* Gain: '<S94>/Gain' */
    Human_in_Loop_B.Gain_m = Human_in_Loop_P.Gain_Gain_n *
      Human_in_Loop_B.DataTypeConversion_f;

    /* Sum: '<S32>/Sum' */
    Human_in_Loop_B.Sum = Human_in_Loop_B.Gain_h - Human_in_Loop_B.Gain_m;

    /* Gain: '<S32>/Gain' */
    Human_in_Loop_B.Gain_e = Human_in_Loop_P.Gain_Gain_m * Human_in_Loop_B.Sum;

    /* DataTypeConversion: '<S91>/Data Type Conversion' */
    Human_in_Loop_B.DataTypeConversion_fh = Human_in_Loop_B.SFunction1_o1_d;

    /* Gain: '<S91>/Gain' */
    Human_in_Loop_B.Gain_kx = Human_in_Loop_P.Gain_Gain_e *
      Human_in_Loop_B.DataTypeConversion_fh;

    /* DataTypeConversion: '<S92>/Data Type Conversion' */
    Human_in_Loop_B.DataTypeConversion_a = Human_in_Loop_B.SFunction1_o2_d;

    /* Gain: '<S92>/Gain' */
    Human_in_Loop_B.Gain_gm = Human_in_Loop_P.Gain_Gain_a *
      Human_in_Loop_B.DataTypeConversion_a;

    /* MATLAB Function: '<S32>/MATLAB Function' incorporates:
     *  Constant: '<S32>/Calibration'
     */
    /* MATLAB Function 'Sensor Data/IMU/MATLAB Function': '<S87>:1' */
    if (Human_in_Loop_P.Calibration_Value == 1.0) {
      /* '<S87>:1:11' */
      /* '<S87>:1:12' */
      Human_in_Loop_DW.angle1_zero = 0.9 * Human_in_Loop_DW.angle1_zero + 0.1 *
        Human_in_Loop_B.Gain_kx;

      /* '<S87>:1:13' */
      Human_in_Loop_DW.angle2_zero = 0.9 * Human_in_Loop_DW.angle2_zero + 0.1 *
        Human_in_Loop_B.Gain_gm;
    }

    /* '<S87>:1:16' */
    Human_in_Loop_B.Angle_Ankle = (Human_in_Loop_B.Gain_kx -
      Human_in_Loop_DW.angle1_zero) - (Human_in_Loop_B.Gain_gm -
      Human_in_Loop_DW.angle2_zero);

    /* '<S87>:1:18' */
    Human_in_Loop_B.Angle_Ankle = -Human_in_Loop_B.Angle_Ankle;

    /* End of MATLAB Function: '<S32>/MATLAB Function' */
    /* S-Function (rti_commonblock): '<S83>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* S-Function (rti_commonblock): '<S84>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* S-Function (rti_commonblock): '<S88>/S-Function1' */
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
          Human_in_Loop_B.SFunction1_o2_e = (real_T)
            can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64]->processed;
          Human_in_Loop_B.SFunction1_o3_j = (real_T)
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
            Human_in_Loop_B.SFunction1_o1_i = ((real_T) CAN_Sgn.IeeeVal32);
          }
        }
      }

      if (!can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64]->processed) {
        /* ... set RX status to 0 because no new message has arrived */
        Human_in_Loop_B.SFunction1_o2_e = 0.0;
      }
    }
  }

  /* Sin: '<S85>/Sine Wave' */
  Human_in_Loop_B.SineWave = sin(Human_in_Loop_P.SineWave_Freq *
    Human_in_Loop_M->Timing.t[0] + Human_in_Loop_P.SineWave_Phase) *
    Human_in_Loop_P.SineWave_Amp + Human_in_Loop_P.SineWave_Bias;
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S89>/S-Function1' */
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
        Human_in_Loop_B.SFunction1_o1_g = (real_T)
          can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->processed;
        Human_in_Loop_B.SFunction1_o2_g = (real_T)
          can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->timestamp;
        Human_in_Loop_B.SFunction1_o3_c = (real_T)
          can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->deltatime;
        Human_in_Loop_B.SFunction1_o4_h = (real_T)
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
    /* Update for DiscreteFilter: '<S33>/0.4low2' */
    Human_in_Loop_DW.u4low2_states[2] = Human_in_Loop_DW.u4low2_states[1];
    Human_in_Loop_DW.u4low2_states[1] = Human_in_Loop_DW.u4low2_states[0];
    Human_in_Loop_DW.u4low2_states[0] = Human_in_Loop_DW.u4low2_tmp;

    /* Update for UnitDelay: '<S95>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE = Human_in_Loop_B.x2k;

    /* Update for UnitDelay: '<S95>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE = Human_in_Loop_B.x1k;

    /* Update for UnitDelay: '<S95>/Unit Delay2' */
    Human_in_Loop_DW.UnitDelay2_DSTATE = Human_in_Loop_B.torque;

    /* Update for UnitDelay: '<S69>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE_h = Human_in_Loop_B.x2k_i;

    /* Update for UnitDelay: '<S69>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_e = Human_in_Loop_B.x1k_i;

    /* Update for UnitDelay: '<S69>/Unit Delay2' */
    Human_in_Loop_DW.UnitDelay2_DSTATE_a = Human_in_Loop_B.Gain3_m;

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

    /* Update for UnitDelay: '<S72>/UD' */
    Human_in_Loop_DW.UD_DSTATE = Human_in_Loop_B.TSamp;

    /* Update for DiscreteFilter: '<S30>/0.4low1' */
    Human_in_Loop_DW.u4low1_states[2] = Human_in_Loop_DW.u4low1_states[1];
    Human_in_Loop_DW.u4low1_states[1] = Human_in_Loop_DW.u4low1_states[0];
    Human_in_Loop_DW.u4low1_states[0] = Human_in_Loop_DW.u4low1_tmp;

    /* Update for UnitDelay: '<S68>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_j = Human_in_Loop_B.Gain3_m;

    /* Update for UnitDelay: '<S68>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE_a = Human_in_Loop_B.Add2_i;
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

  /* Derivatives for StateSpace: '<S50>/low_pass' */
  _rtXdot->low_pass_CSTATE[0] = 0.0;
  _rtXdot->low_pass_CSTATE[1] = 0.0;
  _rtXdot->low_pass_CSTATE[0] += Human_in_Loop_P.low_pass_A[0] *
    Human_in_Loop_X.low_pass_CSTATE[0];
  _rtXdot->low_pass_CSTATE[0] += Human_in_Loop_P.low_pass_A[1] *
    Human_in_Loop_X.low_pass_CSTATE[1];
  _rtXdot->low_pass_CSTATE[1] += Human_in_Loop_P.low_pass_A[2] *
    Human_in_Loop_X.low_pass_CSTATE[0];
  _rtXdot->low_pass_CSTATE[0] += Human_in_Loop_P.low_pass_B *
    Human_in_Loop_B.Abs_h;

  /* Derivatives for StateSpace: '<S51>/low_pass' */
  _rtXdot->low_pass_CSTATE_o[0] = 0.0;
  _rtXdot->low_pass_CSTATE_o[1] = 0.0;
  _rtXdot->low_pass_CSTATE_o[0] += Human_in_Loop_P.low_pass_A_k[0] *
    Human_in_Loop_X.low_pass_CSTATE_o[0];
  _rtXdot->low_pass_CSTATE_o[0] += Human_in_Loop_P.low_pass_A_k[1] *
    Human_in_Loop_X.low_pass_CSTATE_o[1];
  _rtXdot->low_pass_CSTATE_o[1] += Human_in_Loop_P.low_pass_A_k[2] *
    Human_in_Loop_X.low_pass_CSTATE_o[0];
  _rtXdot->low_pass_CSTATE_o[0] += Human_in_Loop_P.low_pass_B_o *
    Human_in_Loop_B.Abs_d;

  /* Derivatives for StateSpace: '<S54>/low_pass' */
  _rtXdot->low_pass_CSTATE_n[0] = 0.0;
  _rtXdot->low_pass_CSTATE_n[1] = 0.0;
  _rtXdot->low_pass_CSTATE_n[0] += Human_in_Loop_P.low_pass_A_f[0] *
    Human_in_Loop_X.low_pass_CSTATE_n[0];
  _rtXdot->low_pass_CSTATE_n[0] += Human_in_Loop_P.low_pass_A_f[1] *
    Human_in_Loop_X.low_pass_CSTATE_n[1];
  _rtXdot->low_pass_CSTATE_n[1] += Human_in_Loop_P.low_pass_A_f[2] *
    Human_in_Loop_X.low_pass_CSTATE_n[0];
  _rtXdot->low_pass_CSTATE_n[0] += Human_in_Loop_P.low_pass_B_j *
    Human_in_Loop_B.Abs_j;

  /* Derivatives for StateSpace: '<S55>/low_pass' */
  _rtXdot->low_pass_CSTATE_p[0] = 0.0;
  _rtXdot->low_pass_CSTATE_p[1] = 0.0;
  _rtXdot->low_pass_CSTATE_p[0] += Human_in_Loop_P.low_pass_A_a[0] *
    Human_in_Loop_X.low_pass_CSTATE_p[0];
  _rtXdot->low_pass_CSTATE_p[0] += Human_in_Loop_P.low_pass_A_a[1] *
    Human_in_Loop_X.low_pass_CSTATE_p[1];
  _rtXdot->low_pass_CSTATE_p[1] += Human_in_Loop_P.low_pass_A_a[2] *
    Human_in_Loop_X.low_pass_CSTATE_p[0];
  _rtXdot->low_pass_CSTATE_p[0] += Human_in_Loop_P.low_pass_B_m *
    Human_in_Loop_B.Abs_a;

  /* Derivatives for StateSpace: '<S61>/low_pass' */
  _rtXdot->low_pass_CSTATE_on[0] = 0.0;
  _rtXdot->low_pass_CSTATE_on[1] = 0.0;
  _rtXdot->low_pass_CSTATE_on[0] += Human_in_Loop_P.low_pass_A_b[0] *
    Human_in_Loop_X.low_pass_CSTATE_on[0];
  _rtXdot->low_pass_CSTATE_on[0] += Human_in_Loop_P.low_pass_A_b[1] *
    Human_in_Loop_X.low_pass_CSTATE_on[1];
  _rtXdot->low_pass_CSTATE_on[1] += Human_in_Loop_P.low_pass_A_b[2] *
    Human_in_Loop_X.low_pass_CSTATE_on[0];
  _rtXdot->low_pass_CSTATE_on[0] += Human_in_Loop_P.low_pass_B_g *
    Human_in_Loop_B.Abs_fy;

  /* Derivatives for StateSpace: '<S60>/low_pass' */
  _rtXdot->low_pass_CSTATE_f[0] = 0.0;
  _rtXdot->low_pass_CSTATE_f[1] = 0.0;
  _rtXdot->low_pass_CSTATE_f[0] += Human_in_Loop_P.low_pass_A_i[0] *
    Human_in_Loop_X.low_pass_CSTATE_f[0];
  _rtXdot->low_pass_CSTATE_f[0] += Human_in_Loop_P.low_pass_A_i[1] *
    Human_in_Loop_X.low_pass_CSTATE_f[1];
  _rtXdot->low_pass_CSTATE_f[1] += Human_in_Loop_P.low_pass_A_i[2] *
    Human_in_Loop_X.low_pass_CSTATE_f[0];
  _rtXdot->low_pass_CSTATE_f[0] += Human_in_Loop_P.low_pass_B_os *
    Human_in_Loop_B.Abs_dj;

  /* Derivatives for StateSpace: '<S56>/high_pass' */
  _rtXdot->high_pass_CSTATE[0] = 0.0;
  _rtXdot->high_pass_CSTATE[1] = 0.0;
  _rtXdot->high_pass_CSTATE[2] = 0.0;
  _rtXdot->high_pass_CSTATE[3] = 0.0;
  _rtXdot->high_pass_CSTATE[0] += Human_in_Loop_P.high_pass_A[0] *
    Human_in_Loop_X.high_pass_CSTATE[0];
  _rtXdot->high_pass_CSTATE[0] += Human_in_Loop_P.high_pass_A[1] *
    Human_in_Loop_X.high_pass_CSTATE[1];
  _rtXdot->high_pass_CSTATE[0] += Human_in_Loop_P.high_pass_A[2] *
    Human_in_Loop_X.high_pass_CSTATE[2];
  _rtXdot->high_pass_CSTATE[1] += Human_in_Loop_P.high_pass_A[3] *
    Human_in_Loop_X.high_pass_CSTATE[0];
  _rtXdot->high_pass_CSTATE[1] += Human_in_Loop_P.high_pass_A[4] *
    Human_in_Loop_X.high_pass_CSTATE[3];
  _rtXdot->high_pass_CSTATE[2] += Human_in_Loop_P.high_pass_A[5] *
    Human_in_Loop_X.high_pass_CSTATE[0];
  _rtXdot->high_pass_CSTATE[3] += Human_in_Loop_P.high_pass_A[6] *
    Human_in_Loop_X.high_pass_CSTATE[1];
  _rtXdot->high_pass_CSTATE[0] += Human_in_Loop_P.high_pass_B *
    Human_in_Loop_B.Gain_k;

  /* Derivatives for StateSpace: '<S57>/high_pass' */
  _rtXdot->high_pass_CSTATE_b[0] = 0.0;
  _rtXdot->high_pass_CSTATE_b[1] = 0.0;
  _rtXdot->high_pass_CSTATE_b[2] = 0.0;
  _rtXdot->high_pass_CSTATE_b[3] = 0.0;
  _rtXdot->high_pass_CSTATE_b[0] += Human_in_Loop_P.high_pass_A_p[0] *
    Human_in_Loop_X.high_pass_CSTATE_b[0];
  _rtXdot->high_pass_CSTATE_b[0] += Human_in_Loop_P.high_pass_A_p[1] *
    Human_in_Loop_X.high_pass_CSTATE_b[1];
  _rtXdot->high_pass_CSTATE_b[0] += Human_in_Loop_P.high_pass_A_p[2] *
    Human_in_Loop_X.high_pass_CSTATE_b[2];
  _rtXdot->high_pass_CSTATE_b[1] += Human_in_Loop_P.high_pass_A_p[3] *
    Human_in_Loop_X.high_pass_CSTATE_b[0];
  _rtXdot->high_pass_CSTATE_b[1] += Human_in_Loop_P.high_pass_A_p[4] *
    Human_in_Loop_X.high_pass_CSTATE_b[3];
  _rtXdot->high_pass_CSTATE_b[2] += Human_in_Loop_P.high_pass_A_p[5] *
    Human_in_Loop_X.high_pass_CSTATE_b[0];
  _rtXdot->high_pass_CSTATE_b[3] += Human_in_Loop_P.high_pass_A_p[6] *
    Human_in_Loop_X.high_pass_CSTATE_b[1];
  _rtXdot->high_pass_CSTATE_b[0] += Human_in_Loop_P.high_pass_B_g *
    Human_in_Loop_B.Gain1_p;

  /* Derivatives for StateSpace: '<S58>/high_pass' */
  _rtXdot->high_pass_CSTATE_c[0] = 0.0;
  _rtXdot->high_pass_CSTATE_c[1] = 0.0;
  _rtXdot->high_pass_CSTATE_c[2] = 0.0;
  _rtXdot->high_pass_CSTATE_c[3] = 0.0;
  _rtXdot->high_pass_CSTATE_c[0] += Human_in_Loop_P.high_pass_A_e[0] *
    Human_in_Loop_X.high_pass_CSTATE_c[0];
  _rtXdot->high_pass_CSTATE_c[0] += Human_in_Loop_P.high_pass_A_e[1] *
    Human_in_Loop_X.high_pass_CSTATE_c[1];
  _rtXdot->high_pass_CSTATE_c[0] += Human_in_Loop_P.high_pass_A_e[2] *
    Human_in_Loop_X.high_pass_CSTATE_c[2];
  _rtXdot->high_pass_CSTATE_c[1] += Human_in_Loop_P.high_pass_A_e[3] *
    Human_in_Loop_X.high_pass_CSTATE_c[0];
  _rtXdot->high_pass_CSTATE_c[1] += Human_in_Loop_P.high_pass_A_e[4] *
    Human_in_Loop_X.high_pass_CSTATE_c[3];
  _rtXdot->high_pass_CSTATE_c[2] += Human_in_Loop_P.high_pass_A_e[5] *
    Human_in_Loop_X.high_pass_CSTATE_c[0];
  _rtXdot->high_pass_CSTATE_c[3] += Human_in_Loop_P.high_pass_A_e[6] *
    Human_in_Loop_X.high_pass_CSTATE_c[1];
  _rtXdot->high_pass_CSTATE_c[0] += Human_in_Loop_P.high_pass_B_i *
    Human_in_Loop_B.Gain2_cz;

  /* Derivatives for StateSpace: '<S59>/high_pass' */
  _rtXdot->high_pass_CSTATE_m[0] = 0.0;
  _rtXdot->high_pass_CSTATE_m[1] = 0.0;
  _rtXdot->high_pass_CSTATE_m[2] = 0.0;
  _rtXdot->high_pass_CSTATE_m[3] = 0.0;
  _rtXdot->high_pass_CSTATE_m[0] += Human_in_Loop_P.high_pass_A_o[0] *
    Human_in_Loop_X.high_pass_CSTATE_m[0];
  _rtXdot->high_pass_CSTATE_m[0] += Human_in_Loop_P.high_pass_A_o[1] *
    Human_in_Loop_X.high_pass_CSTATE_m[1];
  _rtXdot->high_pass_CSTATE_m[0] += Human_in_Loop_P.high_pass_A_o[2] *
    Human_in_Loop_X.high_pass_CSTATE_m[2];
  _rtXdot->high_pass_CSTATE_m[1] += Human_in_Loop_P.high_pass_A_o[3] *
    Human_in_Loop_X.high_pass_CSTATE_m[0];
  _rtXdot->high_pass_CSTATE_m[1] += Human_in_Loop_P.high_pass_A_o[4] *
    Human_in_Loop_X.high_pass_CSTATE_m[3];
  _rtXdot->high_pass_CSTATE_m[2] += Human_in_Loop_P.high_pass_A_o[5] *
    Human_in_Loop_X.high_pass_CSTATE_m[0];
  _rtXdot->high_pass_CSTATE_m[3] += Human_in_Loop_P.high_pass_A_o[6] *
    Human_in_Loop_X.high_pass_CSTATE_m[1];
  _rtXdot->high_pass_CSTATE_m[0] += Human_in_Loop_P.high_pass_B_c *
    Human_in_Loop_B.Gain3_f;

  /* Derivatives for StateSpace: '<S52>/high_pass' */
  _rtXdot->high_pass_CSTATE_a[0] = 0.0;
  _rtXdot->high_pass_CSTATE_a[1] = 0.0;
  _rtXdot->high_pass_CSTATE_a[2] = 0.0;
  _rtXdot->high_pass_CSTATE_a[3] = 0.0;
  _rtXdot->high_pass_CSTATE_a[0] += Human_in_Loop_P.high_pass_A_h[0] *
    Human_in_Loop_X.high_pass_CSTATE_a[0];
  _rtXdot->high_pass_CSTATE_a[0] += Human_in_Loop_P.high_pass_A_h[1] *
    Human_in_Loop_X.high_pass_CSTATE_a[1];
  _rtXdot->high_pass_CSTATE_a[0] += Human_in_Loop_P.high_pass_A_h[2] *
    Human_in_Loop_X.high_pass_CSTATE_a[2];
  _rtXdot->high_pass_CSTATE_a[1] += Human_in_Loop_P.high_pass_A_h[3] *
    Human_in_Loop_X.high_pass_CSTATE_a[0];
  _rtXdot->high_pass_CSTATE_a[1] += Human_in_Loop_P.high_pass_A_h[4] *
    Human_in_Loop_X.high_pass_CSTATE_a[3];
  _rtXdot->high_pass_CSTATE_a[2] += Human_in_Loop_P.high_pass_A_h[5] *
    Human_in_Loop_X.high_pass_CSTATE_a[0];
  _rtXdot->high_pass_CSTATE_a[3] += Human_in_Loop_P.high_pass_A_h[6] *
    Human_in_Loop_X.high_pass_CSTATE_a[1];
  _rtXdot->high_pass_CSTATE_a[0] += Human_in_Loop_P.high_pass_B_d *
    Human_in_Loop_B.Gain4_h;

  /* Derivatives for StateSpace: '<S53>/high_pass' */
  _rtXdot->high_pass_CSTATE_my[0] = 0.0;
  _rtXdot->high_pass_CSTATE_my[1] = 0.0;
  _rtXdot->high_pass_CSTATE_my[2] = 0.0;
  _rtXdot->high_pass_CSTATE_my[3] = 0.0;
  _rtXdot->high_pass_CSTATE_my[0] += Human_in_Loop_P.high_pass_A_i[0] *
    Human_in_Loop_X.high_pass_CSTATE_my[0];
  _rtXdot->high_pass_CSTATE_my[0] += Human_in_Loop_P.high_pass_A_i[1] *
    Human_in_Loop_X.high_pass_CSTATE_my[1];
  _rtXdot->high_pass_CSTATE_my[0] += Human_in_Loop_P.high_pass_A_i[2] *
    Human_in_Loop_X.high_pass_CSTATE_my[2];
  _rtXdot->high_pass_CSTATE_my[1] += Human_in_Loop_P.high_pass_A_i[3] *
    Human_in_Loop_X.high_pass_CSTATE_my[0];
  _rtXdot->high_pass_CSTATE_my[1] += Human_in_Loop_P.high_pass_A_i[4] *
    Human_in_Loop_X.high_pass_CSTATE_my[3];
  _rtXdot->high_pass_CSTATE_my[2] += Human_in_Loop_P.high_pass_A_i[5] *
    Human_in_Loop_X.high_pass_CSTATE_my[0];
  _rtXdot->high_pass_CSTATE_my[3] += Human_in_Loop_P.high_pass_A_i[6] *
    Human_in_Loop_X.high_pass_CSTATE_my[1];
  _rtXdot->high_pass_CSTATE_my[0] += Human_in_Loop_P.high_pass_B_h *
    Human_in_Loop_B.Gain5;

  /* Derivatives for StateSpace: '<S50>/high_pass' */
  _rtXdot->high_pass_CSTATE_p[0] = 0.0;
  _rtXdot->high_pass_CSTATE_p[1] = 0.0;
  _rtXdot->high_pass_CSTATE_p[0] += Human_in_Loop_P.high_pass_A_d[0] *
    Human_in_Loop_X.high_pass_CSTATE_p[1];
  _rtXdot->high_pass_CSTATE_p[1] += Human_in_Loop_P.high_pass_A_d[1] *
    Human_in_Loop_X.high_pass_CSTATE_p[0];
  _rtXdot->high_pass_CSTATE_p[1] += Human_in_Loop_P.high_pass_A_d[2] *
    Human_in_Loop_X.high_pass_CSTATE_p[1];
  _rtXdot->high_pass_CSTATE_p[1] += Human_in_Loop_P.high_pass_B_b *
    Human_in_Loop_B.Gain_k;

  /* Derivatives for StateSpace: '<S51>/high_pass' */
  _rtXdot->high_pass_CSTATE_n[0] = 0.0;
  _rtXdot->high_pass_CSTATE_n[1] = 0.0;
  _rtXdot->high_pass_CSTATE_n[0] += Human_in_Loop_P.high_pass_A_b[0] *
    Human_in_Loop_X.high_pass_CSTATE_n[1];
  _rtXdot->high_pass_CSTATE_n[1] += Human_in_Loop_P.high_pass_A_b[1] *
    Human_in_Loop_X.high_pass_CSTATE_n[0];
  _rtXdot->high_pass_CSTATE_n[1] += Human_in_Loop_P.high_pass_A_b[2] *
    Human_in_Loop_X.high_pass_CSTATE_n[1];
  _rtXdot->high_pass_CSTATE_n[1] += Human_in_Loop_P.high_pass_B_j *
    Human_in_Loop_B.Gain1_p;

  /* Derivatives for StateSpace: '<S54>/high_pass' */
  _rtXdot->high_pass_CSTATE_pw[0] = 0.0;
  _rtXdot->high_pass_CSTATE_pw[1] = 0.0;
  _rtXdot->high_pass_CSTATE_pw[0] += Human_in_Loop_P.high_pass_A_j[0] *
    Human_in_Loop_X.high_pass_CSTATE_pw[1];
  _rtXdot->high_pass_CSTATE_pw[1] += Human_in_Loop_P.high_pass_A_j[1] *
    Human_in_Loop_X.high_pass_CSTATE_pw[0];
  _rtXdot->high_pass_CSTATE_pw[1] += Human_in_Loop_P.high_pass_A_j[2] *
    Human_in_Loop_X.high_pass_CSTATE_pw[1];
  _rtXdot->high_pass_CSTATE_pw[1] += Human_in_Loop_P.high_pass_B_e *
    Human_in_Loop_B.Gain2_cz;

  /* Derivatives for StateSpace: '<S55>/high_pass' */
  _rtXdot->high_pass_CSTATE_k[0] = 0.0;
  _rtXdot->high_pass_CSTATE_k[1] = 0.0;
  _rtXdot->high_pass_CSTATE_k[0] += Human_in_Loop_P.high_pass_A_l[0] *
    Human_in_Loop_X.high_pass_CSTATE_k[1];
  _rtXdot->high_pass_CSTATE_k[1] += Human_in_Loop_P.high_pass_A_l[1] *
    Human_in_Loop_X.high_pass_CSTATE_k[0];
  _rtXdot->high_pass_CSTATE_k[1] += Human_in_Loop_P.high_pass_A_l[2] *
    Human_in_Loop_X.high_pass_CSTATE_k[1];
  _rtXdot->high_pass_CSTATE_k[1] += Human_in_Loop_P.high_pass_B_o *
    Human_in_Loop_B.Gain3_f;

  /* Derivatives for StateSpace: '<S60>/high_pass' */
  _rtXdot->high_pass_CSTATE_l[0] = 0.0;
  _rtXdot->high_pass_CSTATE_l[1] = 0.0;
  _rtXdot->high_pass_CSTATE_l[0] += Human_in_Loop_P.high_pass_A_m[0] *
    Human_in_Loop_X.high_pass_CSTATE_l[1];
  _rtXdot->high_pass_CSTATE_l[1] += Human_in_Loop_P.high_pass_A_m[1] *
    Human_in_Loop_X.high_pass_CSTATE_l[0];
  _rtXdot->high_pass_CSTATE_l[1] += Human_in_Loop_P.high_pass_A_m[2] *
    Human_in_Loop_X.high_pass_CSTATE_l[1];
  _rtXdot->high_pass_CSTATE_l[1] += Human_in_Loop_P.high_pass_B_hu *
    Human_in_Loop_B.Gain5;

  /* Derivatives for StateSpace: '<S61>/high_pass' */
  _rtXdot->high_pass_CSTATE_h[0] = 0.0;
  _rtXdot->high_pass_CSTATE_h[1] = 0.0;
  _rtXdot->high_pass_CSTATE_h[0] += Human_in_Loop_P.high_pass_A_f[0] *
    Human_in_Loop_X.high_pass_CSTATE_h[1];
  _rtXdot->high_pass_CSTATE_h[1] += Human_in_Loop_P.high_pass_A_f[1] *
    Human_in_Loop_X.high_pass_CSTATE_h[0];
  _rtXdot->high_pass_CSTATE_h[1] += Human_in_Loop_P.high_pass_A_f[2] *
    Human_in_Loop_X.high_pass_CSTATE_h[1];
  _rtXdot->high_pass_CSTATE_h[1] += Human_in_Loop_P.high_pass_B_bm *
    Human_in_Loop_B.Gain4_h;
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
    for (i = 0; i < 5; i++) {
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

    /* Start for RateTransition: '<S28>/RT1' */
    Human_in_Loop_B.RT1_c = Human_in_Loop_P.RT1_InitialCondition_c;

    /* Start for RateTransition: '<S28>/RT2' */
    Human_in_Loop_B.RT2_g = Human_in_Loop_P.RT2_InitialCondition_i;

    /* Start for DataStoreMemory: '<Root>/Data Store Memory' */
    memcpy(&Human_in_Loop_DW.TorqueMem[0],
           &Human_in_Loop_P.DataStoreMemory_InitialValue_l[0], 4400U * sizeof
           (real_T));
  }

  Human_in_Loop_PrevZCX.MeanCalculate_Trig_ZCE = UNINITIALIZED_ZCSIG;

  {
    int32_T i;

    /* InitializeConditions for DiscreteFilter: '<S33>/0.4low2' */
    Human_in_Loop_DW.u4low2_states[0] = Human_in_Loop_P.u4low2_InitialStates;
    Human_in_Loop_DW.u4low2_states[1] = Human_in_Loop_P.u4low2_InitialStates;
    Human_in_Loop_DW.u4low2_states[2] = Human_in_Loop_P.u4low2_InitialStates;

    /* InitializeConditions for UnitDelay: '<S95>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE =
      Human_in_Loop_P.UnitDelay1_InitialCondition_i;

    /* InitializeConditions for UnitDelay: '<S95>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE =
      Human_in_Loop_P.UnitDelay_InitialCondition_d;

    /* InitializeConditions for UnitDelay: '<S95>/Unit Delay2' */
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

    /* InitializeConditions for UnitDelay: '<S69>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE_h =
      Human_in_Loop_P.UnitDelay1_InitialCondition_j;

    /* InitializeConditions for UnitDelay: '<S69>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_e =
      Human_in_Loop_P.UnitDelay_InitialCondition_j;

    /* InitializeConditions for UnitDelay: '<S69>/Unit Delay2' */
    Human_in_Loop_DW.UnitDelay2_DSTATE_a =
      Human_in_Loop_P.UnitDelay2_InitialCondition_h;

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
    for (i = 0; i < 5; i++) {
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

    /* InitializeConditions for StateSpace: '<S50>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE[0] =
      Human_in_Loop_P.low_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S51>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_o[0] =
      Human_in_Loop_P.low_pass_InitialCondition_e;

    /* InitializeConditions for StateSpace: '<S54>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_n[0] =
      Human_in_Loop_P.low_pass_InitialCondition_ew;

    /* InitializeConditions for StateSpace: '<S55>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_p[0] =
      Human_in_Loop_P.low_pass_InitialCondition_c;

    /* InitializeConditions for StateSpace: '<S61>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_on[0] =
      Human_in_Loop_P.low_pass_InitialCondition_l;

    /* InitializeConditions for StateSpace: '<S60>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_f[0] =
      Human_in_Loop_P.low_pass_InitialCondition_i;

    /* InitializeConditions for StateSpace: '<S50>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE[1] =
      Human_in_Loop_P.low_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S51>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_o[1] =
      Human_in_Loop_P.low_pass_InitialCondition_e;

    /* InitializeConditions for StateSpace: '<S54>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_n[1] =
      Human_in_Loop_P.low_pass_InitialCondition_ew;

    /* InitializeConditions for StateSpace: '<S55>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_p[1] =
      Human_in_Loop_P.low_pass_InitialCondition_c;

    /* InitializeConditions for StateSpace: '<S61>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_on[1] =
      Human_in_Loop_P.low_pass_InitialCondition_l;

    /* InitializeConditions for StateSpace: '<S60>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_f[1] =
      Human_in_Loop_P.low_pass_InitialCondition_i;
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

    /* InitializeConditions for UnitDelay: '<S72>/UD' */
    Human_in_Loop_DW.UD_DSTATE = Human_in_Loop_P.DiscreteDerivative_ICPrevScaled;

    /* InitializeConditions for DiscreteFilter: '<S30>/0.4low1' */
    Human_in_Loop_DW.u4low1_states[0] = Human_in_Loop_P.u4low1_InitialStates;
    Human_in_Loop_DW.u4low1_states[1] = Human_in_Loop_P.u4low1_InitialStates;
    Human_in_Loop_DW.u4low1_states[2] = Human_in_Loop_P.u4low1_InitialStates;

    /* InitializeConditions for UnitDelay: '<S68>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_j =
      Human_in_Loop_P.UnitDelay_InitialCondition_g;

    /* InitializeConditions for UnitDelay: '<S68>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE_a =
      Human_in_Loop_P.UnitDelay1_InitialCondition_i0;

    /* InitializeConditions for StateSpace: '<S56>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE[0] =
      Human_in_Loop_P.high_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S57>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_b[0] =
      Human_in_Loop_P.high_pass_InitialCondition_g;

    /* InitializeConditions for StateSpace: '<S58>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_c[0] =
      Human_in_Loop_P.high_pass_InitialCondition_n;

    /* InitializeConditions for StateSpace: '<S59>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_m[0] =
      Human_in_Loop_P.high_pass_InitialCondition_b;

    /* InitializeConditions for StateSpace: '<S52>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_a[0] =
      Human_in_Loop_P.high_pass_InitialCondition_bd;

    /* InitializeConditions for StateSpace: '<S53>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_my[0] =
      Human_in_Loop_P.high_pass_InitialCondition_m;

    /* InitializeConditions for StateSpace: '<S56>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE[1] =
      Human_in_Loop_P.high_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S57>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_b[1] =
      Human_in_Loop_P.high_pass_InitialCondition_g;

    /* InitializeConditions for StateSpace: '<S58>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_c[1] =
      Human_in_Loop_P.high_pass_InitialCondition_n;

    /* InitializeConditions for StateSpace: '<S59>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_m[1] =
      Human_in_Loop_P.high_pass_InitialCondition_b;

    /* InitializeConditions for StateSpace: '<S52>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_a[1] =
      Human_in_Loop_P.high_pass_InitialCondition_bd;

    /* InitializeConditions for StateSpace: '<S53>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_my[1] =
      Human_in_Loop_P.high_pass_InitialCondition_m;

    /* InitializeConditions for StateSpace: '<S56>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE[2] =
      Human_in_Loop_P.high_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S57>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_b[2] =
      Human_in_Loop_P.high_pass_InitialCondition_g;

    /* InitializeConditions for StateSpace: '<S58>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_c[2] =
      Human_in_Loop_P.high_pass_InitialCondition_n;

    /* InitializeConditions for StateSpace: '<S59>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_m[2] =
      Human_in_Loop_P.high_pass_InitialCondition_b;

    /* InitializeConditions for StateSpace: '<S52>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_a[2] =
      Human_in_Loop_P.high_pass_InitialCondition_bd;

    /* InitializeConditions for StateSpace: '<S53>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_my[2] =
      Human_in_Loop_P.high_pass_InitialCondition_m;

    /* InitializeConditions for StateSpace: '<S56>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE[3] =
      Human_in_Loop_P.high_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S57>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_b[3] =
      Human_in_Loop_P.high_pass_InitialCondition_g;

    /* InitializeConditions for StateSpace: '<S58>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_c[3] =
      Human_in_Loop_P.high_pass_InitialCondition_n;

    /* InitializeConditions for StateSpace: '<S59>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_m[3] =
      Human_in_Loop_P.high_pass_InitialCondition_b;

    /* InitializeConditions for StateSpace: '<S52>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_a[3] =
      Human_in_Loop_P.high_pass_InitialCondition_bd;

    /* InitializeConditions for StateSpace: '<S53>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_my[3] =
      Human_in_Loop_P.high_pass_InitialCondition_m;

    /* InitializeConditions for StateSpace: '<S50>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_p[0] =
      Human_in_Loop_P.high_pass_InitialCondition_p;

    /* InitializeConditions for StateSpace: '<S51>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_n[0] =
      Human_in_Loop_P.high_pass_InitialCondition_n4;

    /* InitializeConditions for StateSpace: '<S54>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_pw[0] =
      Human_in_Loop_P.high_pass_InitialCondition_m4;

    /* InitializeConditions for StateSpace: '<S55>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_k[0] =
      Human_in_Loop_P.high_pass_InitialCondition_gb;

    /* InitializeConditions for StateSpace: '<S60>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_l[0] =
      Human_in_Loop_P.high_pass_InitialCondition_c;

    /* InitializeConditions for StateSpace: '<S61>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_h[0] =
      Human_in_Loop_P.high_pass_InitialCondition_gw;

    /* InitializeConditions for StateSpace: '<S50>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_p[1] =
      Human_in_Loop_P.high_pass_InitialCondition_p;

    /* InitializeConditions for StateSpace: '<S51>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_n[1] =
      Human_in_Loop_P.high_pass_InitialCondition_n4;

    /* InitializeConditions for StateSpace: '<S54>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_pw[1] =
      Human_in_Loop_P.high_pass_InitialCondition_m4;

    /* InitializeConditions for StateSpace: '<S55>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_k[1] =
      Human_in_Loop_P.high_pass_InitialCondition_gb;

    /* InitializeConditions for StateSpace: '<S60>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_l[1] =
      Human_in_Loop_P.high_pass_InitialCondition_c;

    /* InitializeConditions for StateSpace: '<S61>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_h[1] =
      Human_in_Loop_P.high_pass_InitialCondition_gw;

    /* InitializeConditions for RateTransition: '<S28>/RT1' */
    Human_in_Loop_DW.RT1_Buffer0_e = Human_in_Loop_P.RT1_InitialCondition_c;
    Human_in_Loop_DW.RT1_write_buf_f = -1;
    Human_in_Loop_DW.RT1_read_buf_p = -1;

    /* InitializeConditions for RateTransition: '<S28>/RT2' */
    Human_in_Loop_DW.RT2_Buffer0_i = Human_in_Loop_P.RT2_InitialCondition_i;
    Human_in_Loop_DW.RT2_write_buf_g = -1;
    Human_in_Loop_DW.RT2_read_buf_a = -1;

    /* SystemInitialize for MATLAB Function: '<S33>/Data process' */
    Human_in_Loop_DW.torque_zero = 0.0;

    /* SystemInitialize for MATLAB Function: '<S30>/Data process' */
    Human_in_Loop_DW.angle_zero_f = 0.0;

    /* SystemInitialize for MATLAB Function: '<S30>/Data process1' */
    Human_in_Loop_DW.angle_zero = 0.0;

    /* SystemInitialize for MATLAB Function: '<S31>/FootSwitch Filter' */
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

    /* SystemInitialize for MATLAB Function: '<S33>/MATLAB Function' */
    for (i = 0; i < 15; i++) {
      Human_in_Loop_DW.data[i] = 1.0;
    }

    /* End of SystemInitialize for MATLAB Function: '<S33>/MATLAB Function' */

    /* SystemInitialize for MATLAB Function: '<S56>/MATLAB Function' */
    Human_in_Lo_MATLABFunction_Init(&Human_in_Loop_DW.sf_MATLABFunction_g);

    /* SystemInitialize for MATLAB Function: '<S57>/MATLAB Function' */
    Human_in_Lo_MATLABFunction_Init(&Human_in_Loop_DW.sf_MATLABFunction_m);

    /* SystemInitialize for MATLAB Function: '<S58>/MATLAB Function' */
    Human_in_Lo_MATLABFunction_Init(&Human_in_Loop_DW.sf_MATLABFunction_n);

    /* SystemInitialize for MATLAB Function: '<S59>/MATLAB Function' */
    Human_in_Lo_MATLABFunction_Init(&Human_in_Loop_DW.sf_MATLABFunction_l);

    /* SystemInitialize for MATLAB Function: '<S52>/MATLAB Function' */
    Human_in_Lo_MATLABFunction_Init(&Human_in_Loop_DW.sf_MATLABFunction_g5);

    /* SystemInitialize for MATLAB Function: '<S53>/MATLAB Function' */
    Human_in_Lo_MATLABFunction_Init(&Human_in_Loop_DW.sf_MATLABFunction_p);

    /* SystemInitialize for MATLAB Function: '<S28>/Timer' */
    Human_in_Loop_DW.count = 0.0;
    Human_in_Loop_DW.time = 0.0;

    /* SystemInitialize for S-Function (rti_commonblock): '<S36>/S-Function1' incorporates:
     *  SubSystem: '<S34>/Serial Decoding System'
     */
    Human_SerialDecodingSystem_Init();

    /* End of SystemInitialize for S-Function (rti_commonblock): '<S36>/S-Function1' */

    /* SystemInitialize for MATLAB Function: '<S32>/MATLAB Function' */
    Human_in_Loop_DW.angle1_zero = 0.0;
    Human_in_Loop_DW.angle2_zero = 0.0;
  }
}

/* Model terminate function */
void Human_in_Loop_terminate(void)
{
  /* Terminate for S-Function (rti_commonblock): '<S73>/S-Function1' */

  /* --- Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL1 --- */
  /* --- [RTIEMC, Encoder] - DIO class: 2 - Unit: 1 - Port: 1 - Channel: 1 --- */
  {
    /* Deactivates encoder interface functionality */
    DioCl2EncoderIn_stop(pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1);
  }

  /* Terminate for S-Function (rti_commonblock): '<S74>/S-Function1' */

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

  /* Terminate for S-Function (rti_commonblock): '<S75>/S-Function1' incorporates:
   *  Constant: '<S30>/VCC1'
   */

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup1 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 11 --- */

  /* disable digital output channel 11-11 on port 3 *
   * (set to high-impedance), when the simulation terminates       */
  DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_3_Ch_11,
    DIO_CLASS1_HIGH_Z_ON);
  DioCl1DigOut_write(pRTIDioC1DigOut_Port_3_Ch_11);

  /* Terminate for S-Function (rti_commonblock): '<S76>/S-Function1' incorporates:
   *  Constant: '<S30>/VCC3'
   */

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup3 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 13 --- */

  /* disable digital output channel 13-13 on port 3 *
   * (set to high-impedance), when the simulation terminates       */
  DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_3_Ch_13,
    DIO_CLASS1_HIGH_Z_ON);
  DioCl1DigOut_write(pRTIDioC1DigOut_Port_3_Ch_13);

  /* Terminate for S-Function (rti_commonblock): '<S81>/S-Function1' incorporates:
   *  Constant: '<S31>/Constant'
   */

  /* --- Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_OUT_BL1 --- */
  /* --- [RTI120X, BITOUT] - Port: 1 - Channel: 1 --- */

  /* disable digital output channel 1-1 on port 1 *
   * (set to high-impedance), when the simulation terminates       */
  DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_1_Ch_1,
    DIO_CLASS1_HIGH_Z_ON);
  DioCl1DigOut_write(pRTIDioC1DigOut_Port_1_Ch_1);

  /* Terminate for S-Function (rti_commonblock): '<S37>/S-Function1' */

  /* dSPACE I/O Board DS1202SER #1 Unit:GENSER Group:SETUP */
  dsser_disable(rtiDS1202SER_B1_Ser[0]);
  dsser_fifo_reset(rtiDS1202SER_B1_Ser[0]);

  /* Terminate for S-Function (rti_commonblock): '<S90>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:1 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X1])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S83>/S-Function1' */

  /* dSPACE RTICAN STD Srvc-Message Block */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (CANTP1_RX_SPMSG_M1_C1_STD)) == DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S84>/S-Function1' */

  /* dSPACE RTICAN STD Srvc-Message Block */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (CANTP1_RX_SPMSG_M1_C2_STD)) == DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S88>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:100 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][1] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S89>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:100 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][1] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }
}
