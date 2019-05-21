/*
 * Human_in_Loop.c
 *
 * Code generation for model "Human_in_Loop".
 *
 * Model version              : 1.1216
 * Simulink Coder version : 8.13 (R2017b) 24-Jul-2017
 * C source code generated on : Tue May 21 20:27:12 2019
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
static real_T Human_in_Loo_eml_rand_mt19937ar(uint32_T state[625]);
static void Human_in_Loop_rand(real_T r[2]);
static void nested_function_nested_function(coder_internal_ref_Human_in_L_T
  *varargin_1, const real_T varargin_2[60], real_T varargin_3,
  coder_internal_ref_Human_in_L_T **this_environment_f1, real_T
  this_environment_f2[60], real_T *this_environment_f3);
static void nested_function_nested_functi_f(coder_internal_ref_Human_in_L_T
  *varargin_1_environment_f1, const real_T varargin_1_environment_f2[60], real_T
  varargin_1_environment_f3, const real_T varargin_2[900], real_T varargin_3,
  const real_T varargin_5[30], coder_internal_ref_Human_in_L_T
  **this_environment_f1_environment, real_T this_environment_f1_environme_0[60],
  real_T *this_environment_f1_environme_1, real_T this_environment_f2[900],
  real_T *this_environment_f3, real_T this_environment_f4[30]);
static void nested_function_nested_funct_fs(coder_internal_ref_Human_in_L_T
  *varargin_1, coder_internal_ref_Human_in_L_T *varargin_2_environment_f1, const
  real_T varargin_2_environment_f2[60], real_T varargin_2_environment_f3, const
  real_T varargin_3[900], real_T varargin_4, coder_internal_ref_Human_in_L_T
  **this_environment_f1, coder_internal_ref_Human_in_L_T
  **this_environment_f2_environment, real_T this_environment_f2_environme_0[60],
  real_T *this_environment_f2_environme_1, real_T this_environment_f3[900],
  real_T *this_environment_f4);
static void nested_function_nested_func_fst(coder_internal_ref_Human_in_L_T
  *varargin_1_environment_f1_envir, const real_T
  varargin_1_environment_f1_env_0[60], real_T varargin_1_environment_f1_env_1,
  const real_T varargin_1_environment_f2[900], real_T varargin_1_environment_f3,
  const real_T varargin_1_environment_f4[30], coder_internal_ref_Human_in_L_T
  *varargin_3_environment_f1, coder_internal_ref_Human_in_L_T
  *varargin_3_environment_f2_envir, const real_T
  varargin_3_environment_f2_env_0[60], real_T varargin_3_environment_f2_env_1,
  const real_T varargin_3_environment_f3[900], real_T varargin_3_environment_f4,
  cell_2_Human_in_Loop_T *this_environment_f1_environment,
  cell_3_Human_in_Loop_T *this_environment_f2_environment);
static real_T Human_in_Loop_kernel(const coder_internal_ref_Human_in_L_T
  *b_Theta, const real_T x1[2], const real_T x2[2]);
static void Human_in_Loop_KKernel(const coder_internal_ref_Human_in_L_T *b_Theta,
  const real_T x1[2], const real_T x2_data[], const int32_T x2_size[2], real_T
  k_data[], int32_T k_size[2]);
static void Human_in_Loop_eye(real_T varargin_1, real_T I_data[], int32_T
  I_size[2]);
static void Human__eml_signed_integer_colon(int32_T b, int32_T y_data[], int32_T
  y_size[2]);
static int32_T Human_in_Loop_ixamax(int32_T n, const real_T x_data[], int32_T
  ix0);
static void Human_in_Loop_xswap_a(int32_T n, real_T x_data[], int32_T ix0,
  int32_T iy0);
static real_T Human_in_Loop_xnrm2_m(int32_T n, const real_T x_data[], int32_T
  ix0);
static real_T Human_in_Loop_xzlarfg(int32_T n, real_T *alpha1, real_T x_data[],
  int32_T ix0);
static int32_T Human_in_Loop_ilazlc(int32_T m, int32_T n, const real_T A_data[],
  int32_T ia0, int32_T lda);
static void Human_in_Loop_xgemv(int32_T m, int32_T n, const real_T A_data[],
  int32_T ia0, int32_T lda, const real_T x_data[], int32_T ix0, real_T y_data[]);
static void Human_in_Loop_xgerc(int32_T m, int32_T n, real_T alpha1, int32_T ix0,
  const real_T y_data[], real_T A_data[], int32_T ia0, int32_T lda);
static void Human_in_Loop_xzlarf(int32_T m, int32_T n, int32_T iv0, real_T tau,
  real_T C_data[], int32_T ic0, int32_T ldc, real_T work_data[]);
static void Human_in_Loop_xgeqp3_e(real_T A_data[], int32_T A_size[2], real_T
  tau_data[], int32_T *tau_size, int32_T jpvt_data[], int32_T jpvt_size[2]);
static int32_T Human_in_Loop_rankFromQR(const real_T A_data[], const int32_T
  A_size[2]);
static void Human_in_Loop_LSQFromQR(const real_T A_data[], const int32_T A_size
  [2], const real_T tau_data[], const int32_T jpvt_data[], real_T B_data[],
  int32_T rankA, real_T Y_data[], int32_T *Y_size);
static void Human_in_Loop_qrsolve(const real_T A_data[], const int32_T A_size[2],
  const real_T B_data[], const int32_T *B_size, real_T Y_data[], int32_T *Y_size);
static void Human_in_Loop_xswap(int32_T n, real_T x_data[], int32_T ix0, int32_T
  incx, int32_T iy0, int32_T incy);
static void Human_in_Loop_xgetrf(int32_T m, int32_T n, real_T A_data[], int32_T
  A_size[2], int32_T lda, int32_T ipiv_data[], int32_T ipiv_size[2], int32_T
  *info);
static void Human_in_Loop_xtrsm(int32_T n, const real_T A_data[], int32_T lda,
  real_T B_data[]);
static void Human_in_Loop_xtrsm_e(int32_T n, const real_T A_data[], int32_T lda,
  real_T B_data[]);
static void Human_in_Loop_lusolve(const real_T A_data[], const int32_T A_size[2],
  const real_T B_data[], const int32_T B_size[2], real_T X_data[], int32_T
  X_size[2]);
static void Human_in_Loop_mrdivide(const real_T A_data[], const int32_T A_size[2],
  const real_T B_data[], const int32_T B_size[2], real_T y_data[], int32_T
  y_size[2]);
static real_T Human_in_Loop___anon_fcn(const coder_internal_ref_Human_in_L_T
  *k_environment_f1, const real_T k_environment_f2[60], real_T k_environment_f3,
  const real_T K[900], real_T iter, const real_T Y[30], const real_T x[2]);
static real_T Human_in_Loop___anon_fcn_nf(const coder_internal_ref_Human_in_L_T *
  b_Theta, const coder_internal_ref_Human_in_L_T *k_environment_f1, const real_T
  k_environment_f2[60], real_T k_environment_f3, const real_T K[900], real_T
  iter, const real_T x[2]);
static real_T Human_in_Loop___anon_fcn_n(const coder_internal_ref_Human_in_L_T
  *GP_mu_environment_f1_environmen, const real_T
  GP_mu_environment_f1_environm_0[60], real_T GP_mu_environment_f1_environm_1,
  const real_T GP_mu_environment_f2[900], real_T GP_mu_environment_f3, const
  real_T GP_mu_environment_f4[30], const coder_internal_ref_Human_in_L_T
  *GP_sigma_environment_f1, const coder_internal_ref_Human_in_L_T
  *GP_sigma_environment_f2_environ, const real_T
  GP_sigma_environment_f2_envir_0[60], real_T GP_sigma_environment_f2_envir_1,
  const real_T GP_sigma_environment_f3[900], real_T GP_sigma_environment_f4,
  const real_T x[2]);

/* Forward declaration for local functions */
static void Human_in_Loop_mldivide_g(const real_T A[16], real_T B[4]);
static void Human_in_Loop_power(const real_T a_data[], const int32_T a_size[2],
  real_T y_data[], int32_T y_size[2]);
static void Human_in_Loop_power_a(const real_T a_data[], const int32_T a_size[2],
  real_T y_data[], int32_T y_size[2]);

/* Forward declaration for local functions */
static void Human_in_Loop_eye_n(real_T I_data[], int32_T I_size[2]);
static void Human_in_Loop_log(real_T x_data[], int32_T x_size[2]);
static real_T Human_in_Loop_sum(const real_T x_data[], const int32_T *x_size);
static void Human_in_Loop_power_d(const real_T a_data[], const int32_T *a_size,
  real_T y_data[], int32_T *y_size);
static void Human_in_Loop_sort(real_T x_data[], int32_T idx_data[], int32_T
  idx_size[2]);
static void Human_in_Loop_diag(const real_T v_data[], const int32_T *v_size,
  real_T d_data[], int32_T d_size[2]);
static void Human_in_Loop_triu(real_T x_data[]);
static void Human_in_Loop_eig(const real_T A_data[], creal_T V_data[], int32_T
  V_size[2], creal_T D_data[], int32_T D_size[2]);
static void Human_in_Loop_sqrt(real_T x_data[]);
static void Human_in__genrand_uint32_vector(uint32_T mt[625], uint32_T u[2]);
static real_T Human_in_Loop_genrandu(uint32_T mt[625]);
static real_T Human_in_L_eml_rand_mt19937ar_h(uint32_T state[625]);
static void Human_in_Loop_randn(real_T r_data[], int32_T *r_size);
static boolean_T Human_in_Loop_sortLE(const real_T v_data[], int32_T irow1,
  int32_T irow2);
static void Human_in_Loop_sortrows(real_T y_data[], real_T ndx_data[], int32_T
  *ndx_size);
static void Human_in_Loop_mean(const real_T x_data[], const int32_T x_size[2],
  real_T y[26]);
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

/* Function for MATLAB Function: '<S30>/Estimation' */
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

/* Function for MATLAB Function: '<S30>/Estimation' */
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

/* Function for MATLAB Function: '<S30>/Estimation' */
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

/* System initialize for function-call system: '<S26>/Serial Decoding System' */
void Human_SerialDecodingSystem_Init(void)
{
  int32_T i;

  /* SystemInitialize for MATLAB Function: '<S30>/Estimation' */
  Human_in_Loop_DW.last_count = -1.0;
  for (i = 0; i < 20; i++) {
    Human_in_Loop_DW.x[i] = 10.0 * (real_T)i + 10.0;
    Human_in_Loop_DW.y[i] = 0.0;
  }

  Human_in_Loop_DW.theta[0] = 0.0;
  Human_in_Loop_DW.theta[1] = 0.0;
  Human_in_Loop_DW.num = 1.0;

  /* End of SystemInitialize for MATLAB Function: '<S30>/Estimation' */

  /* SystemInitialize for Outport: '<S30>/E' */
  Human_in_Loop_B.E = Human_in_Loop_P.E_Y0;
}

/* System reset for function-call system: '<S26>/Serial Decoding System' */
void Huma_SerialDecodingSystem_Reset(void)
{
  int32_T i;

  /* SystemReset for MATLAB Function: '<S30>/Estimation' */
  Human_in_Loop_DW.last_count = -1.0;
  for (i = 0; i < 20; i++) {
    Human_in_Loop_DW.x[i] = 10.0 * (real_T)i + 10.0;
    Human_in_Loop_DW.y[i] = 0.0;
  }

  Human_in_Loop_DW.theta[0] = 0.0;
  Human_in_Loop_DW.theta[1] = 0.0;
  Human_in_Loop_DW.num = 1.0;

  /* End of SystemReset for MATLAB Function: '<S30>/Estimation' */
}

/* Output and update for function-call system: '<S26>/Serial Decoding System' */
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

  /* S-Function (rti_commonblock): '<S31>/S-Function1' */
  /* This comment workarounds a code generation problem */

  /* dSPACE I/O Board DS1202SER #1 Unit:GENSER Group:RECEIVE */
  {
    Human_in_Loop_B.SFunction1_o3_ls = dsser_receive(rtiDS1202SER_B1_Ser[0], 16U,
      (UInt8 *) &Human_in_Loop_B.SFunction1_o1_l[0], (UInt32 *)
      &Human_in_Loop_B.SFunction1_o2_bb);
  }

  /* MATLAB Function: '<S30>/MATLAB Function1' */
  /* MATLAB Function 'Opt Module/Cosmed/COSMED/Serial Decoding System/MATLAB Function1': '<S33>:1' */
  /* '<S33>:1:3' */
  Data1 = 0.0;

  /* '<S33>:1:4' */
  Data2 = 0.0;

  /* '<S33>:1:6' */
  for (i = 1.0; Human_in_Loop_B.SFunction1_o1_l[(int32_T)i - 1] != 44; i++) {
    /* '<S33>:1:8' */
    /* '<S33>:1:9' */
    Data1 = ((real_T)Human_in_Loop_B.SFunction1_o1_l[(int32_T)i - 1] - 48.0) +
      Data1 * 10.0;

    /* '<S33>:1:10' */
  }

  /* '<S33>:1:13' */
  for (i++; i < Human_in_Loop_B.SFunction1_o2_bb; i++) {
    /* '<S33>:1:15' */
    /* '<S33>:1:16' */
    Data2 = ((real_T)Human_in_Loop_B.SFunction1_o1_l[(int32_T)i - 1] - 48.0) +
      Data2 * 10.0;

    /* '<S33>:1:17' */
  }

  Human_in_Loop_B.Data1 = Data1;
  Human_in_Loop_B.Data2 = Data2;

  /* End of MATLAB Function: '<S30>/MATLAB Function1' */

  /* MATLAB Function: '<S30>/Estimation' */
  /* MATLAB Function 'Opt Module/Cosmed/COSMED/Serial Decoding System/Estimation': '<S32>:1' */
  /* '<S32>:1:49' */
  if ((Human_in_Loop_DW.last_count != Human_in_Loop_B.RT2_g) ||
      (Human_in_Loop_B.RT2_g == 0.0)) {
    /* '<S32>:1:21' */
    /* '<S32>:1:22' */
    Human_in_Loop_DW.last_count = Human_in_Loop_B.RT2_g;

    /* '<S32>:1:23' */
    /* '<S32>:1:24' */
    for (i_0 = 0; i_0 < 20; i_0++) {
      Human_in_Loop_DW.x[i_0] = 10.0 * (real_T)i_0 + 10.0;
      Human_in_Loop_DW.y[i_0] = 0.0;
    }

    /* '<S32>:1:25' */
    Human_in_Loop_DW.theta[0] = 0.0;
    Human_in_Loop_DW.theta[1] = 0.0;

    /* '<S32>:1:26' */
    Human_in_Loop_DW.num = 1.0;
  } else {
    /* '<S32>:1:28' */
    Human_in_Loop_DW.num++;
  }

  if (Human_in_Loop_B.RT2_g != 0.0) {
    /* '<S32>:1:32' */
    /* '<S32>:1:33' */
    Human_in_Loop_DW.x[(int32_T)Human_in_Loop_DW.num - 1] =
      Human_in_Loop_B.RT1_c;

    /* '<S32>:1:34' */
    Human_in_Loop_DW.y[(int32_T)Human_in_Loop_DW.num - 1] = 0.278 *
      Human_in_Loop_B.Data1 + 0.075 * Human_in_Loop_B.Data2;
    if (Human_in_Loop_DW.num >= 2.0) {
      /* '<S32>:1:37' */
      /* '<S32>:1:38' */
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

      /* '<S32>:1:39' */
      /* '<S32>:1:40' */
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

  /* '<S32>:1:45' */
  /* '<S32>:1:47' */
  memcpy(&Human_in_Loop_B.Time_p[0], &Human_in_Loop_DW.x[0], 12U * sizeof(real_T));

  /* '<S32>:1:48' */
  memcpy(&Human_in_Loop_B.Y[0], &Human_in_Loop_DW.y[0], 12U * sizeof(real_T));

  /* '<S32>:1:49' */
  for (i_0 = 0; i_0 < 12; i_0++) {
    b_x[i_0] = -Human_in_Loop_B.Time_p[i_0] / 42.0;
  }

  for (loop_ub = 0; loop_ub < 12; loop_ub++) {
    Data1 = b_x[loop_ub];
    Data1 = exp(Data1);
    a[loop_ub] = 1.0 - Data1;
    a[loop_ub + 12] = 1.0;
  }

  /* '<S32>:1:50' */
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

  /* End of MATLAB Function: '<S30>/Estimation' */
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static real_T Human_in_Loo_eml_rand_mt19937ar(uint32_T state[625])
{
  real_T r;
  uint32_T u[2];
  uint32_T mti;
  uint32_T y;
  int32_T j;
  int32_T kk;

  /* ========================= COPYRIGHT NOTICE ============================ */
  /*  This is a uniform (0,1) pseudorandom number generator based on:        */
  /*                                                                         */
  /*  A C-program for MT19937, with initialization improved 2002/1/26.       */
  /*  Coded by Takuji Nishimura and Makoto Matsumoto.                        */
  /*                                                                         */
  /*  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,      */
  /*  All rights reserved.                                                   */
  /*                                                                         */
  /*  Redistribution and use in source and binary forms, with or without     */
  /*  modification, are permitted provided that the following conditions     */
  /*  are met:                                                               */
  /*                                                                         */
  /*    1. Redistributions of source code must retain the above copyright    */
  /*       notice, this list of conditions and the following disclaimer.     */
  /*                                                                         */
  /*    2. Redistributions in binary form must reproduce the above copyright */
  /*       notice, this list of conditions and the following disclaimer      */
  /*       in the documentation and/or other materials provided with the     */
  /*       distribution.                                                     */
  /*                                                                         */
  /*    3. The names of its contributors may not be used to endorse or       */
  /*       promote products derived from this software without specific      */
  /*       prior written permission.                                         */
  /*                                                                         */
  /*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS    */
  /*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT      */
  /*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR  */
  /*  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT  */
  /*  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,  */
  /*  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT       */
  /*  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,  */
  /*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY  */
  /*  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT    */
  /*  (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE */
  /*  OF THIS  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  */
  /*                                                                         */
  /* =============================   END   ================================= */
  do {
    for (j = 0; j < 2; j++) {
      mti = state[624] + 1U;
      if (mti >= 625U) {
        for (kk = 0; kk < 227; kk++) {
          y = (state[kk + 1] & 2147483647U) | (state[kk] & 2147483648U);
          if ((int32_T)(y & 1U) == 0) {
            mti = y >> 1U;
          } else {
            mti = y >> 1U ^ 2567483615U;
          }

          state[kk] = state[kk + 397] ^ mti;
        }

        for (kk = 0; kk < 396; kk++) {
          y = (state[kk + 227] & 2147483648U) | (state[kk + 228] & 2147483647U);
          if ((int32_T)(y & 1U) == 0) {
            mti = y >> 1U;
          } else {
            mti = y >> 1U ^ 2567483615U;
          }

          state[kk + 227] = state[kk] ^ mti;
        }

        y = (state[623] & 2147483648U) | (state[0] & 2147483647U);
        if ((int32_T)(y & 1U) == 0) {
          mti = y >> 1U;
        } else {
          mti = y >> 1U ^ 2567483615U;
        }

        state[623] = state[396] ^ mti;
        mti = 1U;
      }

      y = state[(int32_T)mti - 1];
      state[624] = mti;
      y ^= y >> 11U;
      y ^= y << 7U & 2636928640U;
      y ^= y << 15U & 4022730752U;
      y ^= y >> 18U;
      u[j] = y;
    }

    r = ((real_T)(u[0] >> 5U) * 6.7108864E+7 + (real_T)(u[1] >> 6U)) *
      1.1102230246251565E-16;
  } while (r == 0.0);

  return r;
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_rand(real_T r[2])
{
  real_T b;
  b = Human_in_Loo_eml_rand_mt19937ar(Human_in_Loop_DW.state);
  r[0] = b;
  b = Human_in_Loo_eml_rand_mt19937ar(Human_in_Loop_DW.state);
  r[1] = b;
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void nested_function_nested_function(coder_internal_ref_Human_in_L_T
  *varargin_1, const real_T varargin_2[60], real_T varargin_3,
  coder_internal_ref_Human_in_L_T **this_environment_f1, real_T
  this_environment_f2[60], real_T *this_environment_f3)
{
  coder_internal_ref_Human_in_L_T *environment_f1;
  environment_f1 = varargin_1;
  memcpy(&this_environment_f2[0], &varargin_2[0], 60U * sizeof(real_T));
  *this_environment_f1 = environment_f1;
  *this_environment_f3 = varargin_3;
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void nested_function_nested_functi_f(coder_internal_ref_Human_in_L_T
  *varargin_1_environment_f1, const real_T varargin_1_environment_f2[60], real_T
  varargin_1_environment_f3, const real_T varargin_2[900], real_T varargin_3,
  const real_T varargin_5[30], coder_internal_ref_Human_in_L_T
  **this_environment_f1_environment, real_T this_environment_f1_environme_0[60],
  real_T *this_environment_f1_environme_1, real_T this_environment_f2[900],
  real_T *this_environment_f3, real_T this_environment_f4[30])
{
  coder_internal_ref_Human_in_L_T *environment_f1_environment_f1;
  environment_f1_environment_f1 = varargin_1_environment_f1;
  memcpy(&this_environment_f1_environme_0[0], &varargin_1_environment_f2[0], 60U
         * sizeof(real_T));
  memcpy(&this_environment_f2[0], &varargin_2[0], 900U * sizeof(real_T));
  *this_environment_f1_environment = environment_f1_environment_f1;
  *this_environment_f1_environme_1 = varargin_1_environment_f3;
  *this_environment_f3 = varargin_3;
  memcpy(&this_environment_f4[0], &varargin_5[0], 30U * sizeof(real_T));
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void nested_function_nested_funct_fs(coder_internal_ref_Human_in_L_T
  *varargin_1, coder_internal_ref_Human_in_L_T *varargin_2_environment_f1, const
  real_T varargin_2_environment_f2[60], real_T varargin_2_environment_f3, const
  real_T varargin_3[900], real_T varargin_4, coder_internal_ref_Human_in_L_T
  **this_environment_f1, coder_internal_ref_Human_in_L_T
  **this_environment_f2_environment, real_T this_environment_f2_environme_0[60],
  real_T *this_environment_f2_environme_1, real_T this_environment_f3[900],
  real_T *this_environment_f4)
{
  coder_internal_ref_Human_in_L_T *environment_f1;
  coder_internal_ref_Human_in_L_T *environment_f2_environment_f1;
  environment_f1 = varargin_1;
  environment_f2_environment_f1 = varargin_2_environment_f1;
  memcpy(&this_environment_f2_environme_0[0], &varargin_2_environment_f2[0], 60U
         * sizeof(real_T));
  memcpy(&this_environment_f3[0], &varargin_3[0], 900U * sizeof(real_T));
  *this_environment_f1 = environment_f1;
  *this_environment_f2_environment = environment_f2_environment_f1;
  *this_environment_f2_environme_1 = varargin_2_environment_f3;
  *this_environment_f4 = varargin_4;
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void nested_function_nested_func_fst(coder_internal_ref_Human_in_L_T
  *varargin_1_environment_f1_envir, const real_T
  varargin_1_environment_f1_env_0[60], real_T varargin_1_environment_f1_env_1,
  const real_T varargin_1_environment_f2[900], real_T varargin_1_environment_f3,
  const real_T varargin_1_environment_f4[30], coder_internal_ref_Human_in_L_T
  *varargin_3_environment_f1, coder_internal_ref_Human_in_L_T
  *varargin_3_environment_f2_envir, const real_T
  varargin_3_environment_f2_env_0[60], real_T varargin_3_environment_f2_env_1,
  const real_T varargin_3_environment_f3[900], real_T varargin_3_environment_f4,
  cell_2_Human_in_Loop_T *this_environment_f1_environment,
  cell_3_Human_in_Loop_T *this_environment_f2_environment)
{
  coder_internal_ref_Human_in_L_T *environment_f1_environment_f1_e;
  real_T environment_f1_environment_f1_0[60];
  coder_internal_ref_Human_in_L_T *environment_f2_environment_f1;
  coder_internal_ref_Human_in_L_T *environment_f2_environment_f2_e;
  real_T environment_f2_environment_f2_0[60];
  environment_f1_environment_f1_e = varargin_1_environment_f1_envir;
  memcpy(&environment_f1_environment_f1_0[0], &varargin_1_environment_f1_env_0[0],
         60U * sizeof(real_T));
  memcpy(&Human_in_Loop_B.environment_f1_environment_f2[0],
         &varargin_1_environment_f2[0], 900U * sizeof(real_T));
  environment_f2_environment_f1 = varargin_3_environment_f1;
  environment_f2_environment_f2_e = varargin_3_environment_f2_envir;
  memcpy(&environment_f2_environment_f2_0[0], &varargin_3_environment_f2_env_0[0],
         60U * sizeof(real_T));
  memcpy(&Human_in_Loop_B.environment_f2_environment_f3[0],
         &varargin_3_environment_f3[0], 900U * sizeof(real_T));
  this_environment_f1_environment->f1.environment.f1 =
    environment_f1_environment_f1_e;
  memcpy(&this_environment_f1_environment->f1.environment.f2[0],
         &environment_f1_environment_f1_0[0], 60U * sizeof(real_T));
  this_environment_f1_environment->f1.environment.f3 =
    varargin_1_environment_f1_env_1;
  memcpy(&this_environment_f1_environment->f2[0],
         &Human_in_Loop_B.environment_f1_environment_f2[0], 900U * sizeof(real_T));
  this_environment_f1_environment->f3 = varargin_1_environment_f3;
  memcpy(&this_environment_f1_environment->f4[0], &varargin_1_environment_f4[0],
         30U * sizeof(real_T));
  this_environment_f2_environment->f1 = environment_f2_environment_f1;
  this_environment_f2_environment->f2.environment.f1 =
    environment_f2_environment_f2_e;
  memcpy(&this_environment_f2_environment->f2.environment.f2[0],
         &environment_f2_environment_f2_0[0], 60U * sizeof(real_T));
  this_environment_f2_environment->f2.environment.f3 =
    varargin_3_environment_f2_env_1;
  memcpy(&this_environment_f2_environment->f3[0],
         &Human_in_Loop_B.environment_f2_environment_f3[0], 900U * sizeof(real_T));
  this_environment_f2_environment->f4 = varargin_3_environment_f4;
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static real_T Human_in_Loop_kernel(const coder_internal_ref_Human_in_L_T
  *b_Theta, const real_T x1[2], const real_T x2[2])
{
  real_T theta[2];
  real_T d[4];
  real_T a21;
  real_T z1_idx_0;
  real_T z1_idx_1;

  /* '<S36>:1:88' */
  a21 = b_Theta->contents;

  /* '<S36>:1:89' */
  z1_idx_0 = a21 * a21;
  z1_idx_1 = a21 * a21;
  a21 = z1_idx_0;
  d[0] = a21;
  a21 = -(x1[0] - x2[0]);
  z1_idx_0 = a21;
  a21 = z1_idx_1;
  d[3] = a21;
  a21 = -(x1[1] - x2[1]);
  z1_idx_1 = a21;
  a21 = 0.0 / d[0];
  theta[0] = z1_idx_0 / d[0];
  theta[1] = (z1_idx_1 - theta[0] * 0.0) / (d[3] - a21 * 0.0);
  theta[0] -= theta[1] * a21;
  z1_idx_0 = (x1[0] - x2[0]) * theta[0];
  z1_idx_0 += (x1[1] - x2[1]) * theta[1];
  return exp(z1_idx_0 / 2.0);
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_KKernel(const coder_internal_ref_Human_in_L_T *b_Theta,
  const real_T x1[2], const real_T x2_data[], const int32_T x2_size[2], real_T
  k_data[], int32_T k_size[2])
{
  int32_T varargin_2;
  int32_T tmp;

  /* '<S36>:1:97' */
  /* '<S36>:1:95' */
  varargin_2 = x2_size[1];
  k_size[0] = 1;
  k_size[1] = varargin_2;
  if (0 <= varargin_2 - 1) {
    memset(&k_data[0], 0, varargin_2 * sizeof(real_T));
  }

  /* '<S36>:1:96' */
  varargin_2 = 0;
  tmp = x2_size[1] - 1;
  while (varargin_2 <= tmp) {
    /* '<S36>:1:96' */
    /* '<S36>:1:97' */
    k_data[varargin_2] = Human_in_Loop_kernel(b_Theta, x1, &x2_data[x2_size[0] *
      varargin_2]);
    varargin_2++;
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_eye(real_T varargin_1, real_T I_data[], int32_T
  I_size[2])
{
  int32_T loop_ub;
  I_size[0] = (int32_T)varargin_1;
  I_size[1] = (int32_T)varargin_1;
  loop_ub = (int32_T)varargin_1 * (int32_T)varargin_1;
  if (0 <= loop_ub - 1) {
    memset(&I_data[0], 0, loop_ub * sizeof(real_T));
  }

  for (loop_ub = 0; loop_ub + 1 <= (int32_T)varargin_1; loop_ub++) {
    I_data[loop_ub + I_size[0] * loop_ub] = 1.0;
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human__eml_signed_integer_colon(int32_T b, int32_T y_data[], int32_T
  y_size[2])
{
  int32_T yk;
  int32_T k;
  y_size[0] = 1;
  y_size[1] = b;
  y_data[0] = 1;
  yk = 1;
  for (k = 2; k <= b; k++) {
    yk++;
    y_data[k - 1] = yk;
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static int32_T Human_in_Loop_ixamax(int32_T n, const real_T x_data[], int32_T
  ix0)
{
  int32_T idxmax;
  int32_T ix;
  real_T smax;
  real_T s;
  int32_T k;
  if (n < 1) {
    idxmax = 0;
  } else {
    idxmax = 1;
    if (n > 1) {
      ix = ix0 - 1;
      smax = fabs(x_data[ix0 - 1]);
      for (k = 2; k <= n; k++) {
        ix++;
        s = fabs(x_data[ix]);
        if (s > smax) {
          idxmax = k;
          smax = s;
        }
      }
    }
  }

  return idxmax;
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_xswap_a(int32_T n, real_T x_data[], int32_T ix0,
  int32_T iy0)
{
  int32_T ix;
  int32_T iy;
  real_T temp;
  int32_T k;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 1; k <= n; k++) {
    temp = x_data[ix];
    x_data[ix] = x_data[iy];
    x_data[iy] = temp;
    ix++;
    iy++;
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static real_T Human_in_Loop_xnrm2_m(int32_T n, const real_T x_data[], int32_T
  ix0)
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

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static real_T Human_in_Loop_xzlarfg(int32_T n, real_T *alpha1, real_T x_data[],
  int32_T ix0)
{
  real_T tau;
  real_T xnorm;
  int32_T knt;
  int32_T b_k;
  int32_T c_k;
  tau = 0.0;
  if (!(n <= 0)) {
    xnorm = Human_in_Loop_xnrm2_m(n - 1, x_data, ix0);
    if (xnorm != 0.0) {
      xnorm = rt_hypotd_snf(*alpha1, xnorm);
      if (*alpha1 >= 0.0) {
        xnorm = -xnorm;
      }

      if (fabs(xnorm) < 1.0020841800044864E-292) {
        knt = 0;
        b_k = (ix0 + n) - 2;
        do {
          knt++;
          for (c_k = ix0; c_k <= b_k; c_k++) {
            x_data[c_k - 1] *= 9.9792015476736E+291;
          }

          xnorm *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while (!(fabs(xnorm) >= 1.0020841800044864E-292));

        xnorm = rt_hypotd_snf(*alpha1, Human_in_Loop_xnrm2_m(n - 1, x_data, ix0));
        if (*alpha1 >= 0.0) {
          xnorm = -xnorm;
        }

        tau = (xnorm - *alpha1) / xnorm;
        *alpha1 = 1.0 / (*alpha1 - xnorm);
        b_k = (ix0 + n) - 2;
        for (c_k = ix0; c_k <= b_k; c_k++) {
          x_data[c_k - 1] *= *alpha1;
        }

        for (b_k = 1; b_k <= knt; b_k++) {
          xnorm *= 1.0020841800044864E-292;
        }

        *alpha1 = xnorm;
      } else {
        tau = (xnorm - *alpha1) / xnorm;
        *alpha1 = 1.0 / (*alpha1 - xnorm);
        knt = (ix0 + n) - 2;
        for (b_k = ix0; b_k <= knt; b_k++) {
          x_data[b_k - 1] *= *alpha1;
        }

        *alpha1 = xnorm;
      }
    }
  }

  return tau;
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static int32_T Human_in_Loop_ilazlc(int32_T m, int32_T n, const real_T A_data[],
  int32_T ia0, int32_T lda)
{
  int32_T j;
  int32_T coltop;
  int32_T ia;
  int32_T exitg1;
  boolean_T exitg2;
  j = n;
  exitg2 = false;
  while ((!exitg2) && (j > 0)) {
    coltop = (j - 1) * lda + ia0;
    ia = coltop;
    do {
      exitg1 = 0;
      if (ia <= (coltop + m) - 1) {
        if (A_data[ia - 1] != 0.0) {
          exitg1 = 1;
        } else {
          ia++;
        }
      } else {
        j--;
        exitg1 = 2;
      }
    } while (exitg1 == 0);

    if (exitg1 == 1) {
      exitg2 = true;
    }
  }

  return j;
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_xgemv(int32_T m, int32_T n, const real_T A_data[],
  int32_T ia0, int32_T lda, const real_T x_data[], int32_T ix0, real_T y_data[])
{
  int32_T ix;
  real_T c;
  int32_T iy;
  int32_T b;
  int32_T iac;
  int32_T d;
  int32_T ia;
  if (n != 0) {
    for (iy = 1; iy <= n; iy++) {
      y_data[iy - 1] = 0.0;
    }

    iy = 0;
    b = (n - 1) * lda + ia0;
    for (iac = ia0; iac <= b; iac += lda) {
      ix = ix0;
      c = 0.0;
      d = (iac + m) - 1;
      for (ia = iac; ia <= d; ia++) {
        c += A_data[ia - 1] * x_data[ix - 1];
        ix++;
      }

      y_data[iy] += c;
      iy++;
    }
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_xgerc(int32_T m, int32_T n, real_T alpha1, int32_T ix0,
  const real_T y_data[], real_T A_data[], int32_T ia0, int32_T lda)
{
  int32_T jA;
  int32_T jy;
  real_T temp;
  int32_T ix;
  int32_T j;
  int32_T b;
  int32_T ijA;
  if (!(alpha1 == 0.0)) {
    jA = ia0 - 1;
    jy = 0;
    for (j = 1; j <= n; j++) {
      if (y_data[jy] != 0.0) {
        temp = y_data[jy] * alpha1;
        ix = ix0;
        b = m + jA;
        for (ijA = jA; ijA + 1 <= b; ijA++) {
          A_data[ijA] += A_data[ix - 1] * temp;
          ix++;
        }
      }

      jy++;
      jA += lda;
    }
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_xzlarf(int32_T m, int32_T n, int32_T iv0, real_T tau,
  real_T C_data[], int32_T ic0, int32_T ldc, real_T work_data[])
{
  int32_T lastv;
  int32_T lastc;
  if (tau != 0.0) {
    lastv = m;
    lastc = iv0 + m;
    while ((lastv > 0) && (C_data[lastc - 2] == 0.0)) {
      lastv--;
      lastc--;
    }

    lastc = Human_in_Loop_ilazlc(lastv, n, C_data, ic0, ldc);
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    Human_in_Loop_xgemv(lastv, lastc, C_data, ic0, ldc, C_data, iv0, work_data);
    Human_in_Loop_xgerc(lastv, lastc, -tau, iv0, work_data, C_data, ic0, ldc);
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_xgeqp3_e(real_T A_data[], int32_T A_size[2], real_T
  tau_data[], int32_T *tau_size, int32_T jpvt_data[], int32_T jpvt_size[2])
{
  int32_T m;
  int32_T n;
  int32_T mn;
  real_T work_data[29];
  real_T vn1_data[29];
  real_T vn2_data[29];
  int32_T k;
  int32_T i_i;
  int32_T nmi;
  int32_T mmi;
  int32_T pvt;
  int32_T itemp;
  real_T temp2;
  real_T b_atmp;
  int8_T c_idx_0;
  m = A_size[0];
  n = A_size[1];
  k = A_size[0];
  mn = A_size[1];
  if (k < mn) {
    mn = k;
  }

  *tau_size = (int8_T)mn;
  Human__eml_signed_integer_colon(A_size[1], jpvt_data, jpvt_size);
  c_idx_0 = (int8_T)A_size[1];
  k = c_idx_0;
  if (0 <= k - 1) {
    memset(&work_data[0], 0, k * sizeof(real_T));
  }

  k = 1;
  for (mmi = 0; mmi + 1 <= n; mmi++) {
    vn1_data[mmi] = Human_in_Loop_xnrm2_m(m, A_data, k);
    vn2_data[mmi] = vn1_data[mmi];
    k += m;
  }

  for (k = 0; k + 1 <= mn; k++) {
    i_i = k * m + k;
    nmi = (n - k) - 1;
    mmi = m - k;
    pvt = (Human_in_Loop_ixamax(1 + nmi, vn1_data, k + 1) + k) - 1;
    if (pvt + 1 != k + 1) {
      Human_in_Loop_xswap_a(m, A_data, 1 + m * pvt, 1 + m * k);
      itemp = jpvt_data[pvt];
      jpvt_data[pvt] = jpvt_data[k];
      jpvt_data[k] = itemp;
      vn1_data[pvt] = vn1_data[k];
      vn2_data[pvt] = vn2_data[k];
    }

    if (k + 1 < m) {
      b_atmp = A_data[i_i];
      temp2 = Human_in_Loop_xzlarfg(mmi, &b_atmp, A_data, i_i + 2);
      tau_data[k] = temp2;
      A_data[i_i] = b_atmp;
    } else {
      tau_data[k] = 0.0;
    }

    if (k + 1 < n) {
      b_atmp = A_data[i_i];
      A_data[i_i] = 1.0;
      Human_in_Loop_xzlarf(mmi, nmi, i_i + 1, tau_data[k], A_data, (k + (k + 1) *
        m) + 1, m, work_data);
      A_data[i_i] = b_atmp;
    }

    for (i_i = k + 1; i_i + 1 <= n; i_i++) {
      if (vn1_data[i_i] != 0.0) {
        b_atmp = fabs(A_data[A_size[0] * i_i + k]) / vn1_data[i_i];
        b_atmp = 1.0 - b_atmp * b_atmp;
        if (b_atmp < 0.0) {
          b_atmp = 0.0;
        }

        temp2 = vn1_data[i_i] / vn2_data[i_i];
        temp2 = temp2 * temp2 * b_atmp;
        if (temp2 <= 1.4901161193847656E-8) {
          if (k + 1 < m) {
            vn1_data[i_i] = Human_in_Loop_xnrm2_m(mmi - 1, A_data, (m * i_i + k)
              + 2);
            vn2_data[i_i] = vn1_data[i_i];
          } else {
            vn1_data[i_i] = 0.0;
            vn2_data[i_i] = 0.0;
          }
        } else {
          vn1_data[i_i] *= sqrt(b_atmp);
        }
      }
    }
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static int32_T Human_in_Loop_rankFromQR(const real_T A_data[], const int32_T
  A_size[2])
{
  int32_T r;
  real_T tol;
  int32_T minmn;
  int32_T maxmn;
  r = 0;
  if (A_size[0] < A_size[1]) {
    minmn = A_size[0];
    maxmn = A_size[1];
  } else {
    minmn = A_size[1];
    maxmn = A_size[0];
  }

  tol = (real_T)maxmn * fabs(A_data[0]) * 2.2204460492503131E-16;
  while ((r < minmn) && (!(fabs(A_data[A_size[0] * r + r]) <= tol))) {
    r++;
  }

  return r;
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_LSQFromQR(const real_T A_data[], const int32_T A_size
  [2], const real_T tau_data[], const int32_T jpvt_data[], real_T B_data[],
  int32_T rankA, real_T Y_data[], int32_T *Y_size)
{
  int32_T m;
  int32_T mn;
  real_T wj;
  int32_T c_i;
  int8_T b_idx_0;
  int32_T u0;
  b_idx_0 = (int8_T)A_size[1];
  *Y_size = b_idx_0;
  m = b_idx_0;
  if (0 <= m - 1) {
    memset(&Y_data[0], 0, m * sizeof(real_T));
  }

  m = A_size[0];
  u0 = A_size[0];
  mn = A_size[1];
  if (u0 < mn) {
    mn = u0;
  }

  for (u0 = 0; u0 + 1 <= mn; u0++) {
    if (tau_data[u0] != 0.0) {
      wj = B_data[u0];
      for (c_i = u0 + 1; c_i + 1 <= m; c_i++) {
        wj += A_data[A_size[0] * u0 + c_i] * B_data[c_i];
      }

      wj *= tau_data[u0];
      if (wj != 0.0) {
        B_data[u0] -= wj;
        for (c_i = u0 + 1; c_i + 1 <= m; c_i++) {
          B_data[c_i] -= A_data[A_size[0] * u0 + c_i] * wj;
        }
      }
    }
  }

  for (m = 0; m + 1 <= rankA; m++) {
    Y_data[jpvt_data[m] - 1] = B_data[m];
  }

  for (m = rankA - 1; m + 1 > 0; m--) {
    Y_data[jpvt_data[m] - 1] /= A_data[A_size[0] * m + m];
    for (mn = 0; mn + 1 <= m; mn++) {
      Y_data[jpvt_data[mn] - 1] -= A_data[A_size[0] * m + mn] *
        Y_data[jpvt_data[m] - 1];
    }
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_qrsolve(const real_T A_data[], const int32_T A_size[2],
  const real_T B_data[], const int32_T *B_size, real_T Y_data[], int32_T *Y_size)
{
  real_T tau_data[29];
  int32_T jpvt_data[29];
  int32_T rankR;
  real_T B_data_0[29];
  int32_T loop_ub;
  int32_T b_A_size[2];
  int32_T jpvt_size[2];
  b_A_size[0] = A_size[0];
  b_A_size[1] = A_size[1];
  loop_ub = A_size[0] * A_size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&Human_in_Loop_B.b_A_data_b[0], &A_data[0], loop_ub * sizeof(real_T));
  }

  Human_in_Loop_xgeqp3_e(Human_in_Loop_B.b_A_data_b, b_A_size, tau_data, &rankR,
    jpvt_data, jpvt_size);
  rankR = Human_in_Loop_rankFromQR(Human_in_Loop_B.b_A_data_b, b_A_size);
  loop_ub = *B_size;
  if (0 <= loop_ub - 1) {
    memcpy(&B_data_0[0], &B_data[0], loop_ub * sizeof(real_T));
  }

  Human_in_Loop_LSQFromQR(Human_in_Loop_B.b_A_data_b, b_A_size, tau_data,
    jpvt_data, B_data_0, rankR, Y_data, Y_size);
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_xswap(int32_T n, real_T x_data[], int32_T ix0, int32_T
  incx, int32_T iy0, int32_T incy)
{
  int32_T ix;
  int32_T iy;
  real_T temp;
  int32_T k;
  ix = ix0 - 1;
  iy = iy0 - 1;
  for (k = 1; k <= n; k++) {
    temp = x_data[ix];
    x_data[ix] = x_data[iy];
    x_data[iy] = temp;
    ix += incx;
    iy += incy;
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_xgetrf(int32_T m, int32_T n, real_T A_data[], int32_T
  A_size[2], int32_T lda, int32_T ipiv_data[], int32_T ipiv_size[2], int32_T
  *info)
{
  int32_T mmj;
  int32_T j;
  int32_T b_c;
  int32_T idxmax;
  int32_T ix;
  real_T smax;
  real_T s;
  int32_T jy;
  int32_T b_ix;
  int32_T b_j;
  int32_T d;
  int32_T ijA;
  int32_T y;
  if (m < n) {
    y = m;
  } else {
    y = n;
  }

  Human__eml_signed_integer_colon(y, ipiv_data, ipiv_size);
  *info = 0;
  y = m - 1;
  if (!(y < n)) {
    y = n;
  }

  for (j = 1; j <= y; j++) {
    mmj = m - j;
    b_c = (j - 1) * (lda + 1);
    if (mmj + 1 < 1) {
      idxmax = -1;
    } else {
      idxmax = 0;
      if (mmj + 1 > 1) {
        ix = b_c;
        smax = fabs(A_data[b_c]);
        for (jy = 1; jy + 1 <= mmj + 1; jy++) {
          ix++;
          s = fabs(A_data[ix]);
          if (s > smax) {
            idxmax = jy;
            smax = s;
          }
        }
      }
    }

    if (A_data[b_c + idxmax] != 0.0) {
      if (idxmax != 0) {
        ipiv_data[j - 1] = j + idxmax;
        Human_in_Loop_xswap(n, A_data, j, lda, j + idxmax, lda);
      }

      idxmax = (b_c + mmj) + 1;
      for (ix = b_c + 1; ix + 1 <= idxmax; ix++) {
        A_data[ix] /= A_data[b_c];
      }
    } else {
      *info = j;
    }

    idxmax = n - j;
    ix = (b_c + lda) + 1;
    jy = b_c + lda;
    for (b_j = 1; b_j <= idxmax; b_j++) {
      smax = A_data[jy];
      if (A_data[jy] != 0.0) {
        b_ix = b_c + 1;
        d = mmj + ix;
        for (ijA = ix; ijA + 1 <= d; ijA++) {
          A_data[ijA] += A_data[b_ix] * -smax;
          b_ix++;
        }
      }

      jy += lda;
      ix += lda;
    }
  }

  if ((*info == 0) && (m <= n) && (!(A_data[((m - 1) * A_size[0] + m) - 1] !=
        0.0))) {
    *info = m;
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_xtrsm(int32_T n, const real_T A_data[], int32_T lda,
  real_T B_data[])
{
  int32_T jAcol;
  int32_T j;
  int32_T k;
  for (j = 0; j + 1 <= n; j++) {
    jAcol = lda * j;
    for (k = 0; k + 1 <= j; k++) {
      if (A_data[k + jAcol] != 0.0) {
        B_data[j] -= A_data[k + jAcol] * B_data[k];
      }
    }

    B_data[j] *= 1.0 / A_data[j + jAcol];
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_xtrsm_e(int32_T n, const real_T A_data[], int32_T lda,
  real_T B_data[])
{
  int32_T jAcol;
  int32_T j;
  int32_T k;
  for (j = n; j > 0; j--) {
    jAcol = (j - 1) * lda;
    for (k = j; k + 1 <= n; k++) {
      if (A_data[k + jAcol] != 0.0) {
        B_data[j - 1] -= A_data[k + jAcol] * B_data[k];
      }
    }
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_lusolve(const real_T A_data[], const int32_T A_size[2],
  const real_T B_data[], const int32_T B_size[2], real_T X_data[], int32_T
  X_size[2])
{
  int32_T ipiv_data[29];
  real_T temp;
  int32_T loop_ub;
  int32_T b_A_size[2];
  int32_T ipiv_size[2];
  b_A_size[0] = A_size[0];
  b_A_size[1] = A_size[1];
  loop_ub = A_size[0] * A_size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&Human_in_Loop_B.b_A_data[0], &A_data[0], loop_ub * sizeof(real_T));
  }

  Human_in_Loop_xgetrf(A_size[1], A_size[1], Human_in_Loop_B.b_A_data, b_A_size,
                       A_size[1], ipiv_data, ipiv_size, &loop_ub);
  X_size[0] = 1;
  X_size[1] = B_size[1];
  loop_ub = B_size[0] * B_size[1];
  if (0 <= loop_ub - 1) {
    memcpy(&X_data[0], &B_data[0], loop_ub * sizeof(real_T));
  }

  Human_in_Loop_xtrsm(A_size[1], Human_in_Loop_B.b_A_data, A_size[1], X_data);
  Human_in_Loop_xtrsm_e(A_size[1], Human_in_Loop_B.b_A_data, A_size[1], X_data);
  for (loop_ub = A_size[1] - 2; loop_ub + 1 > 0; loop_ub--) {
    if (loop_ub + 1 != ipiv_data[loop_ub]) {
      temp = X_data[X_size[0] * loop_ub];
      X_data[X_size[0] * loop_ub] = X_data[(ipiv_data[loop_ub] - 1) * X_size[0]];
      X_data[X_size[0] * (ipiv_data[loop_ub] - 1)] = temp;
    }
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static void Human_in_Loop_mrdivide(const real_T A_data[], const int32_T A_size[2],
  const real_T B_data[], const int32_T B_size[2], real_T y_data[], int32_T
  y_size[2])
{
  real_T b_data[29];
  real_T A_data_0[29];
  int32_T i;
  int32_T i_0;
  int32_T loop_ub;
  int32_T b_size;
  int32_T B_size_0[2];
  if (B_size[0] == B_size[1]) {
    Human_in_Loop_lusolve(B_data, B_size, A_data, A_size, y_data, y_size);
  } else {
    B_size_0[0] = B_size[1];
    B_size_0[1] = B_size[0];
    b_size = B_size[0];
    for (i_0 = 0; i_0 < b_size; i_0++) {
      loop_ub = B_size[1];
      for (i = 0; i < loop_ub; i++) {
        Human_in_Loop_B.B_data[i + B_size_0[0] * i_0] = B_data[B_size[0] * i +
          i_0];
      }
    }

    loop_ub = A_size[1];
    b_size = A_size[1];
    for (i_0 = 0; i_0 < b_size; i_0++) {
      A_data_0[i_0] = A_data[A_size[0] * i_0];
    }

    Human_in_Loop_qrsolve(Human_in_Loop_B.B_data, B_size_0, A_data_0, &loop_ub,
                          b_data, &b_size);
    y_size[0] = 1;
    y_size[1] = b_size;
    if (0 <= b_size - 1) {
      memcpy(&y_data[0], &b_data[0], b_size * sizeof(real_T));
    }
  }
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static real_T Human_in_Loop___anon_fcn(const coder_internal_ref_Human_in_L_T
  *k_environment_f1, const real_T k_environment_f2[60], real_T k_environment_f3,
  const real_T K[900], real_T iter, const real_T Y[30], const real_T x[2])
{
  real_T varargout_1;
  real_T b_varargout_1_data[29];
  real_T a_data[29];
  real_T b_data[29];
  real_T k_environment_f2_data[58];
  int32_T i;
  int32_T loop_ub;
  int32_T i_0;
  int32_T loop_ub_0;
  int32_T b_varargout_1_size[2];
  int32_T K_size[2];
  int32_T k_environment_f2_size[2];

  /* '<S36>:1:40' */
  /* '<S36>:1:39' */
  k_environment_f2_size[0] = 2;
  k_environment_f2_size[1] = (int32_T)k_environment_f3;
  loop_ub = (int32_T)k_environment_f3;
  for (i = 0; i < loop_ub; i++) {
    k_environment_f2_data[i << 1] = k_environment_f2[i << 1];
    k_environment_f2_data[1 + (i << 1)] = k_environment_f2[(i << 1) + 1];
  }

  Human_in_Loop_KKernel(k_environment_f1, x, k_environment_f2_data,
                        k_environment_f2_size, b_varargout_1_data,
                        b_varargout_1_size);
  Human_in_Loop_eye(iter, Human_in_Loop_B.tmp_data_c, k_environment_f2_size);
  K_size[0] = (int32_T)iter;
  K_size[1] = (int32_T)iter;
  loop_ub = (int32_T)iter;
  for (i = 0; i < loop_ub; i++) {
    loop_ub_0 = (int32_T)iter;
    for (i_0 = 0; i_0 < loop_ub_0; i_0++) {
      Human_in_Loop_B.K_data[i_0 + K_size[0] * i] =
        Human_in_Loop_B.tmp_data_c[k_environment_f2_size[0] * i + i_0] * 2.5E-5
        + K[30 * i + i_0];
    }
  }

  Human_in_Loop_mrdivide(b_varargout_1_data, b_varargout_1_size,
    Human_in_Loop_B.K_data, K_size, a_data, k_environment_f2_size);
  loop_ub = (int32_T)iter;
  if (0 <= loop_ub - 1) {
    memcpy(&b_data[0], &Y[0], loop_ub * sizeof(real_T));
  }

  if ((k_environment_f2_size[1] == 1) || ((int32_T)iter == 1)) {
    varargout_1 = 0.0;
    for (i = 0; i < k_environment_f2_size[1]; i++) {
      varargout_1 += a_data[k_environment_f2_size[0] * i] * b_data[i];
    }
  } else {
    varargout_1 = 0.0;
    for (i = 0; i < k_environment_f2_size[1]; i++) {
      varargout_1 += a_data[k_environment_f2_size[0] * i] * b_data[i];
    }
  }

  return varargout_1;
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static real_T Human_in_Loop___anon_fcn_nf(const coder_internal_ref_Human_in_L_T *
  b_Theta, const coder_internal_ref_Human_in_L_T *k_environment_f1, const real_T
  k_environment_f2[60], real_T k_environment_f3, const real_T K[900], real_T
  iter, const real_T x[2])
{
  real_T varargout_1;
  real_T b_varargout_1_data[29];
  real_T c_varargout_1_data[29];
  real_T y;
  real_T a_data[29];
  real_T b_data[29];
  real_T k_environment_f2_data[58];
  int32_T i;
  int32_T loop_ub;
  int32_T i_0;
  int32_T loop_ub_0;
  int32_T b_varargout_1_size[2];
  int32_T K_size[2];
  int32_T k_environment_f2_size[2];
  int32_T k_environment_f2_size_0[2];

  /* '<S36>:1:41' */
  /* '<S36>:1:39' */
  k_environment_f2_size_0[0] = 2;
  k_environment_f2_size_0[1] = (int32_T)k_environment_f3;
  loop_ub = (int32_T)k_environment_f3;
  for (i = 0; i < loop_ub; i++) {
    k_environment_f2_data[i << 1] = k_environment_f2[i << 1];
    k_environment_f2_data[1 + (i << 1)] = k_environment_f2[(i << 1) + 1];
  }

  Human_in_Loop_KKernel(k_environment_f1, x, k_environment_f2_data,
                        k_environment_f2_size_0, b_varargout_1_data,
                        b_varargout_1_size);

  /* '<S36>:1:39' */
  k_environment_f2_size[0] = 2;
  k_environment_f2_size[1] = (int32_T)k_environment_f3;
  loop_ub = (int32_T)k_environment_f3;
  for (i = 0; i < loop_ub; i++) {
    k_environment_f2_data[i << 1] = k_environment_f2[i << 1];
    k_environment_f2_data[1 + (i << 1)] = k_environment_f2[(i << 1) + 1];
  }

  Human_in_Loop_KKernel(k_environment_f1, x, k_environment_f2_data,
                        k_environment_f2_size, c_varargout_1_data,
                        k_environment_f2_size_0);
  Human_in_Loop_eye(iter, Human_in_Loop_B.tmp_data_p, k_environment_f2_size);
  K_size[0] = (int32_T)iter;
  K_size[1] = (int32_T)iter;
  loop_ub = (int32_T)iter;
  for (i = 0; i < loop_ub; i++) {
    loop_ub_0 = (int32_T)iter;
    for (i_0 = 0; i_0 < loop_ub_0; i_0++) {
      Human_in_Loop_B.K_data_c[i_0 + K_size[0] * i] =
        Human_in_Loop_B.tmp_data_p[k_environment_f2_size[0] * i + i_0] * 2.5E-5
        + K[30 * i + i_0];
    }
  }

  Human_in_Loop_mrdivide(b_varargout_1_data, b_varargout_1_size,
    Human_in_Loop_B.K_data_c, K_size, a_data, k_environment_f2_size);
  loop_ub_0 = k_environment_f2_size_0[1];
  loop_ub = k_environment_f2_size_0[1];
  for (i = 0; i < loop_ub; i++) {
    b_data[i] = c_varargout_1_data[k_environment_f2_size_0[0] * i];
  }

  if ((k_environment_f2_size[1] == 1) || (loop_ub_0 == 1)) {
    y = 0.0;
    for (i = 0; i < k_environment_f2_size[1]; i++) {
      y += a_data[k_environment_f2_size[0] * i] * b_data[i];
    }
  } else {
    y = 0.0;
    for (i = 0; i < k_environment_f2_size[1]; i++) {
      y += a_data[k_environment_f2_size[0] * i] * b_data[i];
    }
  }

  y = Human_in_Loop_kernel(b_Theta, x, x) - y;
  varargout_1 = sqrt(y);
  return varargout_1;
}

/* Function for MATLAB Function: '<S34>/Bayesian Function' */
static real_T Human_in_Loop___anon_fcn_n(const coder_internal_ref_Human_in_L_T
  *GP_mu_environment_f1_environmen, const real_T
  GP_mu_environment_f1_environm_0[60], real_T GP_mu_environment_f1_environm_1,
  const real_T GP_mu_environment_f2[900], real_T GP_mu_environment_f3, const
  real_T GP_mu_environment_f4[30], const coder_internal_ref_Human_in_L_T
  *GP_sigma_environment_f1, const coder_internal_ref_Human_in_L_T
  *GP_sigma_environment_f2_environ, const real_T
  GP_sigma_environment_f2_envir_0[60], real_T GP_sigma_environment_f2_envir_1,
  const real_T GP_sigma_environment_f3[900], real_T GP_sigma_environment_f4,
  const real_T x[2])
{
  real_T varargout_1;
  real_T b_varargout_1;
  real_T c_varargout_1;

  /* '<S36>:1:42' */
  b_varargout_1 = Human_in_Loop___anon_fcn(GP_mu_environment_f1_environmen,
    GP_mu_environment_f1_environm_0, GP_mu_environment_f1_environm_1,
    GP_mu_environment_f2, GP_mu_environment_f3, GP_mu_environment_f4, x);
  c_varargout_1 = Human_in_Loop___anon_fcn_nf(GP_sigma_environment_f1,
    GP_sigma_environment_f2_environ, GP_sigma_environment_f2_envir_0,
    GP_sigma_environment_f2_envir_1, GP_sigma_environment_f3,
    GP_sigma_environment_f4, x);
  varargout_1 = b_varargout_1 - 3.0 * c_varargout_1;
  return varargout_1;
}

/* System initialize for function-call system: '<S25>/Bayeisan Opt' */
void Human_in_Loop_BayeisanOpt_Init(void)
{
  uint32_T r;
  int32_T mti;

  /* SystemInitialize for MATLAB Function: '<S34>/Bayesian Function' */
  Human_in_Loop_DW.iter_not_empty = false;
  memset(&Human_in_Loop_DW.state[0], 0, 625U * sizeof(uint32_T));
  r = 5489U;
  Human_in_Loop_DW.state[0] = 5489U;
  for (mti = 0; mti < 623; mti++) {
    r = ((r >> 30U ^ r) * 1812433253U + mti) + 1U;
    Human_in_Loop_DW.state[mti + 1] = r;
  }

  Human_in_Loop_DW.state[624] = 624U;
  Human_in_Loop_DW.iter = 0.0;
  memset(&Human_in_Loop_DW.X[0], 0, 60U * sizeof(real_T));
  memset(&Human_in_Loop_DW.Y[0], 0, 30U * sizeof(real_T));
  memset(&Human_in_Loop_DW.K[0], 0, 900U * sizeof(real_T));
  memset(&Human_in_Loop_DW.Mu[0], 0, 10000U * sizeof(real_T));
  memset(&Human_in_Loop_DW.LCB[0], 0, 10000U * sizeof(real_T));

  /* End of SystemInitialize for MATLAB Function: '<S34>/Bayesian Function' */

  /* SystemInitialize for Outport: '<S34>/parm' */
  Human_in_Loop_B.parm_acq[0] = Human_in_Loop_P.parm_Y0;
  Human_in_Loop_B.parm_acq[1] = Human_in_Loop_P.parm_Y0;
}

/* System reset for function-call system: '<S25>/Bayeisan Opt' */
void Human_in_Loop_BayeisanOpt_Reset(void)
{
  uint32_T r;
  int32_T mti;

  /* SystemReset for MATLAB Function: '<S34>/Bayesian Function' */
  Human_in_Loop_DW.iter_not_empty = false;
  memset(&Human_in_Loop_DW.state[0], 0, 625U * sizeof(uint32_T));
  r = 5489U;
  Human_in_Loop_DW.state[0] = 5489U;
  for (mti = 0; mti < 623; mti++) {
    r = ((r >> 30U ^ r) * 1812433253U + mti) + 1U;
    Human_in_Loop_DW.state[mti + 1] = r;
  }

  Human_in_Loop_DW.state[624] = 624U;
  Human_in_Loop_DW.iter = 0.0;
  memset(&Human_in_Loop_DW.X[0], 0, 60U * sizeof(real_T));
  memset(&Human_in_Loop_DW.Y[0], 0, 30U * sizeof(real_T));
  memset(&Human_in_Loop_DW.K[0], 0, 900U * sizeof(real_T));
  memset(&Human_in_Loop_DW.Mu[0], 0, 10000U * sizeof(real_T));
  memset(&Human_in_Loop_DW.LCB[0], 0, 10000U * sizeof(real_T));

  /* End of SystemReset for MATLAB Function: '<S34>/Bayesian Function' */
}

/* Output and update for function-call system: '<S25>/Bayeisan Opt' */
void Human_in_Loop_BayeisanOpt(void)
{
  real_T loss;
  boolean_T b_BT_OPT_RESET;
  coder_internal_ref_Human_in_L_T b_Theta;
  coder_internal_ref_Human_in_L_T *k_environment_f1;
  real_T k_environment_f2[60];
  real_T k_environment_f3;
  coder_internal_ref_Human_in_L_T *GP_mu_environment_f1_environmen;
  real_T GP_mu_environment_f1_environm_0[60];
  real_T GP_mu_environment_f3;
  real_T GP_mu_environment_f4[30];
  coder_internal_ref_Human_in_L_T *GP_sigma_environment_f1;
  coder_internal_ref_Human_in_L_T *GP_sigma_environment_f2_environ;
  real_T GP_sigma_environment_f2_envir_0[60];
  real_T GP_sigma_environment_f2_envir_1;
  real_T GP_sigma_environment_f4;
  real_T GP_LCB_environment_f1_environme;
  real_T GP_LCB_environment_f1_environ_0[30];
  real_T x[2];
  real_T kk_data[29];
  real_T extremum[100];
  int32_T iindx[100];
  int32_T ix;
  int32_T cindx;
  int32_T b_ix;
  int32_T b_ixstart;
  real_T varargin_2_data[29];
  boolean_T empty_non_axis_sizes;
  real_T b_varargin_2_data[30];
  real_T varargout_1;
  real_T tmp_data[58];
  int32_T loop_ub;
  int32_T kk_size[2];
  int32_T tmp_size[2];
  int32_T tmp;
  int32_T varargin_2_size_idx_0;
  coder_internal_ref_Human_in_L_T *tmp_0;
  coder_internal_ref_Human_in_L_T *tmp_1;
  coder_internal_ref_Human_in_L_T *tmp_2;
  boolean_T exitg1;

  /* MATLAB Function: '<S34>/Bayesian Function' */
  loss = Human_in_Loop_B.RT2_o;
  b_BT_OPT_RESET = Human_in_Loop_P.BayesianFunction_BT_OPT_RESET;

  /* MATLAB Function 'Opt Module/Subsystem/Bayeisan Opt/Bayesian Function': '<S36>:1' */
  /* '<S36>:1:42' */
  /* '<S36>:1:60' */
  b_Theta.contents = Human_in_Loop_P.BayesianFunction_Theta;
  if (!Human_in_Loop_DW.iter_not_empty) {
    /* '<S36>:1:18' */
    Human_in_Loop_DW.iter_not_empty = true;

    /* '<S36>:1:25' */
    Human_in_Loop_rand(Human_in_Loop_DW.x_opt);
  }

  if (Human_in_Loop_DW.iter == 0.0) {
    /* '<S36>:1:28' */
    /* '<S36>:1:29' */
    Human_in_Loop_rand(Human_in_Loop_B.parm_acq);

    /* '<S36>:1:30' */
    Human_in_Loop_B.parm_opt[0] = Human_in_Loop_B.parm_acq[0];
    Human_in_Loop_B.parm_opt[1] = Human_in_Loop_B.parm_acq[1];

    /* '<S36>:1:32' */
    Human_in_Loop_DW.X[0] = Human_in_Loop_B.parm_acq[0];
    Human_in_Loop_DW.X[1] = Human_in_Loop_B.parm_acq[1];

    /* '<S36>:1:33' */
    Human_in_Loop_DW.K[0] = Human_in_Loop_kernel(&b_Theta,
      Human_in_Loop_B.parm_acq, Human_in_Loop_B.parm_acq);

    /* '<S36>:1:34' */
    Human_in_Loop_DW.iter = 1.0;
  } else if (Human_in_Loop_DW.iter < 30.0) {
    /* '<S36>:1:36' */
    /* '<S36>:1:37' */
    Human_in_Loop_DW.Y[(int32_T)Human_in_Loop_DW.iter - 1] = loss;

    /* '<S36>:1:39' */
    nested_function_nested_function(&b_Theta, Human_in_Loop_DW.X,
      Human_in_Loop_DW.iter, &k_environment_f1, k_environment_f2,
      &k_environment_f3);

    /* '<S36>:1:40' */
    nested_function_nested_functi_f(k_environment_f1, k_environment_f2,
      k_environment_f3, Human_in_Loop_DW.K, Human_in_Loop_DW.iter,
      Human_in_Loop_DW.Y, &tmp_0, GP_mu_environment_f1_environm_0, &loss,
      Human_in_Loop_B.GP_mu_environment_f2, &GP_mu_environment_f3,
      GP_mu_environment_f4);
    GP_mu_environment_f1_environmen = tmp_0;

    /* '<S36>:1:41' */
    nested_function_nested_funct_fs(&b_Theta, k_environment_f1, k_environment_f2,
      k_environment_f3, Human_in_Loop_DW.K, Human_in_Loop_DW.iter, &tmp_1,
      &tmp_2, GP_sigma_environment_f2_envir_0, &GP_sigma_environment_f2_envir_1,
      Human_in_Loop_B.GP_sigma_environment_f3, &GP_sigma_environment_f4);
    GP_sigma_environment_f2_environ = tmp_2;
    GP_sigma_environment_f1 = tmp_1;

    /* '<S36>:1:42' */
    nested_function_nested_func_fst(GP_mu_environment_f1_environmen,
      GP_mu_environment_f1_environm_0, loss,
      Human_in_Loop_B.GP_mu_environment_f2, GP_mu_environment_f3,
      GP_mu_environment_f4, GP_sigma_environment_f1,
      GP_sigma_environment_f2_environ, GP_sigma_environment_f2_envir_0,
      GP_sigma_environment_f2_envir_1, Human_in_Loop_B.GP_sigma_environment_f3,
      GP_sigma_environment_f4, &Human_in_Loop_B.GP_LCB_environment_f1_environme,
      &Human_in_Loop_B.GP_LCB_environment_f2_environme);
    k_environment_f1 = Human_in_Loop_B.GP_LCB_environment_f2_environme.f1;
    GP_sigma_environment_f1 =
      Human_in_Loop_B.GP_LCB_environment_f2_environme.f2.environment.f1;
    memcpy(&k_environment_f2[0],
           &Human_in_Loop_B.GP_LCB_environment_f2_environme.f2.environment.f2[0],
           60U * sizeof(real_T));
    k_environment_f3 =
      Human_in_Loop_B.GP_LCB_environment_f2_environme.f2.environment.f3;
    memcpy(&Human_in_Loop_B.GP_sigma_environment_f3[0],
           &Human_in_Loop_B.GP_LCB_environment_f2_environme.f3[0], 900U * sizeof
           (real_T));
    GP_sigma_environment_f2_envir_1 =
      Human_in_Loop_B.GP_LCB_environment_f2_environme.f4;
    GP_sigma_environment_f2_environ =
      Human_in_Loop_B.GP_LCB_environment_f1_environme.f1.environment.f1;
    memcpy(&GP_sigma_environment_f2_envir_0[0],
           &Human_in_Loop_B.GP_LCB_environment_f1_environme.f1.environment.f2[0],
           60U * sizeof(real_T));
    GP_sigma_environment_f4 =
      Human_in_Loop_B.GP_LCB_environment_f1_environme.f1.environment.f3;
    memcpy(&Human_in_Loop_B.GP_LCB_environment_f1_environ_k[0],
           &Human_in_Loop_B.GP_LCB_environment_f1_environme.f2[0], 900U * sizeof
           (real_T));
    GP_LCB_environment_f1_environme =
      Human_in_Loop_B.GP_LCB_environment_f1_environme.f3;
    memcpy(&GP_LCB_environment_f1_environ_0[0],
           &Human_in_Loop_B.GP_LCB_environment_f1_environme.f4[0], 30U * sizeof
           (real_T));

    /* '<S36>:1:46' */
    /* '<S36>:1:44' */
    for (b_ixstart = 0; b_ixstart < 100; b_ixstart++) {
      /* '<S36>:1:44' */
      /* '<S36>:1:46' */
      /* '<S36>:1:45' */
      x[0] = (1.0 + (real_T)b_ixstart) / 100.0;
      for (ix = 0; ix < 100; ix++) {
        /* '<S36>:1:45' */
        /* '<S36>:1:46' */
        x[1] = (1.0 + (real_T)ix) / 100.0;

        /* '<S36>:1:47' */
        varargout_1 = Human_in_Loop___anon_fcn(GP_mu_environment_f1_environmen,
          GP_mu_environment_f1_environm_0, loss,
          Human_in_Loop_B.GP_mu_environment_f2, GP_mu_environment_f3,
          GP_mu_environment_f4, x);
        Human_in_Loop_DW.Mu[b_ixstart + 100 * ix] = varargout_1;

        /* '<S36>:1:48' */
        varargout_1 = Human_in_Loop___anon_fcn_n(GP_sigma_environment_f2_environ,
          GP_sigma_environment_f2_envir_0, GP_sigma_environment_f4,
          Human_in_Loop_B.GP_LCB_environment_f1_environ_k,
          GP_LCB_environment_f1_environme, GP_LCB_environment_f1_environ_0,
          k_environment_f1, GP_sigma_environment_f1, k_environment_f2,
          k_environment_f3, Human_in_Loop_B.GP_sigma_environment_f3,
          GP_sigma_environment_f2_envir_1, x);
        Human_in_Loop_DW.LCB[b_ixstart + 100 * ix] = varargout_1;
      }
    }

    /* '<S36>:1:52' */
    for (b_ixstart = 0; b_ixstart < 100; b_ixstart++) {
      ix = b_ixstart * 100;
      varargin_2_size_idx_0 = b_ixstart * 100 + 1;
      loss = Human_in_Loop_DW.LCB[ix];
      tmp = 1;
      cindx = 1;
      if (rtIsNaN(Human_in_Loop_DW.LCB[ix])) {
        b_ix = varargin_2_size_idx_0;
        exitg1 = false;
        while ((!exitg1) && (b_ix + 1 <= ix + 100)) {
          cindx++;
          varargin_2_size_idx_0 = b_ix + 1;
          if (!rtIsNaN(Human_in_Loop_DW.LCB[b_ix])) {
            loss = Human_in_Loop_DW.LCB[b_ix];
            tmp = cindx;
            exitg1 = true;
          } else {
            b_ix++;
          }
        }
      }

      if (varargin_2_size_idx_0 < ix + 100) {
        while (varargin_2_size_idx_0 + 1 <= ix + 100) {
          cindx++;
          if (Human_in_Loop_DW.LCB[varargin_2_size_idx_0] < loss) {
            loss = Human_in_Loop_DW.LCB[varargin_2_size_idx_0];
            tmp = cindx;
          }

          varargin_2_size_idx_0++;
        }
      }

      extremum[b_ixstart] = loss;
      iindx[b_ixstart] = tmp;
    }

    /* '<S36>:1:53' */
    b_ixstart = 1;
    loss = extremum[0];
    ix = 0;
    if (rtIsNaN(extremum[0])) {
      varargin_2_size_idx_0 = 2;
      exitg1 = false;
      while ((!exitg1) && (varargin_2_size_idx_0 < 101)) {
        b_ixstart = varargin_2_size_idx_0;
        if (!rtIsNaN(extremum[varargin_2_size_idx_0 - 1])) {
          loss = extremum[varargin_2_size_idx_0 - 1];
          ix = varargin_2_size_idx_0 - 1;
          exitg1 = true;
        } else {
          varargin_2_size_idx_0++;
        }
      }
    }

    if (b_ixstart < 100) {
      while (b_ixstart + 1 < 101) {
        if (extremum[b_ixstart] < loss) {
          loss = extremum[b_ixstart];
          ix = b_ixstart;
        }

        b_ixstart++;
      }
    }

    /* '<S36>:1:54' */
    Human_in_Loop_B.parm_acq[0] = (real_T)iindx[ix] / 100.0;
    Human_in_Loop_B.parm_acq[1] = (real_T)(ix + 1) / 100.0;

    /* '<S36>:1:56' */
    for (b_ixstart = 0; b_ixstart < 100; b_ixstart++) {
      ix = b_ixstart * 100;
      varargin_2_size_idx_0 = b_ixstart * 100 + 1;
      loss = Human_in_Loop_DW.Mu[ix];
      tmp = 1;
      cindx = 1;
      if (rtIsNaN(Human_in_Loop_DW.Mu[ix])) {
        b_ix = varargin_2_size_idx_0;
        exitg1 = false;
        while ((!exitg1) && (b_ix + 1 <= ix + 100)) {
          cindx++;
          varargin_2_size_idx_0 = b_ix + 1;
          if (!rtIsNaN(Human_in_Loop_DW.Mu[b_ix])) {
            loss = Human_in_Loop_DW.Mu[b_ix];
            tmp = cindx;
            exitg1 = true;
          } else {
            b_ix++;
          }
        }
      }

      if (varargin_2_size_idx_0 < ix + 100) {
        while (varargin_2_size_idx_0 + 1 <= ix + 100) {
          cindx++;
          if (Human_in_Loop_DW.Mu[varargin_2_size_idx_0] < loss) {
            loss = Human_in_Loop_DW.Mu[varargin_2_size_idx_0];
            tmp = cindx;
          }

          varargin_2_size_idx_0++;
        }
      }

      extremum[b_ixstart] = loss;
      iindx[b_ixstart] = tmp;
    }

    /* '<S36>:1:57' */
    b_ixstart = 1;
    loss = extremum[0];
    ix = 0;
    if (rtIsNaN(extremum[0])) {
      varargin_2_size_idx_0 = 2;
      exitg1 = false;
      while ((!exitg1) && (varargin_2_size_idx_0 < 101)) {
        b_ixstart = varargin_2_size_idx_0;
        if (!rtIsNaN(extremum[varargin_2_size_idx_0 - 1])) {
          loss = extremum[varargin_2_size_idx_0 - 1];
          ix = varargin_2_size_idx_0 - 1;
          exitg1 = true;
        } else {
          varargin_2_size_idx_0++;
        }
      }
    }

    if (b_ixstart < 100) {
      while (b_ixstart + 1 < 101) {
        if (extremum[b_ixstart] < loss) {
          loss = extremum[b_ixstart];
          ix = b_ixstart;
        }

        b_ixstart++;
      }
    }

    /* '<S36>:1:58' */
    Human_in_Loop_B.parm_opt[0] = (real_T)iindx[ix] / 100.0;
    Human_in_Loop_B.parm_opt[1] = (real_T)(ix + 1) / 100.0;

    /* '<S36>:1:60' */
    tmp_size[0] = 2;
    tmp_size[1] = (int32_T)Human_in_Loop_DW.iter;
    b_ixstart = (int32_T)Human_in_Loop_DW.iter;
    for (ix = 0; ix < b_ixstart; ix++) {
      tmp_data[ix << 1] = Human_in_Loop_DW.X[ix << 1];
      tmp_data[1 + (ix << 1)] = Human_in_Loop_DW.X[(ix << 1) + 1];
    }

    Human_in_Loop_KKernel(&b_Theta, Human_in_Loop_B.parm_acq, tmp_data, tmp_size,
                          kk_data, kk_size);
    varargin_2_size_idx_0 = kk_size[1];
    b_ixstart = kk_size[1];
    for (ix = 0; ix < b_ixstart; ix++) {
      varargin_2_data[ix] = kk_data[kk_size[0] * ix];
    }

    tmp = (int32_T)Human_in_Loop_DW.iter;
    b_ixstart = (int32_T)Human_in_Loop_DW.iter;
    if (!(tmp == 0)) {
      varargin_2_size_idx_0 = b_ixstart;
    }

    empty_non_axis_sizes = (varargin_2_size_idx_0 == 0);
    if (empty_non_axis_sizes) {
      tmp = (int32_T)Human_in_Loop_DW.iter;
    } else {
      tmp = (int32_T)Human_in_Loop_DW.iter;
      if (!(tmp == 0)) {
        tmp = (int32_T)Human_in_Loop_DW.iter;
      } else {
        tmp = 0;
      }
    }

    b_ix = (int32_T)Human_in_Loop_DW.iter;
    b_ixstart = (int32_T)Human_in_Loop_DW.iter;
    for (ix = 0; ix < b_ixstart; ix++) {
      loop_ub = (int32_T)Human_in_Loop_DW.iter;
      for (cindx = 0; cindx < loop_ub; cindx++) {
        Human_in_Loop_B.tmp_data[cindx + b_ix * ix] = Human_in_Loop_DW.K[30 * ix
          + cindx];
      }
    }

    b_ix = tmp + 1;
    for (ix = 0; ix < tmp; ix++) {
      for (cindx = 0; cindx < varargin_2_size_idx_0; cindx++) {
        Human_in_Loop_B.result_data[cindx + varargin_2_size_idx_0 * ix] =
          Human_in_Loop_B.tmp_data[varargin_2_size_idx_0 * ix + cindx];
      }
    }

    for (ix = 0; ix < varargin_2_size_idx_0; ix++) {
      Human_in_Loop_B.result_data[ix + varargin_2_size_idx_0 * tmp] =
        varargin_2_data[ix];
    }

    loss = Human_in_Loop_kernel(&b_Theta, Human_in_Loop_B.parm_acq,
      Human_in_Loop_B.parm_acq);
    b_ixstart = kk_size[0] * kk_size[1];
    if (0 <= b_ixstart - 1) {
      memcpy(&b_varargin_2_data[0], &kk_data[0], b_ixstart * sizeof(real_T));
    }

    b_varargin_2_data[kk_size[0] * kk_size[1]] = loss;

    /* '<S36>:1:61' */
    for (ix = 0; ix < b_ix; ix++) {
      for (cindx = 0; cindx < varargin_2_size_idx_0; cindx++) {
        Human_in_Loop_DW.K[cindx + 30 * ix] =
          Human_in_Loop_B.result_data[varargin_2_size_idx_0 * ix + cindx];
      }
    }

    for (ix = 0; ix < b_ix; ix++) {
      for (cindx = 0; cindx < 1; cindx++) {
        Human_in_Loop_DW.K[varargin_2_size_idx_0 + 30 * ix] =
          b_varargin_2_data[ix];
      }
    }

    /* '<S36>:1:62' */
    ix = (int32_T)(Human_in_Loop_DW.iter + 1.0);
    Human_in_Loop_DW.X[(ix - 1) << 1] = Human_in_Loop_B.parm_acq[0];
    Human_in_Loop_DW.X[1 + ((ix - 1) << 1)] = Human_in_Loop_B.parm_acq[1];

    /* '<S36>:1:63' */
    Human_in_Loop_DW.x_opt[0] = Human_in_Loop_B.parm_opt[0];
    Human_in_Loop_DW.x_opt[1] = Human_in_Loop_B.parm_opt[1];

    /* '<S36>:1:64' */
    Human_in_Loop_DW.iter++;
  } else {
    /* '<S36>:1:67' */
    Human_in_Loop_B.parm_acq[0] = Human_in_Loop_DW.x_opt[0];
    Human_in_Loop_B.parm_acq[1] = Human_in_Loop_DW.x_opt[1];

    /* '<S36>:1:68' */
    Human_in_Loop_B.parm_opt[0] = Human_in_Loop_DW.x_opt[0];
    Human_in_Loop_B.parm_opt[1] = Human_in_Loop_DW.x_opt[1];
  }

  /* '<S36>:1:71' */
  memcpy(&Human_in_Loop_B.target_esitmate[0], &Human_in_Loop_DW.Mu[0], 10000U *
         sizeof(real_T));

  /* '<S36>:1:72' */
  loss = Human_in_Loop_DW.iter;
  if (b_BT_OPT_RESET) {
    /* '<S36>:1:74' */
    /* '<S36>:1:75' */
    Human_in_Loop_DW.iter = 0.0;

    /* '<S36>:1:76' */
    memset(&Human_in_Loop_DW.X[0], 0, 60U * sizeof(real_T));

    /* '<S36>:1:77' */
    memset(&Human_in_Loop_DW.Y[0], 0, 30U * sizeof(real_T));

    /* '<S36>:1:78' */
    memset(&Human_in_Loop_DW.K[0], 0, 900U * sizeof(real_T));

    /* '<S36>:1:79' */
    /* '<S36>:1:80' */
    memset(&Human_in_Loop_DW.Mu[0], 0, 10000U * sizeof(real_T));
    memset(&Human_in_Loop_DW.LCB[0], 0, 10000U * sizeof(real_T));

    /* '<S36>:1:81' */
    Human_in_Loop_B.parm_acq[0] = 0.0;
    Human_in_Loop_B.parm_acq[1] = 0.0;

    /* '<S36>:1:82' */
    Human_in_Loop_B.parm_opt[0] = 0.0;
    Human_in_Loop_B.parm_opt[1] = 0.0;
  }

  Human_in_Loop_B.Iteration = loss;

  /* End of MATLAB Function: '<S34>/Bayesian Function' */
}

/*
 * Output and update for atomic system:
 *    '<S45>/Mux1'
 *    '<S45>/Mux3'
 */
void Human_in_Loop_Mux1(real_T rtu_x1, real_T rtu_x2, real_T rtu_x3, real_T
  rtu_x4, real_T rtu_x5, real_T rtu_x6, B_Mux1_Human_in_Loop_T *localB)
{
  /* MATLAB Function 'Sensor Data/EMG module/Mux1': '<S56>:1' */
  /* '<S56>:1:3' */
  localB->x[0] = rtu_x1;
  localB->x[1] = rtu_x2;
  localB->x[2] = rtu_x3;
  localB->x[3] = rtu_x4;
  localB->x[4] = rtu_x5;
  localB->x[5] = rtu_x6;
}

/*
 * System initialize for atomic system:
 *    '<S60>/MATLAB Function'
 *    '<S61>/MATLAB Function'
 *    '<S64>/MATLAB Function'
 *    '<S65>/MATLAB Function'
 *    '<S66>/MATLAB Function'
 *    '<S67>/MATLAB Function'
 */
void Human_in_Lo_MATLABFunction_Init(DW_MATLABFunction_Human_in_Lo_T *localDW)
{
  memset(&localDW->data_mem[0], 0, 500U * sizeof(real_T));
}

/*
 * Output and update for atomic system:
 *    '<S60>/MATLAB Function'
 *    '<S61>/MATLAB Function'
 *    '<S64>/MATLAB Function'
 *    '<S65>/MATLAB Function'
 *    '<S66>/MATLAB Function'
 *    '<S67>/MATLAB Function'
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

  /* MATLAB Function 'Sensor Data/EMG module/Preprocessing10/MATLAB Function': '<S70>:1' */
  /* '<S70>:1:8' */
  memcpy(&tmp[0], &localDW->data_mem[1], 499U * sizeof(real_T));
  tmp[499] = rtu_data;

  /* '<S70>:1:10' */
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
 *    '<S46>/Mux'
 *    '<S49>/Mux'
 */
void Human_in_Loop_Mux(real_T rtu_x1, real_T rtu_x2, B_Mux_Human_in_Loop_T
  *localB)
{
  /* MATLAB Function 'Sensor Data/Encoder module/Mux': '<S85>:1' */
  /* '<S85>:1:3' */
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
    memcpy(&Human_in_Loop_B.z1_data_p[0], &y_data[0], loop_ub * sizeof(real_T));
  }

  for (loop_ub = 0; loop_ub + 1 <= y_size[1]; loop_ub++) {
    Human_in_Loop_B.z1_data_p[loop_ub] = rt_powd_snf(a_data[loop_ub], 3.0);
  }

  y_size[0] = 1;
  y_size[1] = z1_size_idx_1;
  if (0 <= z1_size_idx_1 - 1) {
    memcpy(&y_data[0], &Human_in_Loop_B.z1_data_p[0], z1_size_idx_1 * sizeof
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
  /* MATLAB Function 'Control Module/Torque track': '<S13>:1' */
  /* '<S13>:1:42' */
  /* '<S13>:1:19' */
  mode = Human_in_Loop_B.RT1[0];

  /* '<S13>:1:20' */
  footstate = Human_in_Loop_B.RT1[1];

  /* '<S13>:1:21' */
  /* '<S13>:1:22' */
  /* '<S13>:1:24' */
  peak_torque = Human_in_Loop_B.RT3[0];

  /* '<S13>:1:25' */
  /* '<S13>:1:26' */
  /* '<S13>:1:27' */
  /* '<S13>:1:29' */
  torque_measure = Human_in_Loop_B.RT4[0];

  /* '<S13>:1:30' */
  troque_delta = Human_in_Loop_B.RT4[1];

  /* '<S13>:1:31' */
  stride_index = Human_in_Loop_B.RT1[3] * 500.0 + 1.0;
  if (stride_index > 750.0) {
    /* '<S13>:1:32' */
    /* '<S13>:1:33' */
    stride_index = 750.0;
  }

  if ((Human_in_Loop_DW.last_footstate == 0.0) && (Human_in_Loop_B.RT1[1] == 1.0)
      && ((Human_in_Loop_B.RT1[0] == 2.0) || (Human_in_Loop_B.RT1[0] == 1.0))) {
    /* '<S13>:1:37' */
    /* '<S13>:1:39' */
    /* '<S13>:1:40' */
    memset(&Human_in_Loop_B.torque_track[0], 0, 750U * sizeof(real_T));
    memset(&Human_in_Loop_B.torque_delta_track[0], 0, 750U * sizeof(real_T));

    /* '<S13>:1:41' */
    /* '<S13>:1:42' */
    /* '<S13>:1:43' */
    index_peak = floor(Human_in_Loop_B.RT3[2] / 100.0 * Human_in_Loop_B.RT1[2] *
                       500.0);

    /* '<S13>:1:44' */
    index_rise = index_peak - floor(Human_in_Loop_B.RT3[1] / 100.0 *
      Human_in_Loop_B.RT1[2] * 500.0);

    /* '<S13>:1:45' */
    index_fall = floor(Human_in_Loop_B.RT3[3] / 100.0 * Human_in_Loop_B.RT1[2] *
                       500.0) + index_peak;

    /* '<S13>:1:48' */
    /* '<S13>:1:52' */
    /* '<S13>:1:53' */
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
      Human_in_Loop_B.tmp_data_k[f] = tmp_data[f];
    }

    Human_in_Loop_power(Human_in_Loop_B.tmp_data_k, tmp_size_2,
                        Human_in_Loop_B.tmp_data_cx, tmp_size_3);
    tmp_size_idx_1 = h - ia;
    loop_ub = h - ia;
    for (f = 0; f < loop_ub; f++) {
      tmp_data[f] = (int16_T)((int16_T)(ia + f) + 1);
    }

    tmp_size_4[0] = 1;
    tmp_size_4[1] = tmp_size_idx_1;
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.tmp_data_k[f] = tmp_data[f];
    }

    Human_in_Loop_power_a(Human_in_Loop_B.tmp_data_k, tmp_size_4,
                          Human_in_Loop_B.tmp_data_b, tmp_size_2);
    for (f = 0; f < b_n; f++) {
      Human_in_Loop_B.tmp_data_cl[f << 2] = 1.0;
    }

    loop_ub = d - br;
    for (f = 0; f <= loop_ub; f++) {
      Human_in_Loop_B.tmp_data_cl[1 + (f << 2)] = (int16_T)((int16_T)((br + f) -
        1) + 1);
    }

    loop_ub = tmp_size_3[1];
    for (f = 0; f < loop_ub; f++) {
      Human_in_Loop_B.tmp_data_cl[2 + (f << 2)] =
        Human_in_Loop_B.tmp_data_cx[tmp_size_3[0] * f];
    }

    loop_ub = tmp_size_2[1];
    for (f = 0; f < loop_ub; f++) {
      Human_in_Loop_B.tmp_data_cl[3 + (f << 2)] =
        Human_in_Loop_B.tmp_data_b[tmp_size_2[0] * f];
    }

    br = b_n;
    for (f = 0; f < b_n; f++) {
      Human_in_Loop_B.result_data_m[f << 2] = Human_in_Loop_B.tmp_data_cl[f << 2];
      Human_in_Loop_B.result_data_m[1 + (f << 2)] = Human_in_Loop_B.tmp_data_cl
        [(f << 2) + 1];
      Human_in_Loop_B.result_data_m[2 + (f << 2)] = Human_in_Loop_B.tmp_data_cl
        [(f << 2) + 2];
      Human_in_Loop_B.result_data_m[3 + (f << 2)] = Human_in_Loop_B.tmp_data_cl
        [(f << 2) + 3];
    }

    tmp_0 = (int16_T)br;
    tmp_size_idx_1 = tmp_0;
    b_n = br - 1;
    if (0 <= tmp_size_idx_1 - 1) {
      memset(&Human_in_Loop_B.tmp_data_k[0], 0, tmp_size_idx_1 * sizeof(real_T));
    }

    if (br != 0) {
      for (br = 1; br - 1 <= b_n; br++) {
        for (d = br; d <= br; d++) {
          Human_in_Loop_B.tmp_data_k[d - 1] = 0.0;
        }
      }

      br = 0;
      for (d = 0; d <= b_n; d++) {
        ar = -1;
        for (f = br; f + 1 <= br + 4; f++) {
          if (Human_in_Loop_B.result_data_m[f] != 0.0) {
            ia = ar;
            for (h = d; h + 1 <= d + 1; h++) {
              ia++;
              Human_in_Loop_B.tmp_data_k[h] += Human_in_Loop_B.result_data_m[f] *
                parm1[ia];
            }
          }

          ar++;
        }

        br += 4;
      }
    }

    /* '<S13>:1:54' */
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.torque_track[y + f] = Human_in_Loop_B.tmp_data_k[f];
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
      Human_in_Loop_B.tmp_data_k[f] = tmp_data[f];
    }

    Human_in_Loop_power(Human_in_Loop_B.tmp_data_k, tmp_size_1,
                        Human_in_Loop_B.tmp_data_cx, tmp_size_2);
    for (f = 0; f < b_n; f++) {
      Human_in_Loop_B.tmp_data_cv[3 * f] = 1.0;
    }

    loop_ub = d - br;
    for (f = 0; f <= loop_ub; f++) {
      Human_in_Loop_B.tmp_data_cv[1 + 3 * f] = (real_T)(int16_T)((int16_T)((br +
        f) - 1) + 1) * 2.0;
    }

    loop_ub = tmp_size_2[1];
    for (f = 0; f < loop_ub; f++) {
      Human_in_Loop_B.tmp_data_cv[2 + 3 * f] =
        Human_in_Loop_B.tmp_data_cx[tmp_size_2[0] * f] * 3.0;
    }

    br = b_n;
    for (f = 0; f < b_n; f++) {
      Human_in_Loop_B.b_result_data[3 * f] = Human_in_Loop_B.tmp_data_cv[3 * f];
      Human_in_Loop_B.b_result_data[1 + 3 * f] = Human_in_Loop_B.tmp_data_cv[3 *
        f + 1];
      Human_in_Loop_B.b_result_data[2 + 3 * f] = Human_in_Loop_B.tmp_data_cv[3 *
        f + 2];
    }

    tmp_0 = (int16_T)br;
    tmp_size_idx_1 = tmp_0;
    b_n = br - 1;
    if (0 <= tmp_size_idx_1 - 1) {
      memset(&Human_in_Loop_B.tmp_data_k[0], 0, tmp_size_idx_1 * sizeof(real_T));
    }

    if (br != 0) {
      for (br = 1; br - 1 <= b_n; br++) {
        for (d = br; d <= br; d++) {
          Human_in_Loop_B.tmp_data_k[d - 1] = 0.0;
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
              Human_in_Loop_B.tmp_data_k[h] += parm1[1 + ia] *
                Human_in_Loop_B.b_result_data[f];
            }
          }

          ar++;
        }

        br += 3;
      }
    }

    /* '<S13>:1:55' */
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.torque_delta_track[y + f] = Human_in_Loop_B.tmp_data_k[f] *
        500.0;
    }

    /* '<S13>:1:57' */
    /* '<S13>:1:61' */
    /* '<S13>:1:62' */
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
      Human_in_Loop_B.tmp_data_k[f] = tmp_data[f];
    }

    Human_in_Loop_power(Human_in_Loop_B.tmp_data_k, tmp_size,
                        Human_in_Loop_B.tmp_data_cx, tmp_size_2);
    tmp_size_idx_1 = h - ia;
    loop_ub = h - ia;
    for (f = 0; f < loop_ub; f++) {
      tmp_data[f] = (int16_T)((int16_T)(ia + f) + 1);
    }

    tmp_size_0[0] = 1;
    tmp_size_0[1] = tmp_size_idx_1;
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.tmp_data_k[f] = tmp_data[f];
    }

    Human_in_Loop_power_a(Human_in_Loop_B.tmp_data_k, tmp_size_0,
                          Human_in_Loop_B.tmp_data_b, tmp_size_3);
    for (f = 0; f < b_n; f++) {
      Human_in_Loop_B.tmp_data_cl[f << 2] = 1.0;
    }

    loop_ub = d - br;
    for (f = 0; f <= loop_ub; f++) {
      Human_in_Loop_B.tmp_data_cl[1 + (f << 2)] = (int16_T)((int16_T)((br + f) -
        1) + 1);
    }

    loop_ub = tmp_size_2[1];
    for (f = 0; f < loop_ub; f++) {
      Human_in_Loop_B.tmp_data_cl[2 + (f << 2)] =
        Human_in_Loop_B.tmp_data_cx[tmp_size_2[0] * f];
    }

    loop_ub = tmp_size_3[1];
    for (f = 0; f < loop_ub; f++) {
      Human_in_Loop_B.tmp_data_cl[3 + (f << 2)] =
        Human_in_Loop_B.tmp_data_b[tmp_size_3[0] * f];
    }

    br = b_n;
    for (f = 0; f < b_n; f++) {
      Human_in_Loop_B.result_data_m[f << 2] = Human_in_Loop_B.tmp_data_cl[f << 2];
      Human_in_Loop_B.result_data_m[1 + (f << 2)] = Human_in_Loop_B.tmp_data_cl
        [(f << 2) + 1];
      Human_in_Loop_B.result_data_m[2 + (f << 2)] = Human_in_Loop_B.tmp_data_cl
        [(f << 2) + 2];
      Human_in_Loop_B.result_data_m[3 + (f << 2)] = Human_in_Loop_B.tmp_data_cl
        [(f << 2) + 3];
    }

    tmp_0 = (int16_T)br;
    tmp_size_idx_1 = tmp_0;
    b_n = br - 1;
    if (0 <= tmp_size_idx_1 - 1) {
      memset(&Human_in_Loop_B.tmp_data_k[0], 0, tmp_size_idx_1 * sizeof(real_T));
    }

    if (br != 0) {
      for (br = 1; br - 1 <= b_n; br++) {
        for (d = br; d <= br; d++) {
          Human_in_Loop_B.tmp_data_k[d - 1] = 0.0;
        }
      }

      br = 0;
      for (d = 0; d <= b_n; d++) {
        ar = -1;
        for (f = br; f + 1 <= br + 4; f++) {
          if (Human_in_Loop_B.result_data_m[f] != 0.0) {
            ia = ar;
            for (h = d; h + 1 <= d + 1; h++) {
              ia++;
              Human_in_Loop_B.tmp_data_k[h] += Human_in_Loop_B.result_data_m[f] *
                parm1[ia];
            }
          }

          ar++;
        }

        br += 4;
      }
    }

    /* '<S13>:1:63' */
    for (f = 0; f < tmp_size_idx_1; f++) {
      Human_in_Loop_B.torque_track[y + f] = Human_in_Loop_B.tmp_data_k[f];
    }

    /* '<S13>:1:67' */
    /* '<S13>:1:69' */
    for (f = 0; f < 750; f++) {
      Human_in_Loop_DW.TorqueMem[f << 2] = Human_in_Loop_B.torque_track[f];
      Human_in_Loop_DW.TorqueMem[2 + (f << 2)] =
        Human_in_Loop_B.torque_delta_track[f];
    }
  }

  /* '<S13>:1:73' */
  Human_in_Loop_DW.last_footstate = footstate;

  /* '<S13>:1:74' */
  Human_in_Loop_DW.TorqueMem[1 + (((int32_T)stride_index - 1) << 2)] =
    torque_measure;

  /* '<S13>:1:75' */
  Human_in_Loop_DW.TorqueMem[3 + (((int32_T)stride_index - 1) << 2)] =
    troque_delta;
  if (mode == 2.0) {
    /* '<S13>:1:77' */
    /* '<S13>:1:78' */
    mode = Human_in_Loop_DW.TorqueMem[((int32_T)stride_index - 1) << 2];

    /* '<S13>:1:79' */
    stride_index = Human_in_Loop_DW.TorqueMem[(((int32_T)stride_index - 1) << 2)
      + 2];
  } else {
    /* '<S13>:1:81' */
    mode = 0.0;

    /* '<S13>:1:82' */
    stride_index = 0.0;
  }

  /* '<S13>:1:85' */
  /* '<S13>:1:86' */
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
  /* MATLAB Function 'Control Module/LRN': '<S11>:1' */
  /* '<S11>:1:15' */
  /* '<S11>:1:19' */
  /* '<S11>:1:22' */
  /* '<S11>:1:23' */
  /* '<S11>:1:25' */
  /* '<S11>:1:27' */
  /* '<S11>:1:29' */
  stride_index = Human_in_Loop_B.RT1[3] * 500.0 + 1.0;
  if (stride_index > 750.0) {
    /* '<S11>:1:30' */
    /* '<S11>:1:31' */
    stride_index = 750.0;
  }

  if ((Human_in_Loop_DW.last_footstate_a == 0.0) && (Human_in_Loop_B.RT1[1] ==
       1.0) && (Human_in_Loop_B.RT1[0] == 2.0) && Human_in_Loop_P.LRN_BT_LRN_ON)
  {
    /* '<S11>:1:34' */
    /* '<S11>:1:35' */
    /* '<S11>:1:36' */
    mode = 1.0 - Human_in_Loop_P.LRN_error_filter_k;

    /* '<S11>:1:37' */
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
    /* '<S11>:1:40' */
    /* '<S11>:1:41' */
    /* '<S11>:1:42' */
    memset(&Human_in_Loop_DW.torque_error_memory[0], 0, 1000U * sizeof(real_T));
    memset(&Human_in_Loop_DW.lrn_cmd_memory[0], 0, 1000U * sizeof(real_T));
  }

  if (Human_in_Loop_B.RT1[0] == 2.0) {
    /* '<S11>:1:45' */
    b_n = Human_in_Loop_P.LRN_time_delay;
    if (b_n < -2147482897) {
      b_n = MAX_int32_T;
    } else {
      b_n = 750 - b_n;
    }

    if (stride_index >= b_n) {
      /* '<S11>:1:46' */
      /* '<S11>:1:47' */
      Human_in_Loop_B.lrn_cmd = 0.0;
    } else {
      /* '<S11>:1:49' */
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
    /* '<S11>:1:52' */
    Human_in_Loop_B.lrn_cmd = 0.0;
  }

  /* '<S11>:1:55' */
  Human_in_Loop_DW.last_footstate_a = Human_in_Loop_B.RT1[1];

  /* '<S11>:1:56' */
  memcpy(&Human_in_Loop_B.lrn_mem[0], &Human_in_Loop_DW.lrn_cmd_memory[0], 750U *
         sizeof(real_T));

  /* End of MATLAB Function: '<S1>/LRN' */

  /* MATLAB Function: '<S1>/Controller' */
  /* MATLAB Function 'Control Module/Controller': '<S10>:1' */
  /* '<S10>:1:21' */
  /* '<S10>:1:22' */
  /* '<S10>:1:26' */
  /* '<S10>:1:27' */
  /* '<S10>:1:31' */
  /* '<S10>:1:33' */
  /* '<S10>:1:36' */
  /* '<S10>:1:37' */
  /* '<S10>:1:41' */
  /* '<S10>:1:42' */
  /* '<S10>:1:43' */
  /* '<S10>:1:44' */
  switch ((int32_T)Human_in_Loop_B.RT1[0]) {
   case 1:
    /* '<S10>:1:48' */
    Human_in_Loop_B.motor_vel_cmd = 0.0;

    /* '<S10>:1:49' */
    Human_in_Loop_DW.calib_state = 0.0;
    break;

   case 3:
    /* '<S10>:1:52' */
    Human_in_Loop_B.motor_vel_cmd = -Human_in_Loop_P.Controller_SLACK_SPEED *
      5.0;
    break;

   case 4:
    if ((Human_in_Loop_DW.calib_state == 0.0) && (Human_in_Loop_B.RT4[0] <
         Human_in_Loop_P.Controller_CALIB_TORQUE)) {
      /* '<S10>:1:55' */
      /* '<S10>:1:56' */
      Human_in_Loop_B.motor_vel_cmd = Human_in_Loop_P.Controller_CALIB_SPEED *
        5.0;
    } else if ((Human_in_Loop_DW.calib_state == 0.0) && (Human_in_Loop_B.RT4[0] >
                Human_in_Loop_P.Controller_CALIB_TORQUE)) {
      /* '<S10>:1:57' */
      /* '<S10>:1:58' */
      Human_in_Loop_DW.calib_state = 1.0;

      /* '<S10>:1:59' */
      Human_in_Loop_B.motor_vel_cmd = 0.0;
    } else {
      /* '<S10>:1:61' */
      Human_in_Loop_B.motor_vel_cmd = 0.0;
    }
    break;

   case 2:
    /* '<S10>:1:65' */
    switch (Human_in_Loop_P.Controller_run_mode) {
     case 1:
      if (Human_in_Loop_B.RT5[0] > 0.0) {
        /* '<S10>:1:68' */
        /* '<S10>:1:69' */
        stride_index = (Human_in_Loop_B.RT5[0] - Human_in_Loop_B.RT5[0] / 50.0 *
                        Human_in_Loop_P.Controller_FOLLOW_SLACK_ANGLE) -
          Human_in_Loop_B.RT6[0] * 0.33333333333333331;
      } else {
        /* '<S10>:1:71' */
        stride_index = Human_in_Loop_B.RT5[0] - Human_in_Loop_B.RT6[0] *
          0.33333333333333331;
      }

      /* '<S10>:1:73' */
      Human_in_Loop_B.motor_vel_cmd = Human_in_Loop_B.RT2[3] * stride_index *
        5.0 / 0.05;
      break;

     case 2:
      if (Human_in_Loop_B.RT1[1] == 1.0) {
        /* '<S10>:1:76' */
        /* '<S10>:1:77' */
        /* '<S10>:1:78' */
        /* '<S10>:1:79' */
        /* '<S10>:1:80' */
        Human_in_Loop_B.motor_vel_cmd = ((((Human_in_Loop_B.torque_des -
          Human_in_Loop_B.RT4[0]) * Human_in_Loop_B.RT2[0] +
          (Human_in_Loop_B.torque_delta_des - Human_in_Loop_B.RT4[1]) *
          Human_in_Loop_B.RT2[1]) + Human_in_Loop_B.lrn_cmd) +
          Human_in_Loop_B.RT2[4] * Human_in_Loop_B.torque_delta_des) * 5.0 /
          0.05;
      } else {
        /* '<S10>:1:82' */
        /* '<S10>:1:83' */
        Human_in_Loop_B.motor_vel_cmd = (0.0 - Human_in_Loop_B.RT6[0] *
          0.33333333333333331) * Human_in_Loop_B.RT2[3] * 5.0 / 0.05;
      }
      break;

     default:
      /* '<S10>:1:86' */
      Human_in_Loop_B.motor_vel_cmd = 0.0;
      break;
    }
    break;

   case 0:
    /* '<S10>:1:90' */
    Human_in_Loop_B.motor_vel_cmd = 0.0;
    break;

   default:
    /* '<S10>:1:93' */
    Human_in_Loop_B.motor_vel_cmd = 0.0;
    break;
  }

  /* End of MATLAB Function: '<S1>/Controller' */

  /* MATLAB Function: '<S12>/MATLAB Function' */
  /* MATLAB Function 'Control Module/Motor/MATLAB Function': '<S15>:1' */
  /* '<S15>:1:3' */
  /* '<S15>:1:4' */
  if (Human_in_Loop_B.RT4[0] > Human_in_Loop_P.MATLABFunction_MAX_TORQUE) {
    /* '<S15>:1:6' */
    /* '<S15>:1:7' */
    Human_in_Loop_B.vel = 0.0;
  } else if (Human_in_Loop_B.RT6[0] >
             Human_in_Loop_P.MATLABFunction_MAX_MOTOR_ANGLE) {
    /* '<S15>:1:8' */
    /* '<S15>:1:9' */
    Human_in_Loop_B.vel = 0.0;
  } else if (Human_in_Loop_B.RT6[0] <
             Human_in_Loop_P.MATLABFunction_MIN_MOTOR_ANGLE) {
    /* '<S15>:1:10' */
    /* '<S15>:1:11' */
    Human_in_Loop_B.vel = 0.0;
  } else if (Human_in_Loop_B.motor_vel_cmd >
             Human_in_Loop_P.MATLABFunction_MAX_SPEED) {
    /* '<S15>:1:12' */
    /* '<S15>:1:13' */
    Human_in_Loop_B.vel = Human_in_Loop_P.MATLABFunction_MAX_SPEED;
  } else if (Human_in_Loop_B.motor_vel_cmd <
             -Human_in_Loop_P.MATLABFunction_MAX_SPEED) {
    /* '<S15>:1:14' */
    /* '<S15>:1:15' */
    Human_in_Loop_B.vel = -Human_in_Loop_P.MATLABFunction_MAX_SPEED;
  } else {
    /* '<S15>:1:17' */
    Human_in_Loop_B.vel = Human_in_Loop_B.motor_vel_cmd;
  }

  /* End of MATLAB Function: '<S12>/MATLAB Function' */

  /* Gain: '<S12>/Gain2' */
  Human_in_Loop_B.Gain2_j = Human_in_Loop_P.Gain2_Gain * Human_in_Loop_B.vel;

  /* Gain: '<S12>/Gain1' */
  Human_in_Loop_B.Gain1_pi = Human_in_Loop_P.Gain1_Gain *
    Human_in_Loop_B.Gain2_j;

  /* S-Function (rti_commonblock): '<S14>/S-Function1' */
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
  /* Terminate for S-Function (rti_commonblock): '<S14>/S-Function1' */

  /* --- Human_in_Loop/Control Module/Motor/DAC_CLASS1_BL1 --- */
  /* --- [RTI120X, DAC C1] - Channel: 16 --- */

  /* All channel outputs are set to high impedance state */
  DacCl1AnalogOut_setOutputHighZ(pRTIDacC1AnalogOut_Ch_16, DAC_CLASS1_HIGH_Z_ON);

  /* Deactivates AnalogOut functionality */
  DacCl1AnalogOut_stop(pRTIDacC1AnalogOut_Ch_16);
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static void Human_in_Loop_eye_n(real_T I_data[], int32_T I_size[2])
{
  I_size[0] = 1;
  I_size[1] = 1;
  I_data[0] = 1.0;
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static void Human_in_Loop_log(real_T x_data[], int32_T x_size[2])
{
  int32_T nx;
  int32_T k;
  nx = x_size[1];
  for (k = 0; k + 1 <= nx; k++) {
    x_data[k] = log(x_data[k]);
  }
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static real_T Human_in_Loop_sum(const real_T x_data[], const int32_T *x_size)
{
  real_T y;
  int32_T k;
  y = x_data[0];
  k = 2;
  while (k <= *x_size) {
    y += x_data[1];
    k = 3;
  }

  return y;
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static void Human_in_Loop_power_d(const real_T a_data[], const int32_T *a_size,
  real_T y_data[], int32_T *y_size)
{
  real_T z1_data[2];
  int32_T loop_ub;
  int8_T a_idx_0;
  int32_T z1_size_idx_0;
  a_idx_0 = (int8_T)*a_size;
  *y_size = a_idx_0;
  z1_size_idx_0 = *y_size;
  loop_ub = *y_size;
  if (0 <= loop_ub - 1) {
    memcpy(&z1_data[0], &y_data[0], loop_ub * sizeof(real_T));
  }

  for (loop_ub = 0; loop_ub + 1 <= *y_size; loop_ub++) {
    z1_data[loop_ub] = a_data[loop_ub] * a_data[loop_ub];
  }

  *y_size = z1_size_idx_0;
  if (0 <= z1_size_idx_0 - 1) {
    memcpy(&y_data[0], &z1_data[0], z1_size_idx_0 * sizeof(real_T));
  }
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static void Human_in_Loop_sort(real_T x_data[], int32_T idx_data[], int32_T
  idx_size[2])
{
  int32_T nNaNs;
  real_T xwork_data[4];
  real_T x4[4];
  int8_T idx4[4];
  int32_T ib;
  int8_T perm[4];
  int32_T i1;
  int32_T i2;
  int32_T i3;
  int32_T i4;
  idx_size[0] = 1;
  idx_size[1] = 4;
  idx_data[0] = 0;
  x4[0] = 0.0;
  idx4[0] = 0;
  xwork_data[0] = 0.0;
  idx_data[1] = 0;
  x4[1] = 0.0;
  idx4[1] = 0;
  xwork_data[1] = 0.0;
  idx_data[2] = 0;
  x4[2] = 0.0;
  idx4[2] = 0;
  xwork_data[2] = 0.0;
  idx_data[3] = 0;
  x4[3] = 0.0;
  idx4[3] = 0;
  xwork_data[3] = 0.0;
  nNaNs = 0;
  ib = 0;
  if (rtIsNaN(x_data[0])) {
    idx_data[3] = 1;
    xwork_data[3] = x_data[0];
    nNaNs = 1;
  } else {
    ib = 1;
    idx4[0] = 1;
    x4[0] = x_data[0];
  }

  if (rtIsNaN(x_data[1])) {
    idx_data[3 - nNaNs] = 2;
    xwork_data[3 - nNaNs] = x_data[1];
    nNaNs++;
  } else {
    ib++;
    idx4[ib - 1] = 2;
    x4[ib - 1] = x_data[1];
  }

  if (rtIsNaN(x_data[2])) {
    idx_data[3 - nNaNs] = 3;
    xwork_data[3 - nNaNs] = x_data[2];
    nNaNs++;
  } else {
    ib++;
    idx4[ib - 1] = 3;
    x4[ib - 1] = x_data[2];
  }

  if (rtIsNaN(x_data[3])) {
    idx_data[3 - nNaNs] = 4;
    xwork_data[3 - nNaNs] = x_data[3];
    nNaNs++;
  } else {
    ib++;
    idx4[ib - 1] = 4;
    x4[ib - 1] = x_data[3];
    if (ib == 4) {
      ib = 3 - nNaNs;
      if (x4[0] <= x4[1]) {
        i1 = 1;
        i2 = 2;
      } else {
        i1 = 2;
        i2 = 1;
      }

      if (x4[2] <= x4[3]) {
        i3 = 3;
        i4 = 4;
      } else {
        i3 = 4;
        i4 = 3;
      }

      if (x4[i1 - 1] <= x4[i3 - 1]) {
        if (x4[i2 - 1] <= x4[i3 - 1]) {
          perm[0] = (int8_T)i1;
          perm[1] = (int8_T)i2;
          perm[2] = (int8_T)i3;
          perm[3] = (int8_T)i4;
        } else if (x4[i2 - 1] <= x4[i4 - 1]) {
          perm[0] = (int8_T)i1;
          perm[1] = (int8_T)i3;
          perm[2] = (int8_T)i2;
          perm[3] = (int8_T)i4;
        } else {
          perm[0] = (int8_T)i1;
          perm[1] = (int8_T)i3;
          perm[2] = (int8_T)i4;
          perm[3] = (int8_T)i2;
        }
      } else if (x4[i1 - 1] <= x4[i4 - 1]) {
        if (x4[i2 - 1] <= x4[i4 - 1]) {
          perm[0] = (int8_T)i3;
          perm[1] = (int8_T)i1;
          perm[2] = (int8_T)i2;
          perm[3] = (int8_T)i4;
        } else {
          perm[0] = (int8_T)i3;
          perm[1] = (int8_T)i1;
          perm[2] = (int8_T)i4;
          perm[3] = (int8_T)i2;
        }
      } else {
        perm[0] = (int8_T)i3;
        perm[1] = (int8_T)i4;
        perm[2] = (int8_T)i1;
        perm[3] = (int8_T)i2;
      }

      idx_data[ib - 3] = idx4[perm[0] - 1];
      idx_data[ib - 2] = idx4[perm[1] - 1];
      idx_data[ib - 1] = idx4[perm[2] - 1];
      idx_data[ib] = idx4[perm[3] - 1];
      x_data[ib - 3] = x4[perm[0] - 1];
      x_data[ib - 2] = x4[perm[1] - 1];
      x_data[ib - 1] = x4[perm[2] - 1];
      x_data[ib] = x4[perm[3] - 1];
      ib = 0;
    }
  }

  i2 = 4;
  if (ib > 0) {
    perm[1] = 0;
    perm[2] = 0;
    perm[3] = 0;
    if (ib == 1) {
      perm[0] = 1;
    } else if (ib == 2) {
      if (x4[0] <= x4[1]) {
        perm[0] = 1;
        perm[1] = 2;
      } else {
        perm[0] = 2;
        perm[1] = 1;
      }
    } else if (x4[0] <= x4[1]) {
      if (x4[1] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 2;
        perm[2] = 3;
      } else if (x4[0] <= x4[2]) {
        perm[0] = 1;
        perm[1] = 3;
        perm[2] = 2;
      } else {
        perm[0] = 3;
        perm[1] = 1;
        perm[2] = 2;
      }
    } else if (x4[0] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 1;
      perm[2] = 3;
    } else if (x4[1] <= x4[2]) {
      perm[0] = 2;
      perm[1] = 3;
      perm[2] = 1;
    } else {
      perm[0] = 3;
      perm[1] = 2;
      perm[2] = 1;
    }

    while (i2 - 3 <= ib) {
      idx_data[(i2 - nNaNs) - ib] = idx4[perm[i2 - 4] - 1];
      x_data[(i2 - nNaNs) - ib] = x4[perm[i2 - 4] - 1];
      i2++;
    }
  }

  i2 = nNaNs >> 1;
  for (ib = 1; ib <= i2; ib++) {
    i1 = idx_data[(ib - nNaNs) + 3];
    idx_data[(ib - nNaNs) + 3] = idx_data[4 - ib];
    idx_data[4 - ib] = i1;
    x_data[(ib - nNaNs) + 3] = xwork_data[4 - ib];
    x_data[4 - ib] = xwork_data[(ib - nNaNs) + 3];
  }

  if ((nNaNs & 1) != 0) {
    x_data[(i2 - nNaNs) + 4] = xwork_data[(i2 - nNaNs) + 4];
  }
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static void Human_in_Loop_diag(const real_T v_data[], const int32_T *v_size,
  real_T d_data[], int32_T d_size[2])
{
  int32_T loop_ub;
  int8_T b_idx_0;
  int8_T b_idx_1;
  b_idx_0 = (int8_T)*v_size;
  b_idx_1 = (int8_T)*v_size;
  d_size[0] = b_idx_0;
  d_size[1] = b_idx_1;
  loop_ub = b_idx_0 * b_idx_1;
  if (0 <= loop_ub - 1) {
    memset(&d_data[0], 0, loop_ub * sizeof(real_T));
  }

  for (loop_ub = 0; loop_ub + 1 <= *v_size; loop_ub++) {
    d_data[loop_ub + d_size[0] * loop_ub] = v_data[loop_ub];
  }
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static void Human_in_Loop_triu(real_T x_data[])
{
  x_data[0] = 0.0;
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static void Human_in_Loop_eig(const real_T A_data[], creal_T V_data[], int32_T
  V_size[2], creal_T D_data[], int32_T D_size[2])
{
  boolean_T p;
  p = ((!rtIsInf(A_data[0])) && (!rtIsNaN(A_data[0])));
  if (!p) {
    V_size[0] = 1;
    V_size[1] = 1;
    V_data[0].re = (rtNaN);
    V_data[0].im = 0.0;
    D_size[0] = 1;
    D_size[1] = 1;
    D_data[0].re = (rtNaN);
    D_data[0].im = 0.0;
  } else {
    V_size[0] = 1;
    V_size[1] = 1;
    V_data[0].re = 1.0;
    V_data[0].im = 0.0;
    D_size[0] = 1;
    D_size[1] = 1;
    D_data[0].re = A_data[0];
    D_data[0].im = 0.0;
  }
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static void Human_in_Loop_sqrt(real_T x_data[])
{
  x_data[0] = sqrt(x_data[0]);
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static void Human_in__genrand_uint32_vector(uint32_T mt[625], uint32_T u[2])
{
  uint32_T mti;
  uint32_T y;
  int32_T j;
  int32_T kk;
  for (j = 0; j < 2; j++) {
    mti = mt[624] + 1U;
    if (mti >= 625U) {
      for (kk = 0; kk < 227; kk++) {
        y = (mt[kk + 1] & 2147483647U) | (mt[kk] & 2147483648U);
        if ((int32_T)(y & 1U) == 0) {
          mti = y >> 1U;
        } else {
          mti = y >> 1U ^ 2567483615U;
        }

        mt[kk] = mt[kk + 397] ^ mti;
      }

      for (kk = 0; kk < 396; kk++) {
        y = (mt[kk + 227] & 2147483648U) | (mt[kk + 228] & 2147483647U);
        if ((int32_T)(y & 1U) == 0) {
          mti = y >> 1U;
        } else {
          mti = y >> 1U ^ 2567483615U;
        }

        mt[kk + 227] = mt[kk] ^ mti;
      }

      y = (mt[623] & 2147483648U) | (mt[0] & 2147483647U);
      if ((int32_T)(y & 1U) == 0) {
        mti = y >> 1U;
      } else {
        mti = y >> 1U ^ 2567483615U;
      }

      mt[623] = mt[396] ^ mti;
      mti = 1U;
    }

    y = mt[(int32_T)mti - 1];
    mt[624] = mti;
    y ^= y >> 11U;
    y ^= y << 7U & 2636928640U;
    y ^= y << 15U & 4022730752U;
    y ^= y >> 18U;
    u[j] = y;
  }
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static real_T Human_in_Loop_genrandu(uint32_T mt[625])
{
  real_T r;
  uint32_T u[2];

  /* ========================= COPYRIGHT NOTICE ============================ */
  /*  This is a uniform (0,1) pseudorandom number generator based on:        */
  /*                                                                         */
  /*  A C-program for MT19937, with initialization improved 2002/1/26.       */
  /*  Coded by Takuji Nishimura and Makoto Matsumoto.                        */
  /*                                                                         */
  /*  Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,      */
  /*  All rights reserved.                                                   */
  /*                                                                         */
  /*  Redistribution and use in source and binary forms, with or without     */
  /*  modification, are permitted provided that the following conditions     */
  /*  are met:                                                               */
  /*                                                                         */
  /*    1. Redistributions of source code must retain the above copyright    */
  /*       notice, this list of conditions and the following disclaimer.     */
  /*                                                                         */
  /*    2. Redistributions in binary form must reproduce the above copyright */
  /*       notice, this list of conditions and the following disclaimer      */
  /*       in the documentation and/or other materials provided with the     */
  /*       distribution.                                                     */
  /*                                                                         */
  /*    3. The names of its contributors may not be used to endorse or       */
  /*       promote products derived from this software without specific      */
  /*       prior written permission.                                         */
  /*                                                                         */
  /*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS    */
  /*  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT      */
  /*  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR  */
  /*  A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT  */
  /*  OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,  */
  /*  SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT       */
  /*  LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,  */
  /*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY  */
  /*  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT    */
  /*  (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE */
  /*  OF THIS  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  */
  /*                                                                         */
  /* =============================   END   ================================= */
  do {
    Human_in__genrand_uint32_vector(mt, u);
    r = ((real_T)(u[0] >> 5U) * 6.7108864E+7 + (real_T)(u[1] >> 6U)) *
      1.1102230246251565E-16;
  } while (r == 0.0);

  return r;
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static real_T Human_in_L_eml_rand_mt19937ar_h(uint32_T state[625])
{
  real_T r;
  int32_T i;
  real_T x;
  uint32_T u32[2];
  real_T c_u;
  static const real_T b[257] = { 0.0, 0.215241895984875, 0.286174591792068,
    0.335737519214422, 0.375121332878378, 0.408389134611989, 0.43751840220787,
    0.46363433679088, 0.487443966139235, 0.50942332960209, 0.529909720661557,
    0.549151702327164, 0.567338257053817, 0.584616766106378, 0.601104617755991,
    0.61689699000775, 0.63207223638606, 0.646695714894993, 0.660822574244419,
    0.674499822837293, 0.687767892795788, 0.700661841106814, 0.713212285190975,
    0.725446140909999, 0.737387211434295, 0.749056662017815, 0.760473406430107,
    0.771654424224568, 0.782615023307232, 0.793369058840623, 0.80392911698997,
    0.814306670135215, 0.824512208752291, 0.834555354086381, 0.844444954909153,
    0.854189171008163, 0.863795545553308, 0.87327106808886, 0.882622229585165,
    0.891855070732941, 0.900975224461221, 0.909987953496718, 0.91889818364959,
    0.927710533401999, 0.936429340286575, 0.945058684468165, 0.953602409881086,
    0.96206414322304, 0.970447311064224, 0.978755155294224, 0.986990747099062,
    0.99515699963509, 1.00325667954467, 1.01129241744, 1.01926671746548,
    1.02718196603564, 1.03504043983344, 1.04284431314415, 1.05059566459093,
    1.05829648333067, 1.06594867476212, 1.07355406579244, 1.0811144097034,
    1.08863139065398, 1.09610662785202, 1.10354167942464, 1.11093804601357,
    1.11829717411934, 1.12562045921553, 1.13290924865253, 1.14016484436815,
    1.14738850542085, 1.15458145035993, 1.16174485944561, 1.16887987673083,
    1.17598761201545, 1.18306914268269, 1.19012551542669, 1.19715774787944,
    1.20416683014438, 1.2111537262437, 1.21811937548548, 1.22506469375653,
    1.23199057474614, 1.23889789110569, 1.24578749554863, 1.2526602218949,
    1.25951688606371, 1.26635828701823, 1.27318520766536, 1.27999841571382,
    1.28679866449324, 1.29358669373695, 1.30036323033084, 1.30712898903073,
    1.31388467315022, 1.32063097522106, 1.32736857762793, 1.33409815321936,
    1.3408203658964, 1.34753587118059, 1.35424531676263, 1.36094934303328,
    1.36764858359748, 1.37434366577317, 1.38103521107586, 1.38772383568998,
    1.39441015092814, 1.40109476367925, 1.4077782768464, 1.41446128977547,
    1.42114439867531, 1.42782819703026, 1.43451327600589, 1.44120022484872,
    1.44788963128058, 1.45458208188841, 1.46127816251028, 1.46797845861808,
    1.47468355569786, 1.48139403962819, 1.48811049705745, 1.49483351578049,
    1.50156368511546, 1.50830159628131, 1.51504784277671, 1.521803020761,
    1.52856772943771, 1.53534257144151, 1.542128153229, 1.54892508547417,
    1.55573398346918, 1.56255546753104, 1.56939016341512, 1.57623870273591,
    1.58310172339603, 1.58997987002419, 1.59687379442279, 1.60378415602609,
    1.61071162236983, 1.61765686957301, 1.62462058283303, 1.63160345693487,
    1.63860619677555, 1.64562951790478, 1.65267414708306, 1.65974082285818,
    1.66683029616166, 1.67394333092612, 1.68108070472517, 1.68824320943719,
    1.69543165193456, 1.70264685479992, 1.7098896570713, 1.71716091501782,
    1.72446150294804, 1.73179231405296, 1.73915426128591, 1.74654827828172,
    1.75397532031767, 1.76143636531891, 1.76893241491127, 1.77646449552452,
    1.78403365954944, 1.79164098655216, 1.79928758454972, 1.80697459135082,
    1.81470317596628, 1.82247454009388, 1.83028991968276, 1.83815058658281,
    1.84605785028518, 1.8540130597602, 1.86201760539967, 1.87007292107127,
    1.878180486293, 1.88634182853678, 1.8945585256707, 1.90283220855043,
    1.91116456377125, 1.91955733659319, 1.92801233405266, 1.93653142827569,
    1.94511656000868, 1.95376974238465, 1.96249306494436, 1.97128869793366,
    1.98015889690048, 1.98910600761744, 1.99813247135842, 2.00724083056053,
    2.0164337349062, 2.02571394786385, 2.03508435372962, 2.04454796521753,
    2.05410793165065, 2.06376754781173, 2.07353026351874, 2.0833996939983,
    2.09337963113879, 2.10347405571488, 2.11368715068665, 2.12402331568952,
    2.13448718284602, 2.14508363404789, 2.15581781987674, 2.16669518035431,
    2.17772146774029, 2.18890277162636, 2.20024554661128, 2.21175664288416,
    2.22344334009251, 2.23531338492992, 2.24737503294739, 2.25963709517379,
    2.27210899022838, 2.28480080272449, 2.29772334890286, 2.31088825060137,
    2.32430801887113, 2.33799614879653, 2.35196722737914, 2.36623705671729,
    2.38082279517208, 2.39574311978193, 2.41101841390112, 2.42667098493715,
    2.44272531820036, 2.4592083743347, 2.47614993967052, 2.49358304127105,
    2.51154444162669, 2.53007523215985, 2.54922155032478, 2.56903545268184,
    2.58957598670829, 2.61091051848882, 2.63311639363158, 2.65628303757674,
    2.68051464328574, 2.70593365612306, 2.73268535904401, 2.76094400527999,
    2.79092117400193, 2.82287739682644, 2.85713873087322, 2.89412105361341,
    2.93436686720889, 2.97860327988184, 3.02783779176959, 3.08352613200214,
    3.147889289518, 3.2245750520478, 3.32024473383983, 3.44927829856143,
    3.65415288536101, 3.91075795952492 };

  static const real_T c[257] = { 1.0, 0.977101701267673, 0.959879091800108,
    0.9451989534423, 0.932060075959231, 0.919991505039348, 0.908726440052131,
    0.898095921898344, 0.887984660755834, 0.878309655808918, 0.869008688036857,
    0.860033621196332, 0.851346258458678, 0.842915653112205, 0.834716292986884,
    0.826726833946222, 0.818929191603703, 0.811307874312656, 0.803849483170964,
    0.796542330422959, 0.789376143566025, 0.782341832654803, 0.775431304981187,
    0.768637315798486, 0.761953346836795, 0.755373506507096, 0.748892447219157,
    0.742505296340151, 0.736207598126863, 0.729995264561476, 0.72386453346863,
    0.717811932630722, 0.711834248878248, 0.705928501332754, 0.700091918136512,
    0.694321916126117, 0.688616083004672, 0.682972161644995, 0.677388036218774,
    0.671861719897082, 0.66639134390875, 0.660975147776663, 0.655611470579697,
    0.650298743110817, 0.645035480820822, 0.639820277453057, 0.634651799287624,
    0.629528779924837, 0.624450015547027, 0.619414360605834, 0.614420723888914,
    0.609468064925773, 0.604555390697468, 0.599681752619125, 0.594846243767987,
    0.590047996332826, 0.585286179263371, 0.580559996100791, 0.575868682972354,
    0.571211506735253, 0.566587763256165, 0.561996775814525, 0.557437893618766,
    0.552910490425833, 0.548413963255266, 0.543947731190026, 0.539511234256952,
    0.535103932380458, 0.530725304403662, 0.526374847171684, 0.522052074672322,
    0.517756517229756, 0.513487720747327, 0.509245245995748, 0.505028667943468,
    0.500837575126149, 0.49667156905249, 0.492530263643869, 0.488413284705458,
    0.484320269426683, 0.480250865909047, 0.476204732719506, 0.47218153846773,
    0.468180961405694, 0.464202689048174, 0.460246417812843, 0.456311852678716,
    0.452398706861849, 0.448506701507203, 0.444635565395739, 0.440785034665804,
    0.436954852547985, 0.433144769112652, 0.429354541029442, 0.425583931338022,
    0.421832709229496, 0.418100649837848, 0.414387534040891, 0.410693148270188,
    0.407017284329473, 0.403359739221114, 0.399720314980197, 0.396098818515832,
    0.392495061459315, 0.388908860018789, 0.385340034840077, 0.381788410873393,
    0.378253817245619, 0.374736087137891, 0.371235057668239, 0.367750569779032,
    0.364282468129004, 0.360830600989648, 0.357394820145781, 0.353974980800077,
    0.350570941481406, 0.347182563956794, 0.343809713146851, 0.340452257044522,
    0.337110066637006, 0.333783015830718, 0.330470981379163, 0.327173842813601,
    0.323891482376391, 0.320623784956905, 0.317370638029914, 0.314131931596337,
    0.310907558126286, 0.307697412504292, 0.30450139197665, 0.301319396100803,
    0.298151326696685, 0.294997087799962, 0.291856585617095, 0.288729728482183,
    0.285616426815502, 0.282516593083708, 0.279430141761638, 0.276356989295668,
    0.273297054068577, 0.270250256365875, 0.267216518343561, 0.264195763997261,
    0.261187919132721, 0.258192911337619, 0.255210669954662, 0.252241126055942,
    0.249284212418529, 0.246339863501264, 0.24340801542275, 0.240488605940501,
    0.237581574431238, 0.23468686187233, 0.231804410824339, 0.228934165414681,
    0.226076071322381, 0.223230075763918, 0.220396127480152, 0.217574176724331,
    0.214764175251174, 0.211966076307031, 0.209179834621125, 0.206405406397881,
    0.203642749310335, 0.200891822494657, 0.198152586545776, 0.195425003514135,
    0.192709036903589, 0.190004651670465, 0.187311814223801, 0.1846304924268,
    0.181960655599523, 0.179302274522848, 0.176655321443735, 0.174019770081839,
    0.171395595637506, 0.168782774801212, 0.166181285764482, 0.163591108232366,
    0.161012223437511, 0.158444614155925, 0.15588826472448, 0.153343161060263,
    0.150809290681846, 0.148286642732575, 0.145775208005994, 0.143274978973514,
    0.140785949814445, 0.138308116448551, 0.135841476571254, 0.133386029691669,
    0.130941777173644, 0.12850872228, 0.126086870220186, 0.123676228201597,
    0.12127680548479, 0.11888861344291, 0.116511665625611, 0.114145977827839,
    0.111791568163838, 0.109448457146812, 0.107116667774684, 0.104796225622487,
    0.102487158941935, 0.10018949876881, 0.0979032790388625, 0.095628536713009,
    0.093365311912691, 0.0911136480663738, 0.0888735920682759,
    0.0866451944505581, 0.0844285095703535, 0.082223595813203,
    0.0800305158146631, 0.0778493367020961, 0.0756801303589272,
    0.0735229737139814, 0.0713779490588905, 0.0692451443970068,
    0.0671246538277886, 0.065016577971243, 0.0629210244377582, 0.06083810834954,
    0.0587679529209339, 0.0567106901062031, 0.0546664613248891,
    0.0526354182767924, 0.0506177238609479, 0.0486135532158687,
    0.0466230949019305, 0.0446465522512946, 0.0426841449164746,
    0.0407361106559411, 0.0388027074045262, 0.0368842156885674,
    0.0349809414617162, 0.0330932194585786, 0.0312214171919203,
    0.0293659397581334, 0.0275272356696031, 0.0257058040085489,
    0.0239022033057959, 0.0221170627073089, 0.0203510962300445,
    0.0186051212757247, 0.0168800831525432, 0.0151770883079353,
    0.0134974506017399, 0.0118427578579079, 0.0102149714397015,
    0.00861658276939875, 0.00705087547137324, 0.00552240329925101,
    0.00403797259336304, 0.00260907274610216, 0.0012602859304986,
    0.000477467764609386 };

  int32_T exitg1;
  do {
    exitg1 = 0;
    Human_in__genrand_uint32_vector(state, u32);
    i = (int32_T)((u32[1] >> 24U) + 1U);
    r = (((real_T)(u32[0] >> 3U) * 1.6777216E+7 + (real_T)((int32_T)u32[1] &
           16777215)) * 2.2204460492503131E-16 - 1.0) * b[i];
    if (fabs(r) <= b[i - 1]) {
      exitg1 = 1;
    } else if (i < 256) {
      x = Human_in_Loop_genrandu(state);
      if ((c[i - 1] - c[i]) * x + c[i] < exp(-0.5 * r * r)) {
        exitg1 = 1;
      }
    } else {
      do {
        x = Human_in_Loop_genrandu(state);
        x = log(x) * 0.273661237329758;
        c_u = Human_in_Loop_genrandu(state);
      } while (!(-2.0 * log(c_u) > x * x));

      if (r < 0.0) {
        r = x - 3.65415288536101;
      } else {
        r = 3.65415288536101 - x;
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return r;
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static void Human_in_Loop_randn(real_T r_data[], int32_T *r_size)
{
  real_T b;
  *r_size = 1;
  b = Human_in_L_eml_rand_mt19937ar_h(Human_in_Loop_DW.state_l);
  r_data[0] = b;
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static boolean_T Human_in_Loop_sortLE(const real_T v_data[], int32_T irow1,
  int32_T irow2)
{
  boolean_T p;
  boolean_T b;
  p = true;
  if ((v_data[irow1 - 1] == v_data[irow2 - 1]) || (rtIsNaN(v_data[irow1 - 1]) &&
       rtIsNaN(v_data[irow2 - 1]))) {
    b = true;
  } else {
    b = false;
  }

  if ((!b) && (!(v_data[irow1 - 1] <= v_data[irow2 - 1])) && (!rtIsNaN
       (v_data[irow2 - 1]))) {
    p = false;
  }

  return p;
}

/* Function for MATLAB Function: '<S4>/CMAES' */
static void Human_in_Loop_sortrows(real_T y_data[], real_T ndx_data[], int32_T
  *ndx_size)
{
  int32_T idx_data[4];
  int32_T iwork_data[4];
  int32_T b_k;
  int32_T i;
  int32_T j;
  int32_T pEnd;
  int32_T p;
  int32_T q;
  real_T ycol_data_idx_0;
  real_T ycol_data_idx_1;
  real_T ycol_data_idx_2;
  real_T ycol_data_idx_3;
  *ndx_size = 4;
  if (Human_in_Loop_sortLE(y_data, 1, 2)) {
    idx_data[0] = 1;
    idx_data[1] = 2;
  } else {
    idx_data[0] = 2;
    idx_data[1] = 1;
  }

  if (Human_in_Loop_sortLE(y_data, 3, 4)) {
    idx_data[2] = 3;
    idx_data[3] = 4;
  } else {
    idx_data[2] = 4;
    idx_data[3] = 3;
  }

  i = 2;
  while (i < 4) {
    j = 1;
    pEnd = 1 + i;
    while (pEnd < 5) {
      p = j;
      q = pEnd;
      b_k = 0;
      while (b_k + 1 <= 5 - j) {
        if (Human_in_Loop_sortLE(y_data, idx_data[p - 1], idx_data[q - 1])) {
          iwork_data[b_k] = idx_data[p - 1];
          p++;
          if (p == pEnd) {
            while (q < 5) {
              b_k++;
              iwork_data[b_k] = idx_data[q - 1];
              q++;
            }
          }
        } else {
          iwork_data[b_k] = idx_data[q - 1];
          q++;
          if (q == 5) {
            while (p < pEnd) {
              b_k++;
              iwork_data[b_k] = idx_data[p - 1];
              p++;
            }
          }
        }

        b_k++;
      }

      for (pEnd = 0; pEnd + 1 <= 5 - j; pEnd++) {
        idx_data[(j + pEnd) - 1] = iwork_data[pEnd];
      }

      j = 5;
      pEnd = 5 + i;
    }

    i <<= 1;
  }

  ycol_data_idx_0 = y_data[idx_data[0] - 1];
  ycol_data_idx_1 = y_data[idx_data[1] - 1];
  ycol_data_idx_2 = y_data[idx_data[2] - 1];
  ycol_data_idx_3 = y_data[idx_data[3] - 1];
  y_data[0] = ycol_data_idx_0;
  ndx_data[0] = idx_data[0];
  y_data[1] = ycol_data_idx_1;
  ndx_data[1] = idx_data[1];
  y_data[2] = ycol_data_idx_2;
  ndx_data[2] = idx_data[2];
  y_data[3] = ycol_data_idx_3;
  ndx_data[3] = idx_data[3];
}

/* Function for MATLAB Function: '<S19>/MATLAB Function1' */
static void Human_in_Loop_mean(const real_T x_data[], const int32_T x_size[2],
  real_T y[26])
{
  int32_T xoffset;
  int32_T j;
  int32_T b_j;
  real_T y_0;
  if (x_size[1] == 0) {
    memset(&y[0], 0, 26U * sizeof(real_T));
  } else {
    memcpy(&y[0], &x_data[0], 26U * sizeof(real_T));
    for (j = 2; j <= x_size[1]; j++) {
      xoffset = (j - 1) * 26;
      for (b_j = 0; b_j < 26; b_j++) {
        y_0 = y[b_j];
        y_0 += x_data[xoffset + b_j];
        y[b_j] = y_0;
      }
    }
  }

  j = x_size[1];
  for (xoffset = 0; xoffset < 26; xoffset++) {
    y_0 = y[xoffset];
    y_0 /= (real_T)j;
    y[xoffset] = y_0;
  }
}

/* Function for MATLAB Function: '<S49>/MATLAB Function' */
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

/* Function for MATLAB Function: '<S49>/MATLAB Function' */
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

/* Function for MATLAB Function: '<S49>/MATLAB Function' */
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
  real_T weights_data[2];
  real_T mueff;
  real_T cc;
  real_T cs;
  real_T c1;
  real_T cmu;
  boolean_T hsig;
  creal_T Bz_data;
  creal_T Dz_data;
  real_T arxx_data[4];
  real_T b_index_data[4];
  real_T b_y_data[2];
  int32_T iidx_data[4];
  int32_T ia;
  real_T x;
  real_T d_y_data[2];
  int32_T br;
  int32_T c_ar;
  int32_T c_ia;
  int32_T c_ib;
  int32_T f_ic;
  real_T zmean_data;
  real_T c_data[2];
  real_T x_0[6];
  real_T x_1[2];
  real_T A[30];
  real_T tol;
  real_T B[15];
  static const int8_T b_A[30] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15 };

  ZCEventType zcEvent;
  int32_T i;
  real_T tmp_data[2];
  real_T tmp;
  real_T tmp_0;
  real_T tmp_data_0[4];
  real_T tmp_data_1[2];
  real_T tmp_data_2[4];
  real_T tmp_data_3[2];
  int32_T weights_size;
  int32_T Dz_size[2];
  int32_T b_y_size[2];
  int32_T tmp_size[2];
  int8_T tmp_1;
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
    /* S-Function (rti_commonblock): '<S103>/S-Function1' */
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

    /* Gain: '<S49>/Gain' */
    Human_in_Loop_B.Gain = Human_in_Loop_P.Gain_Gain *
      Human_in_Loop_B.SFunction1;

    /* DiscreteFilter: '<S49>/0.4low2' */
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

    /* MATLAB Function: '<S49>/Data process' */
    /* MATLAB Function 'Sensor Data/Torque module/Data process': '<S104>:1' */
    if (Human_in_Loop_P.Dataprocess_BT_RESET_TORQUE) {
      /* '<S104>:1:10' */
      /* '<S104>:1:11' */
      Human_in_Loop_DW.torque_zero = Human_in_Loop_B.u4low2 *
        Human_in_Loop_P.Dataprocess_load_vol_gain +
        Human_in_Loop_P.Dataprocess_load_vol_offset;
    }

    /* '<S104>:1:14' */
    Human_in_Loop_B.torque = (Human_in_Loop_B.u4low2 *
      Human_in_Loop_P.Dataprocess_load_vol_gain +
      Human_in_Loop_P.Dataprocess_load_vol_offset) -
      Human_in_Loop_DW.torque_zero;

    /* End of MATLAB Function: '<S49>/Data process' */

    /* UnitDelay: '<S102>/Unit Delay1' */
    Human_in_Loop_B.x2k1 = Human_in_Loop_DW.UnitDelay1_DSTATE;

    /* UnitDelay: '<S102>/Unit Delay' */
    Human_in_Loop_B.x1k1 = Human_in_Loop_DW.UnitDelay_DSTATE;

    /* Gain: '<S102>/Gain1' */
    c1 = Human_in_Loop_P.uOrderTD_T1 * Human_in_Loop_P.uOrderTD_T2;
    tol = 1.0 / c1;
    Human_in_Loop_B.Gain1 = tol * Human_in_Loop_B.x1k1;

    /* Gain: '<S102>/Gain2' */
    tol = Human_in_Loop_P.uOrderTD_T1 + Human_in_Loop_P.uOrderTD_T2;
    c1 = Human_in_Loop_P.uOrderTD_T1 * Human_in_Loop_P.uOrderTD_T2;
    tol /= c1;
    Human_in_Loop_B.Gain2 = tol * Human_in_Loop_B.x2k1;

    /* UnitDelay: '<S102>/Unit Delay2' */
    Human_in_Loop_B.UnitDelay2 = Human_in_Loop_DW.UnitDelay2_DSTATE;

    /* Gain: '<S102>/Gain4' */
    c1 = Human_in_Loop_P.uOrderTD_T1 * Human_in_Loop_P.uOrderTD_T2;
    tol = 1.0 / c1;
    Human_in_Loop_B.Gain4 = tol * Human_in_Loop_B.UnitDelay2;

    /* Sum: '<S102>/Add2' */
    Human_in_Loop_B.Add2 = (Human_in_Loop_B.Gain1 + Human_in_Loop_B.Gain2) -
      Human_in_Loop_B.Gain4;

    /* Gain: '<S102>/Gain3' */
    Human_in_Loop_B.Gain3 = Human_in_Loop_P.uOrderTD_Ts * Human_in_Loop_B.Add2;

    /* Sum: '<S102>/Add1' */
    Human_in_Loop_B.x2k = Human_in_Loop_B.x2k1 - Human_in_Loop_B.Gain3;

    /* MATLAB Function: '<S49>/Mux' */
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

    /* S-Function (rti_commonblock): '<S81>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* Gain: '<S46>/Gain' */
    Human_in_Loop_B.Gain_l = Human_in_Loop_P.Gain_Gain_c *
      Human_in_Loop_B.SFunction1_o1;

    /* MATLAB Function: '<S46>/Data process' */
    /* MATLAB Function 'Sensor Data/Encoder module/Data process': '<S78>:1' */
    if (Human_in_Loop_P.Dataprocess_BT_RESET_ANKLE) {
      /* '<S78>:1:8' */
      /* '<S78>:1:9' */
      Human_in_Loop_DW.angle_zero_f = Human_in_Loop_B.Gain_l;
    }

    /* '<S78>:1:12' */
    Human_in_Loop_B.angle_m = Human_in_Loop_B.Gain_l -
      Human_in_Loop_DW.angle_zero_f;

    /* End of MATLAB Function: '<S46>/Data process' */

    /* Gain: '<S46>/Gain1' */
    Human_in_Loop_B.Gain1_b = Human_in_Loop_P.Gain1_Gain_e *
      Human_in_Loop_B.SFunction1_o2;

    /* MATLAB Function: '<S46>/Mux' */
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

    /* S-Function (rti_commonblock): '<S82>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* Gain: '<S46>/Gain2' */
    Human_in_Loop_B.Gain2_h = Human_in_Loop_P.Gain2_Gain_e *
      Human_in_Loop_B.SFunction1_o1_k;

    /* MATLAB Function: '<S46>/Data process1' */
    /* MATLAB Function 'Sensor Data/Encoder module/Data process1': '<S79>:1' */
    if (Human_in_Loop_P.Dataprocess1_BT_RESET_MOTOR) {
      /* '<S79>:1:8' */
      /* '<S79>:1:9' */
      Human_in_Loop_DW.angle_zero = Human_in_Loop_B.Gain2_h;
    }

    /* '<S79>:1:12' */
    Human_in_Loop_B.angle = Human_in_Loop_B.Gain2_h -
      Human_in_Loop_DW.angle_zero;

    /* End of MATLAB Function: '<S46>/Data process1' */

    /* Gain: '<S46>/Gain3' */
    Human_in_Loop_B.Gain3_m = Human_in_Loop_P.Gain3_Gain *
      Human_in_Loop_B.SFunction1_o2_k;

    /* UnitDelay: '<S77>/Unit Delay1' */
    Human_in_Loop_B.x2k1_k = Human_in_Loop_DW.UnitDelay1_DSTATE_h;

    /* UnitDelay: '<S77>/Unit Delay' */
    Human_in_Loop_B.x1k1_m = Human_in_Loop_DW.UnitDelay_DSTATE_e;

    /* Gain: '<S77>/Gain1' */
    c1 = Human_in_Loop_P.uOrderTD_T1_l * Human_in_Loop_P.uOrderTD_T2_h;
    tol = 1.0 / c1;
    Human_in_Loop_B.Gain1_m = tol * Human_in_Loop_B.x1k1_m;

    /* Gain: '<S77>/Gain2' */
    tol = Human_in_Loop_P.uOrderTD_T1_l + Human_in_Loop_P.uOrderTD_T2_h;
    c1 = Human_in_Loop_P.uOrderTD_T1_l * Human_in_Loop_P.uOrderTD_T2_h;
    tol /= c1;
    Human_in_Loop_B.Gain2_c = tol * Human_in_Loop_B.x2k1_k;

    /* UnitDelay: '<S77>/Unit Delay2' */
    Human_in_Loop_B.UnitDelay2_b = Human_in_Loop_DW.UnitDelay2_DSTATE_a;

    /* Gain: '<S77>/Gain4' */
    c1 = Human_in_Loop_P.uOrderTD_T1_l * Human_in_Loop_P.uOrderTD_T2_h;
    tol = 1.0 / c1;
    Human_in_Loop_B.Gain4_g = tol * Human_in_Loop_B.UnitDelay2_b;

    /* Sum: '<S77>/Add2' */
    Human_in_Loop_B.Add2_m = (Human_in_Loop_B.Gain1_m + Human_in_Loop_B.Gain2_c)
      - Human_in_Loop_B.Gain4_g;

    /* Gain: '<S77>/Gain3' */
    Human_in_Loop_B.Gain3_k = Human_in_Loop_P.uOrderTD_Ts_n *
      Human_in_Loop_B.Add2_m;

    /* Sum: '<S77>/Add1' */
    Human_in_Loop_B.x2k_i = Human_in_Loop_B.x2k1_k - Human_in_Loop_B.Gain3_k;

    /* MATLAB Function: '<S46>/Mux1' */
    /* MATLAB Function 'Sensor Data/Encoder module/Mux1': '<S86>:1' */
    /* '<S86>:1:3' */
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

    /* S-Function (rti_commonblock): '<S88>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* MATLAB Function: '<S47>/FootSwitch Filter' */
    /* MATLAB Function 'Sensor Data/FootSwitch module/FootSwitch Filter': '<S90>:1' */
    /* '<S90>:1:6' */
    if (Human_in_Loop_DW.foot_state == 0.0) {
      /* '<S90>:1:14' */
      if (Human_in_Loop_B.SFunction1_a) {
        /* '<S90>:1:15' */
        /* '<S90>:1:16' */
        Human_in_Loop_DW.foot_state = 1.0;
      } else {
        /* '<S90>:1:18' */
        Human_in_Loop_DW.foot_state = 0.0;
      }
    } else if (Human_in_Loop_B.SFunction1_a) {
      /* '<S90>:1:21' */
      /* '<S90>:1:22' */
      Human_in_Loop_DW.foot_state = 1.0;
    } else {
      /* '<S90>:1:24' */
      Human_in_Loop_DW.filter_time += 0.0002;
      if (Human_in_Loop_DW.filter_time > 0.4) {
        /* '<S90>:1:25' */
        /* '<S90>:1:26' */
        Human_in_Loop_DW.filter_time = 0.0;

        /* '<S90>:1:27' */
        Human_in_Loop_DW.foot_state = 0.0;
      }
    }

    /* '<S90>:1:32' */
    Human_in_Loop_B.state_c = Human_in_Loop_DW.foot_state;

    /* End of MATLAB Function: '<S47>/FootSwitch Filter' */

    /* MATLAB Function: '<S8>/State Machine' */
    /* MATLAB Function 'State Module/State Machine': '<S108>:1' */
    /* '<S108>:1:20' */
    /* '<S108>:1:21' */
    /* '<S108>:1:22' */
    /* '<S108>:1:23' */
    /* '<S108>:1:26' */
    /* '<S108>:1:27' */
    /* '<S108>:1:28' */
    /* '<S108>:1:29' */
    /* '<S108>:1:30' */
    if (Human_in_Loop_P.StateMachine_BT_RUN) {
      /* '<S108>:1:33' */
      /* '<S108>:1:34' */
      Human_in_Loop_DW.bt_run = 1.0;
    }

    if (Human_in_Loop_P.StateMachine_BT_CALIB) {
      /* '<S108>:1:36' */
      /* '<S108>:1:37' */
      Human_in_Loop_DW.reg_mode = 4.0;
    }

    if (Human_in_Loop_P.StateMachine_BT_SLACK) {
      /* '<S108>:1:39' */
      /* '<S108>:1:40' */
      Human_in_Loop_DW.reg_mode = 3.0;
    }

    if (Human_in_Loop_P.StateMachine_BT_IDLE) {
      /* '<S108>:1:42' */
      /* '<S108>:1:43' */
      Human_in_Loop_DW.reg_mode = 1.0;
    }

    if (Human_in_Loop_P.StateMachine_BT_ERROR) {
      /* '<S108>:1:45' */
      /* '<S108>:1:46' */
      Human_in_Loop_DW.reg_mode = 0.0;
    }

    if ((Human_in_Loop_DW.bt_run == 1.0) && (Human_in_Loop_DW.reg_last_switch ==
         0.0) && (Human_in_Loop_B.state_c == 1.0)) {
      /* '<S108>:1:49' */
      /* '<S108>:1:50' */
      Human_in_Loop_DW.reg_mode = 2.0;

      /* '<S108>:1:51' */
      Human_in_Loop_DW.bt_run = 0.0;
    }

    if ((Human_in_Loop_DW.reg_mode == 2.0) || (Human_in_Loop_DW.reg_mode == 1.0))
    {
      /* '<S108>:1:54' */
      if ((Human_in_Loop_DW.reg_last_switch == 0.0) && (Human_in_Loop_B.state_c ==
           1.0)) {
        /* '<S108>:1:55' */
        /* '<S108>:1:56' */
        Human_in_Loop_DW.reg_state = 1.0;

        /* '<S108>:1:57' */
        Human_in_Loop_DW.reg_stride_time = 0.618 *
          Human_in_Loop_DW.reg_stride_time + 0.382 *
          Human_in_Loop_DW.reg_stride_time_count;

        /* '<S108>:1:58' */
        Human_in_Loop_DW.reg_stride_time_count = 0.0;
      } else if ((Human_in_Loop_DW.reg_state == 1.0) &&
                 (Human_in_Loop_DW.reg_stride_time_count > 0.65 *
                  Human_in_Loop_DW.reg_stride_time)) {
        /* '<S108>:1:59' */
        /* '<S108>:1:60' */
        Human_in_Loop_DW.reg_state = 0.0;

        /* '<S108>:1:61' */
        Human_in_Loop_DW.reg_stride_time_count += 0.0002;
      } else {
        /* '<S108>:1:63' */
        Human_in_Loop_DW.reg_stride_time_count += 0.0002;
      }
    }

    /* '<S108>:1:67' */
    Human_in_Loop_DW.reg_last_switch = Human_in_Loop_B.state_c;
    if (Human_in_Loop_DW.reg_stride_time > 1.5) {
      /* '<S108>:1:68' */
      /* '<S108>:1:69' */
      Human_in_Loop_DW.reg_stride_time = 1.5;
    } else {
      if (Human_in_Loop_DW.reg_stride_time < 0.5) {
        /* '<S108>:1:70' */
        /* '<S108>:1:71' */
        Human_in_Loop_DW.reg_stride_time = 0.5;
      }
    }

    /* '<S108>:1:74' */
    Human_in_Loop_B.mode = Human_in_Loop_DW.reg_mode;

    /* '<S108>:1:75' */
    Human_in_Loop_B.state = Human_in_Loop_DW.reg_state;

    /* '<S108>:1:76' */
    Human_in_Loop_B.stride_time = Human_in_Loop_DW.reg_stride_time;

    /* '<S108>:1:77' */
    Human_in_Loop_B.stride_timer = Human_in_Loop_DW.reg_stride_time_count;

    /* End of MATLAB Function: '<S8>/State Machine' */

    /* MATLAB Function: '<S8>/Mux1' */
    /* MATLAB Function 'State Module/Mux1': '<S107>:1' */
    /* '<S107>:1:3' */
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

    /* MATLAB Function: '<S37>/Mux1' incorporates:
     *  Constant: '<S37>/Kd'
     *  Constant: '<S37>/Kl'
     *  Constant: '<S37>/Ko'
     *  Constant: '<S37>/Kp'
     *  Constant: '<S37>/Ks'
     */
    /* MATLAB Function 'Parameter Module/Control Parameter/Mux1': '<S39>:1' */
    /* '<S39>:1:3' */
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

    /* Delay: '<S4>/Delay1' */
    Human_in_Loop_B.Delay1 = Human_in_Loop_DW.Delay1_DSTATE;

    /* MATLAB Function: '<S23>/Timer' */
    /* MATLAB Function 'Opt Module/Cosmed/Timer': '<S27>:1' */
    /* '<S27>:1:7' */
    /* '<S27>:1:8' */
    if (Human_in_Loop_P.Timer_BT_METAB_EST) {
      /* '<S27>:1:16' */
      if (Human_in_Loop_DW.time >= 120.0) {
        /* '<S27>:1:17' */
        /* '<S27>:1:18' */
        Human_in_Loop_DW.time = 0.0;

        /* '<S27>:1:19' */
        Human_in_Loop_DW.count++;

        /* '<S27>:1:20' */
        i = 1;
      } else {
        /* '<S27>:1:22' */
        Human_in_Loop_DW.time += 0.0002;

        /* '<S27>:1:23' */
        i = 0;
      }
    } else {
      /* '<S27>:1:26' */
      i = 0;

      /* '<S27>:1:27' */
      Human_in_Loop_DW.count = 0.0;

      /* '<S27>:1:28' */
      Human_in_Loop_DW.time = 120.0;
    }

    /* '<S27>:1:31' */
    Human_in_Loop_B.Trigger = i;

    /* '<S27>:1:32' */
    Human_in_Loop_B.Time = Human_in_Loop_DW.time;

    /* '<S27>:1:33' */
    Human_in_Loop_B.Count = Human_in_Loop_DW.count;

    /* End of MATLAB Function: '<S23>/Timer' */

    /* RateTransition: '<S23>/RT1' */
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

    /* End of RateTransition: '<S23>/RT1' */

    /* RateTransition: '<S23>/RT2' */
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

    /* End of RateTransition: '<S23>/RT2' */

    /* S-Function (rti_commonblock): '<S28>/S-Function1' */

    /* This comment workarounds a code generation problem */

    /* End of Outputs for S-Function (rti_commonblock): '<S28>/S-Function1' */

    /* RateTransition: '<S4>/RT1' */
    switch (Human_in_Loop_DW.RT1_write_buf_n) {
     case 0:
      Human_in_Loop_DW.RT1_read_buf_i = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT1_read_buf_i = 0;
      break;

     default:
      Human_in_Loop_DW.RT1_read_buf_i = Human_in_Loop_DW.RT1_last_buf_wr_c;
      break;
    }

    if (Human_in_Loop_DW.RT1_read_buf_i != 0) {
      Human_in_Loop_B.RT1_g = Human_in_Loop_DW.RT1_Buffer1_o;
    } else {
      Human_in_Loop_B.RT1_g = Human_in_Loop_DW.RT1_Buffer0_f;
    }

    Human_in_Loop_DW.RT1_read_buf_i = -1;

    /* End of RateTransition: '<S4>/RT1' */

    /* MATLAB Function: '<S4>/CMAES' incorporates:
     *  Constant: '<S4>/Reset'
     *  Constant: '<S4>/sigma_init'
     */
    tol = Human_in_Loop_B.Delay1;

    /* MATLAB Function 'Opt Module/CMAES': '<S22>:1' */
    if (!Human_in_Loop_DW.last_trigger_not_empty) {
      /* '<S22>:1:12' */
      /* '<S22>:1:15' */
      Human_in_Loop_DW.last_trigger = Human_in_Loop_B.Delay1;
      Human_in_Loop_DW.last_trigger_not_empty = true;

      /* '<S22>:1:22' */
      Human_in_Loop_DW.xmean.size = 1;
      Human_in_Loop_DW.xmean.data = 0.0;

      /* '<S22>:1:24' */
      Human_in_Loop_DW.xmean.data = 0.5;

      /* '<S22>:1:25' */
      Human_in_Loop_DW.sigma = Human_in_Loop_P.sigma_init_Value;

      /* '<S22>:1:32' */
      Human_in_Loop_eye_n(&Human_in_Loop_DW.B.data, Human_in_Loop_DW.B.size);

      /* '<S22>:1:34' */
      cmu = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
      tmp_data[0] = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
      Human_in_Loop_DW.C.size[0] = 1;
      Human_in_Loop_DW.C.size[1] = 1;
      tmp = 0.0;
      for (ia = 0; ia < 1; ia++) {
        tmp += tmp_data[0] * cmu;
      }

      Human_in_Loop_DW.C.data = tmp;
    }

    if ((Human_in_Loop_P.Reset_Value == 1.0) && (Human_in_Loop_DW.last_trigger ==
         0.0) && (Human_in_Loop_B.Delay1 == 1.0)) {
      /* '<S22>:1:38' */
      /* '<S22>:1:39' */
      Human_in_Loop_DW.Reset = 1.0;
    }

    if (Human_in_Loop_P.Reset_Value == 0.0) {
      /* '<S22>:1:42' */
      /* '<S22>:1:43' */
      Human_in_Loop_DW.Reset = 0.0;
    }

    if (Human_in_Loop_DW.Reset == 0.0) {
      /* '<S22>:1:46' */
      /* '<S22>:1:47' */
      Human_in_Loop_DW.count_val = 1.0;

      /* '<S22>:1:48' */
      Human_in_Loop_DW.cycle_val = 1.0;

      /* '<S22>:1:50' */
      /* '<S22>:1:51' */
      /* '<S22>:1:52' */
      Human_in_Loop_DW.xmean.size = 1;
      Human_in_Loop_DW.xmean.data = 0.5;

      /* '<S22>:1:53' */
      Human_in_Loop_DW.sigma = Human_in_Loop_P.sigma_init_Value;

      /* '<S22>:1:54' */
      Human_in_Loop_DW.arfiness.size[0] = 1;
      Human_in_Loop_DW.arfiness.size[1] = 4;

      /* '<S22>:1:55' */
      Human_in_Loop_DW.arx.size[0] = 1;
      Human_in_Loop_DW.arx.size[1] = 4;

      /* '<S22>:1:56' */
      Human_in_Loop_DW.arz.size[0] = 1;
      Human_in_Loop_DW.arz.size[1] = 4;
      Human_in_Loop_DW.arfiness.data[0] = 0.0;
      Human_in_Loop_DW.arx.data[0] = 0.0;
      Human_in_Loop_DW.arz.data[0] = 0.0;
      Human_in_Loop_DW.arfiness.data[1] = 0.0;
      Human_in_Loop_DW.arx.data[1] = 0.0;
      Human_in_Loop_DW.arz.data[1] = 0.0;
      Human_in_Loop_DW.arfiness.data[2] = 0.0;
      Human_in_Loop_DW.arx.data[2] = 0.0;
      Human_in_Loop_DW.arz.data[2] = 0.0;
      Human_in_Loop_DW.arfiness.data[3] = 0.0;
      Human_in_Loop_DW.arx.data[3] = 0.0;
      Human_in_Loop_DW.arz.data[3] = 0.0;

      /* '<S22>:1:58' */
      Human_in_Loop_DW.pc.size = 1;
      Human_in_Loop_DW.pc.data = 0.0;

      /* '<S22>:1:59' */
      Human_in_Loop_DW.ps.size = 1;
      Human_in_Loop_DW.ps.data = 0.0;

      /* '<S22>:1:60' */
      Human_in_Loop_eye_n(&Human_in_Loop_DW.B.data, Human_in_Loop_DW.B.size);

      /* '<S22>:1:61' */
      Human_in_Loop_DW.D.size[0] = 1;
      Human_in_Loop_DW.D.size[1] = 1;
      Human_in_Loop_DW.D.data = 1.0;

      /* '<S22>:1:62' */
      cmu = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
      tmp_data[0] = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
      Human_in_Loop_DW.C.size[0] = 1;
      Human_in_Loop_DW.C.size[1] = 1;
      tmp = 0.0;
      for (ia = 0; ia < 1; ia++) {
        tmp += tmp_data[0] * cmu;
      }

      Human_in_Loop_DW.C.data = tmp;

      /* '<S22>:1:64' */
      Human_in_Loop_randn(&zmean_data, &i);

      /* '<S22>:1:65' */
      Human_in_Loop_DW.arz.data[0] = zmean_data;

      /* '<S22>:1:66' */
      tmp_data[0] = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
      tmp = 0.0;
      for (ia = 0; ia < 1; ia++) {
        tmp += tmp_data[0] * Human_in_Loop_DW.arz.data[0];
      }

      tmp *= Human_in_Loop_DW.sigma;
      Human_in_Loop_DW.arx.data[0] = Human_in_Loop_DW.xmean.data + tmp;

      /* '<S22>:1:64' */
      Human_in_Loop_randn(&zmean_data, &i);

      /* '<S22>:1:65' */
      Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0]] = zmean_data;

      /* '<S22>:1:66' */
      tmp_data[0] = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
      tmp = 0.0;
      for (ia = 0; ia < 1; ia++) {
        tmp += tmp_data[0] *
          Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0]];
      }

      tmp *= Human_in_Loop_DW.sigma;
      Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0]] =
        Human_in_Loop_DW.xmean.data + tmp;

      /* '<S22>:1:64' */
      Human_in_Loop_randn(&zmean_data, &i);

      /* '<S22>:1:65' */
      Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0] << 1] = zmean_data;

      /* '<S22>:1:66' */
      tmp_data[0] = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
      tmp = 0.0;
      for (ia = 0; ia < 1; ia++) {
        tmp += Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0] << 1] *
          tmp_data[0];
      }

      tmp *= Human_in_Loop_DW.sigma;
      Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0] << 1] =
        Human_in_Loop_DW.xmean.data + tmp;

      /* '<S22>:1:64' */
      Human_in_Loop_randn(&zmean_data, &i);

      /* '<S22>:1:65' */
      Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0] * 3] = zmean_data;

      /* '<S22>:1:66' */
      tmp_data[0] = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
      tmp = 0.0;
      for (ia = 0; ia < 1; ia++) {
        tmp += Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0] * 3] *
          tmp_data[0];
      }

      tmp *= Human_in_Loop_DW.sigma;
      Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0] * 3] =
        Human_in_Loop_DW.xmean.data + tmp;

      /* '<S22>:1:71' */
      cs = Human_in_Loop_DW.arx.data[0];
      arxx_data[0] = cs;
      cs = Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0]];
      arxx_data[1] = cs;
      cs = Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0] << 1];
      arxx_data[2] = cs;
      cs = Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0] * 3];
      arxx_data[3] = cs;
      Human_in_Loop_sortrows(arxx_data, b_index_data, &weights_size);

      /* '<S22>:1:72' */
      tmp_data_0[0] = Human_in_Loop_DW.arx.data[((int32_T)b_index_data[0] - 1) *
        Human_in_Loop_DW.arx.size[0]];
      tmp_data_0[1] = Human_in_Loop_DW.arx.data[((int32_T)b_index_data[1] - 1) *
        Human_in_Loop_DW.arx.size[0]];
      tmp_data_0[2] = Human_in_Loop_DW.arx.data[((int32_T)b_index_data[2] - 1) *
        Human_in_Loop_DW.arx.size[0]];
      tmp_data_0[3] = Human_in_Loop_DW.arx.data[((int32_T)b_index_data[3] - 1) *
        Human_in_Loop_DW.arx.size[0]];
      Human_in_Loop_DW.arx.size[0] = 1;
      Human_in_Loop_DW.arx.size[1] = 4;

      /* '<S22>:1:73' */
      tmp_data_2[0] = Human_in_Loop_DW.arz.data[((int32_T)b_index_data[0] - 1) *
        Human_in_Loop_DW.arz.size[0]];
      tmp_data_2[1] = Human_in_Loop_DW.arz.data[((int32_T)b_index_data[1] - 1) *
        Human_in_Loop_DW.arz.size[0]];
      tmp_data_2[2] = Human_in_Loop_DW.arz.data[((int32_T)b_index_data[2] - 1) *
        Human_in_Loop_DW.arz.size[0]];
      tmp_data_2[3] = Human_in_Loop_DW.arz.data[((int32_T)b_index_data[3] - 1) *
        Human_in_Loop_DW.arz.size[0]];
      Human_in_Loop_DW.arz.size[0] = 1;
      Human_in_Loop_DW.arz.size[1] = 4;
      for (ia = 0; ia < 4; ia++) {
        for (i = 0; i < 1; i++) {
          Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0] * ia] =
            tmp_data_0[ia];
          Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0] * ia] =
            tmp_data_2[ia];
        }
      }

      /* '<S22>:1:75' */
      for (i = 0; i < 5; i++) {
        Human_in_Loop_B.x_this_cycle[i] = 0.0;
      }

      /* '<S22>:1:76' */
      for (i = 0; i < 5; i++) {
        Human_in_Loop_B.mean_this_generation[i] = 0.0;
      }

      /* '<S22>:1:77' */
      for (i = 0; i < 5; i++) {
        Human_in_Loop_B.variance_this_generation[i] = 0.0;
      }

      /* '<S22>:1:78' */
      memset(&Human_in_Loop_B.x_this_gengration[0], 0, 50U * sizeof(real_T));

      /* '<S22>:1:79' */
      c1 = Human_in_Loop_DW.sigma;

      /* '<S22>:1:80' */
      cmu = 1.0;
    } else {
      if ((Human_in_Loop_DW.last_trigger == 0.0) && (Human_in_Loop_B.Delay1 ==
           1.0)) {
        /* '<S22>:1:84' */
        /* '<S22>:1:86' */
        Human_in_Loop_DW.arfiness.data[(int32_T)Human_in_Loop_DW.cycle_val - 1] =
          Human_in_Loop_B.RT1_g;
        if (Human_in_Loop_DW.cycle_val == 4.0) {
          /* '<S22>:1:88' */
          /* '<S22>:1:91' */
          b_y_size[0] = 1;
          b_y_size[1] = 2;
          b_y_data[0] = 1.0;
          b_y_data[1] = 2.0;
          Human_in_Loop_log(b_y_data, b_y_size);
          weights_size = b_y_size[1];
          c_ar = b_y_size[1];
          for (ia = 0; ia < c_ar; ia++) {
            weights_data[ia] = 0.91629073187415511 - b_y_data[b_y_size[0] * ia];
          }

          /* '<S22>:1:92' */
          c1 = Human_in_Loop_sum(weights_data, &weights_size);
          c_ar = weights_size;
          for (ia = 0; ia < c_ar; ia++) {
            cc = weights_data[ia];
            cc /= c1;
            weights_data[ia] = cc;
          }

          /* '<S22>:1:93' */
          c1 = Human_in_Loop_sum(weights_data, &weights_size);
          Human_in_Loop_power_d(weights_data, &weights_size, tmp_data, &i);
          mueff = c1 * c1 / Human_in_Loop_sum(tmp_data, &i);

          /* '<S22>:1:95' */
          cc = (4.0 + mueff) / (2.0 * mueff + 5.0);

          /* '<S22>:1:96' */
          cs = (mueff + 2.0) / ((1.0 + mueff) + 5.0);

          /* '<S22>:1:97' */
          c1 = 2.0 / (5.2899999999999991 + mueff);

          /* '<S22>:1:99' */
          cmu = ((mueff - 2.0) + 1.0 / mueff) * 2.0 / (2.0 * mueff / 2.0 + 9.0);

          /* '<S22>:1:100' */
          /* '<S22>:1:103' */
          Human_in_Loop_sort(Human_in_Loop_DW.arfiness.data, iidx_data, b_y_size);

          /* '<S22>:1:103' */
          /* '<S22>:1:104' */
          b_y_data[0] = Human_in_Loop_DW.arx.data[(iidx_data[0] - 1) *
            Human_in_Loop_DW.arx.size[0]];
          b_y_data[1] = Human_in_Loop_DW.arx.data[(iidx_data[1] - 1) *
            Human_in_Loop_DW.arx.size[0]];
          if (weights_size == 1) {
            tmp = 0.0;
            for (ia = 0; ia < 2; ia++) {
              tmp += b_y_data[ia] * weights_data[ia];
            }

            Human_in_Loop_DW.xmean.size = 1;
            Human_in_Loop_DW.xmean.data = tmp;
          } else {
            Human_in_Loop_DW.xmean.size = 1;
            Human_in_Loop_DW.xmean.size = 1;
            Human_in_Loop_DW.xmean.data = 0.0;
            Human_in_Loop_DW.xmean.data = 0.0;
            i = 0;
            for (br = 0; br + 1 < 3; br++) {
              if (weights_data[br] != 0.0) {
                ia = i;
                Human_in_Loop_DW.xmean.data += weights_data[br] * b_y_data[ia];
              }

              i++;
            }
          }

          /* '<S22>:1:105' */
          b_y_data[0] = Human_in_Loop_DW.arz.data[(iidx_data[0] - 1) *
            Human_in_Loop_DW.arz.size[0]];
          b_y_data[1] = Human_in_Loop_DW.arz.data[(iidx_data[1] - 1) *
            Human_in_Loop_DW.arz.size[0]];
          if (weights_size == 1) {
            tmp = 0.0;
            for (ia = 0; ia < 2; ia++) {
              tmp += b_y_data[ia] * weights_data[ia];
            }

            zmean_data = tmp;
          } else {
            zmean_data = 0.0;
            i = 0;
            for (br = 0; br + 1 < 3; br++) {
              if (weights_data[br] != 0.0) {
                ia = i;
                zmean_data += weights_data[br] * b_y_data[ia];
              }

              i++;
            }
          }

          x = sqrt((2.0 - cs) * cs * mueff);

          /* '<S22>:1:107' */
          tmp = 1.0 - cs;
          tmp_0 = Human_in_Loop_DW.B.data * zmean_data;
          x *= tmp_0;
          Human_in_Loop_DW.ps.size = 1;
          Human_in_Loop_DW.ps.data = tmp * Human_in_Loop_DW.ps.data + x;

          /* '<S22>:1:108' */
          hsig = (fabs(Human_in_Loop_DW.ps.data) / sqrt(1.0 - rt_powd_snf(1.0 -
                    cs, 2.0 * Human_in_Loop_DW.count_val / 4.0)) /
                  0.79761904761904767 < 2.4);
          x = sqrt((2.0 - cc) * cc * mueff) * (real_T)hsig;

          /* '<S22>:1:109' */
          tmp = 1.0 - cc;
          tmp_data[0] = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
          tmp_0 = 0.0;
          for (ia = 0; ia < 1; ia++) {
            tmp_0 += tmp_data[0] * zmean_data;
          }

          x *= tmp_0;
          Human_in_Loop_DW.pc.size = 1;
          Human_in_Loop_DW.pc.data = tmp * Human_in_Loop_DW.pc.data + x;

          /* '<S22>:1:110' */
          mueff = sqrt((mueff - 1.0) / 2.0) - 1.0;
          if ((0.0 > mueff) || rtIsNaN(mueff)) {
            mueff = 0.0;
          }

          Human_in_Loop_DW.sigma *= exp(cs / ((2.0 * mueff + 1.0) + cs) * (fabs
            (Human_in_Loop_DW.ps.data) / 0.79761904761904767 - 1.0));

          /* '<S22>:1:112' */
          mueff = (1.0 - (real_T)hsig) * cc * (2.0 - cc);
          cc = (1.0 - c1) - cmu;

          /* '<S22>:1:113' */
          tmp_data[0] = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
          tmp_data_0[0] = Human_in_Loop_DW.arz.data[(iidx_data[0] - 1) *
            Human_in_Loop_DW.arz.size[0]];
          tmp_data_0[1] = Human_in_Loop_DW.arz.data[(iidx_data[1] - 1) *
            Human_in_Loop_DW.arz.size[0]];
          tmp_data_1[0] = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
          tmp_data_2[0] = Human_in_Loop_DW.arz.data[(iidx_data[0] - 1) *
            Human_in_Loop_DW.arz.size[0]];
          tmp_data_2[1] = Human_in_Loop_DW.arz.data[(iidx_data[1] - 1) *
            Human_in_Loop_DW.arz.size[0]];
          for (ia = 0; ia < 2; ia++) {
            b_y_data[ia] = 0.0;
            for (i = 0; i < 1; i++) {
              cs = b_y_data[ia];
              cs += tmp_data[0] * tmp_data_0[ia];
              b_y_data[ia] = cs;
            }

            c_data[ia] = b_y_data[ia];
            tmp_data_3[ia] = 0.0;
            for (i = 0; i < 1; i++) {
              cs = tmp_data_3[ia];
              cs += tmp_data_1[0] * tmp_data_2[ia];
              tmp_data_3[ia] = cs;
            }

            d_y_data[ia] = cmu * tmp_data_3[ia];
          }

          Human_in_Loop_diag(weights_data, &weights_size, tmp_data_0, b_y_size);
          if (b_y_size[0] == 1) {
            weights_size = b_y_size[1];
            c_ar = b_y_size[1];
            for (ia = 0; ia < c_ar; ia++) {
              tmp_data[ia] = 0.0;
              for (i = 0; i < 2; i++) {
                cmu = tmp_data[ia];
                cmu += tmp_data_0[i + ia] * d_y_data[i];
                tmp_data[ia] = cmu;
              }
            }
          } else {
            tmp_1 = (int8_T)b_y_size[1];
            weights_size = tmp_1;
            i = b_y_size[1] - 1;
            if (0 <= weights_size - 1) {
              memset(&tmp_data[0], 0, weights_size * sizeof(real_T));
            }

            for (br = 1; br - 1 <= i; br++) {
              for (ia = br; ia <= br; ia++) {
                tmp_data[ia - 1] = 0.0;
              }
            }

            br = 0;
            for (ia = 0; ia <= i; ia++) {
              c_ar = 0;
              for (c_ib = br; c_ib + 1 <= br + 2; c_ib++) {
                if (tmp_data_0[c_ib] != 0.0) {
                  c_ia = c_ar;
                  for (f_ic = ia; f_ic + 1 <= ia + 1; f_ic++) {
                    c_ia++;
                    tmp_data[f_ic] += d_y_data[c_ia - 1] * tmp_data_0[c_ib];
                  }
                }

                c_ar++;
              }

              br += 2;
            }
          }

          for (ia = 0; ia < 2; ia++) {
            d_y_data[ia] = b_y_data[ia];
          }

          for (ia = 0; ia < 2; ia++) {
            tmp_data_0[ia] = d_y_data[ia];
          }

          if (weights_size == 1) {
            for (ia = 0; ia < 1; ia++) {
              b_y_data[0] = 0.0;
              cs = b_y_data[0];
              for (i = 0; i < 1; i++) {
                cs += tmp_data[0] * c_data[0];
                b_y_data[0] = cs;
              }
            }
          } else {
            b_y_data[0] = 0.0;
            i = 0;
            br = 0;
            while (br <= 0) {
              for (ia = i; ia + 1 <= i + 2; ia++) {
                if (tmp_data_0[ia] != 0.0) {
                  c_ar = br;
                  b_y_data[0] += c_data[ia] * tmp_data[c_ar];
                }

                br++;
              }

              i += 2;
              br = 1;
            }
          }

          /* '<S22>:1:112' */
          cmu = Human_in_Loop_DW.pc.data;
          tmp = Human_in_Loop_DW.pc.data * cmu;
          cmu = mueff * Human_in_Loop_DW.C.data + tmp;
          Human_in_Loop_DW.C.size[0] = 1;
          Human_in_Loop_DW.C.size[1] = 1;
          Human_in_Loop_DW.C.data = (cc * Human_in_Loop_DW.C.data + c1 * cmu) +
            b_y_data[0];

          /* '<S22>:1:114' */
          c1 = Human_in_Loop_DW.C.data;
          Human_in_Loop_triu(&c1);
          Human_in_Loop_DW.C.size[0] = 1;
          Human_in_Loop_DW.C.size[1] = 1;
          Human_in_Loop_DW.C.data += c1;

          /* '<S22>:1:115' */
          Human_in_Loop_eig(&Human_in_Loop_DW.C.data, &Bz_data, b_y_size,
                            &Dz_data, Dz_size);

          /* '<S22>:1:116' */
          Human_in_Loop_DW.B.size[0] = 1;
          Human_in_Loop_DW.B.size[1] = 1;
          Human_in_Loop_DW.B.data = Bz_data.re;

          /* '<S22>:1:116' */
          Human_in_Loop_DW.D.size[0] = 1;
          Human_in_Loop_DW.D.size[1] = 1;
          Human_in_Loop_DW.D.data = Dz_data.re;

          /* '<S22>:1:117' */
          i = 1;
          zmean_data = Human_in_Loop_DW.D.data;
          Human_in_Loop_sqrt(&zmean_data);
          Human_in_Loop_diag(&zmean_data, &i, tmp_data_0, b_y_size);
          Human_in_Loop_DW.D.size[0] = b_y_size[0];
          Human_in_Loop_DW.D.size[1] = b_y_size[1];
          c_ar = b_y_size[0] * b_y_size[1];
          if (0 <= c_ar - 1) {
            memcpy(&Human_in_Loop_DW.D.data, &tmp_data_0[0], c_ar * sizeof
                   (real_T));
          }

          /* '<S22>:1:119' */
          Human_in_Loop_randn(&zmean_data, &i);

          /* '<S22>:1:120' */
          Human_in_Loop_DW.arz.data[0] = zmean_data;

          /* '<S22>:1:121' */
          tmp_data[0] = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
          tmp = 0.0;
          for (ia = 0; ia < 1; ia++) {
            tmp += tmp_data[0] * Human_in_Loop_DW.arz.data[0];
          }

          tmp *= Human_in_Loop_DW.sigma;
          Human_in_Loop_DW.arx.data[0] = Human_in_Loop_DW.xmean.data + tmp;

          /* '<S22>:1:119' */
          Human_in_Loop_randn(&zmean_data, &i);

          /* '<S22>:1:120' */
          Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0]] = zmean_data;

          /* '<S22>:1:121' */
          tmp_data[0] = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
          tmp = 0.0;
          for (ia = 0; ia < 1; ia++) {
            tmp += tmp_data[0] *
              Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0]];
          }

          tmp *= Human_in_Loop_DW.sigma;
          Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0]] =
            Human_in_Loop_DW.xmean.data + tmp;

          /* '<S22>:1:119' */
          Human_in_Loop_randn(&zmean_data, &i);

          /* '<S22>:1:120' */
          Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0] << 1] =
            zmean_data;

          /* '<S22>:1:121' */
          tmp_data[0] = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
          tmp = 0.0;
          for (ia = 0; ia < 1; ia++) {
            tmp += Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0] << 1] *
              tmp_data[0];
          }

          tmp *= Human_in_Loop_DW.sigma;
          Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0] << 1] =
            Human_in_Loop_DW.xmean.data + tmp;

          /* '<S22>:1:119' */
          Human_in_Loop_randn(&zmean_data, &i);

          /* '<S22>:1:120' */
          Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0] * 3] =
            zmean_data;

          /* '<S22>:1:121' */
          tmp_data[0] = Human_in_Loop_DW.B.data * Human_in_Loop_DW.D.data;
          tmp = 0.0;
          for (ia = 0; ia < 1; ia++) {
            tmp += Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0] * 3] *
              tmp_data[0];
          }

          tmp *= Human_in_Loop_DW.sigma;
          Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0] * 3] =
            Human_in_Loop_DW.xmean.data + tmp;

          /* '<S22>:1:124' */
          cs = Human_in_Loop_DW.arx.data[0];
          arxx_data[0] = cs;
          cs = Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0]];
          arxx_data[1] = cs;
          cs = Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0] << 1];
          arxx_data[2] = cs;
          cs = Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0] * 3];
          arxx_data[3] = cs;
          Human_in_Loop_sortrows(arxx_data, b_index_data, &weights_size);

          /* '<S22>:1:125' */
          tmp_data_0[0] = Human_in_Loop_DW.arx.data[((int32_T)b_index_data[0] -
            1) * Human_in_Loop_DW.arx.size[0]];
          tmp_data_0[1] = Human_in_Loop_DW.arx.data[((int32_T)b_index_data[1] -
            1) * Human_in_Loop_DW.arx.size[0]];
          tmp_data_0[2] = Human_in_Loop_DW.arx.data[((int32_T)b_index_data[2] -
            1) * Human_in_Loop_DW.arx.size[0]];
          tmp_data_0[3] = Human_in_Loop_DW.arx.data[((int32_T)b_index_data[3] -
            1) * Human_in_Loop_DW.arx.size[0]];
          Human_in_Loop_DW.arx.size[0] = 1;
          Human_in_Loop_DW.arx.size[1] = 4;

          /* '<S22>:1:126' */
          tmp_data_2[0] = Human_in_Loop_DW.arz.data[((int32_T)b_index_data[0] -
            1) * Human_in_Loop_DW.arz.size[0]];
          tmp_data_2[1] = Human_in_Loop_DW.arz.data[((int32_T)b_index_data[1] -
            1) * Human_in_Loop_DW.arz.size[0]];
          tmp_data_2[2] = Human_in_Loop_DW.arz.data[((int32_T)b_index_data[2] -
            1) * Human_in_Loop_DW.arz.size[0]];
          tmp_data_2[3] = Human_in_Loop_DW.arz.data[((int32_T)b_index_data[3] -
            1) * Human_in_Loop_DW.arz.size[0]];
          Human_in_Loop_DW.arz.size[0] = 1;
          Human_in_Loop_DW.arz.size[1] = 4;
          for (ia = 0; ia < 4; ia++) {
            for (i = 0; i < 1; i++) {
              Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0] * ia] =
                tmp_data_0[ia];
              Human_in_Loop_DW.arz.data[Human_in_Loop_DW.arz.size[0] * ia] =
                tmp_data_2[ia];
            }
          }

          /* '<S22>:1:127' */
          Human_in_Loop_DW.cycle_val = 0.0;

          /* '<S22>:1:128' */
          Human_in_Loop_DW.count_val++;
        }

        /* '<S22>:1:130' */
        Human_in_Loop_DW.cycle_val++;
      }

      /* '<S22>:1:133' */
      for (i = 0; i < 5; i++) {
        Human_in_Loop_B.x_this_cycle[i] = 0.0;
      }

      /* '<S22>:1:134' */
      for (i = 0; i < 5; i++) {
        Human_in_Loop_B.mean_this_generation[i] = 0.0;
      }

      /* '<S22>:1:135' */
      for (i = 0; i < 5; i++) {
        Human_in_Loop_B.variance_this_generation[i] = 0.0;
      }

      /* '<S22>:1:136' */
      memset(&Human_in_Loop_B.x_this_gengration[0], 0, 50U * sizeof(real_T));

      /* '<S22>:1:138' */
      Human_in_Loop_B.x_this_gengration[0] = Human_in_Loop_DW.arx.data[0];
      Human_in_Loop_B.x_this_gengration[5] =
        Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0]];
      Human_in_Loop_B.x_this_gengration[10] =
        Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0] << 1];
      Human_in_Loop_B.x_this_gengration[15] =
        Human_in_Loop_DW.arx.data[Human_in_Loop_DW.arx.size[0] * 3];

      /* '<S22>:1:139' */
      ia = (int32_T)Human_in_Loop_DW.cycle_val;
      cmu = Human_in_Loop_DW.arx.data[(ia - 1) * Human_in_Loop_DW.arx.size[0]];
      Human_in_Loop_B.x_this_cycle[0] = cmu;

      /* '<S22>:1:140' */
      Human_in_Loop_B.mean_this_generation[0] = Human_in_Loop_DW.xmean.data;

      /* '<S22>:1:141' */
      c1 = Human_in_Loop_DW.sigma;

      /* '<S22>:1:142' */
      cmu = Human_in_Loop_DW.count_val;

      /* '<S22>:1:143' */
      Human_in_Loop_B.variance_this_generation[0] = Human_in_Loop_DW.D.data;
    }

    /* '<S22>:1:147' */
    Human_in_Loop_DW.last_trigger = tol;
    Human_in_Loop_B.gengeration_count = cmu;
    Human_in_Loop_B.sigma_this_generation = c1;

    /* End of MATLAB Function: '<S4>/CMAES' */

    /* MATLAB Function: '<S38>/MATLAB Function' */
    /* MATLAB Function 'Parameter Module/Torque Parameter/MATLAB Function': '<S40>:1' */
    if (Human_in_Loop_B.x_this_cycle[0] == 0.0) {
      /* '<S40>:1:3' */
      /* '<S40>:1:5' */
      Human_in_Loop_B.parm_torque = 20.0;
    } else {
      /* '<S40>:1:8' */
      Human_in_Loop_B.parm_torque = Human_in_Loop_B.x_this_cycle[0] * 15.0 +
        30.0;
    }

    /* End of MATLAB Function: '<S38>/MATLAB Function' */

    /* MATLAB Function: '<S38>/Mux1' incorporates:
     *  Constant: '<S38>/fall_time'
     *  Constant: '<S38>/peak_time'
     *  Constant: '<S38>/rise_time'
     */
    /* MATLAB Function 'Parameter Module/Torque Parameter/Mux1': '<S41>:1' */
    /* '<S41>:1:3' */
    Human_in_Loop_B.x_k[0] = Human_in_Loop_B.parm_torque;
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

    /* S-Function (rti_commonblock): '<S9>/S-Function1' */

    /* This comment workarounds a code generation problem */

    /* End of Outputs for S-Function (rti_commonblock): '<S9>/S-Function1' */

    /* Gain: '<S3>/Gain' */
    Human_in_Loop_B.Gain_i = Human_in_Loop_P.Gain_Gain_j * Human_in_Loop_B.x[3];

    /* Delay: '<S17>/Delay3' */
    Human_in_Loop_B.Delay3 = Human_in_Loop_DW.Delay3_DSTATE;

    /* Delay: '<S17>/Delay2' */
    Human_in_Loop_B.Delay2 = Human_in_Loop_DW.Delay2_DSTATE;

    /* Delay: '<S17>/Delay7' */
    Human_in_Loop_B.Delay7 = Human_in_Loop_DW.Delay7_DSTATE;

    /* Delay: '<S17>/Delay6' */
    Human_in_Loop_B.Delay6 = Human_in_Loop_DW.Delay6_DSTATE;
  }

  /* StateSpace: '<S58>/low_pass' */
  Human_in_Loop_B.low_pass = 0.0;
  Human_in_Loop_B.low_pass += Human_in_Loop_P.low_pass_C *
    Human_in_Loop_X.low_pass_CSTATE[1];

  /* StateSpace: '<S59>/low_pass' */
  Human_in_Loop_B.low_pass_f = 0.0;
  Human_in_Loop_B.low_pass_f += Human_in_Loop_P.low_pass_C_g *
    Human_in_Loop_X.low_pass_CSTATE_o[1];

  /* StateSpace: '<S62>/low_pass' */
  Human_in_Loop_B.low_pass_l = 0.0;
  Human_in_Loop_B.low_pass_l += Human_in_Loop_P.low_pass_C_d *
    Human_in_Loop_X.low_pass_CSTATE_n[1];

  /* StateSpace: '<S63>/low_pass' */
  Human_in_Loop_B.low_pass_d = 0.0;
  Human_in_Loop_B.low_pass_d += Human_in_Loop_P.low_pass_C_i *
    Human_in_Loop_X.low_pass_CSTATE_p[1];

  /* StateSpace: '<S69>/low_pass' */
  Human_in_Loop_B.low_pass_k = 0.0;
  Human_in_Loop_B.low_pass_k += Human_in_Loop_P.low_pass_C_a *
    Human_in_Loop_X.low_pass_CSTATE_on[1];

  /* StateSpace: '<S68>/low_pass' */
  Human_in_Loop_B.low_pass_fc = 0.0;
  Human_in_Loop_B.low_pass_fc += Human_in_Loop_P.low_pass_C_c *
    Human_in_Loop_X.low_pass_CSTATE_f[1];

  /* MATLAB Function: '<S45>/Mux1' */
  Human_in_Loop_Mux1(Human_in_Loop_B.low_pass, Human_in_Loop_B.low_pass_f,
                     Human_in_Loop_B.low_pass_l, Human_in_Loop_B.low_pass_d,
                     Human_in_Loop_B.low_pass_k, Human_in_Loop_B.low_pass_fc,
                     &Human_in_Loop_B.sf_Mux1_c);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    for (i = 0; i < 6; i++) {
      /* Delay: '<S17>/Delay1' */
      Human_in_Loop_B.Delay1_p[i] = Human_in_Loop_DW.Delay1_DSTATE_c[i];

      /* Delay: '<S17>/Delay5' */
      Human_in_Loop_B.Delay5[i] = Human_in_Loop_DW.Delay5_DSTATE[i];

      /* Delay: '<S17>/Delay4' */
      Human_in_Loop_B.Delay4[i] = Human_in_Loop_DW.Delay4_DSTATE[i];

      /* Delay: '<S17>/Delay8' */
      Human_in_Loop_B.Delay8[i] = Human_in_Loop_DW.Delay8_DSTATE[i];
    }

    /* MATLAB Function: '<S17>/MATLAB Function' incorporates:
     *  Constant: '<S17>/Sample Time'
     */
    /* MATLAB Function 'ObjectiveCala/Single Cycle Analysis1/MATLAB Function': '<S21>:1' */
    /* '<S21>:1:28' */
    /* '<S21>:1:11' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_B.Max[i] = Human_in_Loop_B.Delay1_p[i];
    }

    /* '<S21>:1:12' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_B.Mean[i] = Human_in_Loop_B.Delay5[i];
    }

    /* '<S21>:1:13' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_B.RMS[i] = Human_in_Loop_B.Delay4[i];
    }

    /* '<S21>:1:14' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_B.iEMG[i] = Human_in_Loop_B.Delay8[i];
    }

    /* '<S21>:1:16' */
    tol = Human_in_Loop_B.Delay7;

    /* '<S21>:1:17' */
    Human_in_Loop_B.Cycle_Frequency = Human_in_Loop_B.Delay6;

    /* '<S21>:1:19' */
    for (i = 0; i < 6; i++) {
      /* '<S21>:1:19' */
      /* '<S21>:1:20' */
      b_y_data[0] = Human_in_Loop_DW.EMG_Memory[i];
      b_y_data[1] = Human_in_Loop_B.sf_Mux1_c.x[i];
      br = 1;
      c1 = b_y_data[0];
      if (rtIsNaN(b_y_data[0])) {
        ia = 2;
        exitg1 = false;
        while ((!exitg1) && (ia < 3)) {
          br = 2;
          if (!rtIsNaN(b_y_data[1])) {
            c1 = b_y_data[1];
            exitg1 = true;
          } else {
            ia = 3;
          }
        }
      }

      if ((br < 2) && (b_y_data[1] > c1)) {
        c1 = b_y_data[1];
      }

      Human_in_Loop_DW.EMG_Memory[i] = c1;

      /* '<S21>:1:21' */
      Human_in_Loop_DW.EMG_Memory[6 + i] += Human_in_Loop_B.sf_Mux1_c.x[i];

      /* '<S21>:1:22' */
      Human_in_Loop_DW.EMG_Memory[12 + i] += Human_in_Loop_B.sf_Mux1_c.x[i] *
        Human_in_Loop_B.sf_Mux1_c.x[i];

      /* '<S21>:1:23' */
      Human_in_Loop_DW.EMG_Memory[18 + i] += fabs(Human_in_Loop_B.sf_Mux1_c.x[i]);
    }

    if ((Human_in_Loop_B.x[1] == 1.0) && (Human_in_Loop_B.Delay2 == 0.0)) {
      /* '<S21>:1:26' */
      /* '<S21>:1:27' */
      for (ia = 0; ia < 6; ia++) {
        Human_in_Loop_B.Max[ia] = Human_in_Loop_DW.EMG_Memory[ia];
      }

      /* '<S21>:1:28' */
      for (ia = 0; ia < 6; ia++) {
        Human_in_Loop_B.Mean[ia] = Human_in_Loop_DW.EMG_Memory[6 + ia] /
          Human_in_Loop_B.Delay3;
      }

      /* '<S21>:1:29' */
      for (br = 0; br < 6; br++) {
        x = Human_in_Loop_DW.EMG_Memory[12 + br] / Human_in_Loop_B.Delay3;
        x = sqrt(x);
        x_0[br] = x;
      }

      for (ia = 0; ia < 6; ia++) {
        Human_in_Loop_B.RMS[ia] = x_0[ia];
      }

      /* '<S21>:1:30' */
      for (ia = 0; ia < 6; ia++) {
        Human_in_Loop_B.iEMG[ia] = Human_in_Loop_DW.EMG_Memory[18 + ia] /
          Human_in_Loop_B.Delay3;
      }

      /* '<S21>:1:32' */
      /* '<S21>:1:33' */
      /* '<S21>:1:34' */
      /* '<S21>:1:35' */
      for (ia = 0; ia < 6; ia++) {
        Human_in_Loop_DW.EMG_Memory[ia] = 0.0;
        Human_in_Loop_DW.EMG_Memory[6 + ia] = 0.0;
        Human_in_Loop_DW.EMG_Memory[12 + ia] = 0.0;
        Human_in_Loop_DW.EMG_Memory[18 + ia] = 0.0;
      }

      /* '<S21>:1:37' */
      tol = Human_in_Loop_B.Delay3 * Human_in_Loop_P.SingleCycleAnalysis1_Ts;

      /* '<S21>:1:38' */
      Human_in_Loop_B.Cycle_Frequency = 1.0 / tol;
    }

    Human_in_Loop_B.Cycle_Time = tol;

    /* End of MATLAB Function: '<S17>/MATLAB Function' */

    /* Outputs for Triggered SubSystem: '<S16>/Mean Calculate' incorporates:
     *  TriggerPort: '<S19>/Trigger'
     */
    if (rtmIsMajorTimeStep(Human_in_Loop_M)) {
      zcEvent = rt_ZCFcn(RISING_ZERO_CROSSING,
                         &Human_in_Loop_PrevZCX.MeanCalculate_Trig_ZCE,
                         (Human_in_Loop_B.x[1]));
      if (zcEvent != NO_ZCEVENT) {
        /* SignalConversion: '<S20>/TmpSignal ConversionAt SFunction Inport1' incorporates:
         *  MATLAB Function: '<S19>/MATLAB Function1'
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

        /* End of SignalConversion: '<S20>/TmpSignal ConversionAt SFunction Inport1' */

        /* MATLAB Function: '<S19>/MATLAB Function1' incorporates:
         *  Constant: '<S19>/Constant'
         */
        /* MATLAB Function 'ObjectiveCala/Multi Cycle Analysis1/Mean Calculate/MATLAB Function1': '<S20>:1' */
        /* '<S20>:1:5' */
        if (2.0 > Human_in_Loop_P.MultiCycleAnalysis1_N) {
          i = 0;
          br = 0;
        } else {
          i = 1;
          br = (int32_T)Human_in_Loop_P.MultiCycleAnalysis1_N;
        }

        /* '<S20>:1:5' */
        weights_size = (br - i) + 1;
        c_ar = br - i;
        for (ia = 0; ia < c_ar; ia++) {
          memcpy(&Human_in_Loop_B.tmp_data_m[ia * 26],
                 &Human_in_Loop_DW.SingleCycleData[ia * 26 + i * 26], 26U *
                 sizeof(real_T));
        }

        memcpy(&Human_in_Loop_B.tmp_data_m[(br - i) * 26],
               &Human_in_Loop_B.TmpSignalConversionAtSFunctionI[0], 26U * sizeof
               (real_T));
        for (ia = 0; ia < weights_size; ia++) {
          memcpy(&Human_in_Loop_DW.SingleCycleData[ia * 26],
                 &Human_in_Loop_B.tmp_data_m[ia * 26], 26U * sizeof(real_T));
        }

        if (1.0 > Human_in_Loop_P.MultiCycleAnalysis1_N) {
          i = 0;
        } else {
          i = (int32_T)Human_in_Loop_P.MultiCycleAnalysis1_N;
        }

        /* '<S20>:1:7' */
        tmp_size[0] = 26;
        tmp_size[1] = i;
        for (ia = 0; ia < i; ia++) {
          memcpy(&Human_in_Loop_B.tmp_data_m[ia * 26],
                 &Human_in_Loop_DW.SingleCycleData[ia * 26], 26U * sizeof(real_T));
        }

        Human_in_Loop_mean(Human_in_Loop_B.tmp_data_m, tmp_size,
                           Human_in_Loop_B.Mean_m);
      }
    }

    /* End of Outputs for SubSystem: '<S16>/Mean Calculate' */

    /* MATLAB Function: '<S3>/get loss' */
    /* MATLAB Function 'ObjectiveCala/get loss': '<S18>:1' */
    /* '<S18>:1:3' */
    Human_in_Loop_B.loss = (Human_in_Loop_B.Mean_m[14] + Human_in_Loop_B.Mean_m
      [15]) + Human_in_Loop_B.Mean_m[16];

    /* S-Function (rti_commonblock): '<S29>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* MATLAB Function: '<S4>/Opt Trigger' */
    /* MATLAB Function 'Opt Module/Opt Trigger': '<S24>:1' */
    /* '<S24>:1:11' */
    /* '<S24>:1:12' */
    /* '<S24>:1:15' */
    /* '<S24>:1:16' */
    /* '<S24>:1:17' */
    if (!Human_in_Loop_P.OptTrigger_BT_OPT) {
      /* '<S24>:1:19' */
      /* '<S24>:1:20' */
      Human_in_Loop_DW.reg_stride_num = 20.0;

      /* '<S24>:1:21' */
      Human_in_Loop_B.trig = 0.0;
    } else if ((Human_in_Loop_B.x[0] == 2.0) && (Human_in_Loop_DW.last_state ==
                0.0) && (Human_in_Loop_B.x[1] == 1.0) &&
               (Human_in_Loop_DW.reg_stride_num >=
                Human_in_Loop_P.OptTrigger_STRIDE_NUM)) {
      /* '<S24>:1:22' */
      /* '<S24>:1:23' */
      Human_in_Loop_B.trig = 1.0;

      /* '<S24>:1:24' */
      Human_in_Loop_DW.reg_stride_num = 1.0;
    } else if ((Human_in_Loop_B.x[0] == 2.0) && (Human_in_Loop_DW.last_state ==
                0.0) && (Human_in_Loop_B.x[1] == 1.0)) {
      /* '<S24>:1:25' */
      /* '<S24>:1:26' */
      Human_in_Loop_B.trig = 0.0;

      /* '<S24>:1:27' */
      Human_in_Loop_DW.reg_stride_num++;
    } else {
      /* '<S24>:1:29' */
      Human_in_Loop_B.trig = 0.0;
    }

    /* '<S24>:1:32' */
    Human_in_Loop_DW.last_state = Human_in_Loop_B.x[1];

    /* '<S24>:1:33' */
    Human_in_Loop_B.stride_num = Human_in_Loop_DW.reg_stride_num;

    /* End of MATLAB Function: '<S4>/Opt Trigger' */

    /* Delay: '<S4>/Delay' */
    Human_in_Loop_B.Delay = Human_in_Loop_DW.Delay_DSTATE;

    /* RateTransition: '<S25>/RT2' */
    switch (Human_in_Loop_DW.RT2_read_buf_aw) {
     case 0:
      Human_in_Loop_DW.RT2_write_buf_f = 1;
      break;

     case 1:
      Human_in_Loop_DW.RT2_write_buf_f = 0;
      break;

     default:
      Human_in_Loop_DW.RT2_write_buf_f = (int8_T)
        (Human_in_Loop_DW.RT2_last_buf_wr_b == 0);
      break;
    }

    if (Human_in_Loop_DW.RT2_write_buf_f != 0) {
      Human_in_Loop_DW.RT2_Buffer1_l = 0.0;
    } else {
      Human_in_Loop_DW.RT2_Buffer0_b = 0.0;
    }

    Human_in_Loop_DW.RT2_last_buf_wr_b = Human_in_Loop_DW.RT2_write_buf_f;
    Human_in_Loop_DW.RT2_write_buf_f = -1;

    /* End of RateTransition: '<S25>/RT2' */

    /* Outputs for Triggered SubSystem: '<S25>/Software Interrupt' incorporates:
     *  TriggerPort: '<S35>/Trigger'
     */
    if (rtmIsMajorTimeStep(Human_in_Loop_M)) {
      zcEvent = rt_ZCFcn(RISING_ZERO_CROSSING,
                         &Human_in_Loop_PrevZCX.SoftwareInterrupt_Trig_ZCE,
                         (Human_in_Loop_B.Delay));
      if (zcEvent != NO_ZCEVENT) {
        /* S-Function (rti_commonblock): '<S35>/S-Function1' */

        /* This comment workarounds a code generation problem */

        /* dSPACE Software Interrupt Block: <S25>/Software Interrupt */
        rtk_schedule_task(S_SOFTTASK,0);

        /* End of Outputs for S-Function (rti_commonblock): '<S35>/S-Function1' */
      }
    }

    /* End of Outputs for SubSystem: '<S25>/Software Interrupt' */

    /* RateTransition: '<S25>/RT1' */
    Human_in_Loop_B.RT1_m[0] = Human_in_Loop_B.parm_acq[0];
    Human_in_Loop_B.RT1_m[1] = Human_in_Loop_B.parm_acq[1];

    /* Gain: '<S102>/Gain' */
    Human_in_Loop_B.Gain_g = Human_in_Loop_P.uOrderTD_Ts * Human_in_Loop_B.x2k1;

    /* Sum: '<S102>/Add' */
    Human_in_Loop_B.x1k = Human_in_Loop_B.Gain_g + Human_in_Loop_B.x1k1;

    /* MATLAB Function: '<S49>/MATLAB Function' */
    /* MATLAB Function 'Sensor Data/Torque module/MATLAB Function': '<S105>:1' */
    /* '<S105>:1:13' */
    /* '<S105>:1:8' */
    memcpy(&B[0], &Human_in_Loop_DW.data[1], 14U * sizeof(real_T));
    B[14] = Human_in_Loop_B.torque;
    memcpy(&Human_in_Loop_DW.data[0], &B[0], 15U * sizeof(real_T));

    /* '<S105>:1:13' */
    for (ia = 0; ia < 30; ia++) {
      A[ia] = b_A[ia];
    }

    Human_in_Loop_xgeqp3_n(A, b_y_data, tmp_size);
    br = 0;
    tol = 15.0 * fabs(A[0]) * 2.2204460492503131E-16;
    while ((br < 2) && (!(fabs(A[15 * br + br]) <= tol))) {
      br++;
    }

    x_1[0] = 0.0;
    x_1[1] = 0.0;
    memcpy(&B[0], &Human_in_Loop_DW.data[0], 15U * sizeof(real_T));
    if (b_y_data[0] != 0.0) {
      tol = B[0];
      for (ia = 1; ia + 1 < 16; ia++) {
        tol += A[ia] * B[ia];
      }

      tol *= b_y_data[0];
      if (tol != 0.0) {
        B[0] -= tol;
        for (ia = 1; ia + 1 < 16; ia++) {
          B[ia] -= A[ia] * tol;
        }
      }
    }

    if (b_y_data[1] != 0.0) {
      tol = B[1];
      for (ia = 2; ia + 1 < 16; ia++) {
        tol += A[ia + 15] * B[ia];
      }

      tol *= b_y_data[1];
      if (tol != 0.0) {
        B[1] -= tol;
        for (ia = 2; ia + 1 < 16; ia++) {
          B[ia] -= A[ia + 15] * tol;
        }
      }
    }

    for (i = 0; i + 1 <= br; i++) {
      x_1[tmp_size[i] - 1] = B[i];
    }

    for (i = br - 1; i + 1 > 0; i--) {
      x_1[tmp_size[i] - 1] /= A[15 * i + i];
      br = 1;
      while (br <= i) {
        x_1[tmp_size[0] - 1] -= x_1[tmp_size[i] - 1] * A[15 * i];
        br = 2;
      }
    }

    /* '<S105>:1:15' */
    Human_in_Loop_B.torque_dot = x_1[1] * 5000.0;

    /* End of MATLAB Function: '<S49>/MATLAB Function' */

    /* SampleTimeMath: '<S80>/TSamp'
     *
     * About '<S80>/TSamp':
     *  y = u * K where K = 1 / ( w * Ts )
     */
    Human_in_Loop_B.TSamp = Human_in_Loop_B.Gain3_m * Human_in_Loop_P.TSamp_WtEt;

    /* UnitDelay: '<S80>/UD' */
    Human_in_Loop_B.Uk1 = Human_in_Loop_DW.UD_DSTATE;

    /* Sum: '<S80>/Diff' */
    Human_in_Loop_B.Diff = Human_in_Loop_B.TSamp - Human_in_Loop_B.Uk1;

    /* DiscreteFilter: '<S46>/0.4low1' */
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

    /* UnitDelay: '<S76>/Unit Delay' */
    Human_in_Loop_B.UnitDelay = Human_in_Loop_DW.UnitDelay_DSTATE_j;

    /* UnitDelay: '<S76>/Unit Delay1' */
    Human_in_Loop_B.UnitDelay1 = Human_in_Loop_DW.UnitDelay1_DSTATE_a;

    /* Sum: '<S76>/Add1' */
    Human_in_Loop_B.Add1 = Human_in_Loop_B.UnitDelay -
      Human_in_Loop_B.UnitDelay1;

    /* Gain: '<S76>/Gain' */
    Human_in_Loop_B.Gain_p = Human_in_Loop_P.uOrderTD_Ts_g *
      Human_in_Loop_B.Add1;

    /* Gain: '<S76>/Gain1' */
    c1 = Human_in_Loop_P.uOrderTD_T;
    tol = 1.0 / c1;
    Human_in_Loop_B.Gain1_d = tol * Human_in_Loop_B.Gain_p;

    /* Sum: '<S76>/Add2' */
    Human_in_Loop_B.Add2_i = Human_in_Loop_B.UnitDelay1 +
      Human_in_Loop_B.Gain1_d;

    /* Sum: '<S76>/Add3' */
    Human_in_Loop_B.Add3 = Human_in_Loop_B.Gain3_m - Human_in_Loop_B.Add2_i;

    /* Gain: '<S76>/Gain2' */
    c1 = Human_in_Loop_P.uOrderTD_T;
    tol = 1.0 / c1;
    Human_in_Loop_B.Gain2_i = tol * Human_in_Loop_B.Add3;

    /* Gain: '<S77>/Gain' */
    Human_in_Loop_B.Gain_io = Human_in_Loop_P.uOrderTD_Ts_n *
      Human_in_Loop_B.x2k1_k;

    /* Sum: '<S77>/Add' */
    Human_in_Loop_B.x1k_i = Human_in_Loop_B.Gain_io + Human_in_Loop_B.x1k1_m;

    /* MATLAB Function: '<S46>/Mux2' */
    /* MATLAB Function 'Sensor Data/Encoder module/Mux2': '<S87>:1' */
    /* '<S87>:1:3' */
    Human_in_Loop_B.x_d[0] = Human_in_Loop_B.u4low1;
    Human_in_Loop_B.x_d[1] = Human_in_Loop_B.Gain2_i;
    Human_in_Loop_B.x_d[2] = Human_in_Loop_B.x2k_i;

    /* S-Function (rti_commonblock): '<S83>/S-Function1' incorporates:
     *  Constant: '<S46>/VCC1'
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

    /* S-Function (rti_commonblock): '<S84>/S-Function1' incorporates:
     *  Constant: '<S46>/VCC3'
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

    /* S-Function (rti_commonblock): '<S89>/S-Function1' incorporates:
     *  Constant: '<S47>/Constant'
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

    /* S-Function (rti_commonblock): '<S50>/S-Function1' */
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

    /* S-Function (rti_commonblock): '<S51>/S-Function1' */
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

    /* S-Function (rti_commonblock): '<S52>/S-Function1' */
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

    /* S-Function (rti_commonblock): '<S53>/S-Function1' */
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

    /* S-Function (rti_commonblock): '<S54>/S-Function1' */
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

    /* S-Function (rti_commonblock): '<S55>/S-Function1' */
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

    /* Gain: '<S45>/Gain' */
    Human_in_Loop_B.Gain_k = Human_in_Loop_P.Gain_Gain_h *
      Human_in_Loop_B.SFunction1_o;

    /* Gain: '<S45>/Gain1' */
    Human_in_Loop_B.Gain1_p = Human_in_Loop_P.Gain1_Gain_b *
      Human_in_Loop_B.SFunction1_k;

    /* Gain: '<S45>/Gain2' */
    Human_in_Loop_B.Gain2_cz = Human_in_Loop_P.Gain2_Gain_k *
      Human_in_Loop_B.SFunction1_j;

    /* Gain: '<S45>/Gain3' */
    Human_in_Loop_B.Gain3_f = Human_in_Loop_P.Gain3_Gain_b *
      Human_in_Loop_B.SFunction1_d;

    /* Gain: '<S45>/Gain4' */
    Human_in_Loop_B.Gain4_h = Human_in_Loop_P.Gain4_Gain *
      Human_in_Loop_B.SFunction1_l;

    /* Gain: '<S45>/Gain5' */
    Human_in_Loop_B.Gain5 = Human_in_Loop_P.Gain5_Gain *
      Human_in_Loop_B.SFunction1_i;
  }

  /* StateSpace: '<S64>/high_pass' */
  Human_in_Loop_B.high_pass = 0.0;
  Human_in_Loop_B.high_pass += Human_in_Loop_P.high_pass_C *
    Human_in_Loop_X.high_pass_CSTATE[1];

  /* Abs: '<S64>/Abs' */
  Human_in_Loop_B.Abs = fabs(Human_in_Loop_B.high_pass);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* MATLAB Function: '<S64>/MATLAB Function' */
    Human_in_Loop_MATLABFunction(Human_in_Loop_B.Abs,
      &Human_in_Loop_B.sf_MATLABFunction_g,
      &Human_in_Loop_DW.sf_MATLABFunction_g);
  }

  /* StateSpace: '<S65>/high_pass' */
  Human_in_Loop_B.high_pass_b = 0.0;
  Human_in_Loop_B.high_pass_b += Human_in_Loop_P.high_pass_C_i *
    Human_in_Loop_X.high_pass_CSTATE_b[1];

  /* Abs: '<S65>/Abs' */
  Human_in_Loop_B.Abs_k = fabs(Human_in_Loop_B.high_pass_b);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* MATLAB Function: '<S65>/MATLAB Function' */
    Human_in_Loop_MATLABFunction(Human_in_Loop_B.Abs_k,
      &Human_in_Loop_B.sf_MATLABFunction_m,
      &Human_in_Loop_DW.sf_MATLABFunction_m);
  }

  /* StateSpace: '<S66>/high_pass' */
  Human_in_Loop_B.high_pass_d = 0.0;
  Human_in_Loop_B.high_pass_d += Human_in_Loop_P.high_pass_C_m *
    Human_in_Loop_X.high_pass_CSTATE_c[1];

  /* Abs: '<S66>/Abs' */
  Human_in_Loop_B.Abs_b = fabs(Human_in_Loop_B.high_pass_d);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* MATLAB Function: '<S66>/MATLAB Function' */
    Human_in_Loop_MATLABFunction(Human_in_Loop_B.Abs_b,
      &Human_in_Loop_B.sf_MATLABFunction_n,
      &Human_in_Loop_DW.sf_MATLABFunction_n);
  }

  /* StateSpace: '<S67>/high_pass' */
  Human_in_Loop_B.high_pass_o = 0.0;
  Human_in_Loop_B.high_pass_o += Human_in_Loop_P.high_pass_C_j *
    Human_in_Loop_X.high_pass_CSTATE_m[1];

  /* Abs: '<S67>/Abs' */
  Human_in_Loop_B.Abs_n = fabs(Human_in_Loop_B.high_pass_o);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* MATLAB Function: '<S67>/MATLAB Function' */
    Human_in_Loop_MATLABFunction(Human_in_Loop_B.Abs_n,
      &Human_in_Loop_B.sf_MATLABFunction_l,
      &Human_in_Loop_DW.sf_MATLABFunction_l);
  }

  /* StateSpace: '<S60>/high_pass' */
  Human_in_Loop_B.high_pass_l = 0.0;
  Human_in_Loop_B.high_pass_l += Human_in_Loop_P.high_pass_C_mg *
    Human_in_Loop_X.high_pass_CSTATE_a[1];

  /* Abs: '<S60>/Abs' */
  Human_in_Loop_B.Abs_m = fabs(Human_in_Loop_B.high_pass_l);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* MATLAB Function: '<S60>/MATLAB Function' */
    Human_in_Loop_MATLABFunction(Human_in_Loop_B.Abs_m,
      &Human_in_Loop_B.sf_MATLABFunction_g5,
      &Human_in_Loop_DW.sf_MATLABFunction_g5);
  }

  /* StateSpace: '<S61>/high_pass' */
  Human_in_Loop_B.high_pass_g = 0.0;
  Human_in_Loop_B.high_pass_g += Human_in_Loop_P.high_pass_C_a *
    Human_in_Loop_X.high_pass_CSTATE_my[1];

  /* Abs: '<S61>/Abs' */
  Human_in_Loop_B.Abs_f = fabs(Human_in_Loop_B.high_pass_g);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* MATLAB Function: '<S61>/MATLAB Function' */
    Human_in_Loop_MATLABFunction(Human_in_Loop_B.Abs_f,
      &Human_in_Loop_B.sf_MATLABFunction_p,
      &Human_in_Loop_DW.sf_MATLABFunction_p);

    /* MATLAB Function: '<S45>/Mux3' */
    Human_in_Loop_Mux1(Human_in_Loop_B.sf_MATLABFunction_g.y,
                       Human_in_Loop_B.sf_MATLABFunction_m.y,
                       Human_in_Loop_B.sf_MATLABFunction_n.y,
                       Human_in_Loop_B.sf_MATLABFunction_l.y,
                       Human_in_Loop_B.sf_MATLABFunction_g5.y,
                       Human_in_Loop_B.sf_MATLABFunction_p.y,
                       &Human_in_Loop_B.sf_Mux3);
  }

  /* StateSpace: '<S58>/high_pass' */
  Human_in_Loop_B.high_pass_h = 0.0;
  Human_in_Loop_B.high_pass_h += Human_in_Loop_P.high_pass_C_ja[0] *
    Human_in_Loop_X.high_pass_CSTATE_p[0];
  Human_in_Loop_B.high_pass_h += Human_in_Loop_P.high_pass_C_ja[1] *
    Human_in_Loop_X.high_pass_CSTATE_p[1];
  Human_in_Loop_B.high_pass_h += Human_in_Loop_P.high_pass_D *
    Human_in_Loop_B.Gain_k;

  /* Abs: '<S58>/Abs' */
  Human_in_Loop_B.Abs_h = fabs(Human_in_Loop_B.high_pass_h);

  /* StateSpace: '<S59>/high_pass' */
  Human_in_Loop_B.high_pass_oe = 0.0;
  Human_in_Loop_B.high_pass_oe += Human_in_Loop_P.high_pass_C_d[0] *
    Human_in_Loop_X.high_pass_CSTATE_n[0];
  Human_in_Loop_B.high_pass_oe += Human_in_Loop_P.high_pass_C_d[1] *
    Human_in_Loop_X.high_pass_CSTATE_n[1];
  Human_in_Loop_B.high_pass_oe += Human_in_Loop_P.high_pass_D_m *
    Human_in_Loop_B.Gain1_p;

  /* Abs: '<S59>/Abs' */
  Human_in_Loop_B.Abs_d = fabs(Human_in_Loop_B.high_pass_oe);

  /* StateSpace: '<S62>/high_pass' */
  Human_in_Loop_B.high_pass_ba = 0.0;
  Human_in_Loop_B.high_pass_ba += Human_in_Loop_P.high_pass_C_g[0] *
    Human_in_Loop_X.high_pass_CSTATE_pw[0];
  Human_in_Loop_B.high_pass_ba += Human_in_Loop_P.high_pass_C_g[1] *
    Human_in_Loop_X.high_pass_CSTATE_pw[1];
  Human_in_Loop_B.high_pass_ba += Human_in_Loop_P.high_pass_D_j *
    Human_in_Loop_B.Gain2_cz;

  /* Abs: '<S62>/Abs' */
  Human_in_Loop_B.Abs_j = fabs(Human_in_Loop_B.high_pass_ba);

  /* StateSpace: '<S63>/high_pass' */
  Human_in_Loop_B.high_pass_hl = 0.0;
  Human_in_Loop_B.high_pass_hl += Human_in_Loop_P.high_pass_C_h[0] *
    Human_in_Loop_X.high_pass_CSTATE_k[0];
  Human_in_Loop_B.high_pass_hl += Human_in_Loop_P.high_pass_C_h[1] *
    Human_in_Loop_X.high_pass_CSTATE_k[1];
  Human_in_Loop_B.high_pass_hl += Human_in_Loop_P.high_pass_D_a *
    Human_in_Loop_B.Gain3_f;

  /* Abs: '<S63>/Abs' */
  Human_in_Loop_B.Abs_a = fabs(Human_in_Loop_B.high_pass_hl);

  /* StateSpace: '<S68>/high_pass' */
  Human_in_Loop_B.high_pass_a = 0.0;
  Human_in_Loop_B.high_pass_a += Human_in_Loop_P.high_pass_C_dd[0] *
    Human_in_Loop_X.high_pass_CSTATE_l[0];
  Human_in_Loop_B.high_pass_a += Human_in_Loop_P.high_pass_C_dd[1] *
    Human_in_Loop_X.high_pass_CSTATE_l[1];
  Human_in_Loop_B.high_pass_a += Human_in_Loop_P.high_pass_D_n *
    Human_in_Loop_B.Gain5;

  /* Abs: '<S68>/Abs' */
  Human_in_Loop_B.Abs_dj = fabs(Human_in_Loop_B.high_pass_a);

  /* StateSpace: '<S69>/high_pass' */
  Human_in_Loop_B.high_pass_e = 0.0;
  Human_in_Loop_B.high_pass_e += Human_in_Loop_P.high_pass_C_dc[0] *
    Human_in_Loop_X.high_pass_CSTATE_h[0];
  Human_in_Loop_B.high_pass_e += Human_in_Loop_P.high_pass_C_dc[1] *
    Human_in_Loop_X.high_pass_CSTATE_h[1];
  Human_in_Loop_B.high_pass_e += Human_in_Loop_P.high_pass_D_mh *
    Human_in_Loop_B.Gain4_h;

  /* Abs: '<S69>/Abs' */
  Human_in_Loop_B.Abs_fy = fabs(Human_in_Loop_B.high_pass_e);
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S98>/S-Function1' */
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

            Human_in_Loop_B.SFunction1_o1_e = ((real_T) CAN_Sgn.SignedSgn);

            /* ...... "Accel_y" (16|16, standard signal, signed int, little endian) */
            CAN_Sgn.SgnBytes.Byte0 = CAN_Msg[2];
            CAN_Sgn.SgnBytes.Byte1 = CAN_Msg[3];
            CAN_Sgn.SignedSgn &= 0x0000FFFF;
            if (CAN_Sgn.SignedSgn >> 15) {
              CAN_Sgn.SignedSgn |= 0xFFFF0000;
            }

            Human_in_Loop_B.SFunction1_o2_b = ((real_T) CAN_Sgn.SignedSgn);

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

    /* DataTypeConversion: '<S99>/Data Type Conversion' */
    Human_in_Loop_B.DataTypeConversion = Human_in_Loop_B.SFunction1_o1_e;

    /* Gain: '<S99>/Gain' */
    Human_in_Loop_B.Gain_g5 = Human_in_Loop_P.Gain_Gain_jm *
      Human_in_Loop_B.DataTypeConversion;

    /* DataTypeConversion: '<S100>/Data Type Conversion' */
    Human_in_Loop_B.DataTypeConversion_b = Human_in_Loop_B.SFunction1_o2_b;

    /* Gain: '<S100>/Gain' */
    Human_in_Loop_B.Gain_ik = Human_in_Loop_P.Gain_Gain_i *
      Human_in_Loop_B.DataTypeConversion_b;

    /* DataTypeConversion: '<S101>/Data Type Conversion' */
    Human_in_Loop_B.DataTypeConversion_o = Human_in_Loop_B.SFunction1_o3;

    /* Gain: '<S101>/Gain' */
    Human_in_Loop_B.Gain_lw = Human_in_Loop_P.Gain_Gain_h0 *
      Human_in_Loop_B.DataTypeConversion_o;

    /* S-Function (rti_commonblock): '<S91>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* S-Function (rti_commonblock): '<S92>/S-Function1' */
    /* This comment workarounds a code generation problem */

    /* S-Function (rti_commonblock): '<S96>/S-Function1' */
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
          Human_in_Loop_B.SFunction1_o2_ks = (real_T)
            can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64]->processed;
          Human_in_Loop_B.SFunction1_o3_l = (real_T)
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
        Human_in_Loop_B.SFunction1_o2_ks = 0.0;
      }
    }
  }

  /* Sin: '<S93>/Sine Wave' */
  Human_in_Loop_B.SineWave = sin(Human_in_Loop_P.SineWave_Freq *
    Human_in_Loop_M->Timing.t[0] + Human_in_Loop_P.SineWave_Phase) *
    Human_in_Loop_P.SineWave_Amp + Human_in_Loop_P.SineWave_Bias;
  if (rtmIsMajorTimeStep(Human_in_Loop_M) &&
      Human_in_Loop_M->Timing.TaskCounters.TID[1] == 0) {
    /* S-Function (rti_commonblock): '<S97>/S-Function1' */
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
        Human_in_Loop_B.SFunction1_o1_m = (real_T)
          can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->processed;
        Human_in_Loop_B.SFunction1_o2_n = (real_T)
          can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->timestamp;
        Human_in_Loop_B.SFunction1_o3_b = (real_T)
          can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64]->deltatime;
        Human_in_Loop_B.SFunction1_o4_i = (real_T)
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

    /* S-Function (rti_commonblock): '<S94>/S-Function1' */
    /* This comment workarounds a code generation problem */
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
    /* Update for DiscreteFilter: '<S49>/0.4low2' */
    Human_in_Loop_DW.u4low2_states[2] = Human_in_Loop_DW.u4low2_states[1];
    Human_in_Loop_DW.u4low2_states[1] = Human_in_Loop_DW.u4low2_states[0];
    Human_in_Loop_DW.u4low2_states[0] = Human_in_Loop_DW.u4low2_tmp;

    /* Update for UnitDelay: '<S102>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE = Human_in_Loop_B.x2k;

    /* Update for UnitDelay: '<S102>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE = Human_in_Loop_B.x1k;

    /* Update for UnitDelay: '<S102>/Unit Delay2' */
    Human_in_Loop_DW.UnitDelay2_DSTATE = Human_in_Loop_B.torque;

    /* Update for UnitDelay: '<S77>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE_h = Human_in_Loop_B.x2k_i;

    /* Update for UnitDelay: '<S77>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_e = Human_in_Loop_B.x1k_i;

    /* Update for UnitDelay: '<S77>/Unit Delay2' */
    Human_in_Loop_DW.UnitDelay2_DSTATE_a = Human_in_Loop_B.Gain3_m;

    /* Update for Delay: '<S4>/Delay1' */
    Human_in_Loop_DW.Delay1_DSTATE = Human_in_Loop_B.Trigger;

    /* Update for Delay: '<S17>/Delay3' */
    Human_in_Loop_DW.Delay3_DSTATE = Human_in_Loop_B.Gain_i;

    /* Update for Delay: '<S17>/Delay2' */
    Human_in_Loop_DW.Delay2_DSTATE = Human_in_Loop_B.x[1];

    /* Update for Delay: '<S17>/Delay7' */
    Human_in_Loop_DW.Delay7_DSTATE = Human_in_Loop_B.Cycle_Time;

    /* Update for Delay: '<S17>/Delay6' */
    Human_in_Loop_DW.Delay6_DSTATE = Human_in_Loop_B.Cycle_Frequency;

    /* Update for Delay: '<S17>/Delay1' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_DW.Delay1_DSTATE_c[i] = Human_in_Loop_B.Max[i];
    }

    /* End of Update for Delay: '<S17>/Delay1' */

    /* Update for Delay: '<S17>/Delay5' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_DW.Delay5_DSTATE[i] = Human_in_Loop_B.Mean[i];
    }

    /* End of Update for Delay: '<S17>/Delay5' */

    /* Update for Delay: '<S17>/Delay4' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_DW.Delay4_DSTATE[i] = Human_in_Loop_B.RMS[i];
    }

    /* End of Update for Delay: '<S17>/Delay4' */

    /* Update for Delay: '<S17>/Delay8' */
    for (i = 0; i < 6; i++) {
      Human_in_Loop_DW.Delay8_DSTATE[i] = Human_in_Loop_B.iEMG[i];
    }

    /* End of Update for Delay: '<S17>/Delay8' */

    /* Update for Delay: '<S4>/Delay' */
    Human_in_Loop_DW.Delay_DSTATE = 0.0;

    /* Update for UnitDelay: '<S80>/UD' */
    Human_in_Loop_DW.UD_DSTATE = Human_in_Loop_B.TSamp;

    /* Update for DiscreteFilter: '<S46>/0.4low1' */
    Human_in_Loop_DW.u4low1_states[2] = Human_in_Loop_DW.u4low1_states[1];
    Human_in_Loop_DW.u4low1_states[1] = Human_in_Loop_DW.u4low1_states[0];
    Human_in_Loop_DW.u4low1_states[0] = Human_in_Loop_DW.u4low1_tmp;

    /* Update for UnitDelay: '<S76>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_j = Human_in_Loop_B.Gain3_m;

    /* Update for UnitDelay: '<S76>/Unit Delay1' */
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

  /* Derivatives for StateSpace: '<S58>/low_pass' */
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

  /* Derivatives for StateSpace: '<S59>/low_pass' */
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

  /* Derivatives for StateSpace: '<S62>/low_pass' */
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

  /* Derivatives for StateSpace: '<S63>/low_pass' */
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

  /* Derivatives for StateSpace: '<S69>/low_pass' */
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

  /* Derivatives for StateSpace: '<S68>/low_pass' */
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

  /* Derivatives for StateSpace: '<S64>/high_pass' */
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

  /* Derivatives for StateSpace: '<S65>/high_pass' */
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

  /* Derivatives for StateSpace: '<S66>/high_pass' */
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

  /* Derivatives for StateSpace: '<S67>/high_pass' */
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

  /* Derivatives for StateSpace: '<S60>/high_pass' */
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

  /* Derivatives for StateSpace: '<S61>/high_pass' */
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

  /* Derivatives for StateSpace: '<S58>/high_pass' */
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

  /* Derivatives for StateSpace: '<S59>/high_pass' */
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

  /* Derivatives for StateSpace: '<S62>/high_pass' */
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

  /* Derivatives for StateSpace: '<S63>/high_pass' */
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

  /* Derivatives for StateSpace: '<S68>/high_pass' */
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

  /* Derivatives for StateSpace: '<S69>/high_pass' */
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

    /* Start for RateTransition: '<S23>/RT1' */
    Human_in_Loop_B.RT1_c = Human_in_Loop_P.RT1_InitialCondition_c;

    /* Start for RateTransition: '<S23>/RT2' */
    Human_in_Loop_B.RT2_g = Human_in_Loop_P.RT2_InitialCondition_i;

    /* Start for RateTransition: '<S4>/RT1' */
    Human_in_Loop_B.RT1_g = Human_in_Loop_P.RT1_InitialCondition_cs;

    /* Start for RateTransition: '<Root>/RT3' */
    Human_in_Loop_B.RT3[0] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_B.RT3[1] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_B.RT3[2] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_B.RT3[3] = Human_in_Loop_P.RT3_InitialCondition;

    /* Start for Triggered SubSystem: '<S16>/Mean Calculate' */
    /* Start for DataStoreMemory: '<S19>/Data Store Memory' */
    memcpy(&Human_in_Loop_DW.SingleCycleData[0],
           &Human_in_Loop_P.DataStoreMemory_InitialValue[0], 2600U * sizeof
           (real_T));

    /* End of Start for SubSystem: '<S16>/Mean Calculate' */

    /* Start for DataStoreMemory: '<S17>/Data Store Memory' */
    memcpy(&Human_in_Loop_DW.EMG_Memory[0],
           &Human_in_Loop_P.DataStoreMemory_InitialValue_o[0], 24U * sizeof
           (real_T));

    /* Start for RateTransition: '<S25>/RT2' */
    Human_in_Loop_B.RT2_o = Human_in_Loop_P.RT2_InitialCondition_p;

    /* Start for DataStoreMemory: '<Root>/Data Store Memory' */
    memcpy(&Human_in_Loop_DW.TorqueMem[0],
           &Human_in_Loop_P.DataStoreMemory_InitialValue_l[0], 4400U * sizeof
           (real_T));
  }

  Human_in_Loop_PrevZCX.MeanCalculate_Trig_ZCE = UNINITIALIZED_ZCSIG;
  Human_in_Loop_PrevZCX.SoftwareInterrupt_Trig_ZCE = UNINITIALIZED_ZCSIG;

  {
    uint32_T r;
    int32_T i;

    /* InitializeConditions for DiscreteFilter: '<S49>/0.4low2' */
    Human_in_Loop_DW.u4low2_states[0] = Human_in_Loop_P.u4low2_InitialStates;
    Human_in_Loop_DW.u4low2_states[1] = Human_in_Loop_P.u4low2_InitialStates;
    Human_in_Loop_DW.u4low2_states[2] = Human_in_Loop_P.u4low2_InitialStates;

    /* InitializeConditions for UnitDelay: '<S102>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE =
      Human_in_Loop_P.UnitDelay1_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S102>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE =
      Human_in_Loop_P.UnitDelay_InitialCondition;

    /* InitializeConditions for UnitDelay: '<S102>/Unit Delay2' */
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

    /* InitializeConditions for UnitDelay: '<S77>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE_h =
      Human_in_Loop_P.UnitDelay1_InitialCondition_j;

    /* InitializeConditions for UnitDelay: '<S77>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_e =
      Human_in_Loop_P.UnitDelay_InitialCondition_j;

    /* InitializeConditions for UnitDelay: '<S77>/Unit Delay2' */
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

    /* InitializeConditions for Delay: '<S4>/Delay1' */
    Human_in_Loop_DW.Delay1_DSTATE = Human_in_Loop_P.Delay1_InitialCondition;

    /* InitializeConditions for RateTransition: '<S23>/RT1' */
    Human_in_Loop_DW.RT1_Buffer0_e = Human_in_Loop_P.RT1_InitialCondition_c;
    Human_in_Loop_DW.RT1_write_buf_f = -1;
    Human_in_Loop_DW.RT1_read_buf_p = -1;

    /* InitializeConditions for RateTransition: '<S23>/RT2' */
    Human_in_Loop_DW.RT2_Buffer0_i = Human_in_Loop_P.RT2_InitialCondition_i;
    Human_in_Loop_DW.RT2_write_buf_g = -1;
    Human_in_Loop_DW.RT2_read_buf_a = -1;

    /* InitializeConditions for RateTransition: '<S4>/RT1' */
    Human_in_Loop_DW.RT1_Buffer0_f = Human_in_Loop_P.RT1_InitialCondition_cs;
    Human_in_Loop_DW.RT1_write_buf_n = -1;
    Human_in_Loop_DW.RT1_read_buf_i = -1;

    /* InitializeConditions for RateTransition: '<Root>/RT3' */
    Human_in_Loop_DW.RT3_Buffer0[0] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_DW.RT3_Buffer0[1] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_DW.RT3_Buffer0[2] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_DW.RT3_Buffer0[3] = Human_in_Loop_P.RT3_InitialCondition;
    Human_in_Loop_DW.RT3_write_buf = -1;
    Human_in_Loop_DW.RT3_read_buf = -1;

    /* InitializeConditions for Delay: '<S17>/Delay3' */
    Human_in_Loop_DW.Delay3_DSTATE = Human_in_Loop_P.Delay3_InitialCondition;

    /* InitializeConditions for Delay: '<S17>/Delay2' */
    Human_in_Loop_DW.Delay2_DSTATE = Human_in_Loop_P.Delay2_InitialCondition;

    /* InitializeConditions for Delay: '<S17>/Delay7' */
    Human_in_Loop_DW.Delay7_DSTATE = Human_in_Loop_P.Delay7_InitialCondition;

    /* InitializeConditions for Delay: '<S17>/Delay6' */
    Human_in_Loop_DW.Delay6_DSTATE = Human_in_Loop_P.Delay6_InitialCondition;

    /* InitializeConditions for StateSpace: '<S58>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE[0] =
      Human_in_Loop_P.low_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S59>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_o[0] =
      Human_in_Loop_P.low_pass_InitialCondition_e;

    /* InitializeConditions for StateSpace: '<S62>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_n[0] =
      Human_in_Loop_P.low_pass_InitialCondition_ew;

    /* InitializeConditions for StateSpace: '<S63>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_p[0] =
      Human_in_Loop_P.low_pass_InitialCondition_c;

    /* InitializeConditions for StateSpace: '<S69>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_on[0] =
      Human_in_Loop_P.low_pass_InitialCondition_l;

    /* InitializeConditions for StateSpace: '<S68>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_f[0] =
      Human_in_Loop_P.low_pass_InitialCondition_i;

    /* InitializeConditions for StateSpace: '<S58>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE[1] =
      Human_in_Loop_P.low_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S59>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_o[1] =
      Human_in_Loop_P.low_pass_InitialCondition_e;

    /* InitializeConditions for StateSpace: '<S62>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_n[1] =
      Human_in_Loop_P.low_pass_InitialCondition_ew;

    /* InitializeConditions for StateSpace: '<S63>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_p[1] =
      Human_in_Loop_P.low_pass_InitialCondition_c;

    /* InitializeConditions for StateSpace: '<S69>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_on[1] =
      Human_in_Loop_P.low_pass_InitialCondition_l;

    /* InitializeConditions for StateSpace: '<S68>/low_pass' */
    Human_in_Loop_X.low_pass_CSTATE_f[1] =
      Human_in_Loop_P.low_pass_InitialCondition_i;
    for (i = 0; i < 6; i++) {
      /* InitializeConditions for Delay: '<S17>/Delay1' */
      Human_in_Loop_DW.Delay1_DSTATE_c[i] =
        Human_in_Loop_P.Delay1_InitialCondition_o;

      /* InitializeConditions for Delay: '<S17>/Delay5' */
      Human_in_Loop_DW.Delay5_DSTATE[i] =
        Human_in_Loop_P.Delay5_InitialCondition;

      /* InitializeConditions for Delay: '<S17>/Delay4' */
      Human_in_Loop_DW.Delay4_DSTATE[i] =
        Human_in_Loop_P.Delay4_InitialCondition;

      /* InitializeConditions for Delay: '<S17>/Delay8' */
      Human_in_Loop_DW.Delay8_DSTATE[i] =
        Human_in_Loop_P.Delay8_InitialCondition;
    }

    /* InitializeConditions for Delay: '<S4>/Delay' */
    Human_in_Loop_DW.Delay_DSTATE = Human_in_Loop_P.Delay_InitialCondition;

    /* InitializeConditions for RateTransition: '<S25>/RT2' */
    Human_in_Loop_DW.RT2_Buffer0_b = Human_in_Loop_P.RT2_InitialCondition_p;
    Human_in_Loop_DW.RT2_write_buf_f = -1;
    Human_in_Loop_DW.RT2_read_buf_aw = -1;

    /* InitializeConditions for UnitDelay: '<S80>/UD' */
    Human_in_Loop_DW.UD_DSTATE = Human_in_Loop_P.DiscreteDerivative_ICPrevScaled;

    /* InitializeConditions for DiscreteFilter: '<S46>/0.4low1' */
    Human_in_Loop_DW.u4low1_states[0] = Human_in_Loop_P.u4low1_InitialStates;
    Human_in_Loop_DW.u4low1_states[1] = Human_in_Loop_P.u4low1_InitialStates;
    Human_in_Loop_DW.u4low1_states[2] = Human_in_Loop_P.u4low1_InitialStates;

    /* InitializeConditions for UnitDelay: '<S76>/Unit Delay' */
    Human_in_Loop_DW.UnitDelay_DSTATE_j =
      Human_in_Loop_P.UnitDelay_InitialCondition_g;

    /* InitializeConditions for UnitDelay: '<S76>/Unit Delay1' */
    Human_in_Loop_DW.UnitDelay1_DSTATE_a =
      Human_in_Loop_P.UnitDelay1_InitialCondition_i;

    /* InitializeConditions for StateSpace: '<S64>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE[0] =
      Human_in_Loop_P.high_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S65>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_b[0] =
      Human_in_Loop_P.high_pass_InitialCondition_g;

    /* InitializeConditions for StateSpace: '<S66>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_c[0] =
      Human_in_Loop_P.high_pass_InitialCondition_n;

    /* InitializeConditions for StateSpace: '<S67>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_m[0] =
      Human_in_Loop_P.high_pass_InitialCondition_b;

    /* InitializeConditions for StateSpace: '<S60>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_a[0] =
      Human_in_Loop_P.high_pass_InitialCondition_bd;

    /* InitializeConditions for StateSpace: '<S61>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_my[0] =
      Human_in_Loop_P.high_pass_InitialCondition_m;

    /* InitializeConditions for StateSpace: '<S64>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE[1] =
      Human_in_Loop_P.high_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S65>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_b[1] =
      Human_in_Loop_P.high_pass_InitialCondition_g;

    /* InitializeConditions for StateSpace: '<S66>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_c[1] =
      Human_in_Loop_P.high_pass_InitialCondition_n;

    /* InitializeConditions for StateSpace: '<S67>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_m[1] =
      Human_in_Loop_P.high_pass_InitialCondition_b;

    /* InitializeConditions for StateSpace: '<S60>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_a[1] =
      Human_in_Loop_P.high_pass_InitialCondition_bd;

    /* InitializeConditions for StateSpace: '<S61>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_my[1] =
      Human_in_Loop_P.high_pass_InitialCondition_m;

    /* InitializeConditions for StateSpace: '<S64>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE[2] =
      Human_in_Loop_P.high_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S65>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_b[2] =
      Human_in_Loop_P.high_pass_InitialCondition_g;

    /* InitializeConditions for StateSpace: '<S66>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_c[2] =
      Human_in_Loop_P.high_pass_InitialCondition_n;

    /* InitializeConditions for StateSpace: '<S67>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_m[2] =
      Human_in_Loop_P.high_pass_InitialCondition_b;

    /* InitializeConditions for StateSpace: '<S60>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_a[2] =
      Human_in_Loop_P.high_pass_InitialCondition_bd;

    /* InitializeConditions for StateSpace: '<S61>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_my[2] =
      Human_in_Loop_P.high_pass_InitialCondition_m;

    /* InitializeConditions for StateSpace: '<S64>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE[3] =
      Human_in_Loop_P.high_pass_InitialCondition;

    /* InitializeConditions for StateSpace: '<S65>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_b[3] =
      Human_in_Loop_P.high_pass_InitialCondition_g;

    /* InitializeConditions for StateSpace: '<S66>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_c[3] =
      Human_in_Loop_P.high_pass_InitialCondition_n;

    /* InitializeConditions for StateSpace: '<S67>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_m[3] =
      Human_in_Loop_P.high_pass_InitialCondition_b;

    /* InitializeConditions for StateSpace: '<S60>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_a[3] =
      Human_in_Loop_P.high_pass_InitialCondition_bd;

    /* InitializeConditions for StateSpace: '<S61>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_my[3] =
      Human_in_Loop_P.high_pass_InitialCondition_m;

    /* InitializeConditions for StateSpace: '<S58>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_p[0] =
      Human_in_Loop_P.high_pass_InitialCondition_p;

    /* InitializeConditions for StateSpace: '<S59>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_n[0] =
      Human_in_Loop_P.high_pass_InitialCondition_n4;

    /* InitializeConditions for StateSpace: '<S62>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_pw[0] =
      Human_in_Loop_P.high_pass_InitialCondition_m4;

    /* InitializeConditions for StateSpace: '<S63>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_k[0] =
      Human_in_Loop_P.high_pass_InitialCondition_gb;

    /* InitializeConditions for StateSpace: '<S68>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_l[0] =
      Human_in_Loop_P.high_pass_InitialCondition_c;

    /* InitializeConditions for StateSpace: '<S69>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_h[0] =
      Human_in_Loop_P.high_pass_InitialCondition_gw;

    /* InitializeConditions for StateSpace: '<S58>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_p[1] =
      Human_in_Loop_P.high_pass_InitialCondition_p;

    /* InitializeConditions for StateSpace: '<S59>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_n[1] =
      Human_in_Loop_P.high_pass_InitialCondition_n4;

    /* InitializeConditions for StateSpace: '<S62>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_pw[1] =
      Human_in_Loop_P.high_pass_InitialCondition_m4;

    /* InitializeConditions for StateSpace: '<S63>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_k[1] =
      Human_in_Loop_P.high_pass_InitialCondition_gb;

    /* InitializeConditions for StateSpace: '<S68>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_l[1] =
      Human_in_Loop_P.high_pass_InitialCondition_c;

    /* InitializeConditions for StateSpace: '<S69>/high_pass' */
    Human_in_Loop_X.high_pass_CSTATE_h[1] =
      Human_in_Loop_P.high_pass_InitialCondition_gw;

    /* SystemInitialize for MATLAB Function: '<S49>/Data process' */
    Human_in_Loop_DW.torque_zero = 0.0;

    /* SystemInitialize for MATLAB Function: '<S46>/Data process' */
    Human_in_Loop_DW.angle_zero_f = 0.0;

    /* SystemInitialize for MATLAB Function: '<S46>/Data process1' */
    Human_in_Loop_DW.angle_zero = 0.0;

    /* SystemInitialize for MATLAB Function: '<S47>/FootSwitch Filter' */
    Human_in_Loop_DW.foot_state = 0.0;
    Human_in_Loop_DW.filter_time = 0.0;

    /* SystemInitialize for MATLAB Function: '<S8>/State Machine' */
    Human_in_Loop_DW.reg_stride_time = 1.0;
    Human_in_Loop_DW.reg_stride_time_count = 0.0;
    Human_in_Loop_DW.reg_mode = 1.0;
    Human_in_Loop_DW.reg_state = 1.0;
    Human_in_Loop_DW.bt_run = 0.0;
    Human_in_Loop_DW.reg_last_switch = 1.0;

    /* SystemInitialize for MATLAB Function: '<S23>/Timer' */
    Human_in_Loop_DW.count = 0.0;
    Human_in_Loop_DW.time = 120.0;

    /* SystemInitialize for S-Function (rti_commonblock): '<S28>/S-Function1' incorporates:
     *  SubSystem: '<S26>/Serial Decoding System'
     */
    Human_SerialDecodingSystem_Init();

    /* End of SystemInitialize for S-Function (rti_commonblock): '<S28>/S-Function1' */

    /* SystemInitialize for MATLAB Function: '<S4>/CMAES' */
    Human_in_Loop_DW.last_trigger_not_empty = false;
    memset(&Human_in_Loop_DW.state_l[0], 0, 625U * sizeof(uint32_T));
    r = 5489U;
    Human_in_Loop_DW.state_l[0] = 5489U;
    for (i = 0; i < 623; i++) {
      r = ((r >> 30U ^ r) * 1812433253U + i) + 1U;
      Human_in_Loop_DW.state_l[i + 1] = r;
    }

    Human_in_Loop_DW.state_l[624] = 624U;
    Human_in_Loop_DW.count_val = 1.0;
    Human_in_Loop_DW.cycle_val = 1.0;
    Human_in_Loop_DW.Reset = 0.0;
    Human_in_Loop_DW.arfiness.size[0] = 1;
    Human_in_Loop_DW.arfiness.size[1] = 4;
    Human_in_Loop_DW.arx.size[0] = 1;
    Human_in_Loop_DW.arx.size[1] = 4;
    Human_in_Loop_DW.arz.size[0] = 1;
    Human_in_Loop_DW.arz.size[1] = 4;
    Human_in_Loop_DW.arfiness.data[0] = 0.0;
    Human_in_Loop_DW.arx.data[0] = 0.0;
    Human_in_Loop_DW.arz.data[0] = 0.0;
    Human_in_Loop_DW.arfiness.data[1] = 0.0;
    Human_in_Loop_DW.arx.data[1] = 0.0;
    Human_in_Loop_DW.arz.data[1] = 0.0;
    Human_in_Loop_DW.arfiness.data[2] = 0.0;
    Human_in_Loop_DW.arx.data[2] = 0.0;
    Human_in_Loop_DW.arz.data[2] = 0.0;
    Human_in_Loop_DW.arfiness.data[3] = 0.0;
    Human_in_Loop_DW.arx.data[3] = 0.0;
    Human_in_Loop_DW.arz.data[3] = 0.0;
    Human_in_Loop_DW.pc.size = 1;
    Human_in_Loop_DW.pc.data = 0.0;
    Human_in_Loop_DW.ps.size = 1;
    Human_in_Loop_DW.ps.data = 0.0;
    Human_in_Loop_DW.D.size[0] = 1;
    Human_in_Loop_DW.D.size[1] = 1;
    Human_in_Loop_DW.D.data = 0.0;

    /* End of SystemInitialize for MATLAB Function: '<S4>/CMAES' */

    /* SystemInitialize for S-Function (rti_commonblock): '<S9>/S-Function1' incorporates:
     *  SubSystem: '<Root>/Control Module'
     */
    Human_in_Loo_ControlModule_Init();

    /* End of SystemInitialize for S-Function (rti_commonblock): '<S9>/S-Function1' */

    /* SystemInitialize for Triggered SubSystem: '<S16>/Mean Calculate' */
    /* SystemInitialize for Outport: '<S19>/Mean' */
    for (i = 0; i < 26; i++) {
      Human_in_Loop_B.Mean_m[i] = Human_in_Loop_P.Mean_Y0;
    }

    /* End of SystemInitialize for Outport: '<S19>/Mean' */
    /* End of SystemInitialize for SubSystem: '<S16>/Mean Calculate' */

    /* SystemInitialize for MATLAB Function: '<S4>/Opt Trigger' */
    Human_in_Loop_DW.last_state = 0.0;
    Human_in_Loop_DW.reg_stride_num = 20.0;

    /* SystemInitialize for Triggered SubSystem: '<S25>/Software Interrupt' */

    /* SystemInitialize for S-Function (rti_commonblock): '<S35>/S-Function1' incorporates:
     *  SubSystem: '<S25>/Bayeisan Opt'
     */
    Human_in_Loop_BayeisanOpt_Init();

    /* End of SystemInitialize for S-Function (rti_commonblock): '<S35>/S-Function1' */

    /* End of SystemInitialize for SubSystem: '<S25>/Software Interrupt' */

    /* SystemInitialize for MATLAB Function: '<S49>/MATLAB Function' */
    for (i = 0; i < 15; i++) {
      Human_in_Loop_DW.data[i] = 1.0;
    }

    /* End of SystemInitialize for MATLAB Function: '<S49>/MATLAB Function' */

    /* SystemInitialize for MATLAB Function: '<S64>/MATLAB Function' */
    Human_in_Lo_MATLABFunction_Init(&Human_in_Loop_DW.sf_MATLABFunction_g);

    /* SystemInitialize for MATLAB Function: '<S65>/MATLAB Function' */
    Human_in_Lo_MATLABFunction_Init(&Human_in_Loop_DW.sf_MATLABFunction_m);

    /* SystemInitialize for MATLAB Function: '<S66>/MATLAB Function' */
    Human_in_Lo_MATLABFunction_Init(&Human_in_Loop_DW.sf_MATLABFunction_n);

    /* SystemInitialize for MATLAB Function: '<S67>/MATLAB Function' */
    Human_in_Lo_MATLABFunction_Init(&Human_in_Loop_DW.sf_MATLABFunction_l);

    /* SystemInitialize for MATLAB Function: '<S60>/MATLAB Function' */
    Human_in_Lo_MATLABFunction_Init(&Human_in_Loop_DW.sf_MATLABFunction_g5);

    /* SystemInitialize for MATLAB Function: '<S61>/MATLAB Function' */
    Human_in_Lo_MATLABFunction_Init(&Human_in_Loop_DW.sf_MATLABFunction_p);
  }
}

/* Model terminate function */
void Human_in_Loop_terminate(void)
{
  /* Terminate for S-Function (rti_commonblock): '<S81>/S-Function1' */

  /* --- Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL1 --- */
  /* --- [RTIEMC, Encoder] - DIO class: 2 - Unit: 1 - Port: 1 - Channel: 1 --- */
  {
    /* Deactivates encoder interface functionality */
    DioCl2EncoderIn_stop(pRTIEmcEncoder_Unit_1_DioCl_2_Port_1_Ch1);
  }

  /* Terminate for S-Function (rti_commonblock): '<S82>/S-Function1' */

  /* --- Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL3 --- */
  /* --- [RTIEMC, Encoder] - DIO class: 2 - Unit: 3 - Port: 1 - Channel: 5 --- */
  {
    /* Deactivates encoder interface functionality */
    DioCl2EncoderIn_stop(pRTIEmcEncoder_Unit_3_DioCl_2_Port_1_Ch5);
  }

  /* Terminate for S-Function (rti_commonblock): '<S9>/S-Function1' incorporates:
   *  SubSystem: '<Root>/Control Module'
   */
  Human_in_Loo_ControlModule_Term();

  /* End of Terminate for S-Function (rti_commonblock): '<S9>/S-Function1' */

  /* Terminate for S-Function (rti_commonblock): '<S29>/S-Function1' */

  /* dSPACE I/O Board DS1202SER #1 Unit:GENSER Group:SETUP */
  dsser_disable(rtiDS1202SER_B1_Ser[0]);
  dsser_fifo_reset(rtiDS1202SER_B1_Ser[0]);

  /* Terminate for S-Function (rti_commonblock): '<S83>/S-Function1' incorporates:
   *  Constant: '<S46>/VCC1'
   */

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup1 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 11 --- */

  /* disable digital output channel 11-11 on port 3 *
   * (set to high-impedance), when the simulation terminates       */
  DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_3_Ch_11,
    DIO_CLASS1_HIGH_Z_ON);
  DioCl1DigOut_write(pRTIDioC1DigOut_Port_3_Ch_11);

  /* Terminate for S-Function (rti_commonblock): '<S84>/S-Function1' incorporates:
   *  Constant: '<S46>/VCC3'
   */

  /* --- Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup3 --- */
  /* --- [RTI120X, BITOUT] - Port: 3 - Channel: 13 --- */

  /* disable digital output channel 13-13 on port 3 *
   * (set to high-impedance), when the simulation terminates       */
  DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_3_Ch_13,
    DIO_CLASS1_HIGH_Z_ON);
  DioCl1DigOut_write(pRTIDioC1DigOut_Port_3_Ch_13);

  /* Terminate for S-Function (rti_commonblock): '<S89>/S-Function1' incorporates:
   *  Constant: '<S47>/Constant'
   */

  /* --- Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_OUT_BL1 --- */
  /* --- [RTI120X, BITOUT] - Port: 1 - Channel: 1 --- */

  /* disable digital output channel 1-1 on port 1 *
   * (set to high-impedance), when the simulation terminates       */
  DioCl1DigOut_setChMaskOutHighZ(pRTIDioC1DigOut_Port_1_Ch_1,
    DIO_CLASS1_HIGH_Z_ON);
  DioCl1DigOut_write(pRTIDioC1DigOut_Port_1_Ch_1);

  /* Terminate for S-Function (rti_commonblock): '<S98>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:3 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X3])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S91>/S-Function1' */

  /* dSPACE RTICAN STD Srvc-Message Block */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (CANTP1_RX_SPMSG_M1_C1_STD)) == DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S92>/S-Function1' */

  /* dSPACE RTICAN STD Srvc-Message Block */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][0] = can_tp1_msg_sleep
            (CANTP1_RX_SPMSG_M1_C2_STD)) == DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S96>/S-Function1' */

  /* dSPACE RTICAN RX Message Block: "RX Message" Id:100 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][1] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C1_RX_STD_0X64])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }

  /* Terminate for S-Function (rti_commonblock): '<S97>/S-Function1' */

  /* dSPACE RTICAN TX Message Block: "TX Message" Id:100 */
  {
    /* ... Set the message into sleep mode */
    while ((rtican_type1_tq_error[0][1] = can_tp1_msg_sleep
            (can_type1_msg_M1[CANTP1_M1_C2_TX_STD_0X64])) ==
           DSMCOM_BUFFER_OVERFLOW) ;
  }
}
