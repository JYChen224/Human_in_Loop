/*
 * Human_in_Loop.h
 *
 * Code generation for model "Human_in_Loop".
 *
 * Model version              : 1.1218
 * Simulink Coder version : 8.13 (R2017b) 24-Jul-2017
 * C source code generated on : Wed May 22 17:03:57 2019
 *
 * Target selection: rti1202.tlc
 * Note: GRT includes extra infrastructure and instrumentation for prototyping
 * Embedded hardware selection: Custom Processor->Custom
 * Code generation objectives: Unspecified
 * Validation result: Not run
 */

#ifndef RTW_HEADER_Human_in_Loop_h_
#define RTW_HEADER_Human_in_Loop_h_
#include <math.h>
#include <string.h>
#ifndef Human_in_Loop_COMMON_INCLUDES_
# define Human_in_Loop_COMMON_INCLUDES_
#include <brtenv.h>
#include <rtkernel.h>
#include <rti_assert.h>
#include <rtidefineddatatypes.h>
#include <dsser.h>
#include <rtican_ds1202.h>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* Human_in_Loop_COMMON_INCLUDES_ */

#include "Human_in_Loop_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rtGetNaN.h"
#include "rt_nonfinite.h"
#include "rt_zcfcn.h"
#include "rtGetInf.h"

/* Macros for accessing real-time model data structure */
#ifndef rtmGetContStateDisabled
# define rtmGetContStateDisabled(rtm)  ((rtm)->contStateDisabled)
#endif

#ifndef rtmSetContStateDisabled
# define rtmSetContStateDisabled(rtm, val) ((rtm)->contStateDisabled = (val))
#endif

#ifndef rtmGetContStates
# define rtmGetContStates(rtm)         ((rtm)->contStates)
#endif

#ifndef rtmSetContStates
# define rtmSetContStates(rtm, val)    ((rtm)->contStates = (val))
#endif

#ifndef rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmGetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm) ((rtm)->CTOutputIncnstWithState)
#endif

#ifndef rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag
# define rtmSetContTimeOutputInconsistentWithStateAtMajorStepFlag(rtm, val) ((rtm)->CTOutputIncnstWithState = (val))
#endif

#ifndef rtmGetDerivCacheNeedsReset
# define rtmGetDerivCacheNeedsReset(rtm) ((rtm)->derivCacheNeedsReset)
#endif

#ifndef rtmSetDerivCacheNeedsReset
# define rtmSetDerivCacheNeedsReset(rtm, val) ((rtm)->derivCacheNeedsReset = (val))
#endif

#ifndef rtmGetIntgData
# define rtmGetIntgData(rtm)           ((rtm)->intgData)
#endif

#ifndef rtmSetIntgData
# define rtmSetIntgData(rtm, val)      ((rtm)->intgData = (val))
#endif

#ifndef rtmGetOdeF
# define rtmGetOdeF(rtm)               ((rtm)->odeF)
#endif

#ifndef rtmSetOdeF
# define rtmSetOdeF(rtm, val)          ((rtm)->odeF = (val))
#endif

#ifndef rtmGetPeriodicContStateIndices
# define rtmGetPeriodicContStateIndices(rtm) ((rtm)->periodicContStateIndices)
#endif

#ifndef rtmSetPeriodicContStateIndices
# define rtmSetPeriodicContStateIndices(rtm, val) ((rtm)->periodicContStateIndices = (val))
#endif

#ifndef rtmGetPeriodicContStateRanges
# define rtmGetPeriodicContStateRanges(rtm) ((rtm)->periodicContStateRanges)
#endif

#ifndef rtmSetPeriodicContStateRanges
# define rtmSetPeriodicContStateRanges(rtm, val) ((rtm)->periodicContStateRanges = (val))
#endif

#ifndef rtmGetZCCacheNeedsReset
# define rtmGetZCCacheNeedsReset(rtm)  ((rtm)->zCCacheNeedsReset)
#endif

#ifndef rtmSetZCCacheNeedsReset
# define rtmSetZCCacheNeedsReset(rtm, val) ((rtm)->zCCacheNeedsReset = (val))
#endif

#ifndef rtmGetdX
# define rtmGetdX(rtm)                 ((rtm)->derivs)
#endif

#ifndef rtmSetdX
# define rtmSetdX(rtm, val)            ((rtm)->derivs = (val))
#endif

#ifndef rtmGetErrorStatus
# define rtmGetErrorStatus(rtm)        ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
# define rtmSetErrorStatus(rtm, val)   ((rtm)->errorStatus = (val))
#endif

#ifndef rtmGetStopRequested
# define rtmGetStopRequested(rtm)      ((rtm)->Timing.stopRequestedFlag)
#endif

#ifndef rtmSetStopRequested
# define rtmSetStopRequested(rtm, val) ((rtm)->Timing.stopRequestedFlag = (val))
#endif

#ifndef rtmGetStopRequestedPtr
# define rtmGetStopRequestedPtr(rtm)   (&((rtm)->Timing.stopRequestedFlag))
#endif

#ifndef rtmGetT
# define rtmGetT(rtm)                  (rtmGetTPtr((rtm))[0])
#endif

/* Block signals for system '<S45>/Mux1' */
typedef struct {
  real_T x[6];                         /* '<S45>/Mux1' */
} B_Mux1_Human_in_Loop_T;

/* Block signals for system '<S60>/MATLAB Function' */
typedef struct {
  real_T y;                            /* '<S60>/MATLAB Function' */
} B_MATLABFunction_Human_in_Loo_T;

/* Block states (auto storage) for system '<S60>/MATLAB Function' */
typedef struct {
  real_T data_mem[500];                /* '<S60>/MATLAB Function' */
} DW_MATLABFunction_Human_in_Lo_T;

/* Block signals for system '<S46>/Mux' */
typedef struct {
  real_T x[2];                         /* '<S46>/Mux' */
} B_Mux_Human_in_Loop_T;

/* Block signals (auto storage) */
typedef struct {
  real_T result_data_m[3000];
  real_T tmp_data_cl[3000];
  real_T tmp_data_m[2600];
  real_T b_result_data[2250];
  real_T tmp_data_k[2250];
  cell_2_Human_in_Loop_T GP_LCB_environment_f1_environme;
  cell_3_Human_in_Loop_T GP_LCB_environment_f2_environme;
  real_T GP_mu_environment_f2[900];
  real_T GP_sigma_environment_f3[900];
  real_T GP_LCB_environment_f1_environ_k[900];
  real_T environment_f1_environment_f2[900];
  real_T environment_f2_environment_f3[900];
  real_T result_data[870];
  real_T tmp_data[841];
  real_T tmp_data_c[841];
  real_T K_data[841];
  real_T B_data[841];
  real_T b_A_data[841];
  real_T b_A_data_b[841];
  real_T tmp_data_p[841];
  real_T K_data_c[841];
  real_T torque_track[750];
  real_T torque_delta_track[750];
  real_T tmp_data_cx[750];
  real_T tmp_data_b[750];
  real_T tmp_data_pb[750];
  real_T z1_data[750];
  real_T z1_data_c[750];
  real_T SFunction1;                   /* '<S103>/S-Function1' */
  real_T Gain;                         /* '<S49>/Gain' */
  real_T u4low2;                       /* '<S49>/0.4low2' */
  real_T x2k1;                         /* '<S102>/Unit Delay1' */
  real_T x1k1;                         /* '<S102>/Unit Delay' */
  real_T Gain1;                        /* '<S102>/Gain1' */
  real_T Gain2;                        /* '<S102>/Gain2' */
  real_T UnitDelay2;                   /* '<S102>/Unit Delay2' */
  real_T Gain4;                        /* '<S102>/Gain4' */
  real_T Add2;                         /* '<S102>/Add2' */
  real_T Gain3;                        /* '<S102>/Gain3' */
  real_T x2k;                          /* '<S102>/Add1' */
  real_T RT4[2];                       /* '<Root>/RT4' */
  real_T SFunction1_o1;                /* '<S81>/S-Function1' */
  real_T SFunction1_o2;                /* '<S81>/S-Function1' */
  real_T Gain_l;                       /* '<S46>/Gain' */
  real_T Gain1_b;                      /* '<S46>/Gain1' */
  real_T RT5[2];                       /* '<Root>/RT5' */
  real_T SFunction1_o1_k;              /* '<S82>/S-Function1' */
  real_T SFunction1_o2_k;              /* '<S82>/S-Function1' */
  real_T Gain2_h;                      /* '<S46>/Gain2' */
  real_T Gain3_m;                      /* '<S46>/Gain3' */
  real_T x2k1_k;                       /* '<S77>/Unit Delay1' */
  real_T x1k1_m;                       /* '<S77>/Unit Delay' */
  real_T Gain1_m;                      /* '<S77>/Gain1' */
  real_T Gain2_c;                      /* '<S77>/Gain2' */
  real_T UnitDelay2_b;                 /* '<S77>/Unit Delay2' */
  real_T Gain4_g;                      /* '<S77>/Gain4' */
  real_T Add2_m;                       /* '<S77>/Add2' */
  real_T Gain3_k;                      /* '<S77>/Gain3' */
  real_T x2k_i;                        /* '<S77>/Add1' */
  real_T RT6[3];                       /* '<Root>/RT6' */
  real_T RT1[4];                       /* '<Root>/RT1' */
  real_T RT2[5];                       /* '<Root>/RT2' */
  real_T Delay1;                       /* '<S4>/Delay1' */
  real_T RT1_c;                        /* '<S23>/RT1' */
  real_T RT2_g;                        /* '<S23>/RT2' */
  real_T RT1_g;                        /* '<S4>/RT1' */
  real_T RT3[4];                       /* '<Root>/RT3' */
  real_T Gain_i;                       /* '<S3>/Gain' */
  real_T Delay3;                       /* '<S17>/Delay3' */
  real_T Delay2;                       /* '<S17>/Delay2' */
  real_T Delay7;                       /* '<S17>/Delay7' */
  real_T Delay6;                       /* '<S17>/Delay6' */
  real_T low_pass;                     /* '<S58>/low_pass' */
  real_T low_pass_f;                   /* '<S59>/low_pass' */
  real_T low_pass_l;                   /* '<S62>/low_pass' */
  real_T low_pass_d;                   /* '<S63>/low_pass' */
  real_T low_pass_k;                   /* '<S69>/low_pass' */
  real_T low_pass_fc;                  /* '<S68>/low_pass' */
  real_T Delay1_p[6];                  /* '<S17>/Delay1' */
  real_T Delay5[6];                    /* '<S17>/Delay5' */
  real_T Delay4[6];                    /* '<S17>/Delay4' */
  real_T Delay8[6];                    /* '<S17>/Delay8' */
  real_T Delay;                        /* '<S4>/Delay' */
  real_T RT2_o;                        /* '<S25>/RT2' */
  real_T RT1_m[2];                     /* '<S25>/RT1' */
  real_T Gain_g;                       /* '<S102>/Gain' */
  real_T x1k;                          /* '<S102>/Add' */
  real_T TSamp;                        /* '<S80>/TSamp' */
  real_T Uk1;                          /* '<S80>/UD' */
  real_T Diff;                         /* '<S80>/Diff' */
  real_T u4low1;                       /* '<S46>/0.4low1' */
  real_T UnitDelay;                    /* '<S76>/Unit Delay' */
  real_T UnitDelay1;                   /* '<S76>/Unit Delay1' */
  real_T Add1;                         /* '<S76>/Add1' */
  real_T Gain_p;                       /* '<S76>/Gain' */
  real_T Gain1_d;                      /* '<S76>/Gain1' */
  real_T Add2_i;                       /* '<S76>/Add2' */
  real_T Add3;                         /* '<S76>/Add3' */
  real_T Gain2_i;                      /* '<S76>/Gain2' */
  real_T Gain_io;                      /* '<S77>/Gain' */
  real_T x1k_i;                        /* '<S77>/Add' */
  real_T SFunction1_o;                 /* '<S50>/S-Function1' */
  real_T SFunction1_k;                 /* '<S51>/S-Function1' */
  real_T SFunction1_j;                 /* '<S52>/S-Function1' */
  real_T SFunction1_d;                 /* '<S53>/S-Function1' */
  real_T SFunction1_l;                 /* '<S54>/S-Function1' */
  real_T SFunction1_i;                 /* '<S55>/S-Function1' */
  real_T Gain_k;                       /* '<S45>/Gain' */
  real_T Gain1_p;                      /* '<S45>/Gain1' */
  real_T Gain2_cz;                     /* '<S45>/Gain2' */
  real_T Gain3_f;                      /* '<S45>/Gain3' */
  real_T Gain4_h;                      /* '<S45>/Gain4' */
  real_T Gain5;                        /* '<S45>/Gain5' */
  real_T high_pass;                    /* '<S64>/high_pass' */
  real_T Abs;                          /* '<S64>/Abs' */
  real_T high_pass_b;                  /* '<S65>/high_pass' */
  real_T Abs_k;                        /* '<S65>/Abs' */
  real_T high_pass_d;                  /* '<S66>/high_pass' */
  real_T Abs_b;                        /* '<S66>/Abs' */
  real_T high_pass_o;                  /* '<S67>/high_pass' */
  real_T Abs_n;                        /* '<S67>/Abs' */
  real_T high_pass_l;                  /* '<S60>/high_pass' */
  real_T Abs_m;                        /* '<S60>/Abs' */
  real_T high_pass_g;                  /* '<S61>/high_pass' */
  real_T Abs_f;                        /* '<S61>/Abs' */
  real_T high_pass_h;                  /* '<S58>/high_pass' */
  real_T Abs_h;                        /* '<S58>/Abs' */
  real_T high_pass_oe;                 /* '<S59>/high_pass' */
  real_T Abs_d;                        /* '<S59>/Abs' */
  real_T high_pass_ba;                 /* '<S62>/high_pass' */
  real_T Abs_j;                        /* '<S62>/Abs' */
  real_T high_pass_hl;                 /* '<S63>/high_pass' */
  real_T Abs_a;                        /* '<S63>/Abs' */
  real_T high_pass_a;                  /* '<S68>/high_pass' */
  real_T Abs_dj;                       /* '<S68>/Abs' */
  real_T high_pass_e;                  /* '<S69>/high_pass' */
  real_T Abs_fy;                       /* '<S69>/Abs' */
  real_T SFunction1_o1_e;              /* '<S98>/S-Function1' */
  real_T SFunction1_o2_b;              /* '<S98>/S-Function1' */
  real_T SFunction1_o3;                /* '<S98>/S-Function1' */
  real_T SFunction1_o4;                /* '<S98>/S-Function1' */
  real_T SFunction1_o5;                /* '<S98>/S-Function1' */
  real_T DataTypeConversion;           /* '<S99>/Data Type Conversion' */
  real_T Gain_g5;                      /* '<S99>/Gain' */
  real_T DataTypeConversion_b;         /* '<S100>/Data Type Conversion' */
  real_T Gain_ik;                      /* '<S100>/Gain' */
  real_T DataTypeConversion_o;         /* '<S101>/Data Type Conversion' */
  real_T Gain_lw;                      /* '<S101>/Gain' */
  real_T SFunction1_o1_i;              /* '<S96>/S-Function1' */
  real_T SFunction1_o2_ks;             /* '<S96>/S-Function1' */
  real_T SFunction1_o3_l;              /* '<S96>/S-Function1' */
  real_T SineWave;                     /* '<S93>/Sine Wave' */
  real_T SFunction1_o1_m;              /* '<S97>/S-Function1' */
  real_T SFunction1_o2_n;              /* '<S97>/S-Function1' */
  real_T SFunction1_o3_b;              /* '<S97>/S-Function1' */
  real_T SFunction1_o4_i;              /* '<S97>/S-Function1' */
  real_T Gain2_j;                      /* '<S12>/Gain2' */
  real_T Gain1_pi;                     /* '<S12>/Gain1' */
  real_T torque_des;                   /* '<S1>/Torque track' */
  real_T torque_delta_des;             /* '<S1>/Torque track' */
  real_T torque_trace[1500];           /* '<S1>/Torque track' */
  real_T torque_delta_trace[1500];     /* '<S1>/Torque track' */
  real_T vel;                          /* '<S12>/MATLAB Function' */
  real_T lrn_cmd;                      /* '<S1>/LRN' */
  real_T lrn_mem[750];                 /* '<S1>/LRN' */
  real_T motor_vel_cmd;                /* '<S1>/Controller' */
  real_T mode;                         /* '<S8>/State Machine' */
  real_T state;                        /* '<S8>/State Machine' */
  real_T stride_time;                  /* '<S8>/State Machine' */
  real_T stride_timer;                 /* '<S8>/State Machine' */
  real_T x[4];                         /* '<S8>/Mux1' */
  real_T torque_dot;                   /* '<S49>/MATLAB Function' */
  real_T torque;                       /* '<S49>/Data process' */
  real_T state_c;                      /* '<S47>/FootSwitch Filter' */
  real_T x_d[3];                       /* '<S46>/Mux2' */
  real_T x_i[3];                       /* '<S46>/Mux1' */
  real_T angle;                        /* '<S46>/Data process1' */
  real_T angle_m;                      /* '<S46>/Data process' */
  real_T x_k[4];                       /* '<S38>/Mux1' */
  real_T parm_torque;                  /* '<S38>/MATLAB Function' */
  real_T x_o[5];                       /* '<S37>/Mux1' */
  real_T parm_acq[2];                  /* '<S34>/Bayesian Function' */
  real_T parm_opt[2];                  /* '<S34>/Bayesian Function' */
  real_T target_esitmate[10000];       /* '<S34>/Bayesian Function' */
  real_T Iteration;                    /* '<S34>/Bayesian Function' */
  real_T trig;                         /* '<S4>/Opt Trigger' */
  real_T stride_num;                   /* '<S4>/Opt Trigger' */
  real_T Trigger;                      /* '<S23>/Timer' */
  real_T Time;                         /* '<S23>/Timer' */
  real_T Count;                        /* '<S23>/Timer' */
  real_T Data1;                        /* '<S30>/MATLAB Function1' */
  real_T Data2;                        /* '<S30>/MATLAB Function1' */
  real_T Data3;                        /* '<S30>/MATLAB Function1' */
  real_T E;                            /* '<S30>/Estimation' */
  real_T Time_p[12];                   /* '<S30>/Estimation' */
  real_T Y[12];                        /* '<S30>/Estimation' */
  real_T Y_E[12];                      /* '<S30>/Estimation' */
  real_T est_curve[24];                /* '<S30>/Estimation' */
  real_T x_this_cycle[5];              /* '<S4>/CMAES' */
  real_T gengeration_count;            /* '<S4>/CMAES' */
  real_T x_this_gengration[50];        /* '<S4>/CMAES' */
  real_T mean_this_generation[5];      /* '<S4>/CMAES' */
  real_T variance_this_generation[5];  /* '<S4>/CMAES' */
  real_T sigma_this_generation;        /* '<S4>/CMAES' */
  real_T loss;                         /* '<S3>/get loss' */
  real_T Cycle_Time;                   /* '<S17>/MATLAB Function' */
  real_T Cycle_Frequency;              /* '<S17>/MATLAB Function' */
  real_T Max[6];                       /* '<S17>/MATLAB Function' */
  real_T Mean[6];                      /* '<S17>/MATLAB Function' */
  real_T RMS[6];                       /* '<S17>/MATLAB Function' */
  real_T iEMG[6];                      /* '<S17>/MATLAB Function' */
  real_T TmpSignalConversionAtSFunctionI[26];/* '<S19>/MATLAB Function1' */
  real_T Mean_m[26];                   /* '<S19>/MATLAB Function1' */
  int32_T SFunction1_o3_ls;            /* '<S31>/S-Function1' */
  uint32_T SFunction1_o2_bb;           /* '<S31>/S-Function1' */
  uint8_T SFunction1_o1_l[16];         /* '<S31>/S-Function1' */
  boolean_T SFunction1_a;              /* '<S88>/S-Function1' */
  B_Mux_Human_in_Loop_T sf_Mux;        /* '<S49>/Mux' */
  B_Mux_Human_in_Loop_T sf_Mux_p;      /* '<S46>/Mux' */
  B_MATLABFunction_Human_in_Loo_T sf_MATLABFunction_l;/* '<S67>/MATLAB Function' */
  B_MATLABFunction_Human_in_Loo_T sf_MATLABFunction_n;/* '<S66>/MATLAB Function' */
  B_MATLABFunction_Human_in_Loo_T sf_MATLABFunction_m;/* '<S65>/MATLAB Function' */
  B_MATLABFunction_Human_in_Loo_T sf_MATLABFunction_g;/* '<S64>/MATLAB Function' */
  B_MATLABFunction_Human_in_Loo_T sf_MATLABFunction_p;/* '<S61>/MATLAB Function' */
  B_MATLABFunction_Human_in_Loo_T sf_MATLABFunction_g5;/* '<S60>/MATLAB Function' */
  B_Mux1_Human_in_Loop_T sf_Mux3;      /* '<S45>/Mux3' */
  B_Mux1_Human_in_Loop_T sf_Mux1_c;    /* '<S45>/Mux1' */
} B_Human_in_Loop_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  emxArray_real_T_1x4_Human_in__T arfiness;/* '<S4>/CMAES' */
  emxArray_real_T_1x4_Human_in__T arx; /* '<S4>/CMAES' */
  emxArray_real_T_1x4_Human_in__T arz; /* '<S4>/CMAES' */
  emxArray_real_T_1_Human_in_Lo_T xmean;/* '<S4>/CMAES' */
  emxArray_real_T_1_Human_in_Lo_T pc;  /* '<S4>/CMAES' */
  emxArray_real_T_1_Human_in_Lo_T ps;  /* '<S4>/CMAES' */
  emxArray_real_T_1x1_Human_in__T B;   /* '<S4>/CMAES' */
  emxArray_real_T_1x1_Human_in__T D;   /* '<S4>/CMAES' */
  emxArray_real_T_1x1_Human_in__T C;   /* '<S4>/CMAES' */
  real_T u4low2_states[3];             /* '<S49>/0.4low2' */
  real_T UnitDelay1_DSTATE;            /* '<S102>/Unit Delay1' */
  real_T UnitDelay_DSTATE;             /* '<S102>/Unit Delay' */
  real_T UnitDelay2_DSTATE;            /* '<S102>/Unit Delay2' */
  real_T UnitDelay1_DSTATE_h;          /* '<S77>/Unit Delay1' */
  real_T UnitDelay_DSTATE_e;           /* '<S77>/Unit Delay' */
  real_T UnitDelay2_DSTATE_a;          /* '<S77>/Unit Delay2' */
  real_T Delay1_DSTATE;                /* '<S4>/Delay1' */
  real_T Delay3_DSTATE;                /* '<S17>/Delay3' */
  real_T Delay2_DSTATE;                /* '<S17>/Delay2' */
  real_T Delay7_DSTATE;                /* '<S17>/Delay7' */
  real_T Delay6_DSTATE;                /* '<S17>/Delay6' */
  real_T Delay1_DSTATE_c[6];           /* '<S17>/Delay1' */
  real_T Delay5_DSTATE[6];             /* '<S17>/Delay5' */
  real_T Delay4_DSTATE[6];             /* '<S17>/Delay4' */
  real_T Delay8_DSTATE[6];             /* '<S17>/Delay8' */
  real_T Delay_DSTATE;                 /* '<S4>/Delay' */
  real_T UD_DSTATE;                    /* '<S80>/UD' */
  real_T u4low1_states[3];             /* '<S46>/0.4low1' */
  real_T UnitDelay_DSTATE_j;           /* '<S76>/Unit Delay' */
  real_T UnitDelay1_DSTATE_a;          /* '<S76>/Unit Delay1' */
  real_T u4low2_tmp;                   /* '<S49>/0.4low2' */
  volatile real_T RT4_Buffer0[2];      /* '<Root>/RT4' */
  volatile real_T RT4_Buffer1[2];      /* '<Root>/RT4' */
  volatile real_T RT5_Buffer0[2];      /* '<Root>/RT5' */
  volatile real_T RT5_Buffer1[2];      /* '<Root>/RT5' */
  volatile real_T RT6_Buffer0[3];      /* '<Root>/RT6' */
  volatile real_T RT6_Buffer1[3];      /* '<Root>/RT6' */
  volatile real_T RT1_Buffer0[4];      /* '<Root>/RT1' */
  volatile real_T RT1_Buffer1[4];      /* '<Root>/RT1' */
  volatile real_T RT2_Buffer0[5];      /* '<Root>/RT2' */
  volatile real_T RT2_Buffer1[5];      /* '<Root>/RT2' */
  volatile real_T RT1_Buffer0_e;       /* '<S23>/RT1' */
  volatile real_T RT1_Buffer1_f;       /* '<S23>/RT1' */
  volatile real_T RT2_Buffer0_i;       /* '<S23>/RT2' */
  volatile real_T RT2_Buffer1_o;       /* '<S23>/RT2' */
  volatile real_T RT1_Buffer0_f;       /* '<S4>/RT1' */
  volatile real_T RT1_Buffer1_o;       /* '<S4>/RT1' */
  volatile real_T RT3_Buffer0[4];      /* '<Root>/RT3' */
  volatile real_T RT3_Buffer1[4];      /* '<Root>/RT3' */
  real_T EMG_Memory[24];               /* '<S17>/Data Store Memory' */
  volatile real_T RT2_Buffer0_b;       /* '<S25>/RT2' */
  volatile real_T RT2_Buffer1_l;       /* '<S25>/RT2' */
  real_T u4low1_tmp;                   /* '<S46>/0.4low1' */
  real_T TorqueMem[4400];              /* '<Root>/Data Store Memory' */
  real_T last_footstate;               /* '<S1>/Torque track' */
  real_T last_footstate_a;             /* '<S1>/LRN' */
  real_T torque_error_memory[1000];    /* '<S1>/LRN' */
  real_T lrn_cmd_memory[1000];         /* '<S1>/LRN' */
  real_T calib_state;                  /* '<S1>/Controller' */
  real_T reg_stride_time;              /* '<S8>/State Machine' */
  real_T reg_stride_time_count;        /* '<S8>/State Machine' */
  real_T reg_mode;                     /* '<S8>/State Machine' */
  real_T reg_state;                    /* '<S8>/State Machine' */
  real_T bt_run;                       /* '<S8>/State Machine' */
  real_T reg_last_switch;              /* '<S8>/State Machine' */
  real_T data[15];                     /* '<S49>/MATLAB Function' */
  real_T torque_zero;                  /* '<S49>/Data process' */
  real_T foot_state;                   /* '<S47>/FootSwitch Filter' */
  real_T filter_time;                  /* '<S47>/FootSwitch Filter' */
  real_T angle_zero;                   /* '<S46>/Data process1' */
  real_T angle_zero_f;                 /* '<S46>/Data process' */
  real_T iter;                         /* '<S34>/Bayesian Function' */
  real_T X[60];                        /* '<S34>/Bayesian Function' */
  real_T Y[30];                        /* '<S34>/Bayesian Function' */
  real_T K[900];                       /* '<S34>/Bayesian Function' */
  real_T Mu[10000];                    /* '<S34>/Bayesian Function' */
  real_T LCB[10000];                   /* '<S34>/Bayesian Function' */
  real_T x_opt[2];                     /* '<S34>/Bayesian Function' */
  real_T last_state;                   /* '<S4>/Opt Trigger' */
  real_T reg_stride_num;               /* '<S4>/Opt Trigger' */
  real_T count;                        /* '<S23>/Timer' */
  real_T time;                         /* '<S23>/Timer' */
  real_T last_count;                   /* '<S30>/Estimation' */
  real_T x[20];                        /* '<S30>/Estimation' */
  real_T y[20];                        /* '<S30>/Estimation' */
  real_T theta[2];                     /* '<S30>/Estimation' */
  real_T num;                          /* '<S30>/Estimation' */
  real_T last_trigger;                 /* '<S4>/CMAES' */
  real_T count_val;                    /* '<S4>/CMAES' */
  real_T cycle_val;                    /* '<S4>/CMAES' */
  real_T sigma;                        /* '<S4>/CMAES' */
  real_T Reset;                        /* '<S4>/CMAES' */
  real_T SingleCycleData[2600];        /* '<S19>/Data Store Memory' */
  uint32_T state[625];                 /* '<S34>/Bayesian Function' */
  uint32_T state_l[625];               /* '<S4>/CMAES' */
  int_T SFunction1_IWORK[2];           /* '<S91>/S-Function1' */
  int_T SFunction1_IWORK_i[2];         /* '<S92>/S-Function1' */
  volatile int8_T RT4_write_buf;       /* '<Root>/RT4' */
  volatile int8_T RT4_read_buf;        /* '<Root>/RT4' */
  volatile int8_T RT4_last_buf_wr;     /* '<Root>/RT4' */
  volatile int8_T RT5_write_buf;       /* '<Root>/RT5' */
  volatile int8_T RT5_read_buf;        /* '<Root>/RT5' */
  volatile int8_T RT5_last_buf_wr;     /* '<Root>/RT5' */
  volatile int8_T RT6_write_buf;       /* '<Root>/RT6' */
  volatile int8_T RT6_read_buf;        /* '<Root>/RT6' */
  volatile int8_T RT6_last_buf_wr;     /* '<Root>/RT6' */
  volatile int8_T RT1_write_buf;       /* '<Root>/RT1' */
  volatile int8_T RT1_read_buf;        /* '<Root>/RT1' */
  volatile int8_T RT1_last_buf_wr;     /* '<Root>/RT1' */
  volatile int8_T RT2_write_buf;       /* '<Root>/RT2' */
  volatile int8_T RT2_read_buf;        /* '<Root>/RT2' */
  volatile int8_T RT2_last_buf_wr;     /* '<Root>/RT2' */
  volatile int8_T RT1_write_buf_f;     /* '<S23>/RT1' */
  volatile int8_T RT1_read_buf_p;      /* '<S23>/RT1' */
  volatile int8_T RT1_last_buf_wr_j;   /* '<S23>/RT1' */
  volatile int8_T RT2_write_buf_g;     /* '<S23>/RT2' */
  volatile int8_T RT2_read_buf_a;      /* '<S23>/RT2' */
  volatile int8_T RT2_last_buf_wr_a;   /* '<S23>/RT2' */
  volatile int8_T RT1_write_buf_n;     /* '<S4>/RT1' */
  volatile int8_T RT1_read_buf_i;      /* '<S4>/RT1' */
  volatile int8_T RT1_last_buf_wr_c;   /* '<S4>/RT1' */
  volatile int8_T RT3_write_buf;       /* '<Root>/RT3' */
  volatile int8_T RT3_read_buf;        /* '<Root>/RT3' */
  volatile int8_T RT3_last_buf_wr;     /* '<Root>/RT3' */
  volatile int8_T RT2_write_buf_f;     /* '<S25>/RT2' */
  volatile int8_T RT2_read_buf_aw;     /* '<S25>/RT2' */
  volatile int8_T RT2_last_buf_wr_b;   /* '<S25>/RT2' */
  boolean_T iter_not_empty;            /* '<S34>/Bayesian Function' */
  boolean_T last_trigger_not_empty;    /* '<S4>/CMAES' */
  DW_MATLABFunction_Human_in_Lo_T sf_MATLABFunction_l;/* '<S67>/MATLAB Function' */
  DW_MATLABFunction_Human_in_Lo_T sf_MATLABFunction_n;/* '<S66>/MATLAB Function' */
  DW_MATLABFunction_Human_in_Lo_T sf_MATLABFunction_m;/* '<S65>/MATLAB Function' */
  DW_MATLABFunction_Human_in_Lo_T sf_MATLABFunction_g;/* '<S64>/MATLAB Function' */
  DW_MATLABFunction_Human_in_Lo_T sf_MATLABFunction_p;/* '<S61>/MATLAB Function' */
  DW_MATLABFunction_Human_in_Lo_T sf_MATLABFunction_g5;/* '<S60>/MATLAB Function' */
} DW_Human_in_Loop_T;

/* Continuous states (auto storage) */
typedef struct {
  real_T low_pass_CSTATE[2];           /* '<S58>/low_pass' */
  real_T low_pass_CSTATE_o[2];         /* '<S59>/low_pass' */
  real_T low_pass_CSTATE_n[2];         /* '<S62>/low_pass' */
  real_T low_pass_CSTATE_p[2];         /* '<S63>/low_pass' */
  real_T low_pass_CSTATE_on[2];        /* '<S69>/low_pass' */
  real_T low_pass_CSTATE_f[2];         /* '<S68>/low_pass' */
  real_T high_pass_CSTATE[4];          /* '<S64>/high_pass' */
  real_T high_pass_CSTATE_b[4];        /* '<S65>/high_pass' */
  real_T high_pass_CSTATE_c[4];        /* '<S66>/high_pass' */
  real_T high_pass_CSTATE_m[4];        /* '<S67>/high_pass' */
  real_T high_pass_CSTATE_a[4];        /* '<S60>/high_pass' */
  real_T high_pass_CSTATE_my[4];       /* '<S61>/high_pass' */
  real_T high_pass_CSTATE_p[2];        /* '<S58>/high_pass' */
  real_T high_pass_CSTATE_n[2];        /* '<S59>/high_pass' */
  real_T high_pass_CSTATE_pw[2];       /* '<S62>/high_pass' */
  real_T high_pass_CSTATE_k[2];        /* '<S63>/high_pass' */
  real_T high_pass_CSTATE_l[2];        /* '<S68>/high_pass' */
  real_T high_pass_CSTATE_h[2];        /* '<S69>/high_pass' */
} X_Human_in_Loop_T;

/* State derivatives (auto storage) */
typedef struct {
  real_T low_pass_CSTATE[2];           /* '<S58>/low_pass' */
  real_T low_pass_CSTATE_o[2];         /* '<S59>/low_pass' */
  real_T low_pass_CSTATE_n[2];         /* '<S62>/low_pass' */
  real_T low_pass_CSTATE_p[2];         /* '<S63>/low_pass' */
  real_T low_pass_CSTATE_on[2];        /* '<S69>/low_pass' */
  real_T low_pass_CSTATE_f[2];         /* '<S68>/low_pass' */
  real_T high_pass_CSTATE[4];          /* '<S64>/high_pass' */
  real_T high_pass_CSTATE_b[4];        /* '<S65>/high_pass' */
  real_T high_pass_CSTATE_c[4];        /* '<S66>/high_pass' */
  real_T high_pass_CSTATE_m[4];        /* '<S67>/high_pass' */
  real_T high_pass_CSTATE_a[4];        /* '<S60>/high_pass' */
  real_T high_pass_CSTATE_my[4];       /* '<S61>/high_pass' */
  real_T high_pass_CSTATE_p[2];        /* '<S58>/high_pass' */
  real_T high_pass_CSTATE_n[2];        /* '<S59>/high_pass' */
  real_T high_pass_CSTATE_pw[2];       /* '<S62>/high_pass' */
  real_T high_pass_CSTATE_k[2];        /* '<S63>/high_pass' */
  real_T high_pass_CSTATE_l[2];        /* '<S68>/high_pass' */
  real_T high_pass_CSTATE_h[2];        /* '<S69>/high_pass' */
} XDot_Human_in_Loop_T;

/* State disabled  */
typedef struct {
  boolean_T low_pass_CSTATE[2];        /* '<S58>/low_pass' */
  boolean_T low_pass_CSTATE_o[2];      /* '<S59>/low_pass' */
  boolean_T low_pass_CSTATE_n[2];      /* '<S62>/low_pass' */
  boolean_T low_pass_CSTATE_p[2];      /* '<S63>/low_pass' */
  boolean_T low_pass_CSTATE_on[2];     /* '<S69>/low_pass' */
  boolean_T low_pass_CSTATE_f[2];      /* '<S68>/low_pass' */
  boolean_T high_pass_CSTATE[4];       /* '<S64>/high_pass' */
  boolean_T high_pass_CSTATE_b[4];     /* '<S65>/high_pass' */
  boolean_T high_pass_CSTATE_c[4];     /* '<S66>/high_pass' */
  boolean_T high_pass_CSTATE_m[4];     /* '<S67>/high_pass' */
  boolean_T high_pass_CSTATE_a[4];     /* '<S60>/high_pass' */
  boolean_T high_pass_CSTATE_my[4];    /* '<S61>/high_pass' */
  boolean_T high_pass_CSTATE_p[2];     /* '<S58>/high_pass' */
  boolean_T high_pass_CSTATE_n[2];     /* '<S59>/high_pass' */
  boolean_T high_pass_CSTATE_pw[2];    /* '<S62>/high_pass' */
  boolean_T high_pass_CSTATE_k[2];     /* '<S63>/high_pass' */
  boolean_T high_pass_CSTATE_l[2];     /* '<S68>/high_pass' */
  boolean_T high_pass_CSTATE_h[2];     /* '<S69>/high_pass' */
} XDis_Human_in_Loop_T;

/* Zero-crossing (trigger) state */
typedef struct {
  ZCSigState SoftwareInterrupt_Trig_ZCE;/* '<S25>/Software Interrupt' */
  ZCSigState MeanCalculate_Trig_ZCE;   /* '<S16>/Mean Calculate' */
} PrevZCX_Human_in_Loop_T;

#ifndef ODE1_INTG
#define ODE1_INTG

/* ODE1 Integration Data */
typedef struct {
  real_T *f[1];                        /* derivatives */
} ODE1_IntgData;

#endif

/* Parameters (auto storage) */
struct P_Human_in_Loop_T_ {
  real_T Controller_CALIB_SPEED;       /* Mask Parameter: Controller_CALIB_SPEED
                                        * Referenced by: '<S1>/Controller'
                                        */
  real_T Controller_CALIB_TORQUE;      /* Mask Parameter: Controller_CALIB_TORQUE
                                        * Referenced by: '<S1>/Controller'
                                        */
  real_T Controller_FOLLOW_SLACK_ANGLE;/* Mask Parameter: Controller_FOLLOW_SLACK_ANGLE
                                        * Referenced by: '<S1>/Controller'
                                        */
  real_T DiscreteDerivative_ICPrevScaled;/* Mask Parameter: DiscreteDerivative_ICPrevScaled
                                          * Referenced by: '<S80>/UD'
                                          */
  real_T MATLABFunction_MAX_MOTOR_ANGLE;/* Mask Parameter: MATLABFunction_MAX_MOTOR_ANGLE
                                         * Referenced by: '<S12>/MATLAB Function'
                                         */
  real_T MATLABFunction_MAX_SPEED;     /* Mask Parameter: MATLABFunction_MAX_SPEED
                                        * Referenced by: '<S12>/MATLAB Function'
                                        */
  real_T MATLABFunction_MAX_TORQUE;    /* Mask Parameter: MATLABFunction_MAX_TORQUE
                                        * Referenced by: '<S12>/MATLAB Function'
                                        */
  real_T MATLABFunction_MIN_MOTOR_ANGLE;/* Mask Parameter: MATLABFunction_MIN_MOTOR_ANGLE
                                         * Referenced by: '<S12>/MATLAB Function'
                                         */
  real_T MultiCycleAnalysis1_N;        /* Mask Parameter: MultiCycleAnalysis1_N
                                        * Referenced by: '<S19>/Constant'
                                        */
  real_T Controller_SLACK_SPEED;       /* Mask Parameter: Controller_SLACK_SPEED
                                        * Referenced by: '<S1>/Controller'
                                        */
  real_T Controller_SWING_MAX_ANGLE;   /* Mask Parameter: Controller_SWING_MAX_ANGLE
                                        * Referenced by: '<S1>/Controller'
                                        */
  real_T uOrderTD_T;                   /* Mask Parameter: uOrderTD_T
                                        * Referenced by:
                                        *   '<S76>/Gain1'
                                        *   '<S76>/Gain2'
                                        */
  real_T uOrderTD_T1;                  /* Mask Parameter: uOrderTD_T1
                                        * Referenced by:
                                        *   '<S102>/Gain1'
                                        *   '<S102>/Gain2'
                                        *   '<S102>/Gain4'
                                        */
  real_T uOrderTD_T1_l;                /* Mask Parameter: uOrderTD_T1_l
                                        * Referenced by:
                                        *   '<S77>/Gain1'
                                        *   '<S77>/Gain2'
                                        *   '<S77>/Gain4'
                                        */
  real_T uOrderTD_T2;                  /* Mask Parameter: uOrderTD_T2
                                        * Referenced by:
                                        *   '<S102>/Gain1'
                                        *   '<S102>/Gain2'
                                        *   '<S102>/Gain4'
                                        */
  real_T uOrderTD_T2_h;                /* Mask Parameter: uOrderTD_T2_h
                                        * Referenced by:
                                        *   '<S77>/Gain1'
                                        *   '<S77>/Gain2'
                                        *   '<S77>/Gain4'
                                        */
  real_T BayesianFunction_Theta;       /* Mask Parameter: BayesianFunction_Theta
                                        * Referenced by: '<S34>/Bayesian Function'
                                        */
  real_T uOrderTD_Ts;                  /* Mask Parameter: uOrderTD_Ts
                                        * Referenced by:
                                        *   '<S102>/Gain'
                                        *   '<S102>/Gain3'
                                        */
  real_T uOrderTD_Ts_n;                /* Mask Parameter: uOrderTD_Ts_n
                                        * Referenced by:
                                        *   '<S77>/Gain'
                                        *   '<S77>/Gain3'
                                        */
  real_T SingleCycleAnalysis1_Ts;      /* Mask Parameter: SingleCycleAnalysis1_Ts
                                        * Referenced by: '<S17>/Sample Time'
                                        */
  real_T uOrderTD_Ts_g;                /* Mask Parameter: uOrderTD_Ts_g
                                        * Referenced by: '<S76>/Gain'
                                        */
  real_T LRN_error_filter_k;           /* Mask Parameter: LRN_error_filter_k
                                        * Referenced by: '<S1>/LRN'
                                        */
  real_T Dataprocess_load_vol_gain;    /* Mask Parameter: Dataprocess_load_vol_gain
                                        * Referenced by: '<S49>/Data process'
                                        */
  real_T Dataprocess_load_vol_offset;  /* Mask Parameter: Dataprocess_load_vol_offset
                                        * Referenced by: '<S49>/Data process'
                                        */
  real_T LRN_lrn_shrink;               /* Mask Parameter: LRN_lrn_shrink
                                        * Referenced by: '<S1>/LRN'
                                        */
  int32_T LRN_time_delay;              /* Mask Parameter: LRN_time_delay
                                        * Referenced by: '<S1>/LRN'
                                        */
  int8_T OptTrigger_STRIDE_NUM;        /* Mask Parameter: OptTrigger_STRIDE_NUM
                                        * Referenced by: '<S4>/Opt Trigger'
                                        */
  int8_T Controller_run_mode;          /* Mask Parameter: Controller_run_mode
                                        * Referenced by: '<S1>/Controller'
                                        */
  boolean_T StateMachine_BT_CALIB;     /* Mask Parameter: StateMachine_BT_CALIB
                                        * Referenced by: '<S8>/State Machine'
                                        */
  boolean_T StateMachine_BT_ERROR;     /* Mask Parameter: StateMachine_BT_ERROR
                                        * Referenced by: '<S8>/State Machine'
                                        */
  boolean_T StateMachine_BT_IDLE;      /* Mask Parameter: StateMachine_BT_IDLE
                                        * Referenced by: '<S8>/State Machine'
                                        */
  boolean_T LRN_BT_LRN_CLEAR;          /* Mask Parameter: LRN_BT_LRN_CLEAR
                                        * Referenced by: '<S1>/LRN'
                                        */
  boolean_T LRN_BT_LRN_ON;             /* Mask Parameter: LRN_BT_LRN_ON
                                        * Referenced by: '<S1>/LRN'
                                        */
  boolean_T Timer_BT_METAB_EST;        /* Mask Parameter: Timer_BT_METAB_EST
                                        * Referenced by: '<S23>/Timer'
                                        */
  boolean_T OptTrigger_BT_OPT;         /* Mask Parameter: OptTrigger_BT_OPT
                                        * Referenced by: '<S4>/Opt Trigger'
                                        */
  boolean_T BayesianFunction_BT_OPT_RESET;/* Mask Parameter: BayesianFunction_BT_OPT_RESET
                                           * Referenced by: '<S34>/Bayesian Function'
                                           */
  boolean_T Dataprocess_BT_RESET_ANKLE;/* Mask Parameter: Dataprocess_BT_RESET_ANKLE
                                        * Referenced by: '<S46>/Data process'
                                        */
  boolean_T Dataprocess1_BT_RESET_MOTOR;/* Mask Parameter: Dataprocess1_BT_RESET_MOTOR
                                         * Referenced by: '<S46>/Data process1'
                                         */
  boolean_T Dataprocess_BT_RESET_TORQUE;/* Mask Parameter: Dataprocess_BT_RESET_TORQUE
                                         * Referenced by: '<S49>/Data process'
                                         */
  boolean_T StateMachine_BT_RUN;       /* Mask Parameter: StateMachine_BT_RUN
                                        * Referenced by: '<S8>/State Machine'
                                        */
  boolean_T StateMachine_BT_SLACK;     /* Mask Parameter: StateMachine_BT_SLACK
                                        * Referenced by: '<S8>/State Machine'
                                        */
  real_T Mean_Y0;                      /* Computed Parameter: Mean_Y0
                                        * Referenced by: '<S19>/Mean'
                                        */
  real_T DataStoreMemory_InitialValue[2600];/* Expression: zeros(26,100)
                                             * Referenced by: '<S19>/Data Store Memory'
                                             */
  real_T E_Y0;                         /* Computed Parameter: E_Y0
                                        * Referenced by: '<S30>/E'
                                        */
  real_T parm_Y0;                      /* Computed Parameter: parm_Y0
                                        * Referenced by: '<S34>/parm'
                                        */
  real_T Gain2_Gain;                   /* Expression: 1/360
                                        * Referenced by: '<S12>/Gain2'
                                        */
  real_T Gain1_Gain;                   /* Expression: 1/20
                                        * Referenced by: '<S12>/Gain1'
                                        */
  real_T Gain_Gain;                    /* Expression: 10
                                        * Referenced by: '<S49>/Gain'
                                        */
  real_T u4low2_NumCoef[4];            /* Expression: [0.219606211225382   0.658818633676145   0.658818633676145   0.219606211225382]*1e-3
                                        * Referenced by: '<S49>/0.4low2'
                                        */
  real_T u4low2_DenCoef[4];            /* Expression: [1.000000000000000  -2.748835809214676   2.528231219142559  -0.777638560238080]
                                        * Referenced by: '<S49>/0.4low2'
                                        */
  real_T u4low2_InitialStates;         /* Expression: 0
                                        * Referenced by: '<S49>/0.4low2'
                                        */
  real_T UnitDelay1_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S102>/Unit Delay1'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S102>/Unit Delay'
                                        */
  real_T UnitDelay2_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S102>/Unit Delay2'
                                        */
  real_T RT4_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<Root>/RT4'
                                        */
  real_T Gain_Gain_c;                  /* Expression: 360/1250
                                        * Referenced by: '<S46>/Gain'
                                        */
  real_T Gain1_Gain_e;                 /* Expression: 360/1250
                                        * Referenced by: '<S46>/Gain1'
                                        */
  real_T RT5_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<Root>/RT5'
                                        */
  real_T Gain2_Gain_e;                 /* Expression: -360.0/1250.0
                                        * Referenced by: '<S46>/Gain2'
                                        */
  real_T Gain3_Gain;                   /* Expression: -360.0/1250.0
                                        * Referenced by: '<S46>/Gain3'
                                        */
  real_T UnitDelay1_InitialCondition_j;/* Expression: 0
                                        * Referenced by: '<S77>/Unit Delay1'
                                        */
  real_T UnitDelay_InitialCondition_j; /* Expression: 0
                                        * Referenced by: '<S77>/Unit Delay'
                                        */
  real_T UnitDelay2_InitialCondition_h;/* Expression: 0
                                        * Referenced by: '<S77>/Unit Delay2'
                                        */
  real_T RT6_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<Root>/RT6'
                                        */
  real_T RT1_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<Root>/RT1'
                                        */
  real_T Kp_Value;                     /* Expression: 0
                                        * Referenced by: '<S37>/Kp'
                                        */
  real_T Kd_Value;                     /* Expression: 0
                                        * Referenced by: '<S37>/Kd'
                                        */
  real_T Kl_Value;                     /* Expression: 0
                                        * Referenced by: '<S37>/Kl'
                                        */
  real_T Ks_Value;                     /* Expression: 2
                                        * Referenced by: '<S37>/Ks'
                                        */
  real_T Ko_Value;                     /* Expression: 0
                                        * Referenced by: '<S37>/Ko'
                                        */
  real_T RT2_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<Root>/RT2'
                                        */
  real_T Reset_Value;                  /* Expression: 0
                                        * Referenced by: '<S4>/Reset'
                                        */
  real_T Delay1_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S4>/Delay1'
                                        */
  real_T RT1_InitialCondition_c;       /* Expression: 0
                                        * Referenced by: '<S23>/RT1'
                                        */
  real_T RT2_InitialCondition_i;       /* Expression: 0
                                        * Referenced by: '<S23>/RT2'
                                        */
  real_T RT1_InitialCondition_cs;      /* Expression: 0
                                        * Referenced by: '<S4>/RT1'
                                        */
  real_T sigma_init_Value;             /* Expression: 0.5
                                        * Referenced by: '<S4>/sigma_init'
                                        */
  real_T rise_time_Value;              /* Expression: 25
                                        * Referenced by: '<S38>/rise_time'
                                        */
  real_T peak_time_Value;              /* Expression: 49
                                        * Referenced by: '<S38>/peak_time'
                                        */
  real_T fall_time_Value;              /* Expression: 15
                                        * Referenced by: '<S38>/fall_time'
                                        */
  real_T RT3_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<Root>/RT3'
                                        */
  real_T Gain_Gain_j;                  /* Expression: 5000
                                        * Referenced by: '<S3>/Gain'
                                        */
  real_T Delay3_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S17>/Delay3'
                                        */
  real_T Delay2_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S17>/Delay2'
                                        */
  real_T Delay7_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S17>/Delay7'
                                        */
  real_T Delay6_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S17>/Delay6'
                                        */
  real_T low_pass_A[3];                /* Computed Parameter: low_pass_A
                                        * Referenced by: '<S58>/low_pass'
                                        */
  real_T low_pass_B;                   /* Computed Parameter: low_pass_B
                                        * Referenced by: '<S58>/low_pass'
                                        */
  real_T low_pass_C;                   /* Computed Parameter: low_pass_C
                                        * Referenced by: '<S58>/low_pass'
                                        */
  real_T low_pass_InitialCondition;    /* Expression: 0
                                        * Referenced by: '<S58>/low_pass'
                                        */
  real_T low_pass_A_k[3];              /* Computed Parameter: low_pass_A_k
                                        * Referenced by: '<S59>/low_pass'
                                        */
  real_T low_pass_B_o;                 /* Computed Parameter: low_pass_B_o
                                        * Referenced by: '<S59>/low_pass'
                                        */
  real_T low_pass_C_g;                 /* Computed Parameter: low_pass_C_g
                                        * Referenced by: '<S59>/low_pass'
                                        */
  real_T low_pass_InitialCondition_e;  /* Expression: 0
                                        * Referenced by: '<S59>/low_pass'
                                        */
  real_T low_pass_A_f[3];              /* Computed Parameter: low_pass_A_f
                                        * Referenced by: '<S62>/low_pass'
                                        */
  real_T low_pass_B_j;                 /* Computed Parameter: low_pass_B_j
                                        * Referenced by: '<S62>/low_pass'
                                        */
  real_T low_pass_C_d;                 /* Computed Parameter: low_pass_C_d
                                        * Referenced by: '<S62>/low_pass'
                                        */
  real_T low_pass_InitialCondition_ew; /* Expression: 0
                                        * Referenced by: '<S62>/low_pass'
                                        */
  real_T low_pass_A_a[3];              /* Computed Parameter: low_pass_A_a
                                        * Referenced by: '<S63>/low_pass'
                                        */
  real_T low_pass_B_m;                 /* Computed Parameter: low_pass_B_m
                                        * Referenced by: '<S63>/low_pass'
                                        */
  real_T low_pass_C_i;                 /* Computed Parameter: low_pass_C_i
                                        * Referenced by: '<S63>/low_pass'
                                        */
  real_T low_pass_InitialCondition_c;  /* Expression: 0
                                        * Referenced by: '<S63>/low_pass'
                                        */
  real_T low_pass_A_b[3];              /* Computed Parameter: low_pass_A_b
                                        * Referenced by: '<S69>/low_pass'
                                        */
  real_T low_pass_B_g;                 /* Computed Parameter: low_pass_B_g
                                        * Referenced by: '<S69>/low_pass'
                                        */
  real_T low_pass_C_a;                 /* Computed Parameter: low_pass_C_a
                                        * Referenced by: '<S69>/low_pass'
                                        */
  real_T low_pass_InitialCondition_l;  /* Expression: 0
                                        * Referenced by: '<S69>/low_pass'
                                        */
  real_T low_pass_A_i[3];              /* Computed Parameter: low_pass_A_i
                                        * Referenced by: '<S68>/low_pass'
                                        */
  real_T low_pass_B_os;                /* Computed Parameter: low_pass_B_os
                                        * Referenced by: '<S68>/low_pass'
                                        */
  real_T low_pass_C_c;                 /* Computed Parameter: low_pass_C_c
                                        * Referenced by: '<S68>/low_pass'
                                        */
  real_T low_pass_InitialCondition_i;  /* Expression: 0
                                        * Referenced by: '<S68>/low_pass'
                                        */
  real_T Delay1_InitialCondition_o;    /* Expression: 0
                                        * Referenced by: '<S17>/Delay1'
                                        */
  real_T Delay5_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S17>/Delay5'
                                        */
  real_T Delay4_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S17>/Delay4'
                                        */
  real_T Delay8_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S17>/Delay8'
                                        */
  real_T DataStoreMemory_InitialValue_o[24];/* Expression: zeros(6,4)
                                             * Referenced by: '<S17>/Data Store Memory'
                                             */
  real_T Delay_InitialCondition;       /* Expression: 0
                                        * Referenced by: '<S4>/Delay'
                                        */
  real_T RT2_InitialCondition_p;       /* Expression: 0
                                        * Referenced by: '<S25>/RT2'
                                        */
  real_T TSamp_WtEt;                   /* Computed Parameter: TSamp_WtEt
                                        * Referenced by: '<S80>/TSamp'
                                        */
  real_T u4low1_NumCoef[4];            /* Expression: [   0.001567010350588   0.004701031051765   0.004701031051765   0.001567010350588
                                          ]
                                        * Referenced by: '<S46>/0.4low1'
                                        */
  real_T u4low1_DenCoef[4];            /* Expression: [   1.000000000000000  -2.498608344691178   2.115254127003159  -0.604109699507275
                                          ]
                                        * Referenced by: '<S46>/0.4low1'
                                        */
  real_T u4low1_InitialStates;         /* Expression: 0
                                        * Referenced by: '<S46>/0.4low1'
                                        */
  real_T UnitDelay_InitialCondition_g; /* Expression: 0
                                        * Referenced by: '<S76>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition_i;/* Expression: 0
                                        * Referenced by: '<S76>/Unit Delay1'
                                        */
  real_T Gain_Gain_h;                  /* Expression: 10
                                        * Referenced by: '<S45>/Gain'
                                        */
  real_T Gain1_Gain_b;                 /* Expression: 10
                                        * Referenced by: '<S45>/Gain1'
                                        */
  real_T Gain2_Gain_k;                 /* Expression: 10
                                        * Referenced by: '<S45>/Gain2'
                                        */
  real_T Gain3_Gain_b;                 /* Expression: 10
                                        * Referenced by: '<S45>/Gain3'
                                        */
  real_T Gain4_Gain;                   /* Expression: 10
                                        * Referenced by: '<S45>/Gain4'
                                        */
  real_T Gain5_Gain;                   /* Expression: 10
                                        * Referenced by: '<S45>/Gain5'
                                        */
  real_T high_pass_A[7];               /* Computed Parameter: high_pass_A
                                        * Referenced by: '<S64>/high_pass'
                                        */
  real_T high_pass_B;                  /* Computed Parameter: high_pass_B
                                        * Referenced by: '<S64>/high_pass'
                                        */
  real_T high_pass_C;                  /* Computed Parameter: high_pass_C
                                        * Referenced by: '<S64>/high_pass'
                                        */
  real_T high_pass_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S64>/high_pass'
                                        */
  real_T high_pass_A_p[7];             /* Computed Parameter: high_pass_A_p
                                        * Referenced by: '<S65>/high_pass'
                                        */
  real_T high_pass_B_g;                /* Computed Parameter: high_pass_B_g
                                        * Referenced by: '<S65>/high_pass'
                                        */
  real_T high_pass_C_i;                /* Computed Parameter: high_pass_C_i
                                        * Referenced by: '<S65>/high_pass'
                                        */
  real_T high_pass_InitialCondition_g; /* Expression: 0
                                        * Referenced by: '<S65>/high_pass'
                                        */
  real_T high_pass_A_e[7];             /* Computed Parameter: high_pass_A_e
                                        * Referenced by: '<S66>/high_pass'
                                        */
  real_T high_pass_B_i;                /* Computed Parameter: high_pass_B_i
                                        * Referenced by: '<S66>/high_pass'
                                        */
  real_T high_pass_C_m;                /* Computed Parameter: high_pass_C_m
                                        * Referenced by: '<S66>/high_pass'
                                        */
  real_T high_pass_InitialCondition_n; /* Expression: 0
                                        * Referenced by: '<S66>/high_pass'
                                        */
  real_T high_pass_A_o[7];             /* Computed Parameter: high_pass_A_o
                                        * Referenced by: '<S67>/high_pass'
                                        */
  real_T high_pass_B_c;                /* Computed Parameter: high_pass_B_c
                                        * Referenced by: '<S67>/high_pass'
                                        */
  real_T high_pass_C_j;                /* Computed Parameter: high_pass_C_j
                                        * Referenced by: '<S67>/high_pass'
                                        */
  real_T high_pass_InitialCondition_b; /* Expression: 0
                                        * Referenced by: '<S67>/high_pass'
                                        */
  real_T high_pass_A_h[7];             /* Computed Parameter: high_pass_A_h
                                        * Referenced by: '<S60>/high_pass'
                                        */
  real_T high_pass_B_d;                /* Computed Parameter: high_pass_B_d
                                        * Referenced by: '<S60>/high_pass'
                                        */
  real_T high_pass_C_mg;               /* Computed Parameter: high_pass_C_mg
                                        * Referenced by: '<S60>/high_pass'
                                        */
  real_T high_pass_InitialCondition_bd;/* Expression: 0
                                        * Referenced by: '<S60>/high_pass'
                                        */
  real_T high_pass_A_i[7];             /* Computed Parameter: high_pass_A_i
                                        * Referenced by: '<S61>/high_pass'
                                        */
  real_T high_pass_B_h;                /* Computed Parameter: high_pass_B_h
                                        * Referenced by: '<S61>/high_pass'
                                        */
  real_T high_pass_C_a;                /* Computed Parameter: high_pass_C_a
                                        * Referenced by: '<S61>/high_pass'
                                        */
  real_T high_pass_InitialCondition_m; /* Expression: 0
                                        * Referenced by: '<S61>/high_pass'
                                        */
  real_T high_pass_A_d[3];             /* Computed Parameter: high_pass_A_d
                                        * Referenced by: '<S58>/high_pass'
                                        */
  real_T high_pass_B_b;                /* Computed Parameter: high_pass_B_b
                                        * Referenced by: '<S58>/high_pass'
                                        */
  real_T high_pass_C_ja[2];            /* Computed Parameter: high_pass_C_ja
                                        * Referenced by: '<S58>/high_pass'
                                        */
  real_T high_pass_D;                  /* Computed Parameter: high_pass_D
                                        * Referenced by: '<S58>/high_pass'
                                        */
  real_T high_pass_InitialCondition_p; /* Expression: 0
                                        * Referenced by: '<S58>/high_pass'
                                        */
  real_T high_pass_A_b[3];             /* Computed Parameter: high_pass_A_b
                                        * Referenced by: '<S59>/high_pass'
                                        */
  real_T high_pass_B_j;                /* Computed Parameter: high_pass_B_j
                                        * Referenced by: '<S59>/high_pass'
                                        */
  real_T high_pass_C_d[2];             /* Computed Parameter: high_pass_C_d
                                        * Referenced by: '<S59>/high_pass'
                                        */
  real_T high_pass_D_m;                /* Computed Parameter: high_pass_D_m
                                        * Referenced by: '<S59>/high_pass'
                                        */
  real_T high_pass_InitialCondition_n4;/* Expression: 0
                                        * Referenced by: '<S59>/high_pass'
                                        */
  real_T high_pass_A_j[3];             /* Computed Parameter: high_pass_A_j
                                        * Referenced by: '<S62>/high_pass'
                                        */
  real_T high_pass_B_e;                /* Computed Parameter: high_pass_B_e
                                        * Referenced by: '<S62>/high_pass'
                                        */
  real_T high_pass_C_g[2];             /* Computed Parameter: high_pass_C_g
                                        * Referenced by: '<S62>/high_pass'
                                        */
  real_T high_pass_D_j;                /* Computed Parameter: high_pass_D_j
                                        * Referenced by: '<S62>/high_pass'
                                        */
  real_T high_pass_InitialCondition_m4;/* Expression: 0
                                        * Referenced by: '<S62>/high_pass'
                                        */
  real_T high_pass_A_l[3];             /* Computed Parameter: high_pass_A_l
                                        * Referenced by: '<S63>/high_pass'
                                        */
  real_T high_pass_B_o;                /* Computed Parameter: high_pass_B_o
                                        * Referenced by: '<S63>/high_pass'
                                        */
  real_T high_pass_C_h[2];             /* Computed Parameter: high_pass_C_h
                                        * Referenced by: '<S63>/high_pass'
                                        */
  real_T high_pass_D_a;                /* Computed Parameter: high_pass_D_a
                                        * Referenced by: '<S63>/high_pass'
                                        */
  real_T high_pass_InitialCondition_gb;/* Expression: 0
                                        * Referenced by: '<S63>/high_pass'
                                        */
  real_T high_pass_A_m[3];             /* Computed Parameter: high_pass_A_m
                                        * Referenced by: '<S68>/high_pass'
                                        */
  real_T high_pass_B_hu;               /* Computed Parameter: high_pass_B_hu
                                        * Referenced by: '<S68>/high_pass'
                                        */
  real_T high_pass_C_dd[2];            /* Computed Parameter: high_pass_C_dd
                                        * Referenced by: '<S68>/high_pass'
                                        */
  real_T high_pass_D_n;                /* Computed Parameter: high_pass_D_n
                                        * Referenced by: '<S68>/high_pass'
                                        */
  real_T high_pass_InitialCondition_c; /* Expression: 0
                                        * Referenced by: '<S68>/high_pass'
                                        */
  real_T high_pass_A_f[3];             /* Computed Parameter: high_pass_A_f
                                        * Referenced by: '<S69>/high_pass'
                                        */
  real_T high_pass_B_bm;               /* Computed Parameter: high_pass_B_bm
                                        * Referenced by: '<S69>/high_pass'
                                        */
  real_T high_pass_C_dc[2];            /* Computed Parameter: high_pass_C_dc
                                        * Referenced by: '<S69>/high_pass'
                                        */
  real_T high_pass_D_mh;               /* Computed Parameter: high_pass_D_mh
                                        * Referenced by: '<S69>/high_pass'
                                        */
  real_T high_pass_InitialCondition_gw;/* Expression: 0
                                        * Referenced by: '<S69>/high_pass'
                                        */
  real_T Gain_Gain_jm;                 /* Expression: 4/32767
                                        * Referenced by: '<S99>/Gain'
                                        */
  real_T Gain_Gain_i;                  /* Expression: 4/32767
                                        * Referenced by: '<S100>/Gain'
                                        */
  real_T Gain_Gain_h0;                 /* Expression: 4/32767
                                        * Referenced by: '<S101>/Gain'
                                        */
  real_T SineWave_Amp;                 /* Expression: 100
                                        * Referenced by: '<S93>/Sine Wave'
                                        */
  real_T SineWave_Bias;                /* Expression: 0
                                        * Referenced by: '<S93>/Sine Wave'
                                        */
  real_T SineWave_Freq;                /* Expression: 1
                                        * Referenced by: '<S93>/Sine Wave'
                                        */
  real_T SineWave_Phase;               /* Expression: 0
                                        * Referenced by: '<S93>/Sine Wave'
                                        */
  real_T DataStoreMemory_InitialValue_l[4400];/* Expression: zeros(4,1100)
                                               * Referenced by: '<Root>/Data Store Memory'
                                               */
  uint32_T Delay1_DelayLength;         /* Computed Parameter: Delay1_DelayLength
                                        * Referenced by: '<S4>/Delay1'
                                        */
  uint32_T Delay3_DelayLength;         /* Computed Parameter: Delay3_DelayLength
                                        * Referenced by: '<S17>/Delay3'
                                        */
  uint32_T Delay2_DelayLength;         /* Computed Parameter: Delay2_DelayLength
                                        * Referenced by: '<S17>/Delay2'
                                        */
  uint32_T Delay7_DelayLength;         /* Computed Parameter: Delay7_DelayLength
                                        * Referenced by: '<S17>/Delay7'
                                        */
  uint32_T Delay6_DelayLength;         /* Computed Parameter: Delay6_DelayLength
                                        * Referenced by: '<S17>/Delay6'
                                        */
  uint32_T Delay1_DelayLength_e;       /* Computed Parameter: Delay1_DelayLength_e
                                        * Referenced by: '<S17>/Delay1'
                                        */
  uint32_T Delay5_DelayLength;         /* Computed Parameter: Delay5_DelayLength
                                        * Referenced by: '<S17>/Delay5'
                                        */
  uint32_T Delay4_DelayLength;         /* Computed Parameter: Delay4_DelayLength
                                        * Referenced by: '<S17>/Delay4'
                                        */
  uint32_T Delay8_DelayLength;         /* Computed Parameter: Delay8_DelayLength
                                        * Referenced by: '<S17>/Delay8'
                                        */
  uint32_T Delay_DelayLength;          /* Computed Parameter: Delay_DelayLength
                                        * Referenced by: '<S4>/Delay'
                                        */
  boolean_T VCC1_Value;                /* Computed Parameter: VCC1_Value
                                        * Referenced by: '<S46>/VCC1'
                                        */
  boolean_T VCC3_Value;                /* Computed Parameter: VCC3_Value
                                        * Referenced by: '<S46>/VCC3'
                                        */
  boolean_T Constant_Value;            /* Computed Parameter: Constant_Value
                                        * Referenced by: '<S47>/Constant'
                                        */
};

/* Real-time Model Data Structure */
struct tag_RTM_Human_in_Loop_T {
  const char_T *errorStatus;
  RTWSolverInfo solverInfo;
  X_Human_in_Loop_T *contStates;
  int_T *periodicContStateIndices;
  real_T *periodicContStateRanges;
  real_T *derivs;
  boolean_T *contStateDisabled;
  boolean_T zCCacheNeedsReset;
  boolean_T derivCacheNeedsReset;
  boolean_T CTOutputIncnstWithState;
  real_T odeF[1][48];
  ODE1_IntgData intgData;

  /*
   * Sizes:
   * The following substructure contains sizes information
   * for many of the model attributes such as inputs, outputs,
   * dwork, sample times, etc.
   */
  struct {
    int_T numContStates;
    int_T numPeriodicContStates;
    int_T numSampTimes;
  } Sizes;

  /*
   * Timing:
   * The following substructure contains information regarding
   * the timing information for the model.
   */
  struct {
    uint32_T clockTick0;
    uint32_T clockTickH0;
    time_T stepSize0;
    uint32_T clockTick1;
    uint32_T clockTickH1;
    struct {
      uint8_T TID[3];
    } TaskCounters;

    SimTimeStep simTimeStep;
    boolean_T stopRequestedFlag;
    time_T *t;
    time_T tArray[6];
  } Timing;
};

/* Block parameters (auto storage) */
extern P_Human_in_Loop_T Human_in_Loop_P;

/* Block signals (auto storage) */
extern B_Human_in_Loop_T Human_in_Loop_B;

/* Continuous states (auto storage) */
extern X_Human_in_Loop_T Human_in_Loop_X;

/* Block states (auto storage) */
extern DW_Human_in_Loop_T Human_in_Loop_DW;

/* External data declarations for dependent source files */

/* Zero-crossing (trigger) state */
extern PrevZCX_Human_in_Loop_T Human_in_Loop_PrevZCX;

/* Model entry point functions */
extern void Human_in_Loop_initialize(void);
extern void Human_in_Loop_output(void);
extern void Human_in_Loop_update(void);
extern void Human_in_Loop_terminate(void);

/* Real-time Model object */
extern RT_MODEL_Human_in_Loop_T *const Human_in_Loop_M;

/*-
 * The generated code includes comments that allow you to trace directly
 * back to the appropriate location in the model.  The basic format
 * is <system>/block_name, where system is the system number (uniquely
 * assigned by Simulink) and block_name is the name of the block.
 *
 * Use the MATLAB hilite_system command to trace the generated code back
 * to the model.  For example,
 *
 * hilite_system('<S3>')    - opens system 3
 * hilite_system('<S3>/Kp') - opens and selects block Kp which resides in S3
 *
 * Here is the system hierarchy for this model
 *
 * '<Root>' : 'Human_in_Loop'
 * '<S1>'   : 'Human_in_Loop/Control Module'
 * '<S2>'   : 'Human_in_Loop/Data Capture1'
 * '<S3>'   : 'Human_in_Loop/ObjectiveCala'
 * '<S4>'   : 'Human_in_Loop/Opt Module'
 * '<S5>'   : 'Human_in_Loop/Parameter Module'
 * '<S6>'   : 'Human_in_Loop/RTI Data'
 * '<S7>'   : 'Human_in_Loop/Sensor Data'
 * '<S8>'   : 'Human_in_Loop/State Module'
 * '<S9>'   : 'Human_in_Loop/Timer Interrupt'
 * '<S10>'  : 'Human_in_Loop/Control Module/Controller'
 * '<S11>'  : 'Human_in_Loop/Control Module/LRN'
 * '<S12>'  : 'Human_in_Loop/Control Module/Motor'
 * '<S13>'  : 'Human_in_Loop/Control Module/Torque track'
 * '<S14>'  : 'Human_in_Loop/Control Module/Motor/DAC_CLASS1_BL1'
 * '<S15>'  : 'Human_in_Loop/Control Module/Motor/MATLAB Function'
 * '<S16>'  : 'Human_in_Loop/ObjectiveCala/Multi Cycle Analysis1'
 * '<S17>'  : 'Human_in_Loop/ObjectiveCala/Single Cycle Analysis1'
 * '<S18>'  : 'Human_in_Loop/ObjectiveCala/get loss'
 * '<S19>'  : 'Human_in_Loop/ObjectiveCala/Multi Cycle Analysis1/Mean Calculate'
 * '<S20>'  : 'Human_in_Loop/ObjectiveCala/Multi Cycle Analysis1/Mean Calculate/MATLAB Function1'
 * '<S21>'  : 'Human_in_Loop/ObjectiveCala/Single Cycle Analysis1/MATLAB Function'
 * '<S22>'  : 'Human_in_Loop/Opt Module/CMAES'
 * '<S23>'  : 'Human_in_Loop/Opt Module/Cosmed'
 * '<S24>'  : 'Human_in_Loop/Opt Module/Opt Trigger'
 * '<S25>'  : 'Human_in_Loop/Opt Module/Subsystem'
 * '<S26>'  : 'Human_in_Loop/Opt Module/Cosmed/COSMED'
 * '<S27>'  : 'Human_in_Loop/Opt Module/Cosmed/Timer'
 * '<S28>'  : 'Human_in_Loop/Opt Module/Cosmed/COSMED/DS1202SER_INT_C1_I1'
 * '<S29>'  : 'Human_in_Loop/Opt Module/Cosmed/COSMED/DS1202SER_SETUP_C1'
 * '<S30>'  : 'Human_in_Loop/Opt Module/Cosmed/COSMED/Serial Decoding System'
 * '<S31>'  : 'Human_in_Loop/Opt Module/Cosmed/COSMED/Serial Decoding System/DS1202SER_RX_C1'
 * '<S32>'  : 'Human_in_Loop/Opt Module/Cosmed/COSMED/Serial Decoding System/Estimation'
 * '<S33>'  : 'Human_in_Loop/Opt Module/Cosmed/COSMED/Serial Decoding System/MATLAB Function1'
 * '<S34>'  : 'Human_in_Loop/Opt Module/Subsystem/Bayeisan Opt'
 * '<S35>'  : 'Human_in_Loop/Opt Module/Subsystem/Software Interrupt'
 * '<S36>'  : 'Human_in_Loop/Opt Module/Subsystem/Bayeisan Opt/Bayesian Function'
 * '<S37>'  : 'Human_in_Loop/Parameter Module/Control Parameter'
 * '<S38>'  : 'Human_in_Loop/Parameter Module/Torque Parameter'
 * '<S39>'  : 'Human_in_Loop/Parameter Module/Control Parameter/Mux1'
 * '<S40>'  : 'Human_in_Loop/Parameter Module/Torque Parameter/MATLAB Function'
 * '<S41>'  : 'Human_in_Loop/Parameter Module/Torque Parameter/Mux1'
 * '<S42>'  : 'Human_in_Loop/RTI Data/RTI Data Store'
 * '<S43>'  : 'Human_in_Loop/RTI Data/RTI Data Store/RTI Data Store'
 * '<S44>'  : 'Human_in_Loop/RTI Data/RTI Data Store/RTI Data Store/RTI Data Store'
 * '<S45>'  : 'Human_in_Loop/Sensor Data/EMG module'
 * '<S46>'  : 'Human_in_Loop/Sensor Data/Encoder module'
 * '<S47>'  : 'Human_in_Loop/Sensor Data/FootSwitch module'
 * '<S48>'  : 'Human_in_Loop/Sensor Data/IMU'
 * '<S49>'  : 'Human_in_Loop/Sensor Data/Torque module'
 * '<S50>'  : 'Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL1'
 * '<S51>'  : 'Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL2'
 * '<S52>'  : 'Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL3'
 * '<S53>'  : 'Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL4'
 * '<S54>'  : 'Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL5'
 * '<S55>'  : 'Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL6'
 * '<S56>'  : 'Human_in_Loop/Sensor Data/EMG module/Mux1'
 * '<S57>'  : 'Human_in_Loop/Sensor Data/EMG module/Mux3'
 * '<S58>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing'
 * '<S59>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing1'
 * '<S60>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing10'
 * '<S61>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing11'
 * '<S62>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing2'
 * '<S63>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing3'
 * '<S64>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing4'
 * '<S65>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing5'
 * '<S66>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing6'
 * '<S67>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing7'
 * '<S68>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing8'
 * '<S69>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing9'
 * '<S70>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing10/MATLAB Function'
 * '<S71>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing11/MATLAB Function'
 * '<S72>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing4/MATLAB Function'
 * '<S73>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing5/MATLAB Function'
 * '<S74>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing6/MATLAB Function'
 * '<S75>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing7/MATLAB Function'
 * '<S76>'  : 'Human_in_Loop/Sensor Data/Encoder module/1-Order TD'
 * '<S77>'  : 'Human_in_Loop/Sensor Data/Encoder module/2-Order TD'
 * '<S78>'  : 'Human_in_Loop/Sensor Data/Encoder module/Data process'
 * '<S79>'  : 'Human_in_Loop/Sensor Data/Encoder module/Data process1'
 * '<S80>'  : 'Human_in_Loop/Sensor Data/Encoder module/Discrete Derivative'
 * '<S81>'  : 'Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL1'
 * '<S82>'  : 'Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL3'
 * '<S83>'  : 'Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup1'
 * '<S84>'  : 'Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup3'
 * '<S85>'  : 'Human_in_Loop/Sensor Data/Encoder module/Mux'
 * '<S86>'  : 'Human_in_Loop/Sensor Data/Encoder module/Mux1'
 * '<S87>'  : 'Human_in_Loop/Sensor Data/Encoder module/Mux2'
 * '<S88>'  : 'Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_IN_BL1'
 * '<S89>'  : 'Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_OUT_BL1'
 * '<S90>'  : 'Human_in_Loop/Sensor Data/FootSwitch module/FootSwitch Filter'
 * '<S91>'  : 'Human_in_Loop/Sensor Data/IMU/CAN_TYPE1_SETUP_M1_C1'
 * '<S92>'  : 'Human_in_Loop/Sensor Data/IMU/CAN_TYPE1_SETUP_M1_C2'
 * '<S93>'  : 'Human_in_Loop/Sensor Data/IMU/CAN_Test'
 * '<S94>'  : 'Human_in_Loop/Sensor Data/IMU/DS1202_SENSOR_SUPPLY'
 * '<S95>'  : 'Human_in_Loop/Sensor Data/IMU/IMU'
 * '<S96>'  : 'Human_in_Loop/Sensor Data/IMU/CAN_Test/CAN_TYPE1_RX_M1_C1'
 * '<S97>'  : 'Human_in_Loop/Sensor Data/IMU/CAN_Test/CAN_TYPE1_TX_M1_C2'
 * '<S98>'  : 'Human_in_Loop/Sensor Data/IMU/IMU/CAN_RX_Accel'
 * '<S99>'  : 'Human_in_Loop/Sensor Data/IMU/IMU/Subsystem6'
 * '<S100>' : 'Human_in_Loop/Sensor Data/IMU/IMU/Subsystem7'
 * '<S101>' : 'Human_in_Loop/Sensor Data/IMU/IMU/Subsystem8'
 * '<S102>' : 'Human_in_Loop/Sensor Data/Torque module/2-Order TD'
 * '<S103>' : 'Human_in_Loop/Sensor Data/Torque module/ADC_CLASS1_BL6'
 * '<S104>' : 'Human_in_Loop/Sensor Data/Torque module/Data process'
 * '<S105>' : 'Human_in_Loop/Sensor Data/Torque module/MATLAB Function'
 * '<S106>' : 'Human_in_Loop/Sensor Data/Torque module/Mux'
 * '<S107>' : 'Human_in_Loop/State Module/Mux1'
 * '<S108>' : 'Human_in_Loop/State Module/State Machine'
 */
#endif                                 /* RTW_HEADER_Human_in_Loop_h_ */
