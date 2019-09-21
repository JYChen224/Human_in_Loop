/*
 * Human_in_Loop.h
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
#include <rtican_ds1202.h>
#include "rtwtypes.h"
#include "rtw_continuous.h"
#include "rtw_solver.h"
#endif                                 /* Human_in_Loop_COMMON_INCLUDES_ */

#include "Human_in_Loop_types.h"

/* Shared type includes */
#include "multiword_types.h"
#include "rt_zcfcn.h"
#include "rt_nonfinite.h"
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

/* Block signals for system '<S29>/Mux' */
typedef struct {
  real_T x[2];                         /* '<S29>/Mux' */
} B_Mux_Human_in_Loop_T;

/* Block signals (auto storage) */
typedef struct {
  real_T result_data[3000];
  real_T tmp_data[3000];
  real_T b_result_data[2250];
  real_T tmp_data_m[2250];
  real_T tmp_data_mb[1100];
  real_T torque_track[750];
  real_T torque_delta_track[750];
  real_T tmp_data_c[750];
  real_T tmp_data_k[750];
  real_T tmp_data_cx[750];
  real_T z1_data[750];
  real_T z1_data_b[750];
  real_T SFunction1;                   /* '<S77>/S-Function1' */
  real_T Gain;                         /* '<S32>/Gain' */
  real_T u4low2;                       /* '<S32>/0.4low2' */
  real_T x2k1;                         /* '<S76>/Unit Delay1' */
  real_T x1k1;                         /* '<S76>/Unit Delay' */
  real_T Gain1;                        /* '<S76>/Gain1' */
  real_T Gain2;                        /* '<S76>/Gain2' */
  real_T UnitDelay2;                   /* '<S76>/Unit Delay2' */
  real_T Gain4;                        /* '<S76>/Gain4' */
  real_T Add2;                         /* '<S76>/Add2' */
  real_T Gain3;                        /* '<S76>/Gain3' */
  real_T x2k;                          /* '<S76>/Add1' */
  real_T RT4[2];                       /* '<Root>/RT4' */
  real_T SFunction1_o1;                /* '<S56>/S-Function1' */
  real_T SFunction1_o2;                /* '<S56>/S-Function1' */
  real_T Gain_g;                       /* '<S29>/Gain' */
  real_T Gain1_o;                      /* '<S29>/Gain1' */
  real_T RT5[2];                       /* '<Root>/RT5' */
  real_T SFunction1_o1_j;              /* '<S57>/S-Function1' */
  real_T SFunction1_o2_j;              /* '<S57>/S-Function1' */
  real_T Gain2_c;                      /* '<S29>/Gain2' */
  real_T Gain3_h;                      /* '<S29>/Gain3' */
  real_T x2k1_d;                       /* '<S52>/Unit Delay1' */
  real_T x1k1_e;                       /* '<S52>/Unit Delay' */
  real_T Gain1_p;                      /* '<S52>/Gain1' */
  real_T Gain2_n;                      /* '<S52>/Gain2' */
  real_T UnitDelay2_i;                 /* '<S52>/Unit Delay2' */
  real_T Gain4_o;                      /* '<S52>/Gain4' */
  real_T Add2_f;                       /* '<S52>/Add2' */
  real_T Gain3_hf;                     /* '<S52>/Gain3' */
  real_T x2k_a;                        /* '<S52>/Add1' */
  real_T RT6[3];                       /* '<Root>/RT6' */
  real_T RT1[4];                       /* '<Root>/RT1' */
  real_T RT2[6];                       /* '<Root>/RT2' */
  real_T RT3[4];                       /* '<Root>/RT3' */
  real_T Gain_i;                       /* '<S3>/Gain' */
  real_T Delay3;                       /* '<S16>/Delay3' */
  real_T Delay2;                       /* '<S16>/Delay2' */
  real_T Delay7;                       /* '<S16>/Delay7' */
  real_T Delay6;                       /* '<S16>/Delay6' */
  real_T low_pass;                     /* '<S40>/low_pass' */
  real_T low_pass_e;                   /* '<S41>/low_pass' */
  real_T low_pass_m;                   /* '<S42>/low_pass' */
  real_T low_pass_d;                   /* '<S43>/low_pass' */
  real_T low_pass_p;                   /* '<S45>/low_pass' */
  real_T low_pass_i;                   /* '<S44>/low_pass' */
  real_T Delay1[6];                    /* '<S16>/Delay1' */
  real_T Delay5[6];                    /* '<S16>/Delay5' */
  real_T Delay4[6];                    /* '<S16>/Delay4' */
  real_T Delay8[6];                    /* '<S16>/Delay8' */
  real_T Gain_h;                       /* '<S76>/Gain' */
  real_T x1k;                          /* '<S76>/Add' */
  real_T TSamp;                        /* '<S55>/TSamp' */
  real_T Uk1;                          /* '<S55>/UD' */
  real_T Diff;                         /* '<S55>/Diff' */
  real_T u4low1;                       /* '<S29>/0.4low1' */
  real_T UnitDelay;                    /* '<S51>/Unit Delay' */
  real_T UnitDelay1;                   /* '<S51>/Unit Delay1' */
  real_T Add1;                         /* '<S51>/Add1' */
  real_T Gain_m;                       /* '<S51>/Gain' */
  real_T Gain1_b;                      /* '<S51>/Gain1' */
  real_T Add2_h;                       /* '<S51>/Add2' */
  real_T Add3;                         /* '<S51>/Add3' */
  real_T Gain2_i;                      /* '<S51>/Gain2' */
  real_T Gain_o;                       /* '<S52>/Gain' */
  real_T x1k_k;                        /* '<S52>/Add' */
  real_T SFunction1_g;                 /* '<S33>/S-Function1' */
  real_T SFunction1_i;                 /* '<S34>/S-Function1' */
  real_T SFunction1_m;                 /* '<S35>/S-Function1' */
  real_T SFunction1_j;                 /* '<S36>/S-Function1' */
  real_T SFunction1_l;                 /* '<S37>/S-Function1' */
  real_T SFunction1_k;                 /* '<S38>/S-Function1' */
  real_T Gain_mi;                      /* '<S28>/Gain' */
  real_T Gain1_a;                      /* '<S28>/Gain1' */
  real_T Gain2_p;                      /* '<S28>/Gain2' */
  real_T Gain3_n;                      /* '<S28>/Gain3' */
  real_T Gain4_m;                      /* '<S28>/Gain4' */
  real_T Gain5;                        /* '<S28>/Gain5' */
  real_T high_pass;                    /* '<S40>/high_pass' */
  real_T Abs;                          /* '<S40>/Abs' */
  real_T high_pass_h;                  /* '<S41>/high_pass' */
  real_T Abs_i;                        /* '<S41>/Abs' */
  real_T high_pass_hy;                 /* '<S42>/high_pass' */
  real_T Abs_l;                        /* '<S42>/Abs' */
  real_T high_pass_j;                  /* '<S43>/high_pass' */
  real_T Abs_b;                        /* '<S43>/Abs' */
  real_T high_pass_e;                  /* '<S44>/high_pass' */
  real_T Abs_lh;                       /* '<S44>/Abs' */
  real_T high_pass_b;                  /* '<S45>/high_pass' */
  real_T Abs_p;                        /* '<S45>/Abs' */
  real_T SFunction1_o1_c;              /* '<S72>/S-Function1' */
  real_T SFunction1_o2_l;              /* '<S72>/S-Function1' */
  real_T SFunction1_o3;                /* '<S72>/S-Function1' */
  real_T SFunction1_o4;                /* '<S72>/S-Function1' */
  real_T SFunction1_o5;                /* '<S72>/S-Function1' */
  real_T DataTypeConversion;           /* '<S73>/Data Type Conversion' */
  real_T Gain_gd;                      /* '<S73>/Gain' */
  real_T DataTypeConversion_j;         /* '<S74>/Data Type Conversion' */
  real_T Gain_b;                       /* '<S74>/Gain' */
  real_T DataTypeConversion_m;         /* '<S75>/Data Type Conversion' */
  real_T Gain_n;                       /* '<S75>/Gain' */
  real_T SFunction1_o1_h;              /* '<S70>/S-Function1' */
  real_T SFunction1_o2_f;              /* '<S70>/S-Function1' */
  real_T SFunction1_o3_g;              /* '<S70>/S-Function1' */
  real_T SineWave;                     /* '<S68>/Sine Wave' */
  real_T SFunction1_o1_l;              /* '<S71>/S-Function1' */
  real_T SFunction1_o2_jb;             /* '<S71>/S-Function1' */
  real_T SFunction1_o3_p;              /* '<S71>/S-Function1' */
  real_T SFunction1_o4_k;              /* '<S71>/S-Function1' */
  real_T Gain2_o;                      /* '<S11>/Gain2' */
  real_T Gain1_h;                      /* '<S11>/Gain1' */
  real_T torque_des;                   /* '<S1>/Torque track' */
  real_T torque_delta_des;             /* '<S1>/Torque track' */
  real_T torque_trace[1500];           /* '<S1>/Torque track' */
  real_T torque_delta_trace[1500];     /* '<S1>/Torque track' */
  real_T vel;                          /* '<S11>/MATLAB Function' */
  real_T lrn_cmd;                      /* '<S1>/LRN' */
  real_T lrn_mem[750];                 /* '<S1>/LRN' */
  real_T motor_vel_cmd;                /* '<S1>/Controller' */
  real_T mode;                         /* '<S7>/State Machine' */
  real_T state;                        /* '<S7>/State Machine' */
  real_T stride_time;                  /* '<S7>/State Machine' */
  real_T stride_timer;                 /* '<S7>/State Machine' */
  real_T x[4];                         /* '<S7>/Mux1' */
  real_T torque_dot;                   /* '<S32>/MATLAB Function' */
  real_T torque;                       /* '<S32>/Data process' */
  real_T state_m;                      /* '<S30>/FootSwitch Filter' */
  real_T x_i[3];                       /* '<S29>/Mux2' */
  real_T x_f[3];                       /* '<S29>/Mux1' */
  real_T angle;                        /* '<S29>/Data process1' */
  real_T angle_k;                      /* '<S29>/Data process' */
  real_T y;                            /* '<S45>/MVC' */
  real_T y_e;                          /* '<S43>/MVC' */
  real_T y_j;                          /* '<S42>/MVC' */
  real_T y_g;                          /* '<S41>/MVC' */
  real_T y_b;                          /* '<S40>/MVC' */
  real_T x_c[6];                       /* '<S28>/Mux1' */
  real_T x_k[4];                       /* '<S22>/Mux1' */
  real_T x_h[6];                       /* '<S21>/Mux1' */
  real_T torque_track_rmse;            /* '<S3>/torque_track_loss' */
  real_T torque_track_rmse_mean;       /* '<S3>/torque_track_loss' */
  real_T Cycle_Time;                   /* '<S16>/MATLAB Function' */
  real_T Cycle_Frequency;              /* '<S16>/MATLAB Function' */
  real_T Max[6];                       /* '<S16>/MATLAB Function' */
  real_T Mean[6];                      /* '<S16>/MATLAB Function' */
  real_T RMS[6];                       /* '<S16>/MATLAB Function' */
  real_T iEMG[6];                      /* '<S16>/MATLAB Function' */
  real_T UnitDelay_a;                  /* '<S18>/Unit Delay' */
  real_T UnitDelay1_k[26];             /* '<S18>/Unit Delay1' */
  real_T TmpSignalConversionAtSFunctionI[26];/* '<S18>/MATLAB Function1' */
  real_T count;                        /* '<S18>/MATLAB Function1' */
  real_T Mean_c[26];                   /* '<S18>/MATLAB Function1' */
  boolean_T SFunction1_o3_b;           /* '<S56>/S-Function1' */
  boolean_T SFunction1_it;             /* '<S63>/S-Function1' */
  B_Mux_Human_in_Loop_T sf_Mux;        /* '<S32>/Mux' */
  B_Mux_Human_in_Loop_T sf_Mux_p;      /* '<S29>/Mux' */
} B_Human_in_Loop_T;

/* Block states (auto storage) for system '<Root>' */
typedef struct {
  real_T u4low2_states[3];             /* '<S32>/0.4low2' */
  real_T UnitDelay1_DSTATE;            /* '<S76>/Unit Delay1' */
  real_T UnitDelay_DSTATE;             /* '<S76>/Unit Delay' */
  real_T UnitDelay2_DSTATE;            /* '<S76>/Unit Delay2' */
  real_T UnitDelay1_DSTATE_f;          /* '<S52>/Unit Delay1' */
  real_T UnitDelay_DSTATE_a;           /* '<S52>/Unit Delay' */
  real_T UnitDelay2_DSTATE_e;          /* '<S52>/Unit Delay2' */
  real_T Delay3_DSTATE;                /* '<S16>/Delay3' */
  real_T Delay2_DSTATE;                /* '<S16>/Delay2' */
  real_T Delay7_DSTATE;                /* '<S16>/Delay7' */
  real_T Delay6_DSTATE;                /* '<S16>/Delay6' */
  real_T Delay1_DSTATE[6];             /* '<S16>/Delay1' */
  real_T Delay5_DSTATE[6];             /* '<S16>/Delay5' */
  real_T Delay4_DSTATE[6];             /* '<S16>/Delay4' */
  real_T Delay8_DSTATE[6];             /* '<S16>/Delay8' */
  real_T UD_DSTATE;                    /* '<S55>/UD' */
  real_T u4low1_states[3];             /* '<S29>/0.4low1' */
  real_T UnitDelay_DSTATE_d;           /* '<S51>/Unit Delay' */
  real_T UnitDelay1_DSTATE_f2;         /* '<S51>/Unit Delay1' */
  real_T UnitDelay_DSTATE_g;           /* '<S18>/Unit Delay' */
  real_T UnitDelay1_DSTATE_e[26];      /* '<S18>/Unit Delay1' */
  real_T u4low2_tmp;                   /* '<S32>/0.4low2' */
  volatile real_T RT4_Buffer0[2];      /* '<Root>/RT4' */
  volatile real_T RT4_Buffer1[2];      /* '<Root>/RT4' */
  volatile real_T RT5_Buffer0[2];      /* '<Root>/RT5' */
  volatile real_T RT5_Buffer1[2];      /* '<Root>/RT5' */
  volatile real_T RT6_Buffer0[3];      /* '<Root>/RT6' */
  volatile real_T RT6_Buffer1[3];      /* '<Root>/RT6' */
  volatile real_T RT1_Buffer0[4];      /* '<Root>/RT1' */
  volatile real_T RT1_Buffer1[4];      /* '<Root>/RT1' */
  volatile real_T RT2_Buffer0[6];      /* '<Root>/RT2' */
  volatile real_T RT2_Buffer1[6];      /* '<Root>/RT2' */
  volatile real_T RT3_Buffer0[4];      /* '<Root>/RT3' */
  volatile real_T RT3_Buffer1[4];      /* '<Root>/RT3' */
  real_T EMG_Memory[24];               /* '<S16>/Data Store Memory' */
  real_T u4low1_tmp;                   /* '<S29>/0.4low1' */
  real_T TorqueMem[4400];              /* '<Root>/Data Store Memory' */
  real_T ParmReg[10];                  /* '<Root>/Data Store Memory1' */
  real_T last_footstate;               /* '<S1>/Torque track' */
  real_T last_footstate_a;             /* '<S1>/LRN' */
  real_T torque_error_memory[1000];    /* '<S1>/LRN' */
  real_T lrn_cmd_memory[1000];         /* '<S1>/LRN' */
  real_T last_torque_parm[2];          /* '<S1>/LRN' */
  real_T calib_state;                  /* '<S1>/Controller' */
  real_T reg_stride_time;              /* '<S7>/State Machine' */
  real_T reg_stride_time_count;        /* '<S7>/State Machine' */
  real_T reg_mode;                     /* '<S7>/State Machine' */
  real_T reg_state;                    /* '<S7>/State Machine' */
  real_T bt_run;                       /* '<S7>/State Machine' */
  real_T reg_last_switch;              /* '<S7>/State Machine' */
  real_T data[15];                     /* '<S32>/MATLAB Function' */
  real_T foot_state;                   /* '<S30>/FootSwitch Filter' */
  real_T filter_time;                  /* '<S30>/FootSwitch Filter' */
  real_T last_footstate_p;             /* '<S3>/torque_track_loss' */
  real_T loss_reg;                     /* '<S3>/torque_track_loss' */
  real_T loss_mem[10];                 /* '<S3>/torque_track_loss' */
  real_T SingleCycleData[2600];        /* '<S18>/Data Store Memory' */
  int_T SFunction1_IWORK[2];           /* '<S66>/S-Function1' */
  int_T SFunction1_IWORK_c[2];         /* '<S67>/S-Function1' */
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
  volatile int8_T RT3_write_buf;       /* '<Root>/RT3' */
  volatile int8_T RT3_read_buf;        /* '<Root>/RT3' */
  volatile int8_T RT3_last_buf_wr;     /* '<Root>/RT3' */
  boolean_T last_footstate_not_empty;  /* '<S1>/LRN' */
} DW_Human_in_Loop_T;

/* Continuous states (auto storage) */
typedef struct {
  real_T low_pass_CSTATE[2];           /* '<S40>/low_pass' */
  real_T low_pass_CSTATE_b[2];         /* '<S41>/low_pass' */
  real_T low_pass_CSTATE_m[2];         /* '<S42>/low_pass' */
  real_T low_pass_CSTATE_d[2];         /* '<S43>/low_pass' */
  real_T low_pass_CSTATE_j[2];         /* '<S45>/low_pass' */
  real_T low_pass_CSTATE_a[2];         /* '<S44>/low_pass' */
  real_T high_pass_CSTATE[2];          /* '<S40>/high_pass' */
  real_T high_pass_CSTATE_p[2];        /* '<S41>/high_pass' */
  real_T high_pass_CSTATE_f[2];        /* '<S42>/high_pass' */
  real_T high_pass_CSTATE_a[2];        /* '<S43>/high_pass' */
  real_T high_pass_CSTATE_o[2];        /* '<S44>/high_pass' */
  real_T high_pass_CSTATE_ff[2];       /* '<S45>/high_pass' */
} X_Human_in_Loop_T;

/* State derivatives (auto storage) */
typedef struct {
  real_T low_pass_CSTATE[2];           /* '<S40>/low_pass' */
  real_T low_pass_CSTATE_b[2];         /* '<S41>/low_pass' */
  real_T low_pass_CSTATE_m[2];         /* '<S42>/low_pass' */
  real_T low_pass_CSTATE_d[2];         /* '<S43>/low_pass' */
  real_T low_pass_CSTATE_j[2];         /* '<S45>/low_pass' */
  real_T low_pass_CSTATE_a[2];         /* '<S44>/low_pass' */
  real_T high_pass_CSTATE[2];          /* '<S40>/high_pass' */
  real_T high_pass_CSTATE_p[2];        /* '<S41>/high_pass' */
  real_T high_pass_CSTATE_f[2];        /* '<S42>/high_pass' */
  real_T high_pass_CSTATE_a[2];        /* '<S43>/high_pass' */
  real_T high_pass_CSTATE_o[2];        /* '<S44>/high_pass' */
  real_T high_pass_CSTATE_ff[2];       /* '<S45>/high_pass' */
} XDot_Human_in_Loop_T;

/* State disabled  */
typedef struct {
  boolean_T low_pass_CSTATE[2];        /* '<S40>/low_pass' */
  boolean_T low_pass_CSTATE_b[2];      /* '<S41>/low_pass' */
  boolean_T low_pass_CSTATE_m[2];      /* '<S42>/low_pass' */
  boolean_T low_pass_CSTATE_d[2];      /* '<S43>/low_pass' */
  boolean_T low_pass_CSTATE_j[2];      /* '<S45>/low_pass' */
  boolean_T low_pass_CSTATE_a[2];      /* '<S44>/low_pass' */
  boolean_T high_pass_CSTATE[2];       /* '<S40>/high_pass' */
  boolean_T high_pass_CSTATE_p[2];     /* '<S41>/high_pass' */
  boolean_T high_pass_CSTATE_f[2];     /* '<S42>/high_pass' */
  boolean_T high_pass_CSTATE_a[2];     /* '<S43>/high_pass' */
  boolean_T high_pass_CSTATE_o[2];     /* '<S44>/high_pass' */
  boolean_T high_pass_CSTATE_ff[2];    /* '<S45>/high_pass' */
} XDis_Human_in_Loop_T;

/* Zero-crossing (trigger) state */
typedef struct {
  ZCSigState MeanCalculate_Trig_ZCE;   /* '<S15>/Mean Calculate' */
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
  real_T Controller_FOLLOW_SLACK_ANGLE;/* Mask Parameter: Controller_FOLLOW_SLACK_ANGLE
                                        * Referenced by: '<S1>/Controller'
                                        */
  real_T DiscreteDerivative_ICPrevScaled;/* Mask Parameter: DiscreteDerivative_ICPrevScaled
                                          * Referenced by: '<S55>/UD'
                                          */
  real_T EMGmodule_Kemg;               /* Mask Parameter: EMGmodule_Kemg
                                        * Referenced by:
                                        *   '<S28>/Gain'
                                        *   '<S28>/Gain1'
                                        *   '<S28>/Gain2'
                                        *   '<S28>/Gain3'
                                        *   '<S28>/Gain4'
                                        *   '<S28>/Gain5'
                                        */
  real_T MATLABFunction_MAX_MOTOR_ANGLE;/* Mask Parameter: MATLABFunction_MAX_MOTOR_ANGLE
                                         * Referenced by: '<S11>/MATLAB Function'
                                         */
  real_T MATLABFunction_MAX_SPEED;     /* Mask Parameter: MATLABFunction_MAX_SPEED
                                        * Referenced by: '<S11>/MATLAB Function'
                                        */
  real_T MATLABFunction_MAX_TORQUE;    /* Mask Parameter: MATLABFunction_MAX_TORQUE
                                        * Referenced by: '<S11>/MATLAB Function'
                                        */
  real_T MATLABFunction_MIN_MOTOR_ANGLE;/* Mask Parameter: MATLABFunction_MIN_MOTOR_ANGLE
                                         * Referenced by: '<S11>/MATLAB Function'
                                         */
  real_T MultiCycleAnalysis1_N;        /* Mask Parameter: MultiCycleAnalysis1_N
                                        * Referenced by: '<S18>/Constant'
                                        */
  real_T Controller_SWING_MAX_ANGLE;   /* Mask Parameter: Controller_SWING_MAX_ANGLE
                                        * Referenced by: '<S1>/Controller'
                                        */
  real_T uOrderTD_T;                   /* Mask Parameter: uOrderTD_T
                                        * Referenced by:
                                        *   '<S51>/Gain1'
                                        *   '<S51>/Gain2'
                                        */
  real_T uOrderTD_T1;                  /* Mask Parameter: uOrderTD_T1
                                        * Referenced by:
                                        *   '<S76>/Gain1'
                                        *   '<S76>/Gain2'
                                        *   '<S76>/Gain4'
                                        */
  real_T uOrderTD_T1_p;                /* Mask Parameter: uOrderTD_T1_p
                                        * Referenced by:
                                        *   '<S52>/Gain1'
                                        *   '<S52>/Gain2'
                                        *   '<S52>/Gain4'
                                        */
  real_T uOrderTD_T2;                  /* Mask Parameter: uOrderTD_T2
                                        * Referenced by:
                                        *   '<S76>/Gain1'
                                        *   '<S76>/Gain2'
                                        *   '<S76>/Gain4'
                                        */
  real_T uOrderTD_T2_k;                /* Mask Parameter: uOrderTD_T2_k
                                        * Referenced by:
                                        *   '<S52>/Gain1'
                                        *   '<S52>/Gain2'
                                        *   '<S52>/Gain4'
                                        */
  real_T uOrderTD_Ts;                  /* Mask Parameter: uOrderTD_Ts
                                        * Referenced by:
                                        *   '<S76>/Gain'
                                        *   '<S76>/Gain3'
                                        */
  real_T uOrderTD_Ts_n;                /* Mask Parameter: uOrderTD_Ts_n
                                        * Referenced by:
                                        *   '<S52>/Gain'
                                        *   '<S52>/Gain3'
                                        */
  real_T SingleCycleAnalysis1_Ts;      /* Mask Parameter: SingleCycleAnalysis1_Ts
                                        * Referenced by: '<S16>/Sample Time'
                                        */
  real_T uOrderTD_Ts_d;                /* Mask Parameter: uOrderTD_Ts_d
                                        * Referenced by: '<S51>/Gain'
                                        */
  real_T LRN_error_filter_k;           /* Mask Parameter: LRN_error_filter_k
                                        * Referenced by: '<S1>/LRN'
                                        */
  real_T MVC_l_L_GAS_MAX;              /* Mask Parameter: MVC_l_L_GAS_MAX
                                        * Referenced by: '<S45>/MVC'
                                        */
  real_T MVC_l_L_GAS_MIN;              /* Mask Parameter: MVC_l_L_GAS_MIN
                                        * Referenced by: '<S45>/MVC'
                                        */
  real_T MVC_l_M_GAS_MAX;              /* Mask Parameter: MVC_l_M_GAS_MAX
                                        * Referenced by: '<S43>/MVC'
                                        */
  real_T MVC_l_M_GAS_MIN;              /* Mask Parameter: MVC_l_M_GAS_MIN
                                        * Referenced by: '<S43>/MVC'
                                        */
  real_T Dataprocess_load_vol_gain;    /* Mask Parameter: Dataprocess_load_vol_gain
                                        * Referenced by: '<S32>/Data process'
                                        */
  real_T Dataprocess_load_vol_offset;  /* Mask Parameter: Dataprocess_load_vol_offset
                                        * Referenced by: '<S32>/Data process'
                                        */
  real_T LRN_lrn_shrink;               /* Mask Parameter: LRN_lrn_shrink
                                        * Referenced by: '<S1>/LRN'
                                        */
  real_T MVC_r_L_GAS_MAX;              /* Mask Parameter: MVC_r_L_GAS_MAX
                                        * Referenced by: '<S42>/MVC'
                                        */
  real_T MVC_r_L_GAS_MIN;              /* Mask Parameter: MVC_r_L_GAS_MIN
                                        * Referenced by: '<S42>/MVC'
                                        */
  real_T MVC_r_M_GAS_MAX;              /* Mask Parameter: MVC_r_M_GAS_MAX
                                        * Referenced by: '<S41>/MVC'
                                        */
  real_T MVC_r_M_GAS_MIN;              /* Mask Parameter: MVC_r_M_GAS_MIN
                                        * Referenced by: '<S41>/MVC'
                                        */
  real_T MVC_r_SOL_MAX;                /* Mask Parameter: MVC_r_SOL_MAX
                                        * Referenced by: '<S40>/MVC'
                                        */
  real_T MVC_r_SOL_MIN;                /* Mask Parameter: MVC_r_SOL_MIN
                                        * Referenced by: '<S40>/MVC'
                                        */
  int32_T LRN_time_delay;              /* Mask Parameter: LRN_time_delay
                                        * Referenced by: '<S1>/LRN'
                                        */
  int8_T Controller_MODE;              /* Mask Parameter: Controller_MODE
                                        * Referenced by: '<S1>/Controller'
                                        */
  boolean_T StateMachine_BT_CALIB;     /* Mask Parameter: StateMachine_BT_CALIB
                                        * Referenced by: '<S7>/State Machine'
                                        */
  boolean_T StateMachine_BT_ERROR;     /* Mask Parameter: StateMachine_BT_ERROR
                                        * Referenced by: '<S7>/State Machine'
                                        */
  boolean_T StateMachine_BT_IDLE;      /* Mask Parameter: StateMachine_BT_IDLE
                                        * Referenced by: '<S7>/State Machine'
                                        */
  boolean_T LRN_BT_LRN_CLEAR;          /* Mask Parameter: LRN_BT_LRN_CLEAR
                                        * Referenced by: '<S1>/LRN'
                                        */
  boolean_T LRN_BT_LRN_ON;             /* Mask Parameter: LRN_BT_LRN_ON
                                        * Referenced by: '<S1>/LRN'
                                        */
  boolean_T Dataprocess_BT_RESET_ANKLE;/* Mask Parameter: Dataprocess_BT_RESET_ANKLE
                                        * Referenced by: '<S29>/Data process'
                                        */
  boolean_T Dataprocess1_BT_RESET_MOTOR;/* Mask Parameter: Dataprocess1_BT_RESET_MOTOR
                                         * Referenced by: '<S29>/Data process1'
                                         */
  boolean_T Dataprocess_BT_RESET_TORQUE;/* Mask Parameter: Dataprocess_BT_RESET_TORQUE
                                         * Referenced by: '<S32>/Data process'
                                         */
  boolean_T StateMachine_BT_RUN;       /* Mask Parameter: StateMachine_BT_RUN
                                        * Referenced by: '<S7>/State Machine'
                                        */
  boolean_T StateMachine_BT_SLACK;     /* Mask Parameter: StateMachine_BT_SLACK
                                        * Referenced by: '<S7>/State Machine'
                                        */
  real_T Mean_Y0;                      /* Computed Parameter: Mean_Y0
                                        * Referenced by: '<S18>/Mean'
                                        */
  real_T Count_Y0;                     /* Computed Parameter: Count_Y0
                                        * Referenced by: '<S18>/Count'
                                        */
  real_T UnitDelay_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S18>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S18>/Unit Delay1'
                                        */
  real_T DataStoreMemory_InitialValue[2600];/* Expression: zeros(26,100)
                                             * Referenced by: '<S18>/Data Store Memory'
                                             */
  real_T torque_offset_Value;          /* Expression: 2
                                        * Referenced by: '<S1>/torque_offset'
                                        */
  real_T Gain2_Gain;                   /* Expression: 1/360
                                        * Referenced by: '<S11>/Gain2'
                                        */
  real_T Gain1_Gain;                   /* Expression: 1/100
                                        * Referenced by: '<S11>/Gain1'
                                        */
  real_T Gain_Gain;                    /* Expression: 10
                                        * Referenced by: '<S32>/Gain'
                                        */
  real_T u4low2_NumCoef[4];            /* Expression: [0.219606211225382   0.658818633676145   0.658818633676145   0.219606211225382]*1e-3
                                        * Referenced by: '<S32>/0.4low2'
                                        */
  real_T u4low2_DenCoef[4];            /* Expression: [1.000000000000000  -2.748835809214676   2.528231219142559  -0.777638560238080]
                                        * Referenced by: '<S32>/0.4low2'
                                        */
  real_T u4low2_InitialStates;         /* Expression: 0
                                        * Referenced by: '<S32>/0.4low2'
                                        */
  real_T UnitDelay1_InitialCondition_b;/* Expression: 0
                                        * Referenced by: '<S76>/Unit Delay1'
                                        */
  real_T UnitDelay_InitialCondition_b; /* Expression: 0
                                        * Referenced by: '<S76>/Unit Delay'
                                        */
  real_T UnitDelay2_InitialCondition;  /* Expression: 0
                                        * Referenced by: '<S76>/Unit Delay2'
                                        */
  real_T RT4_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<Root>/RT4'
                                        */
  real_T Gain_Gain_o;                  /* Expression: -360/2000
                                        * Referenced by: '<S29>/Gain'
                                        */
  real_T Gain1_Gain_m;                 /* Expression: -360/2000
                                        * Referenced by: '<S29>/Gain1'
                                        */
  real_T RT5_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<Root>/RT5'
                                        */
  real_T Gain2_Gain_h;                 /* Expression: -360.0/1250.0
                                        * Referenced by: '<S29>/Gain2'
                                        */
  real_T Gain3_Gain;                   /* Expression: -360.0/1250.0
                                        * Referenced by: '<S29>/Gain3'
                                        */
  real_T UnitDelay1_InitialCondition_bo;/* Expression: 0
                                         * Referenced by: '<S52>/Unit Delay1'
                                         */
  real_T UnitDelay_InitialCondition_c; /* Expression: 0
                                        * Referenced by: '<S52>/Unit Delay'
                                        */
  real_T UnitDelay2_InitialCondition_c;/* Expression: 0
                                        * Referenced by: '<S52>/Unit Delay2'
                                        */
  real_T RT6_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<Root>/RT6'
                                        */
  real_T RT1_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<Root>/RT1'
                                        */
  real_T Kp_Value;                     /* Expression: 0
                                        * Referenced by: '<S21>/Kp'
                                        */
  real_T Kd_Value;                     /* Expression: 0
                                        * Referenced by: '<S21>/Kd'
                                        */
  real_T Kl_Value;                     /* Expression: 0
                                        * Referenced by: '<S21>/Kl'
                                        */
  real_T Ksp_Value;                    /* Expression: 1
                                        * Referenced by: '<S21>/Ksp'
                                        */
  real_T Ksd_Value;                    /* Expression: 0.1
                                        * Referenced by: '<S21>/Ksd'
                                        */
  real_T Ko_Value;                     /* Expression: 0
                                        * Referenced by: '<S21>/Ko'
                                        */
  real_T controlreset_Value;           /* Expression: 0
                                        * Referenced by: '<S21>/control reset'
                                        */
  real_T RT2_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<Root>/RT2'
                                        */
  real_T peak_torque_Value;            /* Expression: 0
                                        * Referenced by: '<S22>/peak_torque'
                                        */
  real_T rise_time_Value;              /* Expression: 22.8
                                        * Referenced by: '<S22>/rise_time'
                                        */
  real_T peak_time_Value;              /* Expression: 48
                                        * Referenced by: '<S22>/peak_time'
                                        */
  real_T fall_time_Value;              /* Expression: 11.7
                                        * Referenced by: '<S22>/fall_time'
                                        */
  real_T RT3_InitialCondition;         /* Expression: 0
                                        * Referenced by: '<Root>/RT3'
                                        */
  real_T Gain_Gain_j;                  /* Expression: 5000
                                        * Referenced by: '<S3>/Gain'
                                        */
  real_T Delay3_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S16>/Delay3'
                                        */
  real_T Delay2_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S16>/Delay2'
                                        */
  real_T Delay7_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S16>/Delay7'
                                        */
  real_T Delay6_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S16>/Delay6'
                                        */
  real_T low_pass_A[3];                /* Computed Parameter: low_pass_A
                                        * Referenced by: '<S40>/low_pass'
                                        */
  real_T low_pass_B;                   /* Computed Parameter: low_pass_B
                                        * Referenced by: '<S40>/low_pass'
                                        */
  real_T low_pass_C;                   /* Computed Parameter: low_pass_C
                                        * Referenced by: '<S40>/low_pass'
                                        */
  real_T low_pass_InitialCondition;    /* Expression: 0
                                        * Referenced by: '<S40>/low_pass'
                                        */
  real_T low_pass_A_k[3];              /* Computed Parameter: low_pass_A_k
                                        * Referenced by: '<S41>/low_pass'
                                        */
  real_T low_pass_B_j;                 /* Computed Parameter: low_pass_B_j
                                        * Referenced by: '<S41>/low_pass'
                                        */
  real_T low_pass_C_a;                 /* Computed Parameter: low_pass_C_a
                                        * Referenced by: '<S41>/low_pass'
                                        */
  real_T low_pass_InitialCondition_i;  /* Expression: 0
                                        * Referenced by: '<S41>/low_pass'
                                        */
  real_T low_pass_A_e[3];              /* Computed Parameter: low_pass_A_e
                                        * Referenced by: '<S42>/low_pass'
                                        */
  real_T low_pass_B_f;                 /* Computed Parameter: low_pass_B_f
                                        * Referenced by: '<S42>/low_pass'
                                        */
  real_T low_pass_C_g;                 /* Computed Parameter: low_pass_C_g
                                        * Referenced by: '<S42>/low_pass'
                                        */
  real_T low_pass_InitialCondition_n;  /* Expression: 0
                                        * Referenced by: '<S42>/low_pass'
                                        */
  real_T low_pass_A_l[3];              /* Computed Parameter: low_pass_A_l
                                        * Referenced by: '<S43>/low_pass'
                                        */
  real_T low_pass_B_n;                 /* Computed Parameter: low_pass_B_n
                                        * Referenced by: '<S43>/low_pass'
                                        */
  real_T low_pass_C_e;                 /* Computed Parameter: low_pass_C_e
                                        * Referenced by: '<S43>/low_pass'
                                        */
  real_T low_pass_InitialCondition_o;  /* Expression: 0
                                        * Referenced by: '<S43>/low_pass'
                                        */
  real_T low_pass_A_c[3];              /* Computed Parameter: low_pass_A_c
                                        * Referenced by: '<S45>/low_pass'
                                        */
  real_T low_pass_B_m;                 /* Computed Parameter: low_pass_B_m
                                        * Referenced by: '<S45>/low_pass'
                                        */
  real_T low_pass_C_o;                 /* Computed Parameter: low_pass_C_o
                                        * Referenced by: '<S45>/low_pass'
                                        */
  real_T low_pass_InitialCondition_ir; /* Expression: 0
                                        * Referenced by: '<S45>/low_pass'
                                        */
  real_T low_pass_A_ec[3];             /* Computed Parameter: low_pass_A_ec
                                        * Referenced by: '<S44>/low_pass'
                                        */
  real_T low_pass_B_h;                 /* Computed Parameter: low_pass_B_h
                                        * Referenced by: '<S44>/low_pass'
                                        */
  real_T low_pass_C_p;                 /* Computed Parameter: low_pass_C_p
                                        * Referenced by: '<S44>/low_pass'
                                        */
  real_T low_pass_InitialCondition_h;  /* Expression: 0
                                        * Referenced by: '<S44>/low_pass'
                                        */
  real_T Delay1_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S16>/Delay1'
                                        */
  real_T Delay5_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S16>/Delay5'
                                        */
  real_T Delay4_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S16>/Delay4'
                                        */
  real_T Delay8_InitialCondition;      /* Expression: 0
                                        * Referenced by: '<S16>/Delay8'
                                        */
  real_T DataStoreMemory_InitialValue_o[24];/* Expression: zeros(6,4)
                                             * Referenced by: '<S16>/Data Store Memory'
                                             */
  real_T TSamp_WtEt;                   /* Computed Parameter: TSamp_WtEt
                                        * Referenced by: '<S55>/TSamp'
                                        */
  real_T u4low1_NumCoef[4];            /* Expression: [   0.001567010350588   0.004701031051765   0.004701031051765   0.001567010350588
                                          ]
                                        * Referenced by: '<S29>/0.4low1'
                                        */
  real_T u4low1_DenCoef[4];            /* Expression: [   1.000000000000000  -2.498608344691178   2.115254127003159  -0.604109699507275
                                          ]
                                        * Referenced by: '<S29>/0.4low1'
                                        */
  real_T u4low1_InitialStates;         /* Expression: 0
                                        * Referenced by: '<S29>/0.4low1'
                                        */
  real_T UnitDelay_InitialCondition_i; /* Expression: 0
                                        * Referenced by: '<S51>/Unit Delay'
                                        */
  real_T UnitDelay1_InitialCondition_e;/* Expression: 0
                                        * Referenced by: '<S51>/Unit Delay1'
                                        */
  real_T high_pass_A[3];               /* Computed Parameter: high_pass_A
                                        * Referenced by: '<S40>/high_pass'
                                        */
  real_T high_pass_B;                  /* Computed Parameter: high_pass_B
                                        * Referenced by: '<S40>/high_pass'
                                        */
  real_T high_pass_C[2];               /* Computed Parameter: high_pass_C
                                        * Referenced by: '<S40>/high_pass'
                                        */
  real_T high_pass_D;                  /* Computed Parameter: high_pass_D
                                        * Referenced by: '<S40>/high_pass'
                                        */
  real_T high_pass_InitialCondition;   /* Expression: 0
                                        * Referenced by: '<S40>/high_pass'
                                        */
  real_T high_pass_A_g[3];             /* Computed Parameter: high_pass_A_g
                                        * Referenced by: '<S41>/high_pass'
                                        */
  real_T high_pass_B_i;                /* Computed Parameter: high_pass_B_i
                                        * Referenced by: '<S41>/high_pass'
                                        */
  real_T high_pass_C_b[2];             /* Computed Parameter: high_pass_C_b
                                        * Referenced by: '<S41>/high_pass'
                                        */
  real_T high_pass_D_d;                /* Computed Parameter: high_pass_D_d
                                        * Referenced by: '<S41>/high_pass'
                                        */
  real_T high_pass_InitialCondition_a; /* Expression: 0
                                        * Referenced by: '<S41>/high_pass'
                                        */
  real_T high_pass_A_p[3];             /* Computed Parameter: high_pass_A_p
                                        * Referenced by: '<S42>/high_pass'
                                        */
  real_T high_pass_B_h;                /* Computed Parameter: high_pass_B_h
                                        * Referenced by: '<S42>/high_pass'
                                        */
  real_T high_pass_C_by[2];            /* Computed Parameter: high_pass_C_by
                                        * Referenced by: '<S42>/high_pass'
                                        */
  real_T high_pass_D_o;                /* Computed Parameter: high_pass_D_o
                                        * Referenced by: '<S42>/high_pass'
                                        */
  real_T high_pass_InitialCondition_d; /* Expression: 0
                                        * Referenced by: '<S42>/high_pass'
                                        */
  real_T high_pass_A_a[3];             /* Computed Parameter: high_pass_A_a
                                        * Referenced by: '<S43>/high_pass'
                                        */
  real_T high_pass_B_n;                /* Computed Parameter: high_pass_B_n
                                        * Referenced by: '<S43>/high_pass'
                                        */
  real_T high_pass_C_p[2];             /* Computed Parameter: high_pass_C_p
                                        * Referenced by: '<S43>/high_pass'
                                        */
  real_T high_pass_D_n;                /* Computed Parameter: high_pass_D_n
                                        * Referenced by: '<S43>/high_pass'
                                        */
  real_T high_pass_InitialCondition_m; /* Expression: 0
                                        * Referenced by: '<S43>/high_pass'
                                        */
  real_T high_pass_A_a0[3];            /* Computed Parameter: high_pass_A_a0
                                        * Referenced by: '<S44>/high_pass'
                                        */
  real_T high_pass_B_n3;               /* Computed Parameter: high_pass_B_n3
                                        * Referenced by: '<S44>/high_pass'
                                        */
  real_T high_pass_C_e[2];             /* Computed Parameter: high_pass_C_e
                                        * Referenced by: '<S44>/high_pass'
                                        */
  real_T high_pass_D_l;                /* Computed Parameter: high_pass_D_l
                                        * Referenced by: '<S44>/high_pass'
                                        */
  real_T high_pass_InitialCondition_h; /* Expression: 0
                                        * Referenced by: '<S44>/high_pass'
                                        */
  real_T high_pass_A_f[3];             /* Computed Parameter: high_pass_A_f
                                        * Referenced by: '<S45>/high_pass'
                                        */
  real_T high_pass_B_g;                /* Computed Parameter: high_pass_B_g
                                        * Referenced by: '<S45>/high_pass'
                                        */
  real_T high_pass_C_h[2];             /* Computed Parameter: high_pass_C_h
                                        * Referenced by: '<S45>/high_pass'
                                        */
  real_T high_pass_D_e;                /* Computed Parameter: high_pass_D_e
                                        * Referenced by: '<S45>/high_pass'
                                        */
  real_T high_pass_InitialCondition_e; /* Expression: 0
                                        * Referenced by: '<S45>/high_pass'
                                        */
  real_T Gain_Gain_a;                  /* Expression: 4/32767
                                        * Referenced by: '<S73>/Gain'
                                        */
  real_T Gain_Gain_h;                  /* Expression: 4/32767
                                        * Referenced by: '<S74>/Gain'
                                        */
  real_T Gain_Gain_jd;                 /* Expression: 4/32767
                                        * Referenced by: '<S75>/Gain'
                                        */
  real_T SineWave_Amp;                 /* Expression: 100
                                        * Referenced by: '<S68>/Sine Wave'
                                        */
  real_T SineWave_Bias;                /* Expression: 0
                                        * Referenced by: '<S68>/Sine Wave'
                                        */
  real_T SineWave_Freq;                /* Expression: 1
                                        * Referenced by: '<S68>/Sine Wave'
                                        */
  real_T SineWave_Phase;               /* Expression: 0
                                        * Referenced by: '<S68>/Sine Wave'
                                        */
  real_T DataStoreMemory_InitialValue_l[4400];/* Expression: zeros(4,1100)
                                               * Referenced by: '<Root>/Data Store Memory'
                                               */
  real_T DataStoreMemory1_InitialValue[10];/* Expression: zeros(1,10)
                                            * Referenced by: '<Root>/Data Store Memory1'
                                            */
  uint32_T Delay3_DelayLength;         /* Computed Parameter: Delay3_DelayLength
                                        * Referenced by: '<S16>/Delay3'
                                        */
  uint32_T Delay2_DelayLength;         /* Computed Parameter: Delay2_DelayLength
                                        * Referenced by: '<S16>/Delay2'
                                        */
  uint32_T Delay7_DelayLength;         /* Computed Parameter: Delay7_DelayLength
                                        * Referenced by: '<S16>/Delay7'
                                        */
  uint32_T Delay6_DelayLength;         /* Computed Parameter: Delay6_DelayLength
                                        * Referenced by: '<S16>/Delay6'
                                        */
  uint32_T Delay1_DelayLength;         /* Computed Parameter: Delay1_DelayLength
                                        * Referenced by: '<S16>/Delay1'
                                        */
  uint32_T Delay5_DelayLength;         /* Computed Parameter: Delay5_DelayLength
                                        * Referenced by: '<S16>/Delay5'
                                        */
  uint32_T Delay4_DelayLength;         /* Computed Parameter: Delay4_DelayLength
                                        * Referenced by: '<S16>/Delay4'
                                        */
  uint32_T Delay8_DelayLength;         /* Computed Parameter: Delay8_DelayLength
                                        * Referenced by: '<S16>/Delay8'
                                        */
  boolean_T VCC1_Value;                /* Computed Parameter: VCC1_Value
                                        * Referenced by: '<S29>/VCC1'
                                        */
  boolean_T VCC3_Value;                /* Computed Parameter: VCC3_Value
                                        * Referenced by: '<S29>/VCC3'
                                        */
  boolean_T Constant_Value;            /* Computed Parameter: Constant_Value
                                        * Referenced by: '<S30>/Constant'
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
  real_T odeF[1][24];
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
    time_T tArray[4];
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
 * '<S4>'   : 'Human_in_Loop/Parameter Module'
 * '<S5>'   : 'Human_in_Loop/RTI Data'
 * '<S6>'   : 'Human_in_Loop/Sensor Data'
 * '<S7>'   : 'Human_in_Loop/State Module'
 * '<S8>'   : 'Human_in_Loop/Timer Interrupt'
 * '<S9>'   : 'Human_in_Loop/Control Module/Controller'
 * '<S10>'  : 'Human_in_Loop/Control Module/LRN'
 * '<S11>'  : 'Human_in_Loop/Control Module/Motor'
 * '<S12>'  : 'Human_in_Loop/Control Module/Torque track'
 * '<S13>'  : 'Human_in_Loop/Control Module/Motor/DAC_CLASS1_BL1'
 * '<S14>'  : 'Human_in_Loop/Control Module/Motor/MATLAB Function'
 * '<S15>'  : 'Human_in_Loop/ObjectiveCala/Multi Cycle Analysis1'
 * '<S16>'  : 'Human_in_Loop/ObjectiveCala/Single Cycle Analysis1'
 * '<S17>'  : 'Human_in_Loop/ObjectiveCala/torque_track_loss'
 * '<S18>'  : 'Human_in_Loop/ObjectiveCala/Multi Cycle Analysis1/Mean Calculate'
 * '<S19>'  : 'Human_in_Loop/ObjectiveCala/Multi Cycle Analysis1/Mean Calculate/MATLAB Function1'
 * '<S20>'  : 'Human_in_Loop/ObjectiveCala/Single Cycle Analysis1/MATLAB Function'
 * '<S21>'  : 'Human_in_Loop/Parameter Module/Control Parameter'
 * '<S22>'  : 'Human_in_Loop/Parameter Module/Torque Parameter'
 * '<S23>'  : 'Human_in_Loop/Parameter Module/Control Parameter/Mux1'
 * '<S24>'  : 'Human_in_Loop/Parameter Module/Torque Parameter/Mux1'
 * '<S25>'  : 'Human_in_Loop/RTI Data/RTI Data Store'
 * '<S26>'  : 'Human_in_Loop/RTI Data/RTI Data Store/RTI Data Store'
 * '<S27>'  : 'Human_in_Loop/RTI Data/RTI Data Store/RTI Data Store/RTI Data Store'
 * '<S28>'  : 'Human_in_Loop/Sensor Data/EMG module'
 * '<S29>'  : 'Human_in_Loop/Sensor Data/Encoder module'
 * '<S30>'  : 'Human_in_Loop/Sensor Data/FootSwitch module'
 * '<S31>'  : 'Human_in_Loop/Sensor Data/IMU'
 * '<S32>'  : 'Human_in_Loop/Sensor Data/Torque module'
 * '<S33>'  : 'Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL1'
 * '<S34>'  : 'Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL2'
 * '<S35>'  : 'Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL3'
 * '<S36>'  : 'Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL4'
 * '<S37>'  : 'Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL5'
 * '<S38>'  : 'Human_in_Loop/Sensor Data/EMG module/ADC_CLASS1_BL6'
 * '<S39>'  : 'Human_in_Loop/Sensor Data/EMG module/Mux1'
 * '<S40>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing'
 * '<S41>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing1'
 * '<S42>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing2'
 * '<S43>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing3'
 * '<S44>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing8'
 * '<S45>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing9'
 * '<S46>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing/MVC'
 * '<S47>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing1/MVC'
 * '<S48>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing2/MVC'
 * '<S49>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing3/MVC'
 * '<S50>'  : 'Human_in_Loop/Sensor Data/EMG module/Preprocessing9/MVC'
 * '<S51>'  : 'Human_in_Loop/Sensor Data/Encoder module/1-Order TD'
 * '<S52>'  : 'Human_in_Loop/Sensor Data/Encoder module/2-Order TD'
 * '<S53>'  : 'Human_in_Loop/Sensor Data/Encoder module/Data process'
 * '<S54>'  : 'Human_in_Loop/Sensor Data/Encoder module/Data process1'
 * '<S55>'  : 'Human_in_Loop/Sensor Data/Encoder module/Discrete Derivative'
 * '<S56>'  : 'Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL1'
 * '<S57>'  : 'Human_in_Loop/Sensor Data/Encoder module/EMC_ENCODER_BL3'
 * '<S58>'  : 'Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup1'
 * '<S59>'  : 'Human_in_Loop/Sensor Data/Encoder module/Encoder Power Setup3'
 * '<S60>'  : 'Human_in_Loop/Sensor Data/Encoder module/Mux'
 * '<S61>'  : 'Human_in_Loop/Sensor Data/Encoder module/Mux1'
 * '<S62>'  : 'Human_in_Loop/Sensor Data/Encoder module/Mux2'
 * '<S63>'  : 'Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_IN_BL1'
 * '<S64>'  : 'Human_in_Loop/Sensor Data/FootSwitch module/DIO_CLASS1_BIT_OUT_BL1'
 * '<S65>'  : 'Human_in_Loop/Sensor Data/FootSwitch module/FootSwitch Filter'
 * '<S66>'  : 'Human_in_Loop/Sensor Data/IMU/CAN_TYPE1_SETUP_M1_C1'
 * '<S67>'  : 'Human_in_Loop/Sensor Data/IMU/CAN_TYPE1_SETUP_M1_C2'
 * '<S68>'  : 'Human_in_Loop/Sensor Data/IMU/CAN_Test'
 * '<S69>'  : 'Human_in_Loop/Sensor Data/IMU/IMU'
 * '<S70>'  : 'Human_in_Loop/Sensor Data/IMU/CAN_Test/CAN_TYPE1_RX_M1_C1'
 * '<S71>'  : 'Human_in_Loop/Sensor Data/IMU/CAN_Test/CAN_TYPE1_TX_M1_C2'
 * '<S72>'  : 'Human_in_Loop/Sensor Data/IMU/IMU/CAN_RX_Accel'
 * '<S73>'  : 'Human_in_Loop/Sensor Data/IMU/IMU/Subsystem6'
 * '<S74>'  : 'Human_in_Loop/Sensor Data/IMU/IMU/Subsystem7'
 * '<S75>'  : 'Human_in_Loop/Sensor Data/IMU/IMU/Subsystem8'
 * '<S76>'  : 'Human_in_Loop/Sensor Data/Torque module/2-Order TD'
 * '<S77>'  : 'Human_in_Loop/Sensor Data/Torque module/ADC_CLASS1_BL6'
 * '<S78>'  : 'Human_in_Loop/Sensor Data/Torque module/DS1202_SENSOR_SUPPLY'
 * '<S79>'  : 'Human_in_Loop/Sensor Data/Torque module/Data process'
 * '<S80>'  : 'Human_in_Loop/Sensor Data/Torque module/MATLAB Function'
 * '<S81>'  : 'Human_in_Loop/Sensor Data/Torque module/Mux'
 * '<S82>'  : 'Human_in_Loop/State Module/Mux1'
 * '<S83>'  : 'Human_in_Loop/State Module/State Machine'
 */
#endif                                 /* RTW_HEADER_Human_in_Loop_h_ */
